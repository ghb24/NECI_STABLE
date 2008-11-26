#ifdef PARALLEL
!This is a parallel MPI version of the FciMC code.
!All variables refer to values per processor

!!!!! TO DO !!!!
! Add new test
! Option for calculating non-essential information
! Simulate Magnetization to break sign symmetry

MODULE FciMCParMod
    use SystemData , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,nMsh,Arr
    use CalcData , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
    use CalcData , only : TStartMP1,NEquilSteps,TReadPops,TRegenExcitgens,TFixShiftShell,ShellFix,FixShift
    use CalcData , only : GrowMaxFactor,CullFactor,TStartSinglePart,ScaleWalkers,Lambda,TLocalAnnihilation
    use CalcData , only : NDets,RhoApp,TResumFCIMC,TNoAnnihil,MemoryFac,TAnnihilonproc,MemoryFacExcit
    USE Determinants , only : FDet,GetHElement2,GetHElement4
    USE DetCalc , only : NMRKS,ICILevel
    use IntegralsData , only : fck,NMax,UMat
    USE Logging , only : iWritePopsEvery,TPopsFile,TZeroProjE
    USE Logging , only : TAutoCorr,NoACDets!,iLagMin,iLagMax,iLagStep
    USE global_utilities
    USE HElem
    USE Parallel
    IMPLICIT NONE
    SAVE

    INTEGER , PARAMETER :: Root=0   !This is the rank of the root processor
    INTEGER , PARAMETER :: r2=kind(0.d0)
    INTEGER , PARAMETER :: i2=SELECTED_INT_KIND(18)
    
    TYPE ExcitGenerator
        INTEGER , POINTER :: ExcitData(:)=>null()   !This stores the excitation generator
        INTEGER :: nExcitMemLen                     !This is the length of the excitation generator
        INTEGER :: nPointed                         !This indicates the number of elements in the excitation pointer arrays which are pointing to this position
    END TYPE

    TYPE ExcitPointer
        INTEGER , POINTER :: PointToExcit(:)        !This is a pointer to the excitation generator in ExcitGens
        INTEGER :: IndexinExArr                     !This is the index in ExcitGens which we are pointing at
    END TYPE

    TYPE(ExcitGenerator) , ALLOCATABLE :: ExcitGens(:)  !This is the array to store all the excitation generators
    INTEGER , ALLOCATABLE :: FreeIndArray(:)            !This is a circular list of the free positions to put excitation generators in ExcitGens.
    INTEGER :: BackofList,FrontOfList                   !This indicates where in the list we are.
                                                        !We add indices which become available in ExcitGens to the front of the list, and when we use them up,
                                                        !take off the back of the list. It is circular and so will repeat indefinitly.

    TYPE(ExcitPointer) , ALLOCATABLE , TARGET :: WalkVecExcits(:),WalkVec2Excits(:)   !This will store the excitation generators for the particles on each node

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
    TYPE(ExcitPointer) , POINTER :: CurrentExcits(:), NewExcits(:)
    
    INTEGER , ALLOCATABLE :: HFDet(:)       !This will store the HF determinant
    INTEGER :: HFDetTag=0
    TYPE(ExcitGenerator) :: HFExcit         !This is the excitation generator for the HF determinant
    INTEGER(KIND=i2) :: HFHash               !This is the hash for the HF determinant

    INTEGER :: Seed,MaxWalkers,TotWalkers,TotWalkersOld,PreviousNMCyc,Iter,NoComps,MaxWalkersExcit
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
    REAL*8 :: AvSign           !This is the average sign of the particles on each node
    REAL*8 :: AvSignHFD        !This is the average sign of the particles at HF or Double excitations on each node
    INTEGER :: SumWalkersCyc    !This is the sum of all walkers over an update cycle on each processor
    REAL*8 :: MeanExcitLevel    
    INTEGER :: MinExcitLevel
    INTEGER :: MaxExcitLevel
    INTEGER :: Annihilated      !This is the number annihilated on one processor
    INTEGER :: LocalAnn         !This is the number of locally annihilated particles
    INTEGER :: NoatHF           !This is the instantaneous number of particles at the HF determinant
    INTEGER :: NoatDoubs
    INTEGER :: Acceptances      !This is the number of accepted spawns - this is only calculated per node.
    REAL*8 :: AccRat            !Acceptance ratio for each node over the update cycle
    INTEGER :: PreviousCycles   !This is just for the head node, so that it can store the number of previous cycles when reading from POPSFILE
    INTEGER :: NoBorn,NoDied

!These are the global variables, calculated on the root processor, from the values above
    REAL*8 :: AllGrowRate,AllMeanExcitLevel
    INTEGER :: AllMinExcitLevel,AllMaxExcitLevel,AllTotWalkers,AllTotWalkersOld,AllSumWalkersCyc,AllTotSign,AllTotSignOld
    INTEGER :: AllAnnihilated,AllNoatHF,AllNoatDoubs,AllLocalAnn
    REAL*8 :: AllSumNoatHF,AllSumENum,AllAvSign,AllAvSignHFD
    INTEGER :: AllNoBorn,AllNoDied

    REAL*8 :: MPNorm        !MPNorm is used if TNodalCutoff is set, to indicate the normalisation of the MP Wavevector

    TYPE(HElement) :: rhii
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

!    INTEGER :: mpilongintegertype               !This is used to create an MPI derived type to cope with 8 byte integers

    LOGICAL :: TBalanceNodes                    !This is true when the nodes need to be balanced

    LOGICAL :: TDebug                           !Debugging flag
    INTEGER :: MaxIndex

    LOGICAL :: TTruncSpace=.false.              !This is a flag set as to whether the excitation space should be truncated or not.
    LOGICAL :: TFlippedSign=.false.             !This is to indicate when the sign of the particles have been flipped. This is needed for the calculation of the ACF

    TYPE(timer), save :: Walker_Time,Annihil_Time,ACF_Time

!These are arrays used to store the autocorrelation function
    INTEGER , ALLOCATABLE :: WeightatDets(:)                   !First index - det which is stored, second - weight on proc at that iteration
    INTEGER , ALLOCATABLE :: AutoCorrDets(:,:)                  !(NEl,NoAutoDets)
    INTEGER :: WeightatDetsTag=0,AutoCorrDetsTag=0
    INTEGER :: NoAutoDets

!These are variables relating to calculating the population density of an excitation level for use with TLocalAnnihilation
!Local annihilation is currently commented out
    REAL*8 , ALLOCATABLE :: ApproxExcitDets(:)    !This is the approximate size of each excitation level
    INTEGER , ALLOCATABLE :: PartsInExcitLevel(:) !This is the number of walkers for a given iteration in that excitation level
        
    contains


    SUBROUTINE FciMCPar(Weight,Energyxw)
        use soft_exit, only : test_SOFTEXIT
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: i,j,error
        CHARACTER(len=*), PARAMETER :: this_routine='FciMCPar'
        TYPE(HElement) :: Hamii
        LOGICAL :: TIncrement

        TDebug=.false.  !Set debugging flag

!OpenMPI does not currently support MPI_Comm_set_errhandler - a bug in its F90 interface code.
!Ask Nick McLaren if we need to change the err handler - he has a fix/bypass.
!        CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN,error)
!        CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_ARE_FATAL,error)
        
        CALL InitFCIMCCalcPar()

        CALL WriteFCIMCStats()

!Start MC simulation...
        TIncrement=.true.   !If TIncrement is true, it means that when it comes out of the loop, it wants to subtract 1 from the Iteration count to get the true number of iterations
        do Iter=1,NMCyc

            IF(TBalanceNodes) CALL BalanceWalkersonProcs()      !This routine is call periodically when the nodes need to be balanced.
            
            CALL PerformFCIMCycPar()

            IF(mod(Iter,StepsSft).eq.0) THEN
!This will communicate between all nodes, find the new shift (and other parameters) and broadcast them to the other nodes.
                CALL CalcNewShift()

!Test if the file SOFTEXIT is present. If it is, then exit cleanly.
                IF(test_SOFTEXIT()) THEN
                    TIncrement=.false.  !This is so that when the popsfile is written out before exit, it counteracts the subtraction by one.
                    EXIT
                ENDIF

            ENDIF

            IF(TPopsFile.and.(mod(Iter,iWritePopsEvery).eq.0)) THEN
!This will write out the POPSFILE
                CALL WriteToPopsFilePar()
            ENDIF
            IF(TAutoCorr) CALL WriteHistogrammedDets

!End of MC cycle
        enddo

        IF(TIncrement) Iter=Iter-1     !Reduce the iteration count for the POPSFILE since it is incremented upon leaving the loop (if done naturally)
        IF(TPopsFile) CALL WriteToPopsFilePar()

        Weight=HDElement(0.D0)
        Energyxw=HDElement(ProjectionE)

!Deallocate memory
        CALL DeallocFCIMCMemPar()
        
!        IF(TAutoCorr) CALL CalcAutoCorr()

        IF(iProcIndex.eq.Root) THEN
            CLOSE(15)
            IF(TAutoCorr) CLOSE(44)
        ENDIF
        IF(TDebug) CLOSE(11)

!        CALL MPIEnd(.false.)    !Finalize MPI

        RETURN

    END SUBROUTINE FciMCPar

!This is the heart of FCIMC, where the MC Cycles are performed
    SUBROUTINE PerformFCIMCycPar()
        INTEGER :: VecSlot,i,j,k,l
        INTEGER :: nJ(NEl),ierr,IC,Child,iCount
        REAL*8 :: Prob,rat,HDiag
        INTEGER :: iDie             !Indicated whether a particle should self-destruct on DetCurr
        INTEGER :: ExcitLevel,iGetExcitLevel_2,TotWalkersNew,error,length,temp,Ex(2,2)
        LOGICAL :: WSign,tParity
        INTEGER(KIND=i2) :: HashTemp
        TYPE(HElement) :: HDiagTemp,HOffDiag
        CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: message

        IF(TDebug) THEN
            WRITE(11,*) Iter,TotWalkers,NoatHF,NoatDoubs,MaxIndex
            CALL FLUSH(11)
        ENDIF

        CALL set_timer(Walker_Time,30)
        
!VecSlot indicates the next free position in NewDets
        VecSlot=1
!Reset number at HF and doubles
        NoatHF=0
        NoatDoubs=0
        
        do j=1,TotWalkers
!j runs through all current walkers

!            WRITE(6,*) Iter,j,TotWalkers
!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
            CALL SumEContrib(CurrentDets(:,j),CurrentIC(j),CurrentSign(j))

            IF(TResumFCIMC) THEN

                CALL ResumGraphPar(CurrentDets(:,j),CurrentSign(j),VecSlot,j)

            ELSE

                IF(.not.TRegenExcitgens) THEN
!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
                    CALL SetupExitgenPar(CurrentDets(:,j),CurrentExcits(j))
                    CALL GenRandSymExcitIt4(CurrentDets(:,j),CurrentExcits(j)%PointToExcit,nJ,Seed,IC,0,Prob,iCount,Ex,tParity)
                ELSE
                    CALL GetPartRandExcitPar(CurrentDets(:,j),nJ,Seed,IC,0,Prob,iCount,CurrentIC(j),Ex,tParity)
                ENDIF
!Calculate number of children to spawn

                IF(TTruncSpace) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.

                   IF(CurrentIC(j).eq.(ICILevel-1)) THEN
!The current walker is one below the excitation cutoff - if IC is a double, then could go over

                        IF(IC.eq.2) THEN
!Need to check excitation level of excitation
                            ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,ICILevel)
                            IF(ExcitLevel.gt.ICILevel) THEN
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                                Child=0
                            ELSE
                                Child=AttemptCreatePar(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC,Ex,tParity)
                            ENDIF
                        ELSE
                            Child=AttemptCreatePar(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC,Ex,tParity)
                        ENDIF

                    ELSEIF(CurrentIC(j).eq.ICILevel) THEN
!Walker is at the excitation cutoff level - all possible excitations could be disallowed - check the actual excitation level
                        ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,ICILevel)
                        IF(ExcitLevel.gt.ICILevel) THEN
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                            Child=0
                        ELSE
                            Child=AttemptCreatePar(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC,Ex,tParity)
                        ENDIF
                    ELSE
!Excitation cannot be in a dissallowed excitation level - allow it as normal
                        Child=AttemptCreatePar(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC,Ex,tParity)
                    ENDIF 
                
                ELSE
!SD Space is not truncated - allow attempted spawn as usual

                    Child=AttemptCreatePar(CurrentDets(:,j),CurrentSign(j),nJ,Prob,IC,Ex,tParity)

                ENDIF
                
                IF(Child.ne.0) THEN
!We want to spawn a child - find its information to store

                    NoBorn=NoBorn+abs(Child)     !Update counter about particle birth

                    IF(Child.gt.25) THEN
                        WRITE(6,*) "LARGE PARTICLE BLOOM - ",Child," particles created in one attempt."
                        WRITE(6,*) "BEWARE OF MEMORY PROBLEMS"
!                        WRITE(6,*) "DETERMINANT CREATED IS: ",nJ
                        WRITE(6,*) "PROB IS: ",Prob
                        CALL FLUSH(6)
                    ENDIF

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

                    IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                        HashTemp=CreateHash(nJ)
                    ENDIF

                    do l=1,abs(Child)
!Copy across children - cannot copy excitation generators, as do not know them...unless it is HF (this might save a little time if implimented)
                        NewDets(:,VecSlot)=nJ(:)
                        NewSign(VecSlot)=WSign
                        IF(.not.TRegenExcitgens) NewExcits(VecSlot)%PointToExcit=>null()
                        NewIC(VecSlot)=ExcitLevel
                        NewH(VecSlot)=HDiag                     !Diagonal H-element-Hii
                        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                            Hash2Array(VecSlot)=HashTemp        !Hash put in Hash2Array - no need for pointer since always annihilating if storing hashes
                        ENDIF
                        VecSlot=VecSlot+1
!                        IF(TLocalAnnihilation) THEN
!                            PartsinExcitLevel(ExcitLevel)=PartsinExcitLevel(ExcitLevel)+1
!                        ENDIF
                    enddo

                    Acceptances=Acceptances+ABS(Child)      !Sum the number of created children to use in acceptance ratio
                
                ENDIF   !End if child created

!We now have to decide whether the parent particle (j) wants to self-destruct or not...
                iDie=AttemptDiePar(CurrentDets(:,j),CurrentH(j),CurrentIC(j))
!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births

                NoDied=NoDied+iDie          !Update death counter
                
                IF(iDie.le.0) THEN
!This indicates that the particle is spared and we may want to create more...copy them across to NewDets
!If iDie < 0, then we are creating the same particles multiple times. Copy accross (iDie+1) copies of particle
        
                    do l=1,abs(iDie)+1    !We need to copy accross one more, since we need to include the original spared particle
                        NewDets(:,VecSlot)=CurrentDets(:,j)
!                        NewDets(:,VecSlot)=CurrentDets(:,j)
                        NewSign(VecSlot)=CurrentSign(j)
!Copy excitation generator accross
!                        IF(ASSOCIATED(NewExcits(VecSlot)%PointToExcit)) THEN
!                            WRITE(6,*) "Problem is here",NewExcits(VecSlot)%IndexinExArr
!                        ENDIF
                        IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(j),NewExcits(VecSlot),.true.)
                        NewIC(VecSlot)=CurrentIC(j)
                        NewH(VecSlot)=CurrentH(j)
                        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) Hash2Array(VecSlot)=HashArray(j)
                        VecSlot=VecSlot+1
                    enddo

                ELSEIF(iDie.eq.1) THEN
!The particle wants to be killed. Nothing needs to be copied accross, but we do want to kill the old excitation generator.
                    IF(.not.TRegenExcitgens) CALL DissociateExitgen(CurrentExcits(j))

                ELSE!IF(iDie.gt.1) THEN
!This indicates that extra particles want to be killed.
!Anti-particles will need to be created on the same determinant.

                    do l=1,iDie-1
                        NewDets(:,VecSlot)=CurrentDets(:,j)
!                        NewDets(:,VecSlot)=CurrentDets(:,j)
                        IF(CurrentSign(j)) THEN
!Copy accross new anti-particles
                            NewSign(VecSlot)=.FALSE.
                        ELSE
                            NewSign(VecSlot)=.TRUE.
                        ENDIF
!Copy excitation generator accross
                        IF(l.eq.iDie-1) THEN
!Delete old excitation generator behind us (ie the final one is a move, rather than copy)
                            IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(j),NewExcits(VecSlot),.true.)
                        ELSE
                            IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(j),NewExcits(VecSlot),.false.)
                        ENDIF
                        NewIC(VecSlot)=CurrentIC(j)
                        NewH(VecSlot)=CurrentH(j)
                        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) Hash2Array(VecSlot)=HashArray(j)
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
        rat=(TotWalkersNew+0.D0)/(MaxWalkersExcit+0.D0)
        IF(rat.gt.0.9) THEN
            WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkersExcit"
            CALL FLUSH(6)
        ENDIF
        IF(TotWalkersNew.eq.0) THEN
            CALL Stop_All("PerformFciMCycPar","All walkers have died on a node")
        ENDIF
        
        CALL halt_timer(Walker_Time)
        CALL set_timer(Annihil_Time,30)
        
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

        ELSEIF(TAnnihilonProc) THEN
!In this routine, the particles are just annihilated on their own processors. This means that all simulations are independent and not influenced by each other.
!This means that there is no communication between processors and so should be much faster as the system size increases.

            CALL AnnihilateonProc(TotWalkersNew)
            Annihilated=Annihilated+(TotWalkersNew-TotWalkers)

        ELSE
!This routine now cancels down the particles with opposing sign on each determinant

            CALL AnnihilatePartPar(TotWalkersNew)
            Annihilated=Annihilated+(TotWalkersNew-TotWalkers)

        ENDIF
        
        CALL halt_timer(Annihil_Time)
        
!Find the total number of particles at HF (x sign) across all nodes. If this is negative, flip the sign of all particles.
        AllNoatHF=0
        
!Find sum of noathf, and then use an AllReduce to broadcast it to all nodes
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
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
!            IF(ierr.ne.MPI_SUCCESS) THEN
!                WRITE(6,*) "Problem in broadcasting new TSinglePartPhase"
!                CALL Stop_All("PerformFciMCycPar","Problem in broadcasting new TSinglePartPhase")
!            ENDIF
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
!                CALL FLUSH(6)
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


!A routine to annihilate particles separatly on each node. This should mean less annihilation occurs, but it is effect running nProcessors separate simulations.
!If there are enough particles, then this should be sufficient. Less memory is required, since no hashes need to be stored. Also, no communication is needed,
!so the routine should scale better as the number of walkers grows.
    SUBROUTINE AnnihilateonProc(TotWalkersNew)
        TYPE(ExcitPointer) :: TempExcit
        REAL*8 :: TempH
        INTEGER :: TempIC
        INTEGER :: TotWalkersNew,j,k,l,DetCurr(NEl),VecSlot,TotWalkersDet
        INTEGER :: DetLT

        TempExcit%PointToExcit=>null()
!First, it is necessary to sort the list of determinants
        CALL SortPartsPar(TotWalkersNew,NewDets(:,1:TotWalkersNew),NEl)

!Once ordered, each block of walkers on similar determinants can be analysed, and the residual walker concentration moved to CurrentDets
        j=1
!j is the counter over all uncancelled walkers - it indicates when we have reached the end of the list of total walkers
!DetCurr is the current determinant
        DetCurr(:)=NewDets(:,j)
        TempIC=NewIC(j)
        TempH=NewH(j)
        IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(j),TempExcit,.true.) !This will delete what is behind it - is that ok?
        
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
                    CurrentDets(:,VecSlot)=DetCurr(:)
                    CurrentSign(VecSlot)=.true.
                    CurrentIC(VecSlot)=TempIC
                    CurrentH(VecSlot)=TempH
                    IF(.not.TRegenExcitgens) CALL CopyExitgenPar(TempExcit,CurrentExcits(VecSlot),.false.)
                    VecSlot=VecSlot+1
                enddo
            ELSE
!Negative sign particles want to populate this determinant
                do l=1,abs(TotWalkersDet)
                    CurrentDets(:,VecSlot)=DetCurr(:)
                    CurrentSign(VecSlot)=.false.
                    CurrentIC(VecSlot)=TempIC
                    CurrentH(VecSlot)=TempH
                    IF(.not.TRegenExcitgens) CALL CopyExitgenPar(TempExcit,CurrentExcits(VecSlot),.false.)
                    VecSlot=VecSlot+1
                enddo
            ENDIF
!Now update the current determinant
            DetCurr(:)=NewDets(:,j)
            TempIC=NewIC(j)
            TempH=NewH(j)
            IF(.not.TRegenExcitgens) THEN
                CALL DissociateExitGen(TempExcit)
                CALL CopyExitGenPar(NewExcits(j),TempExcit,.true.)
            ENDIF
        enddo
!The new number of residual cancelled walkers is given by one less that VecSlot again.
        TotWalkers=VecSlot-1

        RETURN

    END SUBROUTINE AnnihilateonProc

!This routine sorts the particles before annihilation. It is identical to the routine in the serial version, but the data which is taken is different.
! Based on SORTI, SORTPARTS sorts arrays of integers, representing the determinant the walkers are on
! It then takes all the corresponding info with it
! Dets is the array (length N) of integers to sort
! NElecs is the length (in numbers of integers) of each element of Dets
! Vectors of NewXXX will be sorted correspondingly
    SUBROUTINE SortPartsPar(N,Dets,NElecs)
        TYPE(ExcitPointer) :: ExcitTemp
        REAL*8 :: HTemp
        INTEGER :: ICTemp
        INTEGER :: TempDet(NElecs)     !This stores a single element of the vector temporarily     
        LOGICAL :: WSignTemp
        INTEGER N,I,L,IR,J,NElecs
        INTEGER Dets(NElecs,N)
        INTEGER DETLT

        ExcitTemp%PointToExcit=>null()
        IF(N.LE.1) RETURN
        L=N/2+1 
        IR=N
10      CONTINUE
        IF(L.GT.1)THEN
            L=L-1
            TempDet(:)=Dets(:,L)
            HTemp=NewH(L)
            ICTemp=NewIC(L)
            WSignTemp=NewSign(L)
            IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(L),ExcitTemp,.true.) !This will delete what is behind it - is this ok?
        ELSE
            TempDet(:)=Dets(:,IR)      !Copy IRth elements to temp
            HTemp=NewH(IR)
            ICTemp=NewIC(IR)
            WSignTemp=NewSign(IR)
            IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(IR),ExcitTemp,.true.)

            Dets(:,IR)=Dets(:,1)    !Copy 1st element to IRth element
            NewH(IR)=NewH(1)
            NewIC(IR)=NewIC(1)
            NewSign(IR)=NewSign(1)
            IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(1),NewExcits(IR),.true.)
            IR=IR-1
            IF(IR.EQ.1)THEN
                Dets(:,1)=TempDet(:)    !Copy temp element to 1st element
                NewH(1)=HTemp
                NewIC(1)=ICTemp
                NewSign(1)=WSignTemp
                IF(.not.TRegenExcitgens) CALL CopyExitgenPar(ExcitTemp,NewExcits(1),.true.)
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
                Dets(:,I)=Dets(:,J)     !Copy Jth element to Ith element
                NewH(I)=NewH(J)
                NewIC(I)=NewIC(J)
                NewSign(I)=NewSign(J)
                IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(J),NewExcits(I),.true.)
                I=J
                J=J+J
            ELSE
                J=IR+1
            ENDIF
            GO TO 20
        ENDIF
        Dets(:,I)=TempDet(:)
        NewH(I)=HTemp
        NewIC(I)=ICTemp
        NewSign(I)=WSignTemp
        IF(.not.TRegenExcitgens) CALL CopyExitgenPar(ExcitTemp,NewExcits(I),.true.)
        GO TO 10

    END SUBROUTINE SortPartsPar


!A routine to annihilate particles in parallel. This involves separating hashes by abs(mod(hash,nProc)) to each node and annihilating there,       
!before sending back the annihilated particles to be removed from their original processors.
    SUBROUTINE AnnihilatePartPar(TotWalkersNew)
        INTEGER :: i,j,k,ToAnnihilateIndex,TotWalkersNew,ierr,error,sendcounts(nProcessors)
        INTEGER :: TotWalkersDet,InitialBlockIndex,FinalBlockIndex,ToAnnihilateOnProc,VecSlot
        INTEGER :: disps(nProcessors),recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER :: Minsendcounts,Maxsendcounts,DebugIter
        REAL*8 :: PopDensity(0:NEl)
        INTEGER , ALLOCATABLE :: TempExcitLevel(:)
        INTEGER(KIND=i2) :: HashCurr!MinBin,RangeofBins,NextBinBound
        CHARACTER(len=*), PARAMETER :: this_routine='AnnihilatePartPar'

!        DebugIter=0
!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "Printing out annihilation debug info for Iteration: ",Iter,DebugIter
!        ENDIF

!        IF(TLocalAnnihilation) THEN
!!We need to calculate the approximate population density of each excitation level
!!ApproxExcitDets contains the approximate number of determinants in each excitation level
!!PartsinExcitlevel is the number of particles in each excitation level for the current iteration
!!PopDensity is simply the approximate population density of particles in a given excitation level
!            do i=0,NEl
!                PopDensity(i)=REAL(PartsinExcitLevel(i),r2)/ApproxExcitDets(i)
!            enddo
!            PartsinExcitLevel(:)=0  !Rezero for the next iteration
!!Allocate memory to hold the excitation levels. This is needed since the amount of local annihilation will be a function of
!!PopDensity which the particle is at. This means it needs to be taken with the hash.
!            ALLOCATE(TempExcitLevel(TotWalkersNew),stat=ierr)
!            TempExcitLevel(1:TotWalkersNew)=NewIC(1:TotWalkersNew)
!        ENDIF

!First, allocate memory to hold the signs and the hashes while we annihilate
        ALLOCATE(TempSign(TotWalkersNew),stat=ierr)
!Comment out the memallocs later
!        CALL LogMemAlloc('TempSign',TotWalkersNew,4,this_routine,TempSignTag,ierr)
        ALLOCATE(TempHash(TotWalkersNew),stat=ierr)
!        CALL LogMemAlloc('TempHash',TotWalkersNew,8,this_routine,TempHashTag,ierr)
        
!Temporary arrays, storing the signs and Hashes need ot be kept, as both these arrays are going to be mixed
        TempSign(1:TotWalkersNew)=NewSign(1:TotWalkersNew)
        TempHash(1:TotWalkersNew)=Hash2Array(1:TotWalkersNew)
    
!Create the arrays for index and process
        do i=1,TotWalkersNew
            IndexTable(i)=i
            ProcessVec(i)=iProcIndex
        enddo

!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) TotWalkersNew
!            do i=1,TotWalkersNew
!                WRITE(6,*) i,Hash2Array(i),IndexTable(i),ProcessVec(i),NewSign(i)
!            enddo
!        ENDIF

!Next, order the hash array, taking the index, CPU and sign with it...
!Order the array by abs(mod(Hash,nProcessors)). This will result in a more load-balanced system
!        IF(TLocalAnnihilation) THEN
!If we are locally annihilating, then we need to take the excitation level of each walker with the hash
!            CALL SortMod4I1LLong(TotWalkersNew,Hash2Array(1:TotWalkersNew),IndexTable(1:TotWalkersNew),ProcessVec(1:TotWalkersNew),TempExcitLevel(1:TotWalkersNew),NewSign(1:TotWalkersNew),nProcessors)
!        ELSE
            CALL SortMod3I1LLong(TotWalkersNew,Hash2Array(1:TotWalkersNew),IndexTable(1:TotWalkersNew),ProcessVec(1:TotWalkersNew),NewSign(1:TotWalkersNew),nProcessors)
!        ENDIF
!        CALL Sort3I1LLong(TotWalkersNew,Hash2Array(1:TotWalkersNew),IndexTable(1:TotWalkersNew),ProcessVec(1:TotWalkersNew),NewSign(1:TotWalkersNew))
        
!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "***************"
!            WRITE(6,*) TotWalkersNew
!            do i=1,TotWalkersNew
!                WRITE(6,*) Hash2Array(i),abs(mod(Hash2Array(i),nProcessors)),IndexTable(i),ProcessVec(i),NewSign(i)
!            enddo
!        ENDIF
        
!        IF(nProcessors.eq.1) THEN
!            CALL STOP_All("AnnihilatePartPar","One processor annihilation not available yet...")
!        ENDIF
 
!Create the send counts and disps for the AlltoAllv. Work out equal ranges of bins for the hashes
!        Rangeofbins=HUGE(Rangeofbins)/(nProcessors/2)
!        MinBin=HUGE(MinBin)*-1
!        NextBinBound=MinBin+Rangeofbins

!Send counts is the size of each block of ordered dets which are going to each processor. This could be binary searched for extra speed
        j=1
        do i=0,nProcessors-1    !Search through all possible values of abs(mod(Hash,nProcessors))

            do while((abs(mod(Hash2Array(j),nProcessors)).eq.i).and.(j.le.TotWalkersNew))
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
            WRITE(6,*) "SENDCOUNTS is: ",sendcounts(:)
            WRITE(6,*) "TOTWALKERSNEW is: ",TotWalkersNew
            CALL FLUSH(6)
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
        
!The disps however do want to be cumulative - this is the array indexing the start of the data block
        disps(1)=0      !Starting element is always the first element
        do i=2,nProcessors
            disps(i)=disps(i-1)+sendcounts(i-1)
        enddo
        
!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "SENDCOUNTS: "
!            WRITE(6,*) sendcounts(:)
!            WRITE(6,*) "DISPS: "
!            WRITE(6,*) disps(:)
!            CALL FLUSH(6)
!        ENDIF

!We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
        recvcounts(1:nProcessors)=0
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

!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "RECVCOUNTS: "
!            WRITE(6,*) recvcounts(:)
!            WRITE(6,*) "RECVDISPS: "
!            WRITE(6,*) recvdisps(:),MaxIndex
!            CALL FLUSH(6)
!        ENDIF

!Insert a load-balance check here...maybe find the s.d. of the sendcounts array - maybe just check the range first.
        IF(TotWalkersNew.gt.200) THEN
            IF((Maxsendcounts-Minsendcounts).gt.(TotWalkersNew/3)) THEN
                WRITE(6,"(A,I12)") "**WARNING** Parallel annihilation not optimally balanced on this node, for iter = ",Iter
                WRITE(6,*) "Sendcounts is: ",sendcounts(:)
!                CALL FLUSH(6)
            ENDIF
        ENDIF
        
!Now send the chunks of hashes to the corresponding processors
        CALL MPI_AlltoAllv(Hash2Array(1:TotWalkersNew),sendcounts,disps,MPI_DOUBLE_PRECISION,HashArray(1:MaxIndex),recvcounts,recvdisps,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,error)        
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in dividing up hashes in annihilation"
!            CALL Stop_All("AnnihilatePartPar","Error in dividing up hashes in annihilation")
!        ENDIF

!The signs of the hashes, index and CPU also need to be taken with them.
        CALL MPI_AlltoAllv(NewSign(1:TotWalkersNew),sendcounts,disps,MPI_LOGICAL,CurrentSign,recvcounts,recvdisps,MPI_LOGICAL,MPI_COMM_WORLD,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in dividing up hashes in annihilation"
!            CALL Stop_All("AnnihilatePartPar","Error in dividing up hashes in annihilation")
!        ENDIF
        CALL MPI_AlltoAllv(IndexTable(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,Index2Table,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in dividing up hashes in annihilation"
!            CALL Stop_All("AnnihilatePartPar","Error in dividing up hashes in annihilation")
!        ENDIF
        CALL MPI_AlltoAllv(ProcessVec(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,Process2Vec,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in dividing up hashes in annihilation"
!            CALL Stop_All("AnnihilatePartPar","Error in dividing up hashes in annihilation")
!        ENDIF
!        IF(TLocalAnnihilation) THEN
!!If we are locally annihilating, then we need to take the excitation level of the particle with us.
!!We can send the information to CurrentIC - this is where the final information will be stored, but currently, it is redundant.
!            CALL MPI_AlltoAllv(TempExcitLevel(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,CurrentIC,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
!        ENDIF
        
!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "AFTER DIVISION:   - No. on processor is: ",MaxIndex
!            do i=1,MaxIndex
!                WRITE(6,*) HashArray(i),abs(mod(HashArray(i),nProcessors)),Index2Table(i),Process2Vec(i),CurrentSign(i)
!            enddo
!            CALL FLUSH(6)
!        ENDIF

!The hashes now need to be sorted again - this time by their number
!        IF(TLocalAnnihilation) THEN
!If we are locally annihilating, then we need to take the excitation level with us.
!            CALL Sort4I1LLong(MaxIndex,HashArray(1:MaxIndex),Index2Table(1:MaxIndex),Process2Vec(1:MaxIndex),CurrentIC(1:MaxIndex),CurrentSign(1:MaxIndex))
!        ELSE
            CALL Sort3I1LLong(MaxIndex,HashArray(1:MaxIndex),Index2Table(1:MaxIndex),Process2Vec(1:MaxIndex),CurrentSign(1:MaxIndex))
!        ENDIF
        
!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "AFTER DIVISION & ORDERING:   - No. on processor is: ",MaxIndex
!            do i=1,MaxIndex
!                WRITE(6,*) HashArray(i),abs(mod(HashArray(i),nProcessors)),Index2Table(i),Process2Vec(i),CurrentSign(i)
!            enddo
!            CALL FLUSH(6)
!        ENDIF

!Work out the index of the particles which want to be annihilated
        j=1
        ToAnnihilateIndex=1
        do while(j.le.MaxIndex)
            TotWalkersDet=0
            InitialBlockIndex=j
            FinalBlockIndex=j-1         !Start at j-1 since we are increasing FinalBlockIndex even with the first det in the next loop
            HashCurr=HashArray(j)
            do while((HashArray(j).eq.HashCurr).and.(j.le.MaxIndex))
!First loop counts walkers in the block - TotWalkersDet is then the residual sign of walkers on that determinant
                IF(CurrentSign(j)) THEN
                    TotWalkersDet=TotWalkersDet+1
                ELSE
                    TotWalkersDet=TotWalkersDet-1
                ENDIF
                FinalBlockIndex=FinalBlockIndex+1
                j=j+1
            enddo

!            IF(Iter.eq.DebugIter) THEN
!                WRITE(6,*) "Common block of dets found from ",InitialBlockIndex," ==> ",FinalBlockIndex
!                WRITE(6,*) "Sum of signs in block is: ",TotWalkersDet
!                CALL FLUSH(6)
!            ENDIF

!            IF(TLocalAnnihilation) THEN
!!This is an attempt to approximate the expected annihilation rates when the occupancy of a determinant is only 1.
!                IF(InitialBlockIndex.eq.FinalBlockIndex) THEN
!!The occupancy of the determinant is only one
!!The walker is at an excitation level of CurrentIC(InitialBlockIndex). The only parameter the local annihilation depends on is the population 
!!density of that excitation level, stored in PopDensity
!                    IF(AttemptLocalAnn(PopDensity(CurrentIC(InitialBlockIndex)))) THEN
!!Particle is killed, even though it is the lone occupier of the determinant
!                        TotWalkersDet=0     !By setting TotWalkersDet to zero, it will kill the particle in the next section
!                        LocalAnn=LocalAnn+1 !Keep a track of the number of particles locally annihilated
!!                        IF(HashCurr.eq.HFHash) THEN
!!                            WRITE(6,*) "HF Determinant particle locally annihilated"
!!                        ENDIF
!                    ENDIF
!                ENDIF
!            ENDIF
            
            do k=InitialBlockIndex,FinalBlockIndex
!Second run through the block of same determinants marks walkers for annihilation
                IF(TotWalkersDet.eq.0) THEN
!All walkers in block want to be annihilated from now on.
                    IndexTable(ToAnnihilateIndex)=Index2Table(k)
                    ProcessVec(ToAnnihilateIndex)=Process2Vec(k)
!                    Hash2Array(ToAnnihilateIndex)=HashArray(k)     !This is not strictly needed - remove after checking
!                    NewSign(ToAnnihilateIndex)=CurrentSign(k)       !This is also need needed, but useful for checking
!                    IF(Iter.eq.DebugIter) WRITE(6,*) "Annihilating from if block 1",j,k
                    ToAnnihilateIndex=ToAnnihilateIndex+1
!                    IF(HashCurr.eq.HFHash) THEN
!                        WRITE(6,*) "HF Determinant particle annihilated"
!                    ENDIF
                ELSEIF((TotWalkersDet.lt.0).and.(CurrentSign(k))) THEN
!Annihilate if block has a net negative walker count, and current walker is positive
                    IndexTable(ToAnnihilateIndex)=Index2Table(k)
                    ProcessVec(ToAnnihilateIndex)=Process2Vec(k)
!                    Hash2Array(ToAnnihilateIndex)=HashArray(k)     !This is not strictly needed - remove after checking
!                    NewSign(ToAnnihilateIndex)=CurrentSign(k)       !This is also need needed, but useful for checking
!                    IF(Iter.eq.DebugIter) WRITE(6,*) "Annihilating from if block 2",j,k
                    ToAnnihilateIndex=ToAnnihilateIndex+1
!                    IF(HashCurr.eq.HFHash) THEN
!                        WRITE(6,*) "HF Determinant particle annihilated"
!                    ENDIF
                ELSEIF((TotWalkersDet.gt.0).and.(.not.CurrentSign(k))) THEN
!Annihilate if block has a net positive walker count, and current walker is negative
                    IndexTable(ToAnnihilateIndex)=Index2Table(k)
                    ProcessVec(ToAnnihilateIndex)=Process2Vec(k)
!                    Hash2Array(ToAnnihilateIndex)=HashArray(k)     !This is not strictly needed - remove after checking
!                    NewSign(ToAnnihilateIndex)=CurrentSign(k)       !This is also need needed, but useful for checking
!                    IF(Iter.eq.DebugIter) WRITE(6,*) "Annihilating from if block 3",j,k
                    ToAnnihilateIndex=ToAnnihilateIndex+1
!                    IF(HashCurr.eq.HFHash) THEN
!                        WRITE(6,*) "HF Determinant particle annihilated"
!                    ENDIF
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
!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "Number of particles to annihilate from hashes on this processor: ",ToAnnihilateIndex
!            CALL FLUSH(6)
!        ENDIF

!The annihilation is complete - particles to be annihilated are stored in IndexTable and need to be sent back to their original processor
!To know which processor that is, we need to order the particles to be annihilated in terms of their CPU, i.e. ProcessVec(1:ToAnnihilateIndex)
!Is the list already ordered according to CPU? Is this further sort even necessary?

        IF(ToAnnihilateIndex.gt.1) THEN
!Do not actually have to take indextable, hash2array or newsign with it...
            CALL Sort2IILongL(ToAnnihilateIndex,ProcessVec(1:ToAnnihilateIndex),IndexTable(1:ToAnnihilateIndex),Hash2Array(1:ToAnnihilateIndex),NewSign(1:ToAnnihilateIndex))
        ENDIF

!We now need to regenerate sendcounts and disps
        sendcounts(1:nProcessors)=0
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
        recvcounts(1:nProcessors)=0
!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        CALL MPI_AlltoAll(sendcounts,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,error)

!We can now get recvdisps from recvcounts in the same way we obtained disps from sendcounts
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo

        ToAnnihilateonProc=recvdisps(nProcessors)+recvcounts(nProcessors)
        
!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "FOR RETURN OF ANNIHILATED PARTICLES, SENDCOUNTS: ",sendcounts(:)
!            WRITE(6,*) "DISPS: ",disps(:)
!            WRITE(6,*) "RECVCOUNTS: ",recvcounts(:)
!            WRITE(6,*) "RECVDISPS: ",recvdisps(:)
!            WRITE(6,*) "ToAnnihilateOnProc: ",ToAnnihilateonProc
!            CALL FLUSH(6)
!        ENDIF

!Perform another matrix transpose of the annihilation data using MPI_AlltoAllv, to send the data back to its correct Processor
!The signs of the hashes, index and CPU also need to be taken with them. (CPU does not need to be taken - every element of CPU should be equal to the rank of the processor+1)
!Hash also does not need to be taken, but will be taken as a precaution
!        CALL MPI_AlltoAllv(Hash2Array(1:TotWalkersNew),sendcounts,disps,MPI_DOUBLE_PRECISION,HashArray,recvcounts,recvdisps,mpilongintegertype,MPI_COMM_WORLD,error)        
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in sending back annihilated particles"
!            CALL Stop_All("AnnihilatePartPar","Error in sending back annihilated particles")
!        ENDIF
!        CALL MPI_AlltoAllv(NewSign(1:TotWalkersNew),sendcounts,disps,MPI_LOGICAL,CurrentSign,recvcounts,recvdisps,MPI_LOGICAL,MPI_COMM_WORLD,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in sending back annihilated particles"
!            CALL Stop_All("AnnihilatePartPar","Error in sending back annihilated particles")
!        ENDIF
!        CALL MPI_AlltoAllv(IndexTable(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,Index2Table,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        CALL MPI_AlltoAllv(IndexTable(1:ToAnnihilateonProc),sendcounts,disps,MPI_INTEGER,Index2Table,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in sending back annihilated particles"
!            CALL Stop_All("AnnihilatePartPar","Error in sending back annihilated particles")
!        ENDIF
!        CALL MPI_AlltoAllv(ProcessVec(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,Process2Vec,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in sending back annihilated particles"
!            CALL Stop_All("AnnihilatePartPar","Error in sending back annihilated particles")
!        ENDIF

!TEST
!        do i=1,ToAnnihilateonProc
!            IF(Process2Vec(i).ne.(iProcIndex)) THEN
!                CALL Stop_All("AnnihilatePartPar","AlltoAllv performed incorrectly")
!            ENDIF
!        enddo

!Index2Table now is a list, of length "ToAnnihilateonProc", of walkers which should NOT be transferred to the next array. 
!Order the list according to this index (Hash and sign does not need to be sorted, but will for debugging purposes)

        CALL SORTIILongL(ToAnnihilateonProc,Index2Table(1:ToAnnihilateonProc),HashArray(1:ToAnnihilateonProc),CurrentSign(1:ToAnnihilateonProc))

!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "Number of hashes originally on processor which need to be removed=",ToAnnihilateonProc
!            WRITE(6,*) "To annihilate from processor: "
!            do i=1,ToAnnihilateonProc
!                WRITE(6,*) Index2Table(i),HashArray(i),CurrentSign(i)
!            enddo
!        ENDIF

!TEST - do the hashes and signs match the ones that are returned?
!        do i=1,ToAnnihilateonProc
!            IF(TempHash(Index2Table(i)).ne.(HashArray(i))) THEN
!                CALL Stop_All("AnnihilatePartPar","Incorrect Hash returned")
!            ENDIF
!            IF(TempSign(Index2Table(i))) THEN
!                IF(.not.CurrentSign(i)) THEN
!                    CALL Stop_All("AnnihilatePartPar","Incorrect Sign returned")
!                ENDIF
!            ELSE
!                IF(CurrentSign(i)) THEN
!                    CALL Stop_All("AnnihilatePartPar","Incorrect Sign returned")
!                ENDIF
!            ENDIF
!        enddo
        

        IF(ToAnnihilateonProc.ne.0) THEN
!Copy across the data, apart from ones which have an index given by the indicies in Index2Table(1:ToAnnihilateonProc)
            VecSlot=1       !VecSlot is the index in the final array of TotWalkers
            i=1             !i is the index in the original array of TotWalkersNew
            do j=1,ToAnnihilateonProc
!Loop over all particles to be annihilated
!                IF(Iter.eq.DebugIter) WRITE(6,*) Index2Table(j)
                do while(i.lt.Index2Table(j))
!Copy accross all particles less than this number
                    CurrentDets(:,VecSlot)=NewDets(:,i)
                    CurrentIC(VecSlot)=NewIC(i)
                    CurrentH(VecSlot)=NewH(i)
                    IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(i),CurrentExcits(VecSlot),.true.)
                    HashArray(VecSlot)=TempHash(i)
                    CurrentSign(VecSlot)=TempSign(i)
                    i=i+1
                    VecSlot=VecSlot+1
                enddo
                IF(.not.TRegenExcitgens) CALL DissociateExitgen(NewExcits(i))    !Destroy particles if not copying accross
                i=i+1
            enddo

!Now need to copy accross the residual - from Index2Table(ToAnnihilateonProc) to TotWalkersNew
            do i=Index2Table(ToAnnihilateonProc)+1,TotWalkersNew
                CurrentDets(:,VecSlot)=NewDets(:,i)
                CurrentIC(VecSlot)=NewIC(i)
                CurrentH(VecSlot)=NewH(i)
                IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(i),CurrentExcits(VecSlot),.true.)
                HashArray(VecSlot)=TempHash(i)
                CurrentSign(VecSlot)=TempSign(i)
                VecSlot=VecSlot+1
            enddo

        ELSE
!No particles annihilated
            VecSlot=1
            do i=1,TotWalkersNew
                CurrentDets(:,VecSlot)=NewDets(:,i)
                CurrentIC(VecSlot)=NewIC(i)
                CurrentH(VecSlot)=NewH(i)
                IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(i),CurrentExcits(VecSlot),.true.)
                HashArray(VecSlot)=TempHash(i)
                CurrentSign(VecSlot)=TempSign(i)
                VecSlot=VecSlot+1
            enddo
        ENDIF
                
        TotWalkers=VecSlot-1

!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "FINAL CONFIGURATION: "
!            do i=1,TotWalkers
!                WRITE(6,*) i,HashArray(i),CurrentSign(i)
!            enddo
!        ENDIF

        IF((TotWalkersNew-TotWalkers).ne.ToAnnihilateonProc) THEN
            WRITE(6,*) TotWalkers,TotWalkersNew,ToAnnihilateonProc,Iter
            CALL FLUSH(6)
            CALL Stop_All("AnnihilatePartPar","Problem with numbers when annihilating")
        ENDIF

        DEALLOCATE(TempSign)
!        CALL LogMemDealloc(this_routine,TempSignTag)
        DEALLOCATE(TempHash)
!        CALL LogMemDealloc(this_routine,TempHashTag)
!        IF(TLocalAnnihilation) THEN
!            DEALLOCATE(TempExcitLevel)
!        ENDIF
        
        RETURN

    END SUBROUTINE AnnihilatePartPar


!This routine sums in the energy contribution from a given walker and updates stats such as mean excit level
    SUBROUTINE SumEContrib(DetCurr,ExcitLevel,WSign)
        INTEGER :: DetCurr(NEl),ExcitLevel,i,HighIndex,LowIndex
        LOGICAL :: WSign,CompiPath
        TYPE(HElement) :: HOffDiag

        MeanExcitLevel=MeanExcitLevel+real(ExcitLevel,r2)
        IF(MinExcitLevel.gt.ExcitLevel) MinExcitLevel=ExcitLevel
        IF(MaxExcitLevel.lt.ExcitLevel) MaxExcitLevel=ExcitLevel
        IF(ExcitLevel.eq.0) THEN
            IF(WSign) THEN
                IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF+1
                NoatHF=NoatHF+1
                AvSign=AvSign+1.D0
                AvSignHFD=AvSignHFD+1.D0
            ELSE
                IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF-1
                NoatHF=NoatHF-1
                AvSign=AvSign-1.D0
                AvSignHFD=AvSignHFD-1.D0
            ENDIF
        ELSEIF(ExcitLevel.eq.2) THEN
            NoatDoubs=NoatDoubs+1
!At double excit - find and sum in energy
            HOffDiag=GetHElement2(HFDet,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
            IF(WSign) THEN
                IF(Iter.gt.NEquilSteps) SumENum=SumENum+REAL(HOffDiag%v,r2)
                AvSign=AvSign+1.D0
                AvSignHFD=AvSignHFD+1.D0
            ELSE
                IF(Iter.gt.NEquilSteps) SumENum=SumENum-REAL(HOffDiag%v,r2)
                AvSign=AvSign-1.D0
                AvSignHFD=AvSignHFD-1.D0
            ENDIF
        ELSE
            IF(WSign) THEN
                AvSign=AvSign+1.D0
            ELSE
                AvSign=AvSign-1.D0
            ENDIF
        ENDIF

!        IF(TLocalAnnihilation) THEN
!We need to count the population of walkers at each excitation level for each iteration
!            PartsinExcitLevel(ExcitLevel)=PartsinExcitLevel(ExcitLevel)+1
!        ENDIF

        IF(TAutoCorr) THEN
!First element is HF Det to histogram. Then come doubles, triples and quads

            IF(ExcitLevel.eq.4) THEN
                LowIndex=2+NoACDets(2)+NoACDets(3)
                HighIndex=NoAutoDets
            ELSEIF(ExcitLevel.eq.3) THEN
                LowIndex=2+NoACDets(2)
                HighIndex=1+NoACDets(2)+NoACDets(3)
            ELSEIF(ExcitLevel.eq.2) THEN
                LowIndex=2
                HighIndex=1+NoACDets(2)
            ELSEIF(ExcitLevel.eq.0) THEN
                LowIndex=1
                HighIndex=1
            ELSE
                LowIndex=0
            ENDIF

            IF(LowIndex.ne.0) THEN

                do i=LowIndex,HighIndex
            
                    IF(CompiPath(DetCurr,AutoCorrDets(:,i),NEl)) THEN
!The walker is at a determinant for which we want to calculate the autocorrelation function
                        IF(WSign) THEN
                            IF(TFlippedSign) THEN
                                WeightatDets(i)=WeightatDets(i)-1
                            ELSE
                                WeightatDets(i)=WeightatDets(i)+1
                            ENDIF
                        ELSE
                            IF(TFlippedSign) THEN
                                WeightatDets(i)=WeightatDets(i)+1
                            ELSE
                                WeightatDets(i)=WeightatDets(i)-1
                            ENDIF
                        ENDIF
                        EXIT
                    ENDIF

                enddo
            ENDIF
        ENDIF
        
        RETURN

    END SUBROUTINE SumEContrib
    
!This is a routine to create a graph from the walker at nI, and ascribe new walkers at each vertex on the graph, according to a number of applications of the rho matrix.
    SUBROUTINE ResumGraphPar(nI,WSign,VecSlot,VecInd)
        INTEGER :: nI(NEl),VecSlot,VecInd
        LOGICAL :: WSign
        REAL*8 :: Prob

!This routine will create the graph. It will calculate it differently depending on the size of the graph, and whether excitation generators are stored. 
        CALL CreateGraphPar(nI,VecInd,Prob)

!Apply the rho matrix successive times. This could be improved if large numbers of applications of rho are needed by diagonalising the rho matrix
        CALL ApplyRhoMatPar()
        
        CALL CreateNewPartsPar(nI,VecInd,WSign,CurrentH(VecInd),VecSlot,Prob)   !Create particles proportionally to the magnitude of the vector elements in GraphVec

        RETURN

    END SUBROUTINE ResumGraphPar

!This routine will create a graph from a given initial determinant
    SUBROUTINE CreateGraphPar(nI,VecInd,Prob)
        TYPE(ExcitPointer) :: TempExcitgen
        INTEGER :: nI(NEl),nJ(NEl),VecInd,i,j,IC,iCount,Attempts
        REAL*8 :: Prob,ExcitProb
        LOGICAL :: SameDet,CompiPath
        TYPE(HElement) :: Hij,Hjj

        TempExcitgen%PointToExcit=>null()

!First, set up the excitation generators for the root determinant
        IF(.not.TRegenExcitgens) THEN
            CALL SetupExitgenPar(nI,CurrentExcits(VecInd))
        ELSE
            CALL SetupExitgenPar(nI,TempExcitgen)
        ENDIF

        IF(NDets.eq.2) THEN
!We know that determinants are not going to be regenerated if NDets=2, so we can do this in a slightly simpler way
            
            IF(TRegenExcitgens) THEN
                CALL GenRandSymExcitIt3(nI,TempExcitgen%PointToExcit,nJ,Seed,IC,0,Prob,iCount)
            ELSE
                CALL GenRandSymExcitIt3(nI,CurrentExcits(VecInd)%PointToExcit,nJ,Seed,IC,0,Prob,iCount)
            ENDIF

            Hij=GetHElement2(nI,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
            Hjj=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)

            GraphRhoMat(1,1)=1.D0-Tau*(CurrentH(VecInd)-DiagSft)
            GraphRhoMat(1,2)=-Tau*real(Hij%v,r2)
            GraphRhoMat(2,1)=GraphRhoMat(1,2)
            GraphRhoMat(2,2)=1.D0-Tau*((real(Hjj%v,r2)-Hii)-DiagSft)
                
            DetsinGraph(:,2)=nJ(:)
            GraphKii(2)=REAL(Hjj%v,r2)-Hii              !store det generated and kii element

        ELSE
!Zero the matrix
            GraphRhoMat(1:NDets,1:NDets)=0.D0
            
            GraphRhoMat(1,1)=1.D0-Tau*(CurrentH(VecInd)-DiagSft)    !This is the first rho-matrix element

            i=2
            Attempts=0
            do while(i.lt.NDets)    !Loop until all determinants found

                IF(TRegenExcitgens) THEN
                    CALL GenRandSymExcitIt3(nI,TempExcitgen%PointToExcit,nJ,Seed,IC,0,Prob,iCount)
                ELSE
                    CALL GenRandSymExcitIt3(nI,CurrentExcits(VecInd)%PointToExcit,nJ,Seed,IC,0,Prob,iCount)
                ENDIF

                SameDet=.false.
                do j=2,(i-1)
                    IF(CompiPath(nJ,DetsinGraph(:,j),NEl)) THEN
!Determinants are the same as already created determinant - ignore it

                        SameDet=.true.
                        Attempts=Attempts+1
                        IF(Attempts.gt.1000) THEN
                            WRITE(6,*) "iCOUNT IS: ", iCount
                            CALL FLUSH(6)
                            CALL Stop_All("CreateGraphPar","More than 1000 attempts needed to grow graph")
                        ENDIF
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
                    DetsinGraph(1:NEl,i)=nJ(1:NEl)
!First find connection to root
                    Hij=GetHElement2(nI,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
                    GraphRhoMat(1,i)=-Tau*REAL(Hij%v,r2)
                    GraphRhoMat(i,1)=GraphRhoMat(1,i)

!Then find connection to other determinants
                    do j=2,(i-1)
                        Hij=GetHElement2(nJ,DetsInGraph(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,-1,ECore)
                        GraphRhoMat(i,j)=-Tau*REAL(Hij%v,r2)
                        GraphRhoMat(j,i)=GraphRhoMat(i,j)
                    enddo

!Find diagonal element - and store it for later on...
                    Hjj=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                    GraphKii(i)=REAL(Hjj%v,r2)-Hii                !Again, the root value is not stored
                    GraphRhoMat(i,i)=1.D0-Tau*(GraphKii(i)-DiagSft)

                    i=i+1

                ENDIF
            enddo

        ENDIF   !End if ndets=2

        IF(TRegenExcitgens) THEN
!Deallocate excitation generator if we are regenerating them
            CALL DissociateExitgen(TempExcitgen)
        ENDIF

        RETURN

    END SUBROUTINE CreateGraphPar

    
!This routine creates new particles from the vector which results from the 
    SUBROUTINE CreateNewPartsPar(nI,VecInd,WSign,Kii,VecSlot,Prob)
        IMPLICIT NONE
        LOGICAL :: WSign,TempSign
        INTEGER :: nI(NEl),VecInd
        INTEGER :: i,j,VecSlot,Create,ExcitLevel,iGetExcitLevel_2
        INTEGER(KIND=i2) :: HashTemp
        REAL*8 :: Ran2,rat,Kii,Kjj,Prob

!First deal with the excitations...
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
                IF(ExcitLevel.eq.0) THEN
                    IF(ABS(GraphKii(i)).gt.1.D-07) THEN
                        CALL Stop_All("ResumGraphPar","Diagonal K-mat element should be zero for HF particles")
                    ENDIF
                ENDIF
                IF(Create.lt.0) THEN
                    TempSign=.false.
                ELSE
                    TempSign=.true.
                ENDIF
                
                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    HashTemp=CreateHash(DetsInGraph(:,i))
                ENDIF
                
!Now actually create the particles in NewDets and NewSign
                do j=1,abs(Create)
                    NewDets(1:NEl,VecSlot)=DetsInGraph(1:NEl,i)
                    NewSign(VecSlot)=TempSign
                    NewIC(VecSlot)=ExcitLevel
                    NewH(VecSlot)=GraphKii(i)       !Diagonal H El previously stored
                    IF(.not.TRegenExcitgens) NewExcits(VecSlot)%PointToExcit=>null()
                    IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                        Hash2Array(VecSlot)=HashTemp
                    ENDIF
                    VecSlot=VecSlot+1
                enddo

            ENDIF

        enddo

!Now deal with root
        Create=INT(abs(GraphVec(1)))

        rat=abs(GraphVec(1))-REAL(Create,r2)    !rat is now the fractional part, to be assigned stochastically
        IF(rat.gt.Ran2(Seed)) Create=Create+1
        
        IF((abs(Create)).gt.0) THEN
        
            IF(.not.WSign) Create=-Create
            IF(GraphVec(1).lt.0.D0) Create=-Create

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

!Now actually create the particles in NewDets and NewSign. It is the same particle as the parent particle.
            do j=1,abs(Create)
                NewDets(:,VecSlot)=nI(:)
                NewSign(VecSlot)=TempSign
!Copy excitation generator accross
                IF(j.eq.abs(Create)) THEN
                    IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(VecInd),NewExcits(VecSlot),.true.)
                ELSE
                    IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(VecInd),NewExcits(VecSlot),.false.)
                ENDIF
                NewIC(VecSlot)=CurrentIC(VecInd)
                NewH(VecSlot)=CurrentH(VecInd)
                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    Hash2Array(VecSlot)=HashArray(VecInd)
                ENDIF
                VecSlot=VecSlot+1
            enddo

        ENDIF

        RETURN

    END SUBROUTINE CreateNewPartsPar

!This applies the rho matrix successive times to a root determinant. From this, GraphVec is filled with the correct probabilities for the determinants in the graph
    SUBROUTINE ApplyRhoMatPar()
        REAL*8 :: TempVec(NDets)
        INTEGER :: i,j,k

!        IF(NDets.eq.2) THEN
!
!            GraphVec(1)=1.D0
!            GraphVec(2)=0.D0
!
!            do i=1,RhoApp
!            
!                TempVec(1)=(GraphRhoMat(1,1)*GraphVec(1))+(GraphRhoMat(1,2)*GraphVec(2))
!                TempVec(2)=(GraphRhoMat(2,1)*GraphVec(1))+(GraphRhoMat(2,2)*GraphVec(2))
!            
!                GraphVec(1)=TempVec(1)
!                GraphVec(2)=TempVec(2)
!
!            enddo
!
!        ELSE
            
            GraphVec(1:NDets)=0.d0
            GraphVec(1)=1.D0        !Set the initial vector to be 1 at the root (i.e. for one walker initially)
        
            do i=1,RhoApp

                CALL DGEMV('n',NDets,NDets,1.D0,GraphRhoMat,NDets,GraphVec,1,0.D0,TempVec,1)
                GraphVec(1:NDets)=TempVec(1:NDets)
                TempVec(1:NDets)=0.d0
            
!            do j=1,NDets
!                TempVec(j)=0.D0
!                do k=1,NDets
!                    TempVec(j)=TempVec(j)+GraphRhoMat(j,k)*GraphVec(k)
!                enddo
!            enddo
!            GraphVec(:)=TempVec(:)

            enddo

!        ENDIF

        RETURN
    END SUBROUTINE ApplyRhoMatPar


!Every StepsSft steps, update the diagonal shift value (the running value for the correlation energy)
!We don't want to do this too often, since we want the population levels to acclimatise between changing the shifts
    SUBROUTINE CalcNewShift()
        INTEGER :: error,rc,MaxWalkersProc,MaxAllowedWalkers
        INTEGER :: inpair(7),outpair(7)
        REAL*8 :: TempSumNoatHF,MeanWalkers

!This first call will calculate the GrowRate for each processor, taking culling into account
        CALL UpdateDiagSftPar()

!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

!We need to collate the information from the different processors
!Inpair and outpair are used to package variables to save on latency time
        inpair(1)=TotWalkers
        inpair(2)=Annihilated
        inpair(3)=NoatDoubs
        inpair(4)=NoBorn
        inpair(5)=NoDied
        inpair(6)=SumWalkersCyc
        inpair(7)=LocalAnn
        outpair(:)=0
!        CALL MPI_Reduce(TotWalkers,AllTotWalkers,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!Find total Annihilated,Total at HF and Total at doubles
!        CALL MPI_Reduce(Annihilated,AllAnnihilated,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!!        CALL MPI_Reduce(NoatHF,AllNoatHF,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)  !This is done every iteration now
!        CALL MPI_Reduce(NoatDoubs,AllNoatDoubs,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(NoBorn,AllNoBorn,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(NoDied,AllNoDied,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(SumWalkersCyc,AllSumWalkersCyc,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(inpair,outpair,7,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
        AllTotWalkers=outpair(1)
        AllAnnihilated=outpair(2)
        AllNoatDoubs=outpair(3)
        AllNoBorn=outpair(4)
        AllNoDied=outpair(5)
        AllSumWalkersCyc=outpair(6)
        AllLocalAnn=outpair(7)


        MeanWalkers=REAL(AllTotWalkers,r2)/REAL(nProcessors,r2)
        MaxAllowedWalkers=NINT((MeanWalkers/12.D0)+MeanWalkers)

!Find the range of walkers on different nodes to see if we need to even up the distribution over nodes
!        inpair(1)=TotWalkers
!        inpair(2)=iProcIndex
        CALL MPI_Reduce(TotWalkers,MaxWalkersProc,1,MPI_INTEGER,MPI_MAX,root,MPI_COMM_WORLD,error)
!        MaxWalkersProc=outpair(1)
!        WRITE(6,*) "***",MaxWalkersProc,MaxAllowedWalkers,MeanWalkers
!        CALL MPI_Reduce(inpair,outpair,1,MPI_2INTEGER,MPI_MINLOC,root,MPI_COMM_WORLD,error)
!        MinWalkersProc=outpair(1)

        IF(iProcIndex.eq.root) THEN
!            RangeWalkers=MaxWalkersProc-MinWalkersProc
!            IF(RangeWalkers.gt.300) THEN
            IF((MaxWalkersProc.gt.MaxAllowedWalkers).and.(AllTotWalkers.gt.(nProcessors*500))) THEN
                TBalanceNodes=.true.
            ENDIF
        ENDIF
!We need to tell all nodes whether to balance the nodes or not...
        CALL MPI_BCast(TBalanceNodes,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
        
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

        AvSign=AvSign/real(SumWalkersCyc,r2)
        AvSignHFD=AvSignHFD/real(SumWalkersCyc,r2)
        CALL MPI_Reduce(AvSign,AllAvSign,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(AvSignHFD,AllAvSignHFD,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        IF(iProcIndex.eq.Root) THEN
            AllAvSign=AllAvSign/real(nProcessors,r2)
            AllAvSignHFD=AllAvSignHFD/real(nProcessors,r2)
        ENDIF

!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        TempSumNoatHF=real(SumNoatHF,r2)
        CALL MPIDSumRoot(TempSumNoatHF,1,AllSumNoatHF,Root)
        CALL MPIDSumRoot(SumENum,1,AllSumENum,Root)

!To find minimum and maximum excitation levels, search for them using MPI_Reduce
!        inpair(1)=MaxExcitLevel
!        inpair(2)=iProcIndex

        CALL MPI_Reduce(MaxExcitLevel,AllMaxExcitLevel,1,MPI_INTEGER,MPI_MAX,Root,MPI_COMM_WORLD,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in finding max excitation level"
!            CALL MPI_ABORT(MPI_COMM_WORLD,rc,error)
!        ENDIF
!Max Excit Level is found on processor outpair(2) and is outpair(1)
!        IF(iProcIndex.eq.Root) THEN
!            AllMaxExcitLevel=outpair(1)
!        ENDIF

!        inpair(1)=MinExcitLevel
!        inpair(2)=iProcIndex
        CALL MPI_Reduce(MinExcitLevel,AllMinExcitLevel,1,MPI_INTEGER,MPI_MIN,Root,MPI_COMM_WORLD,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in finding min excitation level"
!            CALL MPI_ABORT(MPI_COMM_WORLD,rc,error)
!        ENDIF
!        IF(iProcIndex.eq.Root) THEN
!            AllMinExcitLevel=outpair(1)
!        ENDIF

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

        CALL WriteFCIMCStats()

!Now need to reinitialise all variables on all processers
        MinExcitLevel=NEl+10
        MaxExcitLevel=0
        MeanExcitLevel=0.D0
        SumWalkersCyc=0
        AvSign=0.D0        !Rezero this quantity - <s> is now a average over the update cycle
        AvSignHFD=0.D0     !This is the average sign over the HF and doubles
        Annihilated=0
        LocalAnn=0
        Acceptances=0
        NoBorn=0
        NoDied=0
!Reset TotWalkersOld so that it is the number of walkers now
        TotWalkersOld=TotWalkers

!Also reinitialise the global variables - should not necessarily need to do this...
        AllSumENum=0.D0
        AllSumNoatHF=0.D0
        AllTotWalkersOld=AllTotWalkers
        AllTotWalkers=0
        AllGrowRate=0.D0
        AllMeanExcitLevel=0.D0
        AllAvSign=0.D0
        AllAvSignHFD=0.D0
        AllSumWalkersCyc=0
        AllAnnihilated=0
        AllLocalAnn=0
        AllNoatHF=0
        AllNoatDoubs=0
        AllNoBorn=0
        AllNoDied=0

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
                IF(.not.TRegenExcitgens) THEN
                    CALL DissociateExitgen(CurrentExcits(Chosen))    !First, destroy the excitation generator of the chosen particle
                    CALL CopyExitgenPar(CurrentExcits(TotWalkers),CurrentExcits(Chosen),.true.) !Then move the end particle to its place.
                ENDIF

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
                IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(i),CurrentExcits(VecSlot),.false.)

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
            CullInfo(1:10,1:3)=0
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
            CullInfo(1:10,1:3)=0

        ENDIF
        
!        DiagSft=DiagSft-(log(GrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
!        IF((DiagSft).gt.0.D0) THEN
!            WRITE(6,*) "***WARNING*** - DiagSft trying to become positive..."
!            STOP
!        ENDIF

    END SUBROUTINE UpdateDiagSftPar


!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalcPar()
        use CalcData, only : EXCITFUNCS
        use Determinants , only : GetH0Element3
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet
        INTEGER :: DetLT,VecSlot,error,HFConn,MemoryAlloc,iMaxExcit,nStore(6),nJ(Nel)
        TYPE(HElement) :: rh,TempHii
        REAL*8 :: TotDets,SymFactor,Choose
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMCPar'
        CHARACTER(len=12) :: abstr

!        CALL MPIInit(.false.)       !Initialises MPI - now have variables iProcIndex and nProcessors
        
!Set timed routine names
        Walker_Time%timer_name='WalkerTime'
        Annihil_Time%timer_name='AnnihilTime'
        ACF_Time%timer_name='ACFTime'

        IF(TDebug) THEN
!This will open a file called LOCALPOPS-"iprocindex" on unit number 11 on every node.
            abstr=''
            write(abstr,'(I2)') iProcIndex
            abstr='LOCALPOPS-'//adjustl(abstr)
            OPEN(11,FILE=abstr,STATUS='UNKNOWN')
        ENDIF

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

!Setup excitation generator for the HF determinant. If we are using assumed sized excitgens, this will also be assumed size.
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,HFExcit%nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
        ALLOCATE(HFExcit%ExcitData(HFExcit%nExcitMemLen),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("InitFCIMC","Problem allocating excitation generator")
        HFExcit%ExcitData(1)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,HFExcit%ExcitData,nJ,iMaxExcit,0,nStore,3)
        HFExcit%nPointed=0
!Indicate that the excitation generator is now correctly allocated.

!        CALL SetupExitgenPar(HFDet,HFExcit)
        CALL GetSymExcitCount(HFExcit%ExcitData,HFConn)

!Initialise random number seed - since the seeds need to be different on different processors, subract processor rank from random number
        Seed=G_VMC_Seed-iProcIndex

!Calculate Hii
        TempHii=GetHElement2(HFDet,HFDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
        Hii=REAL(TempHii%v,r2)
        TempHii=GetH0Element3(HFDet)
        Fii=REAL(TempHii%v,r2)

        TBalanceNodes=.false.   !Assume that the nodes are initially load-balanced

!Initialise variables for calculation on each node
        ProjectionE=0.D0
        AvSign=0.D0
        AvSignHFD=0.D0
        SumENum=0.D0
        SumNoatHF=0
        NoatHF=0
        NoatDoubs=0
        MeanExcitLevel=0.D0
        MinExcitLevel=NEl+10
        MaxExcitLevel=0
        LocalAnn=0
        Annihilated=0
        Acceptances=0
        PreviousCycles=0
        NoBorn=0
        NoDied=0

!Also reinitialise the global variables - should not necessarily need to do this...
        AllSumENum=0.D0
        AllNoatHF=0
        AllNoatDoubs=0
        AllSumNoatHF=0.D0
        AllGrowRate=0.D0
        AllMeanExcitLevel=0.D0
        AllSumWalkersCyc=0
        AllAvSign=0.D0
        AllAvSignHFD=0.D0
        AllNoBorn=0
        AllNoDied=0
        AllLocalAnn=0
        AllAnnihilated=0

!Need to declare a new MPI type to deal with the long integers we use in the hashing, and when reading in from POPSFILEs
!        CALL MPI_Type_create_f90_integer(18,mpilongintegertype,error)
!        CALL MPI_Type_commit(mpilongintegertype,error)
        
        IF(TPopsFile.and.(mod(iWritePopsEvery,StepsSft).ne.0)) THEN
            CALL Warning("InitFCIMCCalc","POPSFILE writeout should be a multiple of the update cycle length.")
        ENDIF

        IF(TNoAnnihil) THEN
            WRITE(6,*) "No Annihilation to occur. Results are likely not to converge on right value. Proceed with caution. "
        ELSEIF(TAnnihilonproc) THEN
            WRITE(6,*) "Annihilation will occur on each processors' walkers only. This should be faster, but result in less annihilation."
            WRITE(6,*) "This is equivalent to running seperate calculations."
        ENDIF
        IF(TReadPops) THEN
!List of things that readpops can't work with...
            IF(TStartSinglePart.or.TStartMP1) THEN
                CALL Stop_All("InitFCIMCCalcPar","ReadPOPS cannot work with StartSinglePart or StartMP1")
            ENDIF
        ENDIF

        IF(TResumFciMC) THEN
            IF(NDets.ge.2) THEN
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
            ELSEIF(TFixShiftShell) THEN
                CALL Stop_All("InitFCIMCCalcPar","Fixing the shift of the certain excitation levels cannot be used within ResumFCIMC")
            ENDIF
            IF(iProcIndex.eq.root) THEN
                WRITE(6,*) "Resumming in multiple transitions to/from each excitation"
                WRITE(6,"(A,I5,A)") "Graphs to resum will consist of ",NDets," determinants."
            ENDIF
!            IF(.not.TNoAnnihil) THEN
!                CALL Stop_All("InitFCIMCCalcPar","Annihilation is not currently compatable with ResumFCIMC in parallel")
!            ENDIF
        ENDIF
        WRITE(6,*) ""
        WRITE(6,*) "Performing Parallel FCIMC...."
        IF(TStartSinglePart) THEN
            WRITE(6,"(A,F9.3,A,I9)") "Initial number of particles set to 1, and shift will be held at ",DiagSft," until particle number on root node gets to ",InitWalkers
        ELSE
            WRITE(6,*) "Initial number of walkers per processor chosen to be: ", InitWalkers
        ENDIF
        WRITE(6,*) "Maximum connectivity of HF determinant is: ",HFConn
        WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp
        IF(.not.TReadPops) THEN
            WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft
        ENDIF
        IF(TStartSinglePart) THEN
            TSinglePartPhase=.true.
            IF(TReadPops) THEN
                CALL Stop_All("InitFciMCCalcPar","Cannot read in POPSFILE as well as starting with a single particle")
            ENDIF
            IF(TStartMP1) THEN
                CALL Stop_All("InitFciMCCalcPar","Cannot start with a single particle, and as the MP1 wavefunction")
            ENDIF
        ELSE
            TSinglePartPhase=.false.
        ENDIF
        IF(TFixShiftShell) THEN
            IF(iProcIndex.eq.root) THEN
                WRITE(6,"(A,I5,A,F20.10)") "All excitations up to ",ShellFix," will have their shift fixed at ",FixShift
                WRITE(6,*) "With this option, results are going to be non-exact, but can be used to equilibrate calculations"
            ENDIF
        ENDIF
        IF(ICILevel.ne.0) THEN
!We are truncating the excitations at a certain value
            TTruncSpace=.true.
            WRITE(6,'(A,I4)') "Truncating the S.D. space at determinants will an excitation level w.r.t. HF of: ",ICILevel
            IF(TResumFciMC) CALL Stop_All("InitFciMCPar","Space cannot be truncated with ResumFCIMC")
        ENDIF

        IF(TAutoCorr) THEN
!We want to calculate the autocorrelation function over the determinants
!            IF(iLagMin.lt.0) THEN
!                CALL Stop_All("InitFciMCPar","LagMin cannot be less than zero (and when equal 0 should be strictly 1")
!            ELSEIF(iLagMax.gt.NMCyc) THEN
!                CALL Stop_All("InitFciMCPar","LagMax cannot be greater than the number of cycles.")
!            ELSEIF(iLagStep.lt.1) THEN
!                CALL Stop_All("InitFciMCPar","LagStep cannot be less than 1")
!            ENDIF

            CALL ChooseACFDets()

!            WRITE(6,*) "Storing information to calculate the ACF at end of simulation..."
        ENDIF

        IF(TLocalAnnihilation) THEN
!If we are locally annihilating, then we need to know the walker density for a given excitation level, for which we need to approximate number of determinants
!in each excitation level
            ALLOCATE(ApproxExcitDets(0:NEl))
            ALLOCATE(PartsinExcitLevel(0:NEl))
            TotDets=1.D0
            do i=0,NEl
                ApproxExcitDets(i)=Choose(NEl,i)*Choose(nBasis-NEl,i)
!                WRITE(6,*) "Excitation level: ",i,ApproxExcitDets(i)
            enddo
            SymFactor=ApproxExcitDets(2)/(HFConn+0.D0)
            do i=1,NEl
                ApproxExcitDets(i)=ApproxExcitDets(i)/SymFactor
                WRITE(6,*) "Excitation level: ",i,ApproxExcitDets(i)
                TotDets=TotDets+ApproxExcitDets(i)
            enddo
            WRITE(6,*) "Approximate size of determinant space is: ",TotDets
            PartsinExcitLevel(:)=0  !Zero the array to hold the population of walkers in each excitation level
        ELSE
            SymFactor=(Choose(NEl,2)*Choose(nBasis-NEl,2))/(HFConn+0.D0)
            TotDets=1.D0
            do i=1,NEl
                WRITE(6,*) "Approximate excitation level population: ",i,NINT((Choose(NEl,i)*Choose(nBasis-NEl,i))/SymFactor)
                TotDets=TotDets+(Choose(NEl,i)*Choose(nBasis-NEl,i))/SymFactor
            enddo
            WRITE(6,*) "Approximate size of determinant space is: ",NINT(TotDets)
        ENDIF


        IF(TStartMP1) THEN
!Start the initial distribution off at the distribution of the MP1 eigenvector

            WRITE(6,"(A)") "Starting run with particles populating double excitations proportionally to MP1 wavevector..."
            CALL InitWalkersMP1Par()

        ELSEIF(TReadPops) THEN
!Read in particles from multiple POPSFILES for each processor
            WRITE(6,*) "Reading in initial particle configuration from POPSFILES..."

            CALL ReadFromPopsFilePar()

        ELSE
!initialise the particle positions - start at HF with positive sign

!Set the maximum number of walkers allowed
            MaxWalkers=NINT(MemoryFac*InitWalkers)
            MaxWalkersExcit=NINT(MemoryFacExcit*InitWalkers)
            WRITE(6,"(A,F14.5)") "Memory Factor for walkers is: ",MemoryFac
            WRITE(6,"(A,F14.5)") "Memory Factor for excitation generators is: ",MemoryFacExcit
            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node of: ",MaxWalkers
            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for excitgens (+ particle paths) of: ",MaxWalkersExcit

!Put a barrier here so all processes synchronise
            CALL MPI_Barrier(MPI_COMM_WORLD,error)
!Allocate memory to hold walkers
            ALLOCATE(WalkVecDets(NEl,MaxWalkersExcit),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkersExcit*NEl,4,this_routine,WalkVecDetsTag,ierr)
            WalkVecDets(1:NEl,1:MaxWalkersExcit)=0
            ALLOCATE(WalkVec2Dets(NEl,MaxWalkersExcit),stat=ierr)
            CALL LogMemAlloc('WalkVec2Dets',MaxWalkersExcit*NEl,4,this_routine,WalkVec2DetsTag,ierr)
            WalkVec2Dets(1:NEl,1:MaxWalkersExcit)=0

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
            WalkVecH(:)=0.d0
            ALLOCATE(WalkVec2H(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2H',MaxWalkers,8,this_routine,WalkVec2HTag,ierr)
            WalkVec2H(:)=0.d0
            
!            MemoryAlloc=((2*NEl+8)*MaxWalkers)*4    !Memory Allocated in bytes
            MemoryAlloc=((8*MaxWalkers)+(2*NEl*MaxWalkersExcit))*4    !Memory Allocated in bytes

            IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                ALLOCATE(HashArray(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('HashArray',MaxWalkers,8,this_routine,HashArrayTag,ierr)
                HashArray(:)=0
                ALLOCATE(Hash2Array(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('Hash2Array',MaxWalkers,8,this_routine,Hash2ArrayTag,ierr)
                Hash2Array(:)=0
                ALLOCATE(IndexTable(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('IndexTable',MaxWalkers,4,this_routine,IndexTableTag,ierr)
                IndexTable(1:MaxWalkers)=0
                ALLOCATE(Index2Table(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('Index2Table',MaxWalkers,4,this_routine,Index2TableTag,ierr)
                Index2Table(1:MaxWalkers)=0
                ALLOCATE(ProcessVec(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('ProcessVec',MaxWalkers,4,this_routine,ProcessVecTag,ierr)
                ProcessVec(1:MaxWalkers)=0
                ALLOCATE(Process2Vec(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('Process2Vec',MaxWalkers,4,this_routine,Process2VecTag,ierr)
                Process2Vec(1:MaxWalkers)=0

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
                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    HashArray(1)=HFHash
                ENDIF
            ELSE

                do j=1,InitWalkers
                    CurrentDets(:,j)=HFDet(:)
                    CurrentSign(j)=.true.
                    CurrentIC(j)=0
                    CurrentH(j)=0.D0
                    IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                        HashArray(j)=HFHash
                    ENDIF
                enddo
            ENDIF

            WRITE(6,"(A,F14.6,A)") "Initial memory (without excitgens + temp arrays) consists of : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
            WRITE(6,*) "Initial memory allocation sucessful..."
            CALL FLUSH(6)

!Put a barrier here so all processes synchronise
            CALL MPI_Barrier(MPI_COMM_WORLD,error)
            IF(.not.TRegenExcitgens) THEN

                ALLOCATE(WalkVecExcits(MaxWalkersExcit),stat=ierr)
                ALLOCATE(WalkVec2Excits(MaxWalkersExcit),stat=ierr)
                ALLOCATE(ExcitGens(MaxWalkersExcit),stat=ierr)
                ALLOCATE(FreeIndArray(MaxWalkersExcit),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("InitFCIMMCCalcPar","Error in allocating walker excitation generators")

                do j=1,MaxWalkersExcit
                    NULLIFY(WalkVecExcits(j)%PointToExcit)
                    NULLIFY(WalkVec2Excits(j)%PointToExcit)
                    ExcitGens(j)%nPointed=0
                    FreeIndArray(j)=j       !All points are initially free
                enddo

!Allocate pointers to the correct excitation arrays
                CurrentExcits=>WalkVecExcits
                NewExcits=>WalkVec2Excits

!Initialise the first Free Excitgens indices...initially whole list is free
                BackOfList=1    !I.e. first allocation should be at ExcitGens(FreeIndArray(1))
                FrontOfList=1   !i.e. first free index should be put at FreeIndArray(1)

!Setup the first particle on HF...
                CALL SetupExitGenPar(HFDet,CurrentExcits(1))
                
                IF(.not.TStartSinglePart) THEN
                    do j=2,InitWalkers
!Copy the HF excitation generator accross to each initial particle
                        CALL CopyExitGenPar(CurrentExcits(1),CurrentExcits(j),.false.)
                    enddo
                ENDIF
                MemoryAlloc=((HFExcit%nExcitMemLen)+4)*4*MaxWalkersExcit

                WRITE(6,"(A,I12)") "Size of HF excitgen is: ",HFExcit%nExcitMemLen
                WRITE(6,"(A,F14.6,A)") "Probable maximum memory for excitgens is : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
                WRITE(6,*) "Initial allocation of excitation generators successful..."
            ELSE
                WRITE(6,*) "Excitation generators will not be stored, but regenerated each time they are needed..."
            ENDIF
            WRITE(6,"(A,F14.6,A)") "Temp Arrays for annihilation cannot be more than : ",REAL(MaxWalkers*12,r2)/1048576.D0," Mb/Processor"
            CALL FLUSH(6)
        
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

        ENDIF   !End if initial walkers method

        CullInfo(1:10,1:3)=0
        NoCulls=0
        
        IF(iProcIndex.eq.root) THEN
!Print out initial starting configurations
            WRITE(6,*) ""
            IF(TLocalAnnihilation) THEN
                WRITE(6,"(A)") "       Step     Shift      WalkerCng    GrowRate       TotWalkers    LocalAnn   TotAnnihil    NoDied    NoBorn    Proj.E          NoatHF NoatDoubs      AvSign    AvSignHF+D   AccRat       MeanEx     MinEx MaxEx"
                WRITE(15,"(A)") "#       Step     Shift      WalkerCng    GrowRate       TotWalkers   LocalAnn    TotAnnihil    NoDied    NoBorn    Proj.E          NoatHF NoatDoubs       AvSign    AvSignHF+D   AccRat       MeanEx     MinEx MaxEx"

            ELSE
                WRITE(6,"(A)") "       Step     Shift      WalkerCng    GrowRate       TotWalkers    Annihil    NoDied    NoBorn    Proj.E          NoatHF NoatDoubs      AvSign    AvSignHF+D   AccRat       MeanEx     MinEx MaxEx"
                WRITE(15,"(A)") "#       Step     Shift      WalkerCng    GrowRate       TotWalkers    Annihil    NoDied    NoBorn    Proj.E          NoatHF NoatDoubs       AvSign    AvSignHF+D   AccRat       MeanEx     MinEx MaxEx"
            
            ENDIF
        ENDIF
        
!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
        RETURN

    END SUBROUTINE InitFCIMCCalcPar

    SUBROUTINE ChooseACFDets()
        use SystemData , only : tAssumeSizeExcitgen
        INTEGER :: ierr,HFConn,nStore(6),i,nJ(NEl),iExcit,j,MaxIndex,ExcitLevel,iGetExcitLevel_2
        INTEGER :: iMaxExcit,ExcitLength
        REAL*8 , ALLOCATABLE :: TempMax(:)
        LOGICAL :: TurnBackAssumeExGen
        TYPE(HElement) :: Fjj,Hij,Fkk
        REAL*8 :: Compt,MaxWeight
        INTEGER , ALLOCATABLE :: ExcitGenTemp(:)
        CHARACTER(len=*), PARAMETER :: this_routine='ChooseACFDets'

!Commented out code is to just perform the ACF for the HF determinant.
!            NoAutoDets=1    !Initially, we are just after the ACF for one determinant
!Fill the autocorrdets array with the determinants that you want to work out the autocorrelation for.
!            AutoCorrDets(:,1)=HFDet(:)

!The code calculates the ACF for the number of doubles here, or all if it is more
!        NoAutoDets=10
!        WRITE(6,"(A,I5,A)") "Choosing the ",NoAutoDets," highest MP1 weight determinants to calculate the ACF for."
        NoAutoDets=1+NoACDets(2)+NoACDets(3)+NoACDets(4)    !Total number of dets to histogram is sum of the determinants for the individual excitation levels.
        WRITE(6,"(A,I5,A)") "Choosing the ",NoACDets(2)+1," highest MP1 weight determinants to calculate the ACF for."
        WRITE(6,"(A,I5,A,I5,A)") "Also picking ",NoACDets(3), " high weighted triply excited and ",NoACDets(4), " quadruply excited determinants."
        WRITE(6,*) "First determinant will be the HF determinant"
        CALL GetSymExcitCount(HFExcit%ExcitData,HFConn)
        IF(NoAutoDets.gt.HFConn) NoAutoDets=HFConn
        ALLOCATE(AutoCorrDets(NEl,NoAutoDets),stat=ierr)
        CALL LogMemAlloc('AutoCorrDets',NEl*NoAutoDets,4,this_routine,AutoCorrDetsTag,ierr)
        AutoCorrDets(:,:)=0
        
!Set the first determinant to find the ACF of to be the HF determinant
        AutoCorrDets(:,1)=HFDet(:)

!We do not know if tAssumeSizeExcitgen is on - if it is, then we can't enumerate all determinants. Get around this by simply regenerating it anyway.
!First, we need to turn off AssumeSizeExcitgen if it is on.
        IF(tAssumeSizeExcitgen) THEN
            TurnBackAssumeExGen=.true.
            tAssumeSizeExcitgen=.false.
        ELSE
            TurnBackAssumeExGen=.false.
        ENDIF

        ALLOCATE(TempMax(2:NoACDets(2)+1),stat=ierr)   !This will temporarily hold the largest components (+1 since we are no longer considering HF in NoACDets(2))
        IF(ierr.ne.0) THEN
            CALL Stop_All("ChooseACFDets","Problem allocating memory")
        ENDIF
        TempMax(:)=0.D0

!Setup excit generators for HF Determinant
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitLength,nJ,iMaxExcit,0,nStore,2)
        ALLOCATE(ExcitGenTemp(ExcitLength),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("ChooseACFDets","Problem allocating excitation generator")
        ExcitGenTemp(1)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGenTemp,nJ,iMaxExcit,0,nStore,2)

        do while(.true.)
!Generate double excitations
            CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.false.,ExcitGenTemp,nJ,iExcit,0,nStore,2)
            IF(nJ(1).eq.0) EXIT
            Hij=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
            CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fjj)
            Compt=real(Hij%v,r2)/(Fii-(REAL(Fjj%v,r2)))
            do j=2,NoACDets(2)+1!NoAutoDets
                IF(abs(Compt).gt.TempMax(j)) THEN
                    TempMax(j)=abs(Compt)
                    AutoCorrDets(:,j)=nJ(:)
                    EXIT
                ENDIF
            enddo
            
        enddo
        DEALLOCATE(ExcitGenTemp)
!Find largest weight MP1 contribution
        MaxWeight=0.D0
        do i=2,NoACDets(2)+1
            IF(TempMax(i).gt.MaxWeight) THEN
                MaxIndex=i
            ENDIF
        enddo
        DEALLOCATE(TempMax)
        
!        IF(iProcIndex.eq.root) THEN
!            WRITE(6,*) "*** Histogramming the following determinants:"
!            do i=1,NoAutoDets
!                do j=1,NEl
!                    WRITE(6,"(I4)",advance='no') AutoCorrDets(j,i)
!                enddo
!                WRITE(6,*) ""
!            enddo
!        ENDIF

!We have now found the largest weight Doubles. To guess at large weighted triples & quads, we simply take the largest weighted double, and find 
!its large MP1 weighted contributions with it as the reference.
!        WRITE(6,*) "MaxIndex = ", MaxIndex
!Setup excit generators for the highest weighted double excitation
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(AutoCorrDets(:,MaxIndex),NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitLength,nJ,iMaxExcit,0,nStore,3)
        ALLOCATE(ExcitGenTemp(ExcitLength),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("ChooseACFDets","Problem allocating excitation generator")
        ExcitGenTemp(1)=0
        CALL GenSymExcitIt2(AutoCorrDets(:,MaxIndex),NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGenTemp,nJ,iMaxExcit,0,nStore,3)

        CALL GetH0Element(AutoCorrDets(:,MaxIndex),NEl,Arr,nBasis,ECore,Fjj)
        ALLOCATE(TempMax(1:(NoACDets(3)+NoACDets(4))),stat=ierr)
        TempMax(:)=0.D0

        do while(.true.)
            CALL GenSymExcitIt2(AutoCorrDets(:,MaxIndex),NEl,G1,nBasis,nBasisMax,.false.,ExcitGenTemp,nJ,iExcit,0,nStore,3)
            IF(nJ(1).eq.0) EXIT
            Hij=GetHElement2(AutoCorrDets(:,MaxIndex),nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
            CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fkk)
            Compt=real(Hij%v,r2)/(real(Fjj%v,r2)-(real(Fkk%v,r2)))
            ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,4)
            IF(ExcitLevel.eq.3) THEN
!We have generated a triple - try to add it to the list
                do j=1,NoACDets(3)
                    IF(abs(Compt).gt.TempMax(j)) THEN
                        TempMax(j)=abs(Compt)
                        AutoCorrDets(:,(j+1+NoACDets(2)))=nJ(:)
                        EXIT
                    ENDIF
                enddo
            ELSEIF(ExcitLevel.eq.4) THEN
!We have generated a quad - try to add it to the list
                do j=NoACDets(3)+1,(NoACDets(3)+NoACDets(4))
                    IF(abs(Compt).gt.TempMax(j)) THEN
                        TempMax(j)=abs(Compt)
                        AutoCorrDets(:,(j+1+NoACDets(2)))=nJ(:)
                        EXIT
                    ENDIF
                enddo
            ENDIF
        enddo

!Deallocate TempExcitgen
        DEALLOCATE(TempMax)
!        CALL DissociateExitgen(ExcitGenTemp)
        DEALLOCATE(ExcitGenTemp)
        
!Find if any zeros
        do i=1,NoAutoDets
            IF(AutoCorrDets(1,i).eq.0) THEN
                NoAutoDets=NoAutoDets-1
                WRITE(6,*) "Could not find AutoCorrFunc determinant to histogram slot ",i
!Move these zeros to the end of the list
                do j=i,NoAutoDets
                    IF(AutoCorrDets(1,j).ne.0) THEN
                        AutoCorrDets(:,i)=AutoCorrDets(:,j)
                        EXIT
                    ENDIF
                enddo
            ENDIF
        enddo
            
!The number of occurunces of walkers at the determiants selected needs to be stored for all iterations
        ALLOCATE(WeightatDets(NoAutoDets),stat=ierr)
        CALL LogMemAlloc('WeightatDets',NoAutoDets,4,this_routine,WeightatDetsTag,ierr)
        WeightatDets(:)=0
!If we want to calculate the ACF, then we now do this in a seperate step. Now we simply write out the determinant populations at each iteration
!to a file called HFDoublePops
        
        IF(iProcIndex.eq.root) THEN
            WRITE(6,*) "Histogramming the following determinants:"
            do i=1,NoAutoDets
                do j=1,NEl
                    WRITE(6,"(I4)",advance='no') AutoCorrDets(j,i)
                enddo
                WRITE(6,*) ""
            enddo

            OPEN(44,FILE='HFDoublePops',STATUS='UNKNOWN')
        ENDIF

        IF(TurnBackAssumeExGen) THEN
!We turned off assumed sized excitation generators for this routine - turn it back on.
            tAssumeSizeExcitgen=.true.
        ENDIF

        RETURN

    END SUBROUTINE ChooseACFDets

!This routine will write out to a popsfile. It transfers all walkers to the head node sequentially, so does not want to be called too often
!When arriving at this routine, CurrentXXX are arrays with the data, and NewXXX will be used by the root processor to temporarily store the information
    SUBROUTINE WriteToPopsfilePar()
        REAL*8 :: TempSumNoatHF
        INTEGER :: error,WalkersonNodes(0:nProcessors-1)
        INTEGER :: Stat(MPI_STATUS_SIZE),Tag,Total,i,j,k

        CALL MPI_Barrier(MPI_COMM_WORLD,error)  !sync

!First, make sure we have up-to-date information - again collect AllTotWalkers,AllSumNoatHF and AllSumENum...
        CALL MPI_Reduce(TotWalkers,AllTotWalkers,1,MPI_INTEGER,MPI_Sum,root,MPI_COMM_WORLD,error)    
!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        TempSumNoatHF=real(SumNoatHF,r2)
        CALL MPIDSumRoot(TempSumNoatHF,1,AllSumNoatHF,Root)
        CALL MPIDSumRoot(SumENum,1,AllSumENum,Root)

!We also need to tell the root processor how many particles to expect from each node - these are gathered into WalkersonNodes
        CALL MPI_Gather(TotWalkers,1,MPI_INTEGER,WalkersonNodes,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)

        Tag=125

        IF(iProcIndex.eq.root) THEN
!First, check that we are going to receive the correct number of particles...
            Total=0
            do i=0,nProcessors-1
                Total=Total+WalkersonNodes(i)
            enddo
            IF(Total.ne.AllTotWalkers) THEN
                CALL Stop_All("WriteToPopsfilePar","Not all walkers accounted for...")
            ENDIF

!Write header information
            WRITE(6,*) "Writing to POPSFILE..."
            OPEN(17,FILE='POPSFILE',Status='unknown')
            WRITE(17,*) AllTotWalkers,"   TOTWALKERS (all nodes)"
            WRITE(17,*) DiagSft,"   DIAG SHIFT"
            WRITE(17,*) NINT(AllSumNoatHF,i2),"   SUMNOATHF (all nodes)"
            WRITE(17,*) AllSumENum,"   SUMENUM ( \sum<D0|H|Psi> - all nodes)"
            WRITE(17,*) Iter+PreviousCycles,"   PREVIOUS CYCLES"
            do j=1,TotWalkers
!First write out walkers on head node
                do k=1,NEl
                    WRITE(17,"(I5)",advance='no') CurrentDets(k,j)
                enddo
                WRITE(17,*) CurrentSign(j)
            enddo

!Now we need to receive the data from each other processor sequentially
            do i=1,nProcessors-1
!Run through all other processors...receive the data...
                CALL MPI_Recv(NewDets(:,1:WalkersonNodes(i)),WalkersonNodes(i)*NEl,MPI_INTEGER,i,Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(NewSign(1:WalkersonNodes(i)),WalkersonNodes(i),MPI_LOGICAL,i,Tag,MPI_COMM_WORLD,Stat,error)
                
!Then write it out...
                do j=1,WalkersonNodes(i)
                    do k=1,NEl
                        WRITE(17,"(I5)",advance='no') NewDets(k,j)
                    enddo
                    WRITE(17,*) NewSign(j)
                enddo

            enddo

            CLOSE(17)

        ELSE
!All other processors need to send their data to root...
            CALL MPI_Send(CurrentDets(:,1:TotWalkers),TotWalkers*NEl,MPI_INTEGER,root,Tag,MPI_COMM_WORLD,error)
            CALL MPI_Send(CurrentSign(1:TotWalkers),TotWalkers,MPI_LOGICAL,root,Tag,MPI_COMM_WORLD,error)
        ENDIF

!Reset the values of the global variables
        AllSumNoatHF=0.D0
        AllSumENum=0.D0
        AllTotWalkers=0

        RETURN

    END SUBROUTINE WriteToPopsfilePar


!This routine reads in particle configurations from a POPSFILE.
    SUBROUTINE ReadFromPopsfilePar()
        LOGICAL :: exists,First
        INTEGER :: AvWalkers,WalkerstoReceive(nProcessors)
        INTEGER*8 :: NodeSumNoatHF(nProcessors),TempAllSumNoatHF
        INTEGER :: TempInitWalkers,error,i,j,k,l,total,ierr,iGetExcitLevel_2,MemoryAlloc,Tag
        INTEGER :: Stat(MPI_STATUS_SIZE),AvSumNoatHF,VecSlot,IntegerPart,HFPointer
        REAL*8 :: Ran2,FracPart
        TYPE(HElement) :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        
        PreviousCycles=0    !Zero previous cycles
        SumENum=0.D0
        SumNoatHF=0
        DiagSft=0.D0
        Tag=124             !Set Tag

        IF(iProcIndex.eq.root) THEN
            INQUIRE(FILE='POPSFILE',EXIST=exists)
            IF(exists) THEN

                OPEN(17,FILE='POPSFILE',Status='old')
!Read in initial data on processors which have a popsfile
                READ(17,*) AllTotWalkers
                READ(17,*) DiagSft
                READ(17,*) TempAllSumNoatHF     !AllSumNoatHF stored as integer for compatability with serial POPSFILEs
                READ(17,*) AllSumENum
                READ(17,*) PreviousCycles

                AllSumNoatHF=REAL(TempAllSumNoatHF,r2)
                WRITE(6,*) "Number of cycles in previous simulation: ",PreviousCycles
                IF(NEquilSteps.gt.0) THEN
                    WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
                    NEquilSteps=0
                ENDIF
                IF(TZeroProjE) THEN
!Reset energy estimator
                    WRITE(6,*) "Resetting projected energy counters to zero..."
                    AllSumENum=0.D0
                    AllSumNoatHF=0.D0
                ENDIF

!Need to calculate the number of walkers each node will receive...
                AvWalkers=NINT(real(AllTotWalkers,r2)/real(nProcessors,r2))
!Divide up the walkers to receive for each node
                do i=1,nProcessors-1
                    WalkerstoReceive(i)=AvWalkers
                enddo
!The last processor takes the 'remainder'
                WalkerstoReceive(nProcessors)=AllTotWalkers-(AvWalkers*(nProcessors-1))

!Quick check to ensure we have all walkers accounted for
                total=0
                do i=1,nProcessors
                    total=total+WalkerstoReceive(i)
                enddo
                IF(total.ne.AllTotWalkers) THEN
                    CALL Stop_All("ReadFromPopsfilePar","All Walkers not accounted for when reading in from POPSFILE")
                ENDIF
                
!InitWalkers needs to be reset for the culling criteria
                InitWalkers=AvWalkers
                SumENum=AllSumENum/REAL(nProcessors,r2)     !Divide up the SumENum over all processors
                AvSumNoatHF=NINT(AllSumNoatHF/real(nProcessors,r2)) !This is the average Sumnoathf
                do i=1,nProcessors-1
                    NodeSumNoatHF(i)=INT(AvSumNoatHF,i2)
                enddo
                NodeSumNoatHF(nProcessors)=NINT(AllSumNoatHF,i2)-INT((AvSumNoatHF*(nProcessors-1)),i2)
                    
!Reset the global variables
                AllSumENum=0.D0
                AllSumNoatHF=0.D0

            ELSE
                CALL Stop_All("ReadFromPopsfilePar","No POPSFILE found")
            ENDIF

        ENDIF

        CALL MPI_Barrier(MPI_COMM_WORLD,error)  !Sync

!Now we need to scatter the WalkerstoReceive to each node, and allocate the desired memory to each node...
!Broadcast info which needs to go to all processors
        CALL MPI_BCast(DiagSft,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,error)
        CALL MPI_BCast(InitWalkers,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
        CALL MPI_BCast(SumENum,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,error)
!Scatter the number of walkers each node will receive to TempInitWalkers, and the SumNoatHF for each node which is distributed approximatly equally
        CALL MPI_Scatter(WalkerstoReceive,1,MPI_INTEGER,TempInitWalkers,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
        CALL MPI_Scatter(NodeSumNoatHF,1,MPI_DOUBLE_PRECISION,SumNoatHF,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,error)
        
!Now we want to allocate memory on all nodes.
        MaxWalkers=NINT(MemoryFac*InitWalkers)    !All nodes have the same amount of memory allocated
        MaxWalkersExcit=NINT(MemoryFacExcit*InitWalkers)
!Allocate memory to hold walkers at least temporarily
        ALLOCATE(WalkVecDets(NEl,MaxWalkersExcit),stat=ierr)
        CALL LogMemAlloc('WalkVecDets',MaxWalkersExcit*NEl,4,this_routine,WalkVecDetsTag,ierr)
        WalkVecDets(1:NEl,1:MaxWalkersExcit)=0
        ALLOCATE(WalkVecSign(MaxWalkersExcit),stat=ierr)
        CALL LogMemAlloc('WalkVecSign',MaxWalkersExcit*NEl,4,this_routine,WalkVecSignTag,ierr)

        IF(iProcIndex.eq.root) THEN
!Root process reads all walkers in and then sends them to the correct processor
            do i=nProcessors,1,-1
!Read in data for processor i
                do j=1,WalkerstoReceive(i)
                    READ(17,*) WalkVecDets(1:NEl,j),WalkVecSign(j)
                enddo

                IF(i.ne.1) THEN
!Now send data to processor i-1 (Processor rank goes from 0 -> nProcs-1). If i=1, then we want the data so stay at the root processor
                    CALL MPI_Send(WalkVecDets(:,1:WalkerstoReceive(i)),WalkerstoReceive(i)*NEl,MPI_INTEGER,i-1,Tag,MPI_COMM_WORLD,error)
                    CALL MPI_Send(WalkVecSign(1:WalkerstoReceive(i)),WalkerstoReceive(i),MPI_LOGICAL,i-1,Tag,MPI_COMM_WORLD,error)
                ENDIF

            enddo

            CLOSE(17)
        
        ENDIF

        do i=1,nProcessors-1
            IF(iProcIndex.eq.i) THEN
!All other processors want to pick up their data from root
                CALL MPI_Recv(WalkVecDets(:,1:TempInitWalkers),TempInitWalkers*NEl,MPI_INTEGER,0,Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(WalkVecSign(1:TempInitWalkers),TempInitWalkers,MPI_LOGICAL,0,Tag,MPI_COMM_WORLD,Stat,error)
            ENDIF
        enddo

        IF(iProcIndex.eq.root) WRITE(6,*) AllTotWalkers," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers  by a factor of: ",ScaleWalkers
            MaxWalkers=NINT(MemoryFac*(NINT(InitWalkers*ScaleWalkers)))   !InitWalkers here is simply the average number of walkers per node, not actual
            MaxWalkersExcit=NINT(MemoryFacExcit*(NINT(InitWalkers*ScaleWalkers)))   !InitWalkers here is simply the average number of walkers per node, not actual

!Allocate memory for walkvec2, which will temporarily hold walkers
            ALLOCATE(WalkVec2Dets(NEl,MaxWalkersExcit),stat=ierr)
            CALL LogMemAlloc('WalkVec2Dets',MaxWalkersExcit*NEl,4,this_routine,WalkVec2DetsTag,ierr)
            WalkVec2Dets(1:NEl,1:MaxWalkersExcit)=0
            ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)

!Scale up the integer part and fractional part seperately
            IntegerPart=INT(ScaleWalkers)       !Round to zero
            FracPart=ScaleWalkers-REAL(IntegerPart)

            VecSlot=1
            do l=1,TempInitWalkers
                do k=1,IntegerPart
                    WalkVec2Dets(1:NEl,VecSlot)=WalkVecDets(1:NEl,l)
                    WalkVec2Sign(VecSlot)=WalkVecSign(l)
                    VecSlot=VecSlot+1
                enddo
                IF(Ran2(Seed).lt.FracPart) THEN
!Stochastically choose whether to create another particle
                    WalkVec2Dets(1:NEl,VecSlot)=WalkVecDets(1:NEl,l)
                    WalkVec2Sign(VecSlot)=WalkVecSign(l)
                    VecSlot=VecSlot+1
                ENDIF
            enddo

!Redefine new (average) number of initial particles for culling criteria
            InitWalkers=NINT(InitWalkers*ScaleWalkers)
!Define TotWalkers as the number of walkers on that node.
            TotWalkers=VecSlot-1
            TotWalkersOld=TotWalkers
!Collate the information about the new total number of walkers
            CALL MPI_AllReduce(TotWalkers,AllTotWalkers,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,error)
            IF(iProcIndex.eq.root) THEN
                AllTotWalkersOld=AllTotWalkers
                WRITE(6,*) "Total number of initial walkers is now: ", AllTotWalkers
            ENDIF


!Deallocate the arrays used to hold the original particles, and reallocate with correct size.
            DEALLOCATE(WalkVecDets)
            CALL LogMemDealloc(this_routine,WalkVecDetsTag)
            DEALLOCATE(WalkVecSign)
            CALL LogMemDealloc(this_routine,WalkVecSignTag)
            ALLOCATE(WalkVecDets(NEl,MaxWalkersExcit),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkersExcit*NEl,4,this_routine,WalkVecDetsTag,ierr)
            WalkVecDets(1:NEl,1:MaxWalkersExcit)=0
            ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)

!Transfer scaled particles back accross to WalkVecDets
            do l=1,TotWalkers
                WalkVecDets(1:NEl,l)=WalkVec2Dets(1:NEl,l)
                WalkVecSign(l)=WalkVec2Sign(l)
            enddo

!Zero the second array for good measure
            WalkVec2Dets(1:NEl,1:MaxWalkers)=0

        ELSE
!We are not scaling the number of walkers...

            ALLOCATE(WalkVec2Dets(NEl,MaxWalkersExcit),stat=ierr)
            CALL LogMemAlloc('WalkVec2Dets',MaxWalkersExcit*NEl,4,this_routine,WalkVec2DetsTag,ierr)
            WalkVec2Dets(1:NEl,1:MaxWalkersExcit)=0
            ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)

            TotWalkers=TempInitWalkers      !Set the total number of walkers
            TotWalkersOld=TempInitWalkers
            IF(iProcIndex.eq.root) THEN
                AllTotWalkersOld=AllTotWalkers
                WRITE(6,*) "Total number of initial walkers is now: ",AllTotWalkers
            ENDIF

        ENDIF
            
        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft

!Need to now allocate other arrays
        ALLOCATE(WalkVecIC(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVecIC',MaxWalkers,4,this_routine,WalkVecICTag,ierr)
        WalkVecIC(1:MaxWalkers)=0
        ALLOCATE(WalkVec2IC(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec2IC',MaxWalkers,4,this_routine,WalkVec2ICTag,ierr)
        WalkVec2IC(1:MaxWalkers)=0
        ALLOCATE(WalkVecH(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVecH',MaxWalkers,8,this_routine,WalkVecHTag,ierr)
        WalkVecH=0.d0
        ALLOCATE(WalkVec2H(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec2H',MaxWalkers,8,this_routine,WalkVec2HTag,ierr)
        WalkVec2H=0.d0

        MemoryAlloc=(8+(2*NEl))*MaxWalkers*4    !Memory allocated in bytes

        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
            ALLOCATE(HashArray(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('HashArray',MaxWalkers,8,this_routine,HashArrayTag,ierr)
            HashArray(:)=0
            ALLOCATE(Hash2Array(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('Hash2Array',MaxWalkers,8,this_routine,Hash2ArrayTag,ierr)
            Hash2Array(:)=0
            ALLOCATE(IndexTable(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('IndexTable',MaxWalkers,4,this_routine,IndexTableTag,ierr)
            IndexTable(1:MaxWalkers)=0
            ALLOCATE(Index2Table(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('Index2Table',MaxWalkers,4,this_routine,Index2TableTag,ierr)
            Index2Table(1:MaxWalkers)=0
            ALLOCATE(ProcessVec(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('ProcessVec',MaxWalkers,4,this_routine,ProcessVecTag,ierr)
            ProcessVec(1:MaxWalkers)=0
            ALLOCATE(Process2Vec(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('Process2Vec',MaxWalkers,4,this_routine,Process2VecTag,ierr)
            Process2Vec(1:MaxWalkers)=0

            MemoryAlloc=MemoryAlloc+32*MaxWalkers
        ENDIF

        WRITE(6,"(A,F14.6,A)") "Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        CALL FLUSH(6)

!Allocate pointers to the correct walker arrays...
        CurrentDets=>WalkVecDets
        CurrentSign=>WalkVecSign
        CurrentIC=>WalkVecIC
        CurrentH=>WalkVecH
        NewDets=>WalkVec2Dets
        NewSign=>WalkVec2Sign
        NewIC=>WalkVec2IC
        NewH=>WalkVec2H

        IF(.not.TRegenExcitgens) THEN
            ALLOCATE(WalkVecExcits(MaxWalkersExcit),stat=ierr)
            ALLOCATE(WalkVec2Excits(MaxWalkersExcit),stat=ierr)
            ALLOCATE(ExcitGens(MaxWalkersExcit),stat=ierr)
            ALLOCATE(FreeIndArray(MaxWalkersExcit),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("ReadFromPopsfilePar","Error in allocating walker excitation generators")
            do j=1,MaxWalkersExcit
                NULLIFY(WalkVecExcits(j)%PointToExcit)
                NULLIFY(WalkVec2Excits(j)%PointToExcit)
                ExcitGens(j)%nPointed=0
                FreeIndArray(j)=j
            enddo

!Allocate pointers to the correct excitation arrays
            CurrentExcits=>WalkVecExcits
            NewExcits=>WalkVec2Excits

!Initialise the first Free Excitgens indices...initially whole list is free
            BackOfList=1    !I.e. first allocation should be at ExcitGens(FreeIndArray(1))
            FrontOfList=1   !i.e. first free index should be put at FreeIndArray(1)

            MemoryAlloc=((HFExcit%nExcitMemLen)+4)*4*MaxWalkers
            WRITE(6,"(A,F14.6,A)") "Probable maximum memory for excitgens is : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
            WRITE(6,*) "Initial allocation of excitation generators successful..."

        ELSE
            WRITE(6,*) "Excitgens will be regenerated when they are needed..."
        ENDIF

        WRITE(6,"(A,F14.6,A)") "Temp Arrays for annihilation cannot be more than : ",REAL(MaxWalkers*12,r2)/1048576.D0," Mb/Processor"
        CALL FLUSH(6)

!Now find out the data needed for the particles which have been read in...
        First=.true.
        do j=1,TotWalkers
            CurrentIC(j)=iGetExcitLevel_2(HFDet,CurrentDets(:,j),NEl,NEl)
            IF(CurrentIC(j).eq.0) THEN
                CurrentH(j)=0.D0
                IF(First) THEN
!First run - create the excitation.
                    IF(.not.TRegenExcitgens) CALL SetupExitgenPar(HFDet,CurrentExcits(j))
                    HFPointer=j
                    First=.false.
                ELSE
                    IF(.not.TRegenExcitgens) CALL CopyExitGenPar(CurrentExcits(HFPointer),CurrentExcits(j),.false.)
                ENDIF
                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    HashArray(j)=HFHash
                ENDIF
            ELSE
                HElemTemp=GetHElement2(CurrentDets(:,j),CurrentDets(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                CurrentH(j)=REAL(HElemTemp%v,r2)-Hii
                IF(.not.TRegenExcitgens) CurrentExcits(j)%PointToExcit=>null()
                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    HashArray(j)=CreateHash(CurrentDets(:,j))
                ENDIF
                
            ENDIF

        enddo

        RETURN

    END SUBROUTINE ReadFromPopsfilePar

!This will set up the initial walker distribution proportially to the MP1 wavevector.
    SUBROUTINE InitWalkersMP1Par()
        use SystemData , only : tAssumeSizeExcitgen
        INTEGER :: HFConn,error,ierr,MemoryAlloc,VecSlot,nJ(NEl),nStore(6),iExcit,i,j,WalkersonHF,HFPointer
        REAL*8 :: SumMP1Compts,MP2Energy,Compt,Ran2,r
        TYPE(HElement) :: Hij,Hjj,Fjj
        INTEGER , ALLOCATABLE :: MP1Dets(:,:), ExcitgenTemp(:)
        LOGICAL , ALLOCATABLE :: MP1Sign(:)
        REAL*8 , ALLOCATABLE :: MP1Comps(:)
        INTEGER :: MP1DetsTag,MP1SignTag,MP1CompsTag,SumWalkersonHF,ExcitLength,iMaxExcit
        CHARACTER(len=*), PARAMETER :: this_routine='InitWalkersMP1Par'
        LOGICAL :: First,TurnBackAssumeExGen
        
    
!Set the maximum number of walkers allowed
        MaxWalkers=NINT(MemoryFac*InitWalkers)
        MaxWalkersExcit=NINT(MemoryFacExcit*InitWalkers)
        WRITE(6,"(A,F14.5)") "Memory Factor for walkers is: ",MemoryFac
        WRITE(6,"(A,F14.5)") "Memory Factor for excitation generators is: ",MemoryFacExcit
        WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node of: ",MaxWalkers
        WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for excitgens of: ",MaxWalkersExcit
                                            
!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
!Allocate memory to hold walkers
        ALLOCATE(WalkVecDets(NEl,MaxWalkersExcit),stat=ierr)
        CALL LogMemAlloc('WalkVecDets',MaxWalkersExcit*NEl,4,this_routine,WalkVecDetsTag,ierr)
        WalkVecDets(1:NEl,1:MaxWalkersExcit)=0
        ALLOCATE(WalkVec2Dets(NEl,MaxWalkersExcit),stat=ierr)
        CALL LogMemAlloc('WalkVec2Dets',MaxWalkersExcit*NEl,4,this_routine,WalkVec2DetsTag,ierr)
        WalkVec2Dets(1:NEl,1:MaxWalkersExcit)=0
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
        WalkVecH=0.d0
        ALLOCATE(WalkVec2H(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec2H',MaxWalkers,8,this_routine,WalkVec2HTag,ierr)
        WalkVec2H=0.d0
        
        MemoryAlloc=((2*NEl+8)*MaxWalkers)*4    !Memory Allocated in bytes

        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
            ALLOCATE(HashArray(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('HashArray',MaxWalkers,8,this_routine,HashArrayTag,ierr)
            HashArray(:)=0
            ALLOCATE(Hash2Array(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('Hash2Array',MaxWalkers,8,this_routine,Hash2ArrayTag,ierr)
            Hash2Array(:)=0
            ALLOCATE(IndexTable(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('IndexTable',MaxWalkers,4,this_routine,IndexTableTag,ierr)
            IndexTable(1:MaxWalkers)=0
            ALLOCATE(Index2Table(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('Index2Table',MaxWalkers,4,this_routine,Index2TableTag,ierr)
            Index2Table(1:MaxWalkers)=0
            ALLOCATE(ProcessVec(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('ProcessVec',MaxWalkers,4,this_routine,ProcessVecTag,ierr)
            ProcessVec(1:MaxWalkers)=0
            ALLOCATE(Process2Vec(MaxWalkers),stat=ierr)
            CALL LogMemAlloc('Process2Vec',MaxWalkers,4,this_routine,Process2VecTag,ierr)
            Process2Vec(1:MaxWalkers)=0

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

!Now calculate MP1 components - allocate memory for doubles
        CALL GetSymExcitCount(HFExcit%ExcitData,HFConn)
        HFConn=HFConn+1     !Add on one for the HF Det itself

        ALLOCATE(MP1Comps(HFConn),stat=ierr)    !This will store the cumulative absolute values of the mp1 wavevector components
        CALL LogMemAlloc('MP1Comps',HFConn,8,this_routine,MP1CompsTag,ierr)
        MP1Comps=0.d0
        ALLOCATE(MP1Dets(NEl,HFConn),stat=ierr)
        CALL LogMemAlloc('MP1Dets',HFConn*NEl,4,this_routine,MP1DetsTag,ierr)
        MP1Dets(1:NEl,1:HFConn)=0
        ALLOCATE(MP1Sign(HFConn),stat=ierr)
        CALL LogMemAlloc('MP1Sign',HFConn,4,this_routine,MP1SignTag,ierr)
        
!HF Compt. of MP1 is 1
        MP1Dets(1:NEl,1)=HFDet(1:NEl)
        MP1Comps(1)=1.D0
        MP1Sign(1)=.true.

        SumMP1Compts=1.D0   !Initialise the sum of the MP1 wavevector components
        VecSlot=2           !This is the next free slot in the MP1 arrays
        MP2Energy=0.D0      !Calculate the MP2 energy as we go, since the shift will be set to this

!If tAssumeSizeExcitgens is on, then we cannot enumerate all determinants. Regenerate HF excitgen, turning of tAssumeSizeExcitgen if on.
        IF(tAssumeSizeExcitgen) THEN
            TurnBackAssumeExGen=.true.
            tAssumeSizeExcitgen=.false.
        ELSE
            TurnBackAssumeExGen=.false.
        ENDIF

!Setup excit generators for HF Determinant
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitLength,nJ,iMaxExcit,0,nStore,2)
        ALLOCATE(ExcitGenTemp(ExcitLength),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("InitWalkersMP1","Problem allocating excitation generator")
        ExcitGenTemp(1)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGenTemp,nJ,iMaxExcit,0,nStore,2)

        do while(.true.)
!Generate double excitations
            CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.false.,ExcitGenTemp,nJ,iExcit,0,nStore,2)
            IF(nJ(1).eq.0) EXIT
            IF(iExcit.ne.2) THEN
                CALL Stop_All("InitWalkersMP1","Error - excitations other than doubles being generated in MP1 wavevector code")
            ENDIF

            Hij=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
            CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fjj)

            Compt=real(Hij%v,r2)/(Fii-(REAL(Fjj%v,r2)))
            IF(Compt.lt.0.D0) THEN
                MP1Sign(VecSlot)=.false.
            ELSE
                MP1Sign(VecSlot)=.true.
            ENDIF
            MP1Dets(1:NEl,VecSlot)=nJ(:)
            MP1Comps(VecSlot)=MP1Comps(VecSlot-1)+abs(Compt)
            SumMP1Compts=SumMP1Compts+abs(Compt)
            MP2Energy=MP2Energy+((real(Hij%v,r2))**2)/(Fii-(REAL(Fjj%v,r2)))

            VecSlot=VecSlot+1

        enddo

        DEALLOCATE(ExcitGenTemp)

        WRITE(6,"(A,F15.7,A)") "Sum of absolute components of MP1 wavefunction is ",SumMP1Compts," with the HF being 1."

        VecSlot=VecSlot-1

!Total components is VecSlot
        IF(MP1Comps(VecSlot).ne.SumMP1Compts) THEN
            CALL Stop_All("InitWalkersMP1","Error in calculating sum of MP1 components")
        ENDIF

        DiagSft=MP2Energy
        MP2Energy=MP2Energy+Hii
        WRITE(6,"(A,F16.7,A,F16.7)") "MP2 energy is ",MP2Energy," which the initial shift has been set to: ",DiagSft

        WalkersonHF=0       !Find the number of walkers we are assigning to HF

        do j=1,InitWalkers

            i=1

            r=Ran2(Seed)*SumMP1Compts       !Choose the double that this walker wants to be put...
            do while(r.gt.MP1Comps(i))

                i=i+1

                IF(i.gt.VecSlot) THEN
                    CALL Stop_All("InitWalkersMP1","Assigning walkers stochastically has been performed incorrectly")
                ENDIF

            enddo

            IF(i.eq.1) THEN
!If we are at HF, then we do not need to calculate the information for the walker        
                WalkersonHF=WalkersonHF+1
                CurrentDets(1:NEl,j)=HFDet(:)
                CurrentIC(j)=0
                CurrentSign(j)=.true.
                CurrentH(j)=0.D0
                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    HashArray(j)=HFHash
                ENDIF
            ELSE
!We are at a double excitation - we need to calculate most of this information...
                CurrentDets(1:NEl,j)=MP1Dets(1:NEl,i)
                CurrentIC(j)=2
                CurrentSign(j)=MP1Sign(i)
                Hjj=GetHElement2(MP1Dets(1:NEl,i),MP1Dets(1:NEl,i),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)     !Find the diagonal element
                CurrentH(j)=real(Hjj%v,r2)-Hii
                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    HashArray(j)=CreateHash(MP1Dets(1:NEl,i))
                ENDIF
                
            ENDIF

        enddo
        
        CALL MPI_Reduce(WalkersonHF,SumWalkersonHF,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)

        IF(iProcIndex.eq.Root) THEN
            WRITE(6,"(A,I12,A,I12,A)") "Out of ",InitWalkers*nProcessors," initial walkers allocated, ",SumWalkersonHF," of them are situated on the HF determinant."
        ENDIF
        AllNoatHF=SumWalkersonHF
        AllNoatDoubs=(InitWalkers*nProcessors)-SumWalkersonHF

!Deallocate MP1 data
        DEALLOCATE(MP1Comps)
        CALL LogMemDealloc(this_routine,MP1CompsTag)
        DEALLOCATE(MP1Dets)
        CALL LogMemDealloc(this_routine,MP1DetsTag)
        DEALLOCATE(MP1Sign)
        CALL LogMemDealloc(this_routine,MP1SignTag)

!We have finished enumerating all determinants from HF - turn back on assume size excitgen if it was off.
        IF(TurnBackAssumeExGen) THEN
            tAssumeSizeExcitgen=.true.
        ENDIF

        WRITE(6,"(A,F14.6,A)") "Initial memory (without excitgens + temp arrays) consists of : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
        WRITE(6,*) "Initial memory allocation sucessful..."
        CALL FLUSH(6)

!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
        IF(.not.TRegenExcitgens) THEN
            ALLOCATE(WalkVecExcits(MaxWalkersExcit),stat=ierr)
            ALLOCATE(WalkVec2Excits(MaxWalkersExcit),stat=ierr)
            ALLOCATE(ExcitGens(MaxWalkersExcit),stat=ierr)
            ALLOCATE(FreeIndArray(MaxWalkersExcit),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("InitFCIMMCCalcPar","Error in allocating walker excitation generators")

            do j=1,MaxWalkersExcit
                NULLIFY(WalkVecExcits(j)%PointToExcit)
                NULLIFY(WalkVec2Excits(j)%PointToExcit)
                ExcitGens(j)%nPointed=0
                FreeIndArray(j)=j       !All points are initially free
            enddo

!Allocate pointers to the correct excitation arrays
            CurrentExcits=>WalkVecExcits
            NewExcits=>WalkVec2Excits

!Initialise the first Free Excitgens indices...initially whole list is free
            BackOfList=1    !I.e. first allocation should be at ExcitGens(FreeIndArray(1))
            FrontOfList=1   !i.e. first free index should be put at FreeIndArray(1)

            First=.true.
            do j=1,InitWalkers
!Copy the HF excitation generator accross to each initial particle
                IF(CurrentIC(j).eq.0) THEN
!We are at HF - we can save the excitgen
                    IF(First) THEN
                        CALL SetupExitgenPar(HFDet,CurrentExcits(j))
                        HFPointer=j
                        First=.false.
                    ELSE
                        CALL CopyExitGenPar(CurrentExcits(HFPointer),CurrentExcits(j),.false.)
                    ENDIF
                ELSE
                    CurrentExcits(j)%PointToExcit=>null()
                ENDIF
            enddo
            MemoryAlloc=((HFExcit%nExcitMemLen)+2)*4*MaxWalkersExcit

            WRITE(6,"(A,F14.6,A)") "Probable maximum memory for excitgens is : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
            WRITE(6,*) "Initial allocation of excitation generators successful..."
        ELSE
            WRITE(6,*) "Excitation generators will not be stored, but regenerated each time they are needed..."
        ENDIF
        WRITE(6,"(A,F14.6,A)") "Temp Arrays for annihilation cannot be more than : ",REAL(MaxWalkers*12,r2)/1048576.D0," Mb/Processor"
        CALL FLUSH(6)
        
!TotWalkers contains the number of current walkers at each step
        TotWalkers=InitWalkers
        TotWalkersOld=InitWalkers
!Initialise global variables for calculation on the root node
        IF(iProcIndex.eq.root) THEN
            AllTotWalkers=InitWalkers*nProcessors
            AllTotWalkersOld=InitWalkers*nProcessors
        ENDIF

        RETURN

    END SUBROUTINE InitWalkersMP1Par

!This routine flips the sign of all particles on the node
    SUBROUTINE FlipSign()
        INTEGER :: i

        do i=1,TotWalkers
            CurrentSign(i)=.not.CurrentSign(i)
        enddo
!Reverse the flag for whether the sign of the particles has been flipped so the ACF can be correctly calculated
        TFlippedSign=.not.TFlippedSign
        RETURN
    
    END SUBROUTINE FlipSign

!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
    INTEGER FUNCTION AttemptCreatePar(DetCurr,WSign,nJ,Prob,IC,Ex,tParity)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,StoreNumTo,StoreNumFrom,DetLT,i,ExtraCreate,Ex(2,2)
        LOGICAL :: WSign,tParity
        REAL*8 :: Prob,Ran2,rat
        TYPE(HElement) :: rh,rhcheck

!Calculate off diagonal hamiltonian matrix element between determinants
!        rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
        rh=GetHElement4(DetCurr,nJ,IC,Ex,tParity)

!        rhcheck=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!        IF(rh%v.ne.rhcheck%v) THEN
!            WRITE(6,*) "DetCurr: ",DetCurr(:)
!            WRITE(6,*) "nJ: ",nJ(:)
!            WRITE(6,*) "EX: ",Ex(1,:),Ex(2,:)
!            WRITE(6,*) "tParity: ",tParity
!            STOP
!        ENDIF

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

!This is a function which tells us whether to annihilate a particle on a determinant if it is the only one there.
!Hopefully this will simulate the annihilation rate to some extent and stop the shift from becoming too negative.
!The only variable it relies on is the approximate density of particles for the given excitation level - ExcitDensity
    LOGICAL FUNCTION AttemptLocalAnn(ExcitDensity)
        REAL*8 :: ExcitDensity,Ran2,AnnProb

!The function can initially be a simple Tau*EXP(-Lambda*ExcitDensity)
!        AnnProb=Tau*EXP(-Lambda*ExcitDensity)
!        AnnProb=Lambda/ExcitDensity
        AnnProb=Lambda/ExcitDensity
        IF(Ran2(Seed).lt.AnnProb) THEN
!Particle is annihilated
            AttemptLocalAnn=.true.
        ELSE
!Particle survives
!            WRITE(6,*) "Particle survives with prob = ",AnnProb
            AttemptLocalAnn=.false.
        ENDIF
        RETURN

    END FUNCTION AttemptLocalAnn

!This function tells us whether we should kill the particle at determinant DetCurr
!If also diffusing, then we need to know the probability with which we have spawned. This will reduce the death probability.
!The function allows multiple births(if +ve shift) or deaths from the same particle.
!The returned number is the number of deaths if positive, and the number of births if negative.
    INTEGER FUNCTION AttemptDiePar(DetCurr,Kii,IC)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),iKill,IC
!        TYPE(HElement) :: rh,rhij
        REAL*8 :: Ran2,rat,Kii

!Calculate the diagonal hamiltonian matrix element for the determinant
!        rh=GetHElement2(DetCurr,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!Subtract from the diagonal the value of the lowest hamiltonian matrix element
!        rh=rh-Hii

!Subtract the current value of the shift and multiply by tau
        IF(TFixShiftShell.and.(.not.TSinglePartPhase)) THEN
!            IF((IC.eq.0).or.(IC.eq.2)) THEN
            IF(IC.le.ShellFix) THEN
                rat=Tau*(Kii-FixShift)
            ELSE
                rat=Tau*(Kii-DiagSft)
            ENDIF
        ELSE
            rat=Tau*(Kii-DiagSft)
        ENDIF

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
!            CreateHash=13*CreateHash+i*DetCurr(i)
            CreateHash=(1099511628211*CreateHash)+i*DetCurr(i)
            
!            CreateHash=mod(1099511628211*CreateHash,2**64)
!            CreateHash=XOR(CreateHash,DetCurr(i))
        enddo
!        WRITE(6,*) CreateHash
        RETURN

    END FUNCTION CreateHash

!This routine copies an excitation generator from origExcit to NewExit, if the original claims that it is for the correct determinant
    SUBROUTINE CopyExitgenPar(OrigExit,NewExit,DelOldCopy)
        TYPE(ExcitPointer) :: OrigExit,NewExit
        LOGICAL :: DelOldCopy
        INTEGER :: ierr
        
        IF(ASSOCIATED(NewExit%PointToExcit)) THEN
            CALL Stop_All("CopyExitgenPar","Trying to copy an excitation, but new pointer is already associated.")
        ENDIF
        IF(.not.ASSOCIATED(OrigExit%PointToExcit)) THEN
!We have not got a new pointer - it hasn't been created yet.
            RETURN
        ENDIF

        NewExit%PointToExcit => OrigExit%PointToExcit
        NewExit%IndexinExArr=OrigExit%IndexinExArr

        IF(DelOldCopy) THEN
!Delete the old excitation - i.e. we are moving excitation generators, rather than copying them.
            OrigExit%PointToExcit=>null()
        ELSE
! We are copying, so increment the number of objects pointing at the excitgen.
            EXCITGENS(NewExit%IndexinExArr)%nPointed=EXCITGENS(NewExit%IndexinExArr)%nPointed+1
        ENDIF

        RETURN

    END SUBROUTINE CopyExitgenPar

    SUBROUTINE SetupExitgenPar(nI,ExcitGen)
        TYPE(ExcitPointer) :: ExcitGen
        INTEGER :: ierr,iMaxExcit,nExcitMemLen,nJ(NEl),MinIndex,i
        INTEGER :: nI(NEl),nStore(6)

        IF(ASSOCIATED(ExcitGen%PointToExcit)) THEN
!The determinant already has an associated excitation generator set up.
            RETURN
        ELSE

!First, we need to find the next free element in the excitgens array...
!This is simply FreeIndArray(BackOfList)
            MinIndex=FreeIndArray(BackOfList)
!Increment BackOfList in a circular fashion.
            IF(BackOfList.eq.MaxWalkersExcit) THEN
                BackOfList=1
            ELSE
                BackOfList=BackOfList+1
            ENDIF

            IF(associated(ExcitGens(MinIndex)%ExcitData)) THEN
                CALL Stop_All("SetupExitgenPar","Index chosen to create excitation generator is not free.")
            ENDIF

!MinIndex is the array element we want to point our new excitation generator to.
!Setup excit generators for this determinant
            iMaxExcit=0
            nStore(1:6)=0
            CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,EXCITGENS(MinIndex)%nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
            ALLOCATE(EXCITGENS(MinIndex)%ExcitData(EXCITGENS(MinIndex)%nExcitMemLen),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
            EXCITGENS(MinIndex)%ExcitData(1)=0
            CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,EXCITGENS(MinIndex)%ExcitData,nJ,iMaxExcit,0,nStore,3)

!Indicate that the excitation generator is now correctly allocated and pointed to by one particle.
            EXCITGENS(MinIndex)%nPointed=1

!Now point Excitgen to this value
            ExcitGen%PointToExcit=>EXCITGENS(MinIndex)%ExcitData
            ExcitGen%IndexinExArr=MinIndex

        ENDIF

    END SUBROUTINE SetupExitgenPar
                
    SUBROUTINE DissociateExitgen(Exitgen)
        TYPE(ExcitPointer) :: Exitgen
        INTEGER :: ind

        IF(.not.ASSOCIATED(Exitgen%PointToExcit)) THEN
            RETURN
        ENDIF
        Ind=Exitgen%IndexinExArr

        IF(Excitgens(Ind)%nPointed.eq.1) THEN
!We want to delete this excitgen.
            DEALLOCATE(Excitgens(Ind)%ExcitData)
            Excitgens(Ind)%nPointed=0

!Add removed excitgen to the front of the free index list
            FreeIndArray(FrontOfList)=Ind
!Increment frontoflist in a circular fashion.
            IF(FrontOfList.eq.MaxWalkersExcit) THEN
                FrontOfList=1
            ELSE
                FrontOfList=FrontOfList+1
            ENDIF

        ELSE
            Excitgens(Ind)%nPointed=Excitgens(Ind)%nPointed-1
        ENDIF
        Exitgen%PointToExcit=>null()    !Point to null to show that it is now free.

    END SUBROUTINE DissociateExitgen

!This routine gets a random excitation for when we want to generate the excitation generator on the fly, then chuck it.
    SUBROUTINE GetPartRandExcitPar(DetCurr,nJ,Seed,IC,Frz,Prob,iCount,ExcitLevel,Ex,tParity)
        INTEGER :: DetCurr(NEl),nJ(NEl),Seed,IC,Frz,iCount,iMaxExcit,nStore(6),MemLength,ierr
        INTEGER :: Excitlevel,Ex(2,2)
        REAL*8 :: Prob
        LOGICAL :: tParity
        INTEGER , ALLOCATABLE :: ExcitGenTemp(:)

        IF(ExcitLevel.eq.0) THEN
!            CALL GenRandSymExcitIt3(DetCurr,HFExcit%ExcitData,nJ,Seed,IC,Frz,Prob,iCount)
            CALL GenRandSymExcitIt4(DetCurr,HFExcit%ExcitData,nJ,Seed,IC,Frz,Prob,iCount,Ex,tParity)
            RETURN
        ENDIF
            
!Need to generate excitation generator to find excitation.
!Setup excit generators for this determinant 
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,MemLength,nJ,iMaxExcit,0,nStore,3)
        ALLOCATE(ExcitGenTemp(MemLength),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
        ExcitGenTemp(1)=0
        CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGenTemp,nJ,iMaxExcit,0,nStore,3)

!Now generate random excitation
!        CALL GenRandSymExcitIt3(DetCurr,ExcitGenTemp,nJ,Seed,IC,Frz,Prob,iCount)
        CALL GenRandSymExcitIt4(DetCurr,ExcitGenTemp,nJ,Seed,IC,Frz,Prob,iCount,Ex,tParity)

!Deallocate when finished
        DEALLOCATE(ExcitGenTemp)

        RETURN

    END SUBROUTINE GetPartRandExcitPar
    

!This routine will move walkers between the processors, in order to balance the number of walkers on each node.
!This could be made slightly faster by using an MPI_Reduce, rather than searching for the min and max by hand...
    SUBROUTINE BalanceWalkersonProcs()
        INTEGER :: i,ProcWalkers(0:nProcessors-1),error
        INTEGER :: MinWalkers(2),MaxWalkers(2)      !First index gives the number, second the rank of the processor
        REAL*8 :: MeanWalkers,MidWalkers
        INTEGER :: WalktoTransfer(4)                !This gives information about how many walkers to transfer from and to
        INTEGER :: Transfers        !This is the number of transfers of walkers necessary
        INTEGER :: IndexFrom,IndexTo,j,Stat(MPI_STATUS_SIZE),Tag,l

        Tag=123         !Set tag for sends
!        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        IF(iProcIndex.eq.root) THEN
            WRITE(6,*) "Moving walkers between nodes in order to balance the load..."
        ENDIF
        
!First, it is necessary to find the number of walkers on each processor and gather the info to the root.
        CALL MPI_Gather(TotWalkers,1,MPI_INTEGER,ProcWalkers,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
        Transfers=0
        WalktoTransfer(4)=100           !Set to large number to get ball rolling...

        do while(WalktoTransfer(4).gt.1)    !This is the range of walkers accross the processors. Continue reducing it until it is 1.
                
!            do j=1,TotWalkers+3
!                WRITE(6,"(8I4)") CurrentDets(:,j)
!                CALL FLUSH(6)
!            enddo

            IF(iProcIndex.eq.root) THEN
!Let the root processor go through and decide which walkers to move from where...
!Initialise variables with the root processors own values
!                WRITE(6,*) "PROCWALKERS: ",ProcWalkers
!                CALL FLUSH(6)
                MaxWalkers(1)=ProcWalkers(0)   !This should be the same as the totwalkers on the root node.
                MaxWalkers(2)=0                !This indicates that currently the maximum number of walkers resides on the root node
                MinWalkers(1)=ProcWalkers(0)
                MinWalkers(2)=0
                MeanWalkers=REAL(ProcWalkers(0),r2)
                do i=1,nProcessors-1
!Find the minimum, maximum and mean of the walkers accross the nodes, and their node
                    IF(ProcWalkers(i).gt.MaxWalkers(1)) THEN
                        MaxWalkers(1)=ProcWalkers(i)
                        MaxWalkers(2)=i
                    ELSEIF(ProcWalkers(i).lt.MinWalkers(1)) THEN
                        MinWalkers(1)=ProcWalkers(i)
                        MinWalkers(2)=i
                    ENDIF
                    MeanWalkers=MeanWalkers+REAL(ProcWalkers(i),r2)
                enddo
                WalktoTransfer(4)=MaxWalkers(1)-MinWalkers(1)
                IF(Transfers.eq.0) THEN
                    WRITE(6,*) "Initial range of walkers is: ",WalktoTransfer(4)
                ELSE
!                    WRITE(6,*) "After ",Transfers," walker transfers, range of walkers is: ",WalktoTransfer(4)
                    IF(Transfers.gt.(10*(nProcessors**2))) THEN
                        CALL Stop_All("BalanceWalkersonProcs","Too many transfers required to balance nodes - Problem here...")
                    ENDIF
!                    CALL FLUSH(6)
                ENDIF
!                WRITE(6,*) "TOTAL WALKERS = ",MeanWalkers
!                CALL FLUSH(6)
                MeanWalkers=MeanWalkers/REAL(nProcessors,r2)
!                WRITE(6,*) "MEAN WALKERS = ", MeanWalkers
                MidWalkers=REAL((MaxWalkers(1)+MinWalkers(1)),r2)/2.D0  !THis is the mean of the two processors to exchange walkers

! Now it is necessary to decide what walkers go where...want to transfer from maxwalkers to minwalkers
! i.e. want to transfer walkers from MaxWalkers(2) to MinWalkers(2). The number to transfer is given by 
! min[ abs(maxwalkers(1)-MidWalkers), abs(minwalkers(1)-MidWalkers) ]
! Then we want to update ProcWalkers on the root and recalculate Max,Min,Mean and range.
! Broadcast range, since this is the termination criterion
! Repeat this until the range .le. 1
! WalktoTransfer(1) is the number of walkers to transfer
! WalktoTransfer(2) is the rank of the processor to transfer them from
! WalktoTransfer(3) is the rank of the processor to transfer them to
! WalktoTransfer(4) is the range of the different walkers on the processors - this will determine when to stop

                IF(WalktoTransfer(4).le.1) THEN
!This is the termination criteria, and will not transfer any walkers in this iteration.
                    WalktoTransfer(1)=0
                    WalktoTransfer(2)=-1
                    WalktoTransfer(3)=-1
                ELSE
!                    WalktoTransfer(1)=MIN(NINT(ABS(REAL(MaxWalkers(1),r2)-MeanWalkers)),NINT(ABS(REAL(MinWalkers(1),r2)-MeanWalkers)))
                    WalktoTransfer(1)=MIN(NINT(ABS(REAL(MaxWalkers(1),r2)-MidWalkers)),NINT(ABS(REAL(MinWalkers(1),r2)-MidWalkers)))
!                    WRITE(6,*) MaxWalkers(:),MinWalkers(:),MidWalkers
!                    CALL FLUSH(6)
                    WalktoTransfer(2)=MaxWalkers(2)
                    WalktoTransfer(3)=MinWalkers(2)
                ENDIF

            ENDIF

!The information about what to transfer now needs to be broadcast to all processors
!            CALL MPI_Barrier(MPI_COMM_WORLD,error)
            CALL MPI_BCast(WalktoTransfer,4,MPI_INTEGER,root,MPI_COMM_WORLD,error)

!            WRITE(6,*) "WALKTOTRANSFER: ",WalktoTransfer(:)
!            CALL FLUSH(6)

            IF(WalkToTransfer(2).eq.iProcIndex) THEN
!This processor wants to transfer walkers to WalkToTransfer(3).
                IndexFrom=TotWalkers-(WalktoTransfer(1)-1)          !We want to transfer the last WalktoTransfer(1) walkers
!                WRITE(6,*) "TRANSFER: ",Transfers,WalkToTransfer(:),IndexFrom,TotWalkers
!                CALL FLUSH(6)
                
                CALL MPI_Send(CurrentDets(:,IndexFrom:TotWalkers),WalktoTransfer(1)*NEl,MPI_INTEGER,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
                CALL MPI_Send(CurrentSign(IndexFrom:TotWalkers),WalktoTransfer(1),MPI_LOGICAL,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
                CALL MPI_Send(CurrentIC(IndexFrom:TotWalkers),WalktoTransfer(1),MPI_INTEGER,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
                CALL MPI_Send(CurrentH(IndexFrom:TotWalkers),WalktoTransfer(1),MPI_DOUBLE_PRECISION,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
!It seems like too much like hard work to send the excitation generators accross, just let them be regenerated on the other side...
!However, we do need to indicate that these the excitgens are no longer being pointed at.
                IF(.not.TRegenExcitgens) THEN
                    do l=IndexFrom,TotWalkers
!Run through the list of walkers, and make sure that their excitgen is removed.
                        CALL DissociateExitgen(CurrentExcits(l))
                    enddo
                ENDIF

                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) CALL MPI_Send(HashArray(IndexFrom:TotWalkers),WalktoTransfer(1),MPI_DOUBLE_PRECISION,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
!                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) CALL MPI_Send(HashArray(IndexFrom:TotWalkers),WalktoTransfer(1),mpilongintegertype,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
                
                TotWalkers=TotWalkers-WalktoTransfer(1)         !Update TotWalkers on this node to reflect that we have lost some
            
            ELSEIF(WalkToTransfer(3).eq.iProcIndex) THEN
!This processor wants to receive walkers from WalktoTransfer(2)
                IndexFrom=TotWalkers+1
                IndexTo=TotWalkers+WalktoTransfer(1)
!                WRITE(6,*) "RECEIVING: ",Transfers,WalkToTransfer(:),IndexFrom,IndexTo
!                CALL FLUSH(6)

                CALL MPI_Recv(CurrentDets(:,IndexFrom:IndexTo),WalktoTransfer(1)*NEl,MPI_INTEGER,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(CurrentSign(IndexFrom:IndexTo),WalktoTransfer(1),MPI_LOGICAL,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(CurrentIC(IndexFrom:IndexTo),WalktoTransfer(1),MPI_INTEGER,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(CurrentH(IndexFrom:IndexTo),WalktoTransfer(1),MPI_DOUBLE_PRECISION,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) CALL MPI_Recv(HashArray(IndexFrom:IndexTo),WalktoTransfer(1),MPI_DOUBLE_PRECISION,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
!Also need to indicate that the excitation generators are no longer useful...
                IF(.not.TRegenExcitgens) THEN
                    do j=IndexFrom,IndexTo
                        CurrentExcits(j)%PointToExcit=>null()
                    enddo
                ENDIF

                TotWalkers=TotWalkers+WalktoTransfer(1)         !Update to show the new number of walkers

            ENDIF
            
            IF(WalktoTransfer(1).gt.0) THEN
                IF(iProcIndex.eq.root) THEN
!Update the number of transfers and the new ProcWalkers array...
                    Transfers=Transfers+1
                    ProcWalkers(WalktoTransfer(2))=ProcWalkers(WalktoTransfer(2))-WalktoTransfer(1)
                    ProcWalkers(WalktoTransfer(3))=ProcWalkers(WalktoTransfer(3))+WalktoTransfer(1)
                ENDIF
            ENDIF

        enddo       !loop over transfers

        IF(iProcIndex.eq.root) THEN
            WRITE(6,*) "Transfer of walkers finished. Number of transfers needed: ",Transfers
        ENDIF

        TBalanceNodes=.false.

        RETURN

    END SUBROUTINE BalanceWalkersonProcs
    
    SUBROUTINE DeallocFCIMCMemPar()
        INTEGER :: i,error,length,temp
        CHARACTER(len=*), PARAMETER :: this_routine='DeallocFciMCMemPar'
        CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: message
            
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
        IF(TLocalAnnihilation) THEN
            DEALLOCATE(ApproxExcitDets)
            DEALLOCATE(PartsinExcitLevel)
        ENDIF

        DEALLOCATE(HFDet)
        CALL LogMemDealloc(this_routine,HFDetTag)
        DEALLOCATE(HFExcit%ExcitData)
        IF(.not.TRegenExcitgens) THEN
            do i=1,MaxWalkersExcit
                CALL DissociateExitgen(WalkVecExcits(i))

!                IF(Allocated(WalkVecExcits(i)%ExcitData)) THEN
!                    DEALLOCATE(WalkVecExcits(i)%ExcitData)
!                ENDIF
!                IF(Allocated(WalkVec2Excits(i)%ExcitData)) THEN
!                    DEALLOCATE(WalkVec2Excits(i)%ExcitData)
!                ENDIF
            enddo
            DEALLOCATE(WalkVecExcits)
            DEALLOCATE(WalkVec2Excits)
            DEALLOCATE(ExcitGens)
            DEALLOCATE(FreeIndArray)
        ENDIF

!There seems to be some problems freeing the derived mpi type.
!        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
!Free the mpi derived type that we have created for the hashes.
!            CALL MPI_Type_free(mpilongintegertype,error)
!            IF(error.ne.MPI_SUCCESS) THEN
!                CALL MPI_Error_string(error,message,length,temp)
!                IF(temp.ne.MPI_SUCCESS) THEN
!                    WRITE(6,*) "REALLY SERIOUS PROBLEMS HERE!",temp
!                    CALL FLUSH(6)
!                ENDIF
!                WRITE(6,*) message(1:length)
!                WRITE(6,*) "ERROR FOUND"
!                CALL FLUSH(6)
!            ENDIF
!        ENDIF
        RETURN

    END SUBROUTINE DeallocFCIMCMemPar

    SUBROUTINE WriteFCIMCStats()

        IF(iProcIndex.eq.root) THEN

            IF(TLocalAnnihilation) THEN
!TotWalkersOld is the number of walkers last time the shift was changed
                WRITE(15,"(I12,G16.7,I9,G16.7,I12,4I11,G16.7,2I10,2F13.5,2G13.5,2I6)") Iter+PreviousCycles,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,   &
 &                  AllTotWalkers,AllLocalAnn,AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AllNoatHF,AllNoatDoubs,AllAvSign,AllAvSignHFD,AccRat,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
                WRITE(6,"(I12,G16.7,I9,G16.7,I12,4I11,G16.7,2I10,2F13.5,2G13.5,2I6)") Iter+PreviousCycles,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,    &
 &                  AllTotWalkers,AllLocalAnn,AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AllNoatHF,AllNoatDoubs,AllAvSign,AllAvSignHFD,AccRat,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel

            ELSE
                WRITE(15,"(I12,G16.7,I9,G16.7,I12,3I11,G16.7,2I10,2F13.5,2G13.5,2I6)") Iter+PreviousCycles,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,   &
 &                  AllTotWalkers,AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AllNoatHF,AllNoatDoubs,AllAvSign,AllAvSignHFD,AccRat,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
                WRITE(6,"(I12,G16.7,I9,G16.7,I12,3I11,G16.7,2I10,2F13.5,2G13.5,2I6)") Iter+PreviousCycles,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,    &
 &                  AllTotWalkers,AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AllNoatHF,AllNoatDoubs,AllAvSign,AllAvSignHFD,AccRat,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
            ENDIF
            
            CALL FLUSH(6)
            CALL FLUSH(15)
            
        ENDIF

        RETURN

    END SUBROUTINE WriteFCIMCStats

!This routine simply reduces the histogrammed determinants, and prints out their population for the given iteration
    SUBROUTINE WriteHistogrammedDets()
        INTEGER :: AllWeightatDets(NoAutoDets),error,k
        
!First, collate the information onto the root
        AllWeightatDets(:)=0
        CALL MPI_Reduce(WeightatDets(1:NoAutoDets),AllWeightatDets(1:NoAutoDets),NoAutoDets,MPI_INTEGER,MPI_SUM,root,MPI_COMM_WORLD,error)

        IF(iProcIndex.eq.root) THEN
            WRITE(44,"(2I8)",advance='no') Iter,AllWeightatDets(1)
            do k=2,NoAutoDets-1
                WRITE(44,"(I8)",advance='no') AllWeightatDets(k)
            enddo
            WRITE(44,"(I8)") AllWeightatDets(NoAutoDets)
        ENDIF
        WeightatDets(:)=0   !Rezero the array for the next iteration
    END SUBROUTINE WriteHistogrammedDets

    
!This routine calculates the autocorrelation function for the determinants listed in AutoCorrDets array, with lags from iLagMin to iLagMax in steps of iLagStep
!and writes it to a file ACF. The calculation is done in serial. It is not trivial to turn it into an efficient parallel algorithm
!The autocorrelation function is now calculated in a seperate standalone program. The occupations at each point are simply printed out.
!    SUBROUTINE CalcAutoCorr()
!        INTEGER :: i,j,k,error
!        REAL*8 :: Means(NoAutoDets),NVar(NoAutoDets),ACF(NoAutoDets)!,AllACF(NoAutoDets)
!        REAL*8 :: Norm(NoAutoDets)
!!        REAL*8 :: TestVar(NoAutoDets)
!        INTEGER :: SumSquares(NoAutoDets),SumWeights(NoAutoDets),ierr
!        INTEGER , ALLOCATABLE :: AllWeightatDets(:,:)
!        INTEGER :: AllWeightatDetsTag=0,NoCounts
!        CHARACTER(len=*), PARAMETER :: this_routine='CalcAutoCorr'
!        
!        CALL set_timer(ACF_Time,30)
!        
!        IF((iLagMax.lt.0).or.(iLagMax.gt.Iter/2)) THEN
!            iLagMax=Iter/2
!        ENDIF
!        IF(iProcIndex.eq.root) THEN
!            OPEN(43,FILE='AutoCorrFunc',STATUS='UNKNOWN')
!        ENDIF
!
!        WRITE(6,"(A,I8,A,I8,A,I8,A)") "Calculating the ACF with lags from ",iLagMin," to ",iLagMax," in steps of ",iLagStep," and writing it to file 'ACF'"
!        WRITE(6,*) "Calculating the ACF for the following determinants:"
!        do i=1,NoAutoDets
!            do j=1,NEl
!                WRITE(6,"(I4)",advance='no') AutoCorrDets(j,i)
!            enddo
!            WRITE(6,*) ""
!        enddo
!        CALL FLUSH(6)
!        CALL MPI_Barrier(MPI_COMM_WORLD,error)
!
!        ALLOCATE(AllWeightatDets(NoAutoDets,Iter),stat=ierr)
!        CALL LogMemAlloc('AllWeightatDets',NoAutoDets*Iter,4,this_routine,AllWeightatDetsTag,ierr)
!        AllWeightatDets(:,:)=0
!!        WRITE(6,*) Iter
!!        do i=1,Iter
!!            do j=1,NoAutoDets
!!                WRITE(6,'(I10)',advance='no') WeightatDets(j,i)
!!            enddo
!!            WRITE(6,*) ""
!!        enddo
!!        WRITE(6,*) "*****************"
!
!        NoCounts=NoAutoDets*Iter
!!Initially, we will calculate this ACF in serial - this will be easier as there are subtle effects in parallel.
!        CALL MPI_Reduce(WeightatDets(1:NoAutoDets,1:Iter),AllWeightatDets(1:NoAutoDets,1:Iter),NoCounts,MPI_INTEGER,MPI_SUM,root,MPI_COMM_WORLD,error)
!!        do i=1,Iter
!!            do j=1,NoAutoDets
!!                WRITE(6,'(I10)',advance='no') AllWeightatDets(j,i)
!!            enddo
!!            WRITE(6,*) ""
!!        enddo
!!        CALL FLUSH(6)
!
!        IF(iProcIndex.eq.root) THEN
!            OPEN(44,FILE='DoublePops',STATUS='UNKNOWN')
!            do i=1,Iter
!                WRITE(44,"(2I8)",advance='no') i,AllWeightatDets(1,i)
!                do k=2,NoAutoDets-1
!                    WRITE(44,"(I8)",advance='no') AllWeightatDets(k,i)
!                enddo
!                WRITE(44,"(I8)") AllWeightatDets(NoAutoDets,i) 
!            enddo
!            CLOSE(44)
!
!!First, we need to calculate the average value for each of the determinants - this could be calculated during the simulation
!!We also want to divide the ACF components by sum_i ((s_i - av(s))^2) , which is N * Var(s) = sum(s^2) - 1/N (sum(x))^2
!            SumWeights(:)=0
!            SumSquares(:)=0
!            NVar(:)=0.D0
!            Means(:)=0.D0
!            Norm(:)=0.D0
!        
!            do i=1,Iter
!                do j=1,NoAutoDets
!                    SumWeights(j)=SumWeights(j)+AllWeightatDets(j,i)
!                    SumSquares(j)=SumSquares(j)+(AllWeightatDets(j,i)**2)
!!                    IF(i.le.(Iter/2)) Norm(j)=Norm(j)+(AllWeightatDets(j,i)**2)
!                    Norm(j)=Norm(j)+(AllWeightatDets(j,i)**2)
!                enddo
!            enddo
!
!!This then needs to be sent to all nodes
!!            CALL MPI_AllReduce(SumWeights,AllSumWeights,NoAutoDets,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,error)
!!            CALL MPI_AllReduce(SumSquares,AllSumSquares,NoAutoDets,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,error)
!
!            do j=1,NoAutoDets
!                Means(j)=REAL(SumWeights(j),r2)/REAL(Iter,r2)
!                NVar(j)=REAL(SumSquares(j),r2)-((REAL(SumWeights(j),r2)**2)/REAL(Iter,r2))
!                WRITE(6,"(A,I4)",advance='no') "Mean+Var for det: ",AutoCorrDets(1,j)
!                do k=2,NEl
!                    WRITE(6,"(I4)",advance='no') AutoCorrDets(k,j)
!                enddo
!                WRITE(6,"(A,2F20.10)") " is: ", Means(j),NVar(j)/REAL(Iter,r2)
!            enddo
!
!!Alternativly, we can calculate the variance seperatly...
!!            TestVar(:)=0.D0
!!            do i=1,Iter
!!                do j=1,NoAutoDets
!!                    TestVar(j)=TestVar(j)+((REAL(AllWeightatDets(j,i),r2)-Means(j))**2)
!!                enddo
!!            enddo
!!            CALL MPI_AllReduce(TestVar,AllTestVar,NoAutoDets,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
!!            WRITE(6,*) "TESTVAR: ",TestVar/REAL(Iter,r2)
!
!!Now we need to calculate the ACF for the desired values of the lag.
!            do i=iLagMin,iLagMax,iLagStep
!!i is now the value of the lag which we are calculating
!
!                ACF(:)=0.D0
!
!                do j=1,(Iter-i)
!!                do j=1,Iter/2
!!j is the run over the values needed
!
!                    do k=1,NoAutoDets
!!k is the run over the desired determinants which to calculate the ACFs
!!                        ACF(k)=ACF(k)+(REAL(AllWeightatDets(k,j),r2)-Means(k))*(REAL(AllWeightatDets(k,j+i),r2)-Means(k))
!                        ACF(k)=ACF(k)+(AllWeightatDets(k,j)*AllWeightatDets(k,j+i))
!                    enddo
!
!                enddo
!
!!Now we need to collate the information from all processors
!!                CALL MPI_Reduce(ACF,AllACF,NoAutoDets,MPI_DOUBLE_PRECISION,MPI_SUM,root,MPI_COMM_WORLD,error)
!
!                do k=1,NoAutoDets
!!Effectivly 'normalise' the ACF by dividing by the variance
!!                    ACF(k)=(ACF(k)/NVar(k))!*REAL(Iter/(Iter-i+0.D0),r2)
!                    ACF(k)=(ACF(k)/REAL(NORM(k),r2))*REAL(Iter/(Iter-i+0.D0),r2)
!!                    ACF(k)=ACF(k)/REAL(NORM(k),8)
!                enddo
!!Write out the ACF
!                WRITE(43,"(I8,F20.10)",advance='no') i,ACF(1)
!                do k=2,NoAutoDets-1
!                    WRITE(43,"(F20.10)",advance='no') ACF(k)
!                enddo
!                WRITE(43,"(F20.10)") ACF(NoAutoDets)
!
!            enddo
!
!            CLOSE(43)
!
!        ENDIF
!        
!        DEALLOCATE(WeightatDets)
!        CALL LogMemDealloc(this_routine,WeightatDetsTag)
!        DEALLOCATE(AllWeightatDets)
!        CALL LogMemDealloc(this_routine,AllWeightatDetsTag)
!        DEALLOCATE(AutoCorrDets)
!        CALL LogMemDealloc(this_routine,AutoCorrDetsTag)
!
!        CALL halt_timer(ACF_Time)
!        
!        RETURN
!
!    END SUBROUTINE CalcAutoCorr
    

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
    REAL*8 :: AvSign           !This is the fraction of positive particles on each node
    REAL*8 :: AvSignHFD
    INTEGER :: SumWalkersCyc    !This is the sum of all walkers over an update cycle on each processor
    REAL*8 :: MeanExcitLevel    
    INTEGER :: MinExcitLevel
    INTEGER :: MaxExcitLevel

!These are the global variables, calculated on the root processor, from the values above
    REAL*8 :: AllGrowRate,AllMeanExcitLevel
    INTEGER :: AllMinExcitLevel,AllMaxExcitLevel,AllTotWalkers,AllTotWalkersOld,AllSumWalkersCyc
    REAL*8 :: AllSumNoatHF,AllSumENum,AllAvSign,AllAvSignHFD

    REAL*8 :: MPNorm        !MPNorm is used if TNodalCutoff is set, to indicate the normalisation of the MP Wavevector

    TYPE(HElement) :: rhii
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
