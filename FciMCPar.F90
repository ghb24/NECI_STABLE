!This is a parallel MPI version of the FciMC code.
!All variables refer to values per processor

! AJWT
! Bringing you a better FciMCPar.  A vision for the future...
!
!   The module now has the same structure with and without PARALLEL being defined.
!   Some routines require MPI and are enclosed in the #ifdef PARALLEL section.  These
!   should have dummy replacements in the #else of this if required.
!   At the end are functions which do not require parallel directives, and are accessible
!   for both parallel and non-parallel.
MODULE FciMCParMod
    use SystemData , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,nMsh,Arr,LMS,NIfD,tHPHF,tListDets
    use SystemData , only : tHub,tReal,tNonUniRandExcits,tMerTwist,tRotatedOrbs,tImportanceSample,tFindCINatOrbs
    use CalcData , only : InitWalkers,NMCyc,DiagSft,Tau,SftDamp,StepsSft,OccCASorbs,VirtCASorbs,tFindGroundDet,tDirectAnnihil
    use CalcData , only : TStartMP1,NEquilSteps,TReadPops,TRegenExcitgens,TFixShiftShell,ShellFix,FixShift,tMultipleDetsSpawn
    use CalcData , only : tConstructNOs,tAnnihilatebyRange,tRotoAnnihil,MemoryFacSpawn,tRegenDiagHEls,tSpawnAsDet
    use CalcData , only : GrowMaxFactor,CullFactor,TStartSinglePart,ScaleWalkers,Lambda,TLocalAnnihilation,tNoReturnStarDets
    use CalcData , only : NDets,RhoApp,TResumFCIMC,TNoAnnihil,MemoryFacPart,TAnnihilonproc,MemoryFacAnnihil,iStarOrbs,tAllSpawnStarDets
    use CalcData , only : FixedKiiCutoff,tFixShiftKii,tFixCASShift,tMagnetize,BField,NoMagDets,tSymmetricField,tStarOrbs,SinglesBias
    use CalcData , only : tHighExcitsSing,iHighExcitsSing,tFindGuide,iGuideDets,tUseGuide,iInitGuideParts,tNoDomSpinCoup
    use CalcData , only : tPrintDominant,iNoDominantDets,MaxExcDom,MinExcDom,tSpawnDominant,tMinorDetsStar
    use CalcData , only : tCCMC,tTruncCAS,iAnnInterval
    use HPHFRandExcitMod , only : FindExcitBitDetSym,GenRandHPHFExcit,GenRandHPHFExcit2Scratch 
    USE Determinants , only : FDet,GetHElement2,GetHElement4
    USE DetCalc , only : ICILevel,nDet,Det,FCIDetIndex
    use GenRandSymExcitNUMod , only : GenRandSymExcitScratchNU,GenRandSymExcitNU
    use IntegralsData , only : fck,NMax,UMat,tPartFreezeCore,NPartFrozen,NHolesFrozen
    USE UMatCache , only : GTID
    USE Logging , only : iWritePopsEvery,TPopsFile,TZeroProjE,iPopsPartEvery,tBinPops,tHistSpawn,iWriteHistEvery,tHistEnergies
    USE Logging , only : NoACDets,BinRange,iNoBins,OffDiagBinRange,OffDiagMax,tPrintSpinCoupHEl!,iLagMin,iLagMax,iLagStep,tAutoCorr
    USE Logging , only : tPrintTriConnections,tHistTriConHels,tPrintHElAccept,tPrintFCIMCPsi,tCalcFCIMCPsi,NHistEquilSteps
    USE Logging , only : tHFPopStartBlock,tIterStartBlock,IterStartBlocking,HFPopStartBlocking,tInitShiftBlocking
    USE SymData , only : nSymLabels
    USE mt95 , only : genrand_real2
    USE Parallel
    USE FciMCData
    IMPLICIT NONE
    SAVE

    contains

#ifdef PARALLEL

    SUBROUTINE FciMCPar(Weight,Energyxw)
        use soft_exit, only : ChangeVars 
        use CalcData, only : iFullSpaceIter
        use UMatCache, only : UMatInd
        use FciMCLoggingMOD , only : PrintTriConnHist,PrintTriConnHElHist,FinaliseBlocking,FinaliseShiftBlocking
        use RotateOrbsMod , only : RotateOrbs
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: i,j,error,HFConn
        CHARACTER(len=*), PARAMETER :: this_routine='FciMCPar'
        TYPE(HElement) :: Hamii
        LOGICAL :: TIncrement,tWritePopsFound,tSoftExitFound,tSingBiasChange
        REAL(4) :: s,etime,tstart(2),tend(2)

        TDebug=.false.  !Set debugging flag

!OpenMPI does not currently support MPI_Comm_set_errhandler - a bug in its F90 interface code.
!Ask Nick McLaren if we need to change the err handler - he has a fix/bypass.
!        CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_RETURN,error)
!        CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD,MPI_ERRORS_ARE_FATAL,error)
        
        CALL SetupParameters()
        CALL InitFCIMCCalcPar()
        CALL WriteFciMCStatsHeader()

        CALL WriteFCIMCStats()

!Start MC simulation...
        TIncrement=.true.   !If TIncrement is true, it means that when it comes out of the loop, it wants to subtract 1 from the Iteration count to get the true number of iterations
        Iter=1
        do while(Iter.le.NMCyc)   !Iter=1,NMCyc
!Main iteration loop...
!            WRITE(6,*) 'Iter',Iter

            IF(TBalanceNodes) THEN 
                CALL BalanceWalkersonProcs()      !This routine is call periodically when the nodes need to be balanced. However, this will not be called with direct annihilation.
            ENDIF
            
            s=etime(tstart)
            IF(tCCMC) THEN
                CALL PerformCCMCCycPar()
            ELSEIF(tMultipleDetsSpawn) THEN
                CALL MultipleConnFCIMCycPar()
            ELSEIF(tCleanRun) THEN
                CALL PerformCleanFCIMCycPar()
            ELSE
                CALL PerformFCIMCycPar()
            ENDIF
            s=etime(tend)
            IterTime=IterTime+(tend(1)-tstart(1))

            IF(mod(Iter,StepsSft).eq.0) THEN
!This will communicate between all nodes, find the new shift (and other parameters) and broadcast them to the other nodes.
                CALL CalcNewShift()

                IF((tTruncCAS.or.tTruncSpace).and.(Iter.gt.iFullSpaceIter).and.(iFullSpaceIter.ne.0)) THEN
!Test if we want to expand to the full space if an EXPANDSPACE variable has been set
                    IF(tHistSpawn.or.tCalcFCIMCPsi) THEN
                        IF(iProcIndex.eq.0) WRITE(6,*) "Unable to expand space since histgramming the wavefunction..."
                    ELSE
                        ICILevel=0
                        tTruncSpace=.false.
                        tTruncCAS=.false.
                        IF(iProcIndex.eq.0) THEN
                            WRITE(6,*) "Expanding to the full space on iteration ",Iter
                        ENDIF
                    ENDIF
                ENDIF

!This routine will check for a CHANGEVARS file and change the parameters of the calculation accordingly.
                CALL ChangeVars(tSingBiasChange,tSoftExitFound,tWritePopsFound)
                IF(tSoftExitFound) THEN
                    TIncrement=.false.
                    EXIT
                ENDIF
                IF(tWritePopsFound) THEN
!We have explicitly asked to write out the POPSFILE from the CHANGEVARS file.
                    IF(tRotoAnnihil) THEN
                        CALL WriteToPopsFileParOneArr()
                    ELSE
                        CALL WriteToPopsFilePar()
                    ENDIF
                ENDIF
                IF(tSingBiasChange) THEN
                    CALL GetSymExcitCount(HFExcit%ExcitData,HFConn)
                    CALL CalcApproxpDoubles(HFConn)
                ENDIF

                
            ENDIF

            IF(TPopsFile.and.(mod(Iter,iWritePopsEvery).eq.0)) THEN
!This will write out the POPSFILE if wanted
                IF(tRotoAnnihil) THEN
                    CALL WriteToPopsFileParOneArr()
                ELSE
                    CALL WriteToPopsFilePar()
                ENDIF
            ENDIF
!            IF(TAutoCorr) CALL WriteHistogrammedDets()

            IF(tHistSpawn.and.(mod(Iter,iWriteHistEvery).eq.0)) THEN
                CALL WriteHistogram()
            ENDIF

            Iter=Iter+1
!End of MC cycle
        enddo

        IF(TIncrement) Iter=Iter-1     !Reduce the iteration count for the POPSFILE since it is incremented upon leaving the loop (if done naturally)
        IF(TPopsFile) THEN
            IF(tRotoAnnihil) THEN
                CALL WriteToPopsFileParOneArr()
            ELSE
                CALL WriteToPopsFilePar()
            ENDIF
        ENDIF
        IF(tCalcFCIMCPsi) THEN
!This routine will actually only print the matrix if tPrintFCIMCPsi is on
            CALL PrintFCIMCPsi()

            IF(tFindCINatOrbs) THEN
!This routine takes the wavefunction Psi, calculates the one electron density matrix, and rotates the HF orbitals to produce a new ROFCIDUMP file.
                CALL RotateOrbs() 
                CALL MPI_Barrier(MPI_COMM_WORLD,error)
            ENDIF
        ENDIF

        IF(tErrorBlocking) CALL FinaliseBlocking(Iter)
        IF(tShiftBlocking) CALL FinaliseShiftBlocking(Iter)

        IF(tHistSpawn) CALL WriteHistogram()

        Weight=HDElement(0.D0)
        Energyxw=HDElement(ProjectionE)

        IF(tConstructNos) CALL NormandDiagOneRDM()

        IF(tHistEnergies) CALL WriteHistogramEnergies()

!If we are calculating a guiding function, write it out here.
        IF(tFindGuide) CALL WriteGuidingFunc()
        IF(tUseGuide) CALL WriteFinalGuidingFunc()

!If we are writing out the dominant determinats, do this here.
        IF(tPrintDominant) CALL PrintDominantDets()

        IF(tPrintTriConnections) CALL PrintTriConnHist
        IF(tHistTriConHEls) CALL PrintTriConnHElHist

!        IF(TAutoCorr) CALL CalcAutoCorr()

!Deallocate memory
        CALL DeallocFCIMCMemPar()

        IF(iProcIndex.eq.Root) THEN
            CLOSE(15)
!            IF(TAutoCorr) CLOSE(44)
        ENDIF
        IF(TDebug) CLOSE(11)

        RETURN

    END SUBROUTINE FciMCPar

    
!This is the heart of FCIMC, where the MC Cycles are performed. However, this version is clean and does not have unnecessary tests for experimental options.
    SUBROUTINE PerformCleanFCIMCycPar()
        INTEGER :: VecSlot,i,j,k,l,CopySign,iPartBloom
        INTEGER :: nJ(NEl),ierr,IC,Child,DetCurr(NEl),iLutnJ(0:NIfD)
        REAL*8 :: Prob,rat,HDiagCurr
        INTEGER :: iDie,WalkExcitLevel,Proc
        INTEGER :: ExcitLevel,TotWalkersNew,iGetExcitLevel_2,Ex(2,2),WSign,p,Scratch1(2,nSymLabels),Scratch2(2,nSymLabels)
        LOGICAL :: tParity,tFilled

        IF(TDebug.and.(mod(Iter,10).eq.0)) THEN
            WRITE(11,*) Iter,TotWalkers,NoatHF,NoatDoubs,MaxIndex,TotParts
            CALL FLUSH(11)
        ENDIF
        
        CALL set_timer(Walker_Time,30)
        
!VecSlot indicates the next free position in NewDets
        VecSlot=1
!Reset number at HF and doubles
        NoatHF=0
        NoatDoubs=0
        iPartBloom=0
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
        ValidSpawnedList(:)=InitialSpawnedSlots(:)

        do j=1,TotWalkers
!j runs through all determinants in the main list. TotWalkers is the number of unique DETERMINANTS in the main array of this processor, not walkers, despite the name.
!the sign indicates the sum of the signs of the walkers on the determinant, and hence j loops over determinants, not particles.
!            WRITE(6,*) Iter,j,TotWalkers
!            CALL FLUSH(6)

!First, decode the bit-string representation of the determinant the walker is on, into a string of naturally-ordered integers
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j),NEl,NIfD)

!Also, we want to find out the excitation level of the determinant - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.
            IF(tTruncSpace) THEN
!We need to know the exact excitation level for truncated calculations.
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,NEl)
            ELSE
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,2)
            ENDIF
!HDiags are stored.
            HDiagCurr=CurrentH(j)

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
            CALL SumEContrib(DetCurr,WalkExcitLevel,CurrentSign(j),CurrentDets(:,j),HDiagCurr,1.D0)

            tFilled=.false.     !This is for regenerating excitations from the same determinant multiple times. There will be a time saving if we can store the excitation generators temporarily.

            do p=1,abs(CurrentSign(j))
!Here, we spawn each particle on the determinant in a seperate attempt.
!we are simply looping over all the particles on the determinant
!This will only be a help if most determinants are multiply occupied.

                IF(tHPHF) THEN
                    CALL GenRandHPHFExcit2Scratch(DetCurr,CurrentDets(:,j),nJ,iLutnJ,pDoubles,exFlag,Prob,Scratch1,Scratch2,tFilled,tGenMatHEl)
                ELSE
                    CALL GenRandSymExcitScratchNU(DetCurr,CurrentDets(:,j),nJ,pDoubles,IC,Ex,tParity,exFlag,Prob,Scratch1,Scratch2,tFilled)
                ENDIF
!                WRITE(6,"(A,7I5)") "Generated: ",nJ(:)


!Calculate number of children to spawn
                IF(TTruncSpace) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.
                    IF(CheckAllowedTruncSpawn(WalkExcitLevel,nJ,iLutnJ,IC)) THEN
!The excitation is allowed - it is below the ICILevel cutoff.
                        Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,1,.false.)
                    ELSE
                        Child=0
                    ENDIF
                ELSE
!SD Space is not truncated - allow attempted spawn as usual
                    Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,1,.false.)
                ENDIF
                
                IF(Child.ne.0) THEN
!We want to spawn a child - find its information to store

!                    WRITE(6,*) "Spawning particle to:",nJ(:)
!                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,NEl)
!                    WRITE(6,*) "Excitlevel:", ExcitLevel
                    NoBorn=NoBorn+abs(Child)     !Update counter about particle birth
                    IF(IC.eq.1) THEN
                        SpawnFromSing=SpawnFromSing+abs(Child)
                    ENDIF

                    IF(abs(Child).gt.25) THEN
!If more than 25 particles are created in one go, then log this fact and print out later that this has happened.
                        IF(abs(Child).gt.abs(iPartBloom)) THEN
                            IF(IC.eq.1) THEN
                                iPartBloom=-abs(Child)
                            ELSE
                                iPartBloom=abs(Child)
                            ENDIF
                        ENDIF
!                        WRITE(6,"(A,I10,A)") "LARGE PARTICLE BLOOM - ",Child," particles created in one attempt."
!                        WRITE(6,"(A,I5)") "Excitation: ",IC
!                        WRITE(6,"(A,G25.10)") "PROB IS: ",Prob
!                        CALL FLUSH(6)
                    ENDIF

!In direct annihilation, we spawn particles into a seperate array, but we do not store them contiguously in the SpawnedParts/SpawnedSign arrays.
!The processor that the newly-spawned particle is going to be sent to has to be determined, and then it will get put into the the appropriate element determined by ValidSpawnedList.

                    Proc=DetermineDetProc(iLutnJ)   !This wants to return a value between 0 -> nProcessors-1
!                    WRITE(6,*) iLutnJ(:),Proc,ValidSpawnedList(Proc),Child,TotWalkers
!                    CALL FLUSH(6)
                    SpawnedParts(:,ValidSpawnedList(Proc))=iLutnJ(:)
                    SpawnedSign(ValidSpawnedList(Proc))=Child
                    ValidSpawnedList(Proc)=ValidSpawnedList(Proc)+1

                    Acceptances=Acceptances+ABS(Child)      !Sum the number of created children to use in acceptance ratio
                
                ENDIF   !End if child created

            enddo   !End of cycling over mulitple particles on same determinant.

!We now have to decide whether the parent particle (j) wants to die or not...
!we can have multiple particles on the same determinant - these can be stochastically killed at the same time.

            iDie=AttemptDiePar(DetCurr,HDiagCurr,WalkExcitLevel,CurrentSign(j))
            NoDied=NoDied+iDie          !Update death counter

!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births
!We slot the particles back into the same array and position VecSlot if the particle survives. If it dies, then j increases, moving onto the next
!entry, but VecSlot remains where it is, meaning that j should never be less that VecSlot

            IF(CurrentSign(j).le.0) THEN
                CopySign=CurrentSign(j)+iDie    !Copy sign is the total number of particles x sign that we want to copy accross.
                IF(CopySign.gt.0) THEN
!If we are copying to the main array, we have to ensure that we maintain sign-coherence in the array. Therefore, if we are spawning anti-particles,
!it wants to go in the spawning array, rather than the main array, so it has a chance to annihilate. However, since anti-particles should not be created
!in normal circumstances, we will remove this possibility.
                    WRITE(6,*) "***",Iter,CopySign,HDiagCurr,iDie,DetCurr(:)
                    CALL Stop_All("PerformFCIMCyc","Creating anti-particles")
                ENDIF
            ELSE
                CopySign=CurrentSign(j)-iDie
                IF(CopySign.lt.0) THEN
                    WRITE(6,*) "***",Iter,CopySign,HDiagCurr,iDie,DetCurr(:)
                    CALL Stop_All("PerformFCIMCyc","Creating anti-particles")
                ENDIF
            ENDIF

!We only have a single array, therefore surviving particles are simply transferred back into the original array.
!We should not get to the case where we want to overwrite particles that we haven't even got to yet.
            IF(CopySign.ne.0) THEN

                CurrentDets(:,VecSlot)=CurrentDets(:,j)
                CurrentSign(VecSlot)=CopySign
                CurrentH(VecSlot)=HDiagCurr
                VecSlot=VecSlot+1

            ENDIF   !To kill if

!Finish cycling over determinants
        enddo


!***Birth/death processes finished. Tidy up and then annihilate.

!SumWalkersCyc calculates the total number of walkers over an update cycle on each process
        SumWalkersCyc=SumWalkersCyc+(INT(TotParts,i2))
!        WRITE(6,*) "Born, Die: ",NoBorn, NoDied

!Since VecSlot holds the next vacant slot in the main determinant array, TotWalkers will be one less than this.
        TotWalkersNew=VecSlot-1

!Output if there has been a particle bloom this iteration. A negative number indicates that particles were created from a single excitation.
!This is now ONLY written out for blooms on the head node.
        IF((iPartBloom.ne.0).and.(iProcIndex.eq.0)) THEN
            WRITE(6,"(A,I10,A)") "LARGE PARTICLE BLOOMS in iteration ",Iter
            IF(iPartBloom.gt.0) THEN
                WRITE(6,"(A,I10,A)") "A max of ",abs(iPartBloom)," particles created in one attempt from double excit."
            ELSE
                WRITE(6,"(A,I10,A)") "A max of ",abs(iPartBloom)," particles created in one attempt from single excit."
            ENDIF
        ENDIF


        rat=(TotWalkersNew+0.D0)/(MaxWalkersPart+0.D0)
        IF(rat.gt.0.95) THEN
            WRITE(6,*) "*WARNING* - Number of particles/determinants has increased to over 95% of MaxWalkersPart"
            CALL FLUSH(6)
        ENDIF

!Need to test whether any of the sublists in the spawning array are getting to the end of their allotted space.
        IF(nProcessors.gt.1) THEN
            do i=0,nProcessors-1
                rat=(ValidSpawnedList(i)-InitialSpawnedSlots(i))/(InitialSpawnedSlots(1)+0.D0)
                IF(rat.gt.0.95) THEN
                    WRITE(6,*) "*WARNING* - Highest processor spawned particles has reached over 95% of MaxSpawned"
                    CALL FLUSH(6)
                ENDIF
            enddo
        ELSE
            rat=(ValidSpawnedList(0)+0.D0)/(MaxSpawned+0.D0)
            IF(rat.gt.0.9) THEN
                WRITE(6,*) "*WARNING* - Number of spawned particles has reached over 90% of MaxSpawned"
                CALL FLUSH(6)
            ENDIF
        ENDIF
        
        CALL halt_timer(Walker_Time)
        CALL set_timer(Annihil_Time,30)
!        CALL MPI_Barrier(MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get into annihilation"
!        CALL FLUSH(6)

!This is the direct annihilation algorithm. The newly spawned walkers should be in a seperate array (SpawnedParts) and the other list should be ordered.
        CALL DirectAnnihilation(TotWalkersNew)

        CALL halt_timer(Annihil_Time)
        
    END SUBROUTINE PerformCleanFCIMCycPar

!This is the heart of FCIMC, where the MC Cycles are performed.
    SUBROUTINE PerformFCIMCycPar()
!        use HPHFRandExcitMod , only : TestGenRandHPHFExcit 
        USE Determinants , only : GetHElement3
        USE FciMCLoggingMOD , only : FindTriConnections,TrackSpawnAttempts,FindSpinCoupHEl
        INTEGER :: MinorVecSlot,VecSlot,i,j,k,l,MinorValidSpawned,ValidSpawned,CopySign,ParticleWeight,Loop,iPartBloom
        INTEGER :: nJ(NEl),ierr,IC,Child,iCount,DetCurr(NEl),iLutnJ(0:NIfD),NoMinorWalkersNew
        REAL*8 :: Prob,rat,HDiag,HDiagCurr
        INTEGER :: iDie,WalkExcitLevel,Proc             !Indicated whether a particle should self-destruct on DetCurr
        INTEGER :: ExcitLevel,TotWalkersNew,iGetExcitLevel_2,error,length,temp,Ex(2,2),WSign,p,Scratch1(2,nSymLabels),Scratch2(2,nSymLabels),FDetSym,FDetSpin
        LOGICAL :: tParity,DetBitEQ,tMainArr,tFilled,tCheckStarGenDet,tStarDet,tMinorDetList,tAnnihilateMinorTemp,tAnnihilateMinor,TestClosedShellDet
        INTEGER(KIND=i2) :: HashTemp
        TYPE(HElement) :: HDiagTemp
        CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: message
        REAL :: Gap

        IF(TDebug.and.(mod(Iter,10).eq.0)) THEN
            WRITE(11,*) Iter,TotWalkers,NoatHF,NoatDoubs,MaxIndex,TotParts
            CALL FLUSH(11)
        ENDIF
        
!        IF(tRotoAnnihil) THEN
!            CALL CheckOrdering(CurrentDets(:,1:TotWalkers),CurrentSign(1:TotWalkers),TotWalkers,.true.)
!        ENDIF

        CALL set_timer(Walker_Time,30)
        
!VecSlot indicates the next free position in NewDets
        VecSlot=1
!Reset number at HF and doubles
        NoatHF=0
        NoatDoubs=0
!        DetsNorm=0.D0
        iPartBloom=0
        ValidSpawned=1  !This is for rotoannihilation - this is the number of spawned particles (well, one more than this.)
        IF(tDirectAnnihil) THEN
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
            ValidSpawnedList(:)=InitialSpawnedSlots(:)
        ENDIF
        MinorValidSpawned=1
        tMinorDetList=.false.
        IF(.not.tStarOrbs) tStarDet=.false.
        ParticleWeight=1    !This will always be the same unless we are 'spawning as determinants'

        IF(tHistSpawn.or.tCalcFCIMCPsi.and.(Iter.ge.NHistEquilSteps)) THEN
            IF(Iter.eq.NHistEquilSteps) THEN
                IF(iProcIndex.eq.Root) WRITE(6,*) 'The iteration is equal to HISTEQUILSTEPS.  Beginning to histogram.'
            ENDIF
            HistMinInd(1:NEl)=FCIDetIndex(1:NEl)    !This is for the binary search when histogramming
        ENDIF
        
        do j=1,TotWalkers
!j runs through all current walkers
!If we are rotoannihilating/direct annihilating, the sign indicates the sum of the signs on the determinant, and hence j loops over determinants, not particles.
!            WRITE(6,*) Iter,j,TotWalkers
!            CALL FLUSH(6)

!First, decode the bit-string representation of the determinant the walker is on, into a string of naturally-ordered integers
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j),NEl,NIfD)
!            IF((Iter.gt.100)) THEN!.and.(.not.DetBitEQ(CurrentDets(:,j),iLutHF,NIfD))) THEN
!This will test the excitation generator for HPHF wavefunctions
!                IF(.not.(TestClosedShellDet(CurrentDets(:,j),NIfD))) THEN
!                    CALL TestGenRandHPHFExcit(DetCurr,NMCyc,pDoubles)
!                    STOP
!                ENDIF
!            ENDIF
!            FDetSym=0
!            FDetSpin=0
!            do i=1,NEl
!               FDetSym=IEOR(FDetSym,INT(G1(DetCurr(i))%Sym%S,4))
!               FDetSpin=FDetSpin+G1(DetCurr(i))%Ms
!            enddo
!            IF(FDetSym.ne.4) THEN
!                WRITE(6,*) "Symmetry of determinant is: ",FDetSym
!            ENDIF
!            IF(FDetSpin.ne.0) THEN
!                WRITE(6,*) "Spin of determinant is: ",FDetSpin
!            ENDIF

!Also, we want to find out the excitation level - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.
            IF(tTruncSpace.or.tHighExcitsSing.or.tHistSpawn.or.tCalcFCIMCPsi.or.tPrintSpinCoupHEl) THEN
!We need to know the exact excitation level for truncated calculations.
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,NEl)
                IF((WalkExcitLevel.eq.2).and.tPrintSpinCoupHEl) CALL FindSpinCoupHEl(iLutHF,CurrentDets(:,j),NIfD,NEl)
            ELSE
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,2)
            ENDIF
            IF(tRegenDiagHEls) THEN
!We are not storing the diagonal hamiltonian elements for each particle. Therefore, we need to regenerate them.
!Need to find H-element!
                IF(DetBitEQ(CurrentDets(0:NIfD,j),iLutHF,NIfD).and.(.not.(tHub.and.tReal))) THEN
!We know we are at HF - HDiag=0
                    HDiagCurr=0.D0
                ELSE
                    IF(tHPHF) THEN
                        CALL HPHFGetDiagHElement(DetCurr,CurrentDets(:,j),HDiagTemp)
                    ELSE
                        HDiagTemp=GetHElement2(DetCurr,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                    ENDIF
                    HDiagCurr=(REAL(HDiagTemp%v,r2))-Hii
                ENDIF
            ELSE
!HDiags are stored.
                HDiagCurr=CurrentH(j)
            ENDIF
            IF(tFindGroundDet) THEN
                IF(HDiagCurr.lt.0.D0) THEN
!We have found a determinant lower in energy that the "root" determinant.
!This should not happen in a HF basis, but can happen if the orbitals have been rotated.
                    CALL ChangeRefDet(HDiagCurr,DetCurr,CurrentDets(:,j))
                ENDIF
            ENDIF
            IF(tStarOrbs) THEN
!This routine will check to see if any of the orbitals in the determinant are in the orbitals which are only to be attached to HF in a 'star'
                CALL CheckStarOrbs(DetCurr,tStarDet)
            ENDIF

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
            CALL SumEContrib(DetCurr,WalkExcitLevel,CurrentSign(j),CurrentDets(:,j),HDiagCurr,1.D0)

!            IF(TResumFCIMC) CALL ResumGraphPar(DetCurr,CurrentSign(j),VecSlot,j)

            tFilled=.false.     !This is for regenerating excitations from the same determinant multiple times. There will be a time saving if we can store the excitation generators temporarily.
            IF(tSpawnAsDet) THEN
!Here, we spawn all particles on the determinant in one go, by multiplying the probability of spawning by the number of particles on the determinant.
                Loop=1
                ParticleWeight=abs(CurrentSign(j))
            ELSE
!Here, we spawn each particle on the determinant in a seperate attempt.
                Loop=abs(CurrentSign(j))
            ENDIF

            do p=1,Loop
!If rotoannihilating, we are simply looping over all the particles on the determinant

                IF(.not.tImportanceSample) THEN
                    IF(.not.TRegenExcitgens) THEN
!Setup excit generators for this determinant
                        CALL SetupExitgenPar(DetCurr,CurrentExcits(j))
                        CALL GenRandSymExcitIt4(DetCurr,CurrentExcits(j)%PointToExcit,nJ,0,IC,0,Prob,iCount,Ex,tParity)
                    ELSE
                        IF(tStarDet) THEN
                            IF((.not.tNoReturnStarDets).and.(.not.tAllSpawnStarDets)) THEN
!We are at a determinant with high-lying 'star' orbitals - we therefore only want to be able to generate the HF determinant
                                nJ(:)=HFDet(:)
                                IC=WalkExcitLevel
                                Ex(1,1)=WalkExcitLevel
                                CALL GetExcitation(DetCurr,nJ,NEl,Ex,tParity)
                                Prob=1.D0   !This is the only allowed excitation - back to HF
                            ENDIF
                        ELSEIF(tHighExcitsSing) THEN
!Here, we only allow single excitations from iHighExcitsSing if the excitation is going further away from the current excitation.
                            IF(WalkExcitLevel.ge.(iHighExcitsSing+2)) THEN
                                CALL GenRandSymExcitNU(DetCurr,CurrentDets(:,j),nJ,0.D0,IC,Ex,tParity,exFlag,Prob)
                            ELSE
                                CALL GenRandSymExcitNU(DetCurr,CurrentDets(:,j),nJ,pDoubles,IC,Ex,tParity,exFlag,Prob)
                            ENDIF

                        ELSE
                            IF(tNonUniRandExcits) THEN
!This will only be a help if most determinants are multiply occupied.
                                IF(tHPHF) THEN
!                                    CALL GenRandHPHFExcit(DetCurr,CurrentDets(:,j),nJ,iLutnJ,pDoubles,exFlag,Prob)
                                    CALL GenRandHPHFExcit2Scratch(DetCurr,CurrentDets(:,j),nJ,iLutnJ,pDoubles,exFlag,Prob,Scratch1,Scratch2,tFilled,tGenMatHEl)
                                ELSE
                                    CALL GenRandSymExcitScratchNU(DetCurr,CurrentDets(:,j),nJ,pDoubles,IC,Ex,tParity,exFlag,Prob,Scratch1,Scratch2,tFilled)
                                ENDIF
                            ELSE
                                CALL GetPartRandExcitPar(DetCurr,CurrentDets(:,j),nJ,IC,0,Prob,iCount,WalkExcitLevel,Ex,tParity)
                            ENDIF
                        ENDIF
                    ENDIF
                ENDIF

!If we want to look at determinants connected in loops of 3, need to select a third determinant.  Currently have DetCurr exciting to nJ.  Now need DetCurr
!to excite to a different nk. - then add the product of the three connecting Helements to either the sign coherent or sign incoherent loops.
                IF(tPrintTriConnections) CALL FindTriConnections(DetCurr,CurrentDets(:,j),iLutHF,nJ,IC,Ex,pDoubles,tFilled,tParity,Scratch1,Scratch2,exflag)

!Calculate number of children to spawn
                IF(TTruncSpace.or.tTruncCAS.or.tListDets.or.tPartFreezeCore) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.
                    IF(tImportanceSample) CALL Stop_All("PerformFCIMCyc","Truncated calculations not yet working with importance sampling")

                    IF(CheckAllowedTruncSpawn(WalkExcitLevel,nJ,iLutnJ,IC)) THEN
                        Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                    ELSE
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                        Child=0
                    ENDIF
                ELSE
!SD Space is not truncated - allow attempted spawn as usual
                    IF(tStarOrbs) THEN

                        IF(tStarDet) THEN
                            IF(tNoReturnStarDets.or.tAllSpawnStarDets) THEN
!We do not allow a return spawn back to HF, or any other determinant (allSpawnStarDets)
                                Child=0
                            ELSE
!Here, we are at a high-energy det and have generated HF which we will try to spawn to.
                                Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                            ENDIF

                        ELSEIF((WalkExcitLevel.eq.0).or.tAllSpawnStarDets) THEN
!We are at HF - all determinants allowed. No need to check generated excitations
                            Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)

                        ELSE
!We need to check whether the excitation generated is in the allowed or disallowed space. High-lying orbitals cannot be generated from non-HF determinants.
!If tStarDet is true, then we know we are already at one of these determinants, and have alraedy generated HF. If we are at HF, then we are allowed to generate these determinants.
                            CALL CheckStarOrbs(nJ,tCheckStarGenDet)
                            IF(tCheckStarGenDet) THEN
!We have generated a 'star' determinant. We are not at the HF determinant, so disallow it and Child=0
                                Child=0
                            ELSE
!We have not generated a 'star' determinant. Allow as normal.
                                Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                            ENDIF

                        ENDIF

                    ELSEIF(tHighExcitsSing) THEN
!This option only allows singles connections between determinants above iHighExcitsSing threshold.
                        
                        IF((WalkExcitLevel.eq.(iHighExcitsSing+1)).and.(IC.eq.2)) THEN
!Only allow doubles to the cutoff-1
                            ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,iHighExcitsSing-1)
                            IF(ExcitLevel.eq.(iHighExcitsSing-1)) THEN
                                Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                            ELSE
                                Child=0
                            ENDIF
                        ELSEIF((WalkExcitLevel.eq.iHighExcitsSing).and.(IC.eq.2)) THEN
!Only allow doubles to the cutoff, or less
                            ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,iHighExcitsSing-1)
                            IF(ExcitLevel.le.iHighExcitsSing) THEN
                                Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                            ELSE
                                Child=0
                            ENDIF
                        ELSE
                            Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                        ENDIF

                    ELSE

                        Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                    ENDIF

                ENDIF

! Want to put a wee routine in here that monitors the number of accepted vs the number of not accepted attempts at spawns, and the H elements that
! are involved in each.
                IF(tPrintHElAccept) CALL TrackSpawnAttempts(Child,DetCurr,j,nJ,iLutnJ,IC,Ex,tParity)
                
                IF(Child.ne.0) THEN
!We want to spawn a child - find its information to store

!                    WRITE(6,*) "Spawning particle to:",nJ(:)
!                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,NEl)
!                    WRITE(6,*) "Excitlevel:", ExcitLevel
                    NoBorn=NoBorn+abs(Child)     !Update counter about particle birth
                    IF(IC.eq.1) THEN
                        SpawnFromSing=SpawnFromSing+abs(Child)
                    ENDIF

                    IF(abs(Child).gt.25) THEN
!If more than 25 particles are created in one go, then log this fact and print out later that this has happened.
                        IF(abs(Child).gt.abs(iPartBloom)) THEN
                            IF(IC.eq.1) THEN
                                iPartBloom=-abs(Child)
                            ELSE
                                iPartBloom=abs(Child)
                            ENDIF
                        ENDIF
!                        WRITE(6,"(A,I10,A)") "LARGE PARTICLE BLOOM - ",Child," particles created in one attempt."
!                        WRITE(6,"(A,I5)") "Excitation: ",IC
!                        WRITE(6,"(A,G25.10)") "PROB IS: ",Prob
!                        CALL FLUSH(6)
                    ENDIF

                    IF(tRotoAnnihil) THEN
!In the RotoAnnihilation implimentation, we spawn particles into a seperate array - SpawnedParts and SpawnedSign. 
!The excitation level and diagonal matrix element are also found out after the annihilation.
!Cannot use old excitation generators with rotoannihilation.

!In rotoannihilation, we can specify multiple particles on the same entry. 
                        IF(tMinorDetList) THEN
!We want to add the determinants to a seperate list, since they are spawning back from "insignificant" determinants.
                            MinorSpawnDets(0:NIfD,MinorValidSpawned)=iLutnJ(0:NIfD)
                            MinorSpawnParent(0:NIfD,MinorValidSpawned)=CurrentDets(0:NIfD,j) !This is DetCurr in bit form
                            MinorSpawnSign(MinorValidSpawned)=Child
!                            CALL DecodeBitDet(TempDet,iLutnJ(:),NEl,NIfD)
                            HashArray(MinorValidSpawned)=CreateHash(nJ)
                            MinorValidSpawned=MinorValidSpawned+1
                            ! MinorValidSpawned is the number spawned on the minor determinants.
                        ELSE
                            SpawnedParts(:,ValidSpawned)=iLutnJ(:)
                            SpawnedSign(ValidSpawned)=Child
                            ValidSpawned=ValidSpawned+1     !Increase index of spawned particles
                        ENDIF

                    ELSEIF(tDirectAnnihil) THEN
!In direct annihilation, we spawn particles into a seperate array, but we do not store them contiguously in the SpawnedParts/SpawnedSign arrays.
!The processor that the newly-spawned particle is going to be sent to has to be determined, and then it will get put into the the appropriate element determined by ValidSpawnedList.

                        Proc=DetermineDetProc(iLutnJ)   !This wants to return a value between 0 -> nProcessors-1
!                        WRITE(6,*) iLutnJ(:),Proc,ValidSpawnedList(Proc),Child,TotWalkers
!                        CALL FLUSH(6)
                        SpawnedParts(:,ValidSpawnedList(Proc))=iLutnJ(:)
                        SpawnedSign(ValidSpawnedList(Proc))=Child
                        ValidSpawnedList(Proc)=ValidSpawnedList(Proc)+1

                    ELSE
!Calculate diagonal ham element

                        IF(Child.gt.0) THEN
!We have successfully created at least one positive child at nJ
                            WSign=1
                        ELSE
!We have successfully created at least one negative child at nJ
                            WSign=-1
                        ENDIF

                        IF(TMagnetize) THEN
                            CALL FindBitExcitLevel(iLutnJ,iLutHF,NIfD,ExcitLevel,2)
                            CALL FindDiagElwithB(HDiag,ExcitLevel,nJ,WSign)
                        ELSE
                            IF(.not.tRegenDiagHEls) THEN
                                IF(DetBitEQ(iLutnJ,iLutHF,NIfD)) THEN
!We know we are at HF - HDiag=0
                                    HDiag=0.D0
!                                        IF(tHub.and.tReal) THEN
!!Reference determinant is not HF
!                                            HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                                            HDiag=(REAL(HDiagTemp%v,r2))
!                                        ENDIF

                                ELSE
                                    IF(tHPHF) THEN
                                        CALL HPHFGetDiagHElement(nJ,iLutnJ,HDiagTemp)
                                    ELSE
                                        HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                                    ENDIF
                                    HDiag=(REAL(HDiagTemp%v,r2))-Hii
                                ENDIF
                            ENDIF
                        ENDIF

!                            IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                        IF(.not.TNoAnnihil) THEN
                            HashTemp=CreateHash(nJ)
                        ENDIF

                        do l=1,abs(Child)
!Copy across children - cannot copy excitation generators, as do not know them...unless it is HF (this might save a little time if implimented)
                            NewDets(:,VecSlot)=iLutnJ(:)
                            NewSign(VecSlot)=WSign
                            IF(.not.TRegenExcitgens) NewExcits(VecSlot)%PointToExcit=>null()
                            IF(.not.tRegenDiagHEls) NewH(VecSlot)=HDiag                     !Diagonal H-element-Hii
                            IF(.not.TNoAnnihil) THEN
!                                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                                Hash2Array(VecSlot)=HashTemp        !Hash put in Hash2Array - no need for pointer since always annihilating if storing hashes
                            ENDIF
                            VecSlot=VecSlot+1
                        enddo
                    
                    ENDIF   !Endif rotoannihil

                    Acceptances=Acceptances+ABS(Child)      !Sum the number of created children to use in acceptance ratio
                
                ENDIF   !End if child created

            enddo   !End of cycling over mulitple particles on same determinant.

!We now have to decide whether the parent particle (j) wants to self-destruct or not...
!For rotoannihilation, we can have multiple particles on the same determinant - these can be stochastically killed at the same time.

            iDie=AttemptDiePar(DetCurr,HDiagCurr,WalkExcitLevel,CurrentSign(j))
            NoDied=NoDied+iDie          !Update death counter

!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births
            IF(tRotoAnnihil.or.tDirectAnnihil) THEN
!We slot the particles back into the same array and position VecSlot if the particle survives. If it dies, then j increases, moving onto the next
!entry, but VecSlot remains where it is, meaning that j should never be less that VecSlot

                IF(CurrentSign(j).le.0) THEN
                    CopySign=CurrentSign(j)+iDie    !Copy sign is the total number of particles x sign that we want to copy accross.
                    IF(CopySign.le.0) THEN
!If we are copying to the main array, we have to ensure that we maintain sign-coherence in the array. Therefore, if we are spawning anti-particles,
!it wants to go in the spawning array, rather than the main array, so it has a chance to annihilate.
                        tMainArr=.true.
                    ELSE
                        tMainArr=.false.
                    ENDIF
                ELSE
                    CopySign=CurrentSign(j)-iDie
                    IF(CopySign.ge.0) THEN
                        tMainArr=.true.
                    ELSE
                        tMainArr=.false.
                    ENDIF
                ENDIF

                IF(tMainArr.and.(VecSlot.gt.j)) THEN
!We only have a single array, therefore surviving particles are simply transferred back into the original array.
!However, this can not happen if we want to overwrite particles that we haven't even got to yet.
!However, this should not happen under normal circumstances...
                    tMainArr=.false.
                    CALL Stop_All("PerformFCIMCyc","About to overwrite particles which haven't been reached yet. This should not happen in normal situations.")
                ENDIF

                IF(CopySign.ne.0) THEN

                    IF(tMainArr) THEN
!                        NewDets(:,VecSlot)=CurrentDets(:,j)
!                        NewSign(VecSlot)=CopySign
                        CurrentDets(:,VecSlot)=CurrentDets(:,j)
                        CurrentSign(VecSlot)=CopySign
                        IF(.not.tRegenDiagHEls) CurrentH(VecSlot)=HDiagCurr
                        VecSlot=VecSlot+1
                    ELSE
!                        CALL Stop_All("PerformFCIMCyc","Creating anti-particles")
!Here, we are creating anti-particles, and so to keep the sign-coherence of the main array assured, we transfer them to the spawning array.
!This should generally not happen under normal circumstances.
                        do p=1,abs(CopySign)
!In rotoannihilation, we still want to specify the determinants singly - this may change in the future...
                            SpawnedParts(:,ValidSpawned)=CurrentDets(:,j)
                            IF(CopySign.lt.0) THEN
                                SpawnedSign(ValidSpawned)=-1
                            ELSE
                                SpawnedSign(ValidSpawned)=1
                            ENDIF
                            ValidSpawned=ValidSpawned+1     !Increase index of spawned particles
                        enddo
                    ENDIF

                ENDIF

            ELSE    !Not rotoannihilation

                IF(iDie.le.0) THEN
!This indicates that the particle is spared and we may want to create more...copy them across to NewDets
!If iDie < 0, then we are creating the same particles multiple times. Copy accross (iDie+1) copies of particle
        
                    do l=1,abs(iDie)+1    !We need to copy accross one more, since we need to include the original spared particle
                        NewDets(:,VecSlot)=CurrentDets(:,j)
                        IF(CurrentSign(j).lt.0) THEN
                            NewSign(VecSlot)=-1
                        ELSE
                            NewSign(VecSlot)=1
                        ENDIF
!Copy excitation generator accross
!                        IF(ASSOCIATED(NewExcits(VecSlot)%PointToExcit)) THEN
!                            WRITE(6,*) "Problem is here",NewExcits(VecSlot)%IndexinExArr
!                        ENDIF
                        IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(j),NewExcits(VecSlot),.true.)
                        IF(.not.tRegenDiagHEls) NewH(VecSlot)=HDiagCurr
                        IF(.not.TNoAnnihil) Hash2Array(VecSlot)=HashArray(j)
!                        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) Hash2Array(VecSlot)=HashArray(j)
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
                        IF(CurrentSign(j).eq.1) THEN
!Copy accross new anti-particles
                            NewSign(VecSlot)=-1
                        ELSE
                            NewSign(VecSlot)=1
                        ENDIF
!Copy excitation generator accross
                        IF(l.eq.iDie-1) THEN
!Delete old excitation generator behind us (ie the final one is a move, rather than copy)
                            IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(j),NewExcits(VecSlot),.true.)
                        ELSE
                            IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(j),NewExcits(VecSlot),.false.)
                        ENDIF
                        IF(.not.tRegenDiagHEls) NewH(VecSlot)=HDiagCurr
                        IF((.not.TNoAnnihil).and.(.not.tRotoAnnihil)) Hash2Array(VecSlot)=HashArray(j)
!                        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) Hash2Array(VecSlot)=HashArray(j)
                        VecSlot=VecSlot+1
                    enddo

                ENDIF   !Roto Kill If
        
            ENDIF   !To kill if

!Finish cycling over walkers
        enddo

!        IF(MinorValidSpawned.gt.5) THEN
!            WRITE(6,*) 'MinorValidSpawned',MinorValidSpawned
!            WRITE(6,*) 'Det,Parent,Sign,Hash'
!            do l=1,MinorValidSpawned-1
!                WRITE(6,*) MinorSpawnDets(:,l),"**",MinorSpawnParent(:,l),"**",MinorSpawnSign(l),"**",HashArray(l) 
!            enddo
!            stop
!        ENDIF

        MinorValidSpawned=MinorValidSpawned-1

        MinorVecSlot=1
        IF(tMinorDetsStar) THEN
!Run through list of "insignificant" determinants.  Attempt spawning back to parents, and then attempt to die.
            do i=1,NoMinorWalkers
                ! Usually run over all determinants, nJ, and attempt to spawn on these, but we can only spawn back on parent, so only run over this 
                ! with probability 1.D0.
                ! nJ is the determinant we are attempting to spawn on in full (i.e. the parent in full form).
!                CALL DecodeBitDet(nJ,MinorStarParent(0:NIfD,i),NEl,NIfD)

                ! CHECK THIS
                Child=AttemptCreateParBack(MinorStarDets(0:NIfD,i),MinorStarParent(0:NIfD,i),MinorStarSign(i),MinorStarHij(i),abs(MinorStarSign(i)),tMinorDetList)
                ! This will give an integer which is the number of walkers (w sign) being spawned back to the parent (allowed) determinant.

                IF(tMinorDetList) THEN
                    WRITE(6,*) 'Spawing from',MinorStarDets(0:NIfD,i)
                    WRITE(6,*) 'attempting to spawn to',MinorStarParent(0:NIfD,i)
                    CALL Stop_All('PerformFCIMCycPar','ERROR. Attempting to spawn between minor determinants.')
                ENDIF
                ! Check that tMinorDetList is always false (in future change it so that it doesn't search).
                ! If tMinorDetList is true, then nJ is not in the allowed list, and the parent list must be wrong.

!If child.ne.0, then add it, but add it to the normal spawning list (not MinorSpawnDets)
                IF(Child.ne.0) THEN
!                    SpawnedParts(:,ValidSpawned)=iLutnJ(:)
                    SpawnedParts(0:NIfD,ValidSpawned)=MinorStarParent(0:NIfD,i)
                    SpawnedSign(ValidSpawned)=Child
                    ValidSpawned=ValidSpawned+1     !Increase index of spawned particles
                ENDIF
           
!Attempt Die for particles on "insignificant" dets

                ! DetCurr is the current determinant in expanded form, MinorStarDets(:,i) is the bit form.
                CALL DecodeBitDet(DetCurr,MinorStarDets(0:NIfD,i),NEl,NIfD)

                iDie=AttemptDiePar(DetCurr,REAL(MinorStarHii(i)%v,r2),0,abs(MinorStarSign(i)))
                ! Take the ith minor determinant and attempt to die.
                ! iDie gives the number of particles to die on that determinant.

                NoDied=NoDied+iDie
                IF((iDie.lt.0).or.(iDie.gt.ABS(MinorStarSign(i)))) CALL Stop_All('PerformFCIMCycPar','Error in attempting to die, trying to kill or create too many walkers.')
                ! Killing too many, or attempting to create particles.  If this is happening often, need to allow for this in the code.

                IF(MinorStarSign(i).le.0) THEN
                    CopySign=MinorStarSign(i)+iDie    !Copy sign is the total number of particles x sign that we want to copy accross.
                ELSE
                    CopySign=MinorStarSign(i)-iDie
                ENDIF

                IF(CopySign.ne.0) THEN
                    ! If copysign is 0, all particles on that determinant are 0 and this determinant is removed from the list.
                    MinorStarDets(:,MinorVecSlot)=MinorStarDets(:,i)
                    MinorStarSign(MinorVecSlot)=CopySign
                    MinorStarParent(:,MinorVecSlot)=MinorStarParent(:,i)
                    MinorStarHii(MinorVecSlot)=MinorStarHii(i)
                    MinorStarHij(MinorVecSlot)=MinorStarHij(i)
                    MinorVecSlot=MinorVecSlot+1
                ENDIF
            enddo
            NoMinorWalkersNew=MinorVecSlot-1
            ! This is the number of walkers on the minor determinants that have survived, i.e. number in MinorStarDets.  Update this again after annihilation.
            ! This is the number of walkers in the spawned arrays on the minor determinants.
            
            ! Now need to annihilate amongst the minor determinants.
            ! The spawned list may contain the same determinant twice, if the walkers have different parents.
            ! First need to annihilate within the spawned walkers, then rotate around the MinorStarDets and annihilate with these.
            ! Walkers that survive all this are put into MinorStarDets maintaining order (again, determinants may be specified more than once, but these should have different
            ! parents).
            
!            WRITE(6,*) 'MinorValidSpawned before rotoannihilation',MinorValidSpawned
            IF(MinorValidSpawned.gt.0) tAnnihilateMinorTemp=.true.
            CALL MPI_AllReduce(tAnnihilateMinorTemp,tAnnihilateMinor,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,error)

            IF(tAnnihilateMinor) CALL RotoAnnihilateMinorSpawned(MinorValidSpawned,NoMinorWalkersNew)

        ENDIF

                    
        
!SumWalkersCyc calculates the total number of walkers over an update cycle on each process
        SumWalkersCyc=SumWalkersCyc+(INT(TotParts,i2))
!        WRITE(6,*) "Born, Die: ",NoBorn, NoDied

!Since VecSlot holds the next vacant slot in the array, TotWalkers will be one less than this.
        TotWalkersNew=VecSlot-1

!Output if there has been a particle bloom this iteration. A negative number indicates that particles were created from a single excitation.
        IF((iPartBloom.ne.0).and.(iProcIndex.eq.0)) THEN
            WRITE(6,"(A,I10,A)") "LARGE PARTICLE BLOOMS in iteration ",Iter
            IF(iPartBloom.gt.0) THEN
                WRITE(6,"(A,I10,A)") "A max of ",abs(iPartBloom)," particles created in one attempt from double excit."
            ELSE
                WRITE(6,"(A,I10,A)") "A max of ",abs(iPartBloom)," particles created in one attempt from single excit."
            ENDIF
        ENDIF


        rat=(TotWalkersNew+0.D0)/(MaxWalkersPart+0.D0)
        IF(rat.gt.0.95) THEN
            WRITE(6,*) "*WARNING* - Number of particles/determinants has increased to over 95% of MaxWalkersPart"
            CALL FLUSH(6)
        ENDIF

        IF(tRotoAnnihil) THEN
            ValidSpawned=ValidSpawned-1     !For rotoannihilation, this is the number of spawned particles.

            rat=(ValidSpawned+0.D0)/(MaxSpawned+0.D0)
            IF(rat.gt.0.9) THEN
                WRITE(6,*) "*WARNING* - Number of spawned particles has reached over 90% of MaxSpawned"
                CALL FLUSH(6)
            ENDIF
        ELSEIF(tDirectAnnihil) THEN
!Need to test whether any of the sublists are getting to the end of their allotted space.
            IF(nProcessors.gt.1) THEN
                do i=0,nProcessors-1
                    rat=(ValidSpawnedList(i)-InitialSpawnedSlots(i))/(InitialSpawnedSlots(1)+0.D0)
!                    WRITE(6,*) rat,(ValidSpawnedList(i)-InitialSpawnedSlots(i)),InitialSpawnedSlots(1)
                    IF(rat.gt.0.95) THEN
                        WRITE(6,*) "*WARNING* - Highest processor spawned particles has reached over 95% of MaxSpawned"
                        CALL FLUSH(6)
                    ENDIF
                enddo
            ELSE
                rat=(ValidSpawnedList(0)+0.D0)/(MaxSpawned+0.D0)
                IF(rat.gt.0.9) THEN
                    WRITE(6,*) "*WARNING* - Number of spawned particles has reached over 90% of MaxSpawned"
                    CALL FLUSH(6)
                ENDIF
            ENDIF
        ENDIF
        
        CALL halt_timer(Walker_Time)
        CALL set_timer(Annihil_Time,30)
!        CALL MPI_Barrier(MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get into annihilation"
!        CALL FLUSH(6)

        IF(tDirectAnnihil) THEN
!This is the direct annihilation algorithm. The newly spawned walkers should be in a seperate array (SpawnedParts) and the other list should be ordered.
            CALL DirectAnnihilation(TotWalkersNew)

        ELSEIF(tRotoAnnihil) THEN
!This is the rotoannihilation algorithm. The newly spawned walkers should be in a seperate array (SpawnedParts) and the other list should be ordered.

            CALL RotoAnnihilation(ValidSpawned,TotWalkersNew)

        ELSEIF(TNoAnnihil) THEN
!However, we now need to swap around the pointers of CurrentDets and NewDets, since this was done previously explicitly in the annihilation routine
            IF(associated(CurrentDets,target=WalkVecDets)) THEN
                CurrentDets=>WalkVec2Dets
                CurrentSign=>WalkVec2Sign
                CurrentH=>WalkVec2H
                CurrentExcits=>WalkVec2Excits
                NewDets=>WalkVecDets
                NewSign=>WalkVecSign
                NewH=>WalkVecH
                NewExcits=>WalkVecExcits
            ELSE
                CurrentDets=>WalkVecDets
                CurrentSign=>WalkVecSign
                CurrentH=>WalkVecH
                CurrentExcits=>WalkVecExcits
                NewDets=>WalkVec2Dets
                NewSign=>WalkVec2Sign
                NewH=>WalkVec2H
                NewExcits=>WalkVec2Excits
            ENDIF

            TotWalkers=TotWalkersNew
            TotParts=TotWalkers

        ELSEIF(TAnnihilonProc) THEN
!In this routine, the particles are just annihilated on their own processors. This means that all simulations are independent and not influenced by each other.
!This means that there is no communication between processors and so should be much faster as the system size increases.

            CALL Stop_All("PerformFCIMCyc","AnnihilonProc has been disabled")
!            CALL AnnihilateonProc(TotWalkersNew)
!            Annihilated=Annihilated+(TotWalkersNew-TotWalkers)
!            TotParts=TotWalkers

        ELSE
!This routine now cancels down the particles with opposing sign on each determinant

            IF(mod(iAnnInterval,Iter).eq.0) THEN

                CALL AnnihilatePartPar(TotWalkersNew)
                Annihilated=Annihilated+(TotWalkersNew-TotWalkers)
                TotParts=TotWalkers

            ELSE
                TotWalkersNew=TotWalkers
                TotParts=TotWalkers
                Annihilated=0
            ENDIF

        ENDIF

        CALL halt_timer(Annihil_Time)
        

!        IF(TSinglePartPhase) THEN
!!            CALL MPI_Barrier(MPI_COMM_WORLD,error)
!!Do not allow culling if we are still in the single particle phase.
!            IF(iProcIndex.eq.root) THEN     !Only exit phase if particle number is sufficient on head node.
!                IF(TotParts.gt.InitWalkers) THEN
!                    WRITE(6,*) "Exiting the single particle growth phase - shift can now change"
!                    TSinglePartPhase=.false.
!                ENDIF
!            ENDIF
!!Broadcast the fact that TSinglePartPhase may have changed to all processors - unfortunatly, have to do this broadcast every iteration
!            CALL MPI_Bcast(TSinglePartPhase,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
!        ELSE
!
!!Culling (either up or down) has now been disabled...
!
!            IF(TotWalkers.gt.(InitWalkers*GrowMaxFactor)) THEN
!!Particle number is too large - kill them randomly
!                IF(.not.tRotoAnnihil) THEN
!
!!Log the fact that we have made a cull
!                    NoCulls=NoCulls+1
!                    IF(NoCulls.gt.10) THEN
!                        WRITE(6,*) "Too Many Culls"
!                        CALL FLUSH(6)
!                        call Stop_All("PerformFCIMCyc","Too Many Culls")
!                    ENDIF
!!CullInfo(:,1) is walkers before cull
!                    CullInfo(NoCulls,1)=TotParts
!!CullInfo(:,3) is MC Steps into shift cycle before cull
!                    CullInfo(NoCulls,3)=mod(Iter,StepsSft)
!
!                    WRITE(6,"(A,F8.2,A)") "Total number of particles has grown to ",GrowMaxFactor," times initial number on this node..."
!                    WRITE(6,"(A,I12,A)") "Killing randomly selected particles in cycle ", Iter," in order to reduce total number on this node..."
!                    WRITE(6,"(A,F8.2)") "Population on this node will reduce by a factor of ",CullFactor
!                    CALL FLUSH(6)
!                    CALL ThermostatParticlesPar(.true.)
!
!                ENDIF
!
!            ELSEIF(TotParts.lt.(InitWalkers/2)) THEN
!!Particle number is too small - double every particle in its current position
!                IF(.not.tRotoAnnihil) THEN
!!Log the fact that we have made a cull
!                    NoCulls=NoCulls+1
!                    IF(NoCulls.gt.10) CALL Stop_All("PerformFCIMCyc","Too Many Culls")
!!CullInfo(:,1) is walkers before cull
!                    CullInfo(NoCulls,1)=TotParts
!!CullInfo(:,3) is MC Steps into shift cycle before cull
!                    CullInfo(NoCulls,3)=mod(Iter,StepsSft)
!                    
!                    WRITE(6,*) "Doubling particle population on this node to increase total number..."
!                    CALL ThermostatParticlesPar(.false.)
!                ELSE
!!                    WRITE(6,*) "Particle number on this node is less than half InitWalkers value"
!                ENDIF
!            ENDIF
!        
!        ENDIF

    END SUBROUTINE PerformFCIMCycPar

!This is a new annihilation algorithm. In this, determinants are kept on predefined processors, and newlyspawned particles are sent here so that all the annihilations are
!done on a predetermined processor, and not rotated around all of them.
    SUBROUTINE DirectAnnihilation(TotWalkersNew)
        INTEGER :: MaxIndex,TotWalkersNew

!This routine will send all the newly-spawned particles to their correct processor. MaxIndex is returned as the new number of newly-spawned particles on the processor. May have duplicates.
!The particles are now stored in SpawnedParts2/SpawnedSign2.
        CALL SendProcNewParts(MaxIndex)

!        WRITE(6,*) "Sent particles"
!        CALL FLUSH(6)

!CompressSpawnedList works on SpawnedParts arrays, so swap the pointers around.
        IF(associated(SpawnedParts2,target=SpawnVec2)) THEN
            SpawnedParts2 => SpawnVec
            SpawnedSign2 => SpawnSignVec
            SpawnedParts => SpawnVec2
            SpawnedSign => SpawnSignVec2
        ELSE
            SpawnedParts => SpawnVec
            SpawnedSign => SpawnSignVec
            SpawnedParts2 => SpawnVec2
            SpawnedSign2 => SpawnSignVec2
        ENDIF

!Now we want to order and compress the spawned list of particles. This will also annihilate the newly spawned particles amongst themselves.
!MaxIndex will change to reflect the final number of unique determinants in the newly-spawned list, and the particles will end up in the spawnedSign/SpawnedParts lists.
        CALL CompressSpawnedList(MaxIndex)

!        WRITE(6,*) "List compressed",MaxIndex,TotWalkersNew
!        CALL FLUSH(6)

!Binary search the main list and copy accross/annihilate determinants which are found.
!This will also remove the found determinants from the spawnedparts lists.

        CALL AnnihilateSpawnedParts(MaxIndex,TotWalkersNew)

!            WRITE(6,*) "Annihilation finished",MaxIndex,TotWalkersNew
!            CALL FLUSH(6)

!Put the surviving particles in the main list, maintaining order of the main list.
!Now we insert the remaining newly-spawned particles back into the original list (keeping it sorted), and remove the annihilated particles from the main list.
        CALL set_timer(Sort_Time,30)
        CALL InsertRemoveParts(MaxIndex,TotWalkersNew)

!            WRITE(6,*) "Surviving particles merged"
!            CALL FLUSH(6)

        CALL halt_timer(Sort_Time)

    END SUBROUTINE DirectAnnihilation

!This routine is used for sending the determinants to the correct processors. 
    SUBROUTINE SendProcNewParts(MaxIndex)
        REAL :: Gap
        INTEGER :: i,sendcounts(nProcessors),disps(nProcessors),recvcounts(nProcessors),recvdisps(nProcessors),error
        INTEGER :: MaxIndex,MaxSendIndex

!        WRITE(6,*) "ValidSpawnedList ",ValidSpawnedList(:)

        Gap=REAL(MaxSpawned)/REAL(nProcessors)
!        WRITE(6,*) "Gap: ",Gap
        do i=0,nProcessors-1
            sendcounts(i+1)=ValidSpawnedList(i)-(NINT(Gap*i)+1)
            disps(i+1)=NINT(Gap*i)
        enddo

        MaxSendIndex=ValidSpawnedList(nProcessors-1)-1

!We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
        recvcounts(1:nProcessors)=0
        
        CALL set_timer(Comms_Time,30)
        
        CALL MPI_AlltoAll(sendcounts,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,error)

!We can now get recvdisps from recvcounts, since we want the data to be contiguous after the move.
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo

        MaxIndex=recvdisps(nProcessors)+recvcounts(nProcessors)
!Max index is the largest occupied index in the array of hashes to be ordered in each processor 
        IF(MaxIndex.gt.(0.9*MaxSpawned)) THEN
            CALL Warning("SendProcNewParts","Maximum index of newly-spawned array is close to maximum length after annihilation send. Increase MemoryFacSpawn")
        ENDIF

!        WRITE(6,*) "sendcounts: ",sendcounts(:)
!        WRITE(6,*) "disps: ",disps(:)
!        WRITE(6,*) "recvcounts: ", recvcounts(:)
!        WRITE(6,*) "recvdisps: ",recvdisps(:)
!        WRITE(6,*) "Sent Sign: ", NINT(Gap),sendcounts(2)
!        do i=NINT(Gap),NINT(Gap)+sendcounts(2)
!            WRITE(6,*) i,"***",SpawnedSign(i)
!        enddo
        
!This is the main send of newly-spawned particles and signs to each determinants correct processor.
        CALL MPI_AlltoAllv(SpawnedSign(1:MaxSendIndex),sendcounts,disps,MPI_INTEGER,SpawnedSign2(1:MaxIndex),recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        
!        WRITE(6,*) MaxIndex, "Recieved signs: "
!        do i=1,MaxIndex
!            WRITE(6,*) SpawnedSign2(i)
!        enddo

!Update the number of integers we need to send.
        do i=1,nProcessors
            sendcounts(i)=sendcounts(i)*(NIfD+1)
            disps(i)=disps(i)*(NIfD+1)
            recvcounts(i)=recvcounts(i)*(NIfD+1)
            recvdisps(i)=recvdisps(i)*(NIfD+1)
        enddo

!        WRITE(6,*) "Sent Particles: ", NINT(Gap),sendcounts(2)
!        do i=NINT(Gap)+1,NINT(Gap)+sendcounts(2)
!            WRITE(6,*) i,"***",SpawnedParts(:,i)
!        enddo

        CALL MPI_AlltoAllv(SpawnedParts(0:NIfD,1:MaxSendIndex),sendcounts,disps,MPI_INTEGER,SpawnedParts2(0:NIfD,1:MaxIndex),recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)

!        WRITE(6,*) MaxIndex, "Recieved particles: "
!        do i=1,MaxSpawned
!            IF(SpawnedParts2(1,i).ne.0) THEN
!                WRITE(6,*) SpawnedParts2(:,i)
!            ENDIF
!        enddo
        
        CALL halt_timer(Comms_Time)

    END SUBROUTINE SendProcNewParts
    


    INTEGER FUNCTION DetermineDetProc(iLut)
        INTEGER :: iLut(0:NIfD),i,j,Elecs!,TempDet(NEl),MurmurHash2Wrapper
        INTEGER(KIND=i2) :: Summ!,RangeofBins,NextBin

!        CALL DecodeBitDet(TempDet,iLut,NEl,NIfD)
!        i=MurmurHash2Wrapper(TempDet,NEl,13)
!        write(6,*) i
        
        
        
        Summ=0
        Elecs=0
        lp2: do i=0,NIfD
            do j=0,31
                IF(BTEST(iLut(i),j)) THEN
                    Elecs=Elecs+1
                    Summ=(1099511628211_8*Summ)+((i*32)+(j+1))*Elecs
                    IF(Elecs.eq.NEl) EXIT lp2
                ENDIF
            enddo
        enddo lp2
        DetermineDetProc=abs(mod(Summ,INT(nProcessors,8)))
!        WRITE(6,*) DetermineDetProc,Summ,nProcessors

!        RangeofBins=NINT(HUGE(RangeofBins)/(nProcessors/2.D0),8)
!        NextBin=-HUGE(NextBin)+RangeofBins
!        do i=1,nProcessors
!            IF(i.eq.nProcessors) THEN
!!Make sure catch them all...
!                DetermineDetProc=nProcessors-1
!                RETURN
!            ENDIF
!            IF(Summ.gt.NextBin) THEN
!                NextBin=NextBin+RangeofBins
!            ELSE
!                DetermineDetProc=i-1
!                RETURN
!            ENDIF
!        enddo
!!Determine range by simply dividing hash...
!        IF(mod(nProcessors,2).ne.0) THEN
!            CALL Stop_All("DetermineDetProc","Number of processors must be a multiple of two for this hashing algorithm")
!        ENDIF
!        RangeofBins=NINT(HUGE(RangeofBins)/(nProcessors/2.D0),8)
!        IF(Summ.gt.0) THEN
!            DetermineDetProc=INT(((Summ+0.D0)/(RangeofBins+0.D0)),4)
!        ELSE
!            DetermineDetProc=INT(((abs(Summ)+0.D0)/(RangeofBins+0.D0)),4)+nProcessors/2
!        ENDIF
            
            

    END FUNCTION DetermineDetProc

!This sorts and compresses the spawned list to make it easier for the rest of the annihilation process.
!This is not essential, but should proove worthwhile
    SUBROUTINE CompressSpawnedList(ValidSpawned)
        INTEGER :: VecInd,ValidSpawned,DetsMerged,ToRemove,i,SignProd,PartIndex,ExcitLevel
        LOGICAL :: DetBitEQ,tSuc

!We want to sort the list of newly spawned particles, in order for quicker binary searching later on. (this is not essential, but should proove faster)
!They should remain sorted after annihilation between spawned
        CALL SortBitDets(ValidSpawned,SpawnedParts(0:NIfD,1:ValidSpawned),NIfD,SpawnedSign(1:ValidSpawned))
        IF(tHistSpawn) HistMinInd2(1:NEl)=FCIDetIndex(1:NEl)

!First, we compress the list of spawned particles, so that they are only specified at most once in each processors list.
!During this, we transfer the particles over to SpawnedParts2
        IF(ValidSpawned.gt.0) THEN
            SpawnedParts2(0:NIfD,1)=SpawnedParts(0:NIfD,1)
            SpawnedSign2(1)=SpawnedSign(1)
        ENDIF
        VecInd=1
        DetsMerged=0
        ToRemove=0
        do i=2,ValidSpawned
            IF(.not.DetBitEQ(SpawnedParts(0:NIfD,i),SpawnedParts2(0:NIfD,VecInd),NIfD)) THEN
                IF(SpawnedSign2(VecInd).eq.0) ToRemove=ToRemove+1
                VecInd=VecInd+1
                SpawnedParts2(:,VecInd)=SpawnedParts(:,i)
                SpawnedSign2(VecInd)=SpawnedSign(i)
            ELSE
                SignProd=SpawnedSign(i)*SpawnedSign2(VecInd)
                IF(SignProd.lt.0) THEN
!We are actually unwittingly annihilating, but just in serial... we therefore need to count it anyway.
                    Annihilated=Annihilated+2*(MIN(abs(SpawnedSign2(VecInd)),abs(SpawnedSign(i))))

                    IF(tHistSpawn) THEN
!We want to histogram where the particle annihilations are taking place.
                        CALL FindBitExcitLevel(SpawnedParts(:,i),iLutHF,NIfD,ExcitLevel,NEl)
                        IF(ExcitLevel.eq.NEl) THEN
                            CALL BinSearchParts2(SpawnedParts(:,i),HistMinInd2(ExcitLevel),Det,PartIndex,tSuc)
                        ELSEIF(ExcitLevel.eq.0) THEN
                            PartIndex=1
                            tSuc=.true.
                        ELSE
                            CALL BinSearchParts2(SpawnedParts(:,i),HistMinInd2(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartIndex,tSuc)
                        ENDIF
                        HistMinInd2(ExcitLevel)=PartIndex
                        IF(tSuc) THEN
                            AvAnnihil(PartIndex)=AvAnnihil(PartIndex)+REAL(2*(MIN(abs(SpawnedSign2(VecInd)),abs(SpawnedSign(i)))),r2)
                            InstAnnihil(PartIndex)=InstAnnihil(PartIndex)+REAL(2*(MIN(abs(SpawnedSign2(VecInd)),abs(SpawnedSign(i)))),r2)
                        ELSE
!                            WRITE(6,*) "Searching between: ",HistMinInd2(ExcitLevel), " and ",FCIDetIndex(ExcitLevel+1)-1
!                            WRITE(6,*) "***",SpawnedParts(0:NIfD,i)
!                            CALL DecodeBitDet(TempDet,SpawnedParts(0:NIfD,i),NEl,NIfD)
!                            WRITE(6,*) "Full Det is: ",TempDet(:)
!                            IF(tHPHF) THEN
!                                CALL FindExcitBitDetSym(SpawnedParts(0:NIfD,i),iLutSym(:))
!                                WRITE(6,*) "*** Sym: ",iLutSym(:)
!                                CALL DecodeBitDet(TempDet,iLutSym(0:NIfD),NEl,NIfD)
!                                WRITE(6,*) "Full Sym Det is: ",TempDet(:)
!                            ENDIF
                            CALL Stop_All("CompressSpawnedList","Cannot find corresponding FCI determinant when histogramming")
                        ENDIF
                    ENDIF
                ENDIF
                SpawnedSign2(VecInd)=SpawnedSign2(VecInd)+SpawnedSign(i)
                DetsMerged=DetsMerged+1
            ENDIF
        enddo

        ValidSpawned=ValidSpawned-DetsMerged
        IF((ValidSpawned.ne.VecInd).and.(VecInd.ne.1)) THEN
            WRITE(6,*) ValidSpawned,VecInd
            CALL Stop_All("CompressSpawnedList","Error in compression of spawned particle list")
        ENDIF
        IF(SpawnedSign2(ValidSpawned).eq.0.and.(ValidSpawned.gt.0)) ToRemove=ToRemove+1

!Now remove zeros. Not actually necessary, but will be useful I suppose? Shouldn't be too much hassle.
!We can also use it to copy the particles back to SpawnedParts array
        DetsMerged=0
        do i=1,ValidSpawned
            IF(SpawnedSign2(i).eq.0) THEN
                DetsMerged=DetsMerged+1
            ELSE
!We want to move all the elements above this point down to 'fill in' the annihilated determinant.
                SpawnedParts(0:NIfD,i-DetsMerged)=SpawnedParts2(0:NIfD,i)
                SpawnedSign(i-DetsMerged)=SpawnedSign2(i)
            ENDIF
        enddo
        IF(DetsMerged.ne.ToRemove) THEN
            CALL Stop_All("CompressSpawnedList","Wrong number of entries removed from spawned list")
        ENDIF
        ValidSpawned=ValidSpawned-DetsMerged
        
    END SUBROUTINE CompressSpawnedList

    
!This is a new routine to totally annihilate all particles on the same determinant. This is not done using an all-to-all, but rather
!by rotating the newly spawned particles around all determinants and annihilating with the particles on their processor.
!Valid spawned is the number of newly-spawned particles. These 'particles' can be multiply specified on the same determinant.
!Each rotation and annihilation step, the number corresponds to a different processors spawned particles.
!TotWalkersNew indicates the number of particles in NewDets - the list of particles to compare for annihilation.
!Improvements in AnnihilateBetweenSpawned:
!Binary search for sendcounts and others, and only transfer all data when need to.
!Memory improvements
!Call as one array for All-to-alls
!Make sure only sort what need to
    SUBROUTINE RotoAnnihilation(ValidSpawned,TotWalkersNew)
        INTEGER :: ValidSpawned,TotWalkersNew,i
        INTEGER :: ierr,error!,SpawnedBeforeRoto
        CHARACTER , ALLOCATABLE :: mpibuffer(:)

!        InitialSpawned=TotSpawned     !Initial spawned will store the original number of spawned particles, so that we can compare afterwards.
!        InitialSpawned=Annihilated
        
        CALL CompressSpawnedList(ValidSpawned)

!        CALL SortBitDets(ValidSpawned,SpawnedParts(0:NIfD,1:ValidSpawned),NIfD,SpawnedSign(1:ValidSpawned))
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
!        WRITE(6,*) "Entering rotoannilation: ",Iter,InitialSpawned,TotWalkersNew
!        CALL FLUSH(6)

!First, annihilate between newly spawned particles. Memory for this will be allocated dynamically.
!This will be done in the usual fashion using the All-to-All communication and hashes.
        CALL AnnihilateBetweenSpawned(ValidSpawned)
!        CALL AnnihilateBetweenSpawnedOneProc(ValidSpawned)
!        Annihilated=Annihilated+(InitialSpawned-TotSpawned)
!        IF(Annihilated.ne.InitialSpawned) THEN
!            WRITE(6,*) "Have annihilated between newly-spawned...",Annihilated-InitialSpawned,Iter
!        ENDIF
            
!We want to sort the list of newly spawned particles, in order for quicker binary searching later on. (this is not essential, but should proove faster)
!        CALL SortBitDets(ValidSpawned,SpawnedParts(0:NIfD,1:ValidSpawned),NIfD,SpawnedSign(1:ValidSpawned))
!        CALL CheckOrdering(SpawnedParts(:,1:ValidSpawned),SpawnedSign(1:ValidSpawned),ValidSpawned,.true.)
!        do i=1,ValidSpawned
!            WRITE(6,*) 1,i,SpawnedParts(:,i),SpawnedSign(i),Iter
!            CALL FLUSH(6)
!        enddo

!        SpawnedBeforeRoto=ValidSpawned
!        WRITE(6,*) "SpawnedBeforeRoto: ",ValidSpawned

!This RemoveInds is useful scratch space for the removal of particles from lists. It probably isn't essential, but keeps things simpler initially.
!        ALLOCATE(RemoveInds(MaxSpawned),stat=ierr)
        
!This routine annihilates the processors set of newly-spawned particles, with the complete set of particles on the processor.
        CALL AnnihilateSpawnedParts(ValidSpawned,TotWalkersNew)

        CALL MPI_Barrier(MPI_COMM_WORLD,error)

!Allocate a buffer here to hold particles when using a buffered send...
!The buffer wants to be able to hold (MaxSpawned+1)x(NIfD+2) integers (*4 for in bytes). If we could work out the maximum ValidSpawned accross the determinants,
!it could get reduced to this... 
        IF(nProcessors.ne.1) THEN
            ALLOCATE(mpibuffer(8*(MaxSpawned+1)*(NIfD+2)),stat=ierr)
            IF(ierr.ne.0) THEN
                CALL Stop_All("RotoAnnihilation","Error allocating memory for transfer buffers...")
            ENDIF
            CALL MPI_Buffer_attach(mpibuffer,8*(MaxSpawned+1)*(NIfD+2),error)
            IF(error.ne.0) THEN
                CALL Stop_All("RotoAnnihilation","Error allocating memory for transfer buffers...")
            ENDIF
        ENDIF

        do i=1,nProcessors-1
!Move newly-spawned particles which haven't been annihilated around the processors in sequence, annihilating locally each step.
!This moves the set of newly-spawned particles on this processor one to the right, and recieves from the left.
!This also updates the ValidSpawned variable so that it now refers to the new set of spawned-particles.
            CALL RotateParticles(ValidSpawned)
!            WRITE(6,*) "Rotating particles for the ",i," time...",Iter
!            CALL FLUSH(6)

!This routine annihilates the processors set of newly-spawned particles, with the complete set of particles on the processor.
            CALL AnnihilateSpawnedParts(ValidSpawned,TotWalkersNew)
!            CALL CheckOrdering(SpawnedParts(:,1:ValidSpawned),SpawnedSign(1:ValidSpawned),ValidSpawned,.true.)
!            WRITE(6,*) "Annihilated locally....",i
!            CALL FLUSH(6)

        enddo

!One final rotation means that the particles are all on their original processor.
        IF(nProcessors.ne.1) THEN
            CALL RotateParticles(ValidSpawned)

!Detach buffers
            CALL MPI_Buffer_detach(mpibuffer,8*(MaxSpawned+1)*(NIfD+2),error)
            DEALLOCATE(mpibuffer)
        ENDIF
        
!If we are using a guiding function, then we want to attempt to annihilate newly spawned particles with the guiding function here.
        IF(tUseGuide) CALL RotoAnnihilGuidingFunc(ValidSpawned)

!Test that we have annihilated the correct number here (from each lists), and calculate Annihilated for each processor.
!Now we insert the remaining newly-spawned particles back into the original list (keeping it sorted), and remove the annihilated particles from the main list.
        CALL set_timer(Sort_Time,30)
!        WRITE(6,*) "Entering insert/remove..."
!        CALL FLUSH(6)
        CALL InsertRemoveParts(ValidSpawned,TotWalkersNew)
        CALL halt_timer(Sort_Time)

!        DEALLOCATE(RemoveInds)


    END SUBROUTINE RotoAnnihilation


! This routine is based on RotoAnnihilation (with a wee bit of AnnihilatePartPar).
! It first takes MinorSpawnDets and orders the determinants (note, no compression is needed, determinants may be listed more than once, but these will have
! different parents).
! It then runs through these spawned walkers, and annihilates amongst the spawned.
! It then does a rotation around the processors, annihilating with the MinorStarDets.
! Any walkers which survive this are then added to MinorStarDets, maintaining order.
! MinorValidSpawned is the number of newly spawned walkers on the minor determinants whereas NoMinorWalkersNew are the walkers on the minor dets that have survived previous 
! iterations.  NoMinorWalkersNew is the number to be compared to for each processor.
    SUBROUTINE RotoAnnihilateMinorSpawned(MinorValidSpawned,NoMinorWalkersNew)
        INTEGER :: i,j,MinorValidSpawned,NoMinorWalkersNew,n,error,ierr
        CHARACTER , ALLOCATABLE :: mpibuffer(:)

! First order the newly spawned walkers in terms of determinant, then parent, taking the sign and H element information with it.
        CALL Sort2BitDetsPlus3(MinorValidSpawned,MinorSpawnDets(0:NIfD,1:MinorValidSpawned),NIfD,MinorSpawnParent(0:NIfD,1:MinorValidSpawned),&
        &NIfD,MinorSpawnSign(1:MinorValidSpawned))
        
!        IF(Iter.gt.1220) THEN
!            WRITE(6,*) 'sort'
!            CALL FLUSH(6)
!        ENDIF

! Make sure all processors have done this before carrying on.
        CALL MPI_Barrier(MPI_COMM_WORLD,error)


! Run through this list of determinants with walkers on it, and annihilate walkers on the same determinant.  Make sure the correct parent information is kept with
! the walkers that survive.
! Then need to communicate between processors, and annihilate the within the spawned particles, between processors.
! At the end, it should be that walkers on the same determinants have the same sign across all processors.
       
      
        CALL AnnihilateAmongstMinorSpawned(MinorValidSpawned)
 
!        IF(Iter.gt.1220) THEN
!            WRITE(6,*) 'annihilate amongst'
!            CALL FLUSH(6)
!        ENDIF


! Should now have the spawned walkers ordered in terms of determinant then parent, with each determinant/parent combination only 
! specified once.

! At the end, should either have each determinant specified once, or more than once with different parents but the same sign.        

! Annihilate with the MinorStarDets on the original processor        
! If multiple entries of the determinant on which we are annihilating - want to kind of add up these walkers, then randomly select the ones to annihilate.


        CALL AnnihilateMinorSpawnedParts(MinorValidSpawned,NoMinorWalkersNew)
 
!        IF(Iter.gt.1220) THEN
!            WRITE(6,*) 'annihilate minor'
!            CALL FLUSH(6)
!        ENDIF


!Allocate a buffer here to hold particles when using a buffered send...
!The buffer wants to be able to hold (MaxSpawned+1)x(NIfD+2) integers (*4 for in bytes). If we could work out the maximum ValidSpawned accross the determinants,
!it could get reduced to this... 
        IF(nProcessors.ne.1) THEN
            ALLOCATE(mpibuffer(8*(MaxSpawned+1)*(NIfD+3)),stat=ierr)
            IF(ierr.ne.0) THEN
                CALL Stop_All("RotoAnnihilateMinor","Error allocating memory for transfer buffers...")
            ENDIF
            CALL MPI_Buffer_attach(mpibuffer,8*(MaxSpawned+1)*(NIfD+3),error)
            IF(error.ne.0) THEN
                CALL Stop_All("RotoAnnihilateMinor","Error allocating memory for transfer buffers...")
            ENDIF
        ENDIF
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        do n=1,nProcessors-1

! Take the walkers that survive the annihilation amongst spawned particles and rotate them around each processor, annihilating with MinorStarDets etc.
            CALL RotateMinorParticles(MinorValidSpawned)
     
!            IF(Iter.gt.1220) THEN
!                WRITE(6,*) 'rotate'
!                CALL FLUSH(6)
!            ENDIF


            CALL AnnihilateMinorSpawnedParts(MinorValidSpawned,NoMinorWalkersNew)
     
!            IF(Iter.gt.1220) THEN
!                WRITE(6,*) 'annihilate minor'
!                CALL FLUSH(6)
!            ENDIF

        enddo


! Then do one final rotation (if nProcessors.gt.1) to get back to the original processor, and add the survivors into MinorStarDets (MinorStarDets will not have any contribution 
! to the energy - may want to put in clause that we cannot select the dominant 2s).
        IF(nProcessors.gt.1) THEN

            CALL MPI_Barrier(MPI_COMM_WORLD,error)
            CALL RotateMinorParticles(MinorValidSpawned)
 
!            IF(Iter.gt.1220) THEN
!                WRITE(6,*) 'last rotate'
!                CALL FLUSH(6)
!            ENDIF
!Detach buffers
            CALL MPI_Buffer_detach(mpibuffer,8*(MaxSpawned+1)*(NIfD+3),error)
            DEALLOCATE(mpibuffer)

        ENDIF

        CALL InsertRemoveMinorParts(MinorValidSpawned,NoMinorWalkersNew)

!        IF(Iter.gt.1220) THEN
!            WRITE(6,*) 'insert remove'
!            CALL FLUSH(6)
!        ENDIF



    ENDSUBROUTINE RotoAnnihilateMinorSpawned




    SUBROUTINE AnnihilateAmongstMinorSpawned(MinorValidSpawned)
! This routine takes the newly spawned walkers on the minor determinants, and annihilates those on the same determinant    
        INTEGER :: i,j,k,ToAnnihilateIndex,MinorValidSpawned,MinorValidSpawnedNew,ierr,error,sendcounts(nProcessors)
        INTEGER :: TotWalkersDet,InitialBlockIndex,FinalBlockIndex,ToAnnihilateOnProc,VecSlot
        INTEGER :: disps(nProcessors),recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER :: Minsendcounts,Maxsendcounts,DebugIter,SubListInds(2,nProcessors),MinProc,MinInd
        INTEGER :: SumOppSign,WalkertoAnnihil,NoNegWalk,NoPosWalk
        REAL*8 :: r
        INTEGER(KIND=i2) :: HashCurr,MinBin,RangeofBins,NextBinBound,MinHash
        INTEGER(KIND=i2) , ALLOCATABLE :: TempHash(:)
        INTEGER , ALLOCATABLE :: TempSign(:),TempMinorSpawnSign(:)                                                      
        LOGICAL :: tWrite
        CHARACTER(len=*), PARAMETER :: this_routine='AnnihilateAmongstMinorSpawned'


!First, allocate memory to hold the signs and the hashes while we annihilate
        ALLOCATE(TempSign(MinorValidSpawned),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine, 'problem allocating memory to tempsign')
!TempMinorSpawnSign is just MinorSpawnSign, but MinorSpawnSign is used as a kind of AllMinorSpawnSign.
        ALLOCATE(TempMinorSpawnSign(MaxSpawned),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine, 'problem allocating memory to tempminorspawnsign')

!Temporary arrays, storing the signs these are going to be mixed.  The hashes are also mixed, but these are not needed after
!so are not reordered.
        TempSign(1:MinorValidSpawned)=MinorSpawnSign(1:MinorValidSpawned)
        TempMinorSpawnSign(1:MinorValidSpawned)=MinorSpawnSign(1:MinorValidSpawned)
        MinorSpawnSign(:)=0

!Create the arrays for index and process
        do i=1,MinorValidSpawned
            IndexTable(i)=i
        enddo
        ProcessVec(1:MinorValidSpawned)=iProcIndex

!Next, order the hash array, taking the index, CPU and sign with it...
!Order the array by abs(mod(Hash,nProcessors)). This will result in a more load-balanced system

        CALL Sort4ILong(MinorValidSpawned,HashArray(1:MinorValidSpawned),IndexTable(1:MinorValidSpawned),ProcessVec(1:MinorValidSpawned),TempMinorSpawnSign(1:MinorValidSpawned))
!Hash's ordered, taking index, ProcessVec and sign with them.  Forget determinants, they're just determined by their hash now.

        IF(nProcessors.ne.1) THEN
!We also need to know the ranges of the hashes to send to each processor. Each range should be the same.
            Rangeofbins=INT(HUGE(Rangeofbins)/(nProcessors/2),8)
            MinBin=-HUGE(MinBin)
            NextBinBound=MinBin+Rangeofbins

!We need to find the indices for each block of hashes which are to be sent to each processor.
!Sendcounts is the size of each block of ordered dets which are going to each processors. This could be binary searched for extra speed.
            j=1
            do i=1,nProcessors    !Search through all possible values of the hashes
                do while((HashArray(j).le.NextBinBound).and.(j.le.MinorValidSpawned))
                    j=j+1
                enddo
                sendcounts(i)=j-1
                IF(i.eq.nProcessors-1) THEN
!Make sure the final bin catches everything...
                    NextBinBound=HUGE(NextBinBound)
                ELSE
                    NextBinBound=NextBinBound+Rangeofbins
                ENDIF
            enddo
        ELSE
            sendcounts(1)=MinorValidSpawned
        ENDIF

        IF(sendcounts(nProcessors).ne.MinorValidSpawned) THEN
            WRITE(6,*) "SENDCOUNTS is: ",sendcounts(:)
            WRITE(6,*) "TOTWALKERSNEW is: ",MinorValidSpawned
            CALL FLUSH(6)
            CALL Stop_All("RotoAnnihilateMinorSpawned","Incorrect calculation of sendcounts")
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

        CALL MPI_AlltoAll(sendcounts(1:nProcessors),1,MPI_INTEGER,recvcounts(1:nProcessors),1,MPI_INTEGER,MPI_COMM_WORLD,error)

!We can now get recvdisps from recvcounts in the same way we obtained disps from sendcounts
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo

        MaxIndex=recvdisps(nProcessors)+recvcounts(nProcessors)
!Max index is the largest occupied index in the array of hashes to be ordered in each processor 
        IF(MaxIndex.gt.(0.95*MaxSpawned)) THEN
            CALL Warning("AnnihilateAmongstMinorSpawned","Maximum index of annihilation array is close to maximum length. Increase MemoryFacAnnihil")
        ENDIF
!Uncomment this if you want to write out load-balancing statistics.
!        AnnihilPart(:)=0
!        CALL MPI_Gather(MaxIndex,1,MPI_INTEGER,AnnihilPart,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
!        IF(iProcIndex.eq.root) THEN
!            WRITE(13,"(I10)",advance='no') Iter
!            do i=1,nProcessors
!                WRITE(13,"(I10)",advance='no') AnnihilPart(i)
!            enddo
!            WRITE(13,"(A)") ""
!            CALL FLUSH(13)
!        ENDIF

!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "RECVCOUNTS: "
!            WRITE(6,*) recvcounts(:)
!            WRITE(6,*) "RECVDISPS: "
!            WRITE(6,*) recvdisps(:),MaxIndex
!            CALL FLUSH(6)
!        ENDIF

!Insert a load-balance check here...maybe find the s.d. of the sendcounts array - maybe just check the range first.
!        IF(TotWalkersNew.gt.200) THEN
!            IF((Maxsendcounts-Minsendcounts).gt.(TotWalkersNew/3)) THEN
!                WRITE(6,"(A,I12)") "**WARNING** Parallel annihilation not optimally balanced on this node, for iter = ",Iter
!                WRITE(6,*) "Sendcounts is: ",sendcounts(:)
!!                CALL FLUSH(6)
!            ENDIF
!        ENDIF
!
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

!Now send the chunks of hashes to the corresponding processors
!All the '2' arrays are like the 'All' arrays.
!TempMinorSpawnSign is the Signs from each processor, when just MinorSpawnSign is the 'All' array.
        CALL MPI_AlltoAllv(HashArray(1:MinorValidSpawned),sendcounts,disps,MPI_DOUBLE_PRECISION,Hash2Array(1:MaxIndex),recvcounts,recvdisps,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,error)        

!        tWrite=.false.
!        IF(MinorValidSpawned.gt.3) THEN
!            WRITE(6,*) 'TempMinorSpawnSign'
!            do i=1,MinorValidSpawned
!                WRITE(6,*) TempMinorSpawnSign(i)
!            enddo
!            tWrite=.true.
!        ENDIF

!The signs of the hashes, index and CPU also need to be taken with them.
        CALL MPI_AlltoAllv(TempMinorSpawnSign(1:MinorValidSpawned),sendcounts,disps,MPI_INTEGER,MinorSpawnSign(1:MaxIndex),recvcounts,recvdisps,MPI_LOGICAL,MPI_COMM_WORLD,error)
        CALL MPI_AlltoAllv(IndexTable(1:MinorValidSpawned),sendcounts,disps,MPI_INTEGER,Index2Table(1:MaxIndex),recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        CALL MPI_AlltoAllv(ProcessVec(1:MinorValidSpawned),sendcounts,disps,MPI_INTEGER,Process2Vec(1:MaxIndex),recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        
!        IF(tWrite) THEN
!            WRITE(6,*) 'MinorSpawnSign'
!            do i=1,20
!                WRITE(6,*) MinorSpawnSign(i)
!            enddo
!            CALL FLUSH(6)
!            CALL Stop_All('','')
!        ENDIF



!Now we need to perform the actual annihilation, running through all the particles and calculating which ones want to be annihilated.

!Here, because we have ordered the hashes initially numerically, we have a set of ordered lists. It is therefore easier to sort them.
!We have to work out how to run sequentially through the hashes, which are a set of nProc seperate ordered lists.
!We would need to have 2*nProc indices, since we will have a set of nProc disjoint ordered sublists.
!SubListInds(1,iProc)=index of current hash from processor iProc
!SubListInds(2,iProc)=index of final hash from processor iProc
!Indices can be obtained from recvcounts and recvdisps - recvcounts(iProc-1) is number of hashes from iProc
!recvdisps(iProc-1) is the displacement to the start of the hashes from iProc
        do i=1,nProcessors-1
            SubListInds(1,i)=recvdisps(i)+1
            SubListInds(2,i)=recvdisps(i+1)
        enddo
        SubListInds(1,nProcessors)=recvdisps(nProcessors)+1
        SubListInds(2,nProcessors)=MaxIndex

!Reorder the lists so that they are in numerical order.
        j=1
        do while(j.le.MaxIndex)
            do i=1,nProcessors
                IF(SubListInds(1,i).le.SubListInds(2,i)) THEN
!This block still has hashes which want to be sorted
                    MinHash=Hash2Array(SubListInds(1,i))
                    MinProc=i
                    MinInd=SubListInds(1,i)
                    EXIT
                ENDIF
!                    IF(i.eq.nProcessors) THEN
!                        WRITE(6,*) "ERROR HERE!!"
!                        CALL FLUSH(6)
!                    ENDIF
            enddo
            IF(MinHash.ne.HashCurr) THEN
                do i=MinProc+1,nProcessors
                    IF((SubListInds(1,i).le.SubListInds(2,i)).and.(Hash2Array(SubListInds(1,i)).lt.MinHash)) THEN
                        MinHash=Hash2Array(SubListInds(1,i))
                        MinProc=i
                        MinInd=SubListInds(1,i)
                        IF(MinHash.eq.HashCurr) THEN
                            EXIT
                        ENDIF
                    ENDIF
                enddo
            ENDIF
!Next smallest hash is MinHash - move the ordered elements into the other array.
            HashArray(j)=MinHash
            IndexTable(j)=Index2Table(MinInd)
            ProcessVec(j)=Process2Vec(MinInd)
            TempMinorSpawnSign(j)=MinorSpawnSign(MinInd)
            HashCurr=MinHash
!Move through the block
            j=j+1
            SubListInds(1,MinProc)=SubListInds(1,MinProc)+1
        enddo

        IF((j-1).ne.MaxIndex) THEN
            CALL Stop_All(this_routine,"Error here in the merge sort algorithm")
        ENDIF

!Need to copy the lists back to the original array
        do i=1,MaxIndex
            Index2Table(i)=IndexTable(i)
            Process2Vec(i)=ProcessVec(i)
            MinorSpawnSign(i)=TempMinorSpawnSign(i)
            Hash2Array(i)=HashArray(i)
        enddo
 
!            Index2Table(1:MaxIndex)=IndexTable(1:MaxIndex)
!            Process2Vec(1:MaxIndex)=ProcessVec(1:MaxIndex)
!            CurrentSign(1:MaxIndex)=NewSign(1:MaxIndex)
!            HashArray(1:MaxIndex)=Hash2Array(1:MaxIndex)
                
!        WRITE(6,*) 'MinorSpawnSign'
!        do i=1,MaxIndex
!            WRITE(6,*) MinorSpawnSign(i)
!        enddo

!Work out the index of the particles which want to be annihilated
        j=1
        ToAnnihilateIndex=1
        do while(j.le.MaxIndex)
            TotWalkersDet=0
            NoPosWalk=0
            NoNegWalk=0
            InitialBlockIndex=j
            FinalBlockIndex=j-1         !Start at j-1 since we are increasing FinalBlockIndex even with the first det in the next loop
            HashCurr=Hash2Array(j)
            do while((Hash2Array(j).eq.HashCurr).and.(j.le.MaxIndex))
!                WRITE(6,*) 'Hash2Array',Hash2Array(j)
!                WRITE(6,*) 'HashCurr',HashCurr
!                WRITE(6,*) 'MinorSpawnSign',MinorSpawnSign(j)
!First loop counts walkers in the block - TotWalkersDet is then the residual sign of walkers on that determinant
                TotWalkersDet=TotWalkersDet+MinorSpawnSign(j)
                IF(MinorSpawnSign(j).gt.0) NoPosWalk=NoPosWalk+ABS(MinorSpawnSign(j))
                IF(MinorSpawnSign(j).lt.0) NoNegWalk=NoNegWalk+ABS(MinorSpawnSign(j))
! These will just annihilate each other until TotWalkersDet is the Total number of walkers (w sign) that should remain on that determinant.
!                IF(MinorSpawnSign(j).eq.1) THEN
!                    TotWalkersDet=TotWalkersDet+1
!                ELSE
!                    TotWalkersDet=TotWalkersDet-1
!                ENDIF
                FinalBlockIndex=FinalBlockIndex+1
                j=j+1
            enddo

!            IF((NoPosWalk.gt.0).and.(NoNegWalk.gt.0)) THEN
!                WRITE(6,*) 'NoPosWalk gt 0 and NoNegWalk gt 0'
!                WRITE(6,*) 'Index,hash,and Sign'
!                do i=InitialBlockIndex,FinalBlockIndex
!                    WRITE(6,*) i,Index2Table(i),IndexTable(i),hash2array(i),hasharray(i),MinorSpawnSign(i)
!                enddo
!                WRITE(6,*) 'NoPosWalk,',NoPosWalk,'NoNegWalk',NoNegWalk
!                WRITE(6,*) 'TotWalkersDet',TotWalkersDet
!            ENDIF

!Second run through the block of same determinants marks walkers for annihilation
            IF((TotWalkersDet.eq.0).and.(NoPosWalk.gt.0)) THEN

                do k=InitialBlockIndex,FinalBlockIndex
!All walkers in block want to be annihilated from now on.
                    IndexTable(ToAnnihilateIndex)=Index2Table(k)
                    ProcessVec(ToAnnihilateIndex)=Process2Vec(k)
!                    Hash2Array(ToAnnihilateIndex)=HashArray(k)     !This is not strictly needed - remove after checking
!                    NewSign(ToAnnihilateIndex)=CurrentSign(k)       !This is also need needed, but useful for checking
!                    IF(Iter.eq.DebugIter) WRITE(6,*) "Annihilating from if block 1",j,k
                    ToAnnihilateIndex=ToAnnihilateIndex+1
!                    WRITE(6,*) 'adding to annihilate index 01'
!                    CALL FLUSH(6)
!                    IF(HashCurr.eq.HFHash) THEN
!                        WRITE(6,*) "HF Determinant particle annihilated"
!                    ENDIF
                enddo   

!            ELSEIF((TotWalkersDet.lt.0).and.(MinorSpawnSign(k).gt.0)) THEN
            ELSEIF(TotWalkersDet.lt.0) THEN
!Need to run through the determinants, find those with positive walkers, and randomly annihilate these.            
                
                do while (NoPosWalk.gt.0)
!                    WRITE(6,*) 'into this loop'

                    ! call a random number between 1 and 0.
                    IF(tMerTwist) THEN
                        CALL genrand_real2(r) 
                    ELSE
                        CALL RANLUX(r,1)
                    ENDIF

                    ! multiply this by the number we need to annihilate, and the round up to the nearest integer.
                    ! this integer indicates the walker we need to annihilate.
                    WalkertoAnnihil=CEILING(r*NoPosWalk)
                    SumOppSign=0
                    do k=InitialBlockIndex,FinalBlockIndex

                        IF(MinorSpawnSign(k).gt.0) SumOppSign=SumOppSign+ABS(MinorSpawnSign(k))
                        IF(SumOppSign.ge.WalkertoAnnihil) THEN
                            MinorSpawnSign(k)=MinorSpawnSign(k)-1
                            NoPosWalk=NoPosWalk-1
                            IndexTable(ToAnnihilateIndex)=Index2Table(k)
                            ProcessVec(ToAnnihilateIndex)=Process2Vec(k)
!                    Hash2Array(ToAnnihilateIndex)=HashArray(k)     !This is not strictly needed - remove after checking
!                    NewSign(ToAnnihilateIndex)=CurrentSign(k)       !This is also need needed, but useful for checking
!                    IF(Iter.eq.DebugIter) WRITE(6,*) "Annihilating from if block 2",j,k
                            ToAnnihilateIndex=ToAnnihilateIndex+1

!                            WRITE(6,*) 'adding to annihilate index 02'
!                            CALL FLUSH(6)

                            EXIT
                        ENDIF
                    enddo
                enddo

            ELSEIF(TotWalkersDet.gt.0) THEN
!Annihilate if block has a net positive walker count, and current walker is negative
                do while (NoNegWalk.gt.0)

                    ! call a random number between 1 and 0.
                    IF(tMerTwist) THEN
                        CALL genrand_real2(r) 
                    ELSE
                        CALL RANLUX(r,1)
                    ENDIF

                    ! multiply this by the number we need to annihilate, and the round up to the nearest integer.
                    ! this integer indicates the walker we need to annihilate.
                    WalkertoAnnihil=CEILING(r*NoNegWalk)
                    SumOppSign=0
                    do k=InitialBlockIndex,FinalBlockIndex

                        IF(MinorSpawnSign(k).lt.0) SumOppSign=SumOppSign+ABS(MinorSpawnSign(k))
                        IF(SumOppSign.ge.WalkertoAnnihil) THEN
                            MinorSpawnSign(k)=MinorSpawnSign(k)+1
                            NoNegWalk=NoNegWalk-1
                            IndexTable(ToAnnihilateIndex)=Index2Table(k)
                            ProcessVec(ToAnnihilateIndex)=Process2Vec(k)
!                    Hash2Array(ToAnnihilateIndex)=HashArray(k)     !This is not strictly needed - remove after checking
!                    NewSign(ToAnnihilateIndex)=CurrentSign(k)       !This is also need needed, but useful for checking
!                    IF(Iter.eq.DebugIter) WRITE(6,*) "Annihilating from if block 3",j,k
                            ToAnnihilateIndex=ToAnnihilateIndex+1

!                            WRITE(6,*) 'adding to annihilate index 03'
!                            CALL FLUSH(6)


                            EXIT
                        ENDIF
                    enddo
                enddo
            ENDIF

!            IF((ToAnnihilateIndex).gt.1) THEN
!                WRITE(6,*) '** Toannihilateindex gt 1'
!                WRITE(6,*) 'InitialBlockIndex,',InitialBlockIndex,'FinalBlockIndex',FinalBlockIndex
!
!                WRITE(6,*) 'index,hash,sign'
!                do i=InitialBlockIndex,FinalBlockIndex
!                    WRITE(6,*) i,Hash2Array(i),MinorSpawnSign(i)
!                enddo
!                WRITE(6,*) 'ToAnnihilateIndex',ToAnnihilateIndex
!                do i=1,ToAnnihilateIndex-1
!                    WRITE(6,*) IndexTable(i),ProcessVec(i),Hash2Array(i),MinorSpawnSign(i)
!                enddo
!                CALL Stop_All('','')
!            ENDIF


        enddo


        ToAnnihilateIndex=ToAnnihilateIndex-1   !ToAnnihilateIndex now tells us the total number of particles to annihilate from the list on this processor
!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "Number of particles to annihilate from hashes on this processor: ",ToAnnihilateIndex
!            CALL FLUSH(6)
!        ENDIF
        MinorAnnihilated=MinorAnnihilated+ToAnnihilateIndex
        Annihilated=Annihilated+ToAnnihilateIndex

!The annihilation is complete - particles to be annihilated are stored in IndexTable and need to be sent back to their original processor
!To know which processor that is, we need to order the particles to be annihilated in terms of their CPU, i.e. ProcessVec(1:ToAnnihilateIndex)
!Is the list already ordered according to CPU? Is this further sort even necessary?

        IF(ToAnnihilateIndex.gt.1) THEN
!Do not actually have to take indextable, hash2array or newsign with it...
            CALL Sort2IILongI(ToAnnihilateIndex,ProcessVec(1:ToAnnihilateIndex),IndexTable(1:ToAnnihilateIndex),HashArray(1:ToAnnihilateIndex),MinorSpawnSign(1:ToAnnihilateIndex))
        ENDIF

!We now need to regenerate sendcounts and disps
        sendcounts(1:nProcessors)=0
        do i=1,ToAnnihilateIndex
            IF(ProcessVec(i).gt.(nProcessors-1)) THEN
                CALL Stop_All("RotoAnnihilateMinor","Annihilation error")
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

        CALL MPI_Barrier(MPI_COMM_WORLD,error)

!Perform another matrix transpose of the annihilation data using MPI_AlltoAllv, to send the data back to its correct Processor
!The signs of the hashes, index and CPU also need to be taken with them. (CPU does not need to be taken - every element of CPU should be equal to the rank of the processor+1)
!Hash also does not need to be taken, but will be taken as a precaution
        CALL MPI_AlltoAllv(IndexTable(1:ToAnnihilateonProc),sendcounts,disps,MPI_INTEGER,Index2Table,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)


!TEST
!        do i=1,ToAnnihilateonProc
!            IF(Process2Vec(i).ne.(iProcIndex)) THEN
!                CALL Stop_All("AnnihilateAmongstMinorSpawned","AlltoAllv performed incorrectly")
!            ENDIF
!        enddo

!Index2Table now is a list, of length "ToAnnihilateonProc", of walkers which should NOT be transferred to the next array. 
!Order the list according to this index (Hash and sign does not need to be sorted, but will for debugging purposes)
        CALL SORTIILongI(ToAnnihilateonProc,Index2Table(1:ToAnnihilateonProc),Hash2Array(1:ToAnnihilateonProc),MinorSpawnSign(1:ToAnnihilateonProc))

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
!                CALL Stop_All("AnnihilateAmongstMinorSpawned","Incorrect Hash returned")
!            ENDIF
!            IF(TempSign(Index2Table(i))) THEN
!                IF(.not.CurrentSign(i)) THEN
!                    CALL Stop_All("AnnihilateAmongstMinorSpawned","Incorrect Sign returned")
!                ENDIF
!            ELSE
!                IF(CurrentSign(i)) THEN
!                    CALL Stop_All("AnnihilateAmongstMinorSpawned","Incorrect Sign returned")
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
                    MinorSpawnDets(:,VecSlot)=MinorSpawnDets(:,i)
                    MinorSpawnParent(:,VecSlot)=MinorSpawnParent(:,i)
                    MinorSpawnSign(VecSlot)=TempSign(i)
                    i=i+1
                    VecSlot=VecSlot+1
                enddo
                i=i+1
            enddo

!Now need to copy accross the residual - from Index2Table(ToAnnihilateonProc) to TotWalkersNew
            do i=Index2Table(ToAnnihilateonProc)+1,MinorValidSpawned
                MinorSpawnDets(:,VecSlot)=MinorSpawnDets(:,i)
                MinorSpawnParent(:,VecSlot)=MinorSpawnParent(:,i)
                MinorSpawnSign(VecSlot)=TempSign(i)
                VecSlot=VecSlot+1
            enddo

        ELSE
!No particles annihilated
            VecSlot=1
            do i=1,MinorValidSpawned
                MinorSpawnDets(:,VecSlot)=MinorSpawnDets(:,i)
                MinorSpawnParent(:,VecSlot)=MinorSpawnParent(:,i)
                MinorSpawnSign(VecSlot)=TempSign(i)
                VecSlot=VecSlot+1
            enddo
        ENDIF
                
        MinorValidSpawnedNew=VecSlot-1

!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "FINAL CONFIGURATION: "
!            do i=1,TotWalkers
!                WRITE(6,*) i,HashArray(i),CurrentSign(i)
!            enddo
!        ENDIF

        IF((MinorValidSpawned-MinorValidSpawnedNew).ne.ToAnnihilateonProc) THEN
            WRITE(6,*) 'MinorValidSpawnedNew,MinorValidSpawned,ToAnnihilateonProc,Iter'
            WRITE(6,*) MinorValidSpawnedNew,MinorValidSpawned,ToAnnihilateonProc,Iter
            CALL FLUSH(6)
            CALL Stop_All("AnnihilateAmongstMinorSpawned","Problem with numbers when annihilating")
        ENDIF

        MinorValidSpawned=MinorValidSpawnedNew

        ! Don't need these after this, so rather than copying them back in the right order, re-zero to be ready for the 
        ! next set of spawned walkers on the minor determinants.
        HashArray(:)=0
        Hash2Array(:)=0


        DEALLOCATE(TempSign)
        DEALLOCATE(TempMinorSpawnSign)
        


    END SUBROUTINE AnnihilateAmongstMinorSpawned



    SUBROUTINE RotateMinorParticles(MinorValidSpawned)
        INTEGER :: i,MinorValidSpawned,error,Stat(MPI_STATUS_SIZE)

! This is the number of particles spawned (and still alive).  Must be sent with the arrays so the next processor knows the size.        
        MinorSpawnSign(0)=MinorValidSpawned

!Send the signs of the particles (number sent is in the first element)
        CALL MPI_BSend(MinorSpawnSign(0:MinorValidSpawned),MinorValidSpawned+1,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),123,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in sending signs")
        ENDIF

!...then send the particles themselves...
        CALL MPI_BSend(MinorSpawnDets(0:NIfD,1:MinorValidSpawned),MinorValidSpawned*(NIfD+1),MPI_INTEGER,MOD(iProcIndex+1,nProcessors),456,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in sending particles")
        ENDIF


!...and then send the parents of the walkers...
        CALL MPI_BSend(MinorSpawnParent(0:NIfD,1:MinorValidSpawned),MinorValidSpawned*(NIfD+1),MPI_INTEGER,MOD(iProcIndex+1,nProcessors),789,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in sending particle parents")
        ENDIF

        CALL MPI_Barrier(MPI_COMM_WORLD,error)

!Receive signs (let it receive the maximum possible (only the first ValidSpawned will be updated.))
        CALL MPI_Recv(MinorSpawnSign2(0:MaxSpawned),MaxSpawned+1,MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),123,MPI_COMM_WORLD,Stat,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in receiving signs")
        ENDIF

!Update the ValidSpawned variable for this new set of data we are about to receive...
        MinorValidSpawned=MinorSpawnSign2(0)

        CALL MPI_Recv(MinorSpawnDets2(0:NIfD,1:MinorValidSpawned),MinorValidSpawned*(NIfD+1),MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),456,MPI_COMM_WORLD,Stat,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in receiving particles")
        ENDIF

        CALL MPI_Recv(MinorSpawnParent2(0:NIfD,1:MinorValidSpawned),MinorValidSpawned*(NIfD+1),MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),789,MPI_COMM_WORLD,Stat,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in receiving particle parents")
        ENDIF

        do i=1,MinorValidSpawned
            MinorSpawnDets(0:NIfD,i)=MinorSpawnDets2(0:NIfD,i)
            MinorSpawnParent(0:NIfD,i)=MinorSpawnParent2(0:NIfD,i)
            MinorSpawnSign(i)=MinorSpawnSign2(i)
        enddo

!Really need to fix this so that I'm using pointers at some stage...

!We now want to make sure that we are working on the correct array. We have now received particles in SpawnedParts2 - switch it so that we are pointing at the other array.
!We always want to annihilate from the SpawedParts and SpawnedSign arrays.
!        IF(associated(SpawnedParts2,target=SpawnVec2)) THEN
!            SpawnedParts2 => SpawnVec
!            SpawnedSign2 => SpawnSignVec
!            SpawnedParts => SpawnVec2
!            SpawnedSign => SpawnSignVec2
!        ELSE
!            SpawnedParts => SpawnVec
!            SpawnedSign => SpawnSignVec
!            SpawnedParts2 => SpawnVec2
!            SpawnedSign2 => SpawnSignVec2
!        ENDIF



    END SUBROUTINE RotateMinorParticles



    SUBROUTINE AnnihilateMinorSpawnedParts(MinorValidSpawned,NoMinorWalkersNew)
        INTEGER :: i,j,k,MinorValidSpawned,NoMinorWalkersNew,MinInd,MaxDetInd,MinDetInd,SumDetPop,SumMinorDetPop
        INTEGER :: ToRemove,DetsMerged,SignProd,PartInd,EqDetPopsTag,ierr,FinalMinorDet
        LOGICAL :: DetsEq,tSuccess,tAnnihilateOne,DetBitEQ
        REAL*8 :: r,Prob
        INTEGER , ALLOCATABLE :: EqDetPops(:)

        MinInd=1
        ! This is the minimum index to start the search.  We are running through MinorSpawnDets (which is ordered), to find the matching det in 
        ! MinorStarDets (which is also ordered).  So if we find one det at a particular position, we only need to search MinorStarDets at positions
        ! lower than this.  But we start at 1.
        
        ToRemove=0 
        ! This is the number of particles to annihilate.

! Run through the newly spawned walkers
!        WRITE(6,*) 'MinorSpawedDets'
!        do j=1,MinorValidSpawned
!            WRITE(6,*) MinorSpawnDets(:,j),MinorSpawnSign(j)
!        enddo

        i=1
        do while (i.le.MinorValidSpawned)
            DetsEq=.false.
            SumMinorDetPop=MinorSpawnSign(i)
            j=1
            IF((i+j).le.MinorValidSpawned) DetsEq=DetBitEQ(MinorSpawnDets(0:NIfD,i),MinorSpawnDets(0:NIfD,i+j),NIfD)
            do while (DetsEq)
                SumMinorDetPop=SumMinorDetPop+MinorSpawnSign(i+j)
                j=j+1
                IF((i+j).gt.MinorValidSpawned) EXIT
                DetsEq=DetBitEQ(MinorSpawnDets(0:NIfD,i),MinorSpawnDets(0:NIfD,i+j),NIfD)
            enddo
            FinalMinorDet=i+j-1
            IF(FinalMinorDet.gt.MinorValidSpawned) FinalMinorDet=MinorValidSpawned
            ! The current spawned determinants therefore run from i to FinalMinorDet.
            ! These have an overall population of SumMinorDetPop.

!            IF((FinalMinorDet-i).gt.0) THEN
!                WRITE(6,*) 'MinorValidSpawned',MinorValidSpawned
!                WRITE(6,*) 'Starting Determinant',MinorSpawnDets(:,i)
!                WRITE(6,*) 'Determinant',MinorSpawnDets(:,FinalMinorDet)
!                WRITE(6,*) 'i',i
!                WRITE(6,*) 'FinalMinorDet',FinalMinorDet
!                WRITE(6,*) 'SumMinorDetPop',SumMinorDetPop
!                CALL FLUSH(6)
!            ENDIF

! Search for the determinant in the MinorStarDets list.
! tSuccess is true if the particle is found.
! This routine takes the MinorSpawnDets given and searches through MinorStarDets between MinInd and TotWalkersNew to find a match.  The index of this match is
! PartInd.
! In this case we need to check the determinants before and after the one found, to see if these are also equal.
            CALL BinSearchMinorParts(MinorSpawnDets(:,i),MinInd,NoMinorWalkersNew,PartInd,tSuccess)

            IF(tSuccess) THEN
!                WRITE(6,*) 'tSuccess'
                CALL FLUSH(6)
                ! Need to run forwards and backwards in the list of MinorStarDets, finding all the determinants that are equal, summing the particles on this determinant
                ! to find out how many need to be annihilated.  All the equal determinants in MinorStarDets should be the same sign, otherwise they will have annihilated already.

                ! Find out how many walkers are already on this determinant.
                SumDetPop=MinorStarSign(PartInd)
                DetsEq=.false.
                MinDetInd=PartInd
                MaxDetInd=PartInd

                ! First check one below the determinant found.
                j=1
                DetsEq=DetBitEQ(MinorSpawnDets(0:NIfD,i),MinorStarDets(0:NIfD,PartInd-j),NIfD)
                do while (DetsEq)
                    ! If the determinant is still equal, add the walkers on it to SumDetPop, and this index becomes the minimum.
                    SumDetPop=SumDetPop+MinorStarSign(PartInd-j)
                    MinDetInd=PartInd-j
                    j=j+1
                    DetsEq=DetBitEQ(MinorSpawnDets(0:NIfD,i),MinorStarDets(0:NIfD,PartInd-j),NIfD)
                    ! If this is true, the walkers on the next determinant will be added.
                enddo

                ! Now check those above the determinant found.
                j=1
                DetsEq=DetBitEQ(MinorSpawnDets(0:NIfD,i),MinorStarDets(0:NIfD,PartInd+j),NIfD)
                do while (DetsEq)
                    ! If the determinant is still equal, add the walkers on it to SumDetPop, and this index becomes the minimum.
                    SumDetPop=SumDetPop+MinorStarSign(PartInd+j)
                    MaxDetInd=PartInd+j
                    j=j+1
                    DetsEq=DetBitEQ(MinorSpawnDets(0:NIfD,i),MinorStarDets(0:NIfD,PartInd+j),NIfD)
                    ! If this is true, the walkers on the next determinant will be added.
                enddo
                ! SumDetPop now gives the number of walkers (with sign) currently on this determinant, and the Min and Max index of where these lie in MinorStarDets.  

!                SignProd=SumDetPop*MinorSpawnSign(i)
                SignProd=SumDetPop*SumMinorDetPop

!                IF((FinalMinorDet-i).gt.0) THEN
!                    WRITE(6,*) '*** Star stuff'
!                    WRITE(6,*) 'MinDetInd',MinDetInd
!                    WRITE(6,*) 'MaxDetInd',MaxDetInd
!                    WRITE(6,*) 'Determinant',MinorStarDets(:,MinDetInd)
!                    WRITE(6,*) 'SumDetPop',SumDetPop
!                ENDIF
                
                IF(SignProd.lt.0) THEN
                ! This suggests the spawned particles are of opposite sign to those currently on the determinant, and so must undergo annihilation.

                    IF((ABS(SumMinorDetPop)).ge.(ABS(SumDetPop))) THEN

!                        WRITE(6,*) 'In this bit 01'
!                        CALL FLUSH(6)

                        ! i.e. more (or equal) spawned than currently there, all walkers currently on that determinant are annihilated (regardless of parent), and the number
                        ! spawned is accordingly reduced.

                        ! Need to figure out which of the MinorSpawnSigns to annihilate.
!                        MinorSpawnSign(i)=MinorSpawnSign(i)+SumDetPop !!!!!!!!!!!
                        MinorAnnihilated=MinorAnnihilated+2*(ABS(SumDetPop))
                        Annihilated=Annihilated+2*(ABS(SumDetPop))

                        ALLOCATE(EqDetPops(i:FinalMinorDet),stat=ierr)
                        CALL LogMemAlloc('EqDetPops',FinalMinorDet-i+1,4,'AnnihilateMinorSpawnedParts',EqDetPopsTag,ierr)
                        IF(ierr.ne.0) CALL Stop_All('AnnihilateMinorSpawnedParts','Error allocating memory for EqDetPops')

                        do j=i,FinalMinorDet
                            EqDetPops(j)=ABS(MinorSpawnSign(j))
                        enddo

                        ! run through each walker on MinorSpawnSign, annihilating those in MinorStarSign one by one randomly.
                        ! the probability is the population in a particular entry of MinorStarSign / the total population from MinorStarSign on that determinant.
                        ! for each walker that annihilates, a random number is called, and based on these probabilities (which are calculated from the initial
                        ! populations before this annihilation i.e. the probabilities do not change as a walker is annihilated) a walker is annihilated from on of the entries.
!                        do j=1,ABS(SumDetPop)
                        j=1
                        do while (j.le.ABS(SumDetPop))
!                            WRITE(6,*) 'in this loop'
!                            CALL FLUSH(6)

                            ! call a random number
                            IF(tMerTwist) THEN
                                CALL genrand_real2(r) 
                            ELSE
                                CALL RANLUX(r,1)
                            ENDIF

                            tAnnihilateOne=.false.
                            
                            do while(.not.tAnnihilateOne)
                                ! tAnnihilateOne becomes true when a particle is annihilated, otherwise need to run through the probabilities again with a different random number.
!                                WRITE(6,*) 'in this loop 02'
!                                CALL FLUSH(6)
                                
                                Prob=0.D0
                                
                                do k=i,FinalMinorDet
!                                    WRITE(6,*) 'in this loop 03'
!                                    CALL FLUSH(6)
!                                    WRITE(6,*) 'i',i
!                                    WRITE(6,*)'FinalMinorDet',FinalMinorDet
!                                    WRITE(6,*) 'EqDetPops',EqDetPops(k)
!                                    WRITE(6,*) 'MinorSpawnSign',MinorSpawnSign(k)
!                                    WRITE(6,*) 'SumMinorDetPop',SumMinorDetPop

                                    Prob=Prob+ABS(REAL(EqDetPops(k),r2)/REAL(SumMinorDetPop,r2))
!                                    WRITE(6,*) 'Prob',Prob
!                                    WRITE(6,*) 'r',r

                                    IF(r.le.Prob) THEN
                                        IF(MinorSpawnSign(k).gt.0) THEN
                                            MinorSpawnSign(k)=MinorSpawnSign(k)-1
                                            tAnnihilateOne=.true.
                                        ELSEIF(MinorSpawnSign(k).lt.0) THEN
                                            MinorSpawnSign(k)=MinorSpawnSign(k)+1
                                            tAnnihilateOne=.true.
                                        ELSEIF(MinorSpawnSign(k).eq.0) THEN
                                            IF(tMerTwist) THEN
                                                CALL genrand_real2(r) 
                                            ELSE
                                                CALL RANLUX(r,1)
                                            ENDIF
                                            tAnnihilateOne=.false.
                                        ENDIF 
                                        EXIT
                                    ENDIF
                                enddo
                            enddo
                            j=j+1
                        enddo
                        DEALLOCATE(EqDetPops)
                        CALL LogMemDealloc('AnnihilateMinorSpawnedParts',EqDetPopsTag)

!                        WRITE(6,*) 'this o.k'
!                        CALL FLUSH(6)
 
                        do j=MinDetInd,MaxDetInd
                            MinorStarSign(j)=0
                        enddo

!                        WRITE(6,*) 'this o.k too'
!                        CALL FLUSH(6)
 
                        do j=i,FinalMinorDet
                            IF(MinorSpawnSign(j).eq.0) ToRemove=ToRemove+1
                            ! All particles have annihilated each other, there is none left in the spawned array so this determinant can be removed.
                        enddo
                        ! remaining walkers in MinorSpawnSign are not transferred across to MinorStarSign yet, as they need to be rotated, to test for other possible annihilations.

!                        WRITE(6,*) 'this o.k three'
!                        CALL FLUSH(6)

                    ELSE

!                        WRITE(6,*) 'In this bit 02'
!                        CALL FLUSH(6)


                        ! if there are less spawned than are currently on this determinant, the spawned annihilate some but not all.  need to randomly choose which to annihilate.
                        Annihilated=Annihilated+2*(ABS(SumMinorDetPop))
                        MinorAnnihilated=MinorAnnihilated+2*(ABS(SumMinorDetPop))

                        ALLOCATE(EqDetPops(MinDetInd:MaxDetInd),stat=ierr)
                        CALL LogMemAlloc('EqDetPops',MaxDetInd-MinDetInd+1,4,'AnnihilateMinorSpawnedParts',EqDetPopsTag,ierr)
                        IF(ierr.ne.0) CALL Stop_All('AnnihilateMinorSpawnedParts','Error allocating memory for EqDetPops')

                        do j=MinDetInd,MaxDetInd
                            EqDetPops(j)=ABS(MinorStarSign(j))
                        enddo

                        ! run through each walker on MinorSpawnSign, annihilating those in MinorStarSign one by one randomly.
                        ! the probability is the population in a particular entry of MinorStarSign / the total population from MinorStarSign on that determinant.
                        ! for each walker that annihilates, a random number is called, and based on these probabilities (which are calculated from the initial
                        ! populations before this annihilation i.e. the probabilities do not change as a walker is annihilated) a walker is annihilated from on of the entries.
                        j=1
!                        do j=1,ABS(SumMinorDetPop)
                        do while (j.le.ABS(SumMinorDetPop))

                            ! call a random number
                            IF(tMerTwist) THEN
                                CALL genrand_real2(r) 
                            ELSE
                                CALL RANLUX(r,1)
                            ENDIF

                            tAnnihilateOne=.false.
                            
                            do while(.not.tAnnihilateOne)
                                ! tAnnihilateOne becomes true when a particle is annihilated, otherwise need to run through the probabilities again with a different random number.
                                
                                Prob=0.D0
                                
                                do k=MinDetInd,MaxDetInd
                                    Prob=Prob+REAL(EqDetPops(k),r2)/ABS(REAL(SumDetPop,r2))
                                    IF(r.le.Prob) THEN
                                        IF(MinorStarSign(k).gt.0) THEN
                                            MinorStarSign(k)=MinorStarSign(k)-1
                                            tAnnihilateOne=.true.
                                        ELSEIF(MinorStarSign(k).lt.0) THEN
                                            MinorStarSign(k)=MinorStarSign(k)+1
                                            tAnnihilateOne=.true.
                                        ELSEIF(MinorStarSign(k).eq.0) THEN
                                            IF(tMerTwist) THEN
                                                CALL genrand_real2(r) 
                                            ELSE
                                                CALL RANLUX(r,1)
                                            ENDIF
                                            tAnnihilateOne=.false.
                                        ENDIF 
                                        EXIT
                                    ENDIF
                                enddo
                            enddo
                            j=j+1
                        enddo
                        DEALLOCATE(EqDetPops)
                        CALL LogMemDealloc('AnnihilateMinorSpawnedParts',EqDetPopsTag)
 
                        do j=i,FinalMinorDet
                            MinorSpawnSign(j)=0
                            ToRemove=ToRemove+1
                        enddo
                        
                    ENDIF

                ELSEIF(SignProd.gt.0) THEN
                    ! This means that the particle has found other particles on the same determinant with the same sign, therefore it cannot annihilate (as all other 
                    ! walkers on the same sign must be sign-coherent).  Therefore it can be just transferred across now.

                    ! These walkers however, must be added to MinorStarSign walkers with the same parent as these spawned ones.
                    ! I.e run over all the parents of all entries in MinorStarSign with the same determinants, until one is found that is the same as the parent of the spawned.


!                    WRITE(6,*) 'In this bit 03'
!                    CALL FLUSH(6)


                    DetsEq=.false.
                    do j=MinDetInd,MaxDetInd
                        DetsEq=DetBitEQ(MinorSpawnParent(0:NIfD,i),MinorStarParent(0:NIfD,j),NIfD)
                        IF(DetsEq) THEN
                            MinorStarSign(j)=MinorStarSign(j)+MinorSpawnSign(i)
                            MinorSpawnSign(i)=0
                            ToRemove=ToRemove+1
                            EXIT
                        ENDIF
                    enddo
!                    IF(.not.DetsEq) THEN
                        ! This just means the determinant has been spawned on from a different parent.
                        ! Leave this in the spawned list - it will be quicker to just merge them all at once, rather than merging now.
!                        WRITE(6,*) 'determinant then parent of star then spawn'
!                        WRITE(6,*) MinorStarDets(0:NIfD,i),'*',MinorStarParent(0:NIfD,i)
!                        do j=MinDetInd,MaxDetInd
!                            WRITE(6,*) MinorSpawnDets(0:NIfD,j),'*',MinorSpawnParent(0:NIfD,j)
!                        enddo
!                        CALL FLUSH(6)
!                        CALL Stop_All('AnnihilateMinorSpawnedParts','Error adding sign coherent spawned particles to the list of current determinants.')
!                    ENDIF
                ELSEIF(SignProd.eq.0) THEN

                    IF(MinorSpawnSign(i).eq.0) ToRemove=ToRemove+1
                ENDIF
                ! This ENDIF means we have dealt with all the cases where the spawned determinant is found in the MinorStarDets list and the signs are the same/different etc.
                
            ENDIF
            ! If the spawned determinant isn't in the MinorStarDets list there isn't anything else to do.

            MinInd=MaxDetInd
            i=FinalMinorDet+1
        enddo
        ! Do this for all spawned on determinants.

!        WRITE(6,*) 'here o.k'
!        CALL FLUSH(6)

! Now remove all the annihilated particles from the spawned list.  I.e. those which now have 0 particles on that determinant, do not need to be in the list.
        IF(ToRemove.gt.0) THEN
            DetsMerged=0
            do i=1,MinorValidSpawned
                IF(MinorSpawnSign(i).eq.0) THEN
                    DetsMerged=DetsMerged+1
                ELSE
                    MinorSpawnDets2(0:NIfD,i-DetsMerged)=MinorSpawnDets(0:NIfD,i)
                    MinorSpawnSign2(i-DetsMerged)=MinorSpawnSign(i)
                    MinorSpawnParent2(0:NIfD,i-DetsMerged)=MinorSpawnParent(0:NIfD,i)
                ENDIF
            enddo
            MinorValidSpawned=MinorValidSpawned-DetsMerged
            IF(DetsMerged.ne.ToRemove) THEN
                WRITE(6,*) "*** Iteration number", Iter
                WRITE(6,*) 'DetsMerged,',DetsMerged,'ToRemove,',ToRemove
                CALL Stop_All("AnnihilateMinorSpawnedParts","Incorrect number of particles removed from minor spawned list")
            ENDIF

            ! My version of changing the pointers over, need to fix this.
            do i=1,MinorValidSpawned
                MinorSpawnDets(0:NIfD,i)=MinorSpawnDets2(0:NIfD,i)
                MinorSpawnParent(0:NIfD,i)=MinorSpawnParent2(0:NIfD,i)
                MinorSpawnSign(i)=MinorSpawnSign2(i)
            enddo
        ENDIF

!        WRITE(6,*) 'here o.k too'
!        CALL FLUSH(6)


    END SUBROUTINE AnnihilateMinorSpawnedParts


!This routine will run through the total list of minor particles (NoMinorWalkersNew in MinorStarDets with sign MinorStarSign) and the list of newly-spawned but
!surviving particles (MinorValidSpawned in MinorSpawnDets and MinorSpawnSign) and move the new particles into the correct place in the new list,
!while removing the particles with sign = 0 from MinorStarDets. 
!Binary searching can be used to speed up this transfer substantially.
!This needs to be modified slightly compared to InsertRemoveParts, as in this case, it is possible for the same determinant to be specified in both the
!spawned and main list, but these will have different parents, and thus must be kept separate.
    SUBROUTINE InsertRemoveMinorParts(MinorValidSpawned,NoMinorWalkersNew)
        INTEGER :: NoMinorWalkersNew,MinorValidSpawned
        INTEGER :: i,DetsMerged

        
! Remove determinants from the main array which have 0 population.        
        TotParts=0
        DetsMerged=0
        do i=1,NoMinorWalkersNew
            IF(MinorStarSign(i).eq.0) THEN
                DetsMerged=DetsMerged+1
            ELSE
! We want to move all the elements above this point down to 'fill in' the annihilated determinant.
                IF(DetsMerged.ne.0) THEN
                    MinorStarDets(0:NIfD,i-DetsMerged)=MinorStarDets(0:NIfD,i)
                    MinorStarSign(i-DetsMerged)=MinorStarSign(i)
                    MinorStarParent(0:NIfD,i-DetsMerged)=MinorStarParent(0:NIfD,i)
                    MinorStarHii(i-DetsMerged)=MinorStarHii(i)
                    MinorStarHij(i-DetsMerged)=MinorStarHij(i)
                ENDIF
                TotParts=TotParts+abs(MinorStarSign(i))
            ENDIF
        enddo
        NoMinorWalkersNew=NoMinorWalkersNew-DetsMerged
        ! So this is the number of determinants specified in the main list before those spawned and survived have been added.

!We now need to compress the spawned list, so that no particles are specified more than once.
!We also want to find the number of particles we are adding to the list from the spawned list.
!We now calculate the contribution to the total number of particles from the spawned lists.
!The list has previously been compressed before the annihilation began.
        IF(MinorValidSpawned.gt.0) THEN
            TotParts=TotParts+abs(MinorSpawnSign(1))
        ENDIF
        do i=2,MinorValidSpawned
            TotParts=TotParts+abs(MinorSpawnSign(i))
        enddo

!We now want to merge the main list with the spawned list of surviving spawned particles.
!The final list will be of length NoMinorWalkers+MinorValidSpawned. This will be returned in the first element of MergeLists updated.
        
        IF(TotParts.gt.0) THEN

            CALL MergeListswH2(NoMinorWalkersNew,MaxWalkersPart,MinorValidSpawned,MinorSpawnDets(0:NIfD,1:MinorValidSpawned),&
            &MinorSpawnParent(0:NIfD,1:MinorValidSpawned),MinorSpawnSign(1:MinorValidSpawned),NIfD)
            
        ENDIF

        NoMinorWalkers=NoMinorWalkersNew    


    END SUBROUTINE InsertRemoveMinorParts


    
    SUBROUTINE AnnihilateBetweenSpawnedOneProc(ValidSpawned)
        INTEGER :: ValidSpawned,DetCurr(0:NIfD),i,j,k,LowBound,HighBound,WSign
        INTEGER :: VecSlot,TotSign
        LOGICAL :: DetBitEQ

        CALL SortBitDets(ValidSpawned,SpawnedParts(0:NIfD,1:ValidSpawned),NIfD,SpawnedSign(1:ValidSpawned))

        VecSlot=1
        i=1
        do while(i.le.ValidSpawned)
            LowBound=i
            DetCurr(0:NIfD)=SpawnedParts(0:NIfD,i)
            i=i+1
            do while(DetBitEQ(DetCurr(0:NIfD),SpawnedParts(0:NIfD,i),NIfD).and.(i.le.ValidSpawned))
                i=i+1
            enddo
            HighBound=i-1

!Now, run through the block of common particles again, counting the residual sign
            TotSign=0
            do j=LowBound,HighBound
                TotSign=TotSign+SpawnedSign(j)
            enddo

!Now, fill up SpawnedSign2 and SpawnedParts2 with the residual particles
            IF(TotSign.ne.0) THEN
                WSign=INT(TotSign/abs(TotSign))
                do k=1,abs(TotSign)
                    SpawnedParts2(0:NIfD,VecSlot)=DetCurr(0:NIfD)
                    SpawnedSign2(VecSlot)=WSign
                    VecSlot=VecSlot+1
                enddo
            ENDIF

        enddo

        ValidSpawned=VecSlot-1

        do i=1,ValidSpawned
            SpawnedParts(0:NIfD,i)=SpawnedParts2(0:NIfD,i)
            SpawnedSign(i)=SpawnedSign2(i)
        enddo

    END SUBROUTINE AnnihilateBetweenSpawnedOneProc

!This routine will run through the total list of particles (TotWalkersNew in CurrentDets with sign CurrentSign) and the list of newly-spawned but
!non annihilated particles (ValidSpawned in SpawnedParts and SpawnedSign) and move the new particles into the correct place in the new list,
!while removing the particles with sign = 0 from CurrentDets. 
!Binary searching can be used to speed up this transfer substantially.
!The key feature which makes this work, is that it is impossible for the same determinant to be specified in both the spawned and main list at the end of
!the annihilation process. Therefore we will not multiply specify determinants when we merge the lists.
    SUBROUTINE InsertRemoveParts(ValidSpawned,TotWalkersNew)
        USE Determinants , only : GetHElement3
        INTEGER :: TotWalkersNew,ValidSpawned
        INTEGER :: i,DetsMerged,nJ(NEl)
        REAL*8 :: HDiag
        TYPE(HElement) :: HDiagTemp
        LOGICAL :: DetBitEQ

!        IF(Iter.eq.56) THEN
!            WRITE(6,*) "Merging lists, ",TotWalkersNew,Iter
!            do i=1,TotWalkersNew
!                WRITE(6,*) CurrentDets(:,i),CurrentSign(i)
!            enddo
!            WRITE(6,*) "Spawned list: ",ValidSpawned,Iter
!            do i=1,ValidSpawned
!                WRITE(6,*) SpawnedParts(:,i),SpawnedSign(i)
!            enddo
!            WRITE(6,*) "*****"
!            CALL FLUSH(6)
!        ENDIF
        
!If we want to do this while only keeping the data in one array, the first thing which is needed, is for the annihilated
!determinants to be removed from the main array. These are denoted by zeros in the sign array for it.
!Surely we only need to perform this loop if the number of annihilated particles > 0?

        TotParts=0
        DetsMerged=0
        do i=1,TotWalkersNew
            IF(CurrentSign(i).eq.0) THEN
                DetsMerged=DetsMerged+1
            ELSE
!We want to move all the elements above this point down to 'fill in' the annihilated determinant.
                IF(DetsMerged.ne.0) THEN
                    CurrentDets(0:NIfD,i-DetsMerged)=CurrentDets(0:NIfD,i)
                    CurrentSign(i-DetsMerged)=CurrentSign(i)
                    IF(.not.tRegenDiagHEls) THEN
                        CurrentH(i-DetsMerged)=CurrentH(i)
                    ENDIF
                ENDIF
                TotParts=TotParts+abs(CurrentSign(i))
            ENDIF
        enddo
        TotWalkersNew=TotWalkersNew-DetsMerged

!        do i=1,TotWalkersNew
!            IF(CurrentSign(i).eq.0) THEN
!                CALL Stop_All("InsertRemoveParts","Particles not removed correctly")
!            ENDIF
!        enddo
!        CALL CheckOrdering(CurrentDets,CurrentSign(1:TotWalkersNew),TotWalkersNew,.true.)

!We now need to compress the spawned list, so that no particles are specified more than once.
!We also want to find the number of particles we are adding to the list from the spawned list.
!We now calculate the contribution to the total number of particles from the spawned lists.
!The list has previously been compressed before the annihilation began.
        IF(ValidSpawned.gt.0) THEN
            TotParts=TotParts+abs(SpawnedSign(1))
        ENDIF
        do i=2,ValidSpawned
            TotParts=TotParts+abs(SpawnedSign(i))
        enddo

!        CALL CheckOrdering(SpawnedParts,SpawnedSign(1:ValidSpawned),ValidSpawned,.true.)
!
!        WRITE(6,*) "Before merging...",TotWalkersNew,Iter
!        do i=1,TotWalkersNew
!            WRITE(6,*) i,CurrentDets(:,i),CurrentSign(i)
!        enddo
!        WRITE(6,*) "***"
!        do i=1,ValidSpawned
!            WRITE(6,*) i,SpawnedParts(:,i),SpawnedSign(i)
!        enddo
!        WRITE(6,*) "***"
!        CALL FLUSH(6)

!TotWalkersNew is now the number of non-annihilated determinants in the main list left.
!We now want to merge the main list with the spawned list of non-annihilated spawned particles.
!The final list will be of length TotWalkersNew+ValidSpawned. This will be returned in the first element of MergeLists updated.
        
       
        IF(tRegenDiagHEls) THEN

            IF(TotWalkersNew.eq.0) THEN
!Merging algorithm will not work with no determinants in the main list.
                TotWalkersNew=ValidSpawned
                do i=1,ValidSpawned
                    CurrentDets(:,i)=SpawnedParts(:,i)
                    CurrentSign(i)=SpawnedSign(i)
                enddo
            ELSE
                CALL MergeLists(TotWalkersNew,MaxWalkersPart,ValidSpawned,SpawnedParts(0:NIfD,1:ValidSpawned),SpawnedSign(1:ValidSpawned),NIfD)
            ENDIF
        ELSE
            IF(TotWalkersNew.eq.0) THEN
!Merging algorithm will not work with no determinants in the main list.
                TotWalkersNew=ValidSpawned
                do i=1,ValidSpawned
                    CurrentDets(:,i)=SpawnedParts(:,i)
                    CurrentSign(i)=SpawnedSign(i)
!We want to calculate the diagonal hamiltonian matrix element for the new particle to be merged.
                    IF(DetBitEQ(CurrentDets(0:NIfD,i),iLutHF,NIfD)) THEN
!We know we are at HF - HDiag=0
                        HDiag=0.D0
                    ELSE
                        CALL DecodeBitDet(nJ,CurrentDets(0:NIfD,i),NEl,NIfD)
                        IF(tHPHF) THEN
                            CALL HPHFGetDiagHElement(nJ,CurrentDets(0:NIfD,i),HDiagTemp)
                        ELSE
                            HDiagTemp=GetHElement3(nJ,nJ,0)
                        ENDIF
                        HDiag=(REAL(HDiagTemp%v,8))-Hii
                    ENDIF
                    CurrentH(i)=HDiag
                enddo
            ELSE
                CALL MergeListswH(TotWalkersNew,MaxWalkersPart,ValidSpawned,SpawnedParts(0:NIfD,1:ValidSpawned),SpawnedSign(1:ValidSpawned),NIfD)
            ENDIF

        ENDIF
        TotWalkers=TotWalkersNew

!        CALL CheckOrdering(CurrentDets,CurrentSign(1:TotWalkers),TotWalkers,.true.)

!There is now no need to swap the pointers around since we only have one array
!        IF(associated(CurrentDets,target=WalkVecDets)) THEN
!            CurrentDets=>WalkVec2Dets
!            CurrentSign=>WalkVec2Sign
!            CurrentH=>WalkVec2H
!            NewDets=>WalkVecDets
!            NewSign=>WalkVecSign
!            NewH=>WalkVecH
!        ELSE
!            CurrentDets=>WalkVecDets
!            CurrentSign=>WalkVecSign
!            CurrentH=>WalkVecH
!            NewDets=>WalkVec2Dets
!            NewSign=>WalkVec2Sign
!            NewH=>WalkVec2H
!        ENDIF
            
!        WRITE(6,*) "Final Merged List: "
!        do i=1,TotWalkers
!            WRITE(6,*) CurrentDets(:,i),CurrentSign(i)
!        enddo

!Below is a seperate way of doing this, which transfers the particles to a new array.

!        IndSpawned=1
!        IndParts=1
!        VecInd=1
!        TotParts=0
!
!        IF((ValidSpawned.gt.0).and.(TotWalkersNew.gt.0)) THEN
!            do while(IndParts.le.TotWalkersNew)
!                
!                CompParts=DetBitLT(NewDets(0:NIfD,IndParts),SpawnedParts(0:NIfD,IndSpawned),NIfD)
!                IF(CompParts.eq.1) THEN
!!Want to move in the particle from NewDets (unless it wants to be annihilated)
!                    IF(NewSign(IndParts).ne.0) THEN
!!We want to keep this particle
!                        CurrentDets(0:NIfD,VecInd)=NewDets(0:NIfD,IndParts)
!                        CurrentSign(VecInd)=NewSign(IndParts)
!                        IF(.not.tRegenDiagHEls) CurrentH(VecInd)=NewH(IndParts)
!                        VecInd=VecInd+1
!                        TotParts=TotParts+abs(NewSign(IndParts))
!
!                    ENDIF
!                    IndParts=IndParts+1
!!                ELSEIF(IndSpawned.le.ValidSpawned) THEN
!                ELSEIF(CompParts.eq.0) THEN
!!This should be taken out later - the lists will be disjoint.
!!This will add the particles on the same determinant together...
!
!                    CurrentDets(0:NIfD,VecInd)=NewDets(0:NIfD,IndParts)
!                    IF(.not.tRegenDiagHEls) CurrentH(VecInd)=NewH(IndParts)
!                    CurrentSign(VecInd)=NewSign(IndParts)+SpawnedSign(IndSpawned)
!                    IndParts=IndParts+1
!                    IndSpawned=IndSpawned+1
!
!                    do while(DetBitEQ(SpawnedParts(:,IndSpawned-1),SpawnedParts(:,IndSpawned),NIfD).and.(IndSpawned.le.ValidSpawned))
!                        CurrentSign(VecInd)=CurrentSign(VecInd)+SpawnedSign(IndSpawned)
!                        IndSpawned=IndSpawned+1
!                    enddo
!
!                    TotParts=TotParts+abs(CurrentSign(VecInd))
!                    VecInd=VecInd+1
!                    IF(IndSpawned.gt.ValidSpawned) THEN
!                        IndSpawned=IndSpawned+1
!                        EXIT    !We have reached the end of the list of spawned particles
!                    ENDIF
!
!                ELSE
!!Now, we want to transfer a spawned particle, unless we have transferred them all
!                    IF(SpawnedSign(IndSpawned).eq.0) THEN
!                        CALL Stop_All("InsertRemoveParts","Should not have particles marked for annihilation in this array")
!                    ENDIF
!                    CurrentDets(0:NIfD,VecInd)=SpawnedParts(0:NIfD,IndSpawned)
!                    CurrentSign(VecInd)=SpawnedSign(IndSpawned)
!                    IndSpawned=IndSpawned+1
!                    
!                    do while(DetBitEQ(SpawnedParts(:,IndSpawned-1),SpawnedParts(:,IndSpawned),NIfD).and.(IndSpawned.le.ValidSpawned))
!                        CurrentSign(VecInd)=CurrentSign(VecInd)+SpawnedSign(IndSpawned)
!                        IndSpawned=IndSpawned+1
!                    enddo
!                    
!                    TotParts=TotParts+abs(CurrentSign(VecInd))
!
!                    IF(.not.tRegenDiagHEls) THEN
!!Need to find H-element!
!                        IF(DetBitEQ(CurrentDets(0:NIfD,VecInd),iLutHF,NIfD)) THEN
!!We know we are at HF - HDiag=0
!                            HDiag=0.D0
!                            IF(tHub.and.tReal) THEN
!!Reference determinant is not HF
!                                CALL DecodeBitDet(nJ,CurrentDets(0:NIfD,VecInd),NEl,NIfD)
!                                HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                                HDiag=(REAL(HDiagTemp%v,r2))
!                            ENDIF
!                        ELSE
!                            CALL DecodeBitDet(nJ,CurrentDets(0:NIfD,VecInd),NEl,NIfD)
!                            HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                            HDiag=(REAL(HDiagTemp%v,r2))-Hii
!                        ENDIF
!                        CurrentH(VecInd)=HDiag
!
!                    ENDIF
!                    
!                    VecInd=VecInd+1
!                    IF(IndSpawned.gt.ValidSpawned) THEN
!                        IndSpawned=IndSpawned+1
!                        EXIT    !We have reached the end of the list of spawned particles
!                    ENDIF
!                ENDIF
!
!            enddo
!        ENDIF
!
!        IF(IndParts.le.TotWalkersNew) THEN
!!Haven't finished copying rest of original particles
!            do i=IndParts,TotWalkersNew
!                IF(NewSign(i).ne.0) THEN
!                    CurrentDets(0:NIfD,VecInd)=NewDets(0:NIfD,i)
!                    CurrentSign(VecInd)=NewSign(i)
!                    IF(.not.tRegenDiagHEls) CurrentH(VecInd)=NewH(i)
!                    TotParts=TotParts+abs(NewSign(i))
!                    VecInd=VecInd+1
!                ENDIF
!            enddo
!
!        ELSEIF(IndSpawned.le.ValidSpawned) THEN
!            do while(IndSpawned.le.ValidSpawned)
!
!                CurrentDets(0:NIfD,VecInd)=SpawnedParts(0:NIfD,IndSpawned)
!                CurrentSign(VecInd)=SpawnedSign(IndSpawned)
!                IndSpawned=IndSpawned+1
!
!                do while(DetBitEQ(SpawnedParts(:,IndSpawned-1),SpawnedParts(:,IndSpawned),NIfD).and.(IndSpawned.le.ValidSpawned))
!                    CurrentSign(VecInd)=CurrentSign(VecInd)+SpawnedSign(IndSpawned)
!                    IndSpawned=IndSpawned+1
!                enddo
!                    
!                TotParts=TotParts+abs(CurrentSign(VecInd))
!
!                IF(.not.tRegenDiagHEls) THEN
!!Need to find H-element!
!                    IF(DetBitEQ(CurrentDets(0:NIfD,VecInd),iLutHF,NIfD)) THEN
!!We know we are at HF - HDiag=0
!                        HDiag=0.D0
!                        IF(tHub.and.tReal) THEN
!!Reference determinant is not HF
!                            CALL DecodeBitDet(nJ,CurrentDets(0:NIfD,VecInd),NEl,NIfD)
!                            HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                            HDiag=(REAL(HDiagTemp%v,r2))
!                        ENDIF
!                    ELSE
!                        CALL DecodeBitDet(nJ,CurrentDets(0:NIfD,VecInd),NEl,NIfD)
!                        HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                        HDiag=(REAL(HDiagTemp%v,r2))-Hii
!                    ENDIF
!                    CurrentH(VecInd)=HDiag
!
!                ENDIF
!
!                VecInd=VecInd+1
!            enddo
!        ENDIF

!TotParts is now the total number of particles in the system.
!This should match the decrease in size of the SpawnedParts array over the course of the annihilation steps.
!This should also be equal to TotWalkersNew-TotWalkers

!        TotWalkers=VecInd-1     !The new total number of particles after all annihilation steps.
!        Annihilated=Annihilated+(InitialSpawned-ValidSpawned)+OrigPartAnn   !The total number of annihilated particles is simply the number annihilated from spawned
                                                                !list plus the number annihilated from the original list.

!        IF(TotWalkers.ne.(InitialSpawned+TotWalkersNew-Annihilated)) THEN
!            CALL Stop_All("InsertRemoveParts","Error in number of surviving particles")
!        ENDIF

!        AnnFromSpawned=SpawnedBeforeRoto-ValidSpawned
!        TotAnnFromSpawned=0
!        TotAnnFromOrig=0
!        CALL MPI_Reduce(OrigPartAnn,TotAnnFromOrig,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(AnnFromSpawned,TotAnnFromSpawned,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)

!        IF(TotAnnFromOrig.ne.TotAnnFromSpawned) THEN
!            WRITE(6,*) Iter,TotAnnFromOrig,TotAnnFromSpawned,SpawnedBeforeRoto,ValidSpawned
!            CALL Stop_All("InsertRemoveParts","Different numbers of particles annihilated from spawned and original lists.")
!        ENDIF

!        WRITE(6,*) "Annihilated: ",Annihilated,InitialSpawned,ValidSpawned,OrigPartAnn

!        IF(Iter.eq.56) THEN
!            do i=1,VecInd-1
!                WRITE(6,*) i,CurrentDets(:,i),CurrentSign(i)
!            enddo
!        ENDIF

    
    END SUBROUTINE InsertRemoveParts

!This routine wants to take the ValidSpawned particles in the SpawnedParts array and perform All-to-All communication so that 
!we can annihilate all common particles with opposite signs.
!Particles are fed in on the SpawnedParts and SpawnedSign array, and are returned in the same arrays.
!It requires MaxSpawned*36 bytes of memory (on top of the memory of the arrays fed in...)
!Might not need to send hashes in all-to-all - could just use them for determining where they go
!Package up temp arrays?
    SUBROUTINE AnnihilateBetweenSpawned(ValidSpawned)
        INTEGER(KIND=i2) , ALLOCATABLE :: HashArray1(:),HashArray2(:)
        INTEGER , ALLOCATABLE :: IndexTable1(:),IndexTable2(:),ProcessVec1(:),ProcessVec2(:),TempSign(:)
        INTEGER :: i,j,k,ToAnnihilateIndex,ValidSpawned,ierr,error,sendcounts(nProcessors)
        INTEGER :: TotWalkersDet,InitialBlockIndex,FinalBlockIndex,ToAnnihilateOnProc,VecSlot
        INTEGER :: disps(nProcessors),recvcounts(nProcessors),recvdisps(nProcessors),nJ(NEl)
        INTEGER :: Minsendcounts,Maxsendcounts,DebugIter,SubListInds(2,nProcessors),MinProc,MinInd
        INTEGER(KIND=i2) :: HashCurr,MinBin,RangeofBins,NextBinBound,MinHash
        CHARACTER(len=*), PARAMETER :: this_routine='AnnihilateBetweenSpawned'

        CALL set_timer(AnnSpawned_time,30)

!First, we need to allocate memory banks. Each array needs a hash value, a processor value, and an index value.
!We also want to allocate a temporary sign value
        ALLOCATE(TempSign(ValidSpawned),stat=ierr)

!These arrays may as well be kept all the way through the simulation?
        ALLOCATE(HashArray1(MaxSpawned),stat=ierr)
        ALLOCATE(HashArray2(MaxSpawned),stat=ierr)
        ALLOCATE(IndexTable1(MaxSpawned),stat=ierr)
        ALLOCATE(IndexTable2(MaxSpawned),stat=ierr)
        ALLOCATE(ProcessVec1(MaxSpawned),stat=ierr)
        ALLOCATE(ProcessVec2(MaxSpawned),stat=ierr)

        IF(ierr.ne.0) THEN
            CALL Stop_All("AnnihilateBetweenSpawned","Error in allocating initial data")
        ENDIF

        TempSign(1:ValidSpawned)=SpawnedSign(1:ValidSpawned)
        ProcessVec1(1:ValidSpawned)=iProcIndex

!        WRITE(6,*) "***************************************"
        do i=1,ValidSpawned
            IndexTable1(i)=i
            CALL DecodeBitDet(nJ,SpawnedParts(0:NIfD,i),NEl,NIfD)
            HashArray1(i)=CreateHash(nJ)
!            IF(Iter.eq.1346.and.(HashArray1(i).eq.2905380077198165348)) THEN
!                WRITE(6,*) "Hash found, ",i,SpawnedSign(i),HashArray1(i),SpawnedParts(0:NIfD,i)
!            ENDIF
        enddo

!Next, order the hash array, taking the index, CPU and sign with it...
        IF(.not.tAnnihilatebyRange) THEN
!Order the array by abs(mod(Hash,nProcessors)). This will result in a more load-balanced system (no need to actually take ProcessVec1 - this will always be iProcIndex here.
            CALL SortMod4ILong(ValidSpawned,HashArray1(1:ValidSpawned),IndexTable1(1:ValidSpawned),ProcessVec1(1:ValidSpawned),SpawnedSign(1:ValidSpawned),nProcessors)

!Send counts is the size of each block of ordered dets which are going to each processor. This could be binary searched for extra speed
            IF(ValidSpawned.gt.0) THEN
                j=1
                do i=0,nProcessors-1    !Search through all possible values of abs(mod(Hash,nProcessors))
                    do while((abs(mod(HashArray1(j),INT(nProcessors,8))).eq.i).and.(j.le.ValidSpawned))
                        j=j+1
                    enddo
                    sendcounts(i+1)=j-1
                enddo
            ELSE
                sendcounts(1:nProcessors)=0
            ENDIF

        ELSE
!We can try to sort the hashes by range, which may result in worse load-balancing, but will remove the need for a second sort of the hashes once they have been sent to the correct processor.
            CALL Sort4ILong(ValidSpawned,HashArray1(1:ValidSpawned),IndexTable1(1:ValidSpawned),ProcessVec1(1:ValidSpawned),SpawnedSign(1:ValidSpawned))
!We also need to know the ranges of the hashes to send to each processor. Each range should be the same.
            IF(nProcessors.ne.1) THEN
                Rangeofbins=INT(HUGE(Rangeofbins)/(nProcessors/2),8)
                MinBin=-HUGE(MinBin)
                NextBinBound=MinBin+Rangeofbins

!We need to find the indices for each block of hashes which are to be sent to each processor.
!Sendcounts is the size of each block of ordered dets which are going to each processors. This could be binary searched for extra speed.
                j=1
                do i=1,nProcessors    !Search through all possible values of the hashes
                    do while((HashArray1(j).le.NextBinBound).and.(j.le.ValidSpawned))
                        j=j+1
                    enddo
                    sendcounts(i)=j-1
                    IF(i.eq.nProcessors-1) THEN
!Make sure the final bin catches everything...
                        NextBinBound=HUGE(NextBinBound)
                    ELSE
                        NextBinBound=NextBinBound+Rangeofbins
                    ENDIF
                enddo
            ELSE
                sendcounts(1)=ValidSpawned
!                do j=1,ValidSpawned
!                    WRITE(6,*) Iter,j,HashArray1(j),SpawnedSign(j)
!                enddo
                    
            ENDIF
        ENDIF

        IF(sendcounts(nProcessors).ne.ValidSpawned) THEN
            WRITE(6,*) "SENDCOUNTS is: ",sendcounts(:)
            WRITE(6,*) "VALIDSPAWNED is: ",ValidSpawned
            CALL FLUSH(6)
            CALL Stop_All("AnnihilateBetweenSpawned","Incorrect calculation of sendcounts")
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

!We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
        recvcounts(1:nProcessors)=0

        CALL MPI_AlltoAll(sendcounts,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,error)

!We can now get recvdisps from recvcounts in the same way we obtained disps from sendcounts
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo

        MaxIndex=recvdisps(nProcessors)+recvcounts(nProcessors)
!Max index is the largest occupied index in the array of hashes to be ordered in each processor 
        IF(MaxIndex.gt.(0.93*MaxSpawned)) THEN
            CALL Warning("AnnihilateBetweenSpawned","Maximum index of annihilation array is close to maximum length. Increase MemoryFacSpawn")
        ENDIF

!Uncomment this if you want to write out load-balancing statistics.
!        AnnihilPart(:)=0
!        CALL MPI_Gather(MaxIndex,1,MPI_INTEGER,AnnihilPart,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
!        IF(iProcIndex.eq.root) THEN
!            WRITE(13,"(I10)",advance='no') Iter
!            do i=1,nProcessors
!                WRITE(13,"(I10)",advance='no') AnnihilPart(i)
!            enddo
!            WRITE(13,"(A)") ""
!            CALL FLUSH(13)
!        ENDIF

!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "RECVCOUNTS: "
!            WRITE(6,*) recvcounts(:)
!            WRITE(6,*) "RECVDISPS: "
!            WRITE(6,*) recvdisps(:),MaxIndex
!            CALL FLUSH(6)
!        ENDIF

!Insert a load-balance check here...maybe find the s.d. of the sendcounts array - maybe just check the range first.
!        IF(TotWalkersNew.gt.200) THEN
!            IF((Maxsendcounts-Minsendcounts).gt.(TotWalkersNew/3)) THEN
!                WRITE(6,"(A,I12)") "**WARNING** Parallel annihilation not optimally balanced on this node, for iter = ",Iter
!                WRITE(6,*) "Sendcounts is: ",sendcounts(:)
!                CALL FLUSH(6)
!            ENDIF
!        ENDIF

!Now send the chunks of hashes to the corresponding processors
        CALL MPI_AlltoAllv(HashArray1(1:ValidSpawned),sendcounts,disps,MPI_DOUBLE_PRECISION,HashArray2(1:MaxIndex),recvcounts,recvdisps,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,error)

!The signs of the hashes, index and CPU also need to be taken with them.
        CALL MPI_AlltoAllv(SpawnedSign(1:ValidSpawned),sendcounts,disps,MPI_INTEGER,SpawnedSign2(1:MaxIndex),recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        CALL MPI_AlltoAllv(IndexTable1(1:ValidSpawned),sendcounts,disps,MPI_INTEGER,IndexTable2,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        CALL MPI_AlltoAllv(ProcessVec1(1:ValidSpawned),sendcounts,disps,MPI_INTEGER,ProcessVec2,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)

        IF(.not.tAnnihilatebyrange) THEN
!The hashes now need to be sorted again - this time by their number
!This sorting would be redundant if we had initially sorted the hashes by range (ie tAnnihilatebyrange).
            CALL Sort4ILong(MaxIndex,HashArray2(1:MaxIndex),IndexTable2(1:MaxIndex),ProcessVec2(1:MaxIndex),SpawnedSign2(1:MaxIndex))
        ELSE
!Here, because we have ordered the hashes initially numerically, we have a set of ordered lists. It is therefore easier to sort them.
!We have to work out how to run sequentially through the hashes, which are a set of nProc seperate ordered lists.
!We would need to have 2*nProc indices, since we will have a set of nProc disjoint ordered sublists.
!SubListInds(1,iProc)=index of current hash from processor iProc
!SubListInds(2,iProc)=index of final hash from processor iProc
!Indices can be obtained from recvcounts and recvdisps - recvcounts(iProc-1) is number of hashes from iProc
!recvdisps(iProc-1) is the displacement to the start of the hashes from iProc
            do i=1,nProcessors-1
                SubListInds(1,i)=recvdisps(i)+1
                SubListInds(2,i)=recvdisps(i+1)
            enddo
            SubListInds(1,nProcessors)=recvdisps(nProcessors)+1
            SubListInds(2,nProcessors)=MaxIndex
!            WRITE(6,*) "SubListInds(1,:) ", SubListInds(1,:)
!            WRITE(6,*) "SubListInds(2,:) ", SubListInds(2,:)
!            WRITE(6,*) "Original hash list is: "
!Reorder the lists so that they are in numerical order.
            j=1
            do while(j.le.MaxIndex)
                do i=1,nProcessors
                    IF(SubListInds(1,i).le.SubListInds(2,i)) THEN
!This block still has hashes which want to be sorted
                        MinHash=HashArray2(SubListInds(1,i))
                        MinProc=i
                        MinInd=SubListInds(1,i)
                        EXIT
                    ENDIF
!                    IF(i.eq.nProcessors) THEN
!                        WRITE(6,*) "ERROR HERE!!"
!                        CALL FLUSH(6)
!                    ENDIF
                enddo
                IF(MinHash.ne.HashCurr) THEN
                    do i=MinProc+1,nProcessors
                        IF((SubListInds(1,i).le.SubListInds(2,i)).and.(HashArray2(SubListInds(1,i)).lt.MinHash)) THEN
                            MinHash=HashArray2(SubListInds(1,i))
                            MinProc=i
                            MinInd=SubListInds(1,i)
                            IF(MinHash.eq.HashCurr) THEN
                                EXIT
                            ENDIF
                        ENDIF
                    enddo
                ENDIF
!Next smallest hash is MinHash - move the ordered elements into the other array.
                HashArray1(j)=MinHash
                IndexTable1(j)=IndexTable2(MinInd)
                ProcessVec1(j)=ProcessVec2(MinInd)
                SpawnedSign(j)=SpawnedSign2(MinInd)
                HashCurr=MinHash
!Move through the block
                j=j+1
                SubListInds(1,MinProc)=SubListInds(1,MinProc)+1
            enddo

            IF((j-1).ne.MaxIndex) THEN
                CALL Stop_All(this_routine,"Error here in the merge sort algorithm")
            ENDIF

!Need to copy the lists back to the original array to fit in with the rest of the code
            do i=1,MaxIndex
                IndexTable2(i)=IndexTable1(i)
                ProcessVec2(i)=ProcessVec1(i)
                SpawnedSign2(i)=SpawnedSign(i)
                HashArray2(i)=HashArray1(i)
            enddo

        ENDIF

!Work out the index of the particles which want to be annihilated
        j=1
        ToAnnihilateIndex=1
        do while(j.le.MaxIndex)
            TotWalkersDet=0
            InitialBlockIndex=j
            FinalBlockIndex=j-1         !Start at j-1 since we are increasing FinalBlockIndex even with the first det in the next loop
            HashCurr=HashArray2(j)
            do while((HashArray2(j).eq.HashCurr).and.(j.le.MaxIndex))
!First loop counts walkers in the block - TotWalkersDet is then the residual sign of walkers on that determinant
                TotWalkersDet=TotWalkersDet+SpawnedSign2(j)

!                IF(SpawnedSign2(j).eq.1) THEN
!                    TotWalkersDet=TotWalkersDet+1
!                ELSE
!                    TotWalkersDet=TotWalkersDet-1
!                ENDIF
                FinalBlockIndex=FinalBlockIndex+1
                j=j+1
            enddo

!            IF((Iter.eq.1877)) THEN
!                WRITE(6,*) "Common block of dets found from ",InitialBlockIndex," ==> ",FinalBlockIndex
!                WRITE(6,*) "Sum of signs in block is: ",TotWalkersDet,HashCurr
!                do k=InitialBlockIndex,FinalBlockIndex
!                    WRITE(6,*) TotWalkersDet,ToAnnihilateIndex,IndexTable2(k),ProcessVec2(k),SpawnedSign2(k)
!                enddo
!                CALL FLUSH(6)
!            ENDIF
!We need to now run through the block, and count of the same number of surviving particles as given by TotWalkersDet
    ! 1. If particles are of opposite sign, then annihilation
    ! 2. If particles are of same sign, then count out until we have the required number and annihilate the rest.
    ! Now, the sign has to be passed back. This will indicate the sign of the SURVIVING particles on that determinant.
    ! ToAnnihilateIndex now indicates the number of particles who want their sign changed at all...

            do k=InitialBlockIndex,FinalBlockIndex
!Second run through the block of same determinants marks walkers for annihilation
                IF(TotWalkersDet.eq.0) THEN
!All walkers in block want to be annihilated from now on.
                    IndexTable1(ToAnnihilateIndex)=IndexTable2(k)
                    ProcessVec1(ToAnnihilateIndex)=ProcessVec2(k)
                    SpawnedSign(ToAnnihilateIndex)=0   
                    ToAnnihilateIndex=ToAnnihilateIndex+1
                    Annihilated=Annihilated+abs(SpawnedSign2(k))
!                    TotSpawned=TotSpawned-abs(SpawnedSign2(k))
                ELSEIF((TotWalkersDet.lt.0).and.(SpawnedSign2(k).gt.0)) THEN
!Annihilate if block has a net negative walker count, and current walker is positive
                    IndexTable1(ToAnnihilateIndex)=IndexTable2(k)
                    ProcessVec1(ToAnnihilateIndex)=ProcessVec2(k)
                    SpawnedSign(ToAnnihilateIndex)=0
                    ToAnnihilateIndex=ToAnnihilateIndex+1
                    Annihilated=Annihilated+SpawnedSign2(k)
!                    TotSpawned=TotSpawned-SpawnedSign2(k)
                ELSEIF((TotWalkersDet.gt.0).and.(SpawnedSign2(k).lt.0)) THEN
!Annihilate if block has a net positive walker count, and current walker is negative
                    IndexTable1(ToAnnihilateIndex)=IndexTable2(k)
                    ProcessVec1(ToAnnihilateIndex)=ProcessVec2(k)
                    SpawnedSign(ToAnnihilateIndex)=0
                    ToAnnihilateIndex=ToAnnihilateIndex+1
                    Annihilated=Annihilated-SpawnedSign2(k)
!                    TotSpawned=TotSpawned+SpawnedSign2(k)
                ELSE
!If net walkers is positive, and we have a positive walkers, then remove one from the net positive walkers and continue through the block
!Now, we have a particle which is the same sign as the residual sign we want to pass through.
!If the sign on the particle is equal to, or less than the residual sign, then we want to let all particles live.
!Otherwise, we want to annihilate a fraction of them...
                    IF((abs(TotWalkersDet)).ge.(abs(SpawnedSign2(k)))) THEN
!All these particles are ok to be transferred accross...Increase (SpawnedSign2(k) < 0) the sign on totwalkersdet)
                        TotWalkersDet=TotWalkersDet-SpawnedSign2(k)
                    ELSE
!There is a greater number of particles in this entry than the total residual sign. Therefore, this entry want to be PARTIALLY annihilated.
!SpawnedSign will indicate the number of particles we want to remain on this entry.
                        IndexTable1(ToAnnihilateIndex)=IndexTable2(k)
                        ProcessVec1(ToAnnihilateIndex)=ProcessVec2(k)
                        Annihilated=Annihilated+(abs(SpawnedSign2(k))-abs(TotWalkersDet))
!                        TotSpawned=TotSpawned-(abs(SpawnedSign2(k))-abs(TotWalkersDet))
                        SpawnedSign(ToAnnihilateIndex)=TotWalkersDet    !The number of particles that we want left to copy accross is simply the remaining residual sign
                        ToAnnihilateIndex=ToAnnihilateIndex+1
                        TotWalkersDet=0     !All the residual sign has now been compensated for.

                    ENDIF

                ENDIF
            enddo
            IF(TotWalkersDet.ne.0) THEN
                CALL Stop_All("AnnihilateBetweenSpawned","Problem counting residual sign...")
            ENDIF

        enddo

        ToAnnihilateIndex=ToAnnihilateIndex-1   !ToAnnihilateIndex now tells us the total number of particles to annihilate from the list on this processor

!The annihilation is complete - particles to be annihilated are stored in IndexTable and need to be sent back to their original processor
!To know which processor that is, we need to order the particles to be annihilated in terms of their CPU, i.e. ProcessVec(1:ToAnnihilateIndex)
!Is the list already ordered according to CPU? Is this further sort even necessary?

        IF(ToAnnihilateIndex.gt.1) THEN
!Do not actually have to take indextable, hash2array or newsign with it...
            CALL Sort2IILongI(ToAnnihilateIndex,ProcessVec1(1:ToAnnihilateIndex),IndexTable1(1:ToAnnihilateIndex),HashArray1(1:ToAnnihilateIndex),SpawnedSign(1:ToAnnihilateIndex))
        ENDIF

!We now need to regenerate sendcounts and disps
        sendcounts(1:nProcessors)=0
        do i=1,ToAnnihilateIndex
            IF(ProcessVec1(i).gt.(nProcessors-1)) THEN
                WRITE(6,*) i,ToAnnihilateIndex
                WRITE(6,*) "***"
                WRITE(6,*) ProcessVec1(1:ToAnnihilateIndex)
                WRITE(6,*) "***"
                WRITE(6,*) sendcounts(:)
                WRITE(6,*) "***"
                WRITE(6,*) HashArray(1:ToAnnihilateIndex)
                WRITE(6,*) "***"
                WRITE(6,*) IndexTable1(1:ToAnnihilateIndex)

                CALL Stop_All("AnnihilateBetweenSpawned","Annihilation error")
            ENDIF
            sendcounts(ProcessVec1(i)+1)=sendcounts(ProcessVec1(i)+1)+1
        enddo
!The disps however do want to be cumulative
        disps(1)=0      !Starting element is always the first element
        do i=2,nProcessors
            disps(i)=disps(i-1)+sendcounts(i-1)
        enddo

!We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
        recvcounts(1:nProcessors)=0

        CALL MPI_AlltoAll(sendcounts,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,error)

!We can now get recvdisps from recvcounts in the same way we obtained disps from sendcounts
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo

        ToAnnihilateonProc=recvdisps(nProcessors)+recvcounts(nProcessors)

        CALL MPI_AlltoAllv(IndexTable1(1:ToAnnihilateonProc),sendcounts,disps,MPI_INTEGER,IndexTable2,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        CALL MPI_AlltoAllv(SpawnedSign(1:ToAnnihilateonProc),sendcounts,disps,MPI_INTEGER,SpawnedSign2(1:),recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)

!We now need to take with the index, the sign to remain on the entry, as it does not necessarily want to be totally annihilated.
        CALL NECI_SORT2I(ToAnnihilateonProc,IndexTable2(1:ToAnnihilateonProc),SpawnedSign2(1:ToAnnihilateonProc))
!        CALL SORTIILongL(ToAnnihilateonProc,Index2Table(1:ToAnnihilateonProc),HashArray(1:ToAnnihilateonProc),CurrentSign(1:ToAnnihilateonProc))

        IF(ToAnnihilateonProc.ne.0) THEN
!Copy across the data, apart from ones which have an index given by the indicies in Index2Table(1:ToAnnihilateonProc)
            VecSlot=1       !VecSlot is the index in the final array of TotWalkers
            i=1             !i is the index in the original array of TotWalkersNew
            do j=1,ToAnnihilateonProc
!Loop over all particles to be annihilated
                do while(i.lt.IndexTable2(j))
!Copy accross all particles less than this number
                    SpawnedParts2(:,VecSlot)=SpawnedParts(:,i)
                    SpawnedSign(VecSlot)=TempSign(i)
!                    IF(SpawnedSign(VecSlot).eq.0) THEN
!                        CALL Stop_All("AnnihilateBetweenSpawned","Should have non-zero number of particles in this entry")
!                    ENDIF
                    i=i+1
                    VecSlot=VecSlot+1
                enddo
                IF(SpawnedSign2(j).ne.0) THEN
!We want the entry to be partially annihilated. Keep the particle, but change its value to be that given by SpawnedSign2
                    SpawnedParts2(:,VecSlot)=SpawnedParts(:,i)
                    IF(abs(SpawnedSign2(j)).ge.abs(TempSign(i))) THEN
                        WRITE(6,*) "***",Iter,ToannihilateonProc,ValidSpawned
                        CALL Stop_All("AnnihilateBetweenSpawned","Incorrect annihilating here...")
                    ENDIF
                    SpawnedSign(VecSlot)=SpawnedSign2(j)
                    IF(SpawnedSign(VecSlot).eq.0) THEN
                        CALL Stop_All("AnnihilateBetweenSpawned","Should have non-zero number of particles in this entry")
                    ENDIF
                    VecSlot=VecSlot+1
                ENDIF
                i=i+1
            enddo

!Now need to copy accross the residual - from Index2Table(ToAnnihilateonProc) to TotWalkersNew
            do i=IndexTable2(ToAnnihilateonProc)+1,ValidSpawned
                SpawnedParts2(:,VecSlot)=SpawnedParts(:,i)
                SpawnedSign(VecSlot)=TempSign(i)
                VecSlot=VecSlot+1
            enddo

        ELSE
!No particles annihilated
            VecSlot=1
            do i=1,ValidSpawned
                SpawnedParts2(:,VecSlot)=SpawnedParts(:,i)
                SpawnedSign(VecSlot)=TempSign(i)
                VecSlot=VecSlot+1
            enddo
        ENDIF

!Have to swap arrays around here, since the pointers must stay in sync with the arrays they're pointing at.
        ValidSpawned=VecSlot-1
        do i=1,ValidSpawned
            SpawnedSign2(i)=SpawnedSign(i)
        enddo

!        IF((TotWalkersNew-TotWalkers).ne.ToAnnihilateonProc) THEN
!            WRITE(6,*) TotWalkers,TotWalkersNew,ToAnnihilateonProc,Iter
!            CALL FLUSH(6)
!            CALL Stop_All("AnnihilatePartPar","Problem with numbers when annihilating")
!        ENDIF

!Deallocate temp arrays
        DEALLOCATE(TempSign)
        DEALLOCATE(HashArray1)
        DEALLOCATE(HashArray2)
        DEALLOCATE(IndexTable1)
        DEALLOCATE(IndexTable2)
        DEALLOCATE(ProcessVec1)
        DEALLOCATE(ProcessVec2)

!We also need to swap round the pointers to the two arrays, since the next annihilation steps take place on SpawnedParts, not SpawnedParts2 
        IF(associated(SpawnedParts2,target=SpawnVec2)) THEN
            SpawnedParts2 => SpawnVec
            SpawnedSign2 => SpawnSignVec
            SpawnedParts => SpawnVec2
            SpawnedSign => SpawnSignVec2
        ELSE
            SpawnedParts => SpawnVec
            SpawnedSign => SpawnSignVec
            SpawnedParts2 => SpawnVec2
            SpawnedSign2 => SpawnSignVec2
        ENDIF

        CALL halt_timer(AnnSpawned_time)

    END SUBROUTINE AnnihilateBetweenSpawned

    SUBROUTINE LinSearchParts(DetArray,iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: iLut(0:NIfD),MinInd,MaxInd,PartInd,DetArray(0:NIfD,1:MaxInd)
        INTEGER :: i,j,N,Comp,DetBitLT
        LOGICAL :: tSuccess

        N=MinInd
        do while(N.le.MaxInd)
            Comp=DetBitLT(DetArray(:,N),iLut(:),NIfD)
            IF(Comp.eq.1) THEN
                N=N+1
            ELSEIF(Comp.eq.-1) THEN
                PartInd=N-1
                tSuccess=.false.
                RETURN
            ELSE
                tSuccess=.true.
                PartInd=N
                RETURN
            ENDIF
        enddo
        tSuccess=.false.
        PartInd=MaxInd-1

    END SUBROUTINE LinSearchParts

!Do a binary search in CurrentDets, between the indices of MinInd and MaxInd. If successful, tSuccess will be true and 
!PartInd will be a coincident determinant. If there are multiple values, the chosen one may be any of them...
!If failure, then the index will be one less than the index that the particle would be in if it was present in the list.
!(or close enough!)
    SUBROUTINE BinSearchParts(iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: iLut(0:NIfD),MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp,DetBitLT
        LOGICAL :: tSuccess

        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(CurrentDets(:,MaxInd),iLut(:),NIfD)
            IF(Comp.eq.0) THEN
                tSuccess=.true.
                PartInd=MaxInd
                RETURN
            ELSE
                tSuccess=.false.
                PartInd=MinInd
            ENDIF
        ENDIF

        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(CurrentDets(:,N),iLut(:),NIfD)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
                    Comp=DetBitLT(CurrentDets(:,i+1),iLut(:),NIfD)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN

                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchParts

!Do a binary search the guiding function dets, between the indices of MinInd and MaxInd. If successful, tSuccess will be true and 
!PartInd will be a coincident determinant. If there are multiple values, the chosen one may be any of them...
!If failure, then the index will be one less than the index that the particle would be in if it was present in the list.
!(or close enough!)
    SUBROUTINE BinSearchGuideParts(iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: iLut(0:NIfD),MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp,DetBitLT
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
        i=MinInd
        j=MaxInd
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(GuideFuncDets(:,N),iLut(:),NIfD)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
                    Comp=DetBitLT(GuideFuncDets(:,i+1),iLut(:),NIfD)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN

                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchGuideParts


!Do a binary search in MinorStarDets, between the indices of MinInd and MaxInd. If successful, tSuccess will be true and 
!PartInd will be a coincident determinant. If there are multiple values, the chosen one may be any of them...
!If failure, then the index will be one less than the index that the particle would be in if it was present in the list.
!(or close enough!)
    SUBROUTINE BinSearchMinorParts(iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: iLut(0:NIfD),MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp,DetBitLT
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
        i=MinInd
        j=MaxInd
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(MinorStarDets(:,N),iLut(:),NIfD)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
                    Comp=DetBitLT(MinorStarDets(:,i+1),iLut(:),NIfD)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN

                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchMinorParts


!Do a binary search of the DominantDets, between the indices of MinInd and MaxInd. If successful, tSuccess will be true.
    SUBROUTINE BinSearchDomParts(AllExcDets,iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: iLut(0:NIfD),MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp,DetBitLT,AllExcDets(0:NIfD,MinInd:MaxInd)
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
        i=MinInd
        j=MaxInd
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(AllExcDets(0:NIfD,N),iLut(:),NIfD)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
                    Comp=DetBitLT(AllExcDets(0:NIfD,i+1),iLut(:),NIfD)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN

                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchDomParts




!In this routine, we want to search through the list of spawned particles. For each spawned particle, we binary search the list of particles on the processor
!to see if an annihilation event can occur. The annihilated particles are then removed from the spawned list
!to the whole list of spawned particles at the end of the routine.
!In the main list, we change the 'sign' element of the array to zero. These will be deleted at the end of the total annihilation step.
    SUBROUTINE AnnihilateSpawnedParts(ValidSpawned,TotWalkersNew)
        INTEGER :: ValidSpawned,MinInd,TotWalkersNew,PartInd,i,j,k,ToRemove,VecInd,SignProd,DetsMerged,PartIndex!,SearchInd,AnnihilateInd
        INTEGER :: ExcitLevel
        LOGICAL :: DetBitEQ,tSuccess,tSuc!,tSkipSearch

        CALL set_timer(AnnMain_time,30)
!        IF(Iter.eq.1877) THEN
!            WRITE(6,*) "MainList: ",TotWalkersNew
!            do i=1,TotWalkersNew
!                WRITE(6,*) CurrentDets(:,i),CurrentSign(i)
!            enddo
!            WRITE(6,*) "** ** **"
!            do i=1,ValidSpawned
!                WRITE(6,*) SpawnedParts(:,i),SpawnedSign(i)
!            enddo
!            WRITE(6,*) "****************"
!        ENDIF

!MinInd indicates the minimum bound of the main array in which the particle can be found.
!Since the spawnedparts arrays are ordered in the same fashion as the main array, we can find the particle position in the main array by only searching a subset.
        MinInd=1
        IF(tHistSpawn) HistMinInd2(1:NEl)=FCIDetIndex(1:NEl)
        ToRemove=0  !The number of particles to annihilate
!        WRITE(6,*) "Annihilating between ",ValidSpawned, " spawned particles and ",TotWalkersNew," original particles..."
!        WRITE(6,*) "SpawnedParts: "
!        do i=1,ValidSpawned
!            WRITE(6,*) SpawnedParts(:,i),SpawnedSign(i)
!        enddo
!        WRITE(6,*) "Original Parts: "
!        do i=1,TotWalkersNew
!            WRITE(6,*) CurrentDets(:,i),CurrentSign(i)
!        enddo
!        CALL FLUSH(6)
        
        CALL set_timer(BinSearch_time,45)

        do i=1,ValidSpawned

!This will binary search the CurrentDets array to find the desired particle. tSuccess will determine whether the particle has been found or not.
!It will also return the index of the position one below where the particle would be found if was in the list.
!            CALL LinSearchParts(CurrentDets(:,1:TotWalkersNew),SpawnedParts(0:NIfD,i),MinInd,TotWalkersNew,PartInd,tSuccess)
            CALL BinSearchParts(SpawnedParts(:,i),MinInd,TotWalkersNew,PartInd,tSuccess)
!            WRITE(6,*) "Binary search complete: ",i,PartInd,tSuccess
!            CALL FLUSH(6)

            IF(tSuccess) THEN
!A particle on the same list has been found. We now want to search backwards, to find the first particle in this block.
!Actually, this shouldn't be necessary - the CurrentDets array should be sign-coherent. The only time that we need to search forwards/backwards is if we hit upon an already
!annihilated particle, i.e. Sign=0. If we hit upon +-1, then we know that the block is sign coherent.

!                SearchInd=PartInd   !This can actually be min(1,PartInd-1) once we know that the binary search is working, as we know that PartInd is the same particle.
!                MinInd=PartInd      !Make sure we only have a smaller list to search next time since the next particle will not be at an index smaller than PartInd
!                AnnihilateInd=0     !AnnihilateInd indicates the index in CurrentDets of the particle we want to annihilate. It will remain 0 if we find not complimentary particle.
!                tSkipSearch=.false. !This indicates whether we want to continue searching forwards through the list once we exit the loop going backwards.
                
                SignProd=CurrentSign(PartInd)*SpawnedSign(i)
                IF(SignProd.lt.0) THEN
!This indicates that the particle has found the same particle of opposite sign to annihilate with
!Mark these particles for annihilation in both arrays
!If we go to a determinant representation of the spawned particles, then we need to be careful that we can only annihilate against the number of particles on the main list.
!We cannot transfer the rest of the particles across, since we rely on the fact that the main arrays are sign-coherent with each other.
!This means that at the end, only one sign of a determinant will exist, whether on the main array, or spawned array.
!                    AnnihilateInd=SearchInd
                    IF(abs(SpawnedSign(i)).ge.abs(CurrentSign(PartInd))) THEN
!There are more (or equal) numbers of spawned particles to annihilate. We can only annihilate some from the spawned list, but all from main list (or all from both if equal and opposite).
                        SpawnedSign(i)=SpawnedSign(i)+CurrentSign(PartInd)
                        Annihilated=Annihilated+2*(abs(CurrentSign(PartInd)))
                        
                        IF(tHistSpawn) THEN
!We want to histogram where the particle annihilations are taking place.
                            CALL FindBitExcitLevel(SpawnedParts(:,i),iLutHF,NIfD,ExcitLevel,NEl)
                            IF(ExcitLevel.eq.NEl) THEN
                                CALL BinSearchParts2(SpawnedParts(:,i),HistMinInd2(ExcitLevel),Det,PartIndex,tSuc)
                            ELSEIF(ExcitLevel.eq.0) THEN
                                PartIndex=1
                                tSuc=.true.
                            ELSE
                                CALL BinSearchParts2(SpawnedParts(:,i),HistMinInd2(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartIndex,tSuc)
                            ENDIF
                            HistMinInd2(ExcitLevel)=PartIndex
                            IF(tSuc) THEN
                                AvAnnihil(PartIndex)=AvAnnihil(PartIndex)+REAL(2*(abs(CurrentSign(PartInd))),r2)
                                InstAnnihil(PartIndex)=InstAnnihil(PartIndex)+REAL(2*(abs(CurrentSign(PartInd))),r2)
                            ELSE
                                WRITE(6,*) "***",SpawnedParts(0:NIfD,i)
                                Call WriteBitDet(6,SpawnedParts(0:NIfD,i),.true.)
                                CALL Stop_All("AnnihilateSpawnedParts","Cannot find corresponding FCI determinant when histogramming")
                            ENDIF
                        ENDIF

                        CurrentSign(PartInd)=0
                        IF(SpawnedSign(i).eq.0) THEN
!The number of particles were equal and opposite. We want to remove this entry from the spawned list.
                            ToRemove=ToRemove+1
                        ENDIF

                    ELSE
!There are more particles in the main list, than the spawned list. We want to annihilate all particles from the spawned list, but only some from main list.
                        CurrentSign(PartInd)=CurrentSign(PartInd)+SpawnedSign(i)
                        Annihilated=Annihilated+2*(abs(SpawnedSign(i)))
                        
                        IF(tHistSpawn) THEN
!We want to histogram where the particle annihilations are taking place.
                            CALL FindBitExcitLevel(SpawnedParts(:,i),iLutHF,NIfD,ExcitLevel,NEl)
                            IF(ExcitLevel.eq.NEl) THEN
                                CALL BinSearchParts2(SpawnedParts(:,i),HistMinInd2(ExcitLevel),Det,PartIndex,tSuc)
                            ELSEIF(ExcitLevel.eq.0) THEN
                                PartIndex=1
                                tSuc=.true.
                            ELSE
                                CALL BinSearchParts2(SpawnedParts(:,i),HistMinInd2(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartIndex,tSuc)
                            ENDIF
                            HistMinInd2(ExcitLevel)=PartIndex
                            IF(tSuc) THEN
                                AvAnnihil(PartIndex)=AvAnnihil(PartIndex)+REAL(2*(abs(SpawnedSign(i))),r2)
                                InstAnnihil(PartIndex)=InstAnnihil(PartIndex)+REAL(2*(abs(SpawnedSign(i))),r2)
                            ELSE
                                WRITE(6,*) "***",SpawnedParts(0:NIfD,i)
                                CALL Stop_All("AnnihilateSpawnedParts","Cannot find corresponding FCI determinant when histogramming")
                            ENDIF
                        ENDIF

                        SpawnedSign(i)=0
                        ToRemove=ToRemove+1
                    ENDIF

                        
!                    IF(CurrentSign(PartInd).gt.0) THEN
!                        CurrentSign(PartInd)=CurrentSign(PartInd)-1
!                    ELSE
!                        CurrentSign(PartInd)=CurrentSign(PartInd)+1
!                    ENDIF
!                    SpawnedSign(i)=0
!                    ToRemove=ToRemove+1
!!                    RemoveInds(ToRemove)=i  !This is the index of the spawned particle to remove.
!                    Annihilated=Annihilated+2   !Count that we have annihilated two particles

                ELSEIF(SignProd.gt.0) THEN
!This indicates that the particle has found a similar particle of the same sign. It therefore cannot annihilate, since all arrays accross all processors are sign-coherent.
!Therefore, we can just transfer it accross now.
                    CurrentSign(PartInd)=CurrentSign(PartInd)+SpawnedSign(i)
!We have transferred a particle accross between processors. "Annihilate" from the spawned list, but not the main list.
                    SpawnedSign(i)=0
                    ToRemove=ToRemove+1
!                    RemoveInds(ToRemove)=i
!                    AnnihilateInd=-SearchInd
                ELSE
!One of the signs on the list is actually 0. If this zero is on the spawned list, we need to mark it for removal.
                    IF(SpawnedSign(i).eq.0) THEN
                        ToRemove=ToRemove+1
                    ENDIF
                ENDIF


!                do while((DetBitEQ(SpawnedParts(:,i),CurrentDets(:,SearchInd),NIfD)).and.(SearchInd.ge.1))
!!Cycle backwards through the list, checking where the start of this block of determinants starts.
!                    SignProd=CurrentSign(SearchInd)*SpawnedSign(i)
!                    IF(SignProd.lt.0) THEN
!!We have actually found a complimentary particle - mark the index of this particle for annihilation.
!                        AnnihilateInd=SearchInd
!!                        WRITE(6,"(A,2I12,I4,2I12,I4)") "Annihilated from MainList: ",SpawnedParts(:,i),SpawnedSign(i),CurrentDets(:,SearchInd),CurrentSign(SearchInd)
!                        tSkipSearch=.true.
!                        EXIT
!                    ELSEIF(SignProd.gt.0) THEN
!!Since the signs are coherent on CurrentSign, we know that we can not annihilate the SpawnedParts particle if we find a particle of the same sign. 
!!Therefore, we can remove the particle from the list, and instantly transfer the particle to the next array.
!                        tSkipSearch=.true.
!                        CurrentSign(SearchInd)=CurrentSign(SearchInd)+SpawnedSign(i)
!                        AnnihilateInd=-SearchInd    !Give the annihilateInd a negative index to indicate that we only want to annihilate from the spawned list, not main list.
!                        EXIT
!                    ENDIF
!
!                    SearchInd=SearchInd-1
!                    
!                enddo
!
!!                IF((SearchInd.eq.PartInd).and.(AnnihilateInd.eq.0)) THEN
!!!The searchind should not equal partind, since we know that the particles are the same at PartInd, otherwise the binary search should have returned false.
!!!(unless we have already found the particle to annihilate)
!!                    CALL Stop_All("AnnihilateSpawnedParts","Binary search has fatal error")
!!                ENDIF
!
!                IF(.not.tSkipSearch) THEN
!!We have searched from the beginning of the particle block(SearchInd) to PartInd for a complimentary particle, but have not had any success. Now we can search from
!!PartInd+1 to the end of the block for a complimentary particle.
!                    SearchInd=PartInd+1
!                    do while((DetBitEQ(SpawnedParts(0:NIfD,i),CurrentDets(0:NIfD,SearchInd),NIfD)).and.(SearchInd.le.TotWalkersNew))
!                        SignProd=CurrentSign(SearchInd)*SpawnedSign(i)
!                        IF(SignProd.lt.0) THEN
!!We have found a complimentary particle - mark the index of this particle for annihilation.
!!                            WRITE(6,"(A,2I12,I4,2I12,I4)") "Annihilated from MainList: ",SpawnedParts(:,i),SpawnedSign(i),CurrentDets(:,SearchInd),CurrentSign(SearchInd)
!                            AnnihilateInd=SearchInd
!                            EXIT
!                        ELSEIF(SignProd.gt.0) THEN
!                            CurrentSign(SearchInd)=CurrentSign(SearchInd)+SpawnedSign(i)
!                            AnnihilateInd=-SearchInd    !Give the annihilateInd a negative index to indicate that we only want to annihilate from the spawned list, not main list.
!                            EXIT
!                        ENDIF
!
!                        SearchInd=SearchInd+1
!                    enddo
!
!!We can now also move the MinInd to the end of the block if we want
!                    MinInd=SearchInd-1     !This cannot be more than TotWalkersNew
!                ENDIF
!
!                IF(AnnihilateInd.gt.0) THEN
!!We have found a particle to annihilate. Mark the particles for annihilation.
!                    IF(CurrentSign(AnnihilateInd).gt.0) THEN
!                        CurrentSign(AnnihilateInd)=CurrentSign(AnnihilateInd)-1
!                    ELSE
!                        CurrentSign(AnnihilateInd)=CurrentSign(AnnihilateInd)+1
!                    ENDIF
!                    SpawnedSign(i)=0
!                    ToRemove=ToRemove+1
!                    RemoveInds(ToRemove)=i  !This is the index of the spawned particle to remove.
!                    Annihilated=Annihilated+2   !Count that we have annihilated two particles
!
!                ELSEIF(AnnihilateInd.lt.0) THEN
!!We have transferred a particle accross between processors. "Annihilate" from the spawned list, but not the main list.
!                    SpawnedSign(i)=0
!                    ToRemove=ToRemove+1
!                    RemoveInds(ToRemove)=i
!                ENDIF

            ENDIF

!Even if a corresponding particle wasn't found, we can still search a smaller list next time....so not all bad news then...
            MinInd=PartInd

        enddo
        
        CALL halt_timer(BinSearch_time)

!Now we have to remove the annihilated particles from the spawned list. They will be removed from the main list at the end of the annihilation process.
!It may actually be easier to just move the annihilated particles to the end of the list and resort the list?
!Or, the removed indices could be found on the fly? This may have little benefit though if the memory isn't needed.
        IF(ToRemove.gt.0) THEN

!            VecInd=1
!            do i=1,ValidSpawned
!                IF(SpawnedSign(i).ne.0) THEN
!                    SpawnedParts2(:,VecInd)=SpawnedParts(:,i)
!                    SpawnedSign2(VecInd)=SpawnedSign(i)
!                    VecInd=VecInd+1
!                ENDIF
!            enddo
!            ValidSpawned=ValidSpawned-ToRemove
!            IF((VecInd-1).ne.ValidSpawned) THEN
!                CALL Stop_All("AnnihilateSpawnedParts","Not all spawned particles correctly annihilated.")
!            ENDIF
!            do i=1,ValidSpawned
!                SpawnedParts(:,i)=SpawnedParts2(:,i)
!                SpawnedSign(i)=SpawnedSign2(i)
!            enddo

!Since reading and writing from the same array is slow, copy the information accross to the other spawned array, and just swap the pointers around after.
            DetsMerged=0
            do i=1,ValidSpawned
!We want to move all the elements above this point down to 'fill in' the annihilated determinant.
                IF(SpawnedSign(i).eq.0) THEN
                    DetsMerged=DetsMerged+1
                ELSE
                    SpawnedParts2(0:NIfD,i-DetsMerged)=SpawnedParts(0:NIfD,i)
                    SpawnedSign2(i-DetsMerged)=SpawnedSign(i)
                ENDIF
            enddo
            ValidSpawned=ValidSpawned-DetsMerged
            IF(DetsMerged.ne.ToRemove) THEN
                WRITE(6,*) "***", Iter
                CALL Stop_All("AnnihilateSpawnedParts","Incorrect number of particles removed from spawned list")
            ENDIF
!We always want to annihilate from the SpawedParts and SpawnedSign arrays.
            IF(associated(SpawnedParts2,target=SpawnVec2)) THEN
                SpawnedParts2 => SpawnVec
                SpawnedSign2 => SpawnSignVec
                SpawnedParts => SpawnVec2
                SpawnedSign => SpawnSignVec2
            ELSE
                SpawnedParts => SpawnVec
                SpawnedSign => SpawnSignVec
                SpawnedParts2 => SpawnVec2
                SpawnedSign2 => SpawnSignVec2
            ENDIF



!            do i=1,ToRemove
!
!                do j=RemoveInds(i)+1-(i-1),ValidSpawned
!                    SpawnedParts(:,j-1)=SpawnedParts(:,j)
!                    SpawnedSign(j-1)=SpawnedSign(j)
!                enddo
!                ValidSpawned=ValidSpawned-1
!            enddo

        ENDIF

!        do i=1,ValidSpawned
!            WRITE(6,*) SpawnedParts(:,i),SpawnedSign(i)
!            IF(SpawnedSign(i).eq.0) THEN
!                CALL Stop_All("AnnihilateSpawnedParts","Not all spawned particles correctly annihilated")
!            ENDIF
!        enddo
!        CALL CheckOrdering(SpawnedParts,SpawnedSign(1:ValidSpawned),ValidSpawned,.true.)


!            i=1
!            do while(SpawnedSign(i).ne.0)
!                i=i+1
!            enddo
!!i now indicates the index of the first particle to remove.
!            MinInd=i+1    !MinInd now indicates the index of the beginning of the block for which to move 
!            i=MinInd
!            do j=1,ToRemove
!!Run over the number of particles to annihilate. Each particle will mean that a block of particles will be shifted down.
!                do while((SpawnedSign(i).ne.0).and.(i.le.ValidSpawned))
!                    i=i+1
!                enddo
!!We now want to move the block, donoted by MinInd -> i-1 down by j steps. Can we do this as a array operation?
!!                SpawnedParts(:,MinInd-j:i-1-j)=SpawnedParts(:,MinInd:i-1)
!                do k=MinInd,i-1
!                    SpawnedParts(:,k-j)=SpawnedParts(:,k)
!                    SpawnedSign(k-j)=SpawnedSign(k)
!                enddo
!                MinInd=i+1
!                i=MinInd
!            enddo
!
!!Reduce ValidSpawned by the number of newly-spawned particles which have been annihilated.
!            ValidSpawned=ValidSpawned-ToRemove
!
!        ENDIF

        CALL halt_timer(AnnMain_time)

    END SUBROUTINE AnnihilateSpawnedParts
        
!This rotates the spawned (and still alive particles) around the processors. Particles are sent to MOD(iProcIndex+1,nProcessors) and received from MOD(iProcIndex+nProcessors-1,nProcessors).
!Issues here:
!1) Want to avoid deadlock, but also want to avoid having to send data sequentially, therefore blocking is going to be necessary.
!2) This will also mean we have to beware of buffer overflow. Do we need to attach a specific buffer for the particles?
!3) Do we want one of two sets of data? If two, then we need to set up a pointer system. If not, then how do we know how many particles to recieve without
!       extra communication?
    SUBROUTINE RotateParticles(ValidSpawned)
        INTEGER :: error,ValidSpawned
        INTEGER, DIMENSION(MPI_STATUS_SIZE) :: Stat 
        
        CALL set_timer(Comms_Time,30)

!ValidSpawned is the number of particles spawned (and still alive) for this set of particles (index is iProcIndex-no.rotates)
        SpawnedSign(0)=ValidSpawned

!        WRITE(6,*) "Particles to send: ",ValidSpawned
!        CALL FLUSH(6)

!Send the signs of the particles (number sent is in the first element)
        CALL MPI_BSend(SpawnedSign(0:ValidSpawned),ValidSpawned+1,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),123,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in sending signs")
        ENDIF
!        WRITE(6,*) "Sent sign",ValidSpawned+1
!        CALL FLUSH(6)

!...and then send the particles themselves...
        CALL MPI_BSend(SpawnedParts(0:NIfD,1:ValidSpawned),ValidSpawned*(NIfD+1),MPI_INTEGER,MOD(iProcIndex+1,nProcessors),456,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in sending particles")
        ENDIF
!        WRITE(6,*) "Sent particles",ValidSpawned
!        CALL FLUSH(6)

!Receive signs (let it receive the maximum possible (only the first ValidSpawned will be updated.))
        CALL MPI_Recv(SpawnedSign2(0:MaxSpawned),MaxSpawned+1,MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),123,MPI_COMM_WORLD,Stat,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in receiving signs")
        ENDIF
!        WRITE(6,*) "Recieved sign",MaxSpawned+1
!        CALL FLUSH(6)

!Update the ValidSpawned variable for this new set of data we are about to receive...
        ValidSpawned=SpawnedSign2(0)

        CALL MPI_Recv(SpawnedParts2(0:NIfD,1:ValidSpawned),ValidSpawned*(NIfD+1),MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),456,MPI_COMM_WORLD,Stat,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in receiving particles")
        ENDIF
!        WRITE(6,*) "Recieving particles, ",ValidSpawned
!        CALL FLUSH(6)

!We now want to make sure that we are working on the correct array. We have now received particles in SpawnedParts2 - switch it so that we are pointing at the other array.
!We always want to annihilate from the SpawedParts and SpawnedSign arrays.
        IF(associated(SpawnedParts2,target=SpawnVec2)) THEN
            SpawnedParts2 => SpawnVec
            SpawnedSign2 => SpawnSignVec
            SpawnedParts => SpawnVec2
            SpawnedSign => SpawnSignVec2
        ELSE
            SpawnedParts => SpawnVec
            SpawnedSign => SpawnSignVec
            SpawnedParts2 => SpawnVec2
            SpawnedSign2 => SpawnSignVec2
        ENDIF
!        WRITE(6,*) "Switched arrays around..."
!        CALL FLUSH(6)

        CALL halt_timer(Comms_Time)

    END SUBROUTINE RotateParticles

    SUBROUTINE CheckOrdering(DetArray,SignArray,NoDets,tCheckSignCoher)
        INTEGER :: NoDets,DetArray(0:NIfD,1:NoDets),DetBitLT,i,j
        INTEGER :: SignArray(1:NoDets),Comp
        LOGICAL :: tCheckSignCoher

        IF(NoDets.gt.0) THEN
            IF(SignArray(1).eq.0) THEN
                WRITE(6,*) "Iter: ",Iter,1
                CALL Stop_All("CheckOrdering","Array has annihilated particles in it...")
            ENDIF
        ENDIF
        do i=2,NoDets
            IF(SignArray(i).eq.0) THEN
                WRITE(6,*) "Iter: ",Iter,i
                CALL Stop_All("CheckOrdering","Array has annihilated particles in it...")
            ENDIF
            Comp=DetBitLT(DetArray(:,i-1),DetArray(:,i),NIfD)
            IF(Comp.eq.-1) THEN
!Array is in reverse numerical order for these particles
                do j=max(i-5,1),min(i+5,NoDets)
                    WRITE(6,*) Iter,j,DetArray(:,j),SignArray(j)
                enddo
                CALL Stop_All("CheckOrdering","Array not ordered correctly")
            ELSEIF(Comp.eq.0) THEN
!Dets are the same - see if we want to check sign-coherence
                IF(tCheckSignCoher) THEN
!!This bit checks that there is only one copy of the determinants in the list
                    do j=max(i-5,1),min(i+5,NoDets)
                        WRITE(6,*) Iter,j,DetArray(:,j),SignArray(j)
                    enddo
                    CALL Stop_All("CheckOrdering","Determinant same as previous one...")
                ENDIF
                IF(tCheckSignCoher.and.(SignArray(i-1).ne.SignArray(i))) THEN
!This checks that any multple copies in the list are sign-coherent...
                    do j=i-5,i+5
                        WRITE(6,*) Iter,j,DetArray(:,j),SignArray(j)
                    enddo
                    CALL Stop_All("CheckOrdering","Array not sign-coherent")
                ENDIF
            ENDIF
        enddo

    END SUBROUTINE CheckOrdering


!A routine to annihilate particles in parallel. This involves separating hashes by abs(mod(hash,nProc)) to each node and annihilating there,       
!before sending back the annihilated particles to be removed from their original processors.
    SUBROUTINE AnnihilatePartPar(TotWalkersNew)
        INTEGER :: i,j,k,ToAnnihilateIndex,TotWalkersNew,ierr,error,sendcounts(nProcessors)
        INTEGER :: TotWalkersDet,InitialBlockIndex,FinalBlockIndex,ToAnnihilateOnProc,VecSlot
        INTEGER :: disps(nProcessors),recvcounts(nProcessors),recvdisps(nProcessors)!,AnnihilPart(nProcessors)
        INTEGER :: Minsendcounts,Maxsendcounts,DebugIter,SubListInds(2,nProcessors),MinProc,MinInd
        REAL*8 :: PopDensity(0:NEl)
        INTEGER , ALLOCATABLE :: TempExcitLevel(:)
        INTEGER(KIND=i2) :: HashCurr,MinBin,RangeofBins,NextBinBound,MinHash
        CHARACTER(len=*), PARAMETER :: this_routine='AnnihilatePartPar'
        INTEGER , ALLOCATABLE :: TempSign(:)                                                         !Temp array to hold sign of walkers when annihilating
        INTEGER(KIND=i2) , ALLOCATABLE :: TempHash(:)
        INTEGER :: TempSignTag=0,TempHashTag=0

!This is just to see if there are higher-weighted determinants that HF...
!        AllNoatHF=0
!        CALL MPI_AllReduce(NoatHF,AllNoatHF,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)

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
        enddo
        ProcessVec(1:TotWalkersNew)=iProcIndex

!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) TotWalkersNew
!            do i=1,TotWalkersNew
!                WRITE(6,*) i,Hash2Array(i),IndexTable(i),ProcessVec(i),NewSign(i)
!            enddo
!        ENDIF

        CALL set_timer(Sort_Time,30)
!Next, order the hash array, taking the index, CPU and sign with it...
        IF(.not.tAnnihilatebyRange) THEN
!Order the array by abs(mod(Hash,nProcessors)). This will result in a more load-balanced system
!             IF(TLocalAnnihilation) THEN
!If we are locally annihilating, then we need to take the excitation level of each walker with the hash
!                 CALL SortMod4I1LLong(TotWalkersNew,Hash2Array(1:TotWalkersNew),IndexTable(1:TotWalkersNew),ProcessVec(1:TotWalkersNew),TempExcitLevel(1:TotWalkersNew),NewSign(1:TotWalkersNew),nProcessors)
!             ELSE
!            CALL SortMod3I1LLong(TotWalkersNew,Hash2Array(1:TotWalkersNew),IndexTable(1:TotWalkersNew),ProcessVec(1:TotWalkersNew),NewSign(1:TotWalkersNew),nProcessors)
            CALL SortMod4ILong(TotWalkersNew,Hash2Array(1:TotWalkersNew),IndexTable(1:TotWalkersNew),ProcessVec(1:TotWalkersNew),NewSign(1:TotWalkersNew),nProcessors)
!             ENDIF
            CALL halt_timer(Sort_Time)

!Send counts is the size of each block of ordered dets which are going to each processor. This could be binary searched for extra speed
            j=1
            do i=0,nProcessors-1    !Search through all possible values of abs(mod(Hash,nProcessors))
                do while((abs(mod(Hash2Array(j),INT(nProcessors,8))).eq.i).and.(j.le.TotWalkersNew))
                    j=j+1 
                enddo
                sendcounts(i+1)=j-1
            enddo
        
        ELSE
!We can try to sort the hashes by range, which may result in worse load-balancing, but will remove the need for a second sort of the hashes once they have been sent to the correct processor.
!            CALL Sort3I1LLong(TotWalkersNew,Hash2Array(1:TotWalkersNew),IndexTable(1:TotWalkersNew),ProcessVec(1:TotWalkersNew),NewSign(1:TotWalkersNew))
            CALL Sort4ILong(TotWalkersNew,Hash2Array(1:TotWalkersNew),IndexTable(1:TotWalkersNew),ProcessVec(1:TotWalkersNew),NewSign(1:TotWalkersNew))
            CALL halt_timer(Sort_Time)
            IF(nProcessors.ne.1) THEN
!We also need to know the ranges of the hashes to send to each processor. Each range should be the same.
                Rangeofbins=INT(HUGE(Rangeofbins)/(nProcessors/2),8)
                MinBin=-HUGE(MinBin)
                NextBinBound=MinBin+Rangeofbins
!            WRITE(6,*) "Rangeofbins: ",Rangeofbins
!            WRITE(6,*) "MinBin: ",MinBin
!            WRITE(6,*) "NextBinBound: ",NextBinBound

!We need to find the indices for each block of hashes which are to be sent to each processor.
!Sendcounts is the size of each block of ordered dets which are going to each processors. This could be binary searched for extra speed.
                j=1
                do i=1,nProcessors    !Search through all possible values of the hashes
                    do while((Hash2Array(j).le.NextBinBound).and.(j.le.TotWalkersNew))
                        j=j+1
                    enddo
                    sendcounts(i)=j-1
                    IF(i.eq.nProcessors-1) THEN
!Make sure the final bin catches everything...
                        NextBinBound=HUGE(NextBinBound)
                    ELSE
                        NextBinBound=NextBinBound+Rangeofbins
                    ENDIF
                enddo
            ELSE
                sendcounts(1)=TotWalkersNew
            ENDIF

        ENDIF
        
!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "***************"
!            WRITE(6,*) TotWalkersNew
!            do i=1,TotWalkersNew
!                WRITE(6,*) Hash2Array(i),abs(mod(Hash2Array(i),nProcessors)),IndexTable(i),ProcessVec(i),NewSign(i)
!            enddo
!        ENDIF
        
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
!        CALL MPI_Barrier(MPI_COMM_WORLD,error)
        CALL set_timer(Comms_Time,30)

        CALL MPI_AlltoAll(sendcounts,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,MPI_COMM_WORLD,error)

!We can now get recvdisps from recvcounts in the same way we obtained disps from sendcounts
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo

        MaxIndex=recvdisps(nProcessors)+recvcounts(nProcessors)
!Max index is the largest occupied index in the array of hashes to be ordered in each processor 
        IF(MaxIndex.gt.(0.95*MaxWalkersAnnihil)) THEN
            CALL Warning("AnnihilatePartPar","Maximum index of annihilation array is close to maximum length. Increase MemoryFacAnnihil")
        ENDIF
!Uncomment this if you want to write out load-balancing statistics.
!        AnnihilPart(:)=0
!        CALL MPI_Gather(MaxIndex,1,MPI_INTEGER,AnnihilPart,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
!        IF(iProcIndex.eq.root) THEN
!            WRITE(13,"(I10)",advance='no') Iter
!            do i=1,nProcessors
!                WRITE(13,"(I10)",advance='no') AnnihilPart(i)
!            enddo
!            WRITE(13,"(A)") ""
!            CALL FLUSH(13)
!        ENDIF

!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "RECVCOUNTS: "
!            WRITE(6,*) recvcounts(:)
!            WRITE(6,*) "RECVDISPS: "
!            WRITE(6,*) recvdisps(:),MaxIndex
!            CALL FLUSH(6)
!        ENDIF

!Insert a load-balance check here...maybe find the s.d. of the sendcounts array - maybe just check the range first.
!        IF(TotWalkersNew.gt.200) THEN
!            IF((Maxsendcounts-Minsendcounts).gt.(TotWalkersNew/3)) THEN
!                WRITE(6,"(A,I12)") "**WARNING** Parallel annihilation not optimally balanced on this node, for iter = ",Iter
!                WRITE(6,*) "Sendcounts is: ",sendcounts(:)
!!                CALL FLUSH(6)
!            ENDIF
!        ENDIF
        
!Now send the chunks of hashes to the corresponding processors
        CALL MPI_AlltoAllv(Hash2Array(1:TotWalkersNew),sendcounts,disps,MPI_DOUBLE_PRECISION,HashArray(1:MaxIndex),recvcounts,recvdisps,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,error)        

!The signs of the hashes, index and CPU also need to be taken with them.
        CALL MPI_AlltoAllv(NewSign(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,CurrentSign,recvcounts,recvdisps,MPI_LOGICAL,MPI_COMM_WORLD,error)
        CALL MPI_AlltoAllv(IndexTable(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,Index2Table,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
        CALL MPI_AlltoAllv(ProcessVec(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,Process2Vec,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
!        IF(TLocalAnnihilation) THEN
!!If we are locally annihilating, then we need to take the excitation level of the particle with us.
!!We can send the information to CurrentIC - this is where the final information will be stored, but currently, it is redundant.
!            CALL MPI_AlltoAllv(TempExcitLevel(1:TotWalkersNew),sendcounts,disps,MPI_INTEGER,CurrentIC,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
!        ENDIF
        CALL halt_timer(Comms_Time)
        
!        IF(Iter.eq.DebugIter) THEN
!            WRITE(6,*) "AFTER DIVISION:   - No. on processor is: ",MaxIndex
!            do i=1,MaxIndex
!                WRITE(6,*) HashArray(i),abs(mod(HashArray(i),nProcessors)),Index2Table(i),Process2Vec(i),CurrentSign(i)
!            enddo
!            CALL FLUSH(6)
!        ENDIF

!Now we need to perform the actual annihilation, running through all the particles and calculating which ones want to be annihilated.
        CALL set_timer(Sort_Time,30)
        IF(.not.tAnnihilatebyrange) THEN
!The hashes now need to be sorted again - this time by their number
!This sorting would be redundant if we had initially sorted the hashes by range (ie tAnnihilatebyrange).
!            IF(TLocalAnnihilation) THEN
!If we are locally annihilating, then we need to take the excitation level with us.
!                CALL Sort4I1LLong(MaxIndex,HashArray(1:MaxIndex),Index2Table(1:MaxIndex),Process2Vec(1:MaxIndex),CurrentIC(1:MaxIndex),CurrentSign(1:MaxIndex))
!            ELSE
!                CALL Sort3I1LLong(MaxIndex,HashArray(1:MaxIndex),Index2Table(1:MaxIndex),Process2Vec(1:MaxIndex),CurrentSign(1:MaxIndex))
                CALL Sort4ILong(MaxIndex,HashArray(1:MaxIndex),Index2Table(1:MaxIndex),Process2Vec(1:MaxIndex),CurrentSign(1:MaxIndex))
!            ENDIF
        ELSE
!Here, because we have ordered the hashes initially numerically, we have a set of ordered lists. It is therefore easier to sort them.
!We have to work out how to run sequentially through the hashes, which are a set of nProc seperate ordered lists.
!We would need to have 2*nProc indices, since we will have a set of nProc disjoint ordered sublists.
!SubListInds(1,iProc)=index of current hash from processor iProc
!SubListInds(2,iProc)=index of final hash from processor iProc
!Indices can be obtained from recvcounts and recvdisps - recvcounts(iProc-1) is number of hashes from iProc
!recvdisps(iProc-1) is the displacement to the start of the hashes from iProc
            do i=1,nProcessors-1
                SubListInds(1,i)=recvdisps(i)+1
                SubListInds(2,i)=recvdisps(i+1)
            enddo
            SubListInds(1,nProcessors)=recvdisps(nProcessors)+1
            SubListInds(2,nProcessors)=MaxIndex
!            WRITE(6,*) "SubListInds(1,:) ", SubListInds(1,:)
!            WRITE(6,*) "SubListInds(2,:) ", SubListInds(2,:)
!            WRITE(6,*) "Original hash list is: "
!            do i=1,MaxIndex
!                WRITE(6,*) HashArray(i)
!            enddo
!            WRITE(6,*) "**************"
!Reorder the lists so that they are in numerical order.
            j=1
            do while(j.le.MaxIndex)
                do i=1,nProcessors
                    IF(SubListInds(1,i).le.SubListInds(2,i)) THEN
!This block still has hashes which want to be sorted
                        MinHash=HashArray(SubListInds(1,i))
                        MinProc=i
                        MinInd=SubListInds(1,i)
                        EXIT
                    ENDIF
!                    IF(i.eq.nProcessors) THEN
!                        WRITE(6,*) "ERROR HERE!!"
!                        CALL FLUSH(6)
!                    ENDIF
                enddo
                IF(MinHash.ne.HashCurr) THEN
                    do i=MinProc+1,nProcessors
                        IF((SubListInds(1,i).le.SubListInds(2,i)).and.(HashArray(SubListInds(1,i)).lt.MinHash)) THEN
                            MinHash=HashArray(SubListInds(1,i))
                            MinProc=i
                            MinInd=SubListInds(1,i)
                            IF(MinHash.eq.HashCurr) THEN
                                EXIT
                            ENDIF
                        ENDIF
                    enddo
                ENDIF
!Next smallest hash is MinHash - move the ordered elements into the other array.
                Hash2Array(j)=MinHash
                IndexTable(j)=Index2Table(MinInd)
                ProcessVec(j)=Process2Vec(MinInd)
                NewSign(j)=CurrentSign(MinInd)
                HashCurr=MinHash
!Move through the block
                j=j+1
                SubListInds(1,MinProc)=SubListInds(1,MinProc)+1
            enddo

            IF((j-1).ne.MaxIndex) THEN
                CALL Stop_All(this_routine,"Error here in the merge sort algorithm")
            ENDIF

!Need to copy the lists back to the original array
            do i=1,MaxIndex
                Index2Table(i)=IndexTable(i)
                Process2Vec(i)=ProcessVec(i)
                CurrentSign(i)=NewSign(i)
                HashArray(i)=Hash2Array(i)
            enddo
                
!            Index2Table(1:MaxIndex)=IndexTable(1:MaxIndex)
!            Process2Vec(1:MaxIndex)=ProcessVec(1:MaxIndex)
!            CurrentSign(1:MaxIndex)=NewSign(1:MaxIndex)
!            HashArray(1:MaxIndex)=Hash2Array(1:MaxIndex)
                    
        ENDIF

        CALL halt_timer(Sort_Time)

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
                IF(CurrentSign(j).eq.1) THEN
                    TotWalkersDet=TotWalkersDet+1
                ELSE
                    TotWalkersDet=TotWalkersDet-1
                ENDIF
                FinalBlockIndex=FinalBlockIndex+1
                j=j+1
            enddo

!            IF(TotWalkersDet.gt.AllNoatHF) THEN
!                WRITE(6,"(A,I20,2I7)") "High-weighted Det: ",HashCurr,TotWalkersDet,INT(AllSumNoatHF,4)
!            ENDIF
     
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
                ELSEIF((TotWalkersDet.lt.0).and.(CurrentSign(k).eq.1)) THEN
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
                ELSEIF((TotWalkersDet.gt.0).and.(CurrentSign(k).eq.-1)) THEN
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
                    IF(CurrentSign(k).eq.1) THEN
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
            CALL set_timer(Sort_Time,30)
!Do not actually have to take indextable, hash2array or newsign with it...
            CALL Sort2IILongI(ToAnnihilateIndex,ProcessVec(1:ToAnnihilateIndex),IndexTable(1:ToAnnihilateIndex),Hash2Array(1:ToAnnihilateIndex),NewSign(1:ToAnnihilateIndex))
!            CALL Sort2IILongL(ToAnnihilateIndex,ProcessVec(1:ToAnnihilateIndex),IndexTable(1:ToAnnihilateIndex),Hash2Array(1:ToAnnihilateIndex),NewSign(1:ToAnnihilateIndex))
            CALL halt_timer(Sort_Time)
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
!        CALL MPI_Barrier(MPI_COMM_WORLD,error)
        CALL set_timer(Comms_Time,30)

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
        CALL halt_timer(Comms_Time)

!TEST
!        do i=1,ToAnnihilateonProc
!            IF(Process2Vec(i).ne.(iProcIndex)) THEN
!                CALL Stop_All("AnnihilatePartPar","AlltoAllv performed incorrectly")
!            ENDIF
!        enddo

!Index2Table now is a list, of length "ToAnnihilateonProc", of walkers which should NOT be transferred to the next array. 
!Order the list according to this index (Hash and sign does not need to be sorted, but will for debugging purposes)
        CALL set_timer(Sort_Time,30)
        CALL SORTIILongI(ToAnnihilateonProc,Index2Table(1:ToAnnihilateonProc),HashArray(1:ToAnnihilateonProc),CurrentSign(1:ToAnnihilateonProc))
!        CALL SORTIILongL(ToAnnihilateonProc,Index2Table(1:ToAnnihilateonProc),HashArray(1:ToAnnihilateonProc),CurrentSign(1:ToAnnihilateonProc))
        CALL halt_timer(Sort_Time)

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
!                    CurrentIC(VecSlot)=NewIC(i)
                    IF(.not.tRegenDiagHEls) CurrentH(VecSlot)=NewH(i)
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
!                CurrentIC(VecSlot)=NewIC(i)
                IF(.not.tRegenDiagHEls) CurrentH(VecSlot)=NewH(i)
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
!                CurrentIC(VecSlot)=NewIC(i)
                IF(.not.tRegenDiagHEls) CurrentH(VecSlot)=NewH(i)
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
!AJWT added optional argument dProb which is a probability that whatever gave this contribution as generated.
!  It defaults to 1, and weights the contribution of this det. (Only in the projected energy)
    SUBROUTINE SumEContrib(DetCurr,ExcitLevel,WSign,iLutCurr,HDiagCurr,dProbFin)
        use SystemData, only : tNoBrillouin
        INTEGER :: DetCurr(NEl),ExcitLevel,i,HighIndex,LowIndex,iLutCurr(0:NIfD),WSign,Bin
        INTEGER :: PartInd,iLutSym(0:NIfD),OpenOrbs
        LOGICAL :: CompiPath,tSuccess,iLut2(0:NIfD),DetBitEQ
        REAL*8 :: HDiagCurr
        TYPE(HElement) :: HOffDiag
        real*8 dProbFin
!        write(81,*) DetCurr,ExcitLevel,WSign,iLutCurr,HDiagCurr,dProb

!        MeanExcitLevel=MeanExcitLevel+real(ExcitLevel,r2)
!        IF(MinExcitLevel.gt.ExcitLevel) MinExcitLevel=ExcitLevel
!        IF(MaxExcitLevel.lt.ExcitLevel) MaxExcitLevel=ExcitLevel
!        DetsNorm=DetsNorm+REAL((WSign**2),r2)
        IF(ExcitLevel.eq.0) THEN
            IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF+WSign
            NoatHF=NoatHF+WSign
            HFCyc=HFCyc+WSign      !This is simply the number at HF*sign over the course of the update cycle 
!            AvSign=AvSign+REAL(WSign,r2)
!            AvSignHFD=AvSignHFD+REAL(WSign,r2)
            
        ELSEIF(ExcitLevel.eq.2) THEN
            NoatDoubs=NoatDoubs+abs(WSign)
!At double excit - find and sum in energy
            IF(tHPHF) THEN
                CALL HPHFGetOffDiagHElement(HFDet,DetCurr,iLutHF,iLutCurr,HOffDiag)
            ELSE
                HOffDiag=GetHElement2(HFDet,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
            ENDIF
            IF(Iter.gt.NEquilSteps) SumENum=SumENum+(REAL(HOffDiag%v,r2)*WSign/dProbFin)
!            AvSign=AvSign+REAL(WSign,r2)
!            AvSignHFD=AvSignHFD+REAL(WSign,r2)
            ENumCyc=ENumCyc+(REAL(HOffDiag%v,r2)*WSign/dProbFin)     !This is simply the Hij*sign summed over the course of the update cycle
            
            
!        ELSE
!            AvSign=AvSign+REAL(WSign,r2)

        ELSEIF(ExcitLevel.eq.1) THEN
          if(tNoBrillouin.or.(tHub.and.tReal).or.tRotatedOrbs) then
!For the real-space hubbard model, determinants are only connected to excitations one level away, and brillouins theorem can not hold.
!For Rotated orbitals, brillouins theorem also cannot hold, and energy contributions from walkers on singly excited determinants must
!be included in the energy values along with the doubles.
            IF(tHPHF) THEN
                CALL HPHFGetOffDiagHElement(HFDet,DetCurr,iLutHF,iLutCurr,HOffDiag)
            ELSE
                HOffDiag=GetHElement2(HFDet,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
            ENDIF
            IF(Iter.gt.NEquilSteps) SumENum=SumENum+(REAL(HOffDiag%v,r2)*WSign/dProbFin)
!            AvSign=AvSign+REAL(WSign,r2)
!            AvSignHFD=AvSignHFD+REAL(WSign,r2)
            ENumCyc=ENumCyc+(REAL(HOffDiag%v,r2)*WSign/dProbFin)     !This is simply the Hij*sign summed over the course of the update cycle
          endif 

          IF(tConstructNOs) THEN
!Fill the 1-RDM to find natural orbital on-the-fly.
!Find the orbitals that are involved in the excitation (those that differ in occupation to the ref orbital).
                CALL FindSingleOrbs(iLutHF,iLutCurr,NIfD,Orbs)
!Add 1.D0 (or -1.D0) to the off-diagonal element connecting the relevent orbitals.
                IF(Iter.gt.NEquilSteps) THEN
                    OneRDM(Orbs(1),Orbs(2))=OneRDM(Orbs(1),Orbs(2))+REAL(WSign,r2)
                    OneRDM(Orbs(2),Orbs(1))=OneRDM(Orbs(1),Orbs(2))
                ENDIF
!At the end of all iterations, this OneRDM will contain only the unnormalised off-diagonal elements.
          ENDIF
            
        ENDIF

        
!Histogramming diagnostic options...
        IF(tHistSpawn.or.tCalcFCIMCPsi.and.(Iter.ge.NHistEquilSteps)) THEN
            IF(ExcitLevel.eq.NEl) THEN
                CALL BinSearchParts2(iLutCurr,HistMinInd(ExcitLevel),Det,PartInd,tSuccess)
                HistMinInd(ExcitLevel)=PartInd
            ELSEIF(ExcitLevel.eq.0) THEN
                PartInd=1
                tSuccess=.true.
            ELSE
                CALL BinSearchParts2(iLutCurr,HistMinInd(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartInd,tSuccess)
                HistMinInd(ExcitLevel)=PartInd
            ENDIF
            IF(tSuccess) THEN
                IF(tFlippedSign) THEN
                    Histogram(PartInd)=Histogram(PartInd)-REAL(WSign,r2)
                    IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-REAL(WSign,r2)
                ELSE
                    Histogram(PartInd)=Histogram(PartInd)+REAL(WSign,r2)
                    IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+REAL(WSign,r2)
                ENDIF
            ELSE
                WRITE(6,*) DetCurr(:)
                WRITE(6,*) "***",iLutCurr(0:NIfD)
                WRITE(6,*) "***",ExcitLevel,HistMinInd(ExcitLevel),Det
                Call WriteBitDet(6,iLutCurr(0:NIfD),.true.)
                CALL Stop_All("SumEContrib","Cannot find corresponding FCI determinant when histogramming")
            ENDIF
            IF(tHPHF) THEN
!With HPHF space, we need to also include the spin-coupled determinant, which will have the same amplitude as the original determinant, unless it is antisymmetric.
                CALL FindExcitBitDetSym(iLutCurr,iLutSym)
                IF(ExcitLevel.eq.NEl) THEN
                    CALL BinSearchParts2(iLutSym,FCIDetIndex(ExcitLevel),Det,PartInd,tSuccess)
                ELSEIF(ExcitLevel.eq.0) THEN
                    PartInd=1
                    tSuccess=.true.
                ELSE
                    CALL BinSearchParts2(iLutSym,FCIDetIndex(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartInd,tSuccess)
                ENDIF
                IF(tSuccess) THEN
                    CALL CalcOpenOrbs(iLutSym,NIfD,NEl,OpenOrbs)
                    IF(tFlippedSign) THEN
                        IF(mod(OpenOrbs,2).eq.1) THEN
                            Histogram(PartInd)=Histogram(PartInd)+REAL(WSign,r2)
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+REAL(WSign,r2)
                        ELSE
                            Histogram(PartInd)=Histogram(PartInd)-REAL(WSign,r2)
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-REAL(WSign,r2)
                        ENDIF
                    ELSE
                        IF(mod(OpenOrbs,2).eq.1) THEN
                            Histogram(PartInd)=Histogram(PartInd)-REAL(WSign,r2)
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-REAL(WSign,r2)
                        ELSE
                            Histogram(PartInd)=Histogram(PartInd)+REAL(WSign,r2)
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+REAL(WSign,r2)
                        ENDIF
                    ENDIF
                ELSE
                    WRITE(6,*) DetCurr(:)
                    WRITE(6,*) "***",iLutSym(0:NIfD)
                    WRITE(6,*) "***",ExcitLevel,Det
                    CALL Stop_All("SumEContrib","Cannot find corresponding spin-coupled FCI determinant when histogramming")
                ENDIF
            ENDIF
        ELSEIF(tHistEnergies) THEN
!This wil histogramm the energies of the particles, rather than the determinants themselves.
            Bin=INT(HDiagCurr/BinRange)+1
            IF(Bin.gt.iNoBins) THEN
                CALL Stop_All("SumEContrib","Histogramming energies higher than the arrays can cope with. Increase iNoBins or BinRange")
            ENDIF
            Histogram(Bin)=Histogram(Bin)+real(abs(WSign),r2)
        ENDIF


!        IF(TLocalAnnihilation) THEN
!We need to count the population of walkers at each excitation level for each iteration
!            PartsinExcitLevel(ExcitLevel)=PartsinExcitLevel(ExcitLevel)+1
!        ENDIF

!        IF(TAutoCorr) THEN
!!First element is HF Det to histogram. Then come doubles, triples and quads
!
!            IF(ExcitLevel.eq.4) THEN
!                LowIndex=2+NoACDets(2)+NoACDets(3)
!                HighIndex=NoAutoDets
!            ELSEIF(ExcitLevel.eq.3) THEN
!                LowIndex=2+NoACDets(2)
!                HighIndex=1+NoACDets(2)+NoACDets(3)
!            ELSEIF(ExcitLevel.eq.2) THEN
!                LowIndex=2
!                HighIndex=1+NoACDets(2)
!            ELSEIF(ExcitLevel.eq.0) THEN
!                LowIndex=1
!                HighIndex=1
!            ELSE
!                LowIndex=0
!            ENDIF
!
!            IF(LowIndex.ne.0) THEN
!
!                do i=LowIndex,HighIndex
!            
!                    IF(CompiPath(DetCurr,AutoCorrDets(:,i),NEl)) THEN
!!The walker is at a determinant for which we want to calculate the autocorrelation function
!                        IF(TFlippedSign) THEN
!                            WeightatDets(i)=WeightatDets(i)-WSign
!                        ELSE
!                            WeightatDets(i)=WeightatDets(i)+WSign
!                        ENDIF
!                        EXIT
!                    ENDIF
!
!                enddo
!            ENDIF
!        ENDIF
        
        RETURN

    END SUBROUTINE SumEContrib

    
!Every StepsSft steps, update the diagonal shift value (the running value for the correlation energy)
!We don't want to do this too often, since we want the population levels to acclimatise between changing the shifts
    SUBROUTINE CalcNewShift()
        USE FciMCLoggingMOD , only : PrintSpawnAttemptStats,PrintTriConnStats,PrintSpinCoupHEl,InitErrorBlocking,SumInErrorContrib
        USE FciMCLoggingMOD , only : InitShiftErrorBlocking,SumInShiftErrorContrib
        INTEGER :: error,rc,MaxAllowedWalkers,MaxWalkersProc,MinWalkersProc
        INTEGER :: inpair(9),outpair(9)
        REAL*8 :: TempTotWalkers,TempTotParts
        REAL*8 :: TempSumNoatHF,MeanWalkers,TempSumWalkersCyc,TempAllSumWalkersCyc,TempNoMinorWalkers
        REAL*8 :: inpairreal(3),outpairreal(3)
        LOGICAL :: TBalanceNodesTemp
        
        IF(TSinglePartPhase) THEN
!Exit the single particle phase if the number of walkers exceeds the value in the input file.
!            CALL MPI_Barrier(MPI_COMM_WORLD,error)
            IF(iProcIndex.eq.root) THEN     !Only exit phase if particle number is sufficient on head node.
                IF(TotParts.gt.InitWalkers) THEN
                    WRITE(6,*) "Exiting the single particle growth phase - shift can now change"
                    VaryShiftIter=Iter
                    TSinglePartPhase=.false.
                ENDIF
            ENDIF
!Broadcast the fact that TSinglePartPhase may have changed to all processors - unfortunatly, have to do this broadcast every iteration.
            CALL MPI_Bcast(TSinglePartPhase,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
        ENDIF

!This first call will calculate the GrowRate for each processor, taking culling into account
!        WRITE(6,*) "Get Here"
!        CALL FLUSH(6)
        CALL UpdateDiagSftPar()

!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

!Find the total number of particles at HF (x sign) across all nodes. If this is negative, flip the sign of all particles.
        AllNoatHF=0

!Find sum of noathf, and then use an AllReduce to broadcast it to all nodes
        CALL MPI_AllReduce(NoatHF,AllNoatHF,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,error)

        IF(AllNoatHF.lt.0) THEN
!Flip the sign if we're beginning to get a negative population on the HF
            WRITE(6,*) "No. at HF < 0 - flipping sign of entire ensemble of particles..."
            CALL FlipSign()
        ENDIF
                
!IF we're using a guiding function, want to sum in the contributions to the energy from this guiding function.
        IF(tUseGuide) THEN
!These two are not zeroed after each update cycle, so are average energies over the whole simulation
            SumENum=SumENum+GuideFuncDoub
            SumNoatHF=SumNoatHF+GuideFuncHF
!These two variables are zeroed after each update cycle, so are the "instantaneous" energy (averaged only over the update cycle)
            HFCyc=HFCyc+GuideFuncHF
            ENumCyc=ENumCyc+GuideFuncDoub
        ENDIF

!We need to collate the information from the different processors
!Inpair and outpair are used to package variables to save on latency time
!        inpair(1)=TotWalkers
        inpair(1)=Annihilated
        inpair(2)=NoatDoubs
        inpair(3)=NoBorn
        inpair(4)=NoDied
        inpair(5)=HFCyc         !SumWalkersCyc is now an integer*8
        inpair(6)=LocalAnn
        inpair(7)=SpawnFromSing
        inpair(8)=iInitGuideParts
        inpair(9)=MinorAnnihilated
!        inpair(7)=TotParts
!        inpair(9)=iUniqueDets
        outpair(:)=0
!        CALL MPI_Reduce(TotWalkers,AllTotWalkers,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!Find total Annihilated,Total at HF and Total at doubles
!        CALL MPI_Reduce(Annihilated,AllAnnihilated,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!!        CALL MPI_Reduce(NoatHF,AllNoatHF,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)  !This is done every iteration now
!        CALL MPI_Reduce(NoatDoubs,AllNoatDoubs,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(NoBorn,AllNoBorn,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(NoDied,AllNoDied,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(SumWalkersCyc,AllSumWalkersCyc,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get Here 1"
!        CALL FLUSH(6)
        CALL MPI_Reduce(inpair,outpair,9,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get Here 2"
!        CALL FLUSH(6)
!        AllTotWalkers=outpair(1)
        AllAnnihilated=outpair(1)
        AllNoatDoubs=outpair(2)
        AllNoBorn=outpair(3)
        AllNoDied=outpair(4)
        AllHFCyc=REAL(outpair(5),r2)
        AllLocalAnn=outpair(6)
        AllSpawnFromSing=outpair(7)
!        AllTotParts=outpair(7)
!        AlliUniqueDets=REAL(outpair(9),r2)
        AlliInitGuideParts=outpair(8)
        AllMinorAnnihilated=outpair(9)
        TempTotWalkers=REAL(TotWalkers,r2)
        TempTotParts=REAL(TotParts,r2)
        TempNoMinorWalkers=REAL(NoMinorWalkers,r2)
        IF(tMinorDetsStar) THEN
            TempTotParts=TempTotParts+TempNoMinorWalkers
        ENDIF

        CALL MPI_Reduce(TempTotWalkers,AllTotWalkers,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(TempTotParts,AllTotParts,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(TempNoMinorWalkers,AllNoMinorWalkers,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)

        IF(iProcIndex.eq.0) THEN
            IF(AllTotWalkers.le.0.2) THEN
                WRITE(6,*) AllTotWalkers,TotWalkers
                CALL Stop_All("CalcNewShift","All particles have died. Consider choosing new seed, or raising shift value.")
            ENDIF
        ENDIF

!SumWalkersCyc is now an int*8, therefore is needs to be reduced as a real*8
        TempSumWalkersCyc=REAL(SumWalkersCyc,r2)
        TempAllSumWalkersCyc=0.D0
        CALL MPI_Reduce(TempSumWalkersCyc,TempAllSumWalkersCyc,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get Here 3"
!        CALL FLUSH(6)

        IF(.not.tDirectAnnihil) THEN

            MeanWalkers=AllTotWalkers/REAL(nProcessors,r2)
            MaxAllowedWalkers=NINT((MeanWalkers/12.D0)+MeanWalkers)

!Find the range of walkers on different nodes to see if we need to even up the distribution over nodes
!            inpair(1)=TotWalkers
!            inpair(2)=iProcIndex
            CALL MPI_Reduce(TotWalkers,MaxWalkersProc,1,MPI_INTEGER,MPI_MAX,root,MPI_COMM_WORLD,error)
            CALL MPI_Reduce(TotWalkers,MinWalkersProc,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,error)
            IF(iProcIndex.eq.Root) THEN
                WalkersDiffProc=MaxWalkersProc-MinWalkersProc
            ENDIF
!            WRITE(6,*) "Get Here 4"
!            CALL FLUSH(6)
!            MaxWalkersProc=outpair(1)
!            WRITE(6,*) "***",MaxWalkersProc,MaxAllowedWalkers,MeanWalkers
!            CALL MPI_Reduce(inpair,outpair,1,MPI_2INTEGER,MPI_MINLOC,root,MPI_COMM_WORLD,error)
!            MinWalkersProc=outpair(1)

            IF(iProcIndex.eq.root) THEN
!                RangeWalkers=MaxWalkersProc-MinWalkersProc
!                IF(RangeWalkers.gt.300) THEN
                IF((MaxWalkersProc.gt.MaxAllowedWalkers).and.(AllTotWalkers.gt.(REAL(nProcessors*500,r2)))) THEN
                    TBalanceNodesTemp=.true.
                ELSE
                    TBalanceNodesTemp=.false.
                ENDIF
            ENDIF
!Also choose to balance the nodes if all particles have died on one of them
            IF(TotWalkers.eq.0) THEN
                TBalanceNodesTemp=.true.
            ELSE
                IF(iProcIndex.ne.Root) THEN
                    TBalanceNodesTemp=.false.
                ENDIF
            ENDIF
!We need to tell all nodes whether to balance the nodes or not...
            CALL MPI_AllReduce(TBalanceNodesTemp,TBalanceNodes,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,error)
!            WRITE(6,*) "Get Here 5"
!            CALL FLUSH(6)
!            CALL MPI_BCast(TBalanceNodes,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
            IF(iProcIndex.eq.Root) THEN
                IF(TBalanceNodes.and.(.not.TBalanceNodesTemp)) THEN
                    WRITE(6,*) "Balancing nodes since all particles have died on a node..."
                ENDIF
            ENDIF

        ELSE
!Cannot load-balance with direct annihilation, but still want max & min
            CALL MPI_Reduce(TotWalkers,MaxWalkersProc,1,MPI_INTEGER,MPI_MAX,root,MPI_COMM_WORLD,error)
            CALL MPI_Reduce(TotWalkers,MinWalkersProc,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,error)
            IF(iProcIndex.eq.Root) THEN
                WalkersDiffProc=MaxWalkersProc-MinWalkersProc
            ENDIF

!            TBalanceNodes=.false.   !Temporarily turn off node balancing
        ENDIF

!AlliUniqueDets corresponds to the total number of unique determinants, summed over all iterations in the last update cycle, and over all processors.
!Divide by StepsSft to get the average number of unique determinants visited over a single iteration.
!        AlliUniqueDets=AlliUniqueDets/(REAL(StepsSft,r2))
        
        IF(GrowRate.eq.-1.D0) THEN
!tGlobalSftCng is on, and we want to calculate the change in the shift as a global parameter, rather than as a weighted average.
!This will only be a sensible value on the root.
            AllGrowRate=AllTotParts/AllTotPartsOld
        ELSE
!We want to calculate the mean growth rate over the update cycle, weighted by the total number of walkers
            GrowRate=GrowRate*TempSumWalkersCyc                    
            CALL MPIDSumRoot(GrowRate,1,AllGrowRate,Root)   

            IF(iProcIndex.eq.Root) THEN
                AllGrowRate=AllGrowRate/TempAllSumWalkersCyc
            ENDIF
        ENDIF
!        WRITE(6,*) "Get Here 6"
!        CALL FLUSH(6)

        IterTime=IterTime/REAL(StepsSft)    !This is the average time per iteration in the previous update cycle.

!For the unweighted by iterations energy estimator (ProjEIter), we need the sum of the Hij*Sign from all processors over the last update cycle
!        CALL MPIDSumRoot(ENumCyc,1,AllENumCyc,Root)
!        WRITE(6,*) "Get Here 7"
!        CALL FLUSH(6)

!Do the same for the mean excitation level of all walkers, and the total positive particles
!MeanExcitLevel here is just the sum of all the excitation levels - it needs to be divided by the total walkers in the update cycle first.
!        WRITE(6,"(2I10,2G25.16)",advance='no') Iter,TotWalkers,MeanExcitLevel,TempSumWalkersCyc
!        MeanExcitLevel=MeanExcitLevel/TempSumWalkersCyc
!        CALL MPI_Reduce(MeanExcitLevel,AllMeanExcitLevel,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get Here 8"
!        CALL FLUSH(6)
!        CALL MPIDSumRoot(MeanExcitLevel,1,AllMeanExcitLevel,Root)
!        IF(iProcIndex.eq.Root) THEN
!            AllMeanExcitLevel=AllMeanExcitLevel/real(nProcessors,r2)
!        ENDIF

!AvSign no longer calculated (but would be easy to put back in) - ACF much better bet...
!        AvSign=AvSign/real(SumWalkersCyc,r2)
!        AvSignHFD=AvSignHFD/real(SumWalkersCyc,r2)
!        CALL MPI_Reduce(AvSign,AllAvSign,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(AvSignHFD,AllAvSignHFD,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        IF(iProcIndex.eq.Root) THEN
!            AllAvSign=AllAvSign/real(nProcessors,r2)
!            AllAvSignHFD=AllAvSignHFD/real(nProcessors,r2)
!        ENDIF

!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        TempSumNoatHF=real(SumNoatHF,r2)
!        CALL MPIDSumRoot(TempSumNoatHF,1,AllSumNoatHF,Root)
!        WRITE(6,*) "Get Here 9"
!        CALL FLUSH(6)
!        CALL MPIDSumRoot(SumENum,1,AllSumENum,Root)
!        WRITE(6,*) "Get Here 10"
!        CALL FLUSH(6)
        inpairreal(1)=ENumCyc
        inpairreal(2)=TempSumNoatHF
        inpairreal(3)=SumENum
!        inpairreal(4)=DetsNorm
        CALL MPI_Reduce(inpairreal,outpairreal,3,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        AllENumCyc=outpairreal(1)
        AllSumNoatHF=outpairreal(2)
        AllSumENum=outpairreal(3)
!        AllDetsNorm=outpairreal(4)


!To find minimum and maximum excitation levels, search for them using MPI_Reduce
!        inpair(1)=MaxExcitLevel
!        inpair(2)=iProcIndex

!        CALL MPI_Reduce(MaxExcitLevel,AllMaxExcitLevel,1,MPI_INTEGER,MPI_MAX,Root,MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get Here 11"
!        CALL FLUSH(6)
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
!        CALL MPI_Reduce(MinExcitLevel,AllMinExcitLevel,1,MPI_INTEGER,MPI_MIN,Root,MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get Here 12"
!        CALL FLUSH(6)
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in finding min excitation level"
!            CALL MPI_ABORT(MPI_COMM_WORLD,rc,error)
!        ENDIF
!        IF(iProcIndex.eq.Root) THEN
!            AllMinExcitLevel=outpair(1)
!        ENDIF

!We now want to find how the shift should change for the entire ensemble of processors
        IF(iProcIndex.eq.Root) THEN
            IF(.not.TSinglePartPhase) THEN
                DiagSft=DiagSft-(log(AllGrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
                IF((Iter-VaryShiftIter).gt.NEquilSteps) THEN
!                    WRITE(6,*) Iter-VaryShiftIter, NEquilSteps*StepsSft
                    VaryShiftCycles=VaryShiftCycles+1
                    SumDiagSft=SumDiagSft+DiagSft
                    AvDiagSft=SumDiagSft/REAL(VaryShiftCycles,r2)
                ENDIF
            ENDIF

            IF(AllSumNoatHF.ne.0.D0) THEN
!AllSumNoatHF can actually be 0 if we have equilsteps on.
                ProjectionE=AllSumENum/AllSumNoatHF
            ENDIF

!Calculate the projected energy where each update cycle contributes the same weight to the average for its estimator for the energy
            IF(AllHFCyc.ne.0.D0) THEN
                ProjEIterSum=ProjEIterSum+(AllENumCyc/AllHFCyc)
                HFPopCyc=HFPopCyc+1   !This is the number of iterations where we have a non-zero contribution from HF particles
                ProjEIter=ProjEIterSum/REAL(HFPopCyc,r2)
            ENDIF
        ENDIF
!        IF(tHub.and.tReal) THEN
!!Since for the real-space hubbard model the reference is not the HF, it has to be added on to the energy since it is not subtracted from the
!!diagonal hamiltonian elements.
!            ProjectionE=ProjectionE+HubRefEnergy
!            ProjEIter=ProjEIter+HubRefEnergy
!        ENDIF
!We wan to now broadcast this new shift to all processors
        CALL MPI_Bcast(DiagSft,1,MPI_DOUBLE_PRECISION,Root,MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get Here 13"
!        CALL FLUSH(6)
!        IF(error.ne.MPI_SUCCESS) THEN
!            WRITE(6,*) "Error in broadcasting new shift"
!            CALL MPI_ABORT(MPI_COMM_WORLD,rc,error)
!        ENDIF

        AccRat=(REAL(Acceptances,r2))/TempSumWalkersCyc      !The acceptance ratio which is printed is only for the current node - not summed over all nodes

        CALL WriteFCIMCStats()


!This first bit checks if it is time to set up the blocking analysis.  This is obviously only done once, so these logicals become false once it is done. 
        IF(iProcIndex.eq.Root) THEN
            IF(tIterStartBlock) THEN
                IF(Iter.ge.IterStartBlocking) THEN 
                    CALL InitErrorBlocking(Iter)
                    tIterStartBlock=.false.
                    tErrorBlocking=.true.
                ENDIF
            ELSEIF(tHFPopStartBlock) THEN
                IF((AllHFCyc/StepsSft).ge.HFPopStartBlocking) THEN
                    CALL InitErrorBlocking(Iter)
                    tHFPopStartBlock=.false.
                    tErrorBlocking=.true.
                ENDIF
            ENDIF

            IF((.not.TSinglePartPhase).and.tInitShiftBlocking) THEN
                CALL InitShiftErrorBlocking(Iter)
                tInitShiftBlocking=.false.
                tShiftBlocking=.true.
            ENDIF

!Then we perform the blocking at the end of each update cycle.         
            IF(tErrorBlocking) CALL SumInErrorContrib(Iter,AllENumCyc,AllHFCyc)
            IF(tShiftBlocking) CALL SumInShiftErrorContrib(Iter,DiagSft)
        ENDIF

        IF(tPrintTriConnections) CALL PrintTriConnStats(Iter+PreviousCycles)
        IF(tPrintSpinCoupHEl) CALL PrintSpinCoupHEl(Iter+PreviousCycles)

        IF(tPrintHElAccept) CALL PrintSpawnAttemptStats(Iter+PreviousCycles)

!Now need to reinitialise all variables on all processers
        IterTime=0.0
!        MinExcitLevel=NEl+10
!        MaxExcitLevel=0
!        MeanExcitLevel=0.D0
        SumWalkersCyc=0
!        AvSign=0.D0        !Rezero this quantity - <s> is now a average over the update cycle
!        AvSignHFD=0.D0     !This is the average sign over the HF and doubles
!        DetsNorm=0.D0
        Annihilated=0
        MinorAnnihilated=0
        LocalAnn=0
        Acceptances=0
        NoBorn=0
        SpawnFromSing=0
        NoDied=0
        ENumCyc=0.D0
!        ProjEIter=0.D0     Do not want to rezero, since otherwise, if there are no particles at HF in the next update cycle, it will print out zero.
        HFCyc=0
        GuideFuncHF=0
        GuideFuncDoub=0.D0
!Reset TotWalkersOld so that it is the number of walkers now
        TotWalkersOld=TotWalkers
        TotPartsOld=TotParts

!Also reinitialise the global variables - should not necessarily need to do this...
!        AllHFCyc=0.D0
!        AllENumCyc=0.D0
!        AllDetsNorm=0.D0
        AllSumENum=0.D0
        AllSumNoatHF=0.D0
        AllTotWalkersOld=AllTotWalkers
        AllTotPartsOld=AllTotParts
        AllTotWalkers=0.D0
        AllTotParts=0.D0
        AllGrowRate=0.D0
!        AllMeanExcitLevel=0.D0
!        AllAvSign=0.D0
!        AllAvSignHFD=0.D0
        AllSumWalkersCyc=0
        AllAnnihilated=0
        AllMinorAnnihilated=0
        AllLocalAnn=0
        AllNoatHF=0
        AllNoatDoubs=0
        AllNoBorn=0
        AllSpawnFromSing=0
        AllNoDied=0

        RETURN
    END SUBROUTINE CalcNewShift


!This routine acts as a thermostat for the simulation - killing random particles if the population becomes too large, or 
!Doubling them if it gets too low...
    SUBROUTINE ThermostatParticlesPar(HighLow)
        IMPLICIT NONE
        LOGICAL :: HighLow
        INTEGER :: VecSlot,i,j,ToCull,Culled,OrigWalkers,Chosen
        REAL*8 :: r

        IF(HighLow) THEN
!The population is too large - cull TotWalkers/CullFactor randomly selected particles

            OrigWalkers=TotWalkers
            ToCull=TotWalkers-nint((TotWalkers+0.D0)/CullFactor)
            Culled=0

            do while (Culled.lt.ToCull)

!Pick a random walker between 1 and TotWalkers
                IF(tMerTwist) THEN
                    CALL genrand_real2(r) 
                ELSE
                    CALL RANLUX(r,1)
                ENDIF
                Chosen=int((r*TotWalkers)+1.D0)

!Move the Walker at the end of the list to the position of the walker we have chosen to destroy
                CurrentDets(:,Chosen)=CurrentDets(:,TotWalkers)
                CurrentSign(Chosen)=CurrentSign(TotWalkers)
                IF(.not.tRegenDiagHEls) CurrentH(Chosen)=CurrentH(TotWalkers)
!                CurrentIC(Chosen)=CurrentIC(TotWalkers)
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

            IF(tRotoAnnihil) THEN

                do i=1,TotWalkers
                    CurrentSign(i)=CurrentSign(i)*2
                enddo
                TotParts=TotParts*2
                CullInfo(NoCulls,2)=TotParts

            ELSE

                VecSlot=TotWalkers+1
                do i=1,TotWalkers

!Add clone of walker, at the same determinant, to the end of the list
                    CurrentDets(:,VecSlot)=CurrentDets(:,i)
                    CurrentSign(VecSlot)=CurrentSign(i)
                    IF(.not.tRegenDiagHEls) CurrentH(VecSlot)=CurrentH(i)
!                    CurrentIC(VecSlot)=CurrentIC(i)
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

        ENDIF

        RETURN

    END SUBROUTINE ThermostatParticlesPar


!This routine looks at the change in residual particle number over a number of cycles, and adjusts the 
!value of the diagonal shift in the hamiltonian in order to compensate for this
    SUBROUTINE UpdateDiagSftPar()
        USE CalcData , only : tGlobalSftCng
        INTEGER :: j,k,GrowthSteps,MaxCulls,error
        LOGICAL :: Changed

        Changed=.false.
        IF(tGlobalSftCng) THEN
            CALL MPI_AllReduce(NoCulls,MaxCulls,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,error)
            IF(MaxCulls.gt.0) THEN
                IF(iProcIndex.eq.0) WRITE(6,*) "Culling has occurred in this update cycle..."
!At least one of the nodes is culling at least once, therefore every processor has to perform the original grow rate calculation.
                tGlobalSftCng=.false.
                Changed=.true.
            ENDIF
        ENDIF

        IF(NoCulls.eq.0) THEN
            IF(.not.tGlobalSftCng) THEN
                IF(TotWalkersOld.eq.0) THEN
                    GrowRate=0.D0
                ELSE
                    GrowRate=(TotWalkers+0.D0)/(TotWalkersOld+0.D0)
                ENDIF
            ELSE
                GrowRate=-1.D0
            ENDIF
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

        IF(Changed) THEN
!Return the flag for global shift change back to true.
            tGlobalSftCng=.true.
        ENDIF
        
!        DiagSft=DiagSft-(log(GrowRate)*SftDamp)/(Tau*(StepsSft+0.D0))
!        IF((DiagSft).gt.0.D0) THEN
!            WRITE(6,*) "***WARNING*** - DiagSft trying to become positive..."
!            STOP
!        ENDIF

    END SUBROUTINE UpdateDiagSftPar
        
    


!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalcPar()
        use FciMCLoggingMOD , only : InitTriHElStats,InitSpinCoupHel
        use SystemData , only : tRotateOrbs
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet
        INTEGER :: DetLT,VecSlot,error,MemoryAlloc,Proc
        TYPE(HElement) :: rh,TempHii
        LOGICAL :: exists
        REAL*8 :: TotDets
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMCPar'


        IF(TStartMP1) THEN
!Start the initial distribution off at the distribution of the MP1 eigenvector

            WRITE(6,"(A)") "Starting run with particles populating double excitations proportionally to MP1 wavevector..."
            IF(tDirectAnnihil) CALL Stop_All(this_routine,"MP1Start currently disabled with directannihilation")
            CALL InitWalkersMP1Par()

        ELSEIF(TReadPops) THEN
!Read in particles from multiple POPSFILES for each processor
            WRITE(6,*) "Reading in initial particle configuration from POPSFILES..."

            IF(tDirectAnnihil) CALL Stop_All(this_routine,"READPOPS currently disabled with rotoannihilation")
            CALL ReadFromPopsFilePar()

        ELSE
!initialise the particle positions - start at HF with positive sign

!Set the maximum number of walkers allowed
            MaxWalkersPart=NINT(MemoryFacPart*InitWalkers)
            WRITE(6,"(A,F14.5)") "Memory Factor for walkers is: ",MemoryFacPart
            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node of: ",MaxWalkersPart
            IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
                WRITE(6,"(A,F14.5)") "Memory Factor for arrays used for spawning is: ",MemoryFacSpawn
                WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for spawning of: ",MaxSpawned
            ELSE
                MaxWalkersAnnihil=NINT(MemoryFacAnnihil*InitWalkers)
                WRITE(6,"(A,F14.5)") "Memory Factor for arrays used for annihilation is: ",MemoryFacAnnihil
                WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for annihilation of: ",MaxWalkersAnnihil
            ENDIF

!Put a barrier here so all processes synchronise
            CALL MPI_Barrier(MPI_COMM_WORLD,error)
!Allocate memory to hold walkers
            ALLOCATE(WalkVecDets(0:NIfD,MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfD+1),4,this_routine,WalkVecDetsTag,ierr)
            WalkVecDets(0:NIfD,1:MaxWalkersPart)=0
            IF((.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) THEN
!Rotoannihilation only used a single main array. Spawned particles are put into the spawned arrays.
                ALLOCATE(WalkVec2Dets(0:NIfD,MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVec2Dets',MaxWalkersPart*(NIfD+1),4,this_routine,WalkVec2DetsTag,ierr)
                WalkVec2Dets(0:NIfD,1:MaxWalkersPart)=0
            ENDIF

            IF(.not.tRegenDiagHEls) THEN
                ALLOCATE(WalkVecH(MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVecH',MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
                WalkVecH(:)=0.d0
                IF((.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) THEN
                    ALLOCATE(WalkVec2H(MaxWalkersPart),stat=ierr)
                    CALL LogMemAlloc('WalkVec2H',MaxWalkersPart,8,this_routine,WalkVec2HTag,ierr)
                    WalkVec2H(:)=0.d0
                ENDIF
            ELSE
                IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                    WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*8,r2)/1048576.D0," Mb/Processor"
                ELSE
                    WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*16,r2)/1048576.D0," Mb/Processor"
                ENDIF
            ENDIF
            
            IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                ALLOCATE(WalkVecSign(MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVecSign',MaxWalkersPart,4,this_routine,WalkVecSignTag,ierr)
!                ALLOCATE(WalkVec2Sign(MaxWalkersPart),stat=ierr)
!                CALL LogMemAlloc('WalkVec2Sign',MaxWalkersPart,4,this_routine,WalkVec2SignTag,ierr)
                MemoryAlloc=(NIfD+1+3)*MaxWalkersPart*4    !Memory Allocated in bytes
            ELSE
!The sign is sent through when annihilating, so it needs to be longer.
                ALLOCATE(WalkVecSign(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('WalkVecSign',MaxWalkersAnnihil,4,this_routine,WalkVecSignTag,ierr)
                ALLOCATE(WalkVec2Sign(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('WalkVec2Sign',MaxWalkersAnnihil,4,this_routine,WalkVec2SignTag,ierr)
                MemoryAlloc=((2*MaxWalkersAnnihil)+(((2*(NIfD+1))+4)*MaxWalkersPart))*4    !Memory Allocated in bytes
            ENDIF

            IF(tRegenDiagHEls) THEN
                IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                    MemoryAlloc=MemoryAlloc-(MaxWalkersPart*8)
                ELSE
                    MemoryAlloc=MemoryAlloc-(MaxWalkersPart*16)
                ENDIF
            ENDIF
            
            WalkVecSign(:)=0
            IF((.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) WalkVec2Sign(:)=0

            IF(tRotoAnnihil.or.tDirectAnnihil) THEN

                WRITE(6,"(A,I12,A)") "Spawning vectors allowing for a total of ",MaxSpawned," particles to be spawned in any one iteration."
                ALLOCATE(SpawnVec(0:NIfD,MaxSpawned),stat=ierr)
                CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfD+1),4,this_routine,SpawnVecTag,ierr)
                SpawnVec(:,:)=0
                ALLOCATE(SpawnVec2(0:NIfD,MaxSpawned),stat=ierr)
                CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfD+1),4,this_routine,SpawnVec2Tag,ierr)
                SpawnVec2(:,:)=0
                ALLOCATE(SpawnSignVec(0:MaxSpawned),stat=ierr)
                CALL LogMemAlloc('SpawnSignVec',MaxSpawned+1,4,this_routine,SpawnSignVecTag,ierr)
                SpawnSignVec(:)=0
                ALLOCATE(SpawnSignVec2(0:MaxSpawned),stat=ierr)
                CALL LogMemAlloc('SpawnSignVec2',MaxSpawned+1,4,this_routine,SpawnSignVec2Tag,ierr)
                SpawnSignVec2(:)=0

!Point at correct spawning arrays
                SpawnedParts=>SpawnVec
                SpawnedParts2=>SpawnVec2
                SpawnedSign=>SpawnSignVec
                SpawnedSign2=>SpawnSignVec2

                MemoryAlloc=MemoryAlloc+(((MaxSpawned+1)*2)+(2*MaxSpawned*(1+NIfD)))*4

            ELSEIF(.not.TNoAnnihil) THEN
!            IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                ALLOCATE(HashArray(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('HashArray',MaxWalkersAnnihil,8,this_routine,HashArrayTag,ierr)
                HashArray(:)=0
                ALLOCATE(Hash2Array(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('Hash2Array',MaxWalkersAnnihil,8,this_routine,Hash2ArrayTag,ierr)
                Hash2Array(:)=0
                ALLOCATE(IndexTable(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('IndexTable',MaxWalkersAnnihil,4,this_routine,IndexTableTag,ierr)
                IndexTable(1:MaxWalkersAnnihil)=0
                ALLOCATE(Index2Table(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('Index2Table',MaxWalkersAnnihil,4,this_routine,Index2TableTag,ierr)
                Index2Table(1:MaxWalkersAnnihil)=0
                ALLOCATE(ProcessVec(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('ProcessVec',MaxWalkersAnnihil,4,this_routine,ProcessVecTag,ierr)
                ProcessVec(1:MaxWalkersAnnihil)=0
                ALLOCATE(Process2Vec(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('Process2Vec',MaxWalkersAnnihil,4,this_routine,Process2VecTag,ierr)
                Process2Vec(1:MaxWalkersAnnihil)=0

                MemoryAlloc=MemoryAlloc+32*MaxWalkersAnnihil
            ENDIF

!Allocate pointers to the correct walker arrays
            CurrentDets=>WalkVecDets
            CurrentSign=>WalkVecSign
            IF(.not.tRegenDiagHEls) THEN
                CurrentH=>WalkVecH
                IF((.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) NewH=>WalkVec2H
            ENDIF
            IF((.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) THEN
                NewDets=>WalkVec2Dets
                NewSign=>WalkVec2Sign
            ENDIF

            IF(TStartSinglePart) THEN
                IF(tDirectAnnihil) THEN
                    Proc=DetermineDetProc(iLutHF)
                    WRITE(6,*) "HF processor is: ",Proc
                    IF(iProcIndex.eq.Proc) THEN
                        CurrentDets(:,1)=iLutHF(:)
                        CurrentSign(1)=1
                        IF(.not.tRegenDiagHEls) CurrentH(1)=0.D0
                    ENDIF
                ELSE
                    CurrentDets(:,1)=iLutHF(:)
                    CurrentSign(1)=1
                    IF(.not.tRegenDiagHEls) CurrentH(1)=0.D0
!                    IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    IF((.not.TNoAnnihil).and.(.not.tRotoAnnihil)) THEN
                        HashArray(1)=HFHash
                    ENDIF
                ENDIF
            ELSE

                IF(tRotoAnnihil) THEN
                    CurrentDets(:,1)=iLutHF(:)
                    CurrentSign(1)=InitWalkers
                    IF(.not.tRegenDiagHEls) CurrentH(1)=0.D0
                ELSEIF(tDirectAnnihil) THEN
                    Proc=DetermineDetProc(iLutHF)
                    WRITE(6,*) "HF processor is: ",Proc
                    IF(iProcIndex.eq.Proc) THEN
                        CurrentDets(:,1)=iLutHF(:)
                        CurrentSign(1)=InitWalkers
                        IF(.not.tRegenDiagHEls) CurrentH(1)=0.D0
                    ENDIF
                ELSE

                    do j=1,InitWalkers
                        CurrentDets(:,j)=iLutHF(:)
                        CurrentSign(j)=1
                        IF(.not.tRegenDiagHEls) CurrentH(j)=0.D0
                        IF((.not.TNoAnnihil).and.(.not.tRotoAnnihil)) THEN
    !                    IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                            HashArray(j)=HFHash
                        ENDIF
                    enddo
                ENDIF
            ENDIF

            WRITE(6,"(A,F14.6,A)") "Initial memory (without excitgens + temp arrays) consists of : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
            IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                WRITE(6,*) "Only one array of memory to store main particle list allocated..."
            ENDIF
            WRITE(6,*) "Initial memory allocation sucessful..."
            CALL FLUSH(6)

!Put a barrier here so all processes synchronise
            CALL MPI_Barrier(MPI_COMM_WORLD,error)
            IF(.not.TRegenExcitgens) THEN

                ALLOCATE(WalkVecExcits(MaxWalkersPart),stat=ierr)
                ALLOCATE(WalkVec2Excits(MaxWalkersPart),stat=ierr)
                ALLOCATE(ExcitGens(MaxWalkersPart),stat=ierr)
                ALLOCATE(FreeIndArray(MaxWalkersPart),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("InitFCIMMCCalcPar","Error in allocating walker excitation generators")

                do j=1,MaxWalkersPart
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
                    IF(.not.tRotoAnnihil) THEN

                        do j=2,InitWalkers
!Copy the HF excitation generator accross to each initial particle
                            CALL CopyExitGenPar(CurrentExcits(1),CurrentExcits(j),.false.)
                        enddo
                    ENDIF
                ENDIF
                MemoryAlloc=((HFExcit%nExcitMemLen)+8)*4*MaxWalkersPart    !This is the memory needed by all the excitation generator arrays

                WRITE(6,"(A,I12)") "Size of HF excitgen is: ",HFExcit%nExcitMemLen
                WRITE(6,"(A,F14.6,A)") "Probable maximum memory for excitgens is : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
                WRITE(6,*) "Initial allocation of excitation generators successful..."
            ELSE
                WRITE(6,"(A)") "Excitation generators will not be stored, but regenerated each time they are needed..."
            ENDIF
            IF(tRotoAnnihil) THEN
                WRITE(6,"(A,F14.6,A)") "Temp Arrays for annihilation cannot be more than : ",REAL(MaxSpawned*9*4,r2)/1048576.D0," Mb/Processor"
            ELSE
                IF(.not.tDirectAnnihil) THEN
                    WRITE(6,"(A,F14.6,A)") "Temp Arrays for annihilation cannot be more than : ",REAL(MaxWalkersPart*12,r2)/1048576.D0," Mb/Processor"
                ENDIF
            ENDIF
            CALL FLUSH(6)
        
            IF(TStartSinglePart) THEN
                IF(tDirectAnnihil) THEN
                    Proc=DetermineDetProc(iLutHF)
                    IF(iProcIndex.eq.Proc) THEN
                        TotWalkers=1
                        TotWalkersOld=1
                        TotParts=1
                        TotPartsOld=1
                    ELSE
                        TotWalkers=0
                        TotWalkersOld=0
                        TotParts=0
                        TotPartsOld=0
                    ENDIF
!Initialise global variables for calculation on the root node
                    IF(iProcIndex.eq.root) THEN
                        AllTotWalkers=1.D0
                        AllTotWalkersOld=1.D0
                        AllTotParts=1.D0
                        AllTotPartsOld=1.D0
                    ENDIF
                ELSE
                    TotWalkers=1
                    TotWalkersOld=1
                    TotParts=1
                    TotPartsOld=1
!Initialise global variables for calculation on the root node
                    IF(iProcIndex.eq.root) THEN
                        AllTotWalkers=REAL(nProcessors,r2)
                        AllTotWalkersOld=REAL(nProcessors,r2)
                        AllTotParts=REAL(nProcessors,r2)
                        AllTotPartsOld=REAL(nProcessors,r2)
                    ENDIF
                ENDIF
            ELSE
                IF(tDirectAnnihil) THEN
!In this, only one processor has initial particles.
                    Proc=DetermineDetProc(iLutHF)
                    IF(iProcIndex.eq.Proc) THEN
                        TotWalkers=1
                        TotWalkersOld=1
                        TotParts=InitWalkers
                        TotPartsOld=InitWalkers
                    ELSE
                        TotWalkers=0
                        TotWalkersOld=0
                        TotParts=0
                        TotPartsOld=0
                    ENDIF
                    IF(iProcIndex.eq.Root) THEN
                        AllTotWalkers=1.D0
                        AllTotWalkersOld=1.D0
                        AllTotParts=REAL(InitWalkers,r2)
                        AllTotPartsOld=REAL(InitWalkers,r2)
                    ENDIF
                ELSEIF(tRotoAnnihil) THEN
!This is different, since the arrays are stored in a compressed form, but on all processors.
                    TotWalkers=1
                    TotWalkersOld=1
                    TotParts=InitWalkers
                    TotPartsOld=InitWalkers
                    IF(iProcIndex.eq.root) THEN
                        AllTotWalkers=REAL(nProcessors,r2)
                        AllTotWalkersOld=REAL(nProcessors,r2)
                        AllTotParts=REAL(InitWalkers*nProcessors,r2)
                        AllTotPartsOld=REAL(InitWalkers*nProcessors,r2)
                    ENDIF
                ELSE
!TotWalkers contains the number of current walkers at each step
                    TotWalkers=InitWalkers
                    TotWalkersOld=InitWalkers
                    TotParts=TotWalkers
                    TotPartsOld=TotWalkersOld
!Initialise global variables for calculation on the root node
                    IF(iProcIndex.eq.root) THEN
                        AllTotWalkers=REAL(InitWalkers*nProcessors,r2)
                        AllTotWalkersOld=REAL(InitWalkers*nProcessors,r2)
                        AllTotParts=AllTotWalkers
                        AllTotPartsOld=AllTotWalkersOld
                    ENDIF
                ENDIF

            ENDIF

        ENDIF   !End if initial walkers method

        IF(tFindGuide) THEN
            WRITE(6,*) 'Finding the guiding function from approximately ',iGuideDets,' most populated determinants'
            IF((.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) CALL Stop_All("InitFCIMCCalcPar","Cannot use or find the guiding function without using ROTOANNIHILATION")
        ELSEIF(tUseGuide) THEN
            WRITE(6,*) 'Reading in the guiding function and scaling the number of walkers to ',iInitGuideParts
            IF(.not.tRotoAnnihil) CALL Stop_All("InitFCIMCCalcPar","Cannot use or find the guiding function without using ROTOANNIHILATION")
            IF(.not.tStartSinglePart) CALL Stop_All("InitFCIMCCalcPar","Must use single particle start when reading in the guiding function")

            ! Check the guiding function file exists, if so, set up the guiding function arrays based on this file.
            INQUIRE(FILE='GUIDINGFUNC',EXIST=exists)
            IF(exists) THEN
                CALL InitGuidingFunction()
            ELSE
                CALL Stop_All("InitFCIMCCalcPar","The guiding function file (GUIDINGFUNC) does not exist")
            ENDIF
        ENDIF

        IF(tSpawnDominant) THEN
            WRITE(6,*) 'Reading in from DOMINANTDETS, and only allowing spawning on the determinants in this file.'
            INQUIRE(FILE='DOMINANTDETS',EXIST=exists)
            IF(exists) THEN
                CALL InitSpawnDominant()
            ELSE
                CALL Stop_All("InitFCIMCCalcPar","The dominant determinant file (DOMINANTDETS) does not exist")
            ENDIF
        ELSEIF(tMinorDetsStar.and.(.not.tSpawnDominant)) THEN
            CALL Stop_All("InitFCIMCCalcPar","Cannot use the star approximation on the insignificant determinants if the SPAWNDOMINANTONLY option is not on.")
        ENDIF

        IF(tPrintTriConnections.or.tHistTriConHEls.or.tPrintHElAccept) CALL InitTriHElStats()
        IF(tPrintSpinCoupHEl) CALL InitSpinCoupHEl()

        IF((NMCyc.ne.0).and.(tRotateOrbs.and.(.not.tFindCINatOrbs))) CALL Stop_All(this_routine,"Currently not set up to rotate and then go straight into a spawning &
                                                                                    & calculation.  Ordering of orbitals is incorrect.  This may be fixed if needed.")

!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
        RETURN

    END SUBROUTINE InitFCIMCCalcPar


    SUBROUTINE WriteFciMCStatsHeader()

        IF(iProcIndex.eq.root) THEN
!Print out initial starting configurations
            WRITE(6,*) ""
            IF(TLocalAnnihilation) THEN
                WRITE(6,"(A)") "       Step     Shift      WalkerCng    GrowRate       TotWalkers    LocalAnn   TotAnnihil    NoDied    NoBorn    Proj.E          Proj.E.Iter     NoatHF NoatDoubs      AvSign    AvSignHF+D   AccRat       MeanEx     MinEx MaxEx"
                WRITE(15,"(A)") "#       Step     Shift      WalkerCng    GrowRate       TotWalkers   LocalAnn    TotAnnihil    NoDied    NoBorn    Proj.E          Proj.E.Iter     NoatHF NoatDoubs       AvSign    AvSignHF+D   AccRat       MeanEx     MinEx MaxEx"

            ELSEIF(tUseGuide) THEN
                WRITE(6,"(A12,A16,A10,A16,A12,3A11,3A17,2A10,A13,A12,A13)") "Step","Shift","WalkerCng","GrowRate","TotWalkers","Annihil","NoDied","NoBorn","Proj.E","Proj.E.Iter",&
&               "Proj.E.ThisCyc","NoatHF","NoatDoubs","AccRat","UniqueDets","IterTime"

                WRITE(15,"(A12,A16,A10,A16,A12,3A11,3A17,2A10,A13,A12,A13,A18,A14)") "#","Step","Shift","WalkerCng","GrowRate","TotWalkers","Annihil","NoDied","NoBorn","Proj.E","Proj.E.Iter",&
&               "Proj.E.ThisCyc","NoatHF","NoatDoubs","AccRat","UniqueDets","IterTime","FracSpawnFromSing","NoinGuideFunc"

            ELSEIF(tMinorDetsStar) THEN
                WRITE(6,"(A12,A16,A10,A16,A12,3A11,3A17,2A10,A13,A12,A13)") "Step","Shift","WalkerCng","GrowRate","TotWalkers","Annihil","NoDied","NoBorn","Proj.E","Proj.E.Iter",&
&               "Proj.E.ThisCyc","NoatHF","NoatDoubs","AccRat","UniqueDets","IterTime"

                WRITE(15,"(A12,A16,A10,A16,A12,3A11,3A17,2A10,A13,A12,A13,A18,A16)") "#","Step","Shift","WalkerCng","GrowRate","TotWalkers","Annihil","NoDied","NoBorn","Proj.E","Proj.E.Iter",&
&               "Proj.E.ThisCyc","NoatHF","NoatDoubs","AccRat","UniqueDets","IterTime","FracSpawnFromSing","NoMinorWalkers"


            ELSE
                WRITE(6,"(A)") "       Step     Shift      WalkerCng    GrowRate       TotWalkers    Annihil    NoDied    NoBorn    Proj.E          Av.Shift     Proj.E.ThisCyc   NoatHF NoatDoubs      AccRat     UniqueDets     IterTime"
                WRITE(15,"(A)") "#       Step     Shift      WalkerCng    GrowRate       TotWalkers    Annihil    NoDied    NoBorn    Proj.E          Av.Shift",&
&               "Proj.E.ThisCyc   NoatHF NoatDoubs       AccRat     UniqueDets     IterTime    FracSpawnFromSing    WalkersDiffProc"
            
            ENDIF
        ENDIF

    END SUBROUTINE WriteFciMCStatsHeader
        

    
! This is called if tListDets is set, and will read a list of NAllowedDetList determinants in natural order from the SpawnOnlyDets file, 
! store them in AllowedDetList (compressed), and only allow spawning at these determinants.
    SUBROUTINE ReadSpawnListDets()
        use Parallel, only : MPIIBCast_Scal,MPIIBCast
        LOGICAL :: exists
        INTEGER :: i,nI(NEl),ierr

        IF(iProcIndex.eq.Root) THEN

            INQUIRE(FILE='SpawnOnlyDets',EXIST=exists)
            IF(.not.exists) THEN
                CALL Stop_All('ReadSpawnListDets',"No SpawnOnlyDets file to read in allowed determinants...")
            ENDIF

            WRITE(6,*) "Reading in the allowed determinant list from SpawnOnlyDets..."
            OPEN(17,FILE='SpawnOnlyDets',Status='old')

            NAllowedDetList=0
                
            do while(.true.)

                READ(17,*,END=199) nI(1:NEl)
!                WRITE(6,*) nI

                NAllowedDetList=NAllowedDetList+1
                    
            enddo

199         CONTINUE
        ENDIF

        CALL MPIIBCast_Scal(NAllowedDetList,Root)
        WRITE(6,*) NAllowedDetList, " determinants read in from SpawnOnlyDets file..."

        ALLOCATE(AllowedDetList(0:NIfD,NAllowedDetList),stat=ierr)
        IF(ierr.ne.0) THEN
            CALL Stop_All("ReadSpawnListDets","Error allocating AllowedDetList array")
        ENDIF

        IF(iProcIndex.eq.Root) THEN
            REWIND(17)

            do i=1,NAllowedDetList
                READ(17,*) nI(1:NEl)
                CALL EncodeBitDet(nI,AllowedDetList(0:NIfD,i),NEl,NIfD)
!                WRITE(6,*) AllowedDetList(0:NIfD,i)
            enddo

            CLOSE(17)
        ENDIF

        CALL MPIIBCast(AllowedDetList,NAllowedDetList*(NIfD+1),Root)
!        do i=1,NAllowedDetList
!            CALL DecodeBitDet(nI,AllowedDetList(0:NIfD,i),NEl,NIfD)
!            WRITE(6,*) nI(:)
!        enddo

    END SUBROUTINE ReadSpawnListDets


!This routine is the same as WriteToPopsfilePar, but does not require two main arrays to hold the data.
!The root processors data will be stored in a temporary array while it recieves the data from the other processors.
!This routine will write out to a popsfile. It transfers all walkers to the head node sequentially, so does not want to be called too often
    SUBROUTINE WriteToPopsfileParOneArr()
        REAL*8 :: TempSumNoatHF
        INTEGER :: error,WalkersonNodes(0:nProcessors-1)
        INTEGER :: Stat(MPI_STATUS_SIZE),Tag,Total,i,j,k
        INTEGER , ALLOCATABLE :: OrigSign(:), OrigParts(:,:)
        INTEGER :: OrigSignTag=0,OrigPartsTag=0
        CHARACTER(len=*) , PARAMETER :: this_routine='WriteToPopsfileParOneArr'

        CALL MPI_Barrier(MPI_COMM_WORLD,error)  !sync

!First, make sure we have up-to-date information - again collect AllTotWalkers,AllSumNoatHF and AllSumENum...
!        CALL MPI_Reduce(TotWalkers,AllTotWalkers,1,MPI_INTEGER,MPI_Sum,root,MPI_COMM_WORLD,error)    
!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        TempSumNoatHF=real(SumNoatHF,r2)
        CALL MPIDSumRoot(TempSumNoatHF,1,AllSumNoatHF,Root)
        CALL MPIDSumRoot(SumENum,1,AllSumENum,Root)

!We also need to tell the root processor how many particles to expect from each node - these are gathered into WalkersonNodes
        CALL MPI_AllGather(TotWalkers,1,MPI_INTEGER,WalkersonNodes,1,MPI_INTEGER,MPI_COMM_WORLD,error)
        do i=0,nProcessors-1
            IF(INT(WalkersonNodes(i)/iPopsPartEvery).lt.1) THEN
                RETURN
            ENDIF
        enddo

        Tag=125
!        WRITE(6,*) "Get Here"
!        CALL FLUSH(6)

        IF(iProcIndex.eq.root) THEN
!First, check that we are going to receive the correct number of particles...
            Total=0
            do i=0,nProcessors-1
                Total=Total+INT(WalkersonNodes(i)/iPopsPartEvery)
            enddo
            AllTotWalkers=REAL(Total,r2)
!            IF(Total.ne.AllTotWalkers) THEN
!                CALL Stop_All("WriteToPopsfilePar","Not all walkers accounted for...")
!            ENDIF

!Write header information
            IF(iPopsPartEvery.ne.1) THEN
                IF(tBinPops) THEN
                    WRITE(6,"(A,I12,A)") "Writing a binary reduced POPSFILEBIN, printing a total of ",INT(AllTotWalkers,i2), " particles."
                ELSE
                    WRITE(6,"(A,I12,A)") "Writing a reduced POPSFILE, printing a total of ",INT(AllTotWalkers,i2), " particles."
                ENDIF
            ELSE
                IF(tBinPops) THEN
                    WRITE(6,*) "Writing to binary POPSFILEBIN..."
                ELSE
                    WRITE(6,*) "Writing to POPSFILE..."
                ENDIF
            ENDIF
            IF(tBinPops) THEN
                OPEN(17,FILE='POPSFILEHEAD',Status='replace')
            ELSE
                OPEN(17,FILE='POPSFILE',Status='replace')
            ENDIF
            WRITE(17,*) AllTotWalkers,"   TOTWALKERS (all nodes)"
            WRITE(17,*) DiagSft,"   DIAG SHIFT"
            WRITE(17,*) NINT(AllSumNoatHF,i2),"   SUMNOATHF (all nodes)"
            WRITE(17,*) AllSumENum,"   SUMENUM ( \sum<D0|H|Psi> - all nodes)"
            WRITE(17,*) Iter+PreviousCycles,"   PREVIOUS CYCLES"
            IF(tBinPops) THEN
                CLOSE(17)
                OPEN(17,FILE='POPSFILEBIN',Status='replace',form='unformatted')
            ENDIF

            IF(tBinPops) THEN
                do j=1,TotWalkers
!First write out walkers on head node
                    IF(mod(j,iPopsPartEvery).eq.0) THEN
                        WRITE(17) CurrentDets(0:NIfD,j),CurrentSign(j)
                    ENDIF
                enddo
            ELSE
                do j=1,TotWalkers
!First write out walkers on head node
                    IF(mod(j,iPopsPartEvery).eq.0) THEN
                        do k=0,NIfD
                            WRITE(17,"(I20)",advance='no') CurrentDets(k,j)
                        enddo
                        WRITE(17,*) CurrentSign(j)
                    ENDIF
                enddo
            ENDIF
!            WRITE(6,*) "Written out own walkers..."
!            CALL FLUSH(6)

!Now, we copy the head nodes data to a new array...
            ALLOCATE(OrigSign(TotWalkers),stat=error)
            CALL LogMemAlloc('OrigSign',TotWalkers,4,this_routine,OrigSignTag,error)
            ALLOCATE(OrigParts(0:NIfD,TotWalkers),stat=error)
            CALL LogMemAlloc('OrigParts',TotWalkers*(NIfD+1),4,this_routine,OrigPartsTag,error)
            do i=1,TotWalkers
                OrigSign(i)=CurrentSign(i)
                OrigParts(:,i)=CurrentDets(:,i)
            enddo

!Now we need to receive the data from each other processor sequentially
!We can overwrite the head nodes information, since we have now stored it elsewhere.
            do i=1,nProcessors-1
!Run through all other processors...receive the data...
                CALL MPI_Recv(CurrentDets(0:NIfD,1:WalkersonNodes(i)),WalkersonNodes(i)*(NIfD+1),MPI_INTEGER,i,Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(CurrentSign(1:WalkersonNodes(i)),WalkersonNodes(i),MPI_INTEGER,i,Tag,MPI_COMM_WORLD,Stat,error)
!                WRITE(6,*) "Recieved walkers for processor ",i
!                CALL FLUSH(6)
                
!Then write it out...
                IF(tBinPops) THEN
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            WRITE(17) CurrentDets(0:NIfD,j),CurrentSign(j)
                        ENDIF
                    enddo
                ELSE
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            do k=0,NIfD
                                WRITE(17,"(I20)",advance='no') CurrentDets(k,j)
                            enddo
                            WRITE(17,*) CurrentSign(j)
                        ENDIF
                    enddo
                ENDIF
!                WRITE(6,*) "Writted out walkers for processor ",i
!                CALL FLUSH(6)

            enddo

            CLOSE(17)

!Now we need to copy the head processors original information back to itself again.
            do i=1,TotWalkers
                CurrentSign(i)=OrigSign(i)
                CurrentDets(:,i)=OrigParts(:,i)
            enddo
!Deallocate memory for temporary storage of information.
            DEALLOCATE(OrigSign)
            DEALLOCATE(OrigParts)
            CALL LogMemDealloc(this_routine,OrigSignTag)
            CALL LogMemDealloc(this_routine,OrigPartsTag)

        ELSE
!All other processors need to send their data to root...
            CALL MPI_Send(CurrentDets(0:NIfD,1:TotWalkers),TotWalkers*(NIfD+1),MPI_INTEGER,root,Tag,MPI_COMM_WORLD,error)
            CALL MPI_Send(CurrentSign(1:TotWalkers),TotWalkers,MPI_INTEGER,root,Tag,MPI_COMM_WORLD,error)
!            WRITE(6,*) "Have sent info to head node..."
!            CALL FLUSH(6)
        ENDIF

!Reset the values of the global variables
        AllSumNoatHF=0.D0
        AllSumENum=0.D0
        AllTotWalkers=0.D0

        RETURN

    END SUBROUTINE WriteToPopsfileParOneArr

!This routine will write out to a popsfile. It transfers all walkers to the head node sequentially, so does not want to be called too often
!When arriving at this routine, CurrentXXX are arrays with the data, and NewXXX will be used by the root processor to temporarily store the information
    SUBROUTINE WriteToPopsfilePar()
        REAL*8 :: TempSumNoatHF
        INTEGER :: error,WalkersonNodes(0:nProcessors-1)
        INTEGER :: Stat(MPI_STATUS_SIZE),Tag,Total,i,j,k

        CALL MPI_Barrier(MPI_COMM_WORLD,error)  !sync

!First, make sure we have up-to-date information - again collect AllTotWalkers,AllSumNoatHF and AllSumENum...
!        CALL MPI_Reduce(TotWalkers,AllTotWalkers,1,MPI_INTEGER,MPI_Sum,root,MPI_COMM_WORLD,error)    
!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        TempSumNoatHF=real(SumNoatHF,r2)
        CALL MPIDSumRoot(TempSumNoatHF,1,AllSumNoatHF,Root)
        CALL MPIDSumRoot(SumENum,1,AllSumENum,Root)

!We also need to tell the root processor how many particles to expect from each node - these are gathered into WalkersonNodes
        CALL MPI_AllGather(TotWalkers,1,MPI_INTEGER,WalkersonNodes,1,MPI_INTEGER,MPI_COMM_WORLD,error)
        do i=0,nProcessors-1
            IF(INT(WalkersonNodes(i)/iPopsPartEvery).lt.1) THEN
                RETURN
            ENDIF
        enddo

        Tag=125
!        WRITE(6,*) "Get Here"
!        CALL FLUSH(6)

        IF(iProcIndex.eq.root) THEN
!First, check that we are going to receive the correct number of particles...
            Total=0
            do i=0,nProcessors-1
                Total=Total+INT(WalkersonNodes(i)/iPopsPartEvery)
            enddo
            AllTotWalkers=REAL(Total,r2)
!            IF(Total.ne.AllTotWalkers) THEN
!                CALL Stop_All("WriteToPopsfilePar","Not all walkers accounted for...")
!            ENDIF

!Write header information
            IF(iPopsPartEvery.ne.1) THEN
                IF(tBinPops) THEN
                    WRITE(6,"(A,I12,A)") "Writing a binary reduced POPSFILEBIN, printing a total of ",INT(AllTotWalkers,i2), " particles."
                ELSE
                    WRITE(6,"(A,I12,A)") "Writing a reduced POPSFILE, printing a total of ",INT(AllTotWalkers,i2), " particles."
                ENDIF
            ELSE
                IF(tBinPops) THEN
                    WRITE(6,*) "Writing to binary POPSFILEBIN..."
                ELSE
                    WRITE(6,*) "Writing to POPSFILE..."
                ENDIF
            ENDIF
            IF(tBinPops) THEN
                OPEN(17,FILE='POPSFILEHEAD',Status='replace')
            ELSE
                OPEN(17,FILE='POPSFILE',Status='replace')
            ENDIF
            WRITE(17,*) AllTotWalkers,"   TOTWALKERS (all nodes)"
            WRITE(17,*) DiagSft,"   DIAG SHIFT"
            WRITE(17,*) NINT(AllSumNoatHF,i2),"   SUMNOATHF (all nodes)"
            WRITE(17,*) AllSumENum,"   SUMENUM ( \sum<D0|H|Psi> - all nodes)"
            WRITE(17,*) Iter+PreviousCycles,"   PREVIOUS CYCLES"
            IF(tBinPops) THEN
                CLOSE(17)
                OPEN(17,FILE='POPSFILEBIN',Status='replace',form='unformatted')
            ENDIF

            IF(tBinPops) THEN
                do j=1,TotWalkers
!First write out walkers on head node
                    IF(mod(j,iPopsPartEvery).eq.0) THEN
                        WRITE(17) CurrentDets(0:NIfD,j),CurrentSign(j)
                    ENDIF
                enddo
            ELSE
                do j=1,TotWalkers
!First write out walkers on head node
                    IF(mod(j,iPopsPartEvery).eq.0) THEN
                        do k=0,NIfD
                            WRITE(17,"(I20)",advance='no') CurrentDets(k,j)
                        enddo
                        WRITE(17,*) CurrentSign(j)
                    ENDIF
                enddo
            ENDIF
!            WRITE(6,*) "Written out own walkers..."
!            CALL FLUSH(6)

!Now we need to receive the data from each other processor sequentially
            do i=1,nProcessors-1
!Run through all other processors...receive the data...
                CALL MPI_Recv(NewDets(0:NIfD,1:WalkersonNodes(i)),WalkersonNodes(i)*(NIfD+1),MPI_INTEGER,i,Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(NewSign(1:WalkersonNodes(i)),WalkersonNodes(i),MPI_INTEGER,i,Tag,MPI_COMM_WORLD,Stat,error)
!                WRITE(6,*) "Recieved walkers for processor ",i
!                CALL FLUSH(6)
                
!Then write it out...
                IF(tBinPops) THEN
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            WRITE(17) NewDets(0:NIfD,j),NewSign(j)
                        ENDIF
                    enddo
                ELSE
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            do k=0,NIfD
                                WRITE(17,"(I20)",advance='no') NewDets(k,j)
                            enddo
                            WRITE(17,*) NewSign(j)
                        ENDIF
                    enddo
                ENDIF
!                WRITE(6,*) "Writted out walkers for processor ",i
!                CALL FLUSH(6)

            enddo

            CLOSE(17)

        ELSE
!All other processors need to send their data to root...
            CALL MPI_Send(CurrentDets(0:NIfD,1:TotWalkers),TotWalkers*(NIfD+1),MPI_INTEGER,root,Tag,MPI_COMM_WORLD,error)
            CALL MPI_Send(CurrentSign(1:TotWalkers),TotWalkers,MPI_INTEGER,root,Tag,MPI_COMM_WORLD,error)
!            WRITE(6,*) "Have sent info to head node..."
!            CALL FLUSH(6)
        ENDIF

!Reset the values of the global variables
        AllSumNoatHF=0.D0
        AllSumENum=0.D0
        AllTotWalkers=0.D0

        RETURN

    END SUBROUTINE WriteToPopsfilePar


!This routine reads in particle configurations from a POPSFILE.
    SUBROUTINE ReadFromPopsfilePar()
        LOGICAL :: exists,First,tBinRead,DetBitEQ
        INTEGER :: AvWalkers,WalkerstoReceive(nProcessors)
        INTEGER*8 :: NodeSumNoatHF(nProcessors),TempAllSumNoatHF
        REAL*8 :: TempTotParts
        INTEGER :: TempInitWalkers,error,i,j,k,l,total,ierr,MemoryAlloc,Tag
        INTEGER :: Stat(MPI_STATUS_SIZE),AvSumNoatHF,VecSlot,IntegerPart,HFPointer,TempnI(NEl),ExcitLevel,VecInd,DetsMerged
        REAL*8 :: r,FracPart,TempTotWalkers
        TYPE(HElement) :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        
        PreviousCycles=0    !Zero previous cycles
        SumENum=0.D0
        TotParts=0
        SumNoatHF=0
        DiagSft=0.D0
        Tag=124             !Set Tag

        IF(iProcIndex.eq.root) THEN
            INQUIRE(FILE='POPSFILE',EXIST=exists)
            IF(exists) THEN
                OPEN(17,FILE='POPSFILE',Status='old')
                tBinRead=.false.
            ELSE
!Reading in a binary file
                tBinRead=.true.
                INQUIRE(FILE='POPSFILEHEAD',EXIST=exists)
                IF(.not.exists) THEN
                    INQUIRE(FILE='POPSFILEBIN',EXIST=exists)
                    IF(.not.exists) THEN
                        CALL Stop_All(this_routine,"No POPSFILE's of any kind found.")
                    ELSE
                        CALL Stop_All(this_routine,"POPSFILEBIN found, but POPSFILEHEAD also needed for header information")
                    ENDIF
                ELSE
                    INQUIRE(FILE='POPSFILEBIN',EXIST=exists)
                    IF(.not.exists) THEN
                        CALL Stop_All(this_routine,"POPSFILEHEAD found, but no POPSFILEBIN for particle information - this is also needed")
                    ELSE
                        OPEN(17,FILE='POPSFILEHEAD',Status='old')
                    ENDIF
                ENDIF
            ENDIF
                    
!Read in initial data on processors which have a popsfile
            READ(17,*) AllTotWalkers
            READ(17,*) DiagSft
            READ(17,*) TempAllSumNoatHF     !AllSumNoatHF stored as integer for compatability with serial POPSFILEs
            READ(17,*) AllSumENum
            READ(17,*) PreviousCycles

            IF(tBinRead) THEN
                CLOSE(17)
                OPEN(17,FILE='POPSFILEBIN',Status='old',form='unformatted')
            ENDIF

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
            AvWalkers=NINT(AllTotWalkers/real(nProcessors,r2))
!Divide up the walkers to receive for each node
            do i=1,nProcessors-1
                WalkerstoReceive(i)=AvWalkers
            enddo
!The last processor takes the 'remainder'
            WalkerstoReceive(nProcessors)=NINT(AllTotWalkers)-(AvWalkers*(nProcessors-1))

!Quick check to ensure we have all walkers accounted for
            total=0
            do i=1,nProcessors
                total=total+WalkerstoReceive(i)
            enddo
            IF(total.ne.NINT(AllTotWalkers)) THEN
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
        MaxWalkersPart=NINT(MemoryFacPart*InitWalkers)    !All nodes have the same amount of memory allocated
        IF(tRotoAnnihil) THEN
            MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
        ELSE
            MaxWalkersAnnihil=NINT(MemoryFacAnnihil*InitWalkers)
        ENDIF

!Allocate memory to hold walkers at least temporarily
        ALLOCATE(WalkVecDets(0:NIfD,MaxWalkersPart),stat=ierr)
        CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfD+1),4,this_routine,WalkVecDetsTag,ierr)
        WalkVecDets(0:NIfD,1:MaxWalkersPart)=0
        IF(tRotoAnnihil) THEN
            ALLOCATE(WalkVecSign(MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkersPart,4,this_routine,WalkVecSignTag,ierr)
            WalkVecSign(:)=0
        ELSE
            ALLOCATE(WalkVecSign(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkersAnnihil,4,this_routine,WalkVecSignTag,ierr)
        ENDIF

        IF(iProcIndex.eq.root) THEN
!Root process reads all walkers in and then sends them to the correct processor
            do i=nProcessors,1,-1
!Read in data for processor i
                IF(tBinRead) THEN
                    do j=1,WalkerstoReceive(i)
                        READ(17) WalkVecDets(0:NIfD,j),WalkVecSign(j)
                    enddo
                ELSE
                    do j=1,WalkerstoReceive(i)
                        READ(17,*) WalkVecDets(0:NIfD,j),WalkVecSign(j)
                    enddo
                ENDIF

                IF(i.ne.1) THEN
!Now send data to processor i-1 (Processor rank goes from 0 -> nProcs-1). If i=1, then we want the data so stay at the root processor
                    CALL MPI_Send(WalkVecDets(:,1:WalkerstoReceive(i)),WalkerstoReceive(i)*(NIfD+1),MPI_INTEGER,i-1,Tag,MPI_COMM_WORLD,error)
                    CALL MPI_Send(WalkVecSign(1:WalkerstoReceive(i)),WalkerstoReceive(i),MPI_INTEGER,i-1,Tag,MPI_COMM_WORLD,error)
                ENDIF

            enddo

            CLOSE(17)

            IF(tRotoAnnihil) THEN
                WRITE(6,*) "Ordering/compressing all walkers that have been read in..."
                CALL SortBitDets(WalkerstoReceive(1),WalkVecDets(0:NIfD,1:WalkerstoReceive(1)),NIfD,WalkVecSign(1:WalkerstoReceive(1)))
                
                VecInd=1
                DetsMerged=0
                do i=2,TempInitWalkers
                    IF(.not.DetBitEQ(WalkVecDets(0:NIfD,i),WalkVecDets(0:NIfD,VecInd),NIfD)) THEN
                        VecInd=VecInd+1
                        WalkVecDets(:,VecInd)=WalkVecDets(:,i)
                        WalkVecSign(VecInd)=WalkVecSign(i)
                    ELSE
                        WalkVecSign(VecInd)=WalkVecSign(VecInd)+WalkVecSign(i)
                        DetsMerged=DetsMerged+1
                    ENDIF
                enddo

                TempInitWalkers=TempInitWalkers-DetsMerged

            ENDIF

        
        ENDIF

        do i=1,nProcessors-1
            IF(iProcIndex.eq.i) THEN
!All other processors want to pick up their data from root
                CALL MPI_Recv(WalkVecDets(:,1:TempInitWalkers),TempInitWalkers*(NIfD+1),MPI_INTEGER,0,Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(WalkVecSign(1:TempInitWalkers),TempInitWalkers,MPI_INTEGER,0,Tag,MPI_COMM_WORLD,Stat,error)
                IF(tRotoAnnihil) THEN
                    CALL SortBitDets(TempInitWalkers,WalkVecDets(0:NIfD,1:TempInitWalkers),NIfD,WalkVecSign(1:TempInitWalkers))
                    VecInd=1
                    DetsMerged=0
                    do l=2,TempInitWalkers
                        IF(.not.DetBitEQ(WalkVecDets(0:NIfD,l),WalkVecDets(0:NIfD,VecInd),NIfD)) THEN
                            VecInd=VecInd+1
                            WalkVecDets(:,VecInd)=WalkVecDets(:,l)
                            WalkVecSign(VecInd)=WalkVecSign(l)
                        ELSE
                            WalkVecSign(VecInd)=WalkVecSign(VecInd)+WalkVecSign(l)
                            DetsMerged=DetsMerged+1
                        ENDIF
                    enddo

                    TempInitWalkers=TempInitWalkers-DetsMerged

                ENDIF
            ENDIF
        enddo

        IF(iProcIndex.eq.root) WRITE(6,*) INT(AllTotWalkers,i2)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers  by a factor of: ",ScaleWalkers
            MaxWalkersPart=NINT(MemoryFacPart*(NINT(InitWalkers*ScaleWalkers)))   !InitWalkers here is simply the average number of walkers per node, not actual
            IF(tRotoAnnihil) THEN
                MaxSpawned=NINT(MemoryFacSpawn*(NINT(InitWalkers*ScaleWalkers)))
            ELSE
                MaxWalkersAnnihil=NINT(MemoryFacAnnihil*(NINT(InitWalkers*ScaleWalkers)))   !InitWalkers here is simply the average number of walkers per node, not actual
            ENDIF

            IF(tRotoAnnihil) THEN
                IntegerPart=INT(ScaleWalkers)
                FracPart=ScaleWalkers-REAL(IntegerPart)

                do l=1,TempInitWalkers
                    WalkVecSign(l)=WalkVecSign(l)*IntegerPart
                    IF(tMerTwist) THEN
                        CALL genrand_real2(r) 
                    ELSE
                        CALL RANLUX(r,1)
                    ENDIF
                    IF(r.lt.FracPart) THEN
!Stochastically create another particle
                        IF(WalkVecSign(l).lt.0) THEN
                            WalkVecSign(l)=WalkVecSign(l)-1
                        ELSE
                            WalkVecSign(l)=WalkVecSign(l)+1
                        ENDIF
                    ENDIF
                enddo

                InitWalkers=NINT(InitWalkers*ScaleWalkers)  !New (average) number of initial particles for culling criteria
                TotWalkers=TempInitWalkers  !This is now number of determinants, rather than walkers
                TotWalkersOld=TempInitWalkers
!Collate the information about the new total number of walkers
                TempTotWalkers=REAL(TotWalkers,r2)
                CALL MPI_AllReduce(TempTotWalkers,AllTotWalkers,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
                IF(iProcIndex.eq.root) THEN
                    AllTotWalkersOld=AllTotWalkers
                    WRITE(6,*) "Total number of initial determinants occupied is now: ", INT(AllTotWalkers,i2)
                ENDIF

            ELSE

!Allocate memory for walkvec2, which will temporarily hold walkers
                ALLOCATE(WalkVec2Dets(0:NIfD,MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVec2Dets',MaxWalkersPart*(NIfD+1),4,this_routine,WalkVec2DetsTag,ierr)
                WalkVec2Dets(0:NIfD,1:MaxWalkersPart)=0
                ALLOCATE(WalkVec2Sign(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('WalkVec2Sign',MaxWalkersAnnihil,4,this_routine,WalkVec2SignTag,ierr)
                WalkVec2Sign(:)=0

!Scale up the integer part and fractional part seperately
                IntegerPart=INT(ScaleWalkers)       !Round to zero
                FracPart=ScaleWalkers-REAL(IntegerPart)

                VecSlot=1
                do l=1,TempInitWalkers
                    do k=1,IntegerPart
                        WalkVec2Dets(0:NIfD,VecSlot)=WalkVecDets(0:NIfD,l)
                        WalkVec2Sign(VecSlot)=WalkVecSign(l)
                        VecSlot=VecSlot+1
                    enddo
                    IF(tMerTwist) THEN
                        CALL genrand_real2(r) 
                    ELSE
                        CALL RANLUX(r,1)
                    ENDIF
                    IF(r.lt.FracPart) THEN
!Stochastically choose whether to create another particle
                        WalkVec2Dets(0:NIfD,VecSlot)=WalkVecDets(0:NIfD,l)
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
                TempTotWalkers=REAL(TotWalkers,r2)
                CALL MPI_AllReduce(TempTotWalkers,AllTotWalkers,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
                IF(iProcIndex.eq.root) THEN
                    AllTotWalkersOld=AllTotWalkers
                    WRITE(6,*) "Total number of initial walkers is now: ", INT(AllTotWalkers,i2)
                ENDIF


!Deallocate the arrays used to hold the original particles, and reallocate with correct size.
                DEALLOCATE(WalkVecDets)
                CALL LogMemDealloc(this_routine,WalkVecDetsTag)
                DEALLOCATE(WalkVecSign)
                CALL LogMemDealloc(this_routine,WalkVecSignTag)
                ALLOCATE(WalkVecDets(0:NIfD,MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfD+1),4,this_routine,WalkVecDetsTag,ierr)
                WalkVecDets(0:NIfD,1:MaxWalkersPart)=0
                IF(tRotoAnnihil) THEN
                    ALLOCATE(WalkVecSign(MaxWalkersPart),stat=ierr)
                    CALL LogMemAlloc('WalkVecSign',MaxWalkersPart,4,this_routine,WalkVecSignTag,ierr)
                ELSE
                    ALLOCATE(WalkVecSign(MaxWalkersAnnihil),stat=ierr)
                    CALL LogMemAlloc('WalkVecSign',MaxWalkersAnnihil,4,this_routine,WalkVecSignTag,ierr)
                ENDIF
                WalkVecSign(:)=0

!Transfer scaled particles back accross to WalkVecDets
                do l=1,TotWalkers
                    WalkVecDets(0:NIfD,l)=WalkVec2Dets(0:NIfD,l)
                    WalkVecSign(l)=WalkVec2Sign(l)
                enddo

!Zero the second array for good measure
                WalkVec2Dets(0:NIfD,1:MaxWalkersPart)=0

            ENDIF

        ELSE
!We are not scaling the number of walkers...

            IF(.not.tRotoAnnihil) THEN
                ALLOCATE(WalkVec2Dets(0:NIfD,MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVec2Dets',MaxWalkersPart*(NIfD+1),4,this_routine,WalkVec2DetsTag,ierr)
                WalkVec2Dets(0:NIfD,1:MaxWalkersPart)=0
                ALLOCATE(WalkVec2Sign(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('WalkVec2Sign',MaxWalkersAnnihil,4,this_routine,WalkVec2SignTag,ierr)
                WalkVec2Sign(:)=0
            ENDIF
                

            TotWalkers=TempInitWalkers      !Set the total number of walkers
            TotWalkersOld=TempInitWalkers
            IF(iProcIndex.eq.root) THEN
                AllTotWalkersOld=AllTotWalkers
                WRITE(6,*) "Total number of initial walkers is now: ",INT(AllTotWalkers,i2)
            ENDIF

        ENDIF
            
        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft

!Need to now allocate other arrays
        IF(.not.tRegenDiagHEls) THEN
            ALLOCATE(WalkVecH(MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecH',MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
            WalkVecH(:)=0.d0
            IF(.not.tRotoAnnihil) THEN
                ALLOCATE(WalkVec2H(MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVec2H',MaxWalkersPart,8,this_routine,WalkVec2HTag,ierr)
                WalkVec2H(:)=0.d0
            ENDIF
        ELSE
            IF(tRotoAnnihil) THEN
                WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*4*2,r2)/1048576.D0," Mb/Processor"
            ELSE
                WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*4*4,r2)/1048576.D0," Mb/Processor"
            ENDIF
        ENDIF

        IF(tRotoAnnihil) THEN
            MemoryAlloc=((NIfD+1+3)*MaxWalkersPart*4)
        ELSE
            MemoryAlloc=((2*MaxWalkersAnnihil)+(((2*(NIfD+1))+4)*MaxWalkersPart))*4    !Memory Allocated in bytes
        ENDIF
        IF(tRegenDiagHEls) THEN
            IF(tRotoAnnihil) THEN
                MemoryAlloc=MemoryAlloc-(MaxWalkersPart*4*2)
            ELSE
                MemoryAlloc=MemoryAlloc-(MaxWalkersPart*4*4)
            ENDIF
        ENDIF

        IF(tRotoAnnihil) THEN

            WRITE(6,"(A,I12,A)") "Spawning vectors allowing for a total of ",MaxSpawned," particles to be spawned in any one iteration."
            ALLOCATE(SpawnVec(0:NIfD,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfD+1),4,this_routine,SpawnVecTag,ierr)
            SpawnVec(:,:)=0
            ALLOCATE(SpawnVec2(0:NIfD,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfD+1),4,this_routine,SpawnVec2Tag,ierr)
            SpawnVec2(:,:)=0
            ALLOCATE(SpawnSignVec(0:MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnSignVec',MaxSpawned+1,4,this_routine,SpawnSignVecTag,ierr)
            SpawnSignVec(:)=0
            ALLOCATE(SpawnSignVec2(0:MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnSignVec2',MaxSpawned+1,4,this_routine,SpawnSignVec2Tag,ierr)
            SpawnSignVec2(:)=0

!Point at correct spawning arrays
            SpawnedParts=>SpawnVec
            SpawnedParts2=>SpawnVec2
            SpawnedSign=>SpawnSignVec
            SpawnedSign2=>SpawnSignVec2

            MemoryAlloc=MemoryAlloc+(((MaxSpawned+1)*2)+(2*MaxSpawned*(1+NIfD)))*4

!        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
        ELSEIF(.not.TNoAnnihil) THEN
            ALLOCATE(HashArray(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('HashArray',MaxWalkersAnnihil,8,this_routine,HashArrayTag,ierr)
            HashArray(:)=0
            ALLOCATE(Hash2Array(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('Hash2Array',MaxWalkersAnnihil,8,this_routine,Hash2ArrayTag,ierr)
            Hash2Array(:)=0
            ALLOCATE(IndexTable(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('IndexTable',MaxWalkersAnnihil,4,this_routine,IndexTableTag,ierr)
            IndexTable(1:MaxWalkersAnnihil)=0
            ALLOCATE(Index2Table(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('Index2Table',MaxWalkersAnnihil,4,this_routine,Index2TableTag,ierr)
            Index2Table(1:MaxWalkersAnnihil)=0
            ALLOCATE(ProcessVec(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('ProcessVec',MaxWalkersAnnihil,4,this_routine,ProcessVecTag,ierr)
            ProcessVec(1:MaxWalkersAnnihil)=0
            ALLOCATE(Process2Vec(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('Process2Vec',MaxWalkersAnnihil,4,this_routine,Process2VecTag,ierr)
            Process2Vec(1:MaxWalkersAnnihil)=0

            MemoryAlloc=MemoryAlloc+32*MaxWalkersAnnihil
        ENDIF

        WRITE(6,"(A,F14.6,A)") "Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        CALL FLUSH(6)

!Allocate pointers to the correct walker arrays...
        CurrentDets=>WalkVecDets
        CurrentSign=>WalkVecSign
!        CurrentIC=>WalkVecIC
        IF(.not.tRegenDiagHEls) THEN
            CurrentH=>WalkVecH
            NewH=>WalkVec2H
        ENDIF
        NewDets=>WalkVec2Dets
        NewSign=>WalkVec2Sign
!        NewIC=>WalkVec2IC

        IF(.not.TRegenExcitgens) THEN
            ALLOCATE(WalkVecExcits(MaxWalkersPart),stat=ierr)
            ALLOCATE(WalkVec2Excits(MaxWalkersPart),stat=ierr)
            ALLOCATE(ExcitGens(MaxWalkersPart),stat=ierr)
            ALLOCATE(FreeIndArray(MaxWalkersPart),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("ReadFromPopsfilePar","Error in allocating walker excitation generators")
            do j=1,MaxWalkersPart
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

            MemoryAlloc=((HFExcit%nExcitMemLen)+8)*4*MaxWalkersPart
            WRITE(6,"(A,F14.6,A)") "Probable maximum memory for excitgens is : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
            WRITE(6,*) "Initial allocation of excitation generators successful..."

        ELSE
            WRITE(6,*) "Excitgens will be regenerated when they are needed..."
        ENDIF

        IF(tRotoAnnihil) THEN
            WRITE(6,"(A,F14.6,A)") "Temp Arrays for annihilation cannot be more than : ",REAL(MaxSpawned*9*4,r2)/1048576.D0," Mb/Processor"
        ELSE
            WRITE(6,"(A,F14.6,A)") "Temp Arrays for annihilation cannot be more than : ",REAL(MaxWalkersPart*12,r2)/1048576.D0," Mb/Processor"
        ENDIF
        CALL FLUSH(6)

!Now find out the data needed for the particles which have been read in...
        First=.true.
        do j=1,TotWalkers
            CALL DecodeBitDet(TempnI,CurrentDets(:,j),NEl,NIfD)
            CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,Excitlevel,2)
            IF(Excitlevel.eq.0) THEN
                IF(.not.tRegenDiagHEls) CurrentH(j)=0.D0
                IF(First) THEN
!First run - create the excitation.
                    IF(.not.TRegenExcitgens) CALL SetupExitgenPar(HFDet,CurrentExcits(j))
                    HFPointer=j
                    First=.false.
                ELSE
                    IF(.not.TRegenExcitgens) CALL CopyExitGenPar(CurrentExcits(HFPointer),CurrentExcits(j),.false.)
                ENDIF
                IF((.not.TNoAnnihil).and.(.not.tRotoAnnihil)) THEN
!                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    HashArray(j)=HFHash
                ENDIF
            ELSE
                IF(.not.tRegenDiagHEls) THEN
                    IF(tHPHF) THEN
                        CALL HPHFGetDiagHElement(TempnI,CurrentDets(:,j),HElemTemp)
                    ELSE
                        HElemTemp=GetHElement2(TempnI,TempnI,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                    ENDIF
                    CurrentH(j)=REAL(HElemTemp%v,r2)-Hii
                ENDIF
                IF(.not.TRegenExcitgens) CurrentExcits(j)%PointToExcit=>null()
                IF((.not.TNoAnnihil).and.(.not.tRotoAnnihil)) THEN
!                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    HashArray(j)=CreateHash(TempnI)
                ENDIF
                
            ENDIF
            TotParts=TotParts+abs(CurrentSign(j))

        enddo
        TempTotParts=REAL(TotParts,r2)
        CALL MPI_AllReduce(TempTotParts,AllTotParts,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        AllTotPartsOld=AllTotParts

        RETURN

    END SUBROUTINE ReadFromPopsfilePar


!This will set up the initial walker distribution proportially to the MP1 wavevector.
    SUBROUTINE InitWalkersMP1Par()
        use SystemData , only : tAssumeSizeExcitgen
        INTEGER :: HFConn,error,ierr,MemoryAlloc,VecSlot,nJ(NEl),nStore(6),iExcit,i,j,WalkersonHF,HFPointer,ExcitLevel,VecInd
        REAL*8 :: SumMP1Compts,MP2Energy,Compt,r,FracPart,TempTotWalkers,TempTotParts
        TYPE(HElement) :: Hij,Hjj,Fjj
        INTEGER , ALLOCATABLE :: MP1Dets(:,:), ExcitgenTemp(:)
        INTEGER , ALLOCATABLE :: MP1Sign(:)
        REAL*8 , ALLOCATABLE :: MP1Comps(:),MP1CompsNonCum(:)
        INTEGER :: MP1DetsTag,MP1SignTag,MP1CompsTag,SumWalkersonHF,ExcitLength,iMaxExcit,IntParts,MP1CompsNonCumTag
        CHARACTER(len=*), PARAMETER :: this_routine='InitWalkersMP1Par'
        LOGICAL :: First,TurnBackAssumeExGen
        
    
        IF(tHub.and.tReal) THEN
            CALL Stop_All(this_routine,"Cannot initialise walkers in MP1 for real space hubbard calculations.")
        ENDIF
!Set the maximum number of walkers allowed
        MaxWalkersPart=NINT(MemoryFacPart*InitWalkers)
        WRITE(6,"(A,F14.5)") "Memory Factor for walkers is: ",MemoryFacPart
        WRITE(6,"(A,I14)") "Memory allocated for a maximum particle/det number per node of: ",MaxWalkersPart
        IF(tRotoAnnihil) THEN
            MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
            WRITE(6,"(A,F14.5)") "Memory Factor for arrays used for spawning is: ",MemoryFacSpawn
            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for spawning of: ",MaxSpawned
        ELSE
            MaxWalkersAnnihil=NINT(MemoryFacAnnihil*InitWalkers)
            WRITE(6,"(A,F14.5)") "Memory Factor for arrays used for annihilation is: ",MemoryFacAnnihil
            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for annihilation of: ",MaxWalkersAnnihil
        ENDIF

!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
!Allocate memory to hold walkers
        ALLOCATE(WalkVecDets(0:NIfD,MaxWalkersPart),stat=ierr)
        CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfD+1),4,this_routine,WalkVecDetsTag,ierr)
        WalkVecDets(0:NIfD,1:MaxWalkersPart)=0
        IF(tRotoAnnihil) THEN
            ALLOCATE(WalkVecSign(MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkersPart,4,this_routine,WalkVecSignTag,ierr)
        ELSE
            ALLOCATE(WalkVec2Dets(0:NIfD,MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVec2Dets',MaxWalkersPart*(NIfD+1),4,this_routine,WalkVec2DetsTag,ierr)
            WalkVec2Dets(0:NIfD,1:MaxWalkersPart)=0
            ALLOCATE(WalkVecSign(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkersAnnihil,4,this_routine,WalkVecSignTag,ierr)
            ALLOCATE(WalkVec2Sign(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('WalkVec2Sign',MaxWalkersAnnihil,4,this_routine,WalkVec2SignTag,ierr)
            WalkVec2Sign(:)=0
        ENDIF
        WalkVecSign(:)=0

!        ALLOCATE(WalkVecIC(MaxWalkersPart),stat=ierr)
!        CALL LogMemAlloc('WalkVecIC',MaxWalkersPart,4,this_routine,WalkVecICTag,ierr)
!        ALLOCATE(WalkVec2IC(MaxWalkersPart),stat=ierr)
!        CALL LogMemAlloc('WalkVec2IC',MaxWalkersPart,4,this_routine,WalkVec2ICTag,ierr)
        IF(.not.tRegenDiagHEls) THEN
            ALLOCATE(WalkVecH(MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecH',MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
            WalkVecH(:)=0.d0
            IF(.not.tRotoAnnihil) THEN
                ALLOCATE(WalkVec2H(MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVec2H',MaxWalkersPart,8,this_routine,WalkVec2HTag,ierr)
                WalkVec2H(:)=0.d0
            ENDIF
        ELSE
            IF(tRotoAnnihil) THEN
                WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*4*2,r2)/1048576.D0," Mb/Processor"
            ELSE
                WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*4*4,r2)/1048576.D0," Mb/Processor"
            ENDIF
        ENDIF
        
        IF(tRotoAnnihil) THEN
            MemoryAlloc=(NIfD+1+3)*MaxWalkersPart*4    !Memory Allocated in bytes
            IF(tRegenDiagHEls) MemoryAlloc=MemoryAlloc-(MaxWalkersPart*8)
        ELSE
            MemoryAlloc=((2*MaxWalkersAnnihil)+(((2*(NIfD+1))+4)*MaxWalkersPart))*4    !Memory Allocated in bytes
            IF(tRegenDiagHEls) MemoryAlloc=MemoryAlloc-(MaxWalkersPart*16)
        ENDIF


        IF(tRotoAnnihil) THEN
            
            WRITE(6,"(A,I12,A)") "Spawning vectors allowing for a total of ",MaxSpawned," particles to be spawned in any one iteration."
            ALLOCATE(SpawnVec(0:NIfD,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfD+1),4,this_routine,SpawnVecTag,ierr)
            SpawnVec(:,:)=0
            ALLOCATE(SpawnVec2(0:NIfD,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfD+1),4,this_routine,SpawnVec2Tag,ierr)
            SpawnVec2(:,:)=0
            ALLOCATE(SpawnSignVec(0:MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnSignVec',MaxSpawned+1,4,this_routine,SpawnSignVecTag,ierr)
            SpawnSignVec(:)=0
            ALLOCATE(SpawnSignVec2(0:MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnSignVec2',MaxSpawned+1,4,this_routine,SpawnSignVec2Tag,ierr)
            SpawnSignVec2(:)=0
            
            SpawnedParts=>SpawnVec
            SpawnedParts2=>SpawnVec2
            SpawnedSign=>SpawnSignVec
            SpawnedSign2=>SpawnSignVec2

            MemoryAlloc=MemoryAlloc+(((MaxSpawned+1)*2)+(2*MaxSpawned*(1+NIfD)))*4

        ELSEIF(.not.TNoAnnihil) THEN
!        IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
            ALLOCATE(HashArray(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('HashArray',MaxWalkersAnnihil,8,this_routine,HashArrayTag,ierr)
            HashArray(:)=0
            ALLOCATE(Hash2Array(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('Hash2Array',MaxWalkersAnnihil,8,this_routine,Hash2ArrayTag,ierr)
            Hash2Array(:)=0
            ALLOCATE(IndexTable(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('IndexTable',MaxWalkersAnnihil,4,this_routine,IndexTableTag,ierr)
            IndexTable(1:MaxWalkersAnnihil)=0
            ALLOCATE(Index2Table(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('Index2Table',MaxWalkersAnnihil,4,this_routine,Index2TableTag,ierr)
            Index2Table(1:MaxWalkersAnnihil)=0
            ALLOCATE(ProcessVec(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('ProcessVec',MaxWalkersAnnihil,4,this_routine,ProcessVecTag,ierr)
            ProcessVec(1:MaxWalkersAnnihil)=0
            ALLOCATE(Process2Vec(MaxWalkersAnnihil),stat=ierr)
            CALL LogMemAlloc('Process2Vec',MaxWalkersAnnihil,4,this_routine,Process2VecTag,ierr)
            Process2Vec(1:MaxWalkersAnnihil)=0

            MemoryAlloc=MemoryAlloc+32*MaxWalkersAnnihil
        ENDIF

!Allocate pointers to the correct walker arrays
        CurrentDets=>WalkVecDets
        CurrentSign=>WalkVecSign
        IF(.not.tRegenDiagHEls) THEN
            CurrentH=>WalkVecH
            NewH=>WalkVec2H
        ENDIF
        NewDets=>WalkVec2Dets
        NewSign=>WalkVec2Sign

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
        MP1Sign(:)=0
        ALLOCATE(MP1CompsNonCum(HFConn),stat=ierr)
        CALL LogMemAlloc('MP1CompsNonCum',HFConn,8,this_routine,MP1CompsNonCumTag,ierr)
        MP1CompsNonCum(:)=0.D0
        
!HF Compt. of MP1 is 1
        MP1Dets(1:NEl,1)=HFDet(1:NEl)
        MP1Comps(1)=1.D0
        MP1Sign(1)=1
        MP1CompsNonCum(1)=1.D0      !These are the components of the MP1 wavevector - not stored cumulatively as in MP1Comps

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
!            WRITE(6,"(8I5,2G25.10)") nJ(:),real(Hij%v,r2),(Fii-(REAL(Fjj%v,r2)))

            Compt=real(Hij%v,r2)/(Fii-(REAL(Fjj%v,r2)))
            IF(Compt.lt.0.D0) THEN
                MP1Sign(VecSlot)=-1
            ELSE
                MP1Sign(VecSlot)=1
            ENDIF
            MP1Dets(1:NEl,VecSlot)=nJ(:)
            MP1Comps(VecSlot)=MP1Comps(VecSlot-1)+abs(Compt)
            MP1CompsNonCum(VecSlot)=abs(Compt)
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
!        WRITE(6,*) VecSlot,HFConn,SumMP1Compts,ExcitLength

!        DiagSft=MP2Energy
        MP2Energy=MP2Energy+Hii
        WRITE(6,"(A,F16.7,A,F16.7)") "MP2 energy is ",MP2Energy,", but the initial shift has been set to: ",DiagSft

!        do i=1,VecSlot
!            WRITE(6,"(I5,3G20.10)") i,MP1Comps(i),MP1CompsNonCum(i),MP1CompsNonCum(i)/SumMP1Compts*REAL(InitWalkers,r2)
!        enddo

        WalkersonHF=0       !Find the number of walkers we are assigning to HF

        IF(tRotoAnnihil) THEN
!Here, we simply run through the MP1 components and then assign |amp|/sum{amps} x InitWalkers

            VecInd=1
            TotParts=0
            do j=1,VecSlot

                FracPart=MP1CompsNonCum(j)/SumMP1Compts*REAL(InitWalkers,r2)
                IntParts=INT(FracPart)
                FracPart=FracPart-REAL(IntParts)
!Determine whether we want to stochastically create another particle
                IF(tMerTwist) THEN
                    CALL genrand_real2(r) 
                ELSE
                    CALL RANLUX(r,1)
                ENDIF
                IF(FracPart.gt.r) THEN
                    IntParts=IntParts+1
                ENDIF

                IF(IntParts.gt.0) THEN

                    IF(j.eq.1) THEN
                        WalkersonHF=IntParts
                        CurrentDets(0:NIfD,VecInd)=iLutHF(:)
                        CurrentSign(VecInd)=IntParts*MP1Sign(j)
                        TotParts=TotParts+IntParts
                        IF(.not.tRegenDiagHEls) THEN
                            CurrentH(VecInd)=0.D0
                        ENDIF
                    ELSE
!We are at a double excitation - we need to calculate most of this information...
                        CALL EncodeBitDet(MP1Dets(1:NEl,j),CurrentDets(0:NIfD,VecInd),NEl,NIfD)
                        CurrentSign(VecInd)=IntParts*MP1Sign(j)
                        TotParts=TotParts+IntParts
                        IF(.not.tRegenDiagHEls) THEN
                            Hjj=GetHElement2(MP1Dets(1:NEl,j),MP1Dets(1:NEl,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)     !Find the diagonal element
                            CurrentH(VecInd)=real(Hjj%v,r2)-Hii
                        ENDIF
                    ENDIF
                    VecInd=VecInd+1
                ENDIF

            enddo

            TotWalkers=VecInd-1

!A reduced determinant representation could be created more easily by stochastically choosing amplitudes and running over excitations, rather than walkers.
            WRITE(6,*) "Ordering all walkers for rotoannihilation..."
            IF(tRegenDiagHEls) THEN
                CALL SortBitDets(TotWalkers,CurrentDets(0:NIfD,1:TotWalkers),NIfD,CurrentSign(1:TotWalkers))
            ELSE
                CALL SortBitDetswH(TotWalkers,CurrentDets(0:NIfD,1:TotWalkers),NIfD,CurrentSign(1:TotWalkers),CurrentH(1:TotWalkers))
            ENDIF
            
        ELSE
                        
            do j=1,InitWalkers

                i=1

                IF(tMerTwist) THEN
                    CALL genrand_real2(r) 
                ELSE
                    CALL RANLUX(r,1)
                ENDIF
                r=r*SumMP1Compts       !Choose the double that this walker wants to be put...
                do while(r.gt.MP1Comps(i))

                    i=i+1

                    IF(i.gt.VecSlot) THEN
                        CALL Stop_All("InitWalkersMP1","Assigning walkers stochastically has been performed incorrectly")
                    ENDIF

                enddo

                IF(i.eq.1) THEN
!If we are at HF, then we do not need to calculate the information for the walker        
                    WalkersonHF=WalkersonHF+1
                    CurrentDets(0:NIfD,j)=iLutHF(:)
                    CurrentSign(j)=1
                    IF(.not.tRegenDiagHEls) THEN
                        CurrentH(j)=0.D0
                    ENDIF
                    IF(.not.TNoAnnihil) THEN
!                    IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                        HashArray(j)=HFHash
                    ENDIF
                ELSE
!We are at a double excitation - we need to calculate most of this information...
                    CALL EncodeBitDet(MP1Dets(1:NEl,i),CurrentDets(0:NIfD,j),NEl,NIfD)
                    CurrentSign(j)=MP1Sign(i)
                    IF(.not.tRegenDiagHEls) THEN
                        Hjj=GetHElement2(MP1Dets(1:NEl,i),MP1Dets(1:NEl,i),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)     !Find the diagonal element
                        CurrentH(j)=real(Hjj%v,r2)-Hii
                    ENDIF
                    IF(.not.TNoAnnihil) THEN
!                    IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                        HashArray(j)=CreateHash(MP1Dets(1:NEl,i))
                    ENDIF
                    
                ENDIF

            enddo

        ENDIF
        
        CALL MPI_Reduce(WalkersonHF,SumWalkersonHF,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
        TempTotParts=REAL(TotParts,r2)
        CALL MPI_Reduce(TempTotParts,AllTotParts,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)

        IF(iProcIndex.eq.Root) THEN
            IF(tRotoAnnihil) THEN
                WRITE(6,"(A,I12,A,I12,A)") "Out of ",TotParts*nProcessors," initial walkers allocated, ",SumWalkersonHF," of them are situated on the HF determinant."
            ELSE
                WRITE(6,"(A,I12,A,I12,A)") "Out of ",InitWalkers*nProcessors," initial walkers allocated, ",SumWalkersonHF," of them are situated on the HF determinant."
                TotParts=InitWalkers
                AllTotParts=REAL(InitWalkers*nProcessors,r2)
            ENDIF
        ENDIF
        AllNoatHF=SumWalkersonHF
        AllNoatDoubs=INT(AllTotParts)-SumWalkersonHF

!Deallocate MP1 data
        DEALLOCATE(MP1Comps)
        CALL LogMemDealloc(this_routine,MP1CompsTag)
        DEALLOCATE(MP1Dets)
        CALL LogMemDealloc(this_routine,MP1DetsTag)
        DEALLOCATE(MP1Sign)
        CALL LogMemDealloc(this_routine,MP1SignTag)
        DEALLOCATE(MP1CompsNonCum)
        CALL LogMemDealloc(this_routine,MP1CompsNonCumTag)

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
            ALLOCATE(WalkVecExcits(MaxWalkersPart),stat=ierr)
            ALLOCATE(WalkVec2Excits(MaxWalkersPart),stat=ierr)
            ALLOCATE(ExcitGens(MaxWalkersPart),stat=ierr)
            ALLOCATE(FreeIndArray(MaxWalkersPart),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("InitFCIMMCCalcPar","Error in allocating walker excitation generators")

            do j=1,MaxWalkersPart
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
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,ExcitLevel,2)
                IF(ExcitLevel.eq.0) THEN
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
            MemoryAlloc=((HFExcit%nExcitMemLen)+8)*4*MaxWalkersPart

            WRITE(6,"(A,F14.6,A)") "Probable maximum memory for excitgens is : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb/Processor"
            WRITE(6,*) "Initial allocation of excitation generators successful..."
        ELSE
            WRITE(6,"(A)") "Excitation generators will not be stored, but regenerated each time they are needed..."
        ENDIF
        IF(tRotoAnnihil) THEN
            WRITE(6,"(A,F14.6,A)") "Temp Arrays for annihilation cannot be more than : ",REAL(MaxSpawned*9*4,r2)/1048576.D0," Mb/Processor"
        ELSE
            WRITE(6,"(A,F14.6,A)") "Temp Arrays for annihilation cannot be more than : ",REAL(MaxWalkersPart*12,r2)/1048576.D0," Mb/Processor"
        ENDIF
        CALL FLUSH(6)
        
        IF(tRotoAnnihil) THEN
            TotWalkersOld=TotWalkers
            TempTotWalkers=REAL(TotWalkers,r2)
            CALL MPI_Reduce(TempTotWalkers,AllTotWalkers,1,MPI_DOUBLE_PRECISION,MPI_SUM,root,MPI_COMM_WORLD,error)
            AllTotWalkersOld=AllTotWalkers
            AllTotPartsOld=AllTotParts
        ELSE
!TotWalkers contains the number of current walkers at each step
            TotWalkers=InitWalkers
            TotWalkersOld=InitWalkers
!Initialise global variables for calculation on the root node
            IF(iProcIndex.eq.root) THEN
                AllTotWalkers=REAL(InitWalkers*nProcessors,r2)
                AllTotWalkersOld=REAL(InitWalkers*nProcessors,r2)
            ENDIF
        ENDIF

        RETURN

    END SUBROUTINE InitWalkersMP1Par

!This routine flips the sign of all particles on the node
    SUBROUTINE FlipSign()
        INTEGER :: i

        do i=1,TotWalkers
            CurrentSign(i)=-CurrentSign(i)
        enddo
        
        IF(tMinorDetsStar) THEN
            do i=1,NoMinorWalkers
                MinorStarSign(i)=-MinorStarSign(i)
            enddo
        ENDIF

!Reverse the flag for whether the sign of the particles has been flipped so the ACF can be correctly calculated
        TFlippedSign=.not.TFlippedSign
        RETURN
    
    END SUBROUTINE FlipSign


!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
    INTEGER FUNCTION AttemptCreatePar(DetCurr,iLutCurr,WSign,nJ,iLutnJ,Prob,IC,Ex,tParity,nParts,tMinorDetList)
        use GenRandSymExcitNUMod , only : GenRandSymExcitBiased
        use Logging, only : CCMCDebug
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,StoreNumTo,StoreNumFrom,DetLT,i,ExtraCreate,Ex(2,2),WSign,nParts
        INTEGER :: iLutCurr(0:NIfD),Bin,iLutnJ(0:NIfD),PartInd,ExcitLev,iLut(0:NIfD),iLut2(0:NIfD)
        LOGICAL :: tParity,SymAllowed,tSuccess,tMinorDetList,DetBitEQ
        REAL*8 :: Prob,r,rat
        TYPE(HElement) :: rh,rhcheck

        IF(tImportanceSample) THEN
!Here, we generate an excitation and calculate how many to accept all in one go...
!Most of the arguments to this routine are passed out, rather than passed in.
            CALL GenRandSymExcitBiased(DetCurr,iLutCurr,nJ,pDoubles,IC,Ex,tParity,exFlag,nParts,WSign,Tau,AttemptCreatePar)
            RETURN
        ENDIF

        IF(tSpawnDominant) THEN
!We only allow spawning at determinants between iMinDomLev and iMaxDomLev which are in an allowed list of dominant determinants.
!Find the excitation level of the excitation
!We need to calculate the bit-representation of this new child. This can be done easily since the ExcitMat is known.
            tMinorDetList=.false.
            CALL FindExcitBitDet(iLutCurr,iLutnJ,IC,Ex,NIfD)
            CALL FindBitExcitLevel(iLutHF,iLutnJ,NIfD,ExcitLev,iMaxDomLev)
            IF((ExcitLev.le.iMaxDomLev).and.(ExcitLev.ge.iMinDomLev)) THEN
!We now need to binary search the list of allowed determinants to see whether it is in the list or not.
                CALL BinSearchParts3(iLutnJ,DomDets,iNoDomDets,DomExcIndex(ExcitLev),DomExcIndex(ExcitLev+1)-1,PartInd,tSuccess)
                IF(.not.tSuccess) THEN
!If the particle is not in the list, then do not allow a spawning event there, unless tMinorDetsStar is on.
                    IF(tMinorDetsStar) THEN
!Put into the MinorSpawnDets array, with its parent (iLutCurr)
                        tMinorDetList=.true.
                    ELSE
                        AttemptCreatePar=0
                        RETURN
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
                

!Calculate off diagonal hamiltonian matrix element between determinants
!        rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
        IF(tHPHF) THEN
            IF(tGenMatHEl) THEN
!The prob is actually prob/HEl, since the matrix element was generated at the same time as the excitation

                rat=Tau*REAL(nParts,r2)/abs(Prob)

                rh%v=Prob ! to get the signs right for later on.
!                WRITE(6,*) Prob, DetCurr(:),"***",nJ(:)
!                WRITE(6,*) "******"
!                CALL HPHFGetOffDiagHElement(DetCurr,nJ,iLutCurr,iLutnJ,rh)

            ELSE
!The IC given doesn't really matter. It just needs to know whether it is a diagonal or off-diagonal matrix element.
!However, the excitation generator can generate the same HPHF again. If this is done, the routine will send the matrix element back as zero.
                CALL HPHFGetOffDiagHElement(DetCurr,nJ,iLutCurr,iLutnJ,rh)
!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
                rat=Tau*abs(rh%v)*REAL(nParts,r2)/Prob
!                WRITE(6,*) Prob/rh%v, DetCurr(:),"***",nJ(:)
!                WRITE(6,*) "******"

            ENDIF
        ELSE
!Normal determinant spawn

            rh=GetHElement4(DetCurr,nJ,IC,Ex,tParity)
!            WRITE(6,*) rh%v

!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
            rat=Tau*abs(rh%v)*REAL(nParts,r2)/Prob
        ENDIF
        IF(CCMCDebug.gt.5) WRITE(6,*) "Connection H-element to spawnee:",rh
!        CALL IsSymAllowedExcit(DetCurr,nJ,IC,Ex,SymAllowed) 
!        IF((.not.SymAllowed).and.(abs(rh%v).gt.0.D0)) THEN
!            WRITE(17,*) rh%v
!        ENDIF

!        rhcheck=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!        IF(rh%v.ne.rhcheck%v) THEN
!            WRITE(6,*) "DetCurr: ",DetCurr(:)
!            WRITE(6,*) "nJ: ",nJ(:)
!            WRITE(6,*) "EX: ",Ex(1,:),Ex(2,:)
!            WRITE(6,*) "tParity: ",tParity
!            STOP
!        ENDIF

!        IF(abs(rh%v).le.HEpsilon) THEN
!            AttemptCreatePar=0
!            RETURN
!        ENDIF


!If probability is > 1, then we can just create multiple children at the chosen determinant
        ExtraCreate=INT(rat)
        rat=rat-REAL(ExtraCreate)


!Stochastically choose whether to create or not according to ranlux 
        IF(tMerTwist) THEN
            CALL genrand_real2(r) 
        ELSE
            CALL RANLUX(r,1)
        ENDIF
        IF(rat.gt.r) THEN
!            IF(Iter.eq.18925) THEN
!                WRITE(6,*) "Created",rh%v,rat
!            ENDIF

!Child is created - what sign is it?
            IF(WSign.gt.0) THEN
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
!            IF(Iter.eq.18925) THEN
!                WRITE(6,*) "Not Created",rh%v,rat
!            ENDIF
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
                IF(WSign.gt.0) THEN
                    IF(real(rh%v).gt.0.D0) THEN
                        AttemptCreatePar=-ExtraCreate    !Additional particles are negative
                    ELSE
                        AttemptCreatePar=ExtraCreate       !Additional particles are positive
                    ENDIF
                ELSE
                    IF(real(rh%v).gt.0.D0) THEN
                        AttemptCreatePar=ExtraCreate
                    ELSE
                        AttemptCreatePar=-ExtraCreate
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
        
!We know we want to create a particle. Return the bit-representation of this particle (if we have not already got it)
        IF(.not.tHPHF.and.AttemptCreatePar.ne.0) THEN
            CALL FindExcitBitDet(iLutCurr,iLutnJ,IC,Ex,NIfD)
        ENDIF

!        IF(AttemptCreatePar.ne.0) THEN
!            WRITE(6,"(A,F15.5,I5,G25.15,I8,G25.15)") "Outwards ", rat,ExtraCreate,real(rh%v),nParts,Prob
!        ENDIF

        IF(tHistEnergies) THEN
!First histogram off-diagonal matrix elements.
            Bin=INT((real(rh%v,r2)+OffDiagMax)/OffDiagBinRange)+1
            IF(Bin.le.0.or.Bin.gt.iOffDiagNoBins) THEN
                CALL Stop_All("AttemptCreatePar","Trying to histogram off-diagonal matrix elements, but outside histogram array bounds.")
            ENDIF
            IF(IC.eq.1) THEN
                SinglesAttemptHist(Bin)=SinglesAttemptHist(Bin)+(REAL(nParts,r2)*Tau/Prob)
                IF(AttemptCreatePar.ne.0) THEN
                    SinglesHist(Bin)=SinglesHist(Bin)+real(abs(AttemptCreatePar),r2)
                    IF(BRR(Ex(1,1)).le.NEl) THEN
                        IF(BRR(Ex(2,1)).le.NEl) THEN
                            SinglesHistOccOcc(Bin)=SinglesHistOccOcc(Bin)+real(abs(AttemptCreatePar),r2)
                        ELSE
                            SinglesHistOccVirt(Bin)=SinglesHistOccVirt(Bin)+real(abs(AttemptCreatePar),r2)
                        ENDIF
                    ELSE
                        IF(BRR(Ex(2,1)).le.NEl) THEN
                            SinglesHistVirtOcc(Bin)=SinglesHistVirtOcc(Bin)+real(abs(AttemptCreatePar),r2)
                        ELSE
                            SinglesHistVirtVirt(Bin)=SinglesHistVirtVirt(Bin)+real(abs(AttemptCreatePar),r2)
                        ENDIF
                    ENDIF
                ENDIF
            ELSEIF(IC.eq.2) THEN
                DoublesAttemptHist(Bin)=DoublesAttemptHist(Bin)+(REAL(nParts,r2)*Tau/Prob)
                IF(AttemptCreatePar.ne.0) THEN
                    DoublesHist(Bin)=DoublesHist(Bin)+real(abs(AttemptCreatePar),r2)
                ENDIF
            ENDIF

            IF(tHPHF) THEN
                CALL HPHFGetDiagHElement(nJ,iLutnJ,rh)
            ELSE
                rh=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
            ENDIF
            Bin=INT((real(rh%v,r2)-Hii)/BinRange)+1
            IF(Bin.gt.iNoBins) THEN
                CALL Stop_All("AttemptCreatePar","Histogramming energies higher than the arrays can cope with. Increase iNoBins or BinRange")
            ENDIF
            IF(AttemptCreatePar.ne.0) THEN
                SpawnHist(Bin)=SpawnHist(Bin)+real(abs(AttemptCreatePar),r2)
!                WRITE(6,*) "Get Here!", real(abs(AttemptCreatePar),r2),Bin
            ENDIF
            AttemptHist(Bin)=AttemptHist(Bin)+(REAL(nParts,r2)*Tau/Prob)
        ENDIF

    END FUNCTION AttemptCreatePar


! This function is based on attemptcreatepar, however it only attempts to create particles back on a parent determinant from which it
! was spawned.
! It decides whether or not we are going to create a child back on that parent.  It returns 0 for no spawning, and +1/-1 for a child with sign.
    INTEGER FUNCTION AttemptCreateParBack(iLutCurr,iLutParent,WSign,rh,nParts,tMinorDetList)
        use GenRandSymExcitNUMod , only : GenRandSymExcitBiased
        INTEGER :: StoreNumTo,StoreNumFrom,DetLT,i,ExtraCreate,Ex(2,2),WSign,nParts
        INTEGER :: iLutCurr(0:NIfD),Bin,iLutParent(0:NIfD),PartInd,ExcitLev,IC
        LOGICAL :: SymAllowed,tSuccess,tMinorDetList
        REAL*8 :: Prob=1.D0,r,rat
        TYPE(HElement) :: rh,rhcheck

        IF(tSpawnDominant) THEN
!We only allow spawning at determinants between iMinDomLev and iMaxDomLev which are in an allowed list of dominant determinants.
!Find the excitation level of the excitation
!We need to calculate the bit-representation of this new child. This can be done easily since the ExcitMat is known.
            tMinorDetList=.false.
            CALL FindBitExcitLevel(iLutHF,iLutParent,NIfD,ExcitLev,iMaxDomLev)
            IF((ExcitLev.le.iMaxDomLev).and.(ExcitLev.ge.iMinDomLev)) THEN
!We now need to binary search the list of allowed determinants to see whether it is in the list or not.
                CALL BinSearchParts3(iLutParent,DomDets,iNoDomDets,DomExcIndex(ExcitLev),DomExcIndex(ExcitLev+1)-1,PartInd,tSuccess)
                IF(.not.tSuccess) THEN
!If the particle is not in the list, then do not allow a spawning event there, unless tMinorDetsStar is on.
                    IF(tMinorDetsStar) THEN
!Put into the MinorSpawnDets array, with its parent (iLutCurr)
                        tMinorDetList=.true.
                    ELSE
                        AttemptCreateParBack=0
                        RETURN
                    ENDIF
                ENDIF
            ENDIF
        ENDIF
                

!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
        rat=Tau*abs(rh%v)*REAL(nParts,r2)/Prob

!If probability is > 1, then we can just create multiple children at the chosen determinant
        ExtraCreate=INT(rat)
        rat=rat-REAL(ExtraCreate)


!Stochastically choose whether to create or not according to ranlux 
        IF(tMerTwist) THEN
            CALL genrand_real2(r) 
        ELSE
            CALL RANLUX(r,1)
        ENDIF
        IF(rat.gt.r) THEN
!Child is created - what sign is it?
            IF(WSign.gt.0) THEN
!Parent particle is positive
                IF(real(rh%v).gt.0.D0) THEN
                    AttemptCreateParBack=-1     !-ve walker created
                ELSE
                    AttemptCreateParBack=1      !+ve walker created
                ENDIF

            ELSE
!Parent particle is negative
                IF(real(rh%v).gt.0.D0) THEN
                    AttemptCreateParBack=1      !+ve walker created
                ELSE
                    AttemptCreateParBack=-1     !-ve walker created
                ENDIF
            ENDIF

        ELSE
!No child particle created
            AttemptCreateParBack=0
        ENDIF

        IF(ExtraCreate.ne.0) THEN
!Need to include the definitely create additional particles from a initial probability > 1

            IF(AttemptCreateParBack.lt.0) THEN
!In this case particles are negative
                AttemptCreateParBack=AttemptCreateParBack-ExtraCreate
            ELSEIF(AttemptCreateParBack.gt.0) THEN
!Include extra positive particles
                AttemptCreateParBack=AttemptCreateParBack+ExtraCreate
            ELSEIF(AttemptCreateParBack.eq.0) THEN
!No particles were stochastically created, but some particles are still definatly created - we need to determinant their sign...
                IF(WSign.gt.0) THEN
                    IF(real(rh%v).gt.0.D0) THEN
                        AttemptCreateParBack=-1*ExtraCreate    !Additional particles are negative
                    ELSE
                        AttemptCreateParBack=ExtraCreate       !Additional particles are positive
                    ENDIF
                ELSE
                    IF(real(rh%v).gt.0.D0) THEN
                        AttemptCreateParBack=ExtraCreate
                    ELSE
                        AttemptCreateParBack=-1*ExtraCreate
                    ENDIF
                ENDIF
            ENDIF
        ENDIF

        RETURN

    END FUNCTION AttemptCreateParBack



!This is a function which tells us whether to annihilate a particle on a determinant if it is the only one there.
!Hopefully this will simulate the annihilation rate to some extent and stop the shift from becoming too negative.
!The only variable it relies on is the approximate density of particles for the given excitation level - ExcitDensity
    LOGICAL FUNCTION AttemptLocalAnn(ExcitDensity)
        REAL*8 :: ExcitDensity,r,AnnProb

!The function can initially be a simple Tau*EXP(-Lambda*ExcitDensity)
!        AnnProb=Tau*EXP(-Lambda*ExcitDensity)
!        AnnProb=Lambda/ExcitDensity
        AnnProb=Lambda/ExcitDensity
        IF(tMerTwist) THEN
            CALL genrand_real2(r) 
        ELSE
            CALL RANLUX(r,1)
        ENDIF
        IF(r.lt.AnnProb) THEN
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
!Multiple particles can be attempted to die at the same time - here, |WSign| > 1 and the probability of a death will be multiplied by |WSign|
    INTEGER FUNCTION AttemptDiePar(DetCurr,Kii,IC,WSign)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),iKill,IC,WSign
!        TYPE(HElement) :: rh,rhij
        REAL*8 :: r,rat,Kii
        LOGICAL :: tDETinCAS


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
        ELSEIF(tFixShiftKii.and.(.not.TSinglePartPhase)) THEN
            IF(Kii.le.FixedKiiCutoff) THEN
                rat=Tau*(Kii-FixShift)
            ELSE
                rat=Tau*(Kii-DiagSft)
            ENDIF
        ELSEIF(tFixCASShift.and.(.not.TSinglePartPhase)) THEN
! The 'TestifDETinCAS' function finds out if the determinant is in the complete active space or not.  If it is, the shift is fixed at 
! the chosen fixed value (FixShift), if not the shift remains as the changing DiagSft value.
           
            tDETinCAS=TestifDETinCAS(DetCurr)
            
            IF(tDETinCAS) THEN
                rat=Tau*(Kii-FixShift)
            ELSE
                rat=Tau*(Kii-DiagSft)
            ENDIF
        ELSE
            rat=Tau*(Kii-DiagSft)
        ENDIF

!If there are multiple particles, decide how many to kill in total...
        rat=rat*abs(WSign)

        iKill=INT(rat)
        rat=rat-REAL(iKill)

!Stochastically choose whether to die or not
        IF(tMerTwist) THEN
            CALL genrand_real2(r) 
        ELSE
            CALL RANLUX(r,1)
        ENDIF
        IF(abs(rat).gt.r) THEN
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
   
    


!This is the heart of FCIMC, where the MC Cycles are performed. However, this includes the 'inward spawning' attempt.
    SUBROUTINE MultipleConnFCIMCycPar()
!        use CalcData , only : iDetGroup
        use Determinants , only : GetHElement3
!        use HPHFRandExcitMod , only : TestGenRandHPHFExcit 
        INTEGER :: nStore(6),VecSlot,i,j,k,l,ValidSpawned,CopySign,ParticleWeight,Loop,iPartBloom
        INTEGER :: nJ(NEl),ierr,IC,Child,iCount,DetCurr(NEl),iLutnJ(0:NIfD)
        REAL*8 :: Prob,rat,HDiag,HDiagCurr,r,HSum
        INTEGER :: iDie,WalkExcitLevel,iMaxExcit,ExcitLength,PartInd,iExcit
        INTEGER :: ExcitLevel,TotWalkersNew,iGetExcitLevel_2,error,length,temp,Ex(2,2),WSign,p,Scratch1(2,nSymLabels),Scratch2(2,nSymLabels)
        LOGICAL :: tParity,DetBitEQ,tMainArr,tFilled,tSuccess,tMinorDetList
        TYPE(HElement) :: HDiagTemp,HElemTemp
        INTEGER , ALLOCATABLE :: ExcitGenTemp(:)

        IF(TDebug) THEN
            WRITE(11,*) Iter,TotWalkers,NoatHF,NoatDoubs,MaxIndex,TotParts
!            CALL FLUSH(11)
        ENDIF
        
        CALL set_timer(Walker_Time,30)
        
!VecSlot indicates the next free position in NewDets
        VecSlot=1
!Reset number at HF and doubles
        NoatHF=0
        NoatDoubs=0
        iPartBloom=0
        ValidSpawned=1  !This is for rotoannihilation - this is the number of spawned particles (well, one more than this.)
!        Combs=1
!        do j=1,iDetGroup-1
!!Combs is the total number of pair combinations
!            Combs=Combs*(TotWalkers-j)
!        enddo

        do j=1,TotWalkers
!Spawning loop done separately to the death loop
!First, decode the bit-string representation of the determinant the walker is on, into a string of naturally-ordered integers
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j),NEl,NIfD)

!Also, we want to find out the excitation level - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.
            IF(tTruncSpace) THEN
!We need to know the exact excitation level for truncated calculations.
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,NEl)
            ELSE
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,2)
            ENDIF

            tFilled=.false.     !This is for regenerating excitations from the same determinant multiple times. There will be a time saving if we can store the excitation generators temporarily.
            IF(tSpawnAsDet) THEN
!Here, we spawn all particles on the determinant in one go, by multiplying the probability of spawning by the number of particles on the determinant.
                Loop=1
                ParticleWeight=abs(CurrentSign(j))
            ELSE
!Here, we spawn each particle on the determinant in a seperate attempt.
                Loop=abs(CurrentSign(j))
                ParticleWeight=1
            ENDIF

            do p=1,Loop
!If rotoannihilating, we are simply looping over all the particles on the determinant

                IF(.not.tImportanceSample) THEN
                    IF(tNonUniRandExcits) THEN
!This will only be a help if most determinants are multiply occupied.
                        IF(tHPHF) THEN
!                            CALL GenRandHPHFExcit(DetCurr,CurrentDets(:,j),nJ,iLutnJ,pDoubles,exFlag,Prob)
                            CALL GenRandHPHFExcit2Scratch(DetCurr,CurrentDets(:,j),nJ,iLutnJ,pDoubles,exFlag,Prob,Scratch1,Scratch2,tFilled,tGenMatHEl)
                        ELSE
                            CALL GenRandSymExcitScratchNU(DetCurr,CurrentDets(:,j),nJ,pDoubles,IC,Ex,tParity,exFlag,Prob,Scratch1,Scratch2,tFilled)
                        ENDIF
                    ELSE
                        CALL GetPartRandExcitPar(DetCurr,CurrentDets(:,j),nJ,IC,0,Prob,iCount,WalkExcitLevel,Ex,tParity)
                    ENDIF
                ENDIF
                
                ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,2)
!                WRITE(6,*) ExcitLevel

                IF((ExcitLevel.eq.2).or.(ExcitLevel.eq.1)) THEN
!Do not allow spawning at doubles
                    Child=0
                ELSE

!Calculate number of children to spawn
                    IF(TTruncSpace) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.
                        IF(tImportanceSample) CALL Stop_All("PerformFCIMCyc","Truncated calculations not yet working with importance sampling")

                        IF(ExcitLevel.gt.ICILevel) THEN
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                            Child=0
                        ELSE
                            Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                        ENDIF
                    ELSE
!SD Space is not truncated - allow attempted spawn as usual
                        Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                    ENDIF
                ENDIF
                
                IF(Child.ne.0) THEN
!We want to spawn a child - find its information to store

                    NoBorn=NoBorn+abs(Child)     !Update counter about particle birth
                    IF(IC.eq.1) THEN
                        SpawnFromSing=SpawnFromSing+abs(Child)
                    ENDIF

                    IF(abs(Child).gt.25) THEN
!If more than 25 particles are created in one go, then log this fact and print out later that this has happened.
                        IF(abs(Child).gt.abs(iPartBloom)) THEN
                            IF(IC.eq.1) THEN
                                iPartBloom=-abs(Child)
                            ELSE
                                iPartBloom=abs(Child)
                            ENDIF
                        ENDIF
                    ENDIF

!We need to calculate the bit-representation of this new child. This can be done easily since the ExcitMat is known.
                    IF(.not.tHPHF) CALL FindExcitBitDet(CurrentDets(:,j),iLutnJ,IC,Ex,NIfD)

                    IF(tRotoAnnihil) THEN
!In the RotoAnnihilation implimentation, we spawn particles into a seperate array - SpawnedParts and SpawnedSign. 
!The excitation level and diagonal matrix element are also found out after the annihilation.
!Cannot use old excitation generators with rotoannihilation.

!In rotoannihilation, we can specify multiple particles on the same entry. 
                        SpawnedParts(:,ValidSpawned)=iLutnJ(:)
                        SpawnedSign(ValidSpawned)=Child
                        ValidSpawned=ValidSpawned+1     !Increase index of spawned particles

                    ELSE
!Calculate diagonal ham element
                        CALL Stop_All("MultipleConnFcimCycPar","Needs Rotoannihilation")
                    
                    ENDIF   !Endif rotoannihil

                    Acceptances=Acceptances+ABS(Child)      !Sum the number of created children to use in acceptance ratio
                
                ENDIF   !End if child created

            enddo   !End of cycling over mulitple particles on same determinant.

        enddo

        do j=1,NoDoubs

            DetCurr(:)=DoublesDets(:,j)
!Setup excit generators for this double excitation 
            iMaxExcit=0
            nStore(1:6)=0
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitLength,nJ,iMaxExcit,0,nStore,3)
            ALLOCATE(ExcitGenTemp(ExcitLength),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("MultipleConnFcimCycPar","Problem allocating excitation generator")
            ExcitGenTemp(1)=0
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGenTemp,nJ,iMaxExcit,0,nStore,3)

            HSum=0.D0

            do while(.true.)
!Generate double excitations
                CALL GenSymExcitIt2(Detcurr,NEl,G1,nBasis,nBasisMax,.false.,ExcitGenTemp,nJ,iExcit,0,nStore,3)
                IF(nJ(1).eq.0) EXIT

!Find matrix element
                HElemTemp=GetHElement3(DetCurr,nJ,iExcit)
                IF((abs(REAL(HElemTemp%v,r2))).gt.1.D-8) THEN

!Encode this determinant
                    CALL EncodeBitDet(nJ,iLutnJ,NEl,NIfD)

!If matrix element above a certain size, then find whether the determinant is in the main list.
                    CALL BinSearchParts(iLutnJ,1,TotWalkers,PartInd,tSuccess)

                    IF(tSuccess) THEN
!If found, find Ni x Hij and attempt to spawn at j
                        HSum=HSum+REAL(HElemTemp%v,r2)*CurrentSign(PartInd)
!                        WRITE(6,*) REAL(HElemTemp%v,r2),CurrentSign(PartInd)
                    ENDIF
                ENDIF

            enddo

            DEALLOCATE(ExcitGenTemp)

!The sum of Hij*cj is now in HSum.
            rat=Tau*abs(HSum)
            Child=INT(rat)
            rat=rat-REAL(Child)
            CALL genrand_real2(r)
            IF(rat.gt.r) THEN
!Create child at DetCurr 
                Child=Child+1
            ENDIF
            IF(HSum.gt.0.D0) THEN
!Create negative child
                Child=-Child
            ENDIF

            IF(Child.ne.0) THEN
!                WRITE(6,"(A,F15.5,I5,G20.10,I8,G25.15)") "Inwards ", rat,Child,HSum
!We want to spawn a child - find its information to store
                CALL EncodeBitDet(DetCurr,iLutnJ,NEl,NIfD)

                NoBorn=NoBorn+abs(Child)     !Update counter about particle birth

                IF(abs(Child).gt.25) THEN
!If more than 25 particles are created in one go, then log this fact and print out later that this has happened.
                    WRITE(6,"(I9,A)") Child," particles created in one go during inwards spawning step"
                ENDIF

                IF(tRotoAnnihil) THEN
!In the RotoAnnihilation implimentation, we spawn particles into a seperate array - SpawnedParts and SpawnedSign. 
!The excitation level and diagonal matrix element are also found out after the annihilation.
!Cannot use old excitation generators with rotoannihilation.

!In rotoannihilation, we can specify multiple particles on the same entry. 
                    SpawnedParts(:,ValidSpawned)=iLutnJ(:)
                    SpawnedSign(ValidSpawned)=Child
                    ValidSpawned=ValidSpawned+1     !Increase index of spawned particles

                ELSE
                    CALL Stop_All("LIUB","Die horribly")
                ENDIF   !Endif rotoannihil

                Acceptances=Acceptances+ABS(Child)      !Sum the number of created children to use in acceptance ratio
            
            ENDIF   !End if child created

        enddo


        do j=1,TotWalkers
!Death loop done separately to the spawning loop
!We also have to calculate the energy properties here.
!First, decode the bit-string representation of the determinant the walker is on, into a string of naturally-ordered integers
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j),NEl,NIfD)

!Also, we want to find out the excitation level - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
            CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,2)
            HDiagCurr=CurrentH(j)

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
            CALL SumEContrib(DetCurr,WalkExcitLevel,CurrentSign(j),CurrentDets(:,j),HDiagCurr,1.D0)

!We now have to decide whether the parent particle (j) wants to self-destruct or not...
!For rotoannihilation, we can have multiple particles on the same determinant - these can be stochastically killed at the same time.

            iDie=AttemptDiePar(DetCurr,HDiagCurr,WalkExcitLevel,CurrentSign(j))
            NoDied=NoDied+iDie          !Update death counter

!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births
            IF(tRotoAnnihil) THEN
!We slot the particles back into the same array and position VecSlot if the particle survives. If it dies, then j increases, moving onto the next
!entry, but VecSlot remains where it is, meaning that j should never be less that VecSlot

                IF(CurrentSign(j).le.0) THEN
                    CopySign=CurrentSign(j)+iDie    !Copy sign is the total number of particles x sign that we want to copy accross.
                    IF(CopySign.le.0) THEN
!If we are copying to the main array, we have to ensure that we maintain sign-coherence in the array. Therefore, if we are spawning anti-particles,
!it wants to go in the spawning array, rather than the main array, so it has a chance to annihilate.
                        tMainArr=.true.
                    ELSE
                        tMainArr=.false.
                    ENDIF
                ELSE
                    CopySign=CurrentSign(j)-iDie
                    IF(CopySign.ge.0) THEN
                        tMainArr=.true.
                    ELSE
                        tMainArr=.false.
                    ENDIF
                ENDIF

                IF(tMainArr.and.(VecSlot.gt.j)) THEN
!We only have a single array, therefore surviving particles are simply transferred back into the original array.
!However, this can not happen if we want to overwrite particles that we haven't even got to yet.
!However, this should not happen under normal circumstances...
                    tMainArr=.false.
                    CALL Stop_All("PerformFCIMCyc","About to overwrite particles which haven't been reached yet. This should not happen in normal situations.")
                ENDIF

                IF(CopySign.ne.0) THEN

                    IF(tMainArr) THEN
!                        NewDets(:,VecSlot)=CurrentDets(:,j)
!                        NewSign(VecSlot)=CopySign
                        CurrentDets(:,VecSlot)=CurrentDets(:,j)
                        CurrentSign(VecSlot)=CopySign
                        IF(.not.tRegenDiagHEls) CurrentH(VecSlot)=HDiagCurr
                        VecSlot=VecSlot+1
                    ELSE
!                        CALL Stop_All("PerformFCIMCyc","Creating anti-particles")
!Here, we are creating anti-particles, and so to keep the sign-coherence of the main array assured, we transfer them to the spawning array.
!This should generally not happen under normal circumstances.
                        do p=1,abs(CopySign)
!In rotoannihilation, we still want to specify the determinants singly - this may change in the future...
                            SpawnedParts(:,ValidSpawned)=CurrentDets(:,j)
                            IF(CopySign.lt.0) THEN
                                SpawnedSign(ValidSpawned)=-1
                            ELSE
                                SpawnedSign(ValidSpawned)=1
                            ENDIF
                            ValidSpawned=ValidSpawned+1     !Increase index of spawned particles
                        enddo
                    ENDIF

                ENDIF

            ENDIF   !To kill if

!Finish cycling over walkers
        enddo

!SumWalkersCyc calculates the total number of walkers over an update cycle on each process
        SumWalkersCyc=SumWalkersCyc+(INT(TotParts,i2))
!        WRITE(6,*) "Born, Die: ",NoBorn, NoDied

!Since VecSlot holds the next vacant slot in the array, TotWalkers will be one less than this.
        TotWalkersNew=VecSlot-1

!Output if there has been a particle bloom this iteration. A negative number indicates that particles were created from a single excitation.
        IF(iPartBloom.ne.0) THEN
            WRITE(6,"(A,I10,A)") "LARGE PARTICLE BLOOMS in iteration ",Iter
            IF(iPartBloom.gt.0) THEN
                WRITE(6,"(A,I10,A)") "A max of ",abs(iPartBloom)," particles created in one attempt from double excit."
            ELSE
                WRITE(6,"(A,I10,A)") "A max of ",abs(iPartBloom)," particles created in one attempt from single excit."
            ENDIF
        ENDIF


        rat=(TotWalkersNew+0.D0)/(MaxWalkersPart+0.D0)
        IF(rat.gt.0.95) THEN
            WRITE(6,*) "*WARNING* - Number of particles/determinants has increased to over 95% of MaxWalkersPart"
            CALL FLUSH(6)
        ENDIF

        IF(tRotoAnnihil) THEN
            ValidSpawned=ValidSpawned-1     !For rotoannihilation, this is the number of spawned particles.

            rat=(ValidSpawned+0.D0)/(MaxSpawned+0.D0)
            IF(rat.gt.0.9) THEN
                WRITE(6,*) "*WARNING* - Number of spawned particles has reached over 90% of MaxSpawned"
                CALL FLUSH(6)
            ENDIF
        ENDIF
        
        CALL halt_timer(Walker_Time)
        CALL set_timer(Annihil_Time,30)

        IF(tRotoAnnihil) THEN
!This is the rotoannihilation algorithm. The newly spawned walkers should be in a seperate array (SpawnedParts) and the other list should be ordered.

            CALL RotoAnnihilation(ValidSpawned,TotWalkersNew)

        ELSEIF(TNoAnnihil) THEN
!However, we now need to swap around the pointers of CurrentDets and NewDets, since this was done previously explicitly in the annihilation routine
            IF(associated(CurrentDets,target=WalkVecDets)) THEN
                CurrentDets=>WalkVec2Dets
                CurrentSign=>WalkVec2Sign
                CurrentH=>WalkVec2H
                CurrentExcits=>WalkVec2Excits
                NewDets=>WalkVecDets
                NewSign=>WalkVecSign
                NewH=>WalkVecH
                NewExcits=>WalkVecExcits
            ELSE
                CurrentDets=>WalkVecDets
                CurrentSign=>WalkVecSign
                CurrentH=>WalkVecH
                CurrentExcits=>WalkVecExcits
                NewDets=>WalkVec2Dets
                NewSign=>WalkVec2Sign
                NewH=>WalkVec2H
                NewExcits=>WalkVec2Excits
            ENDIF

            TotWalkers=TotWalkersNew
            TotParts=TotWalkers

        ELSEIF(TAnnihilonProc) THEN
!In this routine, the particles are just annihilated on their own processors. This means that all simulations are independent and not influenced by each other.
!This means that there is no communication between processors and so should be much faster as the system size increases.

            CALL Stop_All("PerformFCIMCyc","AnnihilonProc has been disabled")
!            CALL AnnihilateonProc(TotWalkersNew)
!            Annihilated=Annihilated+(TotWalkersNew-TotWalkers)
!            TotParts=TotWalkers

        ELSE
!This routine now cancels down the particles with opposing sign on each determinant

            CALL AnnihilatePartPar(TotWalkersNew)
            Annihilated=Annihilated+(TotWalkersNew-TotWalkers)
            TotParts=TotWalkers

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
                IF(TotParts.gt.InitWalkers) THEN
                    WRITE(6,*) "Exiting the single particle growth phase - shift can now change"
                    TSinglePartPhase=.false.
                    VaryShiftIter=Iter
                ENDIF
            ENDIF
!Broadcast the fact that TSinglePartPhase may have changed to all processors - unfortunatly, have to do this broadcast every iteration
            CALL MPI_Bcast(TSinglePartPhase,1,MPI_LOGICAL,root,MPI_COMM_WORLD,ierr)
        ELSE

            IF(TotWalkers.gt.(InitWalkers*GrowMaxFactor)) THEN
!Particle number is too large - kill them randomly
                IF(.not.tRotoAnnihil) THEN

!Log the fact that we have made a cull
                    NoCulls=NoCulls+1
                    IF(NoCulls.gt.10) THEN
                        WRITE(6,*) "Too Many Culls"
                        CALL FLUSH(6)
                        call Stop_All("PerformFCIMCyc","Too Many Culls")
                    ENDIF
!CullInfo(:,1) is walkers before cull
                    CullInfo(NoCulls,1)=TotParts
!CullInfo(:,3) is MC Steps into shift cycle before cull
                    CullInfo(NoCulls,3)=mod(Iter,StepsSft)

                    WRITE(6,"(A,F8.2,A)") "Total number of particles has grown to ",GrowMaxFactor," times initial number on this node..."
                    WRITE(6,"(A,I12,A)") "Killing randomly selected particles in cycle ", Iter," in order to reduce total number on this node..."
                    WRITE(6,"(A,F8.2)") "Population on this node will reduce by a factor of ",CullFactor
                    CALL FLUSH(6)
                    CALL ThermostatParticlesPar(.true.)

                ENDIF

            ELSEIF(TotParts.lt.(InitWalkers/2)) THEN
!Particle number is too small - double every particle in its current position
                IF(.not.tRotoAnnihil) THEN
!Log the fact that we have made a cull
                    NoCulls=NoCulls+1
                    IF(NoCulls.gt.10) CALL Stop_All("PerformFCIMCyc","Too Many Culls")
!CullInfo(:,1) is walkers before cull
                    CullInfo(NoCulls,1)=TotParts
!CullInfo(:,3) is MC Steps into shift cycle before cull
                    CullInfo(NoCulls,3)=mod(Iter,StepsSft)
                    
                    WRITE(6,*) "Doubling particle population on this node to increase total number..."
                    CALL ThermostatParticlesPar(.false.)
                ELSE
!                    WRITE(6,*) "Particle number on this node is less than half InitWalkers value"
                ENDIF
            ENDIF
        
        ENDIF

    END SUBROUTINE MultipleConnFCIMCycPar

!This routine will check to see if any of the orbitals in the determinant are in the orbitals which are only to be attached to HF in a 'star'
    SUBROUTINE CheckStarOrbs(DetCurr,tStarDet)
        INTEGER :: DetCurr(NEl),i
        LOGICAL :: tStarDet

        tStarDet=.false.
        do i=NEl,1,-1
            IF(SpinInvBrr(DetCurr(i)).gt.(nBasis-iStarOrbs)) THEN
                tStarDet=.true.
                EXIT
            ENDIF
        enddo

    END SUBROUTINE CheckStarOrbs

!This routine will change the reference determinant to DetCurr. It will also re-zero all the energy estimators, since they now correspond to
!projection onto a different determinant.
    SUBROUTINE ChangeRefDet(HDiagCurr,DetCurr,iLutCurr)
        use Determinants , only : GetH0Element3
        INTEGER :: iLutCurr(0:NIfD),DetCurr(NEl),i,nStore(6),ierr,iMaxExcit
        INTEGER :: nJ(NEl)
        TYPE(HElement) :: TempHii
        REAL*8 :: HDiagCurr

        CALL Stop_All("ChangeRefDet","This option does not currently work. Bug ghb24 if its needed")
!Problem is that we need to rerun the simulation from scratch, and particles currently in the simulation will keep on
!changing the reference since their diagonal K element will remain negative.

        do i=1,NEl
            HFDet(i)=DetCurr(i)
        enddo
        Hii=Hii+HDiagCurr
        iLutHF(0:NIfD)=iLutCurr(0:NIfD)
        TempHii=GetH0Element3(HFDet)
        Fii=REAL(TempHii%v,r2)
        HFHash=CreateHash(HFDet)
        DEALLOCATE(HFExcit%ExcitData)
!Setup excitation generator for the HF determinant. If we are using assumed sized excitgens, this will also be assumed size.
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,HFExcit%nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
        ALLOCATE(HFExcit%ExcitData(HFExcit%nExcitMemLen),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("ChangeRefDet","Problem allocating excitation generator")
        HFExcit%ExcitData(1)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,HFExcit%ExcitData,nJ,iMaxExcit,0,nStore,3)
        HFExcit%nPointed=0

        WRITE(6,"(A)") "*** Changing the reference determinant ***"
        WRITE(6,"(A)") "Switching reference and zeroing energy counters"
        WRITE(6,"(A,G25.16)") "New reference energy is: ",Hii
        WRITE(6,"(A)",advance='no') "New Determinant is: "
        do i=1,NEL-1
            WRITE(6,"(I4)",advance='no') HFDet(i)
        enddo
        WRITE(6,"(I4)") HFDet(NEl)
        
        AllSumENum=0.D0
        AllSumNoatHF=0.D0
        AllNoatHF=0
        AllNoatDoubs=0
!        AllDetsNorm=0.D0
        AllHFCyc=0.D0
        AllENumCyc=0.D0
        ProjectionE=0.D0
        SumENum=0.D0
        NoatHF=0
        NoatDoubs=0
!        DetsNorm=0.D0
        HFCyc=0
        HFPopCyc=0
        ENumCyc=0.D0
        ProjEIter=0.D0
        ProjEIterSum=0.D0

    END SUBROUTINE ChangeRefDet
    
! Simulate Magnetization to break sign symmetry. Death probabiliy is now sign-dependent and given by tau*(Kii-S-Bs_i). B on HF is +ve and sign of particle is s_i.
! This routine simply calculates whether the determinant is magnetic, and how the energy is therefore shifted.
    SUBROUTINE FindDiagElwithB(HDiag,ExcitLevel,nJ,WSign)
        REAL*8 :: HDiag
        LOGICAL :: MagDet,CompiPath
        INTEGER :: ExcitLevel,nJ(NEl),i,WSign
        TYPE(HElement) :: HDiagTemp

        IF(ExcitLevel.eq.0) THEN
!The HF determinant is definitely magnetic and we define the sign of the particles on HF to want to be +ve to be parallel with B
            IF(WSign.eq.1) THEN
!Particle is alligned with B - energy is lower by a value B
                IF(tSymmetricField) THEN
                    HDiag=-BField
                ENDIF
            ELSE
!Particle is antiparallel to B - energy is raised by B
                HDiag=BField
            ENDIF
        
        ELSEIF(ExcitLevel.eq.2) THEN
!We can only tell the sign of the doubles or HF with any certainty. This information is calculated in FindMagneticDets.
!First need to find out if the determinant is one which has been selected to be magnetic
            MagDet=.false.
            do i=1,NoMagDets-1 !Is NoMagDets-1 since the HF det is always magnetic
                IF(CompiPath(nJ,MagDets(:,i),NEl)) THEN
!Determinant is magnetic - now need to find if parallel or antiparallel
                    MagDet=.true.
                    EXIT
                ENDIF
            enddo

            IF(MagDet) THEN
!Determinant is magnetic - first find what the unperturbed energy of the determinant is.
                HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                HDiag=(REAL(HDiagTemp%v,r2))-Hii
                
!+ and + wants a +ve number to subtract
!- and - wants a +ve number to subtract
!Else wants a -ve number to subtract (particle is antiparallel)
                IF((WSign*MagDetsSign(i)).eq.1) THEN
!                IF(.not.XOR(WSign,MagDetsSign(i))) THEN
!Particle is correctly alligned. If we have tSymmetricField on, then this lowers the energy. Otherwise, it stays the same.
                    IF(tSymmetricField) THEN
                        HDiag=HDiag-BField
                    ENDIF
                ELSE
!Particle is incorrectly alligned.
                    HDiag=HDiag+BField
                ENDIF

            ELSE
!Double excitation is not magnetic, find diagonal element as normal
                HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
                HDiag=(REAL(HDiagTemp%v,r2))-Hii
            ENDIF

        ELSE
!Give the child the same diagonal K-matrix element it would normally have.
            HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
            HDiag=(REAL(HDiagTemp%v,r2))-Hii
        ENDIF

        RETURN

    END SUBROUTINE FindDiagElwithB


!This routine will find the largest weighted MP1 determinants, from which we can construct energy level splitting dependant on the sign.
    
!This routine will move walkers between the processors, in order to balance the number of walkers on each node.
!This could be made slightly faster by using an MPI_Reduce, rather than searching for the min and max by hand...
!For rotoannihilation, this will balance the *Determinants* over the processors, not the particles themselves.
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
            WRITE(6,*) "Moving particles/determinants between nodes in order to balance the load..."
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
                    WRITE(6,*) "Initial range of walkers/dets is: ",WalktoTransfer(4)
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
                
                CALL MPI_Send(CurrentDets(:,IndexFrom:TotWalkers),WalktoTransfer(1)*(NIfD+1),MPI_INTEGER,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
                CALL MPI_Send(CurrentSign(IndexFrom:TotWalkers),WalktoTransfer(1),MPI_INTEGER,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
!                CALL MPI_Send(CurrentIC(IndexFrom:TotWalkers),WalktoTransfer(1),MPI_INTEGER,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
                IF(.not.tRegenDiagHEls) THEN
                    CALL MPI_Send(CurrentH(IndexFrom:TotWalkers),WalktoTransfer(1),MPI_DOUBLE_PRECISION,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
                ENDIF
!It seems like too much like hard work to send the excitation generators accross, just let them be regenerated on the other side...
!However, we do need to indicate that these the excitgens are no longer being pointed at.
                IF(.not.TRegenExcitgens) THEN
                    do l=IndexFrom,TotWalkers
!Run through the list of walkers, and make sure that their excitgen is removed.
                        CALL DissociateExitgen(CurrentExcits(l))
                    enddo
                ENDIF

                IF((.not.TNoAnnihil).and.(.not.tRotoAnnihil)) CALL MPI_Send(HashArray(IndexFrom:TotWalkers),WalktoTransfer(1),MPI_DOUBLE_PRECISION,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
!                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) CALL MPI_Send(HashArray(IndexFrom:TotWalkers),WalktoTransfer(1),MPI_DOUBLE_PRECISION,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
!                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) CALL MPI_Send(HashArray(IndexFrom:TotWalkers),WalktoTransfer(1),mpilongintegertype,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
                
                TotWalkers=TotWalkers-WalktoTransfer(1)         !Update TotWalkers on this node to reflect that we have lost some
            
            ELSEIF(WalkToTransfer(3).eq.iProcIndex) THEN
!This processor wants to receive walkers from WalktoTransfer(2)
                IndexFrom=TotWalkers+1
                IndexTo=TotWalkers+WalktoTransfer(1)
!                WRITE(6,*) "RECEIVING: ",Transfers,WalkToTransfer(:),IndexFrom,IndexTo
!                CALL FLUSH(6)

                CALL MPI_Recv(CurrentDets(:,IndexFrom:IndexTo),WalktoTransfer(1)*(NIfD+1),MPI_INTEGER,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(CurrentSign(IndexFrom:IndexTo),WalktoTransfer(1),MPI_INTEGER,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
!                CALL MPI_Recv(CurrentIC(IndexFrom:IndexTo),WalktoTransfer(1),MPI_INTEGER,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
                IF(.not.tRegenDiagHEls) THEN
                    CALL MPI_Recv(CurrentH(IndexFrom:IndexTo),WalktoTransfer(1),MPI_DOUBLE_PRECISION,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
                ENDIF
                IF((.not.TNoAnnihil).and.(.not.tRotoAnnihil)) CALL MPI_Recv(HashArray(IndexFrom:IndexTo),WalktoTransfer(1),MPI_DOUBLE_PRECISION,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
!                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) CALL MPI_Recv(HashArray(IndexFrom:IndexTo),WalktoTransfer(1),MPI_DOUBLE_PRECISION,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
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

        IF(tRotoAnnihil) THEN
!If we are using rotoannihilation, then we need to maintain sorted lists. There is a much better way to do it than sorting them all again though!...
!We are also storing the list as determinants, rather than particles. Therefore after sorting, we need to compress the list to remove multiple specifications of the same det.
            IF(.not.tRegenDiagHEls) THEN
                CALL SortCompressListswH(TotWalkers,CurrentDets(0:NIfD,1:TotWalkers),CurrentSign(1:TotWalkers),CurrentH(1:TotWalkers))
            ELSE
                CALL SortCompressLists(TotWalkers,CurrentDets(0:NIfD,1:TotWalkers),CurrentSign(1:TotWalkers))
            ENDIF
        ELSE
            TotParts=TotWalkers
        ENDIF

        IF(iProcIndex.eq.root) THEN
            WRITE(6,*) "Transfer of walkers/dets finished. Number of transfers needed: ",Transfers
        ENDIF

        TBalanceNodes=.false.

        RETURN

    END SUBROUTINE BalanceWalkersonProcs

!This routine will take the particle and sign lists, sort them and then compress them so each determinant is only specified once. This requires sign-coherent lists.
!The 'length' will be returned as the length of the new list.
    SUBROUTINE SortCompressLists(Length,PartList,SignList)
        INTEGER :: Length,PartList(0:NIfD,Length),SignList(Length)
        INTEGER :: i,VecInd,DetsMerged
        LOGICAL :: DetBitEQ

        CALL SortBitDets(Length,PartList,NIfD,SignList)

!Now compress the list.
        VecInd=1
        DetsMerged=0
        TotParts=0
        IF(Length.gt.0) THEN
            TotParts=TotParts+abs(SignList(1))
        ENDIF
        do i=2,Length
            TotParts=TotParts+abs(SignList(i))
            IF(.not.DetBitEQ(PartList(0:NIfD,i),PartList(0:NIfD,VecInd),NIfD)) THEN
                VecInd=VecInd+1
                PartList(:,VecInd)=PartList(:,i)
                SignList(VecInd)=SignList(i)
            ELSE
                SignList(VecInd)=SignList(VecInd)+SignList(i)
                DetsMerged=DetsMerged+1
            ENDIF
        enddo

        Length=Length-DetsMerged

!Now go through the list, finding common determinants and combining their sign. Also, recalculate TotParts
!        TotParts=abs(SignList(1))
!        CurrInd=1
!        DetCurr=PartList(0:NIfD,1)
!        i=2
!        do while(i.le.Length)
!            CompParts=DetBitLT(DetCurr,PartList(0:NIfD,i),NIfD)
!            IF(CompParts.eq.-1) THEN
!                CALL Stop_All("SortCompressLists","Lists not correctly sorted...")
!            ELSEIF(CompParts.eq.0) THEN
!                IF((SignList(CurrInd)*SignList(i)).le.0) THEN
!                    CALL Stop_All("SortCompressLists","Compressing list which is not sign-coherent")
!                ELSE
!                    SignList(CurrInd)=SignList(CurrInd)+SignList(i)
!                    TotParts=TotParts+abs(SignList(i))
!!Now compress the list...
!                    do j=i,Length-1
!                        PartList(:,j)=PartList(:,j+1)
!                        SignList(j)=SignList(j+1)
!                    enddo
!!                    PartList(:,i:Length-1)=PartList(:,i+1:Length)
!!                    SignList(i:Length-1)=SignList(i+1:Length)
!                    Length=Length-1
!                ENDIF
!            ELSE
!                DetCurr=PartList(0:NIfD,i)
!                TotParts=TotParts+abs(SignList(i))
!                CurrInd=i
!                i=i+1
!            ENDIF
!        enddo

    END SUBROUTINE SortCompressLists

!This routine will take the particle and sign lists, sort them and then compress them so each determinant is only specified once. This requires sign-coherent lists.
!The 'length' will be returned as the length of the new list.
!In this version, the hamiltonian matrix elements will be fed through with the rest of the list and taken with the particles.
    SUBROUTINE SortCompressListswH(Length,PartList,SignList,HList)
        INTEGER :: Length,PartList(0:NIfD,Length),SignList(Length),j
        REAL*8 :: HList(Length)
        INTEGER :: i,DetsMerged,VecInd
        LOGICAL :: DetBitEQ

        CALL SortBitDetswH(Length,PartList,NIfD,SignList,HList)
!        CALL CheckOrdering(PartList(:,1:Length),SignList(1:Length),Length,.false.)
  
!Now compress...
        VecInd=1
        DetsMerged=0
        TotParts=0
        IF(Length.gt.0) THEN
            TotParts=TotParts+abs(SignList(1))
        ENDIF
        do i=2,Length
            TotParts=TotParts+abs(SignList(i))
            IF(.not.DetBitEQ(PartList(0:NIfD,i),PartList(0:NIfD,VecInd),NIfD)) THEN
                VecInd=VecInd+1
                PartList(:,VecInd)=PartList(:,i)
                SignList(VecInd)=SignList(i)
                HList(VecInd)=HList(i)
            ELSE
                SignList(VecInd)=SignList(VecInd)+SignList(i)
                DetsMerged=DetsMerged+1
            ENDIF
        enddo

        Length=Length-DetsMerged

!Now go through the list, finding common determinants and combining their sign. Also, recalculate TotParts
!        TotParts=abs(SignList(1))
!        CurrInd=1
!        DetCurr=PartList(0:NIfD,1)
!        i=2
!        do while(i.le.Length)
!            CompParts=DetBitLT(DetCurr,PartList(0:NIfD,i),NIfD)
!            IF(CompParts.eq.-1) THEN
!                CALL Stop_All("SortCompressLists","Lists not correctly sorted...")
!            ELSEIF(CompParts.eq.0) THEN
!                IF((SignList(CurrInd)*SignList(i)).le.0) THEN
!                    WRITE(6,*) "******************************************"
!                    do j=MIN(CurrInd-10,1),Max(i+10,Length)
!                        WRITE(6,*) j,PartList(:,j),SignList(j)
!                    enddo
!                    WRITE(6,*) Iter,CurrInd,Length,i,SignList(CurrInd),SignList(i)
!                    CALL FLUSH(6)
!                    CALL Stop_All("SortCompressListswH","Compressing list which is not sign-coherent")
!                ELSE
!                    SignList(CurrInd)=SignList(CurrInd)+SignList(i)
!                    TotParts=TotParts+abs(SignList(i))
!!Now compress the list...
!                    do j=i,Length-1
!                        PartList(:,j)=PartList(:,j+1)
!                        SignList(j)=SignList(j+1)
!                        HList(j)=HList(j+1)
!                    enddo
!!                    PartList(:,i:Length-1)=PartList(:,i+1:Length)
!!                    SignList(i:Length-1)=SignList(i+1:Length)
!!                    HList(i:Length-1)=HList(i+1:Length)
!                    Length=Length-1
!                ENDIF
!            ELSE
!                DetCurr=PartList(0:NIfD,i)
!                TotParts=TotParts+abs(SignList(i))
!                CurrInd=i
!                i=i+1
!            ENDIF
!        enddo

!        CALL CheckOrdering(PartList(:,1:Length),CurrentSign(1:Length),Length,.true.)

    END SUBROUTINE SortCompressListswH

    
    SUBROUTINE DeallocFCIMCMemPar()
        INTEGER :: i,error,length,temp
        CHARACTER(len=*), PARAMETER :: this_routine='DeallocFciMCMemPar'
        CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: message

        IF(tPrintSpinCoupHEl) CLOSE(87)
            
        IF(tHistSpawn.or.tCalcFCIMCPsi) THEN
            DEALLOCATE(Histogram)
            DEALLOCATE(AllHistogram)
            IF(tHistSpawn) THEN
                DEALLOCATE(InstHist)
                DEALLOCATE(InstAnnihil)
                DEALLOCATE(AvAnnihil)
            ENDIF
            IF(iProcIndex.eq.0) THEN
                IF(tHistSpawn) THEN
                    DEALLOCATE(AllInstHist)
                    DEALLOCATE(AllAvAnnihil)
                    DEALLOCATE(AllInstAnnihil)
                ENDIF
            ENDIF
        ELSEIF(tHistEnergies) THEN
            DEALLOCATE(Histogram)
            DEALLOCATE(AttemptHist)
            DEALLOCATE(SpawnHist)
            DEALLOCATE(SinglesHist)
            DEALLOCATE(DoublesHist)
            DEALLOCATE(DoublesAttemptHist)
            DEALLOCATE(SinglesAttemptHist)
            DEALLOCATE(SinglesHistOccOcc)
            DEALLOCATE(SinglesHistVirtOcc)
            DEALLOCATE(SinglesHistOccVirt)
            DEALLOCATE(SinglesHistVirtVirt)
            IF(iProcIndex.eq.0) THEN
                DEALLOCATE(AllHistogram)
                DEALLOCATE(AllAttemptHist)
                DEALLOCATE(AllSpawnHist)
                DEALLOCATE(AllSinglesAttemptHist)
                DEALLOCATE(AllSinglesHist)
                DEALLOCATE(AllDoublesAttemptHist)
                DEALLOCATE(AllDoublesHist)
                DEALLOCATE(AllSinglesHistOccOcc)
                DEALLOCATE(AllSinglesHistVirtOcc)
                DEALLOCATE(AllSinglesHistOccVirt)
                DEALLOCATE(AllSinglesHistVirtVirt)
            ENDIF
        ENDIF
        DEALLOCATE(WalkVecDets)
        CALL LogMemDealloc(this_routine,WalkVecDetsTag)
        DEALLOCATE(WalkVecSign)
        CALL LogMemDealloc(this_routine,WalkVecSignTag)
        IF((.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) THEN
            DEALLOCATE(WalkVec2Dets)
            CALL LogMemDealloc(this_routine,WalkVec2DetsTag)
            DEALLOCATE(WalkVec2Sign)
            CALL LogMemDealloc(this_routine,WalkVec2SignTag)
        ENDIF
        IF(.not.tRegenDiagHEls) THEN
            DEALLOCATE(WalkVecH)
            CALL LogMemDealloc(this_routine,WalkVecHTag)
            IF((.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) THEN
                DEALLOCATE(WalkVec2H)
                CALL LogMemDealloc(this_routine,WalkVec2HTag)
            ENDIF
        ENDIF
        
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
        DEALLOCATE(iLutHF)
        DEALLOCATE(HFExcit%ExcitData)
        IF(.not.TRegenExcitgens) THEN
            do i=1,MaxWalkersPart
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

        IF(ALLOCATED(SpinInvBrr)) THEN
            CALL LogMemDealloc(this_routine,SpinInvBRRTag)
            DEALLOCATE(SpinInvBRR)
        ENDIF

        IF(tConstructNOs) THEN
            DEALLOCATE(OneRDM)
            CALL LogMemDealloc(this_routine,OneRDMTag)
        ENDIF
        
        IF(tUseGuide) THEN
            DEALLOCATE(GuideFuncDets)
            CALL LogMemDealloc(this_routine,GuideFuncDetsTag)
            DEALLOCATE(GuideFuncSign)
            CALL LogMemDealloc(this_routine,GuideFuncSignTag)

            DEALLOCATE(DetstoRotate)
            CALL LogMemDealloc(this_routine,DetstoRotateTag)
            DEALLOCATE(SigntoRotate)
            CALL LogMemDealloc(this_routine,SigntoRotateTag)
            DEALLOCATE(DetstoRotate2)
            CALL LogMemDealloc(this_routine,DetstoRotate2Tag)
            DEALLOCATE(SigntoRotate2)
            CALL LogMemDealloc(this_routine,SigntoRotate2Tag)
        ENDIF

        IF(tMinorDetsStar) THEN
            DEALLOCATE(MinorStarDets)
            CALL LogMemDealloc(this_routine,MinorStarDetsTag)
            DEALLOCATE(MinorSpawnDets)
            CALL LogMemDealloc(this_routine,MinorSpawnDetsTag)
            DEALLOCATE(MinorSpawnDets2)
            CALL LogMemDealloc(this_routine,MinorSpawnDets2Tag)
            DEALLOCATE(MinorStarSign)
            CALL LogMemDealloc(this_routine,MinorStarSignTag)
            DEALLOCATE(MinorSpawnSign)
            CALL LogMemDealloc(this_routine,MinorSpawnSignTag)
            DEALLOCATE(MinorSpawnSign2)
            CALL LogMemDealloc(this_routine,MinorSpawnSign2Tag)
            DEALLOCATE(MinorStarParent)
            CALL LogMemDealloc(this_routine,MinorStarParentTag)
            DEALLOCATE(MinorSpawnParent)
            CALL LogMemDealloc(this_routine,MinorSpawnParentTag)
            DEALLOCATE(MinorSpawnParent2)
            CALL LogMemDealloc(this_routine,MinorSpawnParent2Tag)
            DEALLOCATE(MinorStarHii)
            CALL LogMemDealloc(this_routine,MinorStarHiiTag)
            DEALLOCATE(MinorStarHij)
            CALL LogMemDealloc(this_routine,MinorStarHijTag)
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
                WRITE(15,"(I12,G16.7,I9,G16.7,I12,4I11,2G17.9,2I10,G13.5)") Iter+PreviousCycles,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,   &
 &                  AllTotWalkers,AllLocalAnn,AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,ProjEIter,AllNoatHF,AllNoatDoubs,AccRat
                WRITE(6,"(I12,G16.7,I9,G16.7,I12,4I11,2G17.9,2I10,G13.5)") Iter+PreviousCycles,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,    &
 &                  AllTotWalkers,AllLocalAnn,AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,ProjEIter,AllNoatHF,AllNoatDoubs,AccRat

            ELSEIF(tUseGuide) THEN
                WRITE(15,"(I12,G16.7,I10,G16.7,I12,3I11,3G17.9,2I10,G13.5,I12,G13.5,G18.5,I14)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,   &
 &                  INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,ProjEIter,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),  &
 &                  IterTime,REAL(AllSpawnFromSing)/REAL(AllNoBorn),AlliInitGuideParts
                WRITE(6,"(I12,G16.7,I10,G16.7,I12,3I11,3G17.9,2I10,G13.5,I12,G13.5)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,    &
 &                  INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,ProjEIter,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),IterTime

            ELSEIF(tMinorDetsStar) THEN
                WRITE(15,"(I12,G16.7,I10,G16.7,I12,3I11,3G17.9,2I10,G13.5,I12,G13.5,G18.5,I16)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,   &
 &                  INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,ProjEIter,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),  &
 &                  IterTime,REAL(AllSpawnFromSing)/REAL(AllNoBorn),INT(AllNoMinorWalkers,i2),AllMinorAnnihilated
                WRITE(6,"(I12,G16.7,I10,G16.7,I12,3I11,3G17.9,2I10,G13.5,I12,G13.5)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,    &
 &                  INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,ProjEIter,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),IterTime

            ELSE
!                WRITE(15,"(I12,G16.7,I9,G16.7,I12,3I11,3G17.9,2I10,2G13.5,2I6)") Iter+PreviousCycles,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,   &
! &                  AllTotWalkers,AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,ProjEIter,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
!                WRITE(6,"(I12,G16.7,I9,G16.7,I12,3I11,3G17.9,2I10,2G13.5,2I6)") Iter+PreviousCycles,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,    &
! &                  AllTotWalkers,AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,ProjEIter,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
                WRITE(15,"(I12,G16.7,I10,G16.7,I12,3I13,3G17.9,2I10,G13.5,I12,G13.5,G13.5,I10)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,   &
 &                  INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AvDiagSft,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),IterTime,REAL(AllSpawnFromSing)/REAL(AllNoBorn),WalkersDiffProc
                WRITE(6,"(I12,G16.7,I10,G16.7,I12,3I11,3G17.9,2I10,G13.5,I12,G13.5)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,    &
 &                  INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AvDiagSft,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),IterTime
            ENDIF
            
            CALL FLUSH(6)
            CALL FLUSH(15)
            
        ENDIF

        RETURN

    END SUBROUTINE WriteFCIMCStats

    SUBROUTINE WriteHistogramEnergies()
        INTEGER :: error,i
        REAL*8 :: Norm,EnergyBin

        IF(iProcIndex.eq.0) THEN
            AllHistogram(:)=0.D0
            AllAttemptHist(:)=0.D0
            AllSpawnHist(:)=0.D0
            AllDoublesHist(:)=0.D0
            AllDoublesAttemptHist(:)=0.D0
            AllSinglesHist(:)=0.D0
            AllSinglesAttemptHist(:)=0.D0
            AllSinglesHistOccOcc(:)=0.D0
            AllSinglesHistOccVirt(:)=0.D0
            AllSinglesHistVirtOcc(:)=0.D0
            AllSinglesHistVirtVirt(:)=0.D0
        ENDIF
        CALL MPI_Reduce(Histogram,AllHistogram,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(AttemptHist,AllAttemptHist,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SpawnHist,AllSpawnHist,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesHist,AllSinglesHist,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesAttemptHist,AllSinglesAttemptHist,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(DoublesHist,AllDoublesHist,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(DoublesAttemptHist,AllDoublesAttemptHist,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesHistOccOcc,AllSinglesHistOccOcc,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesHistOccVirt,AllSinglesHistOccVirt,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesHistVirtOcc,AllSinglesHistVirtOcc,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesHistVirtVirt,AllSinglesHistVirtVirt,iNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
  
        IF(iProcIndex.eq.0) THEN
            Norm=0.D0
            do i=1,iNoBins
                Norm=Norm+AllHistogram(i)
            enddo
            do i=1,iNoBins
                AllHistogram(i)=AllHistogram(i)/Norm
            enddo
            Norm=0.D0
            do i=1,iNoBins
                Norm=Norm+AllAttemptHist(i)
            enddo
!            WRITE(6,*) "AllAttemptHistNorm = ",Norm
            do i=1,iNoBins
                AllAttemptHist(i)=AllAttemptHist(i)/Norm
            enddo
            Norm=0.D0
            do i=1,iNoBins
                Norm=Norm+AllSpawnHist(i)
            enddo
!            WRITE(6,*) "AllSpawnHistNorm = ",Norm
            do i=1,iNoBins
                AllSpawnHist(i)=AllSpawnHist(i)/Norm
            enddo
            Norm=0.D0
            do i=1,iOffDiagNoBins
                Norm=Norm+AllSinglesAttemptHist(i)
            enddo
!            WRITE(6,*) "AllSinglesAttemptHistNorm = ",Norm
            do i=1,iOffDiagNoBins
                AllSinglesAttemptHist(i)=AllSinglesAttemptHist(i)/Norm
            enddo
            Norm=0.D0
            do i=1,iOffDiagNoBins
                Norm=Norm+AllDoublesHist(i)
            enddo
!            WRITE(6,*) "AllDoublesHistNorm = ",Norm
            do i=1,iOffDiagNoBins
                AllDoublesHist(i)=AllDoublesHist(i)/Norm
            enddo
            Norm=0.D0
            do i=1,iOffDiagNoBins
                Norm=Norm+AllDoublesAttemptHist(i)
            enddo
!            WRITE(6,*) "AllDoublesAttemptHist = ",Norm
            do i=1,iOffDiagNoBins
                AllDoublesAttemptHist(i)=AllDoublesAttemptHist(i)/Norm
            enddo
            Norm=0.D0
            do i=1,iOffDiagNoBins
                Norm=Norm+AllSinglesHist(i)
            enddo
!            WRITE(6,*) "AllSinglesHistNorm = ",Norm
            do i=1,iOffDiagNoBins
                AllSinglesHist(i)=AllSinglesHist(i)/Norm
            enddo
 
!            Norm=0.D0
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistOccOcc(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistOccOcc(i)=AllSinglesHistOccOcc(i)/Norm
            enddo
!            Norm=0.D0
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistOccVirt(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistOccVirt(i)=AllSinglesHistOccVirt(i)/Norm
            enddo
!            Norm=0.D0
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistVirtOcc(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistVirtOcc(i)=AllSinglesHistVirtOcc(i)/Norm
            enddo
!            Norm=0.D0
!            do i=1,iOffDiagNoBins
!                Norm=Norm+AllSinglesHistVirtVirt(i)
!            enddo
            do i=1,iOffDiagNoBins
                AllSinglesHistVirtVirt(i)=AllSinglesHistVirtVirt(i)/Norm
            enddo
 
            OPEN(17,FILE='EVERYENERGYHIST',STATUS='UNKNOWN')
            OPEN(18,FILE='ATTEMPTENERGYHIST',STATUS='UNKNOWN')
            OPEN(19,FILE='SPAWNENERGYHIST',STATUS='UNKNOWN')

            EnergyBin=BinRange/2.D0
            do i=1,iNoBins
                WRITE(17,*) EnergyBin, AllHistogram(i)
                WRITE(18,*) EnergyBin, AllAttemptHist(i)
                WRITE(19,*) EnergyBin, AllSpawnHist(i)
                EnergyBin=EnergyBin+BinRange
            enddo
            CLOSE(17)
            CLOSE(18)
            CLOSE(19)
            OPEN(20,FILE='SINGLESHIST',STATUS='UNKNOWN')
            OPEN(21,FILE='ATTEMPTSINGLESHIST',STATUS='UNKNOWN')
            OPEN(22,FILE='DOUBLESHIST',STATUS='UNKNOWN')
            OPEN(23,FILE='ATTEMPTDOUBLESHIST',STATUS='UNKNOWN')
            OPEN(24,FILE='SINGLESHISTOCCOCC',STATUS='UNKNOWN')
            OPEN(25,FILE='SINGLESHISTOCCVIRT',STATUS='UNKNOWN')
            OPEN(26,FILE='SINGLESHISTVIRTOCC',STATUS='UNKNOWN')
            OPEN(27,FILE='SINGLESHISTVIRTVIRT',STATUS='UNKNOWN')

            EnergyBin=-OffDiagMax+OffDiagBinRange/2.D0
            do i=1,iOffDiagNoBins
                WRITE(20,*) EnergyBin, AllSinglesHist(i)
                WRITE(21,*) EnergyBin, AllSinglesAttemptHist(i)
                WRITE(22,*) EnergyBin, AllDoublesHist(i)
                WRITE(23,*) EnergyBin, AllDoublesAttemptHist(i)
                WRITE(24,*) EnergyBin, AllSinglesHistOccOcc(i)
                WRITE(25,*) EnergyBin, AllSinglesHistOccVirt(i)
                WRITE(26,*) EnergyBin, AllSinglesHistVirtOcc(i)
                WRITE(27,*) EnergyBin, AllSinglesHistVirtVirt(i)
                EnergyBin=EnergyBin+OffDiagBinRange
!                WRITE(6,*) i
            enddo

            CLOSE(20)
            CLOSE(21)
            CLOSE(22)
            CLOSE(23)
            CLOSE(24)
            CLOSE(25)
            CLOSE(26)
            CLOSE(27)
        ENDIF

    END SUBROUTINE WriteHistogramEnergies



    SUBROUTINE InitGuidingFunction()
!This routine reads in the guiding function from the GUIDINGFUNC file printed in a previous calculation.
!It then scales the number of walkers on each determinant up so that the HF population is that specified in the input for iInitGuideParts. 
!The result is an array of determinats and a corresponding array of populations (with sign) for the guiding function.
        INTEGER :: i,j,ierr,CurrentGuideParts,NewGuideParts,error,ExcitLevel,DoubDet(NEl),HFPop,PartInd
        CHARACTER(len=*), PARAMETER :: this_routine='InitGuidingFunction'
        TYPE(HElement) :: HDoubTemp
        REAL*8 :: Hdoub
        LOGICAL :: DetsEq,DetBitEQ,tSuccess


        iGuideDets=0
        AlliInitGuideParts=0
        DetsEq=.false.
        IF(iProcIndex.eq.Root) THEN
            OPEN(36,FILE='GUIDINGFUNC',Status='old')
            READ(36,*) iGuideDets 
        ENDIF
 
        CALL MPI_Bcast(iGuideDets,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)

        ALLOCATE(GuideFuncDets(0:NIfD,1:iGuideDets),stat=ierr)
        CALL LogMemAlloc('GuideFuncDets',(NIfD+1)*iGuideDets,4,this_routine,GuideFuncDetsTag,ierr)
        ALLOCATE(GuideFuncSign(0:iGuideDets),stat=ierr)
        CALL LogMemAlloc('GuideFuncSign',iGuideDets+1,4,this_routine,GuideFuncSignTag,ierr)

        ALLOCATE(DetstoRotate(0:NIfD,1:iGuideDets),stat=ierr)
        CALL LogMemAlloc('DetstoRotate',(NIfD+1)*iGuideDets,4,this_routine,DetstoRotateTag,ierr)
        ALLOCATE(SigntoRotate(0:iGuideDets),stat=ierr)
        CALL LogMemAlloc('SigntoRotate',iGuideDets+1,4,this_routine,SigntoRotateTag,ierr)
        ALLOCATE(DetstoRotate2(0:NIfD,1:iGuideDets),stat=ierr)
        CALL LogMemAlloc('DetstoRotate2',(NIfD+1)*iGuideDets,4,this_routine,DetstoRotate2Tag,ierr)
        ALLOCATE(SigntoRotate2(0:iGuideDets),stat=ierr)
        CALL LogMemAlloc('SigntoRotate2',iGuideDets+1,4,this_routine,SigntoRotate2Tag,ierr)


        IF(iProcIndex.eq.Root) THEN
            !Set up the determinant and sign arrays by reading in from the GUIDINGFUNC file.
            j=1
            do while (j.le.iGuideDets)
                READ(36,*) GuideFuncDets(0:NIfD,j),GuideFuncSign(j)
                j=j+1
            enddo
            CLOSE(36)

            !Calculate the total number of particles in the GUIDINGFUNC file.
            !CurrentGuideParts=0
            !do j=1,iGuideDets
                !CurrentGuideParts=CurrentGuideParts+ABS(GuideFuncSign(j))
            !enddo

            !Scale up the populations (sign), by the ratio of the original total to the iInitGuideParts value from the NECI input,
            !and calculate the new sum (should be approximately iInitGuideParts, but maybe not exactly because need to take nearest integer.
            !NewGuideParts=0
            !do j=1,iGuideDets
                !GuideFuncSign(j)=NINT(REAL(GuideFuncSign(j))*(REAL(iInitGuideParts)/REAL(CurrentGuideParts)))
                !NewGuideParts=NewGuideParts+ABS(GuideFuncSign(j))
            !enddo
            !iInitGuideParts=NewGuideParts
            !The total number of particles in the guiding function is now iInitGuideParts
            !WRITE(6,*) 'iInitGuideParts, ',iInitGuideParts

            !Find the HF determinant and the population in the guiding function.  Take the number from the input and find the factor by which
            !the current pop needs to be scaled by to reach the input value.  Scale all the other populations by the same value.
            !First binary search the Guiding function determinants for the HFDet.

            CALL BinSearchGuideParts(iLutHF,1,iGuideDets,PartInd,tSuccess)
            IF(tSuccess) THEN
                GuideFuncSign(0)=PartInd
                HFPop=ABS(GuideFuncSign(PartInd))
            ELSE
                CALL Stop_All(this_routine, 'HF determinant not found in the guiding function.')
            ENDIF
            NewGuideParts=0
            do j=1,iGuideDets
                GuideFuncSign(j)=NINT(REAL(GuideFuncSign(j))*(REAL(iInitGuideParts)/REAL(HFPop)))
                NewGuideParts=NewGuideParts+ABS(GuideFuncSign(j))
            enddo
            iInitGuideParts=NewGuideParts
            

            
            !Maybe put in a test here that checks the determinants are in the right order.  But they should be cos we're printing them out right.

            AlliInitGuideParts=iInitGuideParts*nProcessors

            !Want to know what is happening to the guiding function throughout the spawning, so just print it out at the beginning to compare to the
            !final one.
            !Store the initguidefuncsigns 
            ALLOCATE(InitGuideFuncSign(iGuideDets),stat=ierr)
            CALL LogMemAlloc('InitGuideFuncSign',iGuideDets,4,this_routine,InitGuideFuncSignTag,ierr)
            InitGuideFuncSign(:)=0

!            OPEN(37,file='GUIDINGFUNCinit',status='unknown')
!            WRITE(37,*) iGuideDets,' determinants included in the guiding function.'    
            do j=1,iGuideDets
                InitGuideFuncSign(j)=nProcessors*GuideFuncSign(j)
!                WRITE(37,*) GuideFuncDets(0:NIfD,j),nProcessors*GuideFuncSign(j)
            enddo
!            CLOSE(37)
            
        ENDIF

        !Broadcast the guiding function determinants (and signs) to all processors.
        !The total number of walkers in the guiding function is therefore nProcessors*iInitGuideParts.
        CALL MPI_Bcast(GuideFuncDets(0:NIfD,1:iGuideDets),iGuideDets*(NIfD+1),MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL MPI_Bcast(GuideFuncSign(0:iGuideDets),iGuideDets+1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)

        !Run through the guiding function determinants and find the index that contains the HF.
        !Want this known on all processors, so that we can just look up the sign at this position to get the guiding function HF population.
!        do i=1,iGuideDets
!            DetsEq=DetBitEQ(iLutHF,GuideFuncDets(0:NIfD,i),NIfD)
!            IF(DetsEq) THEN
!                GuideFuncHFIndex=i
!                EXIT
!            ENDIF
!        enddo
!        IF(GuideFuncHFIndex.ne.GuideFuncSign(0)) CALL Stop_All(this_routine,'Wrong HF Index')
        GuideFuncHFIndex=GuideFuncSign(0)

        IF(iProcIndex.eq.Root) THEN
            GuideFuncDoub=0.D0
            !Run through all other determinants in the guiding function.  Find out if they are doubly excited.  Find H elements, and multiply by number on that double.
            do i=1,iGuideDets
                CALL FindBitExcitLevel(GuideFuncDets(0:NIfD,i),iLutHF,NIfD,ExcitLevel,2)
                IF(ExcitLevel.eq.2) THEN
                    DoubDet(:)=0
                    CALL DecodeBitDet(DoubDet,GuideFuncDets(0:NIfD,i),NEl,NIfD)
                    HdoubTemp=GetHElement2(HFDet,DoubDet,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
                    HDoub=REAL(HDoubTemp%v,r2)
                    GuideFuncDoub=GuideFuncDoub+(GuideFuncSign(i)*Hdoub)
                ENDIF
            enddo
            WRITE(6,*) 'The energy of the guiding function alone is ,',GuideFuncDoub/(REAL(GuideFuncSign(GuideFuncHFIndex),r2))
            GuideFuncDoub=0.D0
        ENDIF
            

    ENDSUBROUTINE InitGuidingFunction



    SUBROUTINE RotoAnnihilGuidingFunc(ValidSpawned)
! This routine takes the spawned particles (that have already been annihilated with themselves and with the main wavefunction), and tries to 
! annihilate them with the guiding function.
! The guiding function itself may be annihilated also, but does not spawn or die by itself.
! However, if the guiding function is completely annihilated on one processor, the spawned particles must be rotated around the other processors
! to check for possible annihilations there.
        INTEGER :: i,j,n,ValidSpawned,InitNoDetstoRotate,NoDetstoRotate,CombSign,error
        INTEGER :: Stat(MPI_STATUS_SIZE),ExcitLevel,DoubDet(NEl)
        TYPE(HElement) :: HDoubTemp
        REAL*8 :: Hdoub
        LOGICAL :: tRotateSpawnedTemp,tRotateSpawned,tDetinSpawnList,DetBitEQ,DetsEq


        NoDetstoRotate=0
        CombSign=0
        DetstoRotate(:,:)=0
        SigntoRotate(:)=0
        tRotateSpawnedTemp=.false.
        tRotateSpawned=.false.
        DetsEq=.false.


        !First attempt at annihilation, just on the processor the spawned particles are currently on.

        !Run through the determinats that have been spawned on.
        do i=1,ValidSpawned
            IF(SpawnedSign(i).ne.0) THEN
                !Run through the guiding function, checking if this spawned determinant is in there.
                do j=1,iGuideDets
                    DetsEq=.false.
                    !DetsEq is true if the two determinants are equal
                    DetsEq=DetBitEQ(SpawnedParts(0:NIfD,i),GuideFuncDets(0:NIfD,j),NIfD)
                    IF(DetsEq) THEN
                        CombSign=SpawnedSign(i)*GuideFuncSign(j)
                        !IF this is negative, the guiding function annihilates the spawned particles.
                        IF(CombSign.lt.0) THEN

                            IF(ABS(SpawnedSign(i)).gt.ABS(GuideFuncSign(j))) THEN
                                ! Don't want to change sign of guiding function so if there are more particles in the spawned list, just put the 
                                ! guiding function to 0 and leave the remaining spawned to be rotated to other processors.
                                SpawnedSign(i)=SpawnedSign(i)+GuideFuncSign(j)
                                ! Add because these are opposite signs.
                                GuideFuncSign(j)=0
                                ! Then need to rotate the remaining walkers in Spawned list to see if there are walkers in the guiding function
                                ! to annihilate with on other processors.

                                NoDetstoRotate=NoDetstoRotate+1
                                DetstoRotate(0:NIfD,NoDetstoRotate)=SpawnedParts(0:NIfD,i)
                                SigntoRotate(NoDetstoRotate)=SpawnedSign(i)

                            ELSEIF(ABS(SpawnedSign(i)).eq.ABS(GuideFuncSign(j))) THEN
                                SpawnedSign(i)=0
                                GuideFuncSign(j)=0
                            ELSE
                                ! The spawned are all annihilated, and the guiding function is decreased by that number.
                                GuideFuncSign(j)=GuideFuncSign(j)+SpawnedSign(i)
                                SpawnedSign(i)=0
                            ENDIF

                        !IF the combined sign (CombSign) is positive, signs are the same and the spawned particles remain.
                        !Nothing changes, the guiding function is not annihilated, and the spawned remain to be put into the full list later.

                        !If CombSign is 0, there are no walkers on the guiding function (for that processor).
                        !Need to rotate the spawned walker to see if there are any walkers on this determinant to annihilate with.
                        ELSEIF(CombSign.eq.0) THEN
                            NoDetstoRotate=NoDetstoRotate+1
                            DetstoRotate(0:NIfD,NoDetstoRotate)=SpawnedParts(0:NIfD,i)
                            SigntoRotate(NoDetstoRotate)=SpawnedSign(i)
                        ENDIF

                        !If we have found a determinant in the guiding function that matches that in the spawned, can stop searching the guiding 
                        !function, there will be no more matches.
                        EXIT

                    ENDIF
                enddo
            ENDIF
        enddo

        IF(NoDetstoRotate.gt.0) tRotateSpawnedTemp=.true.
        !If NoDetstoRotate is 0, don't even have to worry about the rotating stuff.

        !If tRotateSpawnedTemp is true on any processor, this routine makes tRotateSpawned true on all processors.
        CALL MPI_AllReduce(tRotateSpawnedTemp,tRotateSpawned,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,error)


!The allocated DetstoRotate arrays are as big as iGuideDets (the number of determinants in the guiding function), but will only need to rotate 
!a portion of these (those determinants for which the guiding function has had all its walkers annihilated).

        IF(tRotateSpawned) THEN !If NoDetstoRotate is greater than 0 on any processor, need to rotate all arrays, otherwise will overwrite each other etc.

            InitNoDetstoRotate=NoDetstoRotate
            !For now, rotate this same sized array each time, even if not completely necessary.
            !We are currently just making the determinant (and its sign) 0 if we no longer want to rotate it, but could in the future remove it from the
            !array.  Probably doesn't make all that much difference because will never find a 0 determinant in the guiding function.

            do n=1,nProcessors-1
            !Rotate the DetstoRotate, SigntoRotate and NoDetstoRotate values.

                DetsEq=.false.
                CombSign=0
                SigntoRotate(0)=InitNoDetstoRotate

                !Send the sign of those we want to rotate to the next processor.
                !Element 0 is the InitNoDetstoRotate value for this processor, send this as well.
                CALL MPI_BSend(SigntoRotate(0:InitNoDetstoRotate),InitNoDetstoRotate+1,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),123,MPI_COMM_WORLD,error)
                IF(error.ne.MPI_SUCCESS) THEN
                    CALL Stop_All("RotoAnnihilGuidingFunc","Error in sending signs")
                ENDIF

                !Then send the determinants
                CALL MPI_BSend(DetstoRotate(0:NIfD,1:InitNoDetstoRotate),(NIfD+1)*InitNoDetstoRotate,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),456,MPI_COMM_WORLD,error)
                IF(error.ne.MPI_SUCCESS) THEN
                    CALL Stop_All("RotoAnnihilGuidingFunc","Error in sending particles")
                ENDIF

                !Receives signs.
                !Receive max possible, will only overwrite those that are actually being sent.
                CALL MPI_Recv(SigntoRotate2(0:iGuideDets),iGuideDets+1,MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),123,MPI_COMM_WORLD,Stat,error)
                IF(error.ne.MPI_SUCCESS) THEN
                    CALL Stop_All("RotoAnnihilGuidingFunc","Error in receiving signs")
                ENDIF

                InitNoDetstoRotate=SigntoRotate2(0)

                !Recieve determinants
                CALL MPI_Recv(DetstoRotate2(0:NIfD,1:InitNoDetstoRotate),InitNoDetstoRotate*(NIfD+1),MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),456,MPI_COMM_WORLD,Stat,error)
                IF(error.ne.MPI_SUCCESS) THEN
                    CALL Stop_All("RotoAnnihilGuidingFunc","Error in receiving particles")
                ENDIF

                do i=1,InitNoDetstoRotate
                    SigntoRotate(i)=SigntoRotate2(i)
                    DetstoRotate(0:NIfD,i)=DetstoRotate2(0:NIfD,i)
                enddo

                !If a determinant has walkers in the guiding function on the same determinant with the same sign, add the spawned (rotate) walkers
                !to the spawned list of that processor (no longer need to rotate).
                !If the walkers have opposite sign, annihilate and rotate any remaining from the spawned list.
                !If there are no walkers on the guiding function determinant, keep rotating.
            
                !Take a rotated determinant and run through to find it in the guiding function of that processor.
                do i=1,InitNoDetstoRotate
                    !If the number of particles being rotated has become zero, don't need to search through the guiding function for particles to annihilate.
                    IF(Signtorotate(i).ne.0) THEN
                        do j=1,iGuideDets
                            DetsEq=.false.
                            DetsEq=DetBitEQ(DetstoRotate(0:NIfD,i),GuideFuncDets(0:NIfD,j),NIfD)

                            IF(DetsEq) THEN
                                CombSign=SigntoRotate(i)*GuideFuncSign(j)
                                !IF this is negative, the guiding function annihilates the spawned particles.
                                IF(CombSign.lt.0) THEN

                                    IF(ABS(SigntoRotate(i)).gt.ABS(GuideFuncSign(j))) THEN
                                        !If there are still too many in the spawned (rotated) array to be all annihilated by the guiding function on that
                                        !processor, annihilate what you can, but leave the determinant in the array to keep rotating.
                                        SigntoRotate(i)=SigntoRotate(i)+GuideFuncSign(j)
                                        ! Add because these are opposite signs.
                                        GuideFuncSign(j)=0
                                        
                                    ELSEIF(ABS(SigntoRotate(i)).eq.ABS(GuideFuncSign(j))) THEN
                                        SigntoRotate(i)=0
                                        GuideFuncSign(j)=0

                                    ELSE
                                        ! The spawned are all annihilated, and the guiding function is decreased by that number.
                                        GuideFuncSign(j)=GuideFuncSign(j)+SigntoRotate(i)
                                        SigntoRotate(i)=0
                                    ENDIF

                                    !IF the combined sign (CombSign) is positive, there are walkers in the guiding function of that processor with the same
                                    !sign. Thus no annihilation occurs and these particles just continue to rotate around (they will just end up back on the
                                    !the original processor where they'll be recombined back into SpawnedPart.

                                    !If CombSign is 0, there are no walkers on the guiding function (for that processor).
                                    !Continue rotating spawned walkers to see if there are any on the next processor to annihilate with.

                                    !If we have found a determinant in the guiding function that matches that in the spawned, can stop searching the guiding 
                                    !function, there will be no more matches.

                                ENDIF

                                EXIT
                            ENDIF

                        enddo
                    ENDIF
                enddo
                    
            enddo

            !If back to original processor and still have walkers in rotating arrays, just add them to the spawned list.
            !Do one final rotation, end up on original processor - add remaining particles to spawned list.

            SigntoRotate(0)=InitNoDetstoRotate

            !Send the sign of those we want to rotate to the next processor.
            CALL MPI_BSend(SigntoRotate(0:InitNoDetstoRotate),InitNoDetstoRotate+1,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),123,MPI_COMM_WORLD,error)
            IF(error.ne.MPI_SUCCESS) THEN
                CALL Stop_All("RotoAnnihilGuidingFunc","Error in sending signs")
            ENDIF

            !Then send the determinants
            CALL MPI_BSend(DetstoRotate(0:NIfD,1:InitNoDetstoRotate),(NIfD+1)*InitNoDetstoRotate,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),456,MPI_COMM_WORLD,error)
            IF(error.ne.MPI_SUCCESS) THEN
                CALL Stop_All("RotoAnnihilGuidingFunc","Error in sending particles")
            ENDIF

            !Receives signs.
            CALL MPI_Recv(SigntoRotate2(0:iGuideDets),iGuideDets+1,MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),123,MPI_COMM_WORLD,Stat,error)
            IF(error.ne.MPI_SUCCESS) THEN
                CALL Stop_All("RotoAnnihilGuidingFunc","Error in receiving signs")
            ENDIF
            
            InitNoDetstoRotate=SigntoRotate2(0)

            !Recieve determinants
            CALL MPI_Recv(DetstoRotate2(0:NIfD,1:InitNoDetstoRotate),InitNoDetstoRotate*(NIfD+1),MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),456,MPI_COMM_WORLD,Stat,error)
            IF(error.ne.MPI_SUCCESS) THEN
                CALL Stop_All("RotoAnnihilGuidingFunc","Error in receiving particles")
            ENDIF

            do i=1,InitNoDetstoRotate
                SigntoRotate(i)=SigntoRotate2(i)
                DetstoRotate(0:NIfD,i)=DetstoRotate2(0:NIfD,i)
            enddo

            !Now add all the remaining DetstoRotate2 and their signs to the SpawnedPart and SpawnedSign lists.
            !Since I just copied the determinants from SpawnedPart to DetstoRotate, any in DetstoRotate will already been in SpawnedPart.
            !Need to just search for it and overwrite its spin value.
            do j=1,InitNoDetstoRotate
                tDetinSpawnList=.false.
                do i=1,ValidSpawned
                    DetsEq=.false.
                    DetsEq=DetBitEQ(SpawnedParts(0:NIfD,i),DetstoRotate(0:NIfD,j),NIfD)
                    IF(DetsEq) THEN
                        SpawnedSign(i)=SigntoRotate(j)
                        tDetinSpawnList=.true.
                        EXIT
                    ENDIF
                enddo
                IF(.not.tDetinSpawnList) THEN
                    WRITE(6,*) 'Determinant from rotate list : ',DetstoRotate(0:NIfD,j)
!                    do i=1,ValidSpawned
!                        WRITE(6,*) SpawnedParts(0:NIfD,i)
!                    enddo
                    CALL FLUSH(6)
                    CALL Stop_All("RotoAnnihilGuidingFunc","Determinant from rotated list cannot be found in SpawnedParts.")
                ENDIF
            enddo

        ENDIF

        !Calculated the number of walkers in the guiding function after this annihilation.
        iInitGuideParts=0
        AlliInitGuideParts=0
        do i=1,iGuideDets
            iInitGuideParts=iInitGuideParts+ABS(GuideFuncSign(i))
        enddo
!        CALL MPI_Reduce(iInitGuideParts,AlliInitGuideParts,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)

        !Need to calculate the contribution to the HF from the guiding function, and then also the contribution from doubles.
        GuideFuncHF=GuideFuncHF+GuideFuncSign(GuideFuncHFIndex)


        !Run through all other determinants in the guiding function.  Find out if they are doubly excited.  Find H elements, and multiply by number on that double.
        do i=1,iGuideDets
            CALL FindBitExcitLevel(GuideFuncDets(0:NIfD,i),iLutHF,NIfD,ExcitLevel,2)
            IF(ExcitLevel.eq.2) THEN
                DoubDet(:)=0
                CALL DecodeBitDet(DoubDet,GuideFuncDets(0:NIfD,i),NEl,NIfD)
                HdoubTemp=GetHElement2(HFDet,DoubDet,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
                HDoub=REAL(HDoubTemp%v,r2)
                GuideFuncDoub=GuideFuncDoub+(GuideFuncSign(i)*Hdoub)
            ENDIF
        enddo


    ENDSUBROUTINE RotoAnnihilGuidingFunc



    SUBROUTINE WriteGuidingFunc()
!This routine writes out the iGuideDets determinants with the largest no of walkers at the end of a calculation.    
!Really it finds the number of walkers on the final determinant in iGuideDets, and writes out any determinants
!with this number or more walkers on them (accounts for degeneracies in the number of walkers on the determinants).
        INTEGER :: i,j,k,MinGuideDetPop,ierr,error,TestSum,OrigiGuideDets
        INTEGER :: AlliGuideDets,CompiGuideDets
        CHARACTER(len=*), PARAMETER :: this_routine='WriteGuidingFunc'
        INTEGER , ALLOCATABLE :: AllCurrentDets(:,:),AllCurrentSign(:)
        INTEGER :: AllCurrentDetsTag,AllCurrentSignTag
        INTEGER :: RecvCounts(nProcessors),Offsets(nProcessors),RecvCounts02(nProcessors),OffSets02(nProcessors)
        LOGICAL :: DetBitEQ

        CALL FLUSH(6)
        IF(iGuideDets.gt.TotWalkers) CALL Stop_All(this_routine,'iGuideDets is greater than the number of populated determinants')

!        WRITE(6,*) 'the determinants and sign, on each processor, before I touched them'
!        do j=1,TotWalkers
!            WRITE(6,*) CurrentDets(0:NIfD,j),CurrentSign(j)
!        enddo
        
! Firstly order CurrentSign in descending absolute value, taking the corresponding CurrentDets with it.
        CALL SortBitSign(TotWalkers,CurrentSign(1:TotWalkers),NIfD,CurrentDets(0:NIfD,1:TotWalkers))

! Then run through CurrentSign, finding out how many walkers are on the iGuideDets most populated, and counting how
! many have this many walkers or more.
    
        OrigiGuideDets=iGuideDets
        MinGuideDetPop=0
        MinGuideDetPop=ABS(CurrentSign(iGuideDets))
        ! This is the minimum number of walkers on the included determinants.
        ! Also want to include determinants with the same number of walkers.
        j=iGuideDets+1
        do while (ABS(CurrentSign(j)).eq.MinGuideDetPop)
            j=j+1
            iGuideDets=iGuideDets+1
        enddo
 
!        WRITE(6,*) 'The most populated determinants on each processor, ordered in terms of sign'
!        WRITE(6,*) iGuideDets,' determinants included in the guiding function.'    
!        do j=1,iGuideDets
!            WRITE(6,*) CurrentDets(0:NIfD,j),CurrentSign(j)
!        enddo
!        CALL FLUSH(6)
        
! Now take the iGuideDets determinants and reorder them back in terms of determinants (taking the sign with them).
        
        CALL SortBitDets(iGuideDets,CurrentDets(0:NIfD,1:iGuideDets),NIfD,CurrentSign(1:iGuideDets))

!        WRITE(6,*) 'The most populated determinants on each processor, ordered by determinant'
!        WRITE(6,*) iGuideDets,' determinants included in the guiding function.'    
!        do j=1,iGuideDets
!            WRITE(6,*) CurrentDets(0:NIfD,j),CurrentSign(j)
!        enddo
        
! Calculate RecvCounts(1:nProcessors), and OffSets for the Gatherv calculation
! Need to gather the iGuideDets values from each processor for this.
        CALL MPI_Gather(iGuideDets,1,MPI_INTEGER,RecvCounts,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        Offsets(:)=0
        do j=1,nProcessors-1
            OffSets(j+1)=RecvCounts(j)+OffSets(j)
        enddo
        do j=1,nProcessors
            RecvCounts02(j)=RecvCounts(j)*(NIfD+1)
        enddo
        OffSets02(:)=0
        do j=1,nProcessors-1
            OffSets02(j+1)=RecvCounts02(j)+OffSets02(j)
        enddo

!        do j=1,nProcessors
!            WRITE(6,*) RecvCounts(j),OffSets(j),RecvCounts02(j),OffSets02(j)
!        enddo


! Find out the total number of iGuideDets on all processors (roughly the input, but not exactly because of 
! of the degeneracies in the populations). 
! This should be the same as the sum of RecvCounts.
        AlliGuideDets=0
        CALL MPI_Reduce(iGuideDets,AlliGuideDets,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
        IF(iProcIndex.eq.Root) THEN
            TestSum=0
            do i=1,nProcessors
                TestSum=TestSum+RecvCounts(i)
            enddo
            IF(AlliGuideDets.ne.TestSum) CALL Stop_All("WriteGuidingFunc","Error in parallel sum")
        ENDIF
        
!        IF(iProcIndex.eq.Root) WRITE(6,*) 'AlliGuideDets ',AlliGuideDets
!        CALL FLUSH(6)


! Make arrays containing the most populated determinants from all processors.        

        IF(iProcIndex.eq.Root) THEN
            ALLOCATE(AllCurrentDets(0:NIfD,1:AlliGuideDets),stat=ierr)
            CALL LogMemAlloc('AllCurrentDets',(NIfD+1)*AlliGuideDets,4,this_routine,AllCurrentDetsTag,ierr)
            ALLOCATE(AllCurrentSign(1:AlliGuideDets),stat=ierr)
            CALL LogMemAlloc('AllCurrentSign',AlliGuideDets,4,this_routine,AllCurrentSignTag,ierr)
        ENDIF

        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        CALL MPI_Gatherv(CurrentSign(1:iGuideDets),iGuideDets,MPI_INTEGER,AllCurrentSign(1:AlliGuideDets),&
        &RecvCounts,Offsets,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL MPI_Gatherv(CurrentDets(0:NIfD,1:iGuideDets),((NIfD+1)*iGuideDets),MPI_INTEGER,&
        &AllCurrentDets(0:NIfD,1:AlliGuideDets),RecvCounts02,Offsets02,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL FLUSH(6)


        IF(iProcIndex.eq.Root) THEN

!            OPEN(31,file='GUIDINGFUNCall-01',status='unknown')
!            WRITE(31,*) 'The determinants from each processor combined'
!            WRITE(31,*) AlliGuideDets,' determinants included in the guiding function.'    
!            do j=1,AlliGuideDets
!                WRITE(31,*) AllCurrentDets(0:NIfD,j),AllCurrentSign(j)
!            enddo
!            CLOSE(31)

! Having taken the largest occupied determinants from each processor, select the iGuideDets most populated from
! this total list.
! I.e. reorder in terms of population, compressing so that each determinant only appears once.  Then select out the top iGuideDets 
! and reorder by determinant once again.

! From the list of determinants from all processors, order in terms of determinant so that this may be compressed to make sure each
! determinant only appears once.
            CALL SortBitDets(AlliGuideDets,AllCurrentDets(0:NIfD,1:AlliGuideDets),NIfD,AllCurrentSign(1:AlliGuideDets))

!            OPEN(32,file='GUIDINGFUNCall-02',status='unknown')
!            WRITE(32,*) 'The determinants from each processor combined'
!            WRITE(32,*) AlliGuideDets,' determinants included in the guiding function.'    
!            do j=1,AlliGuideDets
!                WRITE(32,*) AllCurrentDets(0:NIfD,j),AllCurrentSign(j)
!            enddo
!            CLOSE(32)
            

! Find out if two determinants next to each are the same, if so add their signs and move the ones below them up.
! DetBitEq returns true if two determinants are identical, or false otherwise.
            CompiGuideDets=AlliGuideDets
            do i=1,AlliGuideDets-1

                IF(i.gt.CompiGuideDets) EXIT 
                ! This means we have got to the end of the compressed list, don't want to keep going, as all determinants after this are 0. 

                do while (DetBitEQ(AllCurrentDets(0:NIfD,i),AllCurrentDets(0:NIfD,i+1),NIfD))
                    ! Take a determinant, if the one above it is identical, add its sign to the original and move the others up to overwrite
                    ! the second.
                    ! Repeat this until the i+1 determinant is no longer equal to the i determinant.

                    IF((AllCurrentSign(i)*AllCurrentSign(i+1)).lt.0) THEN
                        WRITE(6,*) 'Determinant populated with opposite signs',AllCurrentDets(0:NIfD,i)
                        WRITE(6,*) 'Identical determinant next to it',AllCurrentDets(0:NIfD,i+1)
                        CALL FLUSH(6)
                        CALL Stop_All("WriteGuidingFunc","Identical determinants populated with opposite sign")
                    ENDIF

                    AllCurrentSign(i)=AllCurrentSign(i)+AllCurrentSign(i+1)

                    CompiGuideDets=CompiGuideDets-1

                    do j=i+2,AlliGuideDets
                        AllCurrentDets(0:NIfD,j-1)=AllCurrentDets(0:NIfD,j)
                        AllCurrentSign(j-1)=AllCurrentSign(j)
                    enddo
                    ! Zero the last determinant
                    AllCurrentSign(AlliGuideDets)=0
                    AllCurrentDets(0:NIfD,AlliGuideDets)=0
                enddo
            enddo
            AlliGuideDets=CompiGuideDets
     
!            OPEN(33,file='GUIDINGFUNCall-03',status='unknown')
!            WRITE(33,*) 'The determinants from each processor combined and compressed, ordered by determinant'
!            WRITE(33,*) AlliGuideDets,' determinants included in the guiding function.'    
!            do j=1,AlliGuideDets
!                WRITE(33,*) AllCurrentDets(0:NIfD,j),AllCurrentSign(j)
!            enddo
!            CLOSE(33)

          
! Reorder the compressed list of determinants by population (i.e. descending according to absolute value of AllCurrentSign, taking determinant with it).        
            CALL SortBitSign(AlliGuideDets,AllCurrentSign(1:AlliGuideDets),NIfD,AllCurrentDets(0:NIfD,1:AlliGuideDets))


! From this total list of most populated determinants, pick out the iGuideDets most populated, along with any with the same number of walkers as the last.        
            iGuideDets=OrigiGuideDets
            MinGuideDetPop=0
            MinGuideDetPop=ABS(AllCurrentSign(iGuideDets))
            ! This is the minimum number of walkers on the included determinants.
            ! Also want to include determinants with the same number of walkers.
            j=iGuideDets+1
            do while (ABS(AllCurrentSign(j)).eq.MinGuideDetPop)
                j=j+1
                iGuideDets=iGuideDets+1
            enddo
      
!            OPEN(34,file='GUIDINGFUNCall-04',status='unknown')
!            WRITE(34,*) 'The determinants from each processor combined and compressed, ordered by sign'
!            WRITE(34,*) iGuideDets,' determinants included in the guiding function.'    
!            do j=1,iGuideDets
!                WRITE(34,*) AllCurrentDets(0:NIfD,j),AllCurrentSign(j)
!            enddo
!            CLOSE(34)


! Now take the iGuideDets determinants and reorder them back in terms of determinants (taking the sign with them).
            CALL SortBitDets(iGuideDets,AllCurrentDets(0:NIfD,1:iGuideDets),NIfD,AllCurrentSign(1:iGuideDets))


! Write the iGuideDets most populated determinants (in order of their bit strings) to a file.
            OPEN(35,file='GUIDINGFUNC',status='unknown')
            WRITE(35,*) iGuideDets,' determinants included in the guiding function.'    
            do j=1,iGuideDets
                WRITE(35,*) AllCurrentDets(0:NIfD,j),AllCurrentSign(j)
            enddo
            CLOSE(35)

            DEALLOCATE(AllCurrentDets)
            CALL LogMemDealloc(this_routine,AllCurrentDetsTag)
            DEALLOCATE(AllCurrentSign)
            CALL LogMemDealloc(this_routine,AllCurrentSignTag)
        
        ENDIF


    END SUBROUTINE WriteGuidingFunc




    SUBROUTINE WriteFinalGuidingFunc()
        INTEGER :: j,AllGuideFuncSignTag,ierr,error
        INTEGER , ALLOCATABLE :: AllGuideFuncSign(:)

        IF(iProcIndex.eq.Root) THEN
            ALLOCATE(AllGuideFuncSign(iGuideDets),stat=ierr)
            CALL LogMemAlloc('AllGuideFuncSign',iGuideDets,4,'WriteFinalGuidingFunc',AllGuideFuncSignTag,ierr)
        ENDIF 

        CALL MPI_Reduce(GuideFuncSign(1:iGuideDets),AllGuideFuncSign,iGuideDets,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
 
        IF(iProcIndex.eq.Root) THEN
            OPEN(38,file='GUIDINGFUNCfinal',status='unknown')
            WRITE(38,*) iGuideDets,' determinants included in the guiding function.'    
            WRITE(38,'(A11,A28,A11,A15)') "Determinant","InitialPop","FinalPop","TotalChange"
            do j=1,iGuideDets
                WRITE(38,*) GuideFuncDets(0:NIfD,j),InitGuideFuncSign(j),AllGuideFuncSign(j),(ABS(AllGuideFuncSign(j))-ABS(InitGuideFuncSign(j)))
            enddo
            CLOSE(38)

            DEALLOCATE(AllGuideFuncSign)
            CALL LogMemDealloc('WriteFinalGuidingFunc',AllGuideFuncSignTag)
            DEALLOCATE(InitGuideFuncSign)
            CALL LogMemDealloc('WriteFinalGuidingFunc',InitGuideFuncSignTag)
 
        ENDIF


    ENDSUBROUTINE WriteFinalGuidingFunc



    SUBROUTINE InitSpawnDominant()
!It then scales the number of walkers on each determinant up so that the total is that specified in the input for iInitGuideParts. 
!The result is an array of determinats and a corresponding array of populations (with sign) for the guiding function.
        INTEGER :: i,j,ierr,error,iNoExcLevels,Temp,iMinDomLevPop,SumExcLevPop,ExcitLevel
        CHARACTER(len=*), PARAMETER :: this_routine='InitSpawnDominant'

        IF(iProcIndex.eq.Root) THEN
            OPEN(41,FILE='DOMINANTDETS',Status='old')
            READ(41,*) iNoDomDets
            READ(41,*) iNoExcLevels 
            READ(41,*) iMinDomLev,SumExcLevPop
            iMaxDomLev=iMinDomLev+iNoExcLevels-1
!            IF(iMinDomLev.le.2) CALL Stop_All('InitSpawnDominant','Code is not set up to deal with removing doubles or lower.')
        ENDIF
        
        CALL MPI_Bcast(iMinDomLev,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL MPI_Bcast(iMaxDomLev,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL MPI_Bcast(iNoDomDets,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL MPI_Bcast(iNoExcLevels,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)

        ALLOCATE(DomExcIndex(iMinDomLev:(iMinDomLev+iNoExcLevels)),stat=ierr)
        CALL LogMemAlloc('DomExcIndex',iNoExcLevels+1,4,this_routine,DomExcIndexTag,ierr)
        DomExcIndex(:)=0
        
        IF(iProcIndex.eq.Root) THEN
            DomExcIndex(iMinDomLev)=1
            i=iMinDomLev+1
            do while(i.le.iMaxDomLev)
                DomExcIndex(i)=SumExcLevPop+DomExcIndex(i-1)
                READ(41,*) ExcitLevel,Temp
                SumExcLevPop=Temp
                i=i+1
            enddo
            DomExcIndex(iMaxDomLev+1)=SumExcLevPop+DomExcIndex(i-1)
            IF(DomExcIndex(iMaxDomLev+1).ne.(iNoDomDets+1)) CALL Stop_All(this_routine,'ERROR in filling the indexing excitation indexing array.')
        ENDIF


        CALL MPI_Bcast(DomExcIndex(iMinDomLev:(iMaxDomLev+1)),iMaxDomLev-iMinDomLev+2,MPI_INTEGER,Root,MPI_COMM_WORLD,error)


        ALLOCATE(DomDets(0:NIfD,1:iNoDomDets),stat=ierr)
        CALL LogMemAlloc('DomDets',(NIfD+1)*iNoDomDets,4,this_routine,DomDetsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory for dominant determinants')

        IF(iProcIndex.eq.Root) THEN
            !Set up the determinant and sign arrays by reading in from the GUIDINGFUNC file.
            j=1
            do while (j.le.iNoDomDets)
                READ(41,*) DomDets(0:NIfD,j)
                j=j+1
            enddo
            CLOSE(41)
        ENDIF

        ! Broadcast this list of DomDets to all processors
        CALL MPI_Bcast(DomDets(0:NIfD,1:iNoDomDets),(NIfD+1)*iNoDomDets,MPI_INTEGER,Root,MPI_COMM_WORLD,ierr)


        IF(tMinorDetsStar) THEN
            IF(.not.tRotoAnnihil) THEN
                CALL Stop_All(this_routine,'STARMINORDETERMINANTS can only be used with rotoannihilation.') 
            ELSEIF(tFixShiftShell) THEN
                CALL Stop_All(this_routine,'STARMINORDETERMINANTS cannot be used with the fixed shift shell approximation.') 
            ELSEIF(iMinDomLev.le.2) THEN
                CALL Stop_All(this_routine,'STARMINORDETERMINANTS is not set up to deal with removing doubles or lower.')
            ELSE
                CALL InitMinorDetsStar()
            ENDIF
        ENDIF


    ENDSUBROUTINE InitSpawnDominant 




    SUBROUTINE PrintDominantDets()
! This routine takes the list of determinants with particles on them, and picks out those with excitation levels between the max and min
! specified in the input file.  It then orders these in terms of population, takes the iNoDominantDets most populated and prints them 
! in order of excitation level and then determinant to a file named DOMINANTDETS.
        INTEGER :: i,j,k,ierr,error,TestSum,ExcitLevel,NoExcDets,AllNoExcDets,ExcDetsTag,ExcSignTag,AllExcDetsTag,AllExcSignTag,ipos
        INTEGER :: ExcDetsIndex,MinDomDetPop,AllExcLevelTag,ExcLevelTag,CurrExcitLevel,NoExcitLevel,OrigiDominantDets,HFPop,CurriNoDominantDets
        CHARACTER(len=*), PARAMETER :: this_routine='PrintDominantDets'
        REAL*8 :: MinRelDomPop,SpinTot,NormDef
        INTEGER , ALLOCATABLE :: ExcDets(:,:),ExcSign(:),AllExcDets(:,:),AllExcSign(:)
        INTEGER :: RecvCounts(nProcessors),Offsets(nProcessors),RecvCounts02(nProcessors),OffSets02(nProcessors),DetCurr(NEl),NormDefTrunc,NormDefTot
        INTEGER :: SpinCoupDetBit(0:NIfD),SpinCoupDet(NEl),OpenShell(2,NEl),UpSpin(NEl),NoOpenShell,NoUpSpin,iRead,PartInd,ID1,ID2,iComb,TempSign
        LOGICAL :: tDoubOcc,tSuccess

        CALL FLUSH(6)
        IF(.not.tRotoAnnihil) CALL Stop_All(this_routine,'PRINTDOMINANTDETS can only be used with rotoannihilation.')

        WRITE(6,'(A13,I10,A53,I3,A5,I3)') 'Printing the ',iNoDominantDets,' dominant determinants with excitation level between ',MinExcDom,' and ',MaxExcDom


!        WRITE(6,*) 'the determinants and sign, on each processor, before I touched them'
!        do j=1,TotWalkers
!            WRITE(6,*) CurrentDets(0:NIfD,j),CurrentSign(j)
!        enddo


! Firstly, copy the determinants with excitation level between the min and max into a separate array (and do the same with sign).  
! Need to first count the number that are going to go in this array.
        NoExcDets=0
        AllNoExcDets=0
        do i=1,TotWalkers
            CALL FindBitExcitLevel(CurrentDets(0:NIfD,i),iLutHF(0:NIfD),NIfD,ExcitLevel,MaxExcDom)
            IF((ExcitLevel.ge.MinExcDom).and.(ExcitLevel.le.MaxExcDom)) THEN
                NoExcDets=NoExcDets+1
            ELSEIF(ExcitLevel.eq.0) THEN
                HFPop=ABS(CurrentSign(i))
            ENDIF
        enddo

        CALL MPI_Reduce(NoExcDets,AllNoExcDets,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)


! Allocate arrays of this size - these are the ones that will be reordered to find the iNoDominantDets most populated etc.        
        ALLOCATE(ExcDets(0:NIfD,1:NoExcDets),stat=ierr)
        CALL LogMemAlloc('ExcDets',(NIfD+1)*NoExcDets,4,this_routine,ExcDetsTag,ierr)
        ALLOCATE(ExcSign(1:NoExcDets),stat=ierr)
        CALL LogMemAlloc('ExcSign',NoExcDets,4,this_routine,ExcSignTag,ierr)
 
        IF(iProcIndex.eq.Root) THEN
            ALLOCATE(AllExcDets(0:NIfD,1:(10*AllNoExcDets)),stat=ierr)
            CALL LogMemAlloc('AllExcDets',(NIfD+1)*10*AllNoExcDets,4,this_routine,AllExcDetsTag,ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'ERROR allocating memory to AllExcDets.')
            ALLOCATE(AllExcSign(1:(10*AllNoExcDets)),stat=ierr)
            CALL LogMemAlloc('AllExcSign',10*AllNoExcDets,4,this_routine,AllExcSignTag,ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'ERROR allocating memory to AllExcSign.')
        ENDIF


! Now run through the occupied determinants.  If the determinant has the correct excitation level, add it to the ExcDets array, and
! add its sign to the ExcSign array.
        ExcDetsIndex=0
        do i=1,TotWalkers
            CALL FindBitExcitLevel(CurrentDets(0:NIfD,i),iLutHF(0:NIfD),NIfD,ExcitLevel,MaxExcDom)
            IF((ExcitLevel.ge.MinExcDom).and.(ExcitLevel.le.MaxExcDom)) THEN
                ExcDetsIndex=ExcDetsIndex+1
                ExcDets(0:NIfD,ExcDetsIndex)=CurrentDets(0:NIfD,i)
                ExcSign(ExcDetsIndex)=CurrentSign(i)
            ENDIF
        enddo


! Now need to collect this array on the root processor and find the most populated determinants.
! First set up the RecvCounts and OffSet arrays.

! Need to gather the NoExcDets values from each processor for this.
        CALL MPI_Gather(NoExcDets,1,MPI_INTEGER,RecvCounts,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        
        Offsets(:)=0
        do j=1,nProcessors-1
            OffSets(j+1)=RecvCounts(j)+OffSets(j)
        enddo
        do j=1,nProcessors
            RecvCounts02(j)=RecvCounts(j)*(NIfD+1)
        enddo
        OffSets02(:)=0
        do j=1,nProcessors-1
            OffSets02(j+1)=RecvCounts02(j)+OffSets02(j)
        enddo

!        do j=1,nProcessors
!            WRITE(6,*) RecvCounts(j),OffSets(j),RecvCounts02(j),OffSets02(j)
!        enddo

        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        CALL MPI_Gatherv(ExcSign(1:NoExcDets),NoExcDets,MPI_INTEGER,AllExcSign(1:AllNoExcDets),&
        &RecvCounts,Offsets,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL MPI_Gatherv(ExcDets(0:NIfD,1:NoExcDets),((NIfD+1)*NoExcDets),MPI_INTEGER,&
        &AllExcDets(0:NIfD,1:AllNoExcDets),RecvCounts02,Offsets02,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL FLUSH(6)
       
! Now that we have arrays on the root processor with the determinants and sign of correct excitation level, need to order these 
! in descending absolute value, taking the corresponding sign with it.
        IF(iProcIndex.eq.Root) THEN
 
            IF(iNoDominantDets.gt.AllNoExcDets) THEN
                WRITE(6,*) 'iNoDominantDets: ',iNoDominantDets
                WRITE(6,*) 'AllNoExcDets: ',AllNoExcDets
                CALL Stop_All(this_routine,'Not enough determinants are occupied to pick out the number of dominant requested.')
            ENDIF


            CALL SortBitSign(AllNoExcDets,AllExcSign(1:AllNoExcDets),NIfD,AllExcDets(0:NIfD,1:AllNoExcDets))

! Then run through AllExcSign, finding out how many walkers are on the iNoDominantDets most populated, and counting how
! many have this many walkers or more.
        
            OrigiDominantDets=iNoDominantDets
            MinDomDetPop=0
            MinDomDetPop=ABS(AllExcSign(iNoDominantDets))
            ! This is the minimum number of walkers on the included determinants.
            ! Also want to include determinants with the same number of walkers.
            j=iNoDominantDets+1
            do while (ABS(AllExcSign(j)).eq.MinDomDetPop)
                j=j+1
                iNoDominantDets=iNoDominantDets+1
            enddo
 
            OPEN(39,file='DOMINANTDETSdescpop',status='unknown')
            WRITE(39,*) AllNoExcDets,' determinants with the right excitation level.'    
            WRITE(39,*) iNoDominantDets,' with population ',MinDomDetPop, ' and above.'
            do j=1,AllNoExcDets
                WRITE(39,*) AllExcDets(0:NIfD,j),AllExcSign(j)
            enddo
            CLOSE(39)
          
            WRITE(6,*) 'This amounts to determinants with population ',MinDomDetPop,' and larger.'
!            WRITE(6,*) 'HFPop',HFPop
!            WRITE(6,*) 'MinDomDetPop',MinDomDetPop
            MinRelDomPop=REAL(MinDomDetPop)/REAL(HFPop)
            WRITE(6,*) 'These determinants have amplitude ; ',MinRelDomPop,' relative to the most populated determinant.' 

! In order to do binary searches for the spin determinants, need to sort the determinants back into order.
! Do this in two separate lots, 1:iNoDominantDets and iNoDominantDets+1:AllNoExcDets

            CALL SortBitDets(iNoDominantDets,AllExcDets(0:NIfD,1:iNoDominantDets),NIfD,AllExcSign(1:iNoDominantDets))
            CALL SortBitDets((AllNoExcDets-iNoDominantDets),AllExcDets(0:NIfD,(iNoDominantDets+1):AllNoExcDets),NIfD,&
                            &AllExcSign((iNoDominantDets+1):AllNoExcDets))
 
            OPEN(47,file='DOMINANTDETSsorted',status='unknown')
            WRITE(47,*) AllNoExcDets,' determinants with the right excitation level.'    
            WRITE(47,*) iNoDominantDets,' with population ',MinDomDetPop, ' and above.'
            do j=1,AllNoExcDets
                WRITE(47,*) AllExcDets(0:NIfD,j),AllExcSign(j)
            enddo
            CLOSE(47)
 

! Run through all the determinants, finding sum_k=1->AllNoExcDets c_k^2 = NormDefTot.  Want to do this before doing the spin
! coupled stuff when some determinants have overwritten others.
            NormDefTot=0
            do i=1,AllNoExcDets
                NormDefTot=NormDefTot+(AllExcSign(i)**2)
            enddo



! Put a bit in here to run through the determinants, construct all the different spin eigenstates from each and make sure 
! they are all included.  - need to add them to the list maintaining order.
! At the moment, these are ordered by determinant, and only the first iNoDominantDets are going to be taken.
! Determinants that are added to the list should be added in order.
            IF(.not.tNoDomSpinCoup) THEN 
                WRITE(6,*) 'Also including determinants that are spin coupled to those in the list of dominant dets.'
                CurriNoDominantDets=iNoDominantDets
                ! The number of DominantDets before including spin coupled ones.
                do j=1,CurriNoDominantDets
                    ! Add the sign from this determinant to the Norm Deficiency calc - this will be in trunc.
                    ! Decode the current determinant
                    CALL DecodeBitDet(DetCurr,AllExcDets(0:NIfD,j),NEl,NIfD)
!                    WRITE(6,*) 'DetCurrBit',AllExcDets(:,j)
!                    WRITE(6,*) 'DetCurr',DetCurr(:)

                    ! Need to find overall spin of DetCurr
                    SpinTot=0.D0
                    do i=1,NEl
                        SpinTot=SpinTot+G1(DetCurr(i))%Ms
                    enddo
                    IF(SpinTot.ne.(LMS/2)) THEN
                        WRITE(6,*) 'SpinTot',SpinTot
                        WRITE(6,*) 'LMS',LMS
                        CALL FLUSH(6)
                        CALL Stop_All(this_routine,'Error, the summed Ms values of each electon does not equal LMS/2.')
                    ENDIF

                    ! Run through orbitals - creating a list of open shell electrons
                    OpenShell(:,:)=0
                    NoOpenShell=0
                    UpSpin(:)=0
                    NoUpSpin=0
                    SpinCoupDet(:)=0
                    do i=1,NEl
                        ! For an electron of the current det, find the spatial orbital.
!                        CALL GTID(nBasisMax,DetCurr(i),ID1)
                        ID1=CEILING(REAL(DetCurr(i))/2.0)
                        ! Now run through all other orbitals finding out if the same spat orb is occupied.
                        do k=1,NEl
                            IF(k.eq.i) CYCLE 
                            tDoubOcc=.false.
!                            CALL GTID(nBasisMax,DetCurr(k),ID2)
                            ID2=CEILING(REAL(DetCurr(k))/2.0)
                            IF(ID2.eq.ID1) THEN
                                ! doubly occupied orbital, these will be the same in the spin coupled det.
                                SpinCoupDet(i)=DetCurr(i)
                                SpinCoupDet(k)=DetCurr(k)
                                tDoubOcc=.true.
                                EXIT
                            ENDIF
                        enddo
                        IF(.not.tDoubOcc) THEN
                            ! Singly occupied orbitals, put these in the openshell array.
                            ! OpenShell(1,:) - the spatial orbital 
                            ! OpenShell(2,:) - the electron number
                            NoOpenShell=NoOpenShell+1
                            OpenShell(1,NoOpenShell)=ID1
                            OpenShell(2,NoOpenShell)=i
                            IF(mod(DetCurr(i),2).eq.0) NoUpSpin=NoUpSpin+1
!                            IF(G1(DetCurr(i))%Ms.gt.0) NoUpSpin=NoUpSpin+1
                        ENDIF
                    enddo
!                    WRITE(6,*) 'SpinCoup with just doub occ',SpinCoupDet(:)
!                    WRITE(6,*) 'OpenShell(1,:) - spat orbs',OpenShell(1,:)
!                    WRITE(6,*) 'OpenShell(2,:) - electron no.',OpenShell(2,:)

                    CALL gennct(NoOpenShell,NoUpSpin,iComb)
                    ! in the future, change this so that it doesn't print then read, but takes array straight away.

                    OPEN(91,FILE='COMBINATIONS',Status='old')

                    iRead=1
                    do while(iRead.le.iComb) 
                        READ(91,*) UpSpin(1:NoUpSpin)
!                        WRITE(6,*) 'NoUpSpin',NoUpSpin
!                        WRITE(6,*) 'UpSpin',UpSpin(1:NoUpSpin)

                        do i=1,NoOpenShell
                            ! Make all the open shell electrons beta, then overwrite the alpha spin.
                            SpinCoupDet(OpenShell(2,i))=(OpenShell(1,i)*2)-1                    
                        enddo

                        do i=1,NoUpSpin
                            SpinCoupDet(OpenShell(2,(UpSpin(i)+1)))=(OpenShell(1,(UpSpin(i)+1)))*2
                        enddo

                        ! Then have to turn SpinCoupDet into the bit string for binary searching...
!                        WRITE(6,*) 'SpinCoup with beta and alpha',SpinCoupDet(:)

                        CALL EncodeBitDet(SpinCoupDet(1:NEl),SpinCoupDetBit(0:NIfD),NEl,NIfD)
!                        WRITE(6,*) 'SpinCoupDetBit',SpinCoupDetBit(:)

                        ! First search through the list of dominant determinants.
                        tSuccess=.false.
                        CALL BinSearchDomParts(AllExcDets(0:NIfD,1:CurriNoDominantDets),SpinCoupDetBit(0:NIfD),1,CurriNoDominantDets,PartInd,tSuccess)
                        IF(tSuccess) THEN
                            ! Determinant found in the dominant list.
                            iRead=iRead+1
                        ELSE
                            ! If not found in the original dominant determinant list then need to search through the determinants that have been added.
                            IF((iNoDominantDets-CurriNoDominantDets).gt.0) THEN
                                CALL BinSearchDomParts(AllExcDets(0:NIfD,(CurriNoDominantDets+1):iNoDominantDets),SpinCoupDetBit(0:NIfD),(CurriNoDominantDets+1),&
                                                    &iNoDominantDets,PartInd,tSuccess)
                                IF((PartInd.le.CurriNoDominantDets).or.(PartInd.gt.iNoDominantDets)) CALL Stop_All(this_routine, '')                                                    
                            ENDIF

                            IF(tSuccess) THEN
                                iRead=iRead+1
                                ! If the determinant has already been added, don't need to do anything.
                            ELSE

                                ! If there are still determinants left in the rest of the list, search through these to check if the determinant is there (to get the sign).
                                IF((AllNoExcDets-iNoDominantDets).gt.0) THEN
                                    CALL BinSearchDomParts(AllExcDets(0:NIfD,(iNoDominantDets+1):AllNoExcDets),SpinCoupDetBit(0:NIfD),(iNoDominantDets+1),&
                                                    &AllNoExcDets,PartInd,tSuccess)
                                ENDIF
                                IF(tSuccess) THEN
                                    ! Determinant found in the rest of the list, put it in dom dets along with its sign.
                                    ! Need to insert the determinant in the correct position, so that these added determinants stay ordered.
                                    ! SearchGen will give back ipos so that the determinant we are inserting is < ipos-1 and ge ipos.
                                    ! I.e this determinant goes in ipos and everything else is moved up 1.
                                    IF((iNoDominantDets-CurriNoDominantDets).gt.0) THEN
                                        CALL SearchGen((iNoDominantDets-CurriNoDominantDets),AllExcDets(0:NIfD,(CurriNoDominantDets+1):iNoDominantDets),&
                                                    &SpinCoupDetBit(0:NIfD),ipos,NIfD)
                                    ELSE
                                        ipos=0
                                    ENDIF
                                    ipos=ipos+CurriNoDominantDets
!                                    CALL SearchGen(iNoDominantDets,AllExcDets(0:NIfD,1:iNoDominantDets),SpinCoupDetBit(0:NIfD),ipos,NIfD)
                                    ! The position the determinant should go is ipos, if the current determinant in ipos is not equal to the SpinCoupDetBit,
                                    ! then we want to insert this determinant here, move all the others up by one, put the determinant that is getting written
                                    ! over at iNoDominantDets+1 in the position of the spincoupdetbit.
                                    
                                    TempSign=AllExcSign(iNoDominantDets+1) 
                                    AllExcDets(0:NIfD,PartInd)=AllExcDets(0:NIfD,iNoDominantDets+1)

                                    do i=iNoDominantDets,ipos,-1                                                    
                                        AllExcDets(0:NIfD,i+1)=AllExcDets(0:NIfD,i)                                                    
                                        AllExcSign(i+1)=AllExcSign(i)
                                    enddo
                                    AllExcDets(0:NIfD,ipos)=SpinCoupDetBit(0:NIfD)
                                    AllExcSign(ipos)=AllExcSign(PartInd)
                                    iRead=iRead+1
                                    iNoDominantDets=iNoDominantDets+1
                                    AllExcSign(PartInd)=TempSign
                                ELSE
                                    ! Determinant not in the list at all, add to added determinants (maintaining order) with a sign of 0.
                                    CALL SearchGen((iNoDominantDets-CurriNoDominantDets),AllExcDets(0:NIfD,(CurriNoDominantDets+1):iNoDominantDets),&
                                                    &SpinCoupDetBit(0:NIfD),ipos,NIfD)
                                    ipos=ipos+CurriNoDominantDets
                                    
                                    AllExcDets(0:NIfD,AllNoExcDets+1)=AllExcDets(0:NIfD,iNoDominantDets+1)
                                    AllExcSign(AllNoExcDets+1)=AllExcSign(iNoDominantDets+1)

                                    do i=iNoDominantDets,ipos,-1                                                    
                                        AllExcDets(0:NIfD,i+1)=AllExcDets(0:NIfD,i)                                                    
                                        AllExcSign(i+1)=AllExcSign(i)
                                    enddo

                                    AllExcDets(0:NIfD,ipos)=SpinCoupDetBit(0:NIfD)
                                    AllExcSign(ipos)=0
                                    iRead=iRead+1
                                    iNoDominantDets=iNoDominantDets+1
                                    AllNoExcDets=AllNoExcDets+1
                                ENDIF
                            ENDIF
                        ENDIF
                        IF(iNoDominantDets.gt.(10*AllNoExcDets)) THEN
                            do i=1,iNoDominantDets
                                WRITE(6,*) AllExcDets(:,i),AllExcSign(i)
                            enddo
                            CALL FLUSH(6)
                            CALL Stop_All(this_routine,'No spin coupled dets has reached larger than AllExcDets array can handle.')
                        ENDIF
                    enddo
                    CLOSE(91)
                    ! do that little bit for all spin coupled determinants and then for all determinants currently in the list.       
!                    CALL FLUSH(6)
!                    CALL Stop_All(this_routine,'cos I meant to.')
                enddo
            ENDIF
            WRITE(6,*) 'By including spin coupling ',iNoDominantDets-CurriNoDominantDets,' more determinants are included.'

! Run through all the determinants we are including in the truncated list (i.e. all those to go in the dominant dets list).            
! Count up the squares of the populations to calculate the normalisation deficiency.
            NormDefTrunc=0
            do i=1,iNoDominantDets
                NormDefTrunc=NormDefTrunc+(AllExcSign(i)**2)
            enddo
            NormDef=1.D0-(REAL(NormDefTrunc)/REAL(NormDefTot))
            WRITE(6,*) 'The NORMALISATION DEFICIENCY of such a truncation is : ',NormDef

 
! Now take the iNoDominantDets determinants and reorder them in terms of excitation level and then determinants.
! SortExcitBitDets orders the determinants first by excitation level and then by determinant, taking the corresponding sign with them.        
! We have just fed the first iNoDominantDets from the list ordered in terms of population.
! Only this amount will be reordered by excitation level then determinant, and then we print out this many only.

            CALL SortExcitBitDets(iNoDominantDets,AllExcDets(0:NIfD,1:iNoDominantDets),NIfD,AllExcSign(1:iNoDominantDets),iLutHF(0:NIfD),NEl)
 
!            OPEN(40,file='DOMINANTDETSexclevelbit',status='unknown')
!            WRITE(40,*) AllNoExcDets,' determinants with the right excitation level.'    
!            do j=1,AllNoExcDets
!                CALL FindBitExcitLevel(AllExcDets(0:NIfD,j),iLutHF(0:NIfD),NIfD,ExcitLevel,MaxExcDom)
!                WRITE(40,*) AllExcDets(0:NIfD,j),AllExcSign(j),ExcitLevel
!            enddo
!            CLOSE(40)


! Write the iGuideDets most populated determinants (in order of their bit strings) to a file.
            OPEN(38,file='DOMINANTDETS',status='unknown')
            WRITE(38,*) iNoDominantDets,' dominant determinants printed here.'    
            WRITE(38,*) (MaxExcDom-MinExcDom+1),' excitation levels included.'
            
! Need to run through these ordered determinants, counting the number in each excitation level.            
            j=1
            do while (j.le.iNoDominantDets)
                NoExcitLevel=0
                i=j 
                IF(ExcitLevel.ne.(CurrExcitLevel+1)) THEN
                    do k=CurrExcitLevel+1,ExcitLevel-1
                        WRITE(38,*) k,0
                    enddo
                ENDIF
                CALL FindBitExcitLevel(AllExcDets(0:NIfD,i),iLutHF(0:NIfD),NIfD,CurrExcitLevel,MaxExcDom)
                ExcitLevel=CurrExcitLevel
                do while (ExcitLevel.eq.CurrExcitLevel)
                    NoExcitLevel=NoExcitLevel+1
                    j=j+1
                    IF(j.gt.iNoDominantDets) EXIT
                    CALL FindBitExcitLevel(AllExcDets(0:NIfD,j),iLutHF(0:NIfD),NIfD,ExcitLevel,MaxExcDom)
                enddo
                WRITE(38,*) CurrExcitLevel,NoExcitLevel
            enddo
            IF(CurrExcitLevel.ne.MaxExcDom) THEN
                do i=CurrExcitLevel+1,MaxExcDom
                    WRITE(38,*) i,0
                enddo
            ENDIF

            do j=1,iNoDominantDets
                WRITE(38,*) AllExcDets(0:NIfD,j),AllExcSign(j)
            enddo
            CLOSE(38)

            DEALLOCATE(AllExcDets)
            CALL LogMemDealloc(this_routine,AllExcDetsTag)
            DEALLOCATE(AllExcSign)
            CALL LogMemDealloc(this_routine,AllExcSignTag)

        ENDIF

        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        DEALLOCATE(ExcDets)
        CALL LogMemDealloc(this_routine,ExcDetsTag)
        DEALLOCATE(ExcSign)
        CALL LogMemDealloc(this_routine,ExcSignTag)

    END SUBROUTINE PrintDominantDets 



    subroutine gennct(n,t,icomb)
       integer n,t
       integer c(t+2)
       integer j
       integer icomb
       icomb=0
       open(90,file='COMBINATIONS',STATUS='UNKNOWN')
!      WRITE(6,*) ' Writing combinations to file COMBINATIONS'
       do j=1,t
           c(j)=j-1
       enddo
       c(t+1)=n
       c(t+2)=0
20     continue
!      visit c(t)
       
       icomb=icomb+1
!      write(90,'(20i3)' ) (c(j)+1,j=1,t)
       
       write(90,'(20i3)' ) (c(j),j=1,t)
       
       do j=1,n
           if(c(j+1).ne.(c(j)+1)) goto 30
           c(j)=j-1
       enddo
30     continue
       if(j.gt.t) then 
!           write(6,*) ' Generated combinations:',ICOMB
           CLOSE(90)
           RETURN
       endif
       c(j)=c(j)+1
       goto 20
       
   end subroutine gennct
       


   SUBROUTINE InitMinorDetsStar()
! This routine simply sets up the arrays etc for the particles spawned into the determinant space that is not in the allowed list.
! MinorStarDets etc are the particles remaining after annihilation, death etc, whereas MinorSpawnDets etc are the walkers newly
! spawned in a particular iteration.
        INTEGER :: ierr
        CHARACTER(len=*), PARAMETER :: this_routine='InitMinorDetsStar'

        ! The actual determinants.
        ALLOCATE(MinorStarDets(0:NIfD,1:MaxWalkersPart),stat=ierr)
        CALL LogMemAlloc('MinorStarDets',(NIfD+1)*MaxWalkersPart,4,this_routine,MinorStarDetsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorStarDets')
        ALLOCATE(MinorSpawnDets(0:NIfD,1:MaxSpawned),stat=ierr)
        CALL LogMemAlloc('MinorSpawnDets',(NIfD+1)*MaxSpawned,4,this_routine,MinorSpawnDetsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorSpawnDets')
        MinorSpawnDets(:,:)=0
        ALLOCATE(MinorSpawnDets2(0:NIfD,1:MaxSpawned),stat=ierr)
        CALL LogMemAlloc('MinorSpawnDets2',(NIfD+1)*MaxSpawned,4,this_routine,MinorSpawnDets2Tag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorSpawnDets2')


        ! The sign (number of walkers) on each determinant.
        ALLOCATE(MinorStarSign(1:MaxWalkersPart),stat=ierr)
        CALL LogMemAlloc('MinorStarSign',MaxWalkersPart,4,this_routine,MinorStarSignTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorStarSign')
        ALLOCATE(MinorSpawnSign(0:MaxSpawned),stat=ierr)
        CALL LogMemAlloc('MinorSpawnSign',MaxSpawned+1,4,this_routine,MinorSpawnSignTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorSpawnSign')
        ALLOCATE(MinorSpawnSign2(0:MaxSpawned),stat=ierr)
        CALL LogMemAlloc('MinorSpawnSign2',MaxSpawned+1,4,this_routine,MinorSpawnSign2Tag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorSpawnSign2')


        ! The parent from which the walker was spawned.
        ALLOCATE(MinorStarParent(0:NIfD,1:MaxWalkersPart),stat=ierr)
        CALL LogMemAlloc('MinorStarParent',(NIfD+1)*MaxWalkersPart,4,this_routine,MinorStarParentTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorStarParent')
        ALLOCATE(MinorSpawnParent(0:NIfD,1:MaxSpawned),stat=ierr)
        CALL LogMemAlloc('MinorSpawnParent',(NIfD+1)*MaxSpawned,4,this_routine,MinorSpawnParentTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorSpawnParent')
        ALLOCATE(MinorSpawnParent2(0:NIfD,1:MaxSpawned),stat=ierr)
        CALL LogMemAlloc('MinorSpawnParent2',(NIfD+1)*MaxSpawned,4,this_routine,MinorSpawnParent2Tag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorSpawnParent2')


        ! The energy of the "forbidden" determinant with a walker on it.
        ALLOCATE(MinorStarHii(1:MaxWalkersPart),stat=ierr)
        CALL LogMemAlloc('MinorStarHii',MaxWalkersPart,8,this_routine,MinorStarHiiTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorStarHii')

        ! The diagonal element connecting the "forbidden" determinant to its allowed parent.
        ALLOCATE(MinorStarHij(1:MaxWalkersPart),stat=ierr)
        CALL LogMemAlloc('MinorStarHij',MaxWalkersPart,8,this_routine,MinorStarHijTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorStarHij')

        ALLOCATE(HashArray(MaxSpawned),stat=ierr)
        CALL LogMemAlloc('HashArray',MaxSpawned,8,this_routine,HashArrayTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to HashArry')
        HashArray(:)=0
        ALLOCATE(Hash2Array(MaxSpawned),stat=ierr)
        CALL LogMemAlloc('Hash2Array',MaxSpawned,8,this_routine,Hash2ArrayTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to Hash2Arry')
        Hash2Array(:)=0
 
        ALLOCATE(IndexTable(MaxSpawned),stat=ierr)
        CALL LogMemAlloc('IndexTable',MaxSpawned,4,this_routine,IndexTableTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to IndexTable')
        IndexTable(:)=0
        ALLOCATE(Index2Table(MaxSpawned),stat=ierr)
        CALL LogMemAlloc('Index2Table',MaxSpawned,4,this_routine,Index2TableTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to Index2Table')
        Index2Table(:)=0

        ALLOCATE(ProcessVec(MaxSpawned),stat=ierr)
        CALL LogMemAlloc('ProcessVec',MaxSpawned,4,this_routine,ProcessVecTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to ProcessVec')
        ProcessVec(:)=0
        ALLOCATE(Process2Vec(MaxSpawned),stat=ierr)
        CALL LogMemAlloc('Process2Vec',MaxSpawned,4,this_routine,Process2VecTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to Process2Vec')
        Process2Vec(:)=0


    END SUBROUTINE InitMinorDetsStar

!Similar to WriteHistogram, but will only print out in order of maximum component, and only the averaged wavefunction
    SUBROUTINE PrintFCIMCPsi()
        use DetCalc , only : FCIDets
        INTEGER :: error,i,nI(NEl),ExcitLevel,j
        REAL*8 :: norm,norm1

        CALL MPI_AllReduce(Histogram,AllHistogram,Det,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
        norm1=0.D0
        do i=1,Det
            norm1=norm1+AllHistogram(i)**2
        enddo
        norm1=SQRT(norm1)
        WRITE(6,*) "Total FCIMC Wavefuction normalisation:",norm1
        do i=1,Det
            AllHistogram(i)=AllHistogram(i)/norm1
        enddo

        IF(tPrintFCIMCPsi) THEN
!Order and print wavefunction

            
            IF(iProcIndex.eq.0) THEN

!We now want to order AllHistogram, taking the corresponding element(s) of FCIDets with it...
                CALL SortRNI(Det,AllHistogram,FCIDets,NIfD+1)
                
                OPEN(17,FILE='FCIMCPsi',STATUS='UNKNOWN')

                norm=0.D0
                do i=1,Det
                    norm=norm+AllHistogram(i)**2
!write out FCIMC Component weight (normalised), current normalisation, excitation level
                    CALL FindBitExcitLevel(iLutHF,FCIDets(0:NIfD,i),NIfD,ExcitLevel,NEl)
                    CALL DecodeBitDet(nI,FCIDets(0:NIfD,i),NEl,NIfD)
                    WRITE(17,"(I13,G25.16,I6,G20.10)",advance='no') i,AllHistogram(i),ExcitLevel,norm
                    do j=1,NEl-1
                        WRITE(17,"(I5)",advance='no') nI(j)
                    enddo
                    WRITE(17,"(I5)") nI(NEl)
                enddo

                CLOSE(17)

            ENDIF
        ENDIF

    END SUBROUTINE PrintFCIMCPsi



!This routine will write out the average wavevector from the spawning run up until now.
    SUBROUTINE WriteHistogram()
        use Determinants , only : GetHElement3
        use SystemData , only : BasisFN
        INTEGER :: i,j,bits,iLut(0:NIfD),error,IterRead
        TYPE(BasisFN) :: ISym
        REAL*8 :: norm,norm1,norm2,norm3,ShiftRead,AllERead,NumParts
        TYPE(HElement) :: HEL
        CHARACTER(len=22) :: abstr,abstr2
        LOGICAL :: exists

!This will open a file called SpawnHist-"Iter" on unit number 17.
        abstr=''
        write(abstr,'(I12)') Iter
        abstr='SpawnHist-'//adjustl(abstr)
        IF(iProcIndex.eq.0) THEN
            WRITE(6,*) "Writing out the average wavevector up to iteration number: ", Iter
            CALL FLUSH(6)
        ENDIF

!        IF(nBasis/32.ne.0) THEN
!            CALL Stop_All("WriteHistogram","System is too large to histogram as it stands...")
!        ENDIF

        IF(iProcIndex.eq.0) THEN
            AllHistogram(:)=0.D0
            AllInstHist(:)=0.D0
            AllAvAnnihil(:)=0.D0
            AllInstAnnihil(:)=0.D0
        ENDIF

        CALL MPI_Reduce(Histogram,AllHistogram,Det,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(InstHist,AllInstHist,Det,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(InstAnnihil,AllInstAnnihil,Det,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(AvAnnihil,AllAvAnnihil,Det,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        
        IF(iProcIndex.eq.0) THEN

!            IF(.not.associated(NMRKS)) THEN
!                CALL Stop_All("WriteHistogram","A Full Diagonalization is required in the same calculation before histogramming can occur.")
!            ENDIF

            norm=0.D0
            norm1=0.D0
            norm2=0.D0
            norm3=0.D0
            do i=1,Det
                norm=norm+AllHistogram(i)**2
                norm1=norm1+AllInstHist(i)**2
                norm2=norm2+AllInstAnnihil(i)**2
                norm3=norm3+AllAvAnnihil(i)**2
            enddo
            norm=SQRT(norm)
            norm1=SQRT(norm1)
            norm2=SQRT(norm2)
            norm3=SQRT(norm3)
!            WRITE(6,*) "NORM",norm
            do i=1,Det
                AllHistogram(i)=AllHistogram(i)/norm
                AllInstHist(i)=AllInstHist(i)/norm1
                IF(norm2.ne.0.D0) THEN
                    AllInstAnnihil(i)=AllInstAnnihil(i)/norm2
                ENDIF
                IF(norm3.ne.0.D0) THEN
                    AllAvAnnihil(i)=AllAvAnnihil(i)/norm3
                ENDIF
            enddo
            
            OPEN(17,FILE=abstr,STATUS='UNKNOWN')

            abstr=''
            write(abstr,'(I12)') Iter-iWriteHistEvery
            abstr='Energies-'//adjustl(abstr)

            abstr2=''
            write(abstr2,'(I12)') Iter
            abstr2='Energies-'//adjustl(abstr2)
            OPEN(29,FILE=abstr2,STATUS='NEW')

            INQUIRE(FILE=abstr,EXIST=exists)
            IF(exists) THEN
                OPEN(28,FILE=abstr,STATUS='OLD',POSITION='REWIND',ACTION='READ')
                do while(.true.)
                    READ(28,"(I13,3G25.16)",END=99) IterRead,ShiftRead,AllERead,NumParts
                    WRITE(29,"(I13,3G25.16)") IterRead,ShiftRead,AllERead,NumParts
                enddo
99              CONTINUE
                IF(AllHFCyc.eq.0) THEN
                    WRITE(19,"(I13,3G25.16)") Iter,DiagSft,AllERead,AllTotPartsOld
                ELSE
                    WRITE(29,"(I13,3G25.16)") Iter,DiagSft,AllENumCyc/AllHFCyc,AllTotPartsOld
                ENDIF
                CLOSE(29)
                CLOSE(28)

            ELSE
                OPEN(29,FILE=abstr2,STATUS='UNKNOWN')
                WRITE(29,"(I13,3G25.16)") Iter,DiagSft,AllENumCyc/AllHFCyc,AllTotPartsOld
                CLOSE(29)
            ENDIF


            norm=0.D0
            norm1=0.D0
            do i=1,Det
                norm=norm+(AllHistogram(i))**2
                norm1=norm1+(AllAvAnnihil(i))**2
                WRITE(17,"(I13,6G25.16)") i,AllHistogram(i),norm,AllInstHist(i),AllInstAnnihil(i),AllAvAnnihil(i),norm1
            enddo

!            do i=1,Maxdet
!                bits=0
!                do j=0,nbasis-1
!                    IF(BTEST(i,j)) THEN
!                        Bits=Bits+1
!                    ENDIF
!                enddo
!                IF(Bits.eq.NEl) THEN
!
!                    do j=1,ndet
!                        CALL EncodeBitDet(NMRKS(:,j),iLut,NEl,nBasis/32)
!                        IF(iLut(0).eq.i) THEN
!                            CALL GETSYM(NMRKS(1,j),NEL,G1,NBASISMAX,ISYM)
!                            IF(ISym%Sym%S.eq.0) THEN
!                                Det=Det+1
!                                WRITE(17,"(3I12)",advance='no') Det,iLut(0)
!                                HEL=GetHElement3(NMRKS(:,j),NMRKS(:,j),0)
!                                norm=norm+(AllHistogram(i))**2
!                                norm1=norm1+(AllInstHist(i))**2
!                                WRITE(17,"(5G25.16)") REAL(HEL%v,8),AllHistogram(i),norm,AllInstHist(i),norm1
!                            ENDIF
!                            EXIT
!                        ENDIF
!                    ENDDO
!                ENDIF
!
!            ENDDO
            CLOSE(17)
        ENDIF
        InstHist(:)=0.D0
        InstAnnihil(:)=0.D0

    END SUBROUTINE WriteHistogram


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
    
    
!    SUBROUTINE ChooseACFDets()
!        use SystemData , only : tAssumeSizeExcitgen
!        INTEGER :: ierr,HFConn,nStore(6),i,nJ(NEl),iExcit,j,MaxIndex,ExcitLevel,iGetExcitLevel_2
!        INTEGER :: iMaxExcit,ExcitLength,MinInd
!        REAL*8 , ALLOCATABLE :: TempMax(:),ACEnergy(:)
!        LOGICAL :: TurnBackAssumeExGen
!        TYPE(HElement) :: Hij,Hjj,Fjj,Fkk
!        REAL*8 :: Compt,MaxWeight,MinValue
!        INTEGER , ALLOCATABLE :: ExcitGenTemp(:),ACExcLevel(:)
!        CHARACTER(len=*), PARAMETER :: this_routine='ChooseACFDets'
!
!!Commented out code is to just perform the ACF for the HF determinant.
!!            NoAutoDets=1    !Initially, we are just after the ACF for one determinant
!!Fill the autocorrdets array with the determinants that you want to work out the autocorrelation for.
!!            AutoCorrDets(:,1)=HFDet(:)
!
!!The code calculates the ACF for the number of doubles here, or all if it is more
!!        NoAutoDets=10
!!        WRITE(6,"(A,I5,A)") "Choosing the ",NoAutoDets," highest MP1 weight determinants to calculate the ACF for."
!        NoAutoDets=1+NoACDets(2)+NoACDets(3)+NoACDets(4)    !Total number of dets to histogram is sum of the determinants for the individual excitation levels.
!        WRITE(6,"(A,I5,A)") "Choosing the ",NoACDets(2)+1," highest MP1 weight determinants to calculate the ACF for."
!        WRITE(6,"(A,I5,A,I5,A)") "Also picking ",NoACDets(3), " high weighted triply excited and ",NoACDets(4), " quadruply excited determinants."
!        WRITE(6,*) "First determinant will be the HF determinant"
!        CALL GetSymExcitCount(HFExcit%ExcitData,HFConn)
!        IF(NoAutoDets.gt.HFConn) NoAutoDets=HFConn
!        ALLOCATE(AutoCorrDets(NEl,NoAutoDets),stat=ierr)
!        CALL LogMemAlloc('AutoCorrDets',NEl*NoAutoDets,4,this_routine,AutoCorrDetsTag,ierr)
!        AutoCorrDets(:,:)=0
!        
!!Set the first determinant to find the ACF of to be the HF determinant
!        AutoCorrDets(:,1)=HFDet(:)
!
!        ALLOCATE(ACExcLevel(NoAutoDets))
!        ACExcLevel(:)=0
!        ALLOCATE(ACEnergy(NoAutoDets))
!        AcEnergy(:)=0
!
!!We do not know if tAssumeSizeExcitgen is on - if it is, then we can't enumerate all determinants. Get around this by simply regenerating it anyway.
!!First, we need to turn off AssumeSizeExcitgen if it is on.
!        IF(tAssumeSizeExcitgen) THEN
!            TurnBackAssumeExGen=.true.
!            tAssumeSizeExcitgen=.false.
!        ELSE
!            TurnBackAssumeExGen=.false.
!        ENDIF
!
!        ALLOCATE(TempMax(2:NoACDets(2)+1),stat=ierr)   !This will temporarily hold the largest components (+1 since we are no longer considering HF in NoACDets(2))
!        IF(ierr.ne.0) THEN
!            CALL Stop_All("ChooseACFDets","Problem allocating memory")
!        ENDIF
!        TempMax(:)=0.D0
!
!!Setup excit generators for HF Determinant
!        iMaxExcit=0
!        nStore(1:6)=0
!        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitLength,nJ,iMaxExcit,0,nStore,2)
!        ALLOCATE(ExcitGenTemp(ExcitLength),stat=ierr)
!        IF(ierr.ne.0) CALL Stop_All("ChooseACFDets","Problem allocating excitation generator")
!        ExcitGenTemp(1)=0
!        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGenTemp,nJ,iMaxExcit,0,nStore,2)
!        ACEnergy(1)=0.D0
!
!        do while(.true.)
!!Generate double excitations
!            CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.false.,ExcitGenTemp,nJ,iExcit,0,nStore,2)
!            IF(nJ(1).eq.0) EXIT
!            Hij=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
!!            CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Hjj)
!!            CALL GetH0Element(HFDet,NEl,Arr,nBasis,ECore,Fii)
!            Compt=real(Hij%v,r2)/(Fii-(REAL(Fjj%v,r2)))
!            MinInd=2
!            MinValue=TempMax(2)
!!First need to find the minimum value to swap out
!            do j=3,NoACDets(2)+1    !NoAutoDets
!                IF(MinValue.gt.TempMax(j)) THEN
!                    MinValue=TempMax(j)
!                    MinInd=j
!                ENDIF
!            enddo
!
!!See if the just calculated value of the MP1 component is larger than the one current minimum.
!            IF(abs(Compt).gt.TempMax(MinInd)) THEN
!                Hjj=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                TempMax(MinInd)=abs(Compt)
!                AutoCorrDets(:,MinInd)=nJ(:)
!                ACExcLevel(MinInd)=2
!                ACEnergy(MinInd)=real(Hjj%v,r2)-Hii
!            ENDIF
!
!!            do j=2,NoACDets(2)+1!NoAutoDets
!!                IF(abs(Compt).gt.TempMax(j)) THEN
!!                    TempMax(j)=abs(Compt)
!!                    AutoCorrDets(:,j)=nJ(:)
!!                    ACExcLevel(j)=2
!!                    ACEnergy(j)=real(Fjj%v,r2)
!!                    EXIT
!!                ENDIF
!!            enddo
!            
!        enddo
!        DEALLOCATE(ExcitGenTemp)
!!Find largest weight MP1 contribution
!        MaxWeight=0.D0
!        do i=2,NoACDets(2)+1
!            IF(TempMax(i).gt.MaxWeight) THEN
!                MaxIndex=i
!            ENDIF
!        enddo
!        DEALLOCATE(TempMax)
!        
!!        IF(iProcIndex.eq.root) THEN
!!            WRITE(6,*) "*** Histogramming the following determinants:"
!!            do i=1,NoAutoDets
!!                do j=1,NEl
!!                    WRITE(6,"(I4)",advance='no') AutoCorrDets(j,i)
!!                enddo
!!                WRITE(6,*) ""
!!            enddo
!!        ENDIF
!
!!We have now found the largest weight Doubles. To guess at large weighted triples & quads, we simply take the largest weighted double, and find 
!!its large MP1 weighted contributions with it as the reference.
!!        WRITE(6,*) "MaxIndex = ", MaxIndex
!!Setup excit generators for the highest weighted double excitation
!        iMaxExcit=0
!        nStore(1:6)=0
!        CALL GenSymExcitIt2(AutoCorrDets(:,MaxIndex),NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitLength,nJ,iMaxExcit,0,nStore,3)
!        ALLOCATE(ExcitGenTemp(ExcitLength),stat=ierr)
!        IF(ierr.ne.0) CALL Stop_All("ChooseACFDets","Problem allocating excitation generator")
!        ExcitGenTemp(1)=0
!        CALL GenSymExcitIt2(AutoCorrDets(:,MaxIndex),NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGenTemp,nJ,iMaxExcit,0,nStore,3)
!
!        CALL GetH0Element(AutoCorrDets(:,MaxIndex),NEl,Arr,nBasis,ECore,Fjj)
!        ALLOCATE(TempMax(1:(NoACDets(3)+NoACDets(4))),stat=ierr)
!        TempMax(:)=0.D0
!
!        do while(.true.)
!            CALL GenSymExcitIt2(AutoCorrDets(:,MaxIndex),NEl,G1,nBasis,nBasisMax,.false.,ExcitGenTemp,nJ,iExcit,0,nStore,3)
!            IF(nJ(1).eq.0) EXIT
!            Hij=GetHElement2(AutoCorrDets(:,MaxIndex),nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
!            CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fkk)
!            Compt=real(Hij%v,r2)/(real(Fjj%v,r2)-(real(Fkk%v,r2)))
!            ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,4)
!            IF(ExcitLevel.eq.3) THEN
!!We have generated a triple - try to add it to the list
!                do j=1,NoACDets(3)
!                    IF(abs(Compt).gt.TempMax(j)) THEN
!                        Hjj=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                        TempMax(j)=abs(Compt)
!                        AutoCorrDets(:,(j+1+NoACDets(2)))=nJ(:)
!                        ACExcLevel(j+1+NoACDets(2))=3
!                        ACEnergy(j+1+NoACDets(2))=real(Hjj%v,r2)-Hii
!                        EXIT
!                    ENDIF
!                enddo
!            ELSEIF(ExcitLevel.eq.4) THEN
!!We have generated a quad - try to add it to the list
!                do j=NoACDets(3)+1,(NoACDets(3)+NoACDets(4))
!                    IF(abs(Compt).gt.TempMax(j)) THEN
!                        Hjj=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                        TempMax(j)=abs(Compt)
!                        AutoCorrDets(:,(j+1+NoACDets(2)))=nJ(:)
!                        ACExcLevel(j+1+NoACDets(2))=4
!                        ACEnergy(j+1+NoACDets(2))=real(Hjj%v,r2)-Hii
!                        EXIT
!                    ENDIF
!                enddo
!            ENDIF
!        enddo
!
!!Deallocate TempExcitgen
!        DEALLOCATE(TempMax)
!!        CALL DissociateExitgen(ExcitGenTemp)
!        DEALLOCATE(ExcitGenTemp)
!        
!!Find if any zeros
!        do i=1,NoAutoDets
!            IF(AutoCorrDets(1,i).eq.0) THEN
!!                NoAutoDets=NoAutoDets-1
!                WRITE(6,*) "Could not find AutoCorrFunc determinant to histogram slot ",i
!!Move these zeros to the end of the list
!                do j=i,NoAutoDets
!                    IF(AutoCorrDets(1,j).ne.0) THEN
!                        AutoCorrDets(:,i)=AutoCorrDets(:,j)
!                        ACExcLevel(i)=ACExcLevel(j)
!                        ACEnergy(i)=ACEnergy(j)
!                        EXIT
!                    ENDIF
!                enddo
!                NoAutoDets=NoAutoDets-1
!            ENDIF
!        enddo
!            
!!The number of occurunces of walkers at the determiants selected needs to be stored for all iterations
!        ALLOCATE(WeightatDets(NoAutoDets),stat=ierr)
!        CALL LogMemAlloc('WeightatDets',NoAutoDets,4,this_routine,WeightatDetsTag,ierr)
!        WeightatDets(:)=0
!!If we want to calculate the ACF, then we now do this in a seperate step. Now we simply write out the determinant populations at each iteration
!!to a file called HFDoublePops
!        
!        IF(iProcIndex.eq.root) THEN
!            WRITE(6,*) "Histogramming the following determinants:"
!            do i=1,NoAutoDets
!                do j=1,NEl
!                    WRITE(6,"(I4)",advance='no') AutoCorrDets(j,i)
!                enddo
!                WRITE(6,*) ""
!            enddo
!
!            OPEN(44,FILE='HFDoublePops',STATUS='UNKNOWN')
!            WRITE(44,*) 'Energy of determinant'
!            WRITE(44,"(A8)",advance='no') "-"
!            do i=1,NoAutoDets
!                WRITE(44,"(F20.10)",advance='no') ACEnergy(i)
!            enddo
!            WRITE(44,*) 'Excitation level'
!            WRITE(44,"(A8)",advance='no') "-"
!            do i=1,NoAutoDets
!                WRITE(44,"(I8)",advance='no') ACExcLevel(i)
!            enddo
!            WRITE(44,"(/,A14,5X,A23)") "Iteration No.","Determinant Populations"
!        ENDIF
!
!        DEALLOCATE(ACExcLevel)
!        DEALLOCATE(ACEnergy)
!
!        IF(TurnBackAssumeExGen) THEN
!!We turned off assumed sized excitation generators for this routine - turn it back on.
!            tAssumeSizeExcitgen=.true.
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE ChooseACFDets
!
!
!    
!!This is a routine to create a graph from the walker at nI, and ascribe new walkers at each vertex on the graph, according to a number of applications of the rho matrix.
!    SUBROUTINE ResumGraphPar(nI,WSign,VecSlot,VecInd)
!        INTEGER :: nI(NEl),VecSlot,VecInd
!        LOGICAL :: WSign
!        REAL*8 :: Prob
!
!!This routine will create the graph. It will calculate it differently depending on the size of the graph, and whether excitation generators are stored. 
!        CALL CreateGraphPar(nI,VecInd,Prob)
!
!!Apply the rho matrix successive times. This could be improved if large numbers of applications of rho are needed by diagonalising the rho matrix
!        CALL ApplyRhoMatPar()
!        
!        IF(.not.tRegenDiagHEls) THEN
!            CALL CreateNewPartsPar(nI,VecInd,WSign,CurrentH(VecInd),VecSlot,Prob)   !Create particles proportionally to the magnitude of the vector elements in GraphVec
!        ELSE
!            CALL Stop_All("ResumGraphPar","ResumGraphPar does not work with regendiaghels.")
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE ResumGraphPar
!
!!This routine will create a graph from a given initial determinant
!    SUBROUTINE CreateGraphPar(nI,VecInd,Prob)
!        TYPE(ExcitPointer) :: TempExcitgen
!        INTEGER :: nI(NEl),nJ(NEl),VecInd,i,j,IC,iCount,Attempts
!        REAL*8 :: Prob,ExcitProb
!        LOGICAL :: SameDet,CompiPath
!        TYPE(HElement) :: Hij,Hjj
!
!        IF(tRegenDiagHEls) THEN
!            CALL Stop_All("CreateGraphPar","CreateGraphPar will not work with RegenDiagHEls")
!        ENDIF
!
!        TempExcitgen%PointToExcit=>null()
!
!!First, set up the excitation generators for the root determinant
!        IF(.not.TRegenExcitgens) THEN
!            CALL SetupExitgenPar(nI,CurrentExcits(VecInd))
!        ELSE
!            CALL SetupExitgenPar(nI,TempExcitgen)
!        ENDIF
!
!        IF(NDets.eq.2) THEN
!!We know that determinants are not going to be regenerated if NDets=2, so we can do this in a slightly simpler way
!            
!            IF(TRegenExcitgens) THEN
!                CALL GenRandSymExcitIt3(nI,TempExcitgen%PointToExcit,nJ,0,IC,0,Prob,iCount)
!            ELSE
!                CALL GenRandSymExcitIt3(nI,CurrentExcits(VecInd)%PointToExcit,nJ,0,IC,0,Prob,iCount)
!            ENDIF
!
!            Hij=GetHElement2(nI,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!            Hjj=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!
!            GraphRhoMat(1,1)=1.D0-Tau*(CurrentH(VecInd)-DiagSft)
!            GraphRhoMat(1,2)=-Tau*real(Hij%v,r2)
!            GraphRhoMat(2,1)=GraphRhoMat(1,2)
!            GraphRhoMat(2,2)=1.D0-Tau*((real(Hjj%v,r2)-Hii)-DiagSft)
!                
!            DetsinGraph(:,2)=nJ(:)
!            GraphKii(2)=REAL(Hjj%v,r2)-Hii              !store det generated and kii element
!
!        ELSE
!!Zero the matrix
!            GraphRhoMat(1:NDets,1:NDets)=0.D0
!            
!            GraphRhoMat(1,1)=1.D0-Tau*(CurrentH(VecInd)-DiagSft)    !This is the first rho-matrix element
!
!            i=2
!            Attempts=0
!            do while(i.lt.NDets)    !Loop until all determinants found
!
!                IF(TRegenExcitgens) THEN
!                    CALL GenRandSymExcitIt3(nI,TempExcitgen%PointToExcit,nJ,0,IC,0,Prob,iCount)
!                ELSE
!                    CALL GenRandSymExcitIt3(nI,CurrentExcits(VecInd)%PointToExcit,nJ,0,IC,0,Prob,iCount)
!                ENDIF
!
!                SameDet=.false.
!                do j=2,(i-1)
!                    IF(CompiPath(nJ,DetsinGraph(:,j),NEl)) THEN
!!Determinants are the same as already created determinant - ignore it
!
!                        SameDet=.true.
!                        Attempts=Attempts+1
!                        IF(Attempts.gt.1000) THEN
!                            WRITE(6,*) "iCOUNT IS: ", iCount
!                            CALL FLUSH(6)
!                            CALL Stop_All("CreateGraphPar","More than 1000 attempts needed to grow graph")
!                        ENDIF
!                        EXIT
!                    ENDIF
!                enddo
!
!                IF(.not.SameDet) THEN
!!Store the unbiased probability of generating excitations from this root - check that it is the same as other excits generated
!                    IF(i.eq.2) THEN
!                        ExcitProb=Prob
!                    ELSE
!                        IF(abs(Prob-ExcitProb).gt.1.D-07) THEN
!                            CALL Stop_All("CreateGraph","Excitation probabilities are not uniform - problem here...")
!                        ENDIF
!                    ENDIF
!
!!Determinant is distinct - add it
!                    DetsinGraph(1:NEl,i)=nJ(1:NEl)
!!First find connection to root
!                    Hij=GetHElement2(nI,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!                    GraphRhoMat(1,i)=-Tau*REAL(Hij%v,r2)
!                    GraphRhoMat(i,1)=GraphRhoMat(1,i)
!
!!Then find connection to other determinants
!                    do j=2,(i-1)
!                        Hij=GetHElement2(nJ,DetsInGraph(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,-1,ECore)
!                        GraphRhoMat(i,j)=-Tau*REAL(Hij%v,r2)
!                        GraphRhoMat(j,i)=GraphRhoMat(i,j)
!                    enddo
!
!!Find diagonal element - and store it for later on...
!                    Hjj=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                    GraphKii(i)=REAL(Hjj%v,r2)-Hii                !Again, the root value is not stored
!                    GraphRhoMat(i,i)=1.D0-Tau*(GraphKii(i)-DiagSft)
!
!                    i=i+1
!
!                ENDIF
!            enddo
!
!        ENDIF   !End if ndets=2
!
!        IF(TRegenExcitgens) THEN
!!Deallocate excitation generator if we are regenerating them
!            CALL DissociateExitgen(TempExcitgen)
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE CreateGraphPar
!
!    
!!This routine creates new particles from the vector which results from the 
!    SUBROUTINE CreateNewPartsPar(nI,VecInd,WSign,Kii,VecSlot,Prob)
!        IMPLICIT NONE
!        LOGICAL :: WSign,TempSign
!        INTEGER :: nI(NEl),VecInd
!        INTEGER :: i,j,VecSlot,Create,ExcitLevel,iGetExcitLevel_2
!        INTEGER(KIND=i2) :: HashTemp
!        REAL*8 :: rat,Kii,Kjj,Prob,r
!
!!First deal with the excitations...
!        do i=2,NDets
!!Now create the new particles according the the final vector GraphVec
!            
!            GraphVec(i)=GraphVec(i)/((NDets-1)*Prob)    !Augment the component by the chances of picking that determinant
!    
!            Create=INT(abs(GraphVec(i)))
!            rat=abs(GraphVec(i))-REAL(Create,r2)    !rat is now the fractional part, to be assigned stochastically
!            IF(tMerTwist) THEN
!                CALL genrand_real2(r) 
!            ELSE
!                CALL RANLUX(r,1)
!            ENDIF
!            IF(rat.gt.r) Create=Create+1
!            IF(abs(Create).gt.0) THEN
!                IF(.not.WSign) Create=-Create
!                IF(GraphVec(i).lt.0.D0) Create=-Create
!!Find needed information out about the new particles
!!Calculate excitation level, connection to HF. Diagonal ham element info is already stored
!                
!                ExcitLevel=iGetExcitLevel_2(HFDet,DetsinGraph(:,i),NEl,NEl)
!                IF(ExcitLevel.eq.0) THEN
!                    IF(ABS(GraphKii(i)).gt.1.D-07) THEN
!                        CALL Stop_All("ResumGraphPar","Diagonal K-mat element should be zero for HF particles")
!                    ENDIF
!                ENDIF
!                IF(Create.lt.0) THEN
!                    TempSign=.false.
!                ELSE
!                    TempSign=.true.
!                ENDIF
!                
!                IF(.not.TNoAnnihil) THEN
!!                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
!                    HashTemp=CreateHash(DetsInGraph(:,i))
!                ENDIF
!                
!!Now actually create the particles in NewDets and NewSign
!                do j=1,abs(Create)
!                    NewDets(1:NEl,VecSlot)=DetsInGraph(1:NEl,i)
!                    NewSign(VecSlot)=TempSign
!!                    NewIC(VecSlot)=ExcitLevel
!                    IF(.not.tRegenDiagHEls) NewH(VecSlot)=GraphKii(i)       !Diagonal H El previously stored
!                    IF(.not.TRegenExcitgens) NewExcits(VecSlot)%PointToExcit=>null()
!                    IF(.not.TNoAnnihil) THEN
!!                    IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
!                        Hash2Array(VecSlot)=HashTemp
!                    ENDIF
!                    VecSlot=VecSlot+1
!                enddo
!
!            ENDIF
!
!        enddo
!
!!Now deal with root
!        Create=INT(abs(GraphVec(1)))
!
!        rat=abs(GraphVec(1))-REAL(Create,r2)    !rat is now the fractional part, to be assigned stochastically
!        IF(tMerTwist) THEN
!            CALL genrand_real2(r) 
!        ELSE
!            CALL RANLUX(r,1)
!        ENDIF
!        IF(rat.gt.r) Create=Create+1
!        
!        IF((abs(Create)).gt.0) THEN
!        
!            IF(.not.WSign) Create=-Create
!            IF(GraphVec(1).lt.0.D0) Create=-Create
!
!!Test since the root should not change sign - comment out later
!            IF(WSign.and.(Create.lt.0)) THEN
!                call Stop_All("CreateNewPartsPar","Root determinant should not change sign")
!            ELSEIF((.not.WSign).and.(Create.gt.0)) THEN
!                call Stop_All("CreateNewPartsPar","Root determinant should not change sign")
!            ENDIF
!            
!            IF(Create.lt.0) THEN
!                TempSign=.false.
!            ELSE
!                TempSign=.true.
!            ENDIF
!
!!Now actually create the particles in NewDets and NewSign. It is the same particle as the parent particle.
!            do j=1,abs(Create)
!                NewDets(:,VecSlot)=nI(:)
!                NewSign(VecSlot)=TempSign
!!Copy excitation generator accross
!                IF(j.eq.abs(Create)) THEN
!                    IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(VecInd),NewExcits(VecSlot),.true.)
!                ELSE
!                    IF(.not.TRegenExcitgens) CALL CopyExitgenPar(CurrentExcits(VecInd),NewExcits(VecSlot),.false.)
!                ENDIF
!!                NewIC(VecSlot)=CurrentIC(VecInd)
!                IF(.not.tRegenDiagHEls) NewH(VecSlot)=CurrentH(VecInd)
!                IF(.not.TNoAnnihil) THEN
!!                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
!                    Hash2Array(VecSlot)=HashArray(VecInd)
!                ENDIF
!                VecSlot=VecSlot+1
!            enddo
!
!        ENDIF
!
!        RETURN
!
!    END SUBROUTINE CreateNewPartsPar
!
!!This applies the rho matrix successive times to a root determinant. From this, GraphVec is filled with the correct probabilities for the determinants in the graph
!    SUBROUTINE ApplyRhoMatPar()
!        REAL*8 :: TempVec(NDets)
!        INTEGER :: i,j,k
!
!!        IF(NDets.eq.2) THEN
!!
!!            GraphVec(1)=1.D0
!!            GraphVec(2)=0.D0
!!
!!            do i=1,RhoApp
!!            
!!                TempVec(1)=(GraphRhoMat(1,1)*GraphVec(1))+(GraphRhoMat(1,2)*GraphVec(2))
!!                TempVec(2)=(GraphRhoMat(2,1)*GraphVec(1))+(GraphRhoMat(2,2)*GraphVec(2))
!!            
!!                GraphVec(1)=TempVec(1)
!!                GraphVec(2)=TempVec(2)
!!
!!            enddo
!!
!!        ELSE
!            
!            GraphVec(1:NDets)=0.d0
!            GraphVec(1)=1.D0        !Set the initial vector to be 1 at the root (i.e. for one walker initially)
!        
!            do i=1,RhoApp
!
!                CALL DGEMV('n',NDets,NDets,1.D0,GraphRhoMat,NDets,GraphVec,1,0.D0,TempVec,1)
!                GraphVec(1:NDets)=TempVec(1:NDets)
!                TempVec(1:NDets)=0.d0
!            
!!            do j=1,NDets
!!                TempVec(j)=0.D0
!!                do k=1,NDets
!!                    TempVec(j)=TempVec(j)+GraphRhoMat(j,k)*GraphVec(k)
!!                enddo
!!            enddo
!!            GraphVec(:)=TempVec(:)
!
!            enddo
!
!!        ENDIF
!
!        RETURN
!    END SUBROUTINE ApplyRhoMatPar
!
!
!
!    
!!A routine to annihilate particles separatly on each node. This should mean less annihilation occurs, but it is effect running nProcessors separate simulations.
!!If there are enough particles, then this should be sufficient. Less memory is required, since no hashes need to be stored. Also, no communication is needed,
!!so the routine should scale better as the number of walkers grows.
!    SUBROUTINE AnnihilateonProc(TotWalkersNew)
!        TYPE(ExcitPointer) :: TempExcit
!        REAL*8 :: TempH
!!        INTEGER :: TempIC
!        INTEGER :: TotWalkersNew,j,k,l,DetCurr(0:NIfD),VecSlot,TotWalkersDet
!        INTEGER :: DetLT
!
!        TempExcit%PointToExcit=>null()
!!First, it is necessary to sort the list of determinants
!        CALL SortPartsPar(TotWalkersNew,NewDets(:,1:TotWalkersNew),NIfD+1)
!
!!Once ordered, each block of walkers on similar determinants can be analysed, and the residual walker concentration moved to CurrentDets
!        j=1
!!j is the counter over all uncancelled walkers - it indicates when we have reached the end of the list of total walkers
!!DetCurr is the current determinant
!        DetCurr(:)=NewDets(:,j)
!!        TempIC=NewIC(j)
!        IF(.not.tRegenDiagHEls) TempH=NewH(j)
!        IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(j),TempExcit,.true.) !This will delete what is behind it - is that ok?
!        
!        VecSlot=1
!
!        do while(j.le.TotWalkersNew)
!!Loop over all walkers
!            TotWalkersDet=0
!            do while ((DetLT(NewDets(:,j),DetCurr,(NIfD+1)).eq.0).and.(j.le.TotWalkersNew))
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
!                    CurrentDets(:,VecSlot)=DetCurr(:)
!                    CurrentSign(VecSlot)=1
!!                    CurrentIC(VecSlot)=TempIC
!                    IF(.not.tRegenDiagHEls) CurrentH(VecSlot)=TempH
!                    IF(.not.TRegenExcitgens) CALL CopyExitgenPar(TempExcit,CurrentExcits(VecSlot),.false.)
!                    VecSlot=VecSlot+1
!                enddo
!            ELSE
!!Negative sign particles want to populate this determinant
!                do l=1,abs(TotWalkersDet)
!                    CurrentDets(:,VecSlot)=DetCurr(:)
!                    CurrentSign(VecSlot)=1
!!                    CurrentIC(VecSlot)=TempIC
!                    IF(.not.tRegenDiagHEls) CurrentH(VecSlot)=TempH
!                    IF(.not.TRegenExcitgens) CALL CopyExitgenPar(TempExcit,CurrentExcits(VecSlot),.false.)
!                    VecSlot=VecSlot+1
!                enddo
!            ENDIF
!!Now update the current determinant
!            DetCurr(:)=NewDets(:,j)
!!            TempIC=NewIC(j)
!            IF(.not.tRegenDiagHEls) TempH=NewH(j)
!            IF(.not.TRegenExcitgens) THEN
!                CALL DissociateExitGen(TempExcit)
!                CALL CopyExitGenPar(NewExcits(j),TempExcit,.true.)
!            ENDIF
!        enddo
!!The new number of residual cancelled walkers is given by one less that VecSlot again.
!        TotWalkers=VecSlot-1
!
!        RETURN
!
!    END SUBROUTINE AnnihilateonProc
!
!!This routine sorts the particles before annihilation. It is identical to the routine in the serial version, but the data which is taken is different.
!! Based on SORTI, SORTPARTS sorts arrays of integers, representing the determinant the walkers are on
!! It then takes all the corresponding info with it
!! Dets is the array (length N) of integers to sort
!! NElecs is the length (in numbers of integers) of each element of Dets
!! Vectors of NewXXX will be sorted correspondingly
!    SUBROUTINE SortPartsPar(N,Dets,NElecs)
!        TYPE(ExcitPointer) :: ExcitTemp
!        REAL*8 :: HTemp
!!        INTEGER :: ICTemp
!        INTEGER :: TempDet(NElecs)     !This stores a single element of the vector temporarily     
!        INTEGER :: WSignTemp
!        INTEGER N,I,L,IR,J,NElecs
!        INTEGER Dets(NElecs,N)
!        INTEGER DETLT
!
!        ExcitTemp%PointToExcit=>null()
!        IF(N.LE.1) RETURN
!        L=N/2+1 
!        IR=N
!10      CONTINUE
!        IF(L.GT.1)THEN
!            L=L-1
!            TempDet(:)=Dets(:,L)
!            IF(.not.tRegenDiagHEls) HTemp=NewH(L)
!!            ICTemp=NewIC(L)
!            WSignTemp=NewSign(L)
!            IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(L),ExcitTemp,.true.) !This will delete what is behind it - is this ok?
!        ELSE
!            TempDet(:)=Dets(:,IR)      !Copy IRth elements to temp
!            IF(.not.tRegenDiagHEls) HTemp=NewH(IR)
!!            ICTemp=NewIC(IR)
!            WSignTemp=NewSign(IR)
!            IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(IR),ExcitTemp,.true.)
!
!            Dets(:,IR)=Dets(:,1)    !Copy 1st element to IRth element
!            IF(.not.tRegenDiagHEls) NewH(IR)=NewH(1)
!!            NewIC(IR)=NewIC(1)
!            NewSign(IR)=NewSign(1)
!            IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(1),NewExcits(IR),.true.)
!            IR=IR-1
!            IF(IR.EQ.1)THEN
!                Dets(:,1)=TempDet(:)    !Copy temp element to 1st element
!                IF(.not.tRegenDiagHEls) NewH(1)=HTemp
!!                NewIC(1)=ICTemp
!                NewSign(1)=WSignTemp
!                IF(.not.TRegenExcitgens) CALL CopyExitgenPar(ExcitTemp,NewExcits(1),.true.)
!                RETURN
!            ENDIF
!        ENDIF
!        I=L
!        J=L+L
!20      IF(J.LE.IR)THEN
!            IF(J.LT.IR)THEN
!                IF((DETLT(Dets(1,J),Dets(1,J+1),NElecs)).eq.-1) J=J+1
!            ENDIF
!            IF((DETLT(TempDet,Dets(1,J),NElecs)).eq.-1)THEN
!                Dets(:,I)=Dets(:,J)     !Copy Jth element to Ith element
!                IF(.not.tRegenDiagHEls) NewH(I)=NewH(J)
!!                NewIC(I)=NewIC(J)
!                NewSign(I)=NewSign(J)
!                IF(.not.TRegenExcitgens) CALL CopyExitgenPar(NewExcits(J),NewExcits(I),.true.)
!                I=J
!                J=J+J
!            ELSE
!                J=IR+1
!            ENDIF
!            GO TO 20
!        ENDIF
!        Dets(:,I)=TempDet(:)
!        IF(.not.tRegenDiagHEls) NewH(I)=HTemp
!!        NewIC(I)=ICTemp
!        NewSign(I)=WSignTemp
!        IF(.not.TRegenExcitgens) CALL CopyExitgenPar(ExcitTemp,NewExcits(I),.true.)
!        GO TO 10
!
!    END SUBROUTINE SortPartsPar



    SUBROUTINE NormandDiagOneRDM()
!This routine takes the OneRDM with non-normalised off-diagonal elements, adds the diagonal elements from the HF determinant
!and normalises it according to the number of walkers on the HF determinant.
!It then diagonalises the 1-RDM to find linear combinations of the HF orbitals that are closer to the natural orbitals,
!and the occupation numbers of these new orbitals (e-values).
        USE NatOrbsMod , only : FindNatOrbsOld
        INTEGER :: i,j,error,ierr
        REAL*8 :: TempSumNoatHF
        REAL*8 , ALLOCATABLE :: Temp1RDM(:,:)

        TempSumNoatHF=real(SumNoatHF,r2)
!Sum TempSumNoatHF over all processors and then send to all processes
        CALL MPI_AllReduce(TempSumNoatHF,AllSumNoatHF,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
        ALLOCATE(Temp1RDM(nBasis,nBasis),stat=ierr)
        IF(ierr.ne.0) THEN
            CALL Stop_All("NormandDiagOneRDM","Could not allocate Temp1RDM")
        ENDIF

        CALL MPI_AllReduce(OneRDM(:,:),Temp1RDM(:,:),nBasis*nBasis,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
        OneRDM(:,:)=Temp1RDM(:,:)
        DEALLOCATE(Temp1RDM)
        
!Normalise the off-diag OneRDM elements using the number of walkers on the HF determinant
        do i=1,nBasis
            do j=i+1,nBasis
                OneRDM(i,j)=OneRDM(i,j)/AllSumNoatHF
                OneRDM(j,i)=OneRDM(i,j)
            enddo
        enddo
!Using the decoded version of the HF determinant, place values of 1.D0 in the diagonal elements
!of the 1-RDM corresponding to the occupied orbitals.
        do i=1,NEl
            OneRDM(HFDet(i),HFDet(i))=1.D0
        enddo
    
        CALL FindNatOrbsOld()           !Diagonalise the 1-RDM

    END SUBROUTINE NormandDiagOneRDM
    
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
    
!This is the same as BinSearchParts1, but this time, it searches though the full list of determinants created by the full diagonalizer when the histogramming option is on.
    SUBROUTINE BinSearchParts2(iLut,MinInd,MaxInd,PartInd,tSuccess)
        use DetCalc , only : FCIDets
        INTEGER :: iLut(0:NIfD),MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp,DetBitLT
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(FCIDets(:,MaxInd),iLut(:),NIfD)
            IF(Comp.eq.0) THEN
                tSuccess=.true.
                PartInd=MaxInd
                RETURN
            ELSE
                tSuccess=.false.
                PartInd=MinInd
            ENDIF
        ENDIF
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(FCIDets(:,N),iLut(:),NIfD)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
                    Comp=DetBitLT(FCIDets(:,i+1),iLut(:),NIfD)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN
                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchParts2

#else
! AJWT
! Bringing you a better FciMCPar.  A vision for the future...
!
!  This section contains parts of FciMCPar which are not dependent on MPI commands.
!  It's not yet complete, but at least compiles and runs

    SUBROUTINE FciMCPar(Weight,Energyxw)
    TYPE(HDElement) :: Weight,Energyxw

        CALL Stop_All("FciMCPar","Entering the wrong FCIMCPar parallel routine")

    END SUBROUTINE FciMCPar
#endif
! AJWT
! Bringing you a better FciMCPar.  A vision for the future...
!
!  This section contains parts of FciMCPar which are not dependent on MPI commands.
!  It's not yet complete, but at least compiles and runs

    SUBROUTINE SetupParameters()
        use SystemData, only : tUseBrillouin,iRanLuxLev,tSpn,tHPHFInts,tRotateOrbs,tNoBrillouin,tROHF,tFindCINatOrbs
        USE mt95 , only : genrand_init
        use CalcData, only : EXCITFUNCS,tFCIMC
        use Calc, only : VirtCASorbs,OccCASorbs,FixShift,G_VMC_Seed
        use Determinants , only : GetH0Element3
        use SymData , only : nSymLabels,SymLabelList,SymLabelCounts
        use Logging , only : tTruncRODump
        use GenRandSymExcitNUMod , only : SpinOrbSymSetup,tNoSingsPossible
        use FciMCLoggingMOD , only : InitTriHElStats,InitSpinCoupHel
        use DetCalc, only : NMRKS,tagNMRKS,FCIDets
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet,HFDetTest(NEl),Seed,alpha,beta,symalpha,symbeta,endsymstate
        INTEGER :: DetLT,VecSlot,error,HFConn,MemoryAlloc,iMaxExcit,nStore(6),nJ(Nel),BRR2(nBasis),LargestOrb,nBits,HighEDet(NEl),iLutTemp(0:NIfD)
        TYPE(HElement) :: rh,TempHii
        REAL*8 :: TotDets,SymFactor,Choose
        CHARACTER(len=*), PARAMETER :: this_routine='SetupParameters'
        CHARACTER(len=12) :: abstr
        LOGICAL :: tSuccess,tFoundOrbs(nBasis),tTurnBackBrillouin
        REAL :: Gap

!        CALL MPIInit(.false.)       !Initialises MPI - now have variables iProcIndex and nProcessors
        WRITE(6,*) ""
        WRITE(6,*) "Performing Parallel FCIMC...."
        
!Set timed routine names
        Walker_Time%timer_name='WalkerTime'
        Annihil_Time%timer_name='AnnihilTime'
        Sort_Time%timer_name='SortTime'
        Comms_Time%timer_name='CommsTime'
        ACF_Time%timer_name='ACFTime'
        AnnSpawned_time%timer_name='AnnSpawnedTime'
        AnnMain_time%timer_name='AnnMainTime'
        BinSearch_time%timer_name='BinSearchTime'

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
        
!test the encoding of the HFdet to bit representation.
        ALLOCATE(iLutHF(0:NIfD),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Cannot allocate memory for iLutHF")
        CALL EncodeBitDet(HFDet,iLutHF,NEl,NIfD)
!Test that the bit operations are working correctly...
        CALL DecodeBitDet(HFDetTest,iLutHF,NEl,NIfD)
        do i=1,NEl
            IF(HFDetTest(i).ne.HFDet(i)) THEN
                WRITE(6,*) "HFDet: ",HFDet(:)
                WRITE(6,*) "HFDetTest: ",HFDetTest(:)
                CALL Stop_All(this_routine,"HF Determinant incorrectly decoded.")
            ENDIF
        enddo
        CALL LargestBitSet(iLutHF,NIfD,LargestOrb)
        IF(LargestOrb.ne.HFDet(NEl)) THEN
            CALL Stop_All(this_routine,"LargestBitSet FAIL")
        ENDIF
        CALL CountBits(iLutHF,NIfD,nBits,NEl)
        IF(nBits.ne.NEl) THEN
            CALL Stop_All(this_routine,"CountBits FAIL")
        ENDIF

!Check that the symmetry routines have set the symmetry up correctly...
        tSuccess=.true.
        tFoundOrbs(:)=.false.

        IF(.not.tHub) THEN
            do i=1,nSymLabels
!                WRITE(6,*) "NSymLabels: ",NSymLabels,i-1
                EndSymState=SymLabelCounts(1,i)+SymLabelCounts(2,i)-1
!                WRITE(6,*) "Number of states: ",SymLabelCounts(2,i)
                do j=SymLabelCounts(1,i),EndSymState

                    Beta=(2*SymLabelList(j))-1
                    Alpha=(2*SymLabelList(j))
                    SymAlpha=INT((G1(Alpha)%Sym%S),4)
                    SymBeta=INT((G1(Beta)%Sym%S),4)
!                    WRITE(6,*) "***",Alpha,Beta

                    IF(.not.tFoundOrbs(Beta)) THEN
                        tFoundOrbs(Beta)=.true.
                    ELSE
                        CALL Stop_All("SetupParameters","Orbital specified twice")
                    ENDIF
                    IF(.not.tFoundOrbs(Alpha)) THEN
                        tFoundOrbs(Alpha)=.true.
                    ELSE
                        CALL Stop_All("SetupParameters","Orbital specified twice")
                    ENDIF

                    IF(G1(Beta)%Ms.ne.-1) THEN
                        tSuccess=.false.
                    ELSEIF(G1(Alpha)%Ms.ne.1) THEN
                        tSuccess=.false.
                    ELSEIF((SymAlpha.ne.(i-1)).or.(SymBeta.ne.(i-1))) THEN
                        tSuccess=.false.
                    ENDIF
                enddo
            enddo
            do i=1,nBasis
                IF(.not.tFoundOrbs(i)) THEN
                    WRITE(6,*) "Orbital: ",i, " not found."
                    CALL Stop_All("SetupParameters","Orbital not found")
                ENDIF
            enddo
        ENDIF
        IF(.not.tSuccess) THEN
            WRITE(6,*) "************************************************"
            WRITE(6,*) "**                 WARNING!!!                 **"
            WRITE(6,*) "************************************************"
            WRITE(6,*) "Symmetry information not set up correctly in NECI initialisation"
            WRITE(6,*) "Will attempt to set up the symmetry again, but now in terms of spin orbitals"
            WRITE(6,*) "I strongly suggest you check that the reference energy is correct."
!            CALL Stop_All("SetupParameters","Error in the setup of the symmetry/spin ordering of the orbitals. This configuration will not work with spawning excitation generators")
            CALL SpinOrbSymSetup(.true.) 
        ELSE
            WRITE(6,*) "Symmetry and spin of orbitals correctly set up for spawning excitation generators."
!            CALL SpinOrbSymSetup(.false.) 
!            CALL Stop_All("SSS","SKCJB")
        ENDIF

!Check whether it is possible to have a determinant where the electrons
!can be arranged in a determinant so that there are no unoccupied
!orbitals with any of the irreps of the occupied orbitals. If this can happen, we need to check for it before generating excitations.
!NEED TO CREATE A TEST HERE
       tNoSingsPossible=.true. 


!Setup excitation generator for the HF determinant. If we are using assumed sized excitgens, this will also be assumed size.
        IF(tUseBrillouin.and.tNonUniRandExcits) THEN
            WRITE(6,*) "Temporarily turning brillouins theorem off in order to calculate pDoubles for non-uniform excitation generators"
            tTurnBackBrillouin=.true.
            tUseBrillouin=.false.
        ELSE
            tTurnBackBrillouin=.false.
        ENDIF
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,HFExcit%nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
        ALLOCATE(HFExcit%ExcitData(HFExcit%nExcitMemLen),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating excitation generator")
        HFExcit%ExcitData(1)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,HFExcit%ExcitData,nJ,iMaxExcit,0,nStore,3)
        HFExcit%nPointed=0
!Indicate that the excitation generator is now correctly allocated.

!        CALL SetupExitgenPar(HFDet,HFExcit)
        CALL GetSymExcitCount(HFExcit%ExcitData,HFConn)

        IF(tTurnBackBrillouin) THEN
            tUseBrillouin=.true.
            WRITE(6,*) "Turning back on brillouins theorem"
        ENDIF

!Initialise random number seed - since the seeds need to be different on different processors, subract processor rank from random number
        Seed=abs(G_VMC_Seed-iProcIndex)
        WRITE(6,*) "Value for seed is: ",Seed
!Initialise...
        IF(tMerTwist) THEN
            CALL genrand_init(Seed)
        ELSE
            CALL RLUXGO(iRanLuxLev,Seed,0,0)
        ENDIF
        
        IF(tHPHF) THEN
            IF(tROHF.or.(LMS.ne.0)) CALL Stop_All("SetupParameters","Cannot use HPHF with high-spin systems.")
            tHPHFInts=.true.
        ENDIF

!Calculate Hii
        IF(tHPHF) THEN
            CALL HPHFGetDiagHElement(HFDet,iLutHF,TempHii)
        ELSE
            TempHii=GetHElement2(HFDet,HFDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)
        ENDIF
        Hii=REAL(TempHii%v,r2)
        WRITE(6,*) "Reference Energy set to: ",Hii
        TempHii=GetH0Element3(HFDet)
        Fii=REAL(TempHii%v,r2)

!Find the highest energy determinant...
        IF(.not.tSpn) THEN
            do i=1,NEl
                HighEDet(i)=Brr(nBasis-(i-1))
            enddo
            IF(tHPHF) THEN
                CALL EncodeBitDet(HighEDet,iLutTemp,NEl,NIfD)
                CALL HPHFGetDiagHElement(HighEDet,iLutTemp,TempHii)
            ELSE
                TempHii=GetHElement2(HighEDet,HighEDet,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
            ENDIF
            WRITE(6,"(A,G25.15)") "Highest energy determinant is (approximately): ",TempHii%v
            WRITE(6,"(A,F25.15)") "This means tau should be no more than about ",-2.D0/TempHii%v
!            WRITE(6,*) "Highest energy determinant is: ", HighEDet(:)
        ENDIF

        IF(tHub) THEN
            IF(tReal) THEN
!We also know that in real-space hubbard calculations, there are only single excitations.
                exFlag=1
            ELSE
!We are doing a momentum space hubbard calculation - set exFlag to 2 since only doubles are connected for momentum conservation.
                exFlag=2
            ENDIF
        ENDIF


!        IF(tSpawnSymDets) THEN
!!This option will spawn on determinants where the alpha and beta strings are swapped for S=0 RHF systems.
!!These determinants should have the same amplitude in the CI wavefunction - see Helgakker for details.
!            IF(tSpn) THEN
!                CALL Stop_All("SetupParameters","SpawnSymDets cannot work with ROHF or UHF systems (currently?)")
!            ENDIF
!            IF(.not.tRotoAnnihil) THEN
!                CALL Stop_All("SetupParameters","SpawnSymDets must be used with RotoAnnihilation currently")
!            ENDIF
!            WRITE(6,*) "Spawning on symmetric determinants for each spawning step"
!        ENDIF

        IF(LMS.ne.0) THEN
            IF(tNoBrillouin.or.(tHub.and.tReal).or.tRotatedOrbs) THEN
                WRITE(6,*) "High spin calculation with single excitations also used to calculate energy."
            ELSE
                WRITE(6,*) "WARNING!! High-spin calculation detected but single excitations will *not* be used to calculate energy."
                WRITE(6,*) "This is ok for UHF, but not ROHF."
            ENDIF
!            tRotatedOrbs=.true.
!        ELSEIF(LMS.ne.0) THEN
!            CALL Stop_All(this_routine,"Ms not equal to zero, but tSpn is false. Error here")
        ENDIF

        TBalanceNodes=.false.   !Assume that the nodes are initially load-balanced

!Initialise variables for calculation on each node
        IterTime=0.0
        ProjectionE=0.D0
        AvSign=0.D0
        AvSignHFD=0.D0
        SumENum=0.D0
        SumNoatHF=0
        NoatHF=0
        NoatDoubs=0
!        DetsNorm=0.D0
!        MeanExcitLevel=0.D0
!        MinExcitLevel=NEl+10
!        MaxExcitLevel=0
        LocalAnn=0
        Annihilated=0
        Annihilated=0
        Acceptances=0
        PreviousCycles=0
        NoBorn=0
        SpawnFromSing=0
        NoDied=0
        HFCyc=0
        HFPopCyc=0
        ENumCyc=0.D0
        ProjEIter=0.D0
        ProjEIterSum=0.D0
        GuideFuncDoub=0.D0
        GuideFuncHF=0
        VaryShiftCycles=0
        AvDiagSft=0.D0
        SumDiagSft=0.D0

!Also reinitialise the global variables - should not necessarily need to do this...
        AllSumENum=0.D0
        AllNoatHF=0
        AllNoatDoubs=0
        AllSumNoatHF=0.D0
        AllGrowRate=0.D0
!        AllMeanExcitLevel=0.D0
        AllSumWalkersCyc=0
        AllAvSign=0.D0
        AllAvSignHFD=0.D0
        AllNoBorn=0
        AllSpawnFromSing=0
        AllNoDied=0
        AllLocalAnn=0
        AllAnnihilated=0
        AllMinorAnnihilated=0
        AllENumCyc=0.D0
        AllHFCyc=0.D0
!        AllDetsNorm=0.D0
        tCleanRun=.false.
        CullInfo(1:10,1:3)=0
        NoCulls=0
        

        IF((tHistSpawn.or.tCalcFCIMCPsi).and.tFCIMC) THEN
            ALLOCATE(HistMinInd(NEl))
            ALLOCATE(HistMinInd2(NEl))
            maxdet=0
            do i=1,nel
                maxdet=maxdet+2**(nbasis-i)
            enddo

            IF(.not.allocated(FCIDets)) THEN
                CALL Stop_All(this_routine,"A Full Diagonalization is required in the same calculation before histogramming can occur.")
            ENDIF

            WRITE(6,*) "Histogramming spawning wavevector, with Dets=", Det
            ALLOCATE(Histogram(1:det),stat=ierr)
            IF(ierr.ne.0) THEN
                CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
            ENDIF
            Histogram(:)=0.D0
            IF(tHistSpawn) THEN
                ALLOCATE(InstHist(1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                ENDIF
                InstHist(:)=0.D0
                ALLOCATE(AvAnnihil(1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                ENDIF
                AvAnnihil(:)=0.D0
                ALLOCATE(InstAnnihil(1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                ENDIF
                InstAnnihil(:)=0.D0
            ENDIF

            ALLOCATE(AllHistogram(1:det),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
            IF(iProcIndex.eq.0) THEN
                IF(tHistSpawn) THEN
                    ALLOCATE(AllInstHist(1:det),stat=ierr)
                    ALLOCATE(AllInstAnnihil(1:det),stat=ierr)
                    ALLOCATE(AllAvAnnihil(1:det),stat=ierr)
                ENDIF
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                ENDIF
            ENDIF
        ELSEIF(tHistEnergies) THEN
            WRITE(6,*) "Histogramming the energies of the particles, with iNoBins=",iNoBins, " and BinRange=", BinRange
            WRITE(6,*) "Histogramming spawning events from ",-OffDiagMax, " with BinRange = ", OffDiagBinRange
            iOffDiagNoBins=INT((2.D0*OffDiagMax)/OffDiagBinRange)+1
            WRITE(6,*) "This gives ",iOffDiagNoBins," bins to histogram the off-diagonal matrix elements."
            ALLOCATE(Histogram(1:iNoBins))
            ALLOCATE(AttemptHist(1:iNoBins))
            ALLOCATE(SpawnHist(1:iNoBins))
            ALLOCATE(SinglesHist(1:iOffDiagNoBins))
            ALLOCATE(SinglesAttemptHist(1:iOffDiagNoBins))
            ALLOCATE(SinglesHistOccOcc(1:iOffDiagNoBins))
            ALLOCATE(SinglesHistOccVirt(1:iOffDiagNoBins))
            ALLOCATE(SinglesHistVirtOcc(1:iOffDiagNoBins))
            ALLOCATE(SinglesHistVirtVirt(1:iOffDiagNoBins))
            ALLOCATE(DoublesHist(1:iOffDiagNoBins))
            ALLOCATE(DoublesAttemptHist(1:iOffDiagNoBins))
            Histogram(:)=0.D0
            AttemptHist(:)=0.D0
            SpawnHist(:)=0.D0
            SinglesHist(:)=0.D0
            SinglesAttemptHist(:)=0.D0
            SinglesHistOccOcc(:)=0.D0
            SinglesHistOccVirt(:)=0.D0
            SinglesHistVirtOcc(:)=0.D0
            SinglesHistVirtVirt(:)=0.D0
            DoublesHist(:)=0.D0
            DoublesAttemptHist(:)=0.D0
            IF(iProcIndex.eq.0) THEN
                ALLOCATE(AllHistogram(1:iNoBins))
                ALLOCATE(AllAttemptHist(1:iNoBins))
                ALLOCATE(AllSpawnHist(1:iNoBins))
                ALLOCATE(AllSinglesHist(1:iNoBins))
                ALLOCATE(AllDoublesHist(1:iNoBins))
                ALLOCATE(AllSinglesAttemptHist(1:iNoBins))
                ALLOCATE(AllDoublesAttemptHist(1:iNoBins))
                ALLOCATE(AllSinglesHistOccOcc(1:iNoBins))
                ALLOCATE(AllSinglesHistOccVirt(1:iNoBins))
                ALLOCATE(AllSinglesHistVirtOcc(1:iNoBins))
                ALLOCATE(AllSinglesHistVirtVirt(1:iNoBins))
            ENDIF
        ENDIF

!Need to declare a new MPI type to deal with the long integers we use in the hashing, and when reading in from POPSFILEs
!        CALL MPI_Type_create_f90_integer(18,mpilongintegertype,error)
!        CALL MPI_Type_commit(mpilongintegertype,error)
        IF(tNonUniRandExcits) THEN
!These are the new random excitations which give non-uniform generation probabilities.
            IF(tUseBrillouin) THEN
                WRITE(6,*) "Brillouin theorem specified, but this will not be in use with the non-uniform excitation generators."
            ENDIF
            WRITE(6,*) "Non-uniform excitation generators in use."
            IF(.not.tRegenExcitgens) THEN
                WRITE(6,*) "Currently, we can not store the excitation generators with non-uniform excitations. Instead, they will be regenerated."
                tRegenExcitgens=.true.
            ENDIF
            CALL CalcApproxpDoubles(HFConn)
        ENDIF
        IF((.not.tRegenExcitgens).and.(tRotoAnnihil)) THEN
            CALL Stop_All("SetupParameters","Storage of excitation generators is incompatable with RotoAnnihilation. Regenerate excitation generators.")
        ENDIF

!This is a list of options which cannot be used with the stripped-down spawning routine. New options not added to this routine should be put in this list.
        IF(tHighExcitsSing.or.tHistSpawn.or.tRegenDiagHEls.or.tFindGroundDet.or.tStarOrbs.or.tResumFCIMC.or.tSpawnAsDet.or.tImportanceSample    &
     &      .or.(.not.tRegenExcitgens).or.(.not.tNonUniRandExcits).or.(.not.tDirectAnnihil).or.tMinorDetsStar.or.tSpawnDominant.or.(DiagSft.gt.0.D0).or.   &
     &      tPrintTriConnections.or.tHistTriConHEls.or.tCalcFCIMCPsi.or.tTruncCAS.or.tListDets.or.tPartFreezeCore.or.(iAnnInterval.ne.1)) THEN
            WRITE(6,*) "It is not possible to use to clean spawning routine..."
        ELSE
            WRITE(6,*) "Clean spawning routine in use..."
            tCleanRun=.true.
        ENDIF

        IF(iAnnInterval.ne.1) THEN
            WRITE(6,*) "Only annihilating every",iAnnInterval," steps. This will only be possible with default annihilation."
            WRITE(6,*) "Warning - this is experimental!"
            IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                CALL Stop_All(this_routine,"ANNIHILATEEVERY will only work with default annihilation...")
            ENDIF
        ENDIF

        IF(tListDets) THEN
! When this is on, we have to read the list of determinants which we are allowed to spawn at from a file called SpawnOnlyDets.
#ifdef PARALLEL
            CALL ReadSpawnListDets()
#else
            Call Stop_All("FciMCPar/SetupParameters","AJWT disabled call to ReadSpawnList in Serial version.")
#endif
        ENDIF

        IF(tConstructNOs) THEN
! This is the option for constructing the natural orbitals actually during a NECI calculation.  This is different (and probably a lot more complicated and doesn't 
! currently work) from the FINDCINATORBS option which finds the natural orbitals given a final wavefunction.
            ALLOCATE(OneRDM(nBasis,nBasis),stat=ierr)
            CALL LogMemAlloc('OneRDM',nBasis*nBasis,8,this_routine,OneRDMTag,ierr)
            OneRDM(:,:)=0.D0
        ENDIF

        IF(TPopsFile.and.(mod(iWritePopsEvery,StepsSft).ne.0)) THEN
            CALL Warning(this_routine,"POPSFILE writeout should be a multiple of the update cycle length.")
        ENDIF

        IF(TNoAnnihil) THEN
            WRITE(6,*) "No Annihilation to occur. Results are likely not to converge on right value. Proceed with caution. "
        ELSEIF(TAnnihilonproc) THEN
            CALL Stop_All("SetupParameters","Annihilonproc feature is currently disabled")
            WRITE(6,*) "Annihilation will occur on each processors' walkers only. This should be faster, but result in less annihilation."
            WRITE(6,*) "This is equivalent to running seperate calculations."
        ELSEIF(tRotoAnnihil) THEN
            IF(tDirectAnnihil) THEN
                CALL Stop_All("SetupParameters","Cannot specify both direct annihilation and rotoannihilation.")
            ENDIF
            WRITE(6,*) "RotoAnnihilation in use...!"
        ELSEIF(tDirectAnnihil) THEN
            WRITE(6,*) "Direct Annihilation in use...Explicit load-balancing disabled."
            ALLOCATE(ValidSpawnedList(0:nProcessors-1),stat=ierr)
            ALLOCATE(InitialSpawnedSlots(0:nProcessors-1),stat=ierr)
!InitialSpawnedSlots now holds the first free position in the newly-spawned list for each processor, so it does not need to be reevaluated each iteration.
            MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
            Gap=REAL(MaxSpawned)/REAL(nProcessors)
            do j=0,nProcessors-1
                InitialSpawnedSlots(j)=NINT(Gap*j)+1
            enddo
        ENDIF
        IF(TReadPops) THEN
!List of things that readpops can't work with...
            IF(TStartSinglePart.or.TStartMP1) THEN
                CALL Stop_All("SetupParameters","ReadPOPS cannot work with StartSinglePart or StartMP1")
            ENDIF
        ENDIF

        IF(TResumFciMC) THEN
            CALL Stop_All("SetupParameters","ResumFCIMC code is currently disabled")
            IF(NDets.ge.2) THEN
                IF(.not.EXCITFUNCS(10)) THEN
                    WRITE(6,*) "Cannot have an excitation bias with multiple determinant graphs...exiting."
                    CALL Stop_All("SetupParameters","Cannot have biasing with Graphsizes > 2")
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
                CALL Stop_All("SetupParameters","Graphs cannot be smaller than two vertices")
            ELSEIF(TFixShiftShell.or.tFixShiftKii.or.tFixCASShift) THEN
                CALL Stop_All("SetupParameters","Fixing the shift of the certain excitation levels cannot be used within ResumFCIMC")
            ENDIF
            IF(iProcIndex.eq.root) THEN
                WRITE(6,*) "Resumming in multiple transitions to/from each excitation"
                WRITE(6,"(A,I5,A)") "Graphs to resum will consist of ",NDets," determinants."
            ENDIF
        ENDIF
        
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
                CALL Stop_All("SetupParameters","Cannot read in POPSFILE as well as starting with a single particle")
            ENDIF
            IF(TStartMP1) THEN
                CALL Stop_All("SetupParameters","Cannot start with a single particle, and as the MP1 wavefunction")
            ENDIF
        ELSE
            TSinglePartPhase=.false.
        ENDIF
        IF(TFixShiftShell) THEN
            IF(iProcIndex.eq.root) THEN
                WRITE(6,"(A,I5,A,F20.10)") "All excitations up to ",ShellFix," will have their shift fixed at ",FixShift
                WRITE(6,*) "With this option, results are going to be non-exact"
            ENDIF
        ELSEIF(tFixShiftKii) THEN
            IF(iProcIndex.eq.root) THEN
                WRITE(6,"(A,G25.16,A,F20.10)") "All excitations with Kii values less than ",FixedKiiCutoff," will have their shift fixed at ",FixShift
                WRITE(6,*) "With this option, results are going to be non-exact"
            ENDIF
        ELSEIF(tFixCASShift) THEN
            IF(iProcIndex.eq.root) THEN
                WRITE(6,"(A,I5)") "All determinants containing excitations within the active space of ",OccCASorbs
                WRITE(6,"(A,I5,A)")" highest energy, occupied spin orbitals, and ",VirtCASorbs," lowest energy, "
                WRITE(6,"(A,F20.10)") "virtual spin orbitals, will have their shift fixed at ",FixShift
                WRITE(6,*) "With this option, results are going to be non-exact"
            ENDIF
!The SpinInvBRR array is required for the FixCASShift option. Its properties are explained more fully in the subroutine. 

            CALL CreateSpinInvBRR()

            CASmax=NEl+VirtCASorbs
! CASmax is the max spin orbital number (when ordered energetically) within the chosen active space.
! Spin orbitals with energies larger than this maximum value must be unoccupied for the determinant
! to be in the active space.
            CASmin=NEl-OccCASorbs
! CASmin is the max spin orbital number below the active space.  As well as the above criteria, spin 
! orbitals with energies equal to, or below that of the CASmin orbital must be completely occupied for 
! the determinant to be in the active space.

        ENDIF
        IF(tStarOrbs) THEN
            IF(tImportanceSample.or.(ICILevel.ne.0).or.(.not.tRegenExcitgens)) THEN
                CALL Stop_All("SetupParameters","Cannot use star orbs while storing excitation generators, or truncation or importance sampling...")
            ENDIF
            CALL CreateSpinInvBrr()
        ENDIF
        IF(ICILevel.ne.0) THEN
!We are truncating the excitations at a certain value
            TTruncSpace=.true.
            WRITE(6,'(A,I4)') "Truncating the S.D. space at determinants will an excitation level w.r.t. HF of: ",ICILevel
            IF(TResumFciMC) CALL Stop_All("SetupParameters","Space cannot be truncated with ResumFCIMC")
        ENDIF
        IF(tTruncCAS) THEN
!We are truncating the FCI space by only allowing excitations in a predetermined CAS space.
            WRITE(6,'(A,I4,A,I5)') "Truncating the S.D. space as determinants must be within a CAS of ",OccCASOrbs," , ",VirtCASOrbs
!The SpinInvBRR array is required for the tTruncCAS option. Its properties are explained more fully in the subroutine. 

            CALL CreateSpinInvBRR()

            CASmax=NEl+VirtCASorbs
! CASmax is the max spin orbital number (when ordered energetically) within the chosen active space.
! Spin orbitals with energies larger than this maximum value must be unoccupied for the determinant
! to be in the active space.
            CASmin=NEl-OccCASorbs
! CASmin is the max spin orbital number below the active space.  As well as the above criteria, spin 
! orbitals with energies equal to, or below that of the CASmin orbital must be completely occupied for 
! the determinant to be in the active space.

            IF(OccCASOrbs.gt.NEl) CALL Stop_All("SetupParameters","Occupied orbitals in CAS space specified is greater than number of electrons")
            IF(VirtCASOrbs.gt.(nBasis-NEl)) CALL Stop_All("SetupParameters","Virtual orbitals in CAS space specified is greater than number of unoccupied orbitals")
        ENDIF
        IF(tPartFreezeCore) THEN
            WRITE(6,'(A,I4,A,I5)') 'Partially freezing the lowest ',NPartFrozen,' spin orbitals so that no more than ',NHolesFrozen,' holes exist within this core.'
            CALL CreateSpinInvBRR()
        ENDIF

!        IF(TAutoCorr) THEN
!!We want to calculate the autocorrelation function over the determinants
!            IF(iLagMin.lt.0) THEN
!                CALL Stop_All("SetupParameters","LagMin cannot be less than zero (and when equal 0 should be strictly 1")
!            ELSEIF(iLagMax.gt.NMCyc) THEN
!                CALL Stop_All("SetupParameters","LagMax cannot be greater than the number of cycles.")
!            ELSEIF(iLagStep.lt.1) THEN
!                CALL Stop_All("SetupParameters","LagStep cannot be less than 1")
!            ENDIF
!
!            CALL ChooseACFDets()
!            CALL Stop_All("ChooseACFDets","This code has been commented out...")
!
!!            WRITE(6,*) "Storing information to calculate the ACF at end of simulation..."
!        ENDIF

        IF(tMagnetize) THEN

            IF(tRotoAnnihil) THEN
                CALL Stop_All("SetupParameters","Rotoannihilation not currently supporting Magnetization")
            ENDIF
            CALL FindMagneticDets()
        ENDIF

        IF(TLocalAnnihilation) THEN
!If we are locally annihilating, then we need to know the walker density for a given excitation level, for which we need to approximate number of determinants
!in each excitation level
            CALL Stop_All("SetupParameters","LocalAnnihilation is currently disabled")
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
        IF(tHighExcitsSing) THEN
            WRITE(6,*) "Only allowing single excitations between determinants where one of them has an excitation level w.r.t. HF of more than ",iHighExcitsSing
            IF(iHighExcitsSing.gt.NEl) THEN
                CALL Stop_All("SetupParameters","iHighExcitsSing.ge.NEl")
            ELSEIF(iHighExcitsSing.eq.NEl) THEN
                CALL Warning("SetupParameters","iHighExcitsSing = NEl - this will no longer have any effect.")
            ENDIF
            IF((.not.tNonUniRandExcits).or.tStarOrbs.or.tTruncSpace.or.tTruncCAS.or.tListDets.or.tPartFreezeCore) THEN
                CALL Stop_All("SetupParameters","Cannot use HighExcitsSing without Nonuniformrandexcits, or with starorbs or truncated spaces...")
            ENDIF
        ENDIF

        IF(tMultipleDetsSpawn) THEN
!We need to store a list of all double excitations of HF.
            CALL StoreDoubs()
        ENDIF

    END SUBROUTINE SetupParameters
    LOGICAL FUNCTION TestifDETinCAS(DetCurr)
        INTEGER :: k,z,DetCurr(NEl)

!        CASmax=NEl+VirtCASorbs
! CASmax is the max spin orbital number (when ordered energetically) within the chosen active space.
! Spin orbitals with energies larger than this maximum value must be unoccupied for the determinant
! to be in the active space.
!        CASmin=NEl-OccCASorbs   (These have been moved to the InitCalc subroutine so they're not calculated
! each time.
! CASmin is the max spin orbital number below the active space.  As well as the above criteria, spin 
! orbitals with energies equal to, or below that of the CASmin orbital must be completely occupied for 
! the determinant to be in the active space.

        z=0
        do k=1,NEl      ! running over all electrons
            if (SpinInvBRR(DetCurr(k)).gt.CASmax) THEN
                TestifDETinCAS=.false.
                EXIT            
! if at any stage an electron has an energy greater than the CASmax value, the determinant can be ruled out
! of the active space.  Upon identifying this, it is not necessary to check the remaining electrons.
            else
                if (SpinInvBRR(DetCurr(k)).le.CASmin) THEN
                    z=z+1
                endif
! while running over all electrons, the number that occupy orbitals equal to or below the CASmin cutoff are
! counted.
            endif
        enddo
! the final number of electrons in this low energy region must be equal to the number of orbitals (CASmin), 
! otherwise an orbital is unoccupied, and the determinant cannot be part of the active space.
        if (z.eq.CASmin) THEN
            TestifDETinCAS=.true.
        else
            TestifDETinCAS=.false.
        endif
        
        RETURN

    END FUNCTION TestifDETinCAS

!This function will tell us whether we should allow attempted spawning at an excitation when we are truncating the space.
!We pass in the excitation level of the original particle, the two representations of the excitation (we only need the bit-representation of the excitation
!for HPHF) and the magnitude of the excitation (for determinant representation).
    LOGICAL FUNCTION CheckAllowedTruncSpawn(WalkExcitLevel,nJ,iLutnJ,IC)
        INTEGER :: nJ(NEl),WalkExcitLevel,iLutnJ(0:NIfD),ExcitLevel,IC,iGetExcitLevel_2,i,NoInFrozenCore
        LOGICAL :: DetBitEQ

        CheckAllowedTruncSpawn=.true.

        IF(tTruncSpace) THEN
!We are truncating the space by excitation level
            IF(tHPHF) THEN
!With HPHF, we can't rely on this, since one excitation could be a single, and one a double. Also, IC is not returned.
                CALL FindBitExcitLevel(iLutHF,iLutnJ,NIfD,ExcitLevel,ICILevel)
                IF(ExcitLevel.gt.ICILevel) THEN
                    CheckAllowedTruncSpawn=.false.
                ENDIF

            ELSE
!Determinant representation.

                IF(WalkExcitLevel.eq.(ICILevel-1)) THEN
!The current walker is one below the excitation cutoff - if IC is a double, then could go over - we need to check
                    
                    IF(IC.eq.2) THEN
                        ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,ICILevel)
                    ELSE
!Always allow this - a single cannot put us over the truncated excitation level
                        ExcitLevel=0
!                        CheckAllowedTruncSpawn=.true.
!                        RETURN
                    ENDIF
                    IF(ExcitLevel.gt.ICILevel) THEN
                        CheckAllowedTruncSpawn=.false.
!                    ELSE
!                        CheckAllowedTruncSpawn=.true.
                    ENDIF

                ELSEIF(WalkExcitLevel.ge.ICILevel) THEN
!Walker is at the excitation cutoff level - all possible excitations could be disallowed - check the actual excitation level
                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,ICILevel)
                    IF(ExcitLevel.gt.ICILevel) THEN
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                        CheckAllowedTruncSpawn=.false.
!                    ELSE
!                        CheckAllowedTruncSpawn=.true.
                    ENDIF
!                ELSE
!!Excitation cannot be in a dissallowed excitation level - allow it as normal
!                    CheckAllowedTruncSpawn=.true.
                ENDIF 

            ENDIF   !endif tHPHF

        ENDIF   !endif tTruncSpace

        IF(tTruncCAS.and.CheckAllowedTruncSpawn) THEN
!This flag determines if the FCI space is restricted by whether the determinants are in the predescribed CAS.
            IF(.not.TestifDetinCAS(nJ)) THEN
!Excitation not in allowed CAS space.
!                WRITE(6,*) "Not in CAS:",nJ(:)
                CheckAllowedTruncSpawn=.false.
!            ELSE
!                WRITE(6,*) "In Cas:",nJ(:)
            ENDIF

        ENDIF

        IF(tListDets.and.CheckAllowedTruncSpawn) THEN
!This will check to see if the determinants are in a list of determinants in AllowedDetList
            CheckAllowedTruncSpawn=.false.
            IF(.not.tHPHF) THEN
                CALL EncodeBitDet(nJ,iLutnJ,NEl,NIfD)
            ENDIF
            do i=1,NAllowedDetList
!                WRITE(6,*) ILutnJ,AllowedDetList(0:NIfD,i)
                IF(DetBitEQ(iLutnJ,AllowedDetList(0:NIfD,i),NIfD)) THEN
!                    WRITE(6,*) "Allowed Det"
                    CheckAllowedTruncSpawn=.true.
                    EXIT
                ENDIF
            enddo
        ENDIF

        IF(tPartFreezeCore) THEN
!Want to check if the determinant we're about to spawn on has more than the restricted number of holes in the partially frozen core.            

!Run through the electrons in nJ, count the number in the partially frozen core - ie those occupying orbitals with energy (from BRR) 
!less than that of the partially frozen core limit.
!If this is less than NPartFrozen-NHolesFrozen then spawning is forbidden.
            NoInFrozenCore=0
!BRR(i)=j: orbital i is the j-th lowest in energy  
            do i=1,NEl
                IF(SpinInvBRR(nJ(i)).le.NPartFrozen) NoInFrozenCore=NoInFrozenCore+1
                IF(NoInFrozenCore.eq.(NPartFrozen-NHolesFrozen)) EXIT   ! Can exit out of the loop if this is satisfied, since excitation will definitely be accepted.
            enddo
            IF(NoInFrozenCore.lt.(NPartFrozen-NHolesFrozen)) THEN
!There are more holes in the partially frozen core than has been specified as allowed.
                CheckAllowedTruncSpawn=.false.
            ELSE
!Either the 'partially frozen core' is completely full, or it has the allowed number of holes or less.                
!Allowed to spawn, CheckAllowedTruncSpawn=.true.
                CheckAllowedTruncSpawn=.true.
            ENDIF

        ENDIF

    END FUNCTION CheckAllowedTruncSpawn
!This is the same as BinSearchParts1, but this time, the list to search is passed in as an argument. The list goes from 1 to Length, but only between MinInd and MaxInd is actually searched.
    SUBROUTINE BinSearchParts3(iLut,List,Length,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: iLut(0:NIfD),MinInd,MaxInd,PartInd
        INTEGER :: List(0:NIfD,Length),Length
        INTEGER :: i,j,N,Comp,DetBitLT
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(List(:,MaxInd),iLut(:),NIfD)
            IF(Comp.eq.0) THEN
                tSuccess=.true.
                PartInd=MaxInd
                RETURN
            ELSE
                tSuccess=.false.
                PartInd=MinInd
            ENDIF
        ENDIF
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(List(:,N),iLut(:),NIfD)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
                    Comp=DetBitLT(List(:,i+1),iLut(:),NIfD)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN
                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchParts3
    FUNCTION CreateHash(DetCurr)
        INTEGER :: DetCurr(NEl),i
        INTEGER(KIND=i2) :: CreateHash

        CreateHash=0
        do i=1,NEl
!            CreateHash=13*CreateHash+i*DetCurr(i)
            CreateHash=(1099511628211_8*CreateHash)+i*DetCurr(i)
            
!            CreateHash=mod(1099511628211*CreateHash,2**64)
!            CreateHash=XOR(CreateHash,DetCurr(i))
        enddo
!        WRITE(6,*) CreateHash
        RETURN

    END FUNCTION CreateHash
    SUBROUTINE CalcApproxpDoubles(HFConn)
        use SystemData , only : tAssumeSizeExcitgen
        use SymData , only : SymClassSize
        INTEGER :: HFConn,PosExcittypes,iTotal,i
        INTEGER :: nSing,nDoub,ExcitInd

        IF(tHub) THEN
            IF(tReal) THEN
                WRITE(6,*) "Since we are using a real-space hubbard model, only single excitations are connected."
                WRITE(6,*) "Setting pDoub to 0.D0"
                pDoubles=0.D0
                RETURN
            ELSE
                WRITE(6,*) "Since we are using a momentum-space hubbard model, only double excitaitons are connected."
                WRITE(6,*) "Setting pDoub to 1.D0"
                pDoubles=1.D0
                RETURN
            ENDIF
        ENDIF

        WRITE(6,"(A)") "Calculating approximate pDoubles for use with excitation generator by looking a excitations from HF."
        IF(tAssumeSizeExcitgen) THEN
            PosExcittypes=SymClassSize*NEL+NBASIS/32+4
            iTotal=HFExcit%ExcitData(1)
        ELSE
            PosExcittypes=HFExcit%ExcitData(2)
            iTotal=HFExcit%ExcitData(23)
        ENDIF
        IF(iTotal.ne.HFConn) THEN
            CALL Stop_All("CalcApproxpDoubles","Number of excitations from HF determinant has been confused somewhere...")
        ENDIF
!        WRITE(6,*) "**********"
!        do i=0,100
!            WRITE(6,*) i,HFExcit%ExcitData(PosExcittypes+i)
!        enddo
        nSing=0
        nDoub=0
        i=iTotal
!ExcitTypes is normally a (5,NExcitTypes) array, so we have to be a little careful with the indexing.
        ExcitInd=4+PosExcittypes

        do while((i.gt.0).and.HFExcit%ExcitData(ExcitInd).le.i)
            IF(HFExcit%ExcitData(ExcitInd-4).eq.1) THEN
!We are counting single excitations
                NSing=NSing+HFExcit%ExcitData(ExcitInd)
            ELSEIF(HFExcit%ExcitData(ExcitInd-4).eq.2) THEN
                NDoub=NDoub+HFExcit%ExcitData(ExcitInd)
            ELSE
                CALL Stop_All("CalcApproxpDoubles","Cannot read excittypes")
            ENDIF
            i=i-HFExcit%ExcitData(ExcitInd)
            ExcitInd=ExcitInd+5
        enddo

        WRITE(6,"(I7,A,I7,A)") NDoub, " double excitations, and ",NSing," single excitations found from HF. This will be used to calculate pDoubles."

        IF(SinglesBias.ne.1.D0) THEN
            WRITE(6,*) "Singles Bias detected. Multiplying single excitation connectivity of HF determinant by ",SinglesBias," to determine pDoubles."
        ENDIF

        IF((NSing+nDoub).ne.iTotal) THEN
            CALL Stop_All("CalcApproxpDoubles","Sum of number of singles and number of doubles does not equal total number of excitations")
        ENDIF
        IF((NSing.eq.0).or.(NDoub.eq.0)) THEN
            WRITE(6,*) "Number of singles or doubles found equals zero. pDoubles will be set to 0.95. Is this correct?"
            pDoubles=0.95
            RETURN
        ELSEIF((NSing.lt.0).or.(NDoub.lt.0)) THEN
            CALL Stop_All("CalcApproxpDoubles","Number of singles or doubles found to be a negative number. Error here.")
        ENDIF

!Set pDoubles to be the fraction of double excitations.
        pDoubles=real(nDoub,r2)/((real(NSing,r2)*SinglesBias)+real(NDoub,r2))

        IF(SinglesBias.ne.1.D0) THEN
            WRITE(6,"(A,F14.6,A,F14.6)") "pDoubles set to: ",pDoubles, " rather than (without bias): ",real(nDoub,r2)/real(iTotal,r2)
        ELSE
            WRITE(6,"(A,F14.6)") "pDoubles set to: ",pDoubles
        ENDIF

    END SUBROUTINE CalcApproxpDoubles
    SUBROUTINE CreateSpinInvBRR()
    ! Create an SpinInvBRR containing spin orbitals, 
    ! unlike 'createInvBRR' which only has spatial orbitals.
    ! This is used for the FixCASshift option in establishing whether or not
    ! a determinant is in the complete active space.
    ! In:
    !    BRR(i)=j: orbital i is the j-th lowest in energy.
    !    nBasis: size of basis
    ! SpinInvBRR is the inverse of BRR.  SpinInvBRR(j)=i: the j-th lowest energy
    ! orbital corresponds to the i-th orbital in the original basis.
    ! i.e the position in SpinInvBRR now corresponds to the orbital number and 
    ! the value to the relative energy of this orbital. 
    
        IMPLICIT NONE
        INTEGER :: I,t,ierr
        CHARACTER(len=*), PARAMETER :: this_routine='CreateSpinInvBrr'
            
        ALLOCATE(SpinInvBRR(NBASIS),STAT=ierr)
        CALL LogMemAlloc('SpinInvBRR',NBASIS,4,this_routine,SpinInvBRRTag,ierr)
            

!        IF(iProcIndex.eq.root) THEN
!            WRITE(6,*) "================================"
!            WRITE(6,*) "BRR is "
!            WRITE(6,*) BRR(:)
!        ENDIF
        
        SpinInvBRR(:)=0
        
        t=0
        DO I=1,NBASIS
            t=t+1
            SpinInvBRR(BRR(I))=t
        ENDDO

!        IF(iProcIndex.eq.root) THEN
!            WRITE(6,*) "================================"
!            WRITE(6,*) "SpinInvBRR is "
!            WRITE(6,*) SpinInvBRR(:)
!        ENDIF
        
        RETURN
        
    END SUBROUTINE CreateSpinInvBRR




! This creates a hash based not only on one current determinant, but is also dependent on the 
! determinant from which the walkers on this determinant came.
    FUNCTION CreateHashBit(DetCurr)
        INTEGER :: DetCurr(0:NIfD),i,Elecs,j
        INTEGER(KIND=i2) :: CreateHashBit

        CreateHashBit=0
        Elecs=0
        do i=0,NIfD
            do j=0,31
                IF(BTEST(DetCurr(i),j)) THEN
                    CreateHashBit=(1099511628211_8*CreateHashBit)+((i*32)+(j+1))*j
                    Elecs=Elecs+1
                    IF(Elecs.eq.NEl) RETURN
                ENDIF
            enddo
        enddo

    END FUNCTION CreateHashBit


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
            IF(BackOfList.eq.MaxWalkersPart) THEN
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
            IF(FrontOfList.eq.MaxWalkersPart) THEN
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
    SUBROUTINE GetPartRandExcitPar(DetCurr,iLutDet,nJ,IC,Frz,Prob,iCount,ExcitLevel,Ex,tParity)
        use GenRandSymExcitNUMod , only : GenRandSymExcitNU
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,Frz,iCount,iMaxExcit,nStore(6),MemLength,ierr
        INTEGER :: Excitlevel,Ex(2,2),iLutDet(0:NIfD)
        REAL*8 :: Prob
        LOGICAL :: tParity
        INTEGER , ALLOCATABLE :: ExcitGenTemp(:)

        IF(tNonUniRandExcits) THEN
!Generate non-uniform random excitations
            CALL GenRandSymExcitNU(DetCurr,iLutDet,nJ,pDoubles,IC,Ex,tParity,exFlag,Prob)

        ELSE
            IF(ExcitLevel.eq.0) THEN
!            CALL GenRandSymExcitIt3(DetCurr,HFExcit%ExcitData,nJ,Seed,IC,Frz,Prob,iCount)
                CALL GenRandSymExcitIt4(DetCurr,HFExcit%ExcitData,nJ,0,IC,Frz,Prob,iCount,Ex,tParity)
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
            CALL GenRandSymExcitIt4(DetCurr,ExcitGenTemp,nJ,0,IC,Frz,Prob,iCount,Ex,tParity)

!Deallocate when finished
            DEALLOCATE(ExcitGenTemp)
        ENDIF

        RETURN

    END SUBROUTINE GetPartRandExcitPar
    SUBROUTINE FindMagneticDets()
        use SystemData , only : tAssumeSizeExcitgen
        INTEGER :: j,i,ierr,iMaxExcit,nStore(6),MinIndex,ExcitLength,nJ(NEl),iExcit
        TYPE(HElement) :: Hij,Fjj,Kiitemp
        LOGICAL :: TurnBackAssumeExGen
        INTEGER , ALLOCATABLE :: ExcitGenTemp(:)
        REAL*8 , ALLOCATABLE :: TempMax(:)
        REAL*8 :: MP1Energy,Compt,Kii,MinCompt
        CHARACTER(len=*), PARAMETER :: this_routine='FindMagneticDets'

        WRITE(6,"(A,I8,A)") "Finding the sign of the ",NoMagDets," largest weighted MP1 determinants..."
        IF(tSymmetricField) THEN
            WRITE(6,"(A,F14.6)") "Magnetized determinants will raise/lower the energy of anti-parallel/parallel particles by ",BField
        ELSE
            WRITE(6,"(A,F14.6)") "Magnetized determinants will only raise the energy of anti-parallel particles by ",BField
        ENDIF
        IF(NoMagDets.lt.1) THEN
            CALL Stop_All("FindMagneticDets","Number of determinant signs to find < 1 - exiting...")
        ENDIF
        IF(NoMagDets.eq.1) THEN
            WRITE(6,*) "Only fixing the sign of HF determinant"
        ENDIF
        IF(BField.lt.0.D0) THEN
            CALL Stop_All("FindMagneticDets","Magnetic field cannot be negative...")
        ENDIF

!First allocate memory for chosen determinants. HF path is already stored and has a positive sign by definition, so only store NoMagDets-1 of them
        ALLOCATE(MagDets(NEl,NoMagDets-1),stat=ierr)
        CALL LogMemAlloc('MagDets',NEl*(NoMagDets-1),4,this_routine,MagDetsTag,ierr)
        ALLOCATE(MagDetsSign(NoMagDets-1),stat=ierr)
        CALL LogMemAlloc('MagDetsSign',NoMagDets-1,4,this_routine,MagDetsSignTag,ierr)

!Do an MP1 calculation to find determinants to fix the sign of...
!We do not know if tAssumeSizeExcitgen is on - if it is, then we can't enumerate all determinants. Get around this by simply regenerating it anyway.
!First, we need to turn off AssumeSizeExcitgen if it is on.
        IF(tAssumeSizeExcitgen) THEN
            TurnBackAssumeExGen=.true.
            tAssumeSizeExcitgen=.false.
        ELSE
            TurnBackAssumeExGen=.false.
        ENDIF

        ALLOCATE(TempMax(NoMagDets-1),stat=ierr)   !This will temporarily hold the largest components
        IF(ierr.ne.0) THEN
            CALL Stop_All(this_routine,"Problem allocating memory")
        ENDIF
        TempMax(:)=0.D0
        MP1Energy=0.D0

!Setup excit generators for HF Determinant
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitLength,nJ,iMaxExcit,0,nStore,2)
        ALLOCATE(ExcitGenTemp(ExcitLength),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating excitation generator")
        ExcitGenTemp(1)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGenTemp,nJ,iMaxExcit,0,nStore,2)

        i=0
        do while(.true.)
!Generate double excitations
            CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.false.,ExcitGenTemp,nJ,iExcit,0,nStore,2)
            IF(nJ(1).eq.0) EXIT
            i=i+1
            Hij=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
            CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fjj)
            Compt=real(Hij%v,r2)/(Fii-(REAL(Fjj%v,r2)))
            MP1Energy=MP1Energy+((real(Hij%v,r2)**2)/(Fii-(REAL(Fjj%v,r2))))
!Find position of minimum MP1 component stored
            MinCompt=abs(TempMax(1))
            MinIndex=1
            do j=2,NoMagDets-1
                IF(MinCompt.gt.abs(TempMax(j))) THEN
                    MinIndex=j
                    MinCompt=abs(TempMax(j))
                ENDIF
            enddo

!Compare the minimum index MP1 component to the one found from this excitation generated.
            IF(abs(Compt).gt.MinCompt) THEN
                TempMax(MinIndex)=Compt
                MagDets(:,MinIndex)=nJ(:)
                IF(Compt.lt.0) THEN
                    MagDetsSign(MinIndex)=-1
                ELSE
                    MagDetsSign(MinIndex)=1
                ENDIF
            ENDIF
            
        enddo
        WRITE(6,*) i," double excitations found from HF"
        IF(i.lt.NoMagDets-1) THEN
            CALL Stop_All(this_routine,"Not enough double excitations to satisfy number of magnetic determinants requested.")
        ENDIF
        DEALLOCATE(ExcitGenTemp)

        WRITE(6,*) "Determinants picked for magnetisation are (Det   MP1Comp   OrigKii   Kij) :"
        CALL WRITEDET(6,HFDet,NEl,.false.)
        WRITE(6,"(2F14.6)") 1.D0,0.D0
        do j=1,NoMagDets-1
            CALL WRITEDET(6,MagDets(:,j),NEl,.false.)
            Kiitemp=GetHElement2(MagDets(:,j),MagDets(:,j),NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
            Kii=REAL(Kiitemp%v,r2)-Hii
            Hij=GetHElement2(MagDets(:,j),HFDet,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,2,ECore)
            WRITE(6,"(3F14.6)") TempMax(j),Kii,REAL(Hij%v,r2)
        enddo

        WRITE(6,*) "MP1 ENERGY is: ", MP1Energy

        DEALLOCATE(TempMax)
        
        IF(TurnBackAssumeExGen) THEN
!We turned off assumed sized excitation generators for this routine - turn it back on.
            tAssumeSizeExcitgen=.true.
        ENDIF

        RETURN

    END SUBROUTINE FindMagneticDets
!This will store all the double excitations.
    SUBROUTINE StoreDoubs()
        use SystemData , only : tAssumeSizeExcitgen,tUseBrillouin
        INTEGER :: iMaxExcit,nStore(6),ExcitLength,nJ(NEl),ierr,iExcit,VecSlot
        INTEGER , ALLOCATABLE :: ExcitGenTemp(:)

        IF(tAssumeSizeExcitgen) THEN
            CALL Stop_All("StoreDoubs","Cannot have assumed sized excitation generators for full enumeration of determinants")
        ENDIF
        IF(tUseBrillouin) THEN
            CALL Stop_All("StoreDoubs","Cannot have Brillouin theorem as now storing singles too...")
        ENDIF

        
!NoDoubs here is actually the singles + doubles of HF
        CALL GetSymExcitCount(HFExcit%ExcitData,NoDoubs)

        ALLOCATE(DoublesDets(NEl,NoDoubs),stat=ierr)
        CALL LogMemAlloc('DoublesDets',NoDoubs*NEl,4,"StoreDoubs",DoublesDetsTag,ierr)
        DoublesDets(1:NEl,1:NoDoubs)=0
        
        VecSlot=1           !This is the next free slot in the DoublesDets array

!Setup excit generators for HF Determinant
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitLength,nJ,iMaxExcit,0,nStore,3)
        ALLOCATE(ExcitGenTemp(ExcitLength),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("InitWalkersMP1","Problem allocating excitation generator")
        ExcitGenTemp(1)=0
        CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.TRUE.,ExcitGenTemp,nJ,iMaxExcit,0,nStore,3)

        do while(.true.)
!Generate double excitations
            CALL GenSymExcitIt2(HFDet,NEl,G1,nBasis,nBasisMax,.false.,ExcitGenTemp,nJ,iExcit,0,nStore,3)
            IF(nJ(1).eq.0) EXIT
!            IF(iExcit.ne.2) THEN
!                CALL Stop_All("StoreDoubles","Error - excitations other than doubles being generated in DoublesDets wavevector code")
!            ENDIF

            DoublesDets(1:NEl,VecSlot)=nJ(:)
            VecSlot=VecSlot+1

        enddo

!This means that now NoDoubs is double excitations AND singles
!        NoDoubs=VecSlot-1

        IF(VecSlot.ne.(NoDoubs+1)) THEN
            WRITE(6,*) VecSlot,NoDoubs
            CALL Stop_All("StoreDoubs","Problem enumerating all double excitations")
        ENDIF

        DEALLOCATE(ExcitGenTemp)

    END SUBROUTINE StoreDoubs

END MODULE FciMCParMod
    
