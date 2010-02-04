#include "macros.h"
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
    use SystemData , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,nMsh,Arr,LMS,tHPHF,tListDets,NIfD,NIfTot,NIfDBO,NIfY
    use SystemData , only : tHub,tReal,tNonUniRandExcits,tMerTwist,tRotatedOrbs,tImportanceSample,tFindCINatOrbs,tFixLz,LzTot,tUEG, tLatticeGens,tCSF
    use CalcData , only : InitWalkers,NMCyc,DiagSft,Tau,SftDamp,StepsSft,OccCASorbs,VirtCASorbs,tFindGroundDet,tDirectAnnihil
    use CalcData , only : TStartMP1,NEquilSteps,TReadPops,TRegenExcitgens,TFixShiftShell,ShellFix,FixShift,tMultipleDetsSpawn
    use CalcData , only : tConstructNOs,tAnnihilatebyRange,tRotoAnnihil,MemoryFacSpawn,tRegenDiagHEls,tSpawnAsDet
    use CalcData , only : GrowMaxFactor,CullFactor,TStartSinglePart,ScaleWalkers,Lambda,TLocalAnnihilation,tNoReturnStarDets
    use CalcData , only : NDets,RhoApp,TResumFCIMC,TNoAnnihil,MemoryFacPart,TAnnihilonproc,MemoryFacAnnihil,iStarOrbs,tAllSpawnStarDets
    use CalcData , only : FixedKiiCutoff,tFixShiftKii,tFixCASShift,tMagnetize,BField,NoMagDets,tSymmetricField,tStarOrbs,SinglesBias
    use CalcData , only : tHighExcitsSing,iHighExcitsSing,tFindGuide,iGuideDets,tUseGuide,iInitGuideParts,tNoDomSpinCoup
    use CalcData , only : tPrintDominant,iNoDominantDets,MaxExcDom,MinExcDom,tSpawnDominant,tMinorDetsStar,MaxNoatHF,HFPopThresh
    use CalcData , only : tCCMC,tTruncCAS,tTruncInitiator,tDelayTruncInit,IterTruncInit,NShiftEquilSteps,tWalkContGrow
    use HPHFRandExcitMod , only : FindExcitBitDetSym,GenRandHPHFExcit,GenRandHPHFExcit2Scratch
    USE Determinants , only : FDet,GetHElement2,GetHElement4
    USE DetCalc , only : ICILevel,nDet,Det,FCIDetIndex
    use GenRandSymExcitNUMod , only : GenRandSymExcitScratchNU,GenRandSymExcitNU,ScratchSize
    use IntegralsData , only : fck,NMax,UMat,tPartFreezeCore,NPartFrozen,NHolesFrozen,tPartFreezeVirt,NVirtPartFrozen,NElVirtFrozen
    USE UMatCache , only : GTID
    USE Logging , only : iWritePopsEvery,TPopsFile,TZeroProjE,iPopsPartEvery,tBinPops,tHistSpawn,iWriteHistEvery,tHistEnergies,IterShiftBlock
    USE Logging , only : NoACDets,BinRange,iNoBins,OffDiagBinRange,OffDiagMax,tPrintSpinCoupHEl!,iLagMin,iLagMax,iLagStep,tAutoCorr
    USE Logging , only : tPrintTriConnections,tHistTriConHels,tPrintHElAccept,tPrintFCIMCPsi,tCalcFCIMCPsi,NHistEquilSteps,tPrintOrbOcc,StartPrintOrbOcc
    USE Logging , only : tHFPopStartBlock,tIterStartBlock,IterStartBlocking,HFPopStartBlocking,tInitShiftBlocking,tHistHamil,iWriteHamilEvery
    USE Logging , only : OrbOccs,OrbOccsTag,tPrintPopsDefault 
    USE SymData , only : nSymLabels
    USE mt95 , only : genrand_real2
    USE Parallel
    USE FciMCData
    USE AnnihilationMod
    use DetBitops, only: EncodeBitDet, DecodeBitDet, DetBitEQ, DetBitLT
    use DetBitOps, only: FindExcitBitDet, FindBitExcitLevel
    use csf, only: get_csf_bit_yama
    IMPLICIT NONE
    integer, parameter :: dp = selected_real_kind(15,307)
    SAVE

    contains

#ifdef PARALLEL

    SUBROUTINE FciMCPar(Weight,Energyxw)
        use soft_exit, only : ChangeVars 
        use CalcData, only : iFullSpaceIter
        use UMatCache, only : UMatInd
        use FciMCLoggingMOD , only : PrintTriConnHist,PrintTriConnHElHist,FinaliseBlocking,FinaliseShiftBlocking,PrintShiftBlocking,PrintBlocking
        use RotateOrbsMod , only : RotateOrbs
        use NatOrbsMod , only : PrintOrbOccs
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
                    IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                        CALL WriteToPopsfileParOneArr()
                    ELSE
                        CALL WriteToPopsFilePar()
                    ENDIF
                ENDIF
                IF(tSingBiasChange) THEN
                    IF(.not.tNoSpinSymExcitgens) CALL GetSymExcitCount(HFExcit%ExcitData,HFConn)
                    CALL CalcApproxpDoubles(HFConn)
                ENDIF
            
                IF(mod(Iter,StepsSft*100).eq.0) THEN
                    !Every 100 update cycles, write out a new blocking file.
                    IF(tErrorBlocking.and.(Iter.gt.IterStartBlocking)) CALL PrintBlocking(Iter) 
                    IF(tShiftBlocking.and.(Iter.gt.(VaryShiftIter+IterShiftBlock))) CALL PrintShiftBlocking(Iter)
                ENDIF

            ENDIF

            IF(TPopsFile.and.(.not.tPrintPopsDefault).and.(mod(Iter,iWritePopsEvery).eq.0)) THEN
!This will write out the POPSFILE if wanted
                IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                    CALL WriteToPopsfileParOneArr()
                ELSE
                    CALL WriteToPopsfilePar()
                ENDIF
            ENDIF
!            IF(TAutoCorr) CALL WriteHistogrammedDets()

            IF(tHistSpawn.and.(mod(Iter,iWriteHistEvery).eq.0)) THEN
                CALL WriteHistogram()
            ENDIF
            IF(tHistHamil.and.(mod(Iter,iWriteHamilEvery).eq.0)) THEN
                CALL WriteHamilHistogram()
            ENDIF

            Iter=Iter+1
!End of MC cycle
        enddo

        IF(TIncrement) Iter=Iter-1     !Reduce the iteration count for the POPSFILE since it is incremented upon leaving the loop (if done naturally)
        IF(TPopsFile) THEN
            IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                CALL WriteToPopsfileParOneArr()
            ELSE
                CALL WriteToPopsfilePar()
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

        IF(tHistHamil) CALL WriteHamilHistogram()

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
        IF(tPrintOrbOcc) THEN
            CALL PrintOrbOccs(OrbOccs)
            DEALLOCATE(OrbOccs)
            CALL LogMemDeAlloc(this_routine,OrbOccsTag)
        ENDIF
     

!        IF(TAutoCorr) CALL CalcAutoCorr()

!Deallocate memory
        CALL DeallocFCIMCMemPar()

        IF(iProcIndex.eq.Root) THEN
            CLOSE(15)
            IF(tTruncInitiator.or.tDelayTruncInit) CLOSE(16)
!            IF(TAutoCorr) CLOSE(44)
        ENDIF
        IF(TDebug) CLOSE(11)

        RETURN

    END SUBROUTINE FciMCPar

    
!This is the heart of FCIMC, where the MC Cycles are performed. However, this version is clean and does not have unnecessary tests for experimental options.
    SUBROUTINE PerformCleanFCIMCycPar()
        INTEGER :: VecSlot,i,j,k,l,CopySign,iPartBloom
        INTEGER :: nJ(NEl),ierr,IC,Child,DetCurr(NEl),iLutnJ(0:NIfTot)
        REAL*8 :: Prob,rat,HDiagCurr
        INTEGER :: iDie,WalkExcitLevel,Proc
        INTEGER :: ExcitLevel,TotWalkersNew,iGetExcitLevel_2,Ex(2,2),WSign,p,Scratch1(ScratchSize),Scratch2(ScratchSize)
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
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j))

!Also, we want to find out the excitation level of the determinant - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.
            IF(tTruncSpace) THEN
!We need to know the exact excitation level for truncated calculations.
                WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j),&
                                                   nel)
            ELSE
                WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j),2)
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

                IF(IsNullDet(nJ)) THEN
                    Child=0
!Calculate number of children to spawn
                ELSEIF(TTruncSpace) THEN
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
!        use GenRandSymExcitCSF, only: TestCSF123
        use GenRandSymExcitCSF, only: GenRandSymCSFExcit
        USE CalcData , only : tAddtoInitiator,InitiatorWalkNo,tInitIncDoubs
        use detbitops, only : countbits
        INTEGER :: MinorVecSlot,VecSlot,i,j,k,l,MinorValidSpawned,ValidSpawned,CopySign,ParticleWeight,Loop,iPartBloom
        INTEGER :: nJ(NEl),ierr,IC,Child,iCount,DetCurr(NEl),iLutnJ(0:NIfTot),NoMinorWalkersNew
        REAL*8 :: Prob,rat,HDiag,HDiagCurr
        INTEGER :: iDie,WalkExcitLevel,Proc             !Indicated whether a particle should self-destruct on DetCurr
        INTEGER :: ExcitLevel,TotWalkersNew,iGetExcitLevel_2,error,length,temp,Ex(2,2),WSign,p,Scratch1(ScratchSize),Scratch2(ScratchSize),Scratch3(Scratchsize),FDetSym,FDetSpin
        LOGICAL :: tParity,tMainArr,tFilled,tCheckStarGenDet,tStarDet,tMinorDetList,tAnnihilateMinorTemp,tAnnihilateMinor,TestClosedShellDet,tParentInCAS
        INTEGER(KIND=i2) :: HashTemp
        TYPE(HElement) :: HDiagTemp
        CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: message
        REAL :: Gap

        IF(TDebug.and.(mod(Iter,10).eq.0)) THEN
            WRITE(11,*) Iter,TotWalkers,NoatHF,NoatDoubs,MaxIndex,TotParts
            CALL FLUSH(11)
        ENDIF

        IF(tDelayTruncInit.and.(Iter.ge.IterTruncInit)) THEN 
            IF(Iter.eq.IterTruncInit) THEN
                IF(iProcIndex.eq.root) THEN
                    Tau=Tau/10.D0
                    WRITE(6,'(A,F10.5)') 'Beginning truncated initiator calculation and reducing tau by a factor of 10. New tau is : ',Tau
                ENDIF
                CALL MPI_BCast(Tau,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,error)
            ENDIF
            tTruncInitiator=.true.
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

        CALL InitHistMin()

        do j=1,TotWalkers
!j runs through all current walkers
!If we are rotoannihilating/direct annihilating, the sign indicates the sum of the signs on the determinant, and hence j loops over determinants, not particles.
            !WRITE(6,*) Iter,j,TotWalkers
            CALL FLUSH(6)

!First, decode the bit-string representation of the determinant the walker is on, into a string of naturally-ordered integers
!            WRITE(6,*) 'CurrentDet (bit)',CurrentDets(:,j)
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j))

!            IF((Iter.gt.100)) THEN!.and.(.not.DetBitEQ(CurrentDets(:,j),iLutHF))) THEN
!This will test the excitation generator for HPHF wavefunctions
!                IF(.not.(TestClosedShellDet(CurrentDets(:,j)))) THEN
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

            ! This is where you can call the CSF testing routine
            ! call TestCSF123 (DetCurr)

            IF(tTruncSpace.or.tHighExcitsSing.or.tHistSpawn.or.tCalcFCIMCPsi.or.tPrintSpinCoupHEl.or.tHistHamil) THEN
!We need to know the exact excitation level for truncated calculations.
                WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j),&
                                                   nel)
                IF((WalkExcitLevel.eq.2).and.tPrintSpinCoupHEl) CALL FindSpinCoupHEl(iLutHF,CurrentDets(:,j))
            ELSE
                WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j),2)
            ENDIF

            IF(tTruncInitiator) THEN
                IF(tTruncCAS) THEN
                    tParentInCAS=.true.
                    tParentInCAS=TestIfDetInCAS(DetCurr)
                    IF(tParentInCAS) THEN
                        ParentInitiator=0
                        NoInitDets=NoInitDets+1
                        NoInitWalk=NoInitWalk+(ABS(CurrentSign(j)))
!The parent walker from which we are attempting to spawn is in the active space - all children will carry this flag, and these spawn like usual.
                    ELSEIF(tInitIncDoubs.and.(WalkExcitLevel.eq.2)) THEN
                        ParentInitiator=0
                        NoInitDets=NoInitDets+1
                        NoInitWalk=NoInitWalk+(ABS(CurrentSign(j)))
                        NoExtraInitDoubs=NoExtraInitDoubs+1
                    ELSEIF(tAddtoInitiator.and.(ABS(CurrentSign(j)).gt.InitiatorWalkNo)) THEN
                        ParentInitiator=0
                        IF(mod(Iter,StepsSft).eq.0) NoAddedInitiators=NoAddedInitiators+1
                        NoInitDets=NoInitDets+1
                        NoInitWalk=NoInitWalk+(ABS(CurrentSign(j)))
                    ELSE
                        ParentInitiator=1
                        NoNonInitDets=NoNonInitDets+1
                        NoNonInitWalk=NoNonInitWalk+(ABS(CurrentSign(j)))
!The parent from which we are attempting to spawn is outside the active space - children spawned on unoccupied determinants with this flag will be killed.
                    ENDIF
                ELSEIF(tTruncSpace) THEN
                    IF(WalkExcitLevel.le.ICILevel) THEN
                        ParentInitiator=0
                        NoInitDets=NoInitDets+1
                        NoInitWalk=NoInitWalk+(ABS(CurrentSign(j)))
!Parent in allowed space.                        
                    ELSEIF(tInitIncDoubs.and.(WalkExcitLevel.eq.2)) THEN
                        ParentInitiator=0
                        NoInitDets=NoInitDets+1
                        NoInitWalk=NoInitWalk+(ABS(CurrentSign(j)))
                        NoExtraInitDoubs=NoExtraInitDoubs+1
                    ELSEIF(tAddtoInitiator.and.(ABS(CurrentSign(j)).gt.InitiatorWalkNo)) THEN
                        ParentInitiator=0
                        IF(mod(Iter,StepsSft).eq.0) NoAddedInitiators=NoAddedInitiators+1
                        NoInitDets=NoInitDets+1
                        NoInitWalk=NoInitWalk+(ABS(CurrentSign(j)))
                    ELSE
                        ParentInitiator=1
                        NoNonInitDets=NoNonInitDets+1
                        NoNonInitWalk=NoNonInitWalk+(ABS(CurrentSign(j)))
!Parent outside allowed space.                        
                    ENDIF
                ENDIF
            ENDIF

            IF(tRegenDiagHEls) THEN
!We are not storing the diagonal hamiltonian elements for each particle. Therefore, we need to regenerate them.
!Need to find H-element!
                IF(DetBitEQ(CurrentDets(0:NIfTot,j),iLutHF,NIfDBO).and.(.not.(tHub.and.tReal))) THEN
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
!TODO: This is where projected energy calculated - make sure that it can call
!the right HElement call
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

!Ali wanted this debug line left in: this is an appropriate place to call the histogramming of the excitation generator
!at least for UEG and Hubbard model. - jjs
!        write(6,*) "***** DEBUG *****"
!        CALL TestGenRandSymExcitNU(DetCurr,10000000,0.D0,2,1000000)
!        STOP
!            write(6,*) DetCurr

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
                                elseif (tCSF) then
                                    ! TODO: fix this exFlag
                                    exFlag = 7
                                    call GenRandSymCSFExcit (DetCurr, CurrentDets(:,j), nJ, pSingles, pDoubles, IC, Ex, exFlag, Prob, Scratch1, Scratch2, Scratch3, tFilled, tParity)
                                else
                                    CALL GenRandSymExcitScratchNU(DetCurr,CurrentDets(:,j),nJ,pDoubles,IC,Ex,tParity,exFlag,Prob,Scratch1,Scratch2,tFilled)
!                                    WRITE(6,'(A,8I3)') 'determinant generated for spawning',nJ
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
                IF(IsNullDet(nJ)) THEN
                    Child=0
                ELSEIF(((TTruncSpace.or.tTruncCAS).and.(.not.tTruncInitiator)).or.tListDets.or.tPartFreezeCore.or.tPartFreezeVirt.or.tFixLz.or.tUEG) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.
!If we are using the CASStar - all spawns are allowed so no need to check.
                    IF(tImportanceSample) CALL Stop_All("PerformFCIMCyc","Truncated calculations not yet working with importance sampling")
!                    WRITE(6,*) 'cheking if a spawn is allowed'
!                    WRITE(6,*) 'tTruncSpace',tTruncSpace
!                    WRITE(6,*) 'tTruncCAS',tTruncCAS
!                    WRITE(6,*) 'tTruncInitiator',tTruncInitiator

                    IF(CheckAllowedTruncSpawn(WalkExcitLevel,nJ,iLutnJ,IC)) THEN
!                        do i=1,NEl
!                            IF(nJ(i).gt.36) CALL Stop_All('attempt spawning','det not in case space when it says it is')
!                        enddo
                        Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                    ELSE
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
!                        k=0
!                        do i=1,NEl
!                            IF(nJ(i).le.36) THEN
!                                k=k+1
!                            ENDIF
!                        enddo
!                        IF(k.eq.NEl) THEN
!                            WRITE(6,*) 'DetCurr',nJ
!                            WRITE(6,*) 'TestIfDetinCas',tParentInCAS
!                            CALL FLUSH(6)
!                            CALL Stop_All('attempt spawning','det in cas space when it says its not')
!                        ENDIF
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

!                        WRITE(6,'(A,3I20,A,8I3)') 'Child to be created from:',CurrentDets(:,j),' which is ',DetCurr(:)
!                        WRITE(6,'(A,3I20,A,8I3)') 'Child to be created on:',iLutnJ(:),' which is ',nJ(:)
                        !print*, 'attempt create'
                        Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity,ParticleWeight,tMinorDetList)
                        !print*, 'attempted'
                        !if (child /= 0) then
                            !WRITE(6,'(A,3I20)') 'Child to be created on:',iLutnJ(:)
                            !WRITE(6,*) 'Child',Child
                        !endif
                    ENDIF

                ENDIF

! Want to put a wee routine in here that monitors the number of accepted vs the number of not accepted attempts at spawns, and the H elements that
! are involved in each.
                IF(tPrintHElAccept) CALL TrackSpawnAttempts(Child,DetCurr,j,nJ,iLutnJ,IC,Ex,tParity)
                
                IF(Child.ne.0) THEN
!We want to spawn a child - find its information to store
                    IF(tHistHamil) THEN
                        CALL AddHistHamilEl(CurrentDets(:,j),iLutnJ,WalkExcitLevel,Child,1)   !Histogram the hamiltonian - iLutnI,iLutnJ,Excitlevel of parent,spawning indicator
                    ENDIF

!                    WRITE(6,'("Spawning particle to:")',advance='no')
!                    call writedet(6,nJ,nel,.true.)
!                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,NEl)
!                    WRITE(6,*) "Excitlevel:", ExcitLevel
                    NoBorn=NoBorn+abs(Child)     !Update counter about particle birth
                    IF(IC.eq.1) THEN
                        SpawnFromSing=SpawnFromSing+abs(Child)
                    ENDIF

                    IF(tAddtoInitiator.and.(abs(Child).gt.InitiatorWalkNo)) THEN
                        IF(abs(Child).gt.abs(iPartBloom)) THEN
                            IF(IC.eq.1) THEN
                                iPartBloom=-abs(Child)
                            ELSE
                                iPartBloom=abs(Child)
                            ENDIF
                        ENDIF
                    ELSEIF(abs(Child).gt.25) THEN
!If more than 25 particles are created in one go, then log this fact and print out later that this has happened.
                        IF(abs(Child).gt.abs(iPartBloom)) THEN
                            IF(IC.eq.1) THEN
                                iPartBloom=-abs(Child)
                            ELSE
                                iPartBloom=abs(Child)
                            ENDIF
                        ENDIF
!                        WRITE(6,"(A,I10,A)") "LARGE PARTICLE BLOOM - ",Child," particles created in one attempt."
                    ENDIF

                    IF(tRotoAnnihil) THEN
!In the RotoAnnihilation implimentation, we spawn particles into a seperate array - SpawnedParts and SpawnedSign. 
!The excitation level and diagonal matrix element are also found out after the annihilation.
!Cannot use old excitation generators with rotoannihilation.

!In rotoannihilation, we can specify multiple particles on the same entry. 
                        IF(tMinorDetList) THEN
!We want to add the determinants to a seperate list, since they are spawning back from "insignificant" determinants.
                            MinorSpawnDets(0:NIfTot,MinorValidSpawned)=iLutnJ(0:NIfTot)
                            MinorSpawnParent(0:NIfTot,MinorValidSpawned)=CurrentDets(0:NIfTot,j) !This is DetCurr in bit form
                            MinorSpawnSign(MinorValidSpawned)=Child
!                            CALL DecodeBitDet(TempDet,iLutnJ(:))
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

                        !WRITE(6,'(A,3I20)') 'when dealing with directannihilation',iLutnJ(:)
                        Proc=DetermineDetProc(iLutnJ)   !This wants to return a value between 0 -> nProcessors-1
                        !WRITE(6,*) iLutnJ(:),Proc,ValidSpawnedList(Proc),Child,TotWalkers
                        !CALL FLUSH(6)
                        SpawnedParts(:,ValidSpawnedList(Proc))=iLutnJ(:)
                        SpawnedSign(ValidSpawnedList(Proc))=Child
                        IF(tTruncInitiator) SpawnedParts(NIfTot,ValidSpawnedList(Proc))=ParentInitiator
!                        WRITE(6,*) 'SpawnedParts',SpawnedParts(:,ValidSpawnedList(Proc))
                        ValidSpawnedList(Proc)=ValidSpawnedList(Proc)+1
!Set the last integer of the determinant in SpawnedParts to be either 0 or 1 according to whether it's parent is inside or outside the active space.

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
                            ExcitLevel = FindBitExcitLevel(iLutnJ, iLutHF, 2)
                            CALL FindDiagElwithB(HDiag,ExcitLevel,nJ,WSign)
                        ELSE
                            IF(.not.tRegenDiagHEls) THEN
                                IF(DetBitEQ(iLutnJ,iLutHF,NIfDBO)) THEN
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
                        IF(tRotoAnnihil) THEN
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
                        ELSE    !Direct annihilation
                            SpawnedParts(:,ValidSpawnedList(iProcIndex))=CurrentDets(:,j)
                            SpawnedSign(ValidSpawnedList(iProcIndex))=CopySign
                            ValidSpawnedList(iProcIndex)=ValidSpawnedList(iProcIndex)+1
                        ENDIF
                    ENDIF
                ENDIF
            ELSE    !Not rotoannihilation or direct annihilation

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
!                CALL DecodeBitDet(nJ,MinorStarParent(0:NIfTot,i))

                ! CHECK THIS
                Child=AttemptCreateParBack(MinorStarDets(0:NIfTot,i),MinorStarParent(0:NIfTot,i),MinorStarSign(i),MinorStarHij(i),abs(MinorStarSign(i)),tMinorDetList)
                ! This will give an integer which is the number of walkers (w sign) being spawned back to the parent (allowed) determinant.

                IF(tMinorDetList) THEN
                    WRITE(6,*) 'Spawing from',MinorStarDets(0:NIfTot,i)
                    WRITE(6,*) 'attempting to spawn to',MinorStarParent(0:NIfTot,i)
                    CALL Stop_All('PerformFCIMCycPar','ERROR. Attempting to spawn between minor determinants.')
                ENDIF
                ! Check that tMinorDetList is always false (in future change it so that it doesn't search).
                ! If tMinorDetList is true, then nJ is not in the allowed list, and the parent list must be wrong.

!If child.ne.0, then add it, but add it to the normal spawning list (not MinorSpawnDets)
                IF(Child.ne.0) THEN
!                    SpawnedParts(:,ValidSpawned)=iLutnJ(:)
                    SpawnedParts(0:NIfTot,ValidSpawned)=MinorStarParent(0:NIfTot,i)
                    SpawnedSign(ValidSpawned)=Child
                    ValidSpawned=ValidSpawned+1     !Increase index of spawned particles
                ENDIF
           
!Attempt Die for particles on "insignificant" dets

                ! DetCurr is the current determinant in expanded form, MinorStarDets(:,i) is the bit form.
                CALL DecodeBitDet(DetCurr,MinorStarDets(0:NIfTot,i))

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
            IF(tAddtoInitiator) THEN
                WRITE(6,"(A,I10,A)") "Particle Blooms of more than 'n_add' in iteration ",Iter
            ELSE
                WRITE(6,"(A,I10,A)") "LARGE Particle Blooms in iteration ",Iter
            ENDIF
            IF(iPartBloom.gt.0) THEN
                WRITE(6,"(A,I10,A)") "A max of ",abs(iPartBloom)," particles created in one attempt from double excit."
            ELSE
                WRITE(6,"(A,I10,A)") "A max of ",abs(iPartBloom)," particles created in one attempt from single excit."
            ENDIF
        ENDIF


        rat=REAL(TotWalkersNew,r2)/REAL(MaxWalkersPart,r2)
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

            CALL AnnihilatePartPar(TotWalkersNew)
            Annihilated=Annihilated+(TotWalkersNew-TotWalkers)
            TotParts=TotWalkers

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
                        
    
    
!Every StepsSft steps, update the diagonal shift value (the running value for the correlation energy)
!We don't want to do this too often, since we want the population levels to acclimatise between changing the shifts
    SUBROUTINE CalcNewShift()
        USE FciMCLoggingMOD , only : PrintSpawnAttemptStats,PrintTriConnStats,PrintSpinCoupHEl,InitErrorBlocking,SumInErrorContrib
        USE FciMCLoggingMOD , only : InitShiftErrorBlocking,SumInShiftErrorContrib
        INTEGER :: error,rc,MaxAllowedWalkers,MaxWalkersProc,MinWalkersProc
        INTEGER :: inpair(9),outpair(9),inpairInit(8),outpairInit(8)
        REAL*8 :: TempTotWalkers,TempTotParts
        REAL*8 :: TempSumNoatHF,MeanWalkers,TempSumWalkersCyc,TempAllSumWalkersCyc,TempNoMinorWalkers
        REAL*8 :: inpairreal(3),outpairreal(3)
        LOGICAL :: TBalanceNodesTemp,tReZeroShift

        TotImagTime=TotImagTime+StepsSft*Tau

!Find the total number of particles at HF (x sign) across all nodes. If this is negative, flip the sign of all particles.
        AllNoatHF=0
        tReZeroShift=.false.

!Find sum of noathf, and then use an AllReduce to broadcast it to all nodes
        CALL MPIISum(NoatHF,1,AllNoatHF)

        IF(AllNoatHF.lt.0) THEN
!Flip the sign if we're beginning to get a negative population on the HF
            WRITE(6,*) "No. at HF < 0 - flipping sign of entire ensemble of particles..."
            CALL FlipSign()
        ENDIF
        
        IF(TSinglePartPhase) THEN
!Exit the single particle phase if the number of walkers exceeds the value in the input file.
!            CALL MPI_Barrier(MPI_COMM_WORLD,error)
            IF(iProcIndex.eq.root) THEN     !Only exit phase if particle number is sufficient on head node.
                IF((TotParts.gt.InitWalkers).or.(ABS(AllNoatHF).gt.MaxNoatHF)) THEN
                    WRITE(6,*) "Exiting the single particle growth phase - shift can now change"
                    VaryShiftIter=Iter
                    TSinglePartPhase=.false.
                ENDIF
            ENDIF
!Broadcast the fact that TSinglePartPhase may have changed to all processors - unfortunatly, have to do this broadcast every iteration.
            CALL MPILBcast(TSinglePartPhase,1,root)
        ELSE
!Exit the single particle phase if the number of walkers exceeds the value in the input file.
!            CALL MPI_Barrier(MPI_COMM_WORLD,error)
            IF(iProcIndex.eq.root) THEN     !Only exit phase if particle number is sufficient on head node.
                IF((ABS(AllNoatHF).lt.(MaxNoatHF-HFPopThresh))) THEN
                    WRITE(6,'(A)') "No at HF has fallen too low - reentering the single particle growth phase - particle number may grow again."
!                    VaryShiftIter=Iter
                    TSinglePartPhase=.true.
                    tReZeroShift=.true.
                ENDIF
            ENDIF
!Broadcast the fact that TSinglePartPhase may have changed to all processors - unfortunatly, have to do this broadcast every iteration.
            CALL MPILBcast(TSinglePartPhase,1,root)
        ENDIF


!This first call will calculate the GrowRate for each processor, taking culling into account
!        WRITE(6,*) "Get Here"
!        CALL FLUSH(6)
        CALL UpdateDiagSftPar()

!Put a barrier here so all processes synchronise
        CALL MPIBarrier(error)

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

        IF(tTruncInitiator) THEN
            inpairInit(1)=NoAborted
            inpairInit(2)=NoAddedInitiators
            inpairInit(3)=NoInitDets
            inpairInit(4)=NoNonInitDets
            inpairInit(5)=NoInitWalk
            inpairInit(6)=NoNonInitWalk
            inpairInit(7)=NoDoubSpawns
            inpairInit(8)=NoExtraInitDoubs
 
            CALL MPI_Reduce(inpairInit,outpairInit,8,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)

            AllNoAborted=outpairInit(1)
            AllNoAddedInitiators=outpairInit(2)
            AllNoInitDets=outpairInit(3)
            AllNoNonInitDets=outpairInit(4)
            AllNoInitWalk=outpairInit(5)
            AllNoNonInitWalk=outpairInit(6)
            AllNoDoubSpawns=outpairInit(7)
            AllNoExtraInitDoubs=outpairInit(8)
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
            IF(tTruncInitiator) AllGrowRateAbort=(AllTotParts+REAL(AllNoAborted))/(AllTotPartsOld+REAL(AllNoAbortedOld))
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
                IF((Iter-VaryShiftIter).ge.NShiftEquilSteps) THEN
!                    WRITE(6,*) Iter-VaryShiftIter, NEquilSteps*StepsSft
                    IF((Iter-VaryShiftIter).eq.NShiftEquilSteps) WRITE(6,*) 'Beginning to average shift value.'
                    VaryShiftCycles=VaryShiftCycles+1
                    SumDiagSft=SumDiagSft+DiagSft
                    AvDiagSft=SumDiagSft/REAL(VaryShiftCycles,r2)
                ENDIF

                IF(tTruncInitiator) THEN
                    DiagSftAbort=DiagSftAbort-(log(AllGrowRateAbort)*SftDamp)/(Tau*(StepsSft+0.D0))
                    IF((Iter-VaryShiftIter).ge.NShiftEquilSteps) THEN
                        SumDiagSftAbort=SumDiagSftAbort+DiagSftAbort
                        AvDiagSftAbort=SumDiagSftAbort/REAL(VaryShiftCycles,r2)
                    ENDIF
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
        
        IF(tReZeroShift) THEN
            DiagSft=0.D0
            VaryShiftCycles=0
            SumDiagSft=0.D0
            AvDiagSft=0.D0
        ENDIF

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

            IF((.not.TSinglePartPhase).and.tInitShiftBlocking.and.(Iter.eq.(VaryShiftIter+IterShiftBlock))) THEN
                CALL InitShiftErrorBlocking(Iter)
                tInitShiftBlocking=.false.
                tShiftBlocking=.true.
            ENDIF

!Then we perform the blocking at the end of each update cycle.         
            IF(tErrorBlocking) CALL SumInErrorContrib(Iter,AllENumCyc,AllHFCyc)
            IF(tShiftBlocking.and.(Iter.ge.(VaryShiftIter+IterShiftBlock))) CALL SumInShiftErrorContrib(Iter,DiagSft)
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
        NoAborted=0
        NoAddedInitiators=0
        NoInitDets=0
        NoNonInitDets=0
        NoInitWalk=0
        NoNonInitWalk=0
        NoDoubSpawns=0
        NoExtraInitDoubs=0

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
        AllNoAbortedOld=AllNoAborted
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
        AllNoAborted=0
        AllNoAddedInitiators=0
        AllNoInitDets=0
        AllNoNonInitDets=0
        AllNoInitWalk=0
        AllNoNonInitWalk=0
        AllNoDoubSpawns=0
        AllNoExtraInitDoubs=0



        RETURN
    END SUBROUTINE CalcNewShift

    
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

!            IF(tDirectAnnihil) CALL Stop_All(this_routine,"READPOPS currently disabled with directannihilation")
            CALL ReadFromPopsFilePar()

        ELSE
!initialise the particle positions - start at HF with positive sign

!Set the maximum number of walkers allowed
            MaxWalkersPart=NINT(MemoryFacPart*InitWalkers)
!            WRITE(6,"(A,F14.5)") "Memory Factor for walkers is: ",MemoryFacPart
            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node of: ",MaxWalkersPart
            IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
!                WRITE(6,"(A,F14.5)") "Memory Factor for arrays used for spawning is: ",MemoryFacSpawn
!                WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for spawning of: ",MaxSpawned
            ELSE
                MaxWalkersAnnihil=NINT(MemoryFacAnnihil*InitWalkers)
!                WRITE(6,"(A,F14.5)") "Memory Factor for arrays used for annihilation is: ",MemoryFacAnnihil
                WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for annihilation of: ",MaxWalkersAnnihil
            ENDIF

!Put a barrier here so all processes synchronise
            CALL MPI_Barrier(MPI_COMM_WORLD,error)
!Allocate memory to hold walkers
            ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),4,this_routine,WalkVecDetsTag,ierr)
            WalkVecDets(0:NIfTot,1:MaxWalkersPart)=0
            IF((.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) THEN
!Rotoannihilation only used a single main array. Spawned particles are put into the spawned arrays.
                ALLOCATE(WalkVec2Dets(0:NIfTot,MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVec2Dets',MaxWalkersPart*(NIfTot+1),4,this_routine,WalkVec2DetsTag,ierr)
                WalkVec2Dets(0:NIfTot,1:MaxWalkersPart)=0
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
                MemoryAlloc=(NIfTot+1+3)*MaxWalkersPart*4    !Memory Allocated in bytes
            ELSE
!The sign is sent through when annihilating, so it needs to be longer.
                ALLOCATE(WalkVecSign(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('WalkVecSign',MaxWalkersAnnihil,4,this_routine,WalkVecSignTag,ierr)
                ALLOCATE(WalkVec2Sign(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('WalkVec2Sign',MaxWalkersAnnihil,4,this_routine,WalkVec2SignTag,ierr)
                MemoryAlloc=((2*MaxWalkersAnnihil)+(((2*(NIfTot+1))+4)*MaxWalkersPart))*4    !Memory Allocated in bytes
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
                ALLOCATE(SpawnVec(0:NIftot,MaxSpawned),stat=ierr)
                CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),4,this_routine,SpawnVecTag,ierr)
                SpawnVec(:,:)=0
                ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
                CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),4,this_routine,SpawnVec2Tag,ierr)
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

                MemoryAlloc=MemoryAlloc+(((MaxSpawned+1)*2)+(2*MaxSpawned*(1+NIfTot)))*4

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
                        AllNoAbortedOld=0
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
                        AllNoAbortedOld=0
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

        IF(tTruncInitiator.or.tDelayTruncInit) THEN
            IF(.not.tDirectAnnihil) THEN
                CALL Stop_All(this_routine,'TRUNCINITIATOR can only be used with direct annihilation.') 
            ELSEIF(tFixShiftShell) THEN
                CALL Stop_All(this_routine,'TRUNCINITIATOR cannot be used with the fixed shift shell approximation.') 
            ELSEIF(tTruncInitiator.and.(.not.tTruncCAS).and.(.not.tTruncSpace)) THEN
                CALL Stop_All(this_routine,'The initiator space has not been defined - need to use EXCITE X or TRUNCATECAS X Y.')
            ELSEIF(tDelayTruncInit) THEN
                tTruncInitiator=.false.
            ENDIF
        ENDIF

        IF(tPrintTriConnections.or.tHistTriConHEls.or.tPrintHElAccept) CALL InitTriHElStats()
        IF(tPrintSpinCoupHEl) CALL InitSpinCoupHEl()

        IF(tPrintOrbOcc) THEN
            ALLOCATE(OrbOccs(nBasis),stat=ierr)
            CALL LogMemAlloc('OrbOccs',nBasis,8,this_routine,OrbOccsTag,ierr)
            OrbOccs(:)=0.D0
        ENDIF

        IF(MaxNoatHF.eq.0) THEN
            MaxNoatHF=InitWalkers*nProcessors
            HFPopThresh=MaxNoatHF
        ENDIF

        IF((NMCyc.ne.0).and.(tRotateOrbs.and.(.not.tFindCINatOrbs))) CALL Stop_All(this_routine,"Currently not set up to rotate and then go straight into a spawning &
                                                                                    & calculation.  Ordering of orbitals is incorrect.  This may be fixed if needed.")

!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
        RETURN

    END SUBROUTINE InitFCIMCCalcPar


        



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

    SUBROUTINE CheckOrdering(DetArray,SignArray,NoDets,tCheckSignCoher)
        INTEGER :: NoDets,DetArray(0:NIfTot,1:NoDets),i,j
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
            Comp=DetBitLT(DetArray(:,i-1),DetArray(:,i),NIfDBO)
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

        ALLOCATE(AllowedDetList(0:NIfTot,NAllowedDetList),stat=ierr)
        IF(ierr.ne.0) THEN
            CALL Stop_All("ReadSpawnListDets","Error allocating AllowedDetList array")
        ENDIF

        IF(iProcIndex.eq.Root) THEN
            REWIND(17)

            do i=1,NAllowedDetList
                READ(17,*) nI(1:NEl)
                CALL EncodeBitDet(nI,AllowedDetList(0:NIfTot,i))
!                WRITE(6,*) AllowedDetList(0:NIfTot,i)
            enddo

            CLOSE(17)
        ENDIF

        CALL MPIIBCast(AllowedDetList,NAllowedDetList*(NIfTot+1),Root)
!        do i=1,NAllowedDetList
!            CALL DecodeBitDet(nI,AllowedDetList(0:NIfTot,i))
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
                        WRITE(17) CurrentDets(0:NIfTot,j),CurrentSign(j)
                    ENDIF
                enddo
            ELSE
                do j=1,TotWalkers
!First write out walkers on head node
                    IF(mod(j,iPopsPartEvery).eq.0) THEN
                        do k=0,NIfTot
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
            ALLOCATE(OrigParts(0:NIfTot,TotWalkers),stat=error)
            CALL LogMemAlloc('OrigParts',TotWalkers*(NIfTot+1),4,this_routine,OrigPartsTag,error)
            do i=1,TotWalkers
                OrigSign(i)=CurrentSign(i)
                OrigParts(:,i)=CurrentDets(:,i)
            enddo

!Now we need to receive the data from each other processor sequentially
!We can overwrite the head nodes information, since we have now stored it elsewhere.
            do i=1,nProcessors-1
!Run through all other processors...receive the data...
                CALL MPI_Recv(CurrentDets(0:NIfTot,1:WalkersonNodes(i)),WalkersonNodes(i)*(NIfTot+1),MPI_INTEGER,i,Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(CurrentSign(1:WalkersonNodes(i)),WalkersonNodes(i),MPI_INTEGER,i,Tag,MPI_COMM_WORLD,Stat,error)
!                WRITE(6,*) "Recieved walkers for processor ",i
!                CALL FLUSH(6)
                
!Then write it out...
                IF(tBinPops) THEN
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            WRITE(17) CurrentDets(0:NIfTot,j),CurrentSign(j)
                        ENDIF
                    enddo
                ELSE
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            do k=0,NIfTot
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
            CALL MPI_Send(CurrentDets(0:NIfTot,1:TotWalkers),TotWalkers*(NIfTot+1),MPI_INTEGER,root,Tag,MPI_COMM_WORLD,error)
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
                        WRITE(17) CurrentDets(0:NIfTot,j),CurrentSign(j)
                    ENDIF
                enddo
            ELSE
                do j=1,TotWalkers
!First write out walkers on head node
                    IF(mod(j,iPopsPartEvery).eq.0) THEN
                        do k=0,NIfTot
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
                CALL MPI_Recv(NewDets(0:NIfTot,1:WalkersonNodes(i)),WalkersonNodes(i)*(NIfTot+1),MPI_INTEGER,i,Tag,MPI_COMM_WORLD,Stat,error)
                CALL MPI_Recv(NewSign(1:WalkersonNodes(i)),WalkersonNodes(i),MPI_INTEGER,i,Tag,MPI_COMM_WORLD,Stat,error)
!                WRITE(6,*) "Recieved walkers for processor ",i
!                CALL FLUSH(6)
                
!Then write it out...
                IF(tBinPops) THEN
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            WRITE(17) NewDets(0:NIfTot,j),NewSign(j)
                        ENDIF
                    enddo
                ELSE
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            do k=0,NIfTot
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
            CALL MPI_Send(CurrentDets(0:NIfTot,1:TotWalkers),TotWalkers*(NIfTot+1),MPI_INTEGER,root,Tag,MPI_COMM_WORLD,error)
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
        LOGICAL :: exists,First,tBinRead
        INTEGER :: AvWalkers,WalkerstoReceive(nProcessors)
        INTEGER*8 :: NodeSumNoatHF(nProcessors),TempAllSumNoatHF
        REAL*8 :: TempTotParts
        INTEGER :: TempInitWalkers,error,i,j,k,l,total,ierr,MemoryAlloc,Tag,iLutTemp(0:NIfTot),TempSign,Proc,CurrWalkers
        INTEGER :: Stat(MPI_STATUS_SIZE),AvSumNoatHF,VecSlot,IntegerPart,HFPointer,TempnI(NEl),ExcitLevel,VecInd,DetsMerged
        REAL*8 :: r,FracPart,TempTotWalkers,Gap
        TYPE(HElement) :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        
        PreviousCycles=0    !Zero previous cycles
        SumENum=0.D0
        TotParts=0
        SumNoatHF=0
        DiagSft=0.D0
        Tag=124             !Set Tag
        
        IF(tDirectAnnihil) THEN
            INQUIRE(FILE='POPSFILE',EXIST=exists)
            IF(exists) THEN
                OPEN(17,FILE='POPSFILE',Status='old')
                tBinRead=.false.
            ELSE
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

            IF(DiagSft.eq.0.D0) tWalkContGrow=.true.

            IF(tBinRead) THEN
                CLOSE(17)
                OPEN(17,FILE='POPSFILEBIN',Status='old',form='unformatted')
            ENDIF

        ELSE
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

                IF(DiagSft.eq.0.D0) tWalkContGrow=.true.

                IF(tBinRead) THEN
                    CLOSE(17)
                    OPEN(17,FILE='POPSFILEBIN',Status='old',form='unformatted')
                ENDIF
            ENDIF
        ENDIF

        IF(iProcIndex.eq.Root) THEN

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
            IF(.not.tWalkContGrow) THEN
                InitWalkers=AvWalkers
            ELSE
                TSinglePartPhase=.true.
            ENDIF
            SumENum=AllSumENum/REAL(nProcessors,r2)     !Divide up the SumENum over all processors
            AvSumNoatHF=NINT(AllSumNoatHF/real(nProcessors,r2)) !This is the average Sumnoathf
            do i=1,nProcessors-1
                NodeSumNoatHF(i)=INT(AvSumNoatHF,i2)
            enddo
            NodeSumNoatHF(nProcessors)=NINT(AllSumNoatHF,i2)-INT((AvSumNoatHF*(nProcessors-1)),i2)

            ProjectionE=AllSumENum/AllSumNoatHF
                
!Reset the global variables
            AllSumENum=0.D0
            AllSumNoatHF=0.D0

        ENDIF

        CALL MPI_Barrier(MPI_COMM_WORLD,error)  !Sync

!Now we need to scatter the WalkerstoReceive to each node, and allocate the desired memory to each node...
!Broadcast info which needs to go to all processors
        CALL MPI_BCast(DiagSft,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,error)
        CALL MPI_BCast(InitWalkers,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
        CALL MPI_BCast(NEquilSteps,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
        CALL MPI_BCast(NShiftEquilSteps,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
        CALL MPI_BCast(SumENum,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,error)
        CALL MPI_BCast(TSinglePartPhase,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
!        CALL MPI_BCast(tChangenProcessors,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
!Scatter the number of walkers each node will receive to TempInitWalkers, and the SumNoatHF for each node which is distributed approximatly equally
        CALL MPI_Scatter(WalkerstoReceive,1,MPI_INTEGER,TempInitWalkers,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
        CALL MPI_Scatter(NodeSumNoatHF,1,MPI_DOUBLE_PRECISION,SumNoatHF,1,MPI_DOUBLE_PRECISION,root,MPI_COMM_WORLD,error)

        IF(tDirectAnnihil) THEN
            IF(ScaleWalkers.ne.1) CALL Stop_All(this_routine,'Sorry, direct annihilation cannot cope with scaling just yet.')
            IF(MemoryFacPart.le.1.D0) THEN
                WRITE(6,*) 'MemoryFacPart must be larger than 1.0 when reading in a POPSFILE - increasing it to 1.50.'
                MemoryFacPart=1.50
            ENDIF
        ENDIF
        
!Now we want to allocate memory on all nodes.
        MaxWalkersPart=NINT(MemoryFacPart*(NINT(InitWalkers*ScaleWalkers)))   !InitWalkers here is simply the average number of walkers per node, not actual
        IF(tRotoAnnihil.or.tDirectAnnihil) THEN
            MaxSpawned=NINT(MemoryFacSpawn*(NINT(InitWalkers*ScaleWalkers)))
        ELSE
            MaxWalkersAnnihil=NINT(MemoryFacAnnihil*(NINT(InitWalkers*ScaleWalkers)))   !InitWalkers here is simply the average number of walkers per node, not actual
        ENDIF

        Gap=REAL(MaxSpawned)/REAL(nProcessors)
        do i=0,nProcessors-1
            InitialSpawnedSlots(i)=NINT(Gap*i)+1
        enddo
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
        ValidSpawnedList(:)=InitialSpawnedSlots(:)

        CALL MPI_Barrier(MPI_COMM_WORLD,error)  !Sync

!Allocate memory to hold walkers at least temporarily
        ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating WalkVecDets array.')
        CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),4,this_routine,WalkVecDetsTag,ierr)
!        WalkVecDets(0:NIfTot,1:MaxWalkersPart)=0
        WalkVecDets(:,:)=0
        IF(tRotoAnnihil.or.tDirectAnnihil) THEN
            ALLOCATE(WalkVecSign(MaxWalkersPart),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating WalkVecSign array.')
            CALL LogMemAlloc('WalkVecSign',MaxWalkersPart,4,this_routine,WalkVecSignTag,ierr)
            WalkVecSign(:)=0
        ELSE
            ALLOCATE(WalkVecSign(MaxWalkersAnnihil),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating WalkVecSign array.')
            CALL LogMemAlloc('WalkVecSign',MaxWalkersAnnihil,4,this_routine,WalkVecSignTag,ierr)
        ENDIF

!Just allocating this here, so that the SpawnParts arrays can be used for sorting the determinants when using direct annihilation.
        IF(tRotoAnnihil.or.tDirectAnnihil) THEN

            WRITE(6,"(A,I12,A)") " Spawning vectors allowing for a total of ",MaxSpawned," particles to be spawned in any one iteration."
            ALLOCATE(SpawnVec(0:NIfTot,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),4,this_routine,SpawnVecTag,ierr)
            SpawnVec(:,:)=0
            ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),4,this_routine,SpawnVec2Tag,ierr)
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

            MemoryAlloc=MemoryAlloc+(((MaxSpawned+1)*2)+(2*MaxSpawned*(1+NIfTot)))*4
!            IF(tRotoAnnihil) MemoryAlloc=MemoryAlloc+(((MaxSpawned+1)*2)+(2*MaxSpawned*(1+NIfTot)))*4

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

!Allocate pointers to the correct walker arrays...
        CurrentDets=>WalkVecDets
        CurrentSign=>WalkVecSign
!        CurrentIC=>WalkVecIC

!Need to now allocate other arrays
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
                WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*4*2,r2)/1048576.D0," Mb/Processor"
            ELSE
                WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*4*4,r2)/1048576.D0," Mb/Processor"
            ENDIF
        ENDIF

        IF(.not.tRegenDiagHEls) THEN
            CurrentH=>WalkVecH
            NewH=>WalkVec2H
        ENDIF
        NewDets=>WalkVec2Dets
        NewSign=>WalkVec2Sign
!        NewIC=>WalkVec2IC

        IF(tDirectAnnihil) THEN

! The hashing will be different in the new calculation from the one where the POPSFILE was produced, this means we must recalculate the processor each determinant wants to go to.                
! This is done by reading in all walkers to the root and then distributing them in the same way as the spawning steps are done - by finding the determinant and sending it there.

!            IF(iProcIndex.eq.root) THEN

!                do i=1,AllTotWalkers
!                    iLutTemp(:)=0
!                    READ(17,*) iLutTemp(0:NIfTot),TempSign

                    !WRITE(6,'(A,3I20)') 'when dealing with directannihilation',iLutnJ(:)
!                    Proc=DetermineDetProc(iLutTemp)   !This wants to return a value between 0 -> nProcessors-1
                    !WRITE(6,*) iLutnJ(:),Proc,ValidSpawnedList(Proc),Child,TotWalkers
                    !CALL FLUSH(6)
!                    SpawnedParts(:,ValidSpawnedList(Proc))=iLutTemp(:)
!                    SpawnedSign(ValidSpawnedList(Proc))=TempSign
!                    IF(tTruncInitiator) THEN
!                        IF(SpawnedParts(NIfTot,ValidSpawnedList(Proc)).ne.0) CALL Stop_All(this_routine,'The parent initiator flags should all be 0 by the time             & 
!&                                                                                                        the POPSFILE is written, however for some reason this is not so.')
!                    ENDIF
!                    ValidSpawnedList(Proc)=ValidSpawnedList(Proc)+1
!                enddo
!                CLOSE(17)
! the spawnedparts and spawnedsign are then sent to the currentdets and currentsign arrays of the relevant processors in the same way as 
! they are in each iteration - however obviously no annihilation needs to be done.

!                IF(iProcIndex.ne.root) AllTotWalkers=0
!Walkers are only read in from the root processor - this is just making sure it is clear the other processers don't currently have any walkers to distribute.

!                CALL MPI_Barrier(MPI_COMM_WORLD,error)  !Sync

!                TempInitWalkers=AllTotWalkers

!                IF((CurrentDets(0,1).ne.0).or.(CurrentDets(NIfD,1).ne.0)) CALL Stop_All(this_routine,'Memory issue. Need to increase MemoryFacSpawn.')

!                CurrWalkers=0
! This is just a parameter to say that there are no walkers in the CurrentDets array yet - these all come from the POPSFILE as though they are all spawned at once.            
!                IF(tTruncInitiator) THEN
!                    WRITE(6,*) 'TRUNCINITIATOR on.  Turning it off temporarily.'
!                    CALL FLUSH(6)
!                    tTruncInitiator=.false.
!                    CALL DirectAnnihilation(CurrWalkers)
!                    tTruncInitiator=.true.
!                    WRITE(6,*) 'Turning TRUNCINITIATOR back on.'
!                ELSE
!                    CALL DirectAnnihilation(CurrWalkers)
!                ENDIF
!                CALL MPI_Barrier(MPI_COMM_WORLD,error)  !Sync
!            ENDIF


            CurrWalkers=0
            do i=1,AllTotWalkers
                iLutTemp(:)=0
                IF(tBinRead) THEN
                    READ(17) iLutTemp(0:NIfTot),TempSign
                ELSE
                    READ(17,*) iLutTemp(0:NIfTot),TempSign
                ENDIF
                Proc=DetermineDetProc(iLutTemp)   !This wants to return a value between 0 -> nProcessors-1
                IF(Proc.eq.iProcIndex) THEN
                    CurrWalkers=CurrWalkers+1
                    CurrentDets(0:NIfTot,CurrWalkers)=iLutTemp(0:NIfTot)
                    CurrentSign(CurrWalkers)=TempSign
                ENDIF
            enddo
            CLOSE(17)
            
!            DEALLOCATE(SpawnVec)
!            CALL LogMemDeAlloc(this_routine,SpawnVecTag)
!            DEALLOCATE(SpawnVec2)
!            CALL LogMemDeAlloc(this_routine,SpawnVec2Tag)
!            DEALLOCATE(SpawnSignVec)
!            CALL LogMemDeAlloc(this_routine,SpawnSignVecTag)
!            DEALLOCATE(SpawnSignVec2)
!            CALL LogMemDeAlloc(this_routine,SpawnSignVec2Tag)

!            MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
!            Gap=REAL(MaxSpawned)/REAL(nProcessors)
!            do i=0,nProcessors-1
!                InitialSpawnedSlots(i)=NINT(Gap*i)+1
!            enddo
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
!            ValidSpawnedList(:)=InitialSpawnedSlots(:)

!            ALLOCATE(SpawnVec(0:NIfTot,MaxSpawned),stat=ierr)
!            CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),4,this_routine,SpawnVecTag,ierr)
!            SpawnVec(:,:)=0
!            ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
!            CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),4,this_routine,SpawnVec2Tag,ierr)
!            SpawnVec2(:,:)=0
!            ALLOCATE(SpawnSignVec(0:MaxSpawned),stat=ierr)
!            CALL LogMemAlloc('SpawnSignVec',MaxSpawned+1,4,this_routine,SpawnSignVecTag,ierr)
!            SpawnSignVec(:)=0
!            ALLOCATE(SpawnSignVec2(0:MaxSpawned),stat=ierr)
!            CALL LogMemAlloc('SpawnSignVec2',MaxSpawned+1,4,this_routine,SpawnSignVec2Tag,ierr)
!            SpawnSignVec2(:)=0

!Point at correct spawning arrays
!            SpawnedParts=>SpawnVec
!            SpawnedParts2=>SpawnVec2
!            SpawnedSign=>SpawnSignVec
!            SpawnedSign2=>SpawnSignVec2

!            MemoryAlloc=MemoryAlloc+(((MaxSpawned+1)*2)+(2*MaxSpawned*(1+NIfTot)))*4

        ELSE

            IF(iProcIndex.eq.root) THEN

!Root process reads all walkers in and then sends them to the correct processor
                do i=nProcessors,1,-1
!Read in data for processor i
                    IF(tBinRead) THEN
                        do j=1,WalkerstoReceive(i)
                            READ(17) WalkVecDets(0:NIfTot,j),WalkVecSign(j)
                        enddo
                    ELSE
                        do j=1,WalkerstoReceive(i)
                            READ(17,*) WalkVecDets(0:NIfTot,j),WalkVecSign(j)
                        enddo
                    ENDIF

                    IF(i.ne.1) THEN
!Now send data to processor i-1 (Processor rank goes from 0 -> nProcs-1). If i=1, then we want the data so stay at the root processor.
                        CALL MPI_Send(WalkVecDets(:,1:WalkerstoReceive(i)),WalkerstoReceive(i)*(NIfTot+1),MPI_INTEGER,i-1,Tag,MPI_COMM_WORLD,error)
                        CALL MPI_Send(WalkVecSign(1:WalkerstoReceive(i)),WalkerstoReceive(i),MPI_INTEGER,i-1,Tag,MPI_COMM_WORLD,error)
                    ENDIF

                enddo

                CLOSE(17)

                IF(tRotoAnnihil) THEN
                    WRITE(6,*) "Ordering/compressing all walkers that have been read in..."
                    CALL SortBitDets(WalkerstoReceive(1), &
                                     WalkVecDets(:,1:WalkerstoReceive(1)), &
                                     WalkVecSign(1:WalkerstoReceive(1)))
                    
                    VecInd=1
                    DetsMerged=0
                    do i=2,TempInitWalkers
                        IF(.not.DetBitEQ(WalkVecDets(0:NIfTot,i),WalkVecDets(0:NIfTot,VecInd),NIfDBO)) THEN
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

            CALL MPI_Barrier(MPI_COMM_WORLD,error)  !Sync
            

            do i=1,nProcessors-1
                IF(iProcIndex.eq.i) THEN
!All other processors want to pick up their data from root
                    CALL MPI_Recv(WalkVecDets(:,1:TempInitWalkers),TempInitWalkers*(NIfTot+1),MPI_INTEGER,0,Tag,MPI_COMM_WORLD,Stat,error)
                    CALL MPI_Recv(WalkVecSign(1:TempInitWalkers),TempInitWalkers,MPI_INTEGER,0,Tag,MPI_COMM_WORLD,Stat,error)
                    IF(tRotoAnnihil) THEN
                        CALL SortBitDets(TempInitWalkers, &
                                         WalkVecDets(:,1:TempInitWalkers), &
                                         WalkVecSign(1:TempInitWalkers))
                        VecInd=1
                        DetsMerged=0
                        do l=2,TempInitWalkers
                            IF(.not.DetBitEQ(WalkVecDets(0:NIfTot,l),WalkVecDets(0:NIfTot,VecInd),NIfDBO)) THEN
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

        ENDIF


        IF(iProcIndex.eq.root) WRITE(6,'(I10,A)') INT(AllTotWalkers,i2)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers  by a factor of: ",ScaleWalkers

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
                ALLOCATE(WalkVec2Dets(0:NIfTot,MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVec2Dets',MaxWalkersPart*(NIfTot+1),4,this_routine,WalkVec2DetsTag,ierr)
                WalkVec2Dets(0:NIfTot,1:MaxWalkersPart)=0
                ALLOCATE(WalkVec2Sign(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('WalkVec2Sign',MaxWalkersAnnihil,4,this_routine,WalkVec2SignTag,ierr)
                WalkVec2Sign(:)=0

!Scale up the integer part and fractional part seperately
                IntegerPart=INT(ScaleWalkers)       !Round to zero
                FracPart=ScaleWalkers-REAL(IntegerPart)

                VecSlot=1
                do l=1,TempInitWalkers
                    do k=1,IntegerPart
                        WalkVec2Dets(0:NIfTot,VecSlot)=WalkVecDets(0:NIfTot,l)
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
                        WalkVec2Dets(0:NIfTot,VecSlot)=WalkVecDets(0:NIfTot,l)
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
!This is no longer necessary when the scaled size is calculated at the very beginning.
!                DEALLOCATE(WalkVecDets)
!                CALL LogMemDealloc(this_routine,WalkVecDetsTag)
!                DEALLOCATE(WalkVecSign)
!                CALL LogMemDealloc(this_routine,WalkVecSignTag)
!                ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
!                CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),4,this_routine,WalkVecDetsTag,ierr)
                WalkVecDets(0:NIfTot,1:MaxWalkersPart)=0
!                IF(tRotoAnnihil) THEN
!                    ALLOCATE(WalkVecSign(MaxWalkersPart),stat=ierr)
!                    CALL LogMemAlloc('WalkVecSign',MaxWalkersPart,4,this_routine,WalkVecSignTag,ierr)
!                ELSE
!                    ALLOCATE(WalkVecSign(MaxWalkersAnnihil),stat=ierr)
!                    CALL LogMemAlloc('WalkVecSign',MaxWalkersAnnihil,4,this_routine,WalkVecSignTag,ierr)
!                ENDIF
                WalkVecSign(:)=0

!Transfer scaled particles back accross to WalkVecDets
                do l=1,TotWalkers
                    WalkVecDets(0:NIfTot,l)=WalkVec2Dets(0:NIfTot,l)
                    WalkVecSign(l)=WalkVec2Sign(l)
                enddo

!Zero the second array for good measure
                WalkVec2Dets(0:NIfTot,1:MaxWalkersPart)=0

            ENDIF

        ELSE
!We are not scaling the number of walkers...

            IF((.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) THEN
                ALLOCATE(WalkVec2Dets(0:NIfTot,MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVec2Dets',MaxWalkersPart*(NIfTot+1),4,this_routine,WalkVec2DetsTag,ierr)
                WalkVec2Dets(0:NIfTot,1:MaxWalkersPart)=0
                ALLOCATE(WalkVec2Sign(MaxWalkersAnnihil),stat=ierr)
                CALL LogMemAlloc('WalkVec2Sign',MaxWalkersAnnihil,4,this_routine,WalkVec2SignTag,ierr)
                WalkVec2Sign(:)=0
            ENDIF
                
            IF(tDirectAnnihil) THEN
                IF(iProcIndex.eq.root) THEN
                    AllTotWalkers=TotWalkers
                    AllTotWalkersOld=AllTotWalkers
                    WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(AllTotWalkers,i2)
                ENDIF
                TotWalkers=CurrWalkers
                TotWalkersOld=CurrWalkers
            ELSE

                TotWalkers=TempInitWalkers      !Set the total number of walkers
                TotWalkersOld=TempInitWalkers
                IF(iProcIndex.eq.root) THEN
                    AllTotWalkersOld=AllTotWalkers
                    WRITE(6,*) "Total number of initial walkers is now: ",INT(AllTotWalkers,i2)
                ENDIF
            ENDIF

        ENDIF
            
        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft

        IF(tRotoAnnihil.or.tDirectAnnihil) THEN
            MemoryAlloc=((NIfTot+1+3)*MaxWalkersPart*4)
        ELSE
            MemoryAlloc=((2*MaxWalkersAnnihil)+(((2*(NIfTot+1))+4)*MaxWalkersPart))*4    !Memory Allocated in bytes
        ENDIF
        IF(tRegenDiagHEls) THEN
            IF(tRotoAnnihil.or.tDirectAnnihil) THEN
                MemoryAlloc=MemoryAlloc-(MaxWalkersPart*4*2)
            ELSE
                MemoryAlloc=MemoryAlloc-(MaxWalkersPart*4*4)
            ENDIF
        ENDIF

        WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,r2)/1048576.D0," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        CALL FLUSH(6)


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
            IF(.not.tDirectAnnihil) THEN
                WRITE(6,"(A,F14.6,A)") "Temp Arrays for annihilation cannot be more than : ",REAL(MaxWalkersPart*12,r2)/1048576.D0," Mb/Processor"
            ENDIF
        ENDIF
        CALL FLUSH(6)

!Now find out the data needed for the particles which have been read in...
        First=.true.
        TotParts=0
        do j=1,TotWalkers
            CALL DecodeBitDet(TempnI,CurrentDets(:,j))
            Excitlevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j), 2)
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

                IF((.not.TNoAnnihil).and.(.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) THEN
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

                IF((.not.TNoAnnihil).and.(.not.tRotoAnnihil).and.(.not.tDirectAnnihil)) THEN
!                IF((.not.TNoAnnihil).and.(.not.TAnnihilonproc)) THEN
                    HashArray(j)=CreateHash(TempnI)
                ENDIF
            ENDIF
            TotParts=TotParts+abs(CurrentSign(j))

        enddo

        TempTotParts=REAL(TotParts,r2)

        CALL MPI_Barrier(MPI_COMM_WORLD,error)  !Sync
        CALL MPI_Reduce(TempTotParts,AllTotParts,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)

        IF(iProcIndex.eq.root) AllTotPartsOld=AllTotParts

        WRITE(6,'(A,F20.1)') ' The total number of particles read from the POPSFILE is: ',AllTotParts

        RETURN

    END SUBROUTINE ReadFromPopsfilePar

!This will set up the initial walker distribution proportially to the MP1 wavevector.
    SUBROUTINE InitWalkersMP1Par()
        use SystemData , only : tAssumeSizeExcitgen
        use symexcit3 , only : GenExcitations3,CountExcitations3 
        INTEGER :: HFConn,error,ierr,MemoryAlloc,VecSlot,nJ(NEl),nStore(6),iExcit,i,j,WalkersonHF,HFPointer,ExcitLevel,VecInd
        REAL*8 :: SumMP1Compts,MP2Energy,Compt,r,FracPart,TempTotWalkers,TempTotParts
        TYPE(HElement) :: Hij,Hjj,Fjj
        INTEGER , ALLOCATABLE :: MP1Dets(:,:), ExcitgenTemp(:)
        INTEGER , ALLOCATABLE :: MP1Sign(:)
        REAL*8 , ALLOCATABLE :: MP1Comps(:),MP1CompsNonCum(:)
        INTEGER :: MP1DetsTag,MP1SignTag,MP1CompsTag,SumWalkersonHF,ExcitLength,iMaxExcit,IntParts,MP1CompsNonCumTag
        CHARACTER(len=*), PARAMETER :: this_routine='InitWalkersMP1Par'
        LOGICAL :: First,TurnBackAssumeExGen
        INTEGER :: ExcitMat3(2,2),nSingles,nDoubles
        LOGICAL :: tParity,tAllExcitFound

    
        IF(tHub.and.tReal) THEN
            CALL Stop_All(this_routine,"Cannot initialise walkers in MP1 for real space hubbard calculations.")
        ENDIF
!Set the maximum number of walkers allowed
        MaxWalkersPart=NINT(MemoryFacPart*InitWalkers)
!        WRITE(6,"(A,F14.5)") "Memory Factor for walkers is: ",MemoryFacPart
        WRITE(6,"(A,I14)") "Memory allocated for a maximum particle/det number per node of: ",MaxWalkersPart
        IF(tRotoAnnihil) THEN
            MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
!            WRITE(6,"(A,F14.5)") "Memory Factor for arrays used for spawning is: ",MemoryFacSpawn
            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for spawning of: ",MaxSpawned
        ELSE
            MaxWalkersAnnihil=NINT(MemoryFacAnnihil*InitWalkers)
!            WRITE(6,"(A,F14.5)") "Memory Factor for arrays used for annihilation is: ",MemoryFacAnnihil
            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for annihilation of: ",MaxWalkersAnnihil
        ENDIF

!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
!Allocate memory to hold walkers
        ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
        CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),4,this_routine,WalkVecDetsTag,ierr)
        WalkVecDets(0:NIfTot,1:MaxWalkersPart)=0
        IF(tRotoAnnihil) THEN
            ALLOCATE(WalkVecSign(MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkersPart,4,this_routine,WalkVecSignTag,ierr)
        ELSE
            ALLOCATE(WalkVec2Dets(0:NIfTot,MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVec2Dets',MaxWalkersPart*(NIfTot+1),4,this_routine,WalkVec2DetsTag,ierr)
            WalkVec2Dets(0:NIfTot,1:MaxWalkersPart)=0
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
            MemoryAlloc=(NIfTot+1+3)*MaxWalkersPart*4    !Memory Allocated in bytes
            IF(tRegenDiagHEls) MemoryAlloc=MemoryAlloc-(MaxWalkersPart*8)
        ELSE
            MemoryAlloc=((2*MaxWalkersAnnihil)+(((2*(NIfTot+1))+4)*MaxWalkersPart))*4    !Memory Allocated in bytes
            IF(tRegenDiagHEls) MemoryAlloc=MemoryAlloc-(MaxWalkersPart*16)
        ENDIF


        IF(tRotoAnnihil) THEN
            
            WRITE(6,"(A,I12,A)") " Spawning vectors allowing for a total of ",MaxSpawned," particles to be spawned in any one iteration."
            ALLOCATE(SpawnVec(0:NIfTot,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),4,this_routine,SpawnVecTag,ierr)
            SpawnVec(:,:)=0
            ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
            CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),4,this_routine,SpawnVec2Tag,ierr)
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

            MemoryAlloc=MemoryAlloc+(((MaxSpawned+1)*2)+(2*MaxSpawned*(1+NIfTot)))*4

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
        IF(tNoSpinSymExcitgens) THEN
            exflag=3
            CALL CountExcitations3(HFDet,exflag,nSingles,nDoubles)
            HFConn=nSingles+nDoubles+1
        ELSE
            CALL GetSymExcitCount(HFExcit%ExcitData,HFConn)
            HFConn=HFConn+1     !Add on one for the HF Det itself
        ENDIF


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
        IF(tNoSpinSymExcitgens) THEN

            tAllExcitFound=.false.
            ExcitMat3(:,:)=0
! exflag of 2 means only generate double excitations from the HF determinant.
            exflag=2

            do while (.not.tAllExcitFound)
                CALL GenExcitations3(HFDet,iLutHF,nJ,exflag,ExcitMat3,tParity,tAllExcitFound)
                IF(tAllExcitFound) EXIT

                Hij=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,2,ECore)
                CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fjj)
!                WRITE(6,"(4I5,2G25.10)") nJ(:),real(Hij%v,r2),(Fii-(REAL(Fjj%v,r2)))

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

        ELSE

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
                IF(IsNullDet(nJ)) EXIT
                IF(iExcit.ne.2) THEN
                    CALL Stop_All("InitWalkersMP1","Error - excitations other than doubles being generated in MP1 wavevector code")
                ENDIF

                Hij=GetHElement2(HFDet,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,iExcit,ECore)
                CALL GetH0Element(nJ,NEl,Arr,nBasis,ECore,Fjj)
!                WRITE(6,"(8I5,2G25.10)") nJ(:),real(Hij%v,r2),(Fii-(REAL(Fjj%v,r2)))

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

        ENDIF


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
                        CurrentDets(0:NIfTot,VecInd)=iLutHF(:)
                        CurrentSign(VecInd)=IntParts*MP1Sign(j)
                        TotParts=TotParts+IntParts
                        IF(.not.tRegenDiagHEls) THEN
                            CurrentH(VecInd)=0.D0
                        ENDIF
                    ELSE
!We are at a double excitation - we need to calculate most of this information...
                        CALL EncodeBitDet(MP1Dets(1:NEl,j),CurrentDets(0:NIfTot,VecInd))
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
                CALL SortBitDets(TotWalkers,CurrentDets(:,1:TotWalkers), &
                                 CurrentSign(1:TotWalkers))
            ELSE
                CALL SortBitDetswH(TotWalkers,CurrentDets(:,1:TotWalkers), &
                         CurrentSign(1:TotWalkers),CurrentH(1:TotWalkers))
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
                    CurrentDets(0:NIfTot,j)=iLutHF(:)
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
                    CALL EncodeBitDet(MP1Dets(1:NEl,i),CurrentDets(0:NIfTot,j))
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
                ExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j), 2)
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


!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
    INTEGER FUNCTION AttemptCreatePar(DetCurr,iLutCurr,WSign,nJ,iLutnJ,Prob,IC,Ex,tParity,nParts,tMinorDetList)
        use GenRandSymExcitNUMod , only : GenRandSymExcitBiased
        use Logging, only : CCMCDebug
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,StoreNumTo,StoreNumFrom,DetLT,i,ExtraCreate,Ex(2,2),WSign,nParts
        INTEGER :: iLutCurr(0:NIfTot),Bin,iLutnJ(0:NIfTot),PartInd,ExcitLev,iLut(0:NIfTot),iLut2(0:NIfTot)
        LOGICAL :: tParity,SymAllowed,tSuccess,tMinorDetList
        integer :: yama(NIfY)
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
            CALL FindExcitBitDet(iLutCurr,iLutnJ,IC,Ex)
            ExcitLev = FindBitExcitLevel(iLutHF, iLutnJ, iMaxDomLev)
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

            ! TODO: for csfs, have ilutnJ by now. Useful to eval. energy.
            rh=GetHElement4(DetCurr,nJ,IC,Ex,tParity)
            !WRITE(6,*) rh%v

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

!        IF(Iter.eq.463) THEN
!            WRITE(6,"(10I5)") DetCurr(:)
!            WRITE(6,"(10I5,2F15.5,I5)") nJ(:),REAL(rh%v,8),Prob,AttemptCreatePar
!        ENDIF

        
!We know we want to create a particle. Return the bit-representation of this particle (if we have not already got it)
        IF(.not.tHPHF.and.AttemptCreatePar.ne.0) THEN
            if (tCSF) then
                ! This makes sure that the Yamanouchi symbol is correct. It
                ! also makes it work if we have tTruncateCSF on, and ex would
                ! therefore leave all singles as beta, when we switch to dets.
                call EncodeBitDet (nJ, iLutnJ)
            else
                call FindExcitBitDet(iLutCurr,iLutnJ,IC,Ex,yama)
            endif
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
        INTEGER :: iLutCurr(0:NIfTot),Bin,iLutParent(0:NIfTot),PartInd,ExcitLev,IC
        LOGICAL :: SymAllowed,tSuccess,tMinorDetList
        REAL*8 :: Prob=1.D0,r,rat
        TYPE(HElement) :: rh,rhcheck

        IF(tSpawnDominant) THEN
!We only allow spawning at determinants between iMinDomLev and iMaxDomLev which are in an allowed list of dominant determinants.
!Find the excitation level of the excitation
!We need to calculate the bit-representation of this new child. This can be done easily since the ExcitMat is known.
            tMinorDetList=.false.
            ExcitLev = FindBitExcitLevel(iLutHF, iLutParent, iMaxDomLev)
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
        INTEGER :: nJ(NEl),ierr,IC,Child,iCount,DetCurr(NEl),iLutnJ(0:NIfTot)
        REAL*8 :: Prob,rat,HDiag,HDiagCurr,r,HSum
        INTEGER :: iDie,WalkExcitLevel,iMaxExcit,ExcitLength,PartInd,iExcit
        INTEGER :: ExcitLevel,TotWalkersNew,iGetExcitLevel_2,error,length,temp,Ex(2,2),WSign,p,Scratch1(ScratchSize),Scratch2(ScratchSize)
        LOGICAL :: tParity,tMainArr,tFilled,tSuccess,tMinorDetList
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
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j))

!Also, we want to find out the excitation level - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.
            IF(tTruncSpace) THEN
!We need to know the exact excitation level for truncated calculations.
                WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j),&
                                                   nel)
            ELSE
                WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j),2)
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
                    IF(.not.tHPHF) CALL FindExcitBitDet(CurrentDets(:,j),iLutnJ,IC,Ex)

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
                IF(IsNullDet(nJ)) EXIT

!Find matrix element
                HElemTemp=GetHElement3(DetCurr,nJ,iExcit)
                IF((abs(REAL(HElemTemp%v,r2))).gt.1.D-8) THEN

!Encode this determinant
                    CALL EncodeBitDet(nJ,iLutnJ)

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
                CALL EncodeBitDet(DetCurr,iLutnJ)

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
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j))

!Also, we want to find out the excitation level - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
            WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j), 2)
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
        INTEGER :: iLutCurr(0:NIfTot),DetCurr(NEl),i,nStore(6),ierr,iMaxExcit
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
        iLutHF(0:NIfTot)=iLutCurr(0:NIfTot)
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
                
                CALL MPI_Send(CurrentDets(:,IndexFrom:TotWalkers),WalktoTransfer(1)*(NIfTot+1),MPI_INTEGER,WalktoTransfer(3),Tag,MPI_COMM_WORLD,error)
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

                CALL MPI_Recv(CurrentDets(:,IndexFrom:IndexTo),WalktoTransfer(1)*(NIfTot+1),MPI_INTEGER,WalktoTransfer(2),Tag,MPI_COMM_WORLD,Stat,error)
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
                CALL SortCompressListswH(TotWalkers,CurrentDets(0:NIfTot,1:TotWalkers),CurrentSign(1:TotWalkers),CurrentH(1:TotWalkers))
            ELSE
                CALL SortCompressLists(TotWalkers,CurrentDets(0:NIfTot,1:TotWalkers),CurrentSign(1:TotWalkers))
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
        INTEGER :: Length,PartList(0:NIfTot,Length),SignList(Length)
        INTEGER :: i,VecInd,DetsMerged

        CALL SortBitDets(Length,PartList,SignList)

!Now compress the list.
        VecInd=1
        DetsMerged=0
        TotParts=0
        IF(Length.gt.0) THEN
            TotParts=TotParts+abs(SignList(1))
        ENDIF
        do i=2,Length
            TotParts=TotParts+abs(SignList(i))
            IF(.not.DetBitEQ(PartList(0:NIfTot,i),PartList(0:NIfTot,VecInd),NIfDBO)) THEN
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
!        DetCurr=PartList(:,1)
!        i=2
!        do while(i.le.Length)
!            CompParts=DetBitLT(DetCurr,PartList(0:NIfTot,i))
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
!                DetCurr=PartList(:,i)
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
        INTEGER :: Length,PartList(0:NIfTot,Length),SignList(Length),j
        REAL*8 :: HList(Length)
        INTEGER :: i,DetsMerged,VecInd

        CALL SortBitDetswH(Length,PartList,SignList,HList)
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
            IF(.not.DetBitEQ(PartList(0:NIfTot,i),PartList(0:NIfTot,VecInd),NIfDBO)) THEN
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
!        DetCurr=PartList(:,1)
!        i=2
!        do while(i.le.Length)
!            CompParts=DetBitLT(DetCurr,PartList(0:NIfTot,i))
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
!                DetCurr=PartList(:,i)
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
            IF(iProcIndex.eq.Root) THEN
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
        IF(tHistHamil) THEN
            DEALLOCATE(HistHamil)
            DEALLOCATE(AvHistHamil)
            IF(iProcIndex.eq.0) THEN
                DEALLOCATE(AllHistHamil)
                DEALLOCATE(AllAvHistHamil)
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

        IF(.not.tNoSpinSymExcitgens) DEALLOCATE(HFExcit%ExcitData)
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


    SUBROUTINE WriteHistogramEnergies()
        INTEGER :: error,i
        REAL*8 :: Norm,EnergyBin

        IF(iProcIndex.eq.Root) THEN
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
        CALL MPI_Reduce(SinglesHist,AllSinglesHist,iOffDiagNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesAttemptHist,AllSinglesAttemptHist,iOffDiagNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(DoublesHist,AllDoublesHist,iOffDiagNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(DoublesAttemptHist,AllDoublesAttemptHist,iOffDiagNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesHistOccOcc,AllSinglesHistOccOcc,iOffDiagNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesHistOccVirt,AllSinglesHistOccVirt,iOffDiagNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesHistVirtOcc,AllSinglesHistVirtOcc,iOffDiagNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SinglesHistVirtVirt,AllSinglesHistVirtVirt,iOffDiagNoBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
  

        IF(iProcIndex.eq.Root) THEN
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
                IF(AllHistogram(i).gt.0.D0) WRITE(17,*) EnergyBin, AllHistogram(i)
                IF(AllAttemptHist(i).gt.0.D0) WRITE(18,*) EnergyBin, AllAttemptHist(i)
                IF(AllSpawnHist(i).gt.0.D0) WRITE(19,*) EnergyBin, AllSpawnHist(i)
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
                IF(AllSinglesHist(i).gt.0.D0) WRITE(20,*) EnergyBin, AllSinglesHist(i)
                IF(AllSinglesAttemptHist(i).gt.0.D0) WRITE(21,*) EnergyBin, AllSinglesAttemptHist(i)
                IF(AllDoublesHist(i).gt.0.D0) WRITE(22,*) EnergyBin, AllDoublesHist(i)
                IF(AllDoublesAttemptHist(i).gt.0.D0) WRITE(23,*) EnergyBin, AllDoublesAttemptHist(i)
                IF(AllSinglesHistOccOcc(i).gt.0.D0) WRITE(24,*) EnergyBin, AllSinglesHistOccOcc(i)
                IF(AllSinglesHistOccVirt(i).gt.0.D0) WRITE(25,*) EnergyBin, AllSinglesHistOccVirt(i)
                IF(AllSinglesHistVirtOcc(i).gt.0.D0) WRITE(26,*) EnergyBin, AllSinglesHistVirtOcc(i)
                IF(AllSinglesHistVirtVirt(i).gt.0.D0) WRITE(27,*) EnergyBin, AllSinglesHistVirtVirt(i)
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
        LOGICAL :: DetsEq,tSuccess


        iGuideDets=0
        AlliInitGuideParts=0
        DetsEq=.false.
        IF(iProcIndex.eq.Root) THEN
            OPEN(36,FILE='GUIDINGFUNC',Status='old')
            READ(36,*) iGuideDets 
        ENDIF
 
        CALL MPI_Bcast(iGuideDets,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)

        ALLOCATE(GuideFuncDets(0:NIfTot,1:iGuideDets),stat=ierr)
        CALL LogMemAlloc('GuideFuncDets',(NIfTot+1)*iGuideDets,4,this_routine,GuideFuncDetsTag,ierr)
        ALLOCATE(GuideFuncSign(0:iGuideDets),stat=ierr)
        CALL LogMemAlloc('GuideFuncSign',iGuideDets+1,4,this_routine,GuideFuncSignTag,ierr)

        ALLOCATE(DetstoRotate(0:NIfTot,1:iGuideDets),stat=ierr)
        CALL LogMemAlloc('DetstoRotate',(NIfTot+1)*iGuideDets,4,this_routine,DetstoRotateTag,ierr)
        ALLOCATE(SigntoRotate(0:iGuideDets),stat=ierr)
        CALL LogMemAlloc('SigntoRotate',iGuideDets+1,4,this_routine,SigntoRotateTag,ierr)
        ALLOCATE(DetstoRotate2(0:NIfTot,1:iGuideDets),stat=ierr)
        CALL LogMemAlloc('DetstoRotate2',(NIfTot+1)*iGuideDets,4,this_routine,DetstoRotate2Tag,ierr)
        ALLOCATE(SigntoRotate2(0:iGuideDets),stat=ierr)
        CALL LogMemAlloc('SigntoRotate2',iGuideDets+1,4,this_routine,SigntoRotate2Tag,ierr)


        IF(iProcIndex.eq.Root) THEN
            !Set up the determinant and sign arrays by reading in from the GUIDINGFUNC file.
            j=1
            do while (j.le.iGuideDets)
                READ(36,*) GuideFuncDets(0:NIfTot,j),GuideFuncSign(j)
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
!                WRITE(37,*) GuideFuncDets(0:NIfTot,j),nProcessors*GuideFuncSign(j)
            enddo
!            CLOSE(37)
            
        ENDIF

        !Broadcast the guiding function determinants (and signs) to all processors.
        !The total number of walkers in the guiding function is therefore nProcessors*iInitGuideParts.
        CALL MPI_Bcast(GuideFuncDets(0:NIfTot,1:iGuideDets),iGuideDets*(NIfTot+1),MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL MPI_Bcast(GuideFuncSign(0:iGuideDets),iGuideDets+1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)

        !Run through the guiding function determinants and find the index that contains the HF.
        !Want this known on all processors, so that we can just look up the sign at this position to get the guiding function HF population.
!        do i=1,iGuideDets
!            DetsEq=DetBitEQ(iLutHF,GuideFuncDets(0:NIfTot,i))
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
                ExcitLevel = FindBitExcitLevel(GuideFuncDets(:,i), iLutHF, 2)
                IF(ExcitLevel.eq.2) THEN
                    DoubDet(:)=0
                    CALL DecodeBitDet(DoubDet,GuideFuncDets(0:NIfTot,i))
                    HdoubTemp=GetHElement2(HFDet,DoubDet,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
                    HDoub=REAL(HDoubTemp%v,r2)
                    GuideFuncDoub=GuideFuncDoub+(GuideFuncSign(i)*Hdoub)
                ENDIF
            enddo
            WRITE(6,*) 'The energy of the guiding function alone is ,',GuideFuncDoub/(REAL(GuideFuncSign(GuideFuncHFIndex),r2))
            GuideFuncDoub=0.D0
        ENDIF
            

    ENDSUBROUTINE InitGuidingFunction




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

        CALL FLUSH(6)
        IF(iGuideDets.gt.TotWalkers) CALL Stop_All(this_routine,'iGuideDets is greater than the number of populated determinants')

!        WRITE(6,*) 'the determinants and sign, on each processor, before I touched them'
!        do j=1,TotWalkers
!            WRITE(6,*) CurrentDets(0:NIfTot,j),CurrentSign(j)
!        enddo
        
! Firstly order CurrentSign in descending absolute value, taking the corresponding CurrentDets with it.
        CALL SortBitSign(TotWalkers,CurrentSign(1:TotWalkers), &
                         CurrentDets(:,1:TotWalkers))

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
!            WRITE(6,*) CurrentDets(0:NIfTot,j),CurrentSign(j)
!        enddo
!        CALL FLUSH(6)
        
! Now take the iGuideDets determinants and reorder them back in terms of determinants (taking the sign with them).
        
        CALL SortBitDets(iGuideDets,CurrentDets(:,1:iGuideDets), &
                         CurrentSign(1:iGuideDets))

!        WRITE(6,*) 'The most populated determinants on each processor, ordered by determinant'
!        WRITE(6,*) iGuideDets,' determinants included in the guiding function.'    
!        do j=1,iGuideDets
!            WRITE(6,*) CurrentDets(0:NIfTot,j),CurrentSign(j)
!        enddo
        
! Calculate RecvCounts(1:nProcessors), and OffSets for the Gatherv calculation
! Need to gather the iGuideDets values from each processor for this.
        CALL MPI_Gather(iGuideDets,1,MPI_INTEGER,RecvCounts,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        Offsets(:)=0
        do j=1,nProcessors-1
            OffSets(j+1)=RecvCounts(j)+OffSets(j)
        enddo
        do j=1,nProcessors
            RecvCounts02(j)=RecvCounts(j)*(NIfTot+1)
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
            ALLOCATE(AllCurrentDets(0:NIfTot,1:AlliGuideDets),stat=ierr)
            CALL LogMemAlloc('AllCurrentDets',(NIfTot+1)*AlliGuideDets,4,this_routine,AllCurrentDetsTag,ierr)
            ALLOCATE(AllCurrentSign(1:AlliGuideDets),stat=ierr)
            CALL LogMemAlloc('AllCurrentSign',AlliGuideDets,4,this_routine,AllCurrentSignTag,ierr)
        ENDIF

        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        CALL MPI_Gatherv(CurrentSign(1:iGuideDets),iGuideDets,MPI_INTEGER,AllCurrentSign(1:AlliGuideDets),&
        &RecvCounts,Offsets,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL MPI_Gatherv(CurrentDets(0:NIfTot,1:iGuideDets),((NIfTot+1)*iGuideDets),MPI_INTEGER,&
        &AllCurrentDets(0:NIfTot,1:AlliGuideDets),RecvCounts02,Offsets02,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL FLUSH(6)


        IF(iProcIndex.eq.Root) THEN

!            OPEN(31,file='GUIDINGFUNCall-01',status='unknown')
!            WRITE(31,*) 'The determinants from each processor combined'
!            WRITE(31,*) AlliGuideDets,' determinants included in the guiding function.'    
!            do j=1,AlliGuideDets
!                WRITE(31,*) AllCurrentDets(0:NIfTot,j),AllCurrentSign(j)
!            enddo
!            CLOSE(31)

! Having taken the largest occupied determinants from each processor, select the iGuideDets most populated from
! this total list.
! I.e. reorder in terms of population, compressing so that each determinant only appears once.  Then select out the top iGuideDets 
! and reorder by determinant once again.

! From the list of determinants from all processors, order in terms of determinant so that this may be compressed to make sure each
! determinant only appears once.
            CALL SortBitDets(AlliGuideDets,AllCurrentDets(:,1:AlliGuideDets),&
                             AllCurrentSign(1:AlliGuideDets))

!            OPEN(32,file='GUIDINGFUNCall-02',status='unknown')
!            WRITE(32,*) 'The determinants from each processor combined'
!            WRITE(32,*) AlliGuideDets,' determinants included in the guiding function.'    
!            do j=1,AlliGuideDets
!                WRITE(32,*) AllCurrentDets(0:NIfTot,j),AllCurrentSign(j)
!            enddo
!            CLOSE(32)
            

! Find out if two determinants next to each are the same, if so add their signs and move the ones below them up.
! DetBitEq returns true if two determinants are identical, or false otherwise.
            CompiGuideDets=AlliGuideDets
            do i=1,AlliGuideDets-1

                IF(i.gt.CompiGuideDets) EXIT 
                ! This means we have got to the end of the compressed list, don't want to keep going, as all determinants after this are 0. 

                do while (DetBitEQ(AllCurrentDets(0:NIfTot,i),AllCurrentDets(0:NIfTot,i+1),NIfDBO))
                    ! Take a determinant, if the one above it is identical, add its sign to the original and move the others up to overwrite
                    ! the second.
                    ! Repeat this until the i+1 determinant is no longer equal to the i determinant.

                    IF((AllCurrentSign(i)*AllCurrentSign(i+1)).lt.0) THEN
                        WRITE(6,*) 'Determinant populated with opposite signs',AllCurrentDets(0:NIfTot,i)
                        WRITE(6,*) 'Identical determinant next to it',AllCurrentDets(0:NIfTot,i+1)
                        CALL FLUSH(6)
                        CALL Stop_All("WriteGuidingFunc","Identical determinants populated with opposite sign")
                    ENDIF

                    AllCurrentSign(i)=AllCurrentSign(i)+AllCurrentSign(i+1)

                    CompiGuideDets=CompiGuideDets-1

                    do j=i+2,AlliGuideDets
                        AllCurrentDets(0:NIfTot,j-1)=AllCurrentDets(0:NIfTot,j)
                        AllCurrentSign(j-1)=AllCurrentSign(j)
                    enddo
                    ! Zero the last determinant
                    AllCurrentSign(AlliGuideDets)=0
                    AllCurrentDets(0:NIfTot,AlliGuideDets)=0
                enddo
            enddo
            AlliGuideDets=CompiGuideDets
     
!            OPEN(33,file='GUIDINGFUNCall-03',status='unknown')
!            WRITE(33,*) 'The determinants from each processor combined and compressed, ordered by determinant'
!            WRITE(33,*) AlliGuideDets,' determinants included in the guiding function.'    
!            do j=1,AlliGuideDets
!                WRITE(33,*) AllCurrentDets(0:NIfTot,j),AllCurrentSign(j)
!            enddo
!            CLOSE(33)

          
! Reorder the compressed list of determinants by population (i.e. descending according to absolute value of AllCurrentSign, taking determinant with it).        
            CALL SortBitSign(AlliGuideDets,AllCurrentSign(1:AlliGuideDets), &
                             AllCurrentDets(:,1:AlliGuideDets))

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
!                WRITE(34,*) AllCurrentDets(0:NIfTot,j),AllCurrentSign(j)
!            enddo
!            CLOSE(34)


! Now take the iGuideDets determinants and reorder them back in terms of determinants (taking the sign with them).
            CALL SortBitDets(iGuideDets,AllCurrentDets(:,1:iGuideDets), &
                             AllCurrentSign(1:iGuideDets))


! Write the iGuideDets most populated determinants (in order of their bit strings) to a file.
            OPEN(35,file='GUIDINGFUNC',status='unknown')
            WRITE(35,*) iGuideDets,' determinants included in the guiding function.'    
            do j=1,iGuideDets
                WRITE(35,*) AllCurrentDets(0:NIfTot,j),AllCurrentSign(j)
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
                WRITE(38,*) GuideFuncDets(0:NIfTot,j),InitGuideFuncSign(j),AllGuideFuncSign(j),(ABS(AllGuideFuncSign(j))-ABS(InitGuideFuncSign(j)))
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


        ALLOCATE(DomDets(0:NIfTot,1:iNoDomDets),stat=ierr)
        CALL LogMemAlloc('DomDets',(NIfTot+1)*iNoDomDets,4,this_routine,DomDetsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory for dominant determinants')

        IF(iProcIndex.eq.Root) THEN
            !Set up the determinant and sign arrays by reading in from the GUIDINGFUNC file.
            j=1
            do while (j.le.iNoDomDets)
                READ(41,*) DomDets(0:NIfTot,j)
                j=j+1
            enddo
            CLOSE(41)
        ENDIF

        ! Broadcast this list of DomDets to all processors
        CALL MPI_Bcast(DomDets(0:NIfTot,1:iNoDomDets),(NIfTot+1)*iNoDomDets,MPI_INTEGER,Root,MPI_COMM_WORLD,ierr)


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
        INTEGER :: SpinCoupDetBit(0:NIfTot),SpinCoupDet(NEl),OpenShell(2,NEl),UpSpin(NEl),NoOpenShell,NoUpSpin,iRead,PartInd,ID1,ID2,iComb,TempSign
        LOGICAL :: tDoubOcc,tSuccess

        CALL FLUSH(6)
        IF(.not.tRotoAnnihil) CALL Stop_All(this_routine,'PRINTDOMINANTDETS can only be used with rotoannihilation.')

        WRITE(6,'(A13,I10,A53,I3,A5,I3)') 'Printing the ',iNoDominantDets,' dominant determinants with excitation level between ',MinExcDom,' and ',MaxExcDom


!        WRITE(6,*) 'the determinants and sign, on each processor, before I touched them'
!        do j=1,TotWalkers
!            WRITE(6,*) CurrentDets(0:NIfTot,j),CurrentSign(j)
!        enddo


! Firstly, copy the determinants with excitation level between the min and max into a separate array (and do the same with sign).  
! Need to first count the number that are going to go in this array.
        NoExcDets=0
        AllNoExcDets=0
        do i=1,TotWalkers
            ExcitLevel = FindBitExcitLevel(CurrentDets(:,i), iLutHF, &
                                           MaxExcDom)
            IF((ExcitLevel.ge.MinExcDom).and.(ExcitLevel.le.MaxExcDom)) THEN
                NoExcDets=NoExcDets+1
            ELSEIF(ExcitLevel.eq.0) THEN
                HFPop=ABS(CurrentSign(i))
            ENDIF
        enddo

        CALL MPI_Reduce(NoExcDets,AllNoExcDets,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)


! Allocate arrays of this size - these are the ones that will be reordered to find the iNoDominantDets most populated etc.        
        ALLOCATE(ExcDets(0:NIfTot,1:NoExcDets),stat=ierr)
        CALL LogMemAlloc('ExcDets',(NIfTot+1)*NoExcDets,4,this_routine,ExcDetsTag,ierr)
        ALLOCATE(ExcSign(1:NoExcDets),stat=ierr)
        CALL LogMemAlloc('ExcSign',NoExcDets,4,this_routine,ExcSignTag,ierr)
 
        IF(iProcIndex.eq.Root) THEN
            ALLOCATE(AllExcDets(0:NIfTot,1:(10*AllNoExcDets)),stat=ierr)
            CALL LogMemAlloc('AllExcDets',(NIfTot+1)*10*AllNoExcDets,4,this_routine,AllExcDetsTag,ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'ERROR allocating memory to AllExcDets.')
            ALLOCATE(AllExcSign(1:(10*AllNoExcDets)),stat=ierr)
            CALL LogMemAlloc('AllExcSign',10*AllNoExcDets,4,this_routine,AllExcSignTag,ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'ERROR allocating memory to AllExcSign.')
        ENDIF


! Now run through the occupied determinants.  If the determinant has the correct excitation level, add it to the ExcDets array, and
! add its sign to the ExcSign array.
        ExcDetsIndex=0
        do i=1,TotWalkers
            ExcitLevel = FindBitExcitLevel(CurrentDets(:,i), iLutHF(:), &
                                           MaxExcDom)
            IF((ExcitLevel.ge.MinExcDom).and.(ExcitLevel.le.MaxExcDom)) THEN
                ExcDetsIndex=ExcDetsIndex+1
                ExcDets(0:NIfTot,ExcDetsIndex)=CurrentDets(0:NIfTot,i)
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
            RecvCounts02(j)=RecvCounts(j)*(NIfTot+1)
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
        CALL MPI_Gatherv(ExcDets(0:NIfTot,1:NoExcDets),((NIfTot+1)*NoExcDets),MPI_INTEGER,&
        &AllExcDets(0:NIfTot,1:AllNoExcDets),RecvCounts02,Offsets02,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        CALL FLUSH(6)
       
! Now that we have arrays on the root processor with the determinants and sign of correct excitation level, need to order these 
! in descending absolute value, taking the corresponding sign with it.
        IF(iProcIndex.eq.Root) THEN
 
            IF(iNoDominantDets.gt.AllNoExcDets) THEN
                WRITE(6,*) 'iNoDominantDets: ',iNoDominantDets
                WRITE(6,*) 'AllNoExcDets: ',AllNoExcDets
                CALL Stop_All(this_routine,'Not enough determinants are occupied to pick out the number of dominant requested.')
            ENDIF


            CALL SortBitSign(AllNoExcDets,AllExcSign(1:AllNoExcDets), &
                             AllExcDets(:,1:AllNoExcDets))

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
                WRITE(39,*) AllExcDets(0:NIfTot,j),AllExcSign(j)
            enddo
            CLOSE(39)
          
            WRITE(6,*) 'This amounts to determinants with population ',MinDomDetPop,' and larger.'
!            WRITE(6,*) 'HFPop',HFPop
!            WRITE(6,*) 'MinDomDetPop',MinDomDetPop
            MinRelDomPop=REAL(MinDomDetPop)/REAL(HFPop)
            WRITE(6,*) 'These determinants have amplitude ; ',MinRelDomPop,' relative to the most populated determinant.' 

! In order to do binary searches for the spin determinants, need to sort the determinants back into order.
! Do this in two separate lots, 1:iNoDominantDets and iNoDominantDets+1:AllNoExcDets

            CALL SortBitDets(iNoDominantDets,AllExcDets(:,1:iNoDominantDets),&
                             AllExcSign(1:iNoDominantDets))
            CALL SortBitDets((AllNoExcDets-iNoDominantDets), &
                             AllExcDets(:,(iNoDominantDets+1):AllNoExcDets), &
                             AllExcSign((iNoDominantDets+1):AllNoExcDets))
 
            OPEN(47,file='DOMINANTDETSsorted',status='unknown')
            WRITE(47,*) AllNoExcDets,' determinants with the right excitation level.'    
            WRITE(47,*) iNoDominantDets,' with population ',MinDomDetPop, ' and above.'
            do j=1,AllNoExcDets
                WRITE(47,*) AllExcDets(0:NIfTot,j),AllExcSign(j)
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
                    CALL DecodeBitDet(DetCurr,AllExcDets(0:NIfTot,j))
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
!                        ID1 = GTID(DetCurr(i))
                        ID1=CEILING(REAL(DetCurr(i))/2.0)
                        ! Now run through all other orbitals finding out if the same spat orb is occupied.
                        do k=1,NEl
                            IF(k.eq.i) CYCLE 
                            tDoubOcc=.false.
!                            ID2 = GTID(DetCurr(k))
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

                        CALL EncodeBitDet(SpinCoupDet(1:NEl),SpinCoupDetBit(0:NIfTot))
!                        WRITE(6,*) 'SpinCoupDetBit',SpinCoupDetBit(:)

                        ! First search through the list of dominant determinants.
                        tSuccess=.false.
                        CALL BinSearchDomParts(AllExcDets(0:NIfTot,1:CurriNoDominantDets),SpinCoupDetBit(0:NIfTot),1,CurriNoDominantDets,PartInd,tSuccess)
                        IF(tSuccess) THEN
                            ! Determinant found in the dominant list.
                            iRead=iRead+1
                        ELSE
                            ! If not found in the original dominant determinant list then need to search through the determinants that have been added.
                            IF((iNoDominantDets-CurriNoDominantDets).gt.0) THEN
                                CALL BinSearchDomParts(AllExcDets(0:NIfTot,(CurriNoDominantDets+1):iNoDominantDets),SpinCoupDetBit(0:NIfTot),(CurriNoDominantDets+1),&
                                                    &iNoDominantDets,PartInd,tSuccess)
                                IF((PartInd.le.CurriNoDominantDets).or.(PartInd.gt.iNoDominantDets)) CALL Stop_All(this_routine, '')                                                    
                            ENDIF

                            IF(tSuccess) THEN
                                iRead=iRead+1
                                ! If the determinant has already been added, don't need to do anything.
                            ELSE

                                ! If there are still determinants left in the rest of the list, search through these to check if the determinant is there (to get the sign).
                                IF((AllNoExcDets-iNoDominantDets).gt.0) THEN
                                    CALL BinSearchDomParts(AllExcDets(0:NIfTot,(iNoDominantDets+1):AllNoExcDets),SpinCoupDetBit(0:NIfTot),(iNoDominantDets+1),&
                                                    &AllNoExcDets,PartInd,tSuccess)
                                ENDIF
                                IF(tSuccess) THEN
                                    ! Determinant found in the rest of the list, put it in dom dets along with its sign.
                                    ! Need to insert the determinant in the correct position, so that these added determinants stay ordered.
                                    ! SearchGen will give back ipos so that the determinant we are inserting is < ipos-1 and ge ipos.
                                    ! I.e this determinant goes in ipos and everything else is moved up 1.
                                    IF((iNoDominantDets-CurriNoDominantDets).gt.0) THEN
                                        CALL SearchGen((iNoDominantDets-CurriNoDominantDets),AllExcDets(0:NIfTot,(CurriNoDominantDets+1):iNoDominantDets),&
                                                    &SpinCoupDetBit(0:NIfTot),ipos)
                                    ELSE
                                        ipos=0
                                    ENDIF
                                    ipos=ipos+CurriNoDominantDets
!                                    CALL SearchGen(iNoDominantDets,AllExcDets(0:NIfTot,1:iNoDominantDets),SpinCoupDetBit(0:NIfTot),ipos)
                                    ! The position the determinant should go is ipos, if the current determinant in ipos is not equal to the SpinCoupDetBit,
                                    ! then we want to insert this determinant here, move all the others up by one, put the determinant that is getting written
                                    ! over at iNoDominantDets+1 in the position of the spincoupdetbit.
                                    
                                    TempSign=AllExcSign(iNoDominantDets+1) 
                                    AllExcDets(0:NIfTot,PartInd)=AllExcDets(0:NIfTot,iNoDominantDets+1)

                                    do i=iNoDominantDets,ipos,-1                                                    
                                        AllExcDets(0:NIfTot,i+1)=AllExcDets(0:NIfTot,i)                                                    
                                        AllExcSign(i+1)=AllExcSign(i)
                                    enddo
                                    AllExcDets(0:NIfTot,ipos)=SpinCoupDetBit(0:NIfTot)
                                    AllExcSign(ipos)=AllExcSign(PartInd)
                                    iRead=iRead+1
                                    iNoDominantDets=iNoDominantDets+1
                                    AllExcSign(PartInd)=TempSign
                                ELSE
                                    ! Determinant not in the list at all, add to added determinants (maintaining order) with a sign of 0.
                                    CALL SearchGen((iNoDominantDets-CurriNoDominantDets),AllExcDets(0:NIfTot,(CurriNoDominantDets+1):iNoDominantDets),&
                                                    &SpinCoupDetBit(0:NIfTot),ipos)
                                    ipos=ipos+CurriNoDominantDets
                                    
                                    AllExcDets(0:NIfTot,AllNoExcDets+1)=AllExcDets(0:NIfTot,iNoDominantDets+1)
                                    AllExcSign(AllNoExcDets+1)=AllExcSign(iNoDominantDets+1)

                                    do i=iNoDominantDets,ipos,-1                                                    
                                        AllExcDets(0:NIfTot,i+1)=AllExcDets(0:NIfTot,i)                                                    
                                        AllExcSign(i+1)=AllExcSign(i)
                                    enddo

                                    AllExcDets(0:NIfTot,ipos)=SpinCoupDetBit(0:NIfTot)
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

            CALL SortExcitBitDets(iNoDominantDets,AllExcDets(:,1:iNoDominantDets),AllExcSign(1:iNoDominantDets),iLutHF)
 
!            OPEN(40,file='DOMINANTDETSexclevelbit',status='unknown')
!            WRITE(40,*) AllNoExcDets,' determinants with the right excitation level.'    
!            do j=1,AllNoExcDets
!                ExcitLevel = FindBitExcitLevel(AllExcDets(:,j), iLutHF, &
!                                               MaxExcDom)
!                WRITE(40,*) AllExcDets(0:NIfTot,j),AllExcSign(j),ExcitLevel
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
                CurrExcitLevel = FindBitExcitLevel(AllExcDets(:,i), iLutHF, &
                                                   MaxExcDom)
                ExcitLevel=CurrExcitLevel
                do while (ExcitLevel.eq.CurrExcitLevel)
                    NoExcitLevel=NoExcitLevel+1
                    j=j+1
                    IF(j.gt.iNoDominantDets) EXIT
                    ExcitLevel = FindBitExcitLevel(AllExcDets(:,j), iLutHF, &
                                                   MaxExcDom)
                enddo
                WRITE(38,*) CurrExcitLevel,NoExcitLevel
            enddo
            IF(CurrExcitLevel.ne.MaxExcDom) THEN
                do i=CurrExcitLevel+1,MaxExcDom
                    WRITE(38,*) i,0
                enddo
            ENDIF

            do j=1,iNoDominantDets
                WRITE(38,*) AllExcDets(0:NIfTot,j),AllExcSign(j)
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
        ALLOCATE(MinorStarDets(0:NIfTot,1:MaxWalkersPart),stat=ierr)
        CALL LogMemAlloc('MinorStarDets',(NIfTot+1)*MaxWalkersPart,4,this_routine,MinorStarDetsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorStarDets')
        ALLOCATE(MinorSpawnDets(0:NIfTot,1:MaxSpawned),stat=ierr)
        CALL LogMemAlloc('MinorSpawnDets',(NIfTot+1)*MaxSpawned,4,this_routine,MinorSpawnDetsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorSpawnDets')
        MinorSpawnDets(:,:)=0
        ALLOCATE(MinorSpawnDets2(0:NIfTot,1:MaxSpawned),stat=ierr)
        CALL LogMemAlloc('MinorSpawnDets2',(NIftot+1)*MaxSpawned,4,this_routine,MinorSpawnDets2Tag,ierr)
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
        ALLOCATE(MinorStarParent(0:NIfTot,1:MaxWalkersPart),stat=ierr)
        CALL LogMemAlloc('MinorStarParent',(NIfTot+1)*MaxWalkersPart,4,this_routine,MinorStarParentTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorStarParent')
        ALLOCATE(MinorSpawnParent(0:NIfTot,1:MaxSpawned),stat=ierr)
        CALL LogMemAlloc('MinorSpawnParent',(NIfTot+1)*MaxSpawned,4,this_routine,MinorSpawnParentTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating memory to MinorSpawnParent')
        ALLOCATE(MinorSpawnParent2(0:NIfTot,1:MaxSpawned),stat=ierr)
        CALL LogMemAlloc('MinorSpawnParent2',(NIfTot+1)*MaxSpawned,4,this_routine,MinorSpawnParent2Tag,ierr)
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
                CALL SortRNI(Det,AllHistogram,FCIDets,NIfTot+1)
                
                OPEN(17,FILE='FCIMCPsi',STATUS='UNKNOWN')

                norm=0.D0
                do i=1,Det
                    norm=norm+AllHistogram(i)**2
!write out FCIMC Component weight (normalised), current normalisation, excitation level
                    ExcitLevel = FindBitExcitLevel(iLutHF, FCIDets(:,i), nel)
                    CALL DecodeBitDet(nI,FCIDets(0:NIfTot,i))
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
        INTEGER :: i,j,bits,iLut(0:NIfTot),error,IterRead
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

!        IF(NIfTot.ne.0) THEN
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
!                        CALL EncodeBitDet(NMRKS(:,j),iLut)
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

!This routine will write out the average hamiltonian from the spawning run up until now.
    SUBROUTINE WriteHamilHistogram()
        INTEGER :: i,j,error
        CHARACTER(len=22) :: abstr

!This will open a file called HamilHist-"Iter" on unit number 17.
        abstr=''
        write(abstr,'(I12)') Iter
        abstr='HamilHist-'//adjustl(abstr)
        IF(iProcIndex.eq.0) THEN
            WRITE(6,*) "Writing out the average hamiltonian up to iteration number: ", Iter
            CALL FLUSH(6)
        ENDIF

        IF(iProcIndex.eq.0) THEN
            AllHistHamil(:,:)=0.D0
            AllAvHistHamil(:,:)=0.D0
        ENDIF

        CALL MPI_Reduce(HistHamil,AllHistHamil,Det**2,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(AvHistHamil,AllAvHistHamil,Det**2,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        
        IF(iProcIndex.eq.0) THEN
!How do we normalise this!
            OPEN(17,FILE=abstr,STATUS='UNKNOWN')
            do i=1,Det
                do j=1,Det
                    WRITE(17,*) j,i,AllAvHistHamil(j,i),AllHistHamil(j,i)
                enddo
                WRITE(17,*) ""
            enddo
            CLOSE(17)
        ENDIF
        HistHamil(:,:)=0.D0

    END SUBROUTINE WriteHamilHistogram


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
!        INTEGER :: TotWalkersNew,j,k,l,DetCurr(0:NIfTot),VecSlot,TotWalkersDet
!        INTEGER :: DetLT
!
!        TempExcit%PointToExcit=>null()
!!First, it is necessary to sort the list of determinants
!        CALL SortPartsPar(TotWalkersNew,NewDets(:,1:TotWalkersNew),NIfTot+1)
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
!            do while ((DetLT(NewDets(:,j),DetCurr,(NIfTot+1)).eq.0).and.(j.le.TotWalkersNew))
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

    SUBROUTINE AddHistHamilEl(iLutnI,iLutnJ,WalkExcitLevel,Child,iTypeMatEl)
        INTEGER :: iTypeMatEl,WalkExcitLevel,iLutnI(0:NIfD),iLutnJ(0:NIfD),Child
        LOGICAL :: tSuccess
        INTEGER :: PartInd,PartIndChild,ChildExcitLevel

        IF(iTypeMatEl.eq.1) THEN
        !This is a spawning event
        !Need to histogram the hamiltonian - find the correct indicies of parent and child.
            IF(WalkExcitLevel.eq.NEl) THEN
                CALL BinSearchParts2(iLutnI,FCIDetIndex(WalkExcitLevel),Det,PartInd,tSuccess)
            ELSEIF(WalkExcitLevel.eq.0) THEN
                PartInd=1
                tSuccess=.true.
            ELSE
                CALL BinSearchParts2(iLutnI,FCIDetIndex(WalkExcitLevel),FCIDetIndex(WalkExcitLevel+1)-1,PartInd,tSuccess)
            ENDIF
            IF(.not.tSuccess) THEN
                CALL Stop_All("AddHistHamil","Cannot find determinant nI in list")
            ENDIF
            ChildExcitLevel = FindBitExcitLevel(iLutHF, iLutnJ, nel)
            IF(ChildExcitLevel.eq.NEl) THEN
                CALL BinSearchParts2(iLutnJ,FCIDetIndex(ChildExcitLevel),Det,PartIndChild,tSuccess)
            ELSEIF(ChildExcitLevel.eq.0) THEN
                PartIndChild=1
                tSuccess=.true.
            ELSE
                CALL BinSearchParts2(iLutnJ,FCIDetIndex(ChildExcitLevel),FCIDetIndex(ChildExcitLevel+1)-1,PartIndChild,tSuccess)
            ENDIF
            IF(.not.tSuccess) THEN
                CALL Stop_All("AddHistHamil","Cannot find determinant nJ in list")
            ENDIF
            
            HistHamil(PartIndChild,PartInd)=HistHamil(PartIndChild,PartInd)+(1.D0*Child)
            HistHamil(PartInd,PartIndChild)=HistHamil(PartInd,PartIndChild)+(1.D0*Child)
            AvHistHamil(PartIndChild,PartInd)=AvHistHamil(PartIndChild,PartInd)+(1.D0*Child)
            AvHistHamil(PartInd,PartIndChild)=AvHistHamil(PartInd,PartIndChild)+(1.D0*Child)
        ENDIF

    END SUBROUTINE AddHistHamilEl




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
!Very crudely hacked from the parallel version.  MPI calls commented out with !!


!Every StepsSft steps, update the diagonal shift value (the running value for the correlation energy)
!We don't want to do this too often, since we want the population levels to acclimatise between changing the shifts
    SUBROUTINE CalcNewShift()
        USE FciMCLoggingMOD , only : PrintSpawnAttemptStats,PrintTriConnStats,PrintSpinCoupHEl,InitErrorBlocking,SumInErrorContrib
        USE FciMCLoggingMOD , only : InitShiftErrorBlocking,SumInShiftErrorContrib
        INTEGER :: error,rc,MaxAllowedWalkers,MaxWalkersProc,MinWalkersProc
        INTEGER :: inpair(9),outpair(9),inpairInit(8),outpairInit(8)
        REAL*8 :: TempTotWalkers,TempTotParts
        REAL*8 :: TempSumNoatHF,MeanWalkers,TempSumWalkersCyc,TempAllSumWalkersCyc,TempNoMinorWalkers
        REAL*8 :: inpairreal(3),outpairreal(3)
        LOGICAL :: TBalanceNodesTemp

        TotImagTime=TotImagTime+StepsSft*Tau
        
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
            CALL MPILBcast(TSinglePartPhase,1,root)
        ENDIF

!This first call will calculate the GrowRate for each processor, taking culling into account
!        WRITE(6,*) "Get Here"
!        CALL FLUSH(6)
        CALL UpdateDiagSftPar()

!Put a barrier here so all processes synchronise
!!        CALL MPI_Barrier(MPI_COMM_WORLD,error)
        CALL MPIBarrier(error)

!Find the total number of particles at HF (x sign) across all nodes. If this is negative, flip the sign of all particles.
        AllNoatHF=0

!Find sum of noathf, and then use an AllReduce to broadcast it to all nodes
!!        CALL MPI_AllReduce(NoatHF,AllNoatHF,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,error)
        CALL MPIISum(NoatHF,1,AllNoatHF)
      

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
!        CALL MPI_Reduce(NoatHF,AllNoatHF,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)  !This is done every iteration now
!        CALL MPI_Reduce(NoatDoubs,AllNoatDoubs,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(NoBorn,AllNoBorn,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(NoDied,AllNoDied,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(SumWalkersCyc,AllSumWalkersCyc,1,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get Here 1"
!        CALL FLUSH(6)
!!        CALL MPI_Reduce(inpair,outpair,9,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
         call MPIISumArr(inpair,9,outpair)
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

        IF(tTruncInitiator) THEN
            inpairInit(1)=NoAborted
            inpairInit(2)=NoAddedInitiators
            inpairInit(3)=NoInitDets
            inpairInit(4)=NoNonInitDets
            inpairInit(5)=NoInitWalk
            inpairInit(6)=NoNonInitWalk
            inpairInit(7)=NoDoubSpawns
            inpairInit(8)=NoExtraInitDoubs
 
!!            CALL MPI_Reduce(inpairInit,outpairInit,8,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
            Call MPIISumArr(inpairInit,8,outpairInit)

            AllNoAborted=outpairInit(1)
            AllNoAddedInitiators=outpairInit(2)
            AllNoInitDets=outpairInit(3)
            AllNoNonInitDets=outpairInit(4)
            AllNoInitWalk=outpairInit(5)
            AllNoNonInitWalk=outpairInit(6)
            AllNoDoubSpawns=outpairInit(7)
            AllNoExtraInitDoubs=outpairInit(8)
        ENDIF

!!        CALL MPI_Reduce(TempTotWalkers,AllTotWalkers,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
!!        CALL MPI_Reduce(TempTotParts,AllTotParts,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
!!        CALL MPI_Reduce(TempNoMinorWalkers,AllNoMinorWalkers,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        Call MPIDSum(TempTotWalkers,1,AllTotWalkers)
        Call MPIDSum(TempTotParts,1,AllTotParts)
        Call MPIDSum(TempNoMinorWalkers,1,AllNoMinorWalkers)
      

        IF(iProcIndex.eq.0) THEN
            IF(AllTotWalkers.le.0.2) THEN
                WRITE(6,*) AllTotWalkers,TotWalkers
                CALL Stop_All("CalcNewShift","All particles have died. Consider choosing new seed, or raising shift value.")
            ENDIF
        ENDIF

!SumWalkersCyc is now an int*8, therefore is needs to be reduced as a real*8
        TempSumWalkersCyc=REAL(SumWalkersCyc,r2)
        TempAllSumWalkersCyc=0.D0
!!        CALL MPI_Reduce(TempSumWalkersCyc,TempAllSumWalkersCyc,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        Call MPIDSum(TempSumWalkersCyc,1,TempAllSumWalkersCyc)
!        WRITE(6,*) "Get Here 3"
!        CALL FLUSH(6)

        IF(.not.tDirectAnnihil) THEN

            MeanWalkers=AllTotWalkers/REAL(nProcessors,r2)
            MaxAllowedWalkers=NINT((MeanWalkers/12.D0)+MeanWalkers)

!Find the range of walkers on different nodes to see if we need to even up the distribution over nodes
!            inpair(1)=TotWalkers
!            inpair(2)=iProcIndex
!!            CALL MPI_Reduce(TotWalkers,MaxWalkersProc,1,MPI_INTEGER,MPI_MAX,root,MPI_COMM_WORLD,error)
!!            CALL MPI_Reduce(TotWalkers,MinWalkersProc,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,error)
            CALL MPIIReduce(TotWalkers,1,MPI_MAX,MaxWalkersProc)
            CALL MPIIReduce(TotWalkers,1,MPI_MIN,MinWalkersProc)

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
!!            CALL MPI_AllReduce(TBalanceNodesTemp,TBalanceNodes,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,error)
            CALL MPIAllReduceLORScal(TBalanceNodesTemp,TBalanceNodes,error)
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
!!            CALL MPI_Reduce(TotWalkers,MaxWalkersProc,1,MPI_INTEGER,MPI_MAX,root,MPI_COMM_WORLD,error)
!!            CALL MPI_Reduce(TotWalkers,MinWalkersProc,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,error)
            CALL MPIIReduce(TotWalkers,1,MPI_MAX,MaxWalkersProc)
            CALL MPIIReduce(TotWalkers,1,MPI_MIN,MinWalkersProc)
            
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
            IF(tTruncInitiator) AllGrowRateAbort=(AllTotParts+REAL(AllNoAborted))/(AllTotPartsOld+REAL(AllNoAbortedOld))
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
        CALL MPIDSumArr(inpairreal,3,outpairreal)
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
                IF((Iter-VaryShiftIter).ge.NShiftEquilSteps) THEN
!                    WRITE(6,*) Iter-VaryShiftIter, NEquilSteps*StepsSft
                    IF((Iter-VaryShiftIter).eq.NShiftEquilSteps) WRITE(6,*) 'Beginning to average shift value.'
                    VaryShiftCycles=VaryShiftCycles+1
                    SumDiagSft=SumDiagSft+DiagSft
                    AvDiagSft=SumDiagSft/REAL(VaryShiftCycles,r2)
                ENDIF

                IF(tTruncInitiator) THEN
                    DiagSftAbort=DiagSftAbort-(log(AllGrowRateAbort)*SftDamp)/(Tau*(StepsSft+0.D0))
                    IF((Iter-VaryShiftIter).ge.NShiftEquilSteps) THEN
                        SumDiagSftAbort=SumDiagSftAbort+DiagSftAbort
                        AvDiagSftAbort=SumDiagSftAbort/REAL(VaryShiftCycles,r2)
                    ENDIF
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
        CALL MPIDBcast(DiagSft,1,Root)
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

            IF((.not.TSinglePartPhase).and.tInitShiftBlocking.and.(Iter.eq.(VaryShiftIter+IterShiftBlock))) THEN
                CALL InitShiftErrorBlocking(Iter)
                tInitShiftBlocking=.false.
                tShiftBlocking=.true.
            ENDIF

!Then we perform the blocking at the end of each update cycle.         
            IF(tErrorBlocking) CALL SumInErrorContrib(Iter,AllENumCyc,AllHFCyc)
            IF(tShiftBlocking.and.(Iter.ge.(VaryShiftIter+IterShiftBlock))) CALL SumInShiftErrorContrib(Iter,DiagSft)
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
        NoAborted=0
        NoAddedInitiators=0
        NoInitDets=0
        NoNonInitDets=0
        NoInitWalk=0
        NoNonInitWalk=0
        NoDoubSpawns=0
        NoExtraInitDoubs=0

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
        AllNoAbortedOld=AllNoAborted
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
        AllNoAborted=0
        AllNoAddedInitiators=0
        AllNoInitDets=0
        AllNoNonInitDets=0
        AllNoInitWalk=0
        AllNoNonInitWalk=0
        AllNoDoubSpawns=0
        AllNoExtraInitDoubs=0



        RETURN
    END SUBROUTINE CalcNewShift
    SUBROUTINE WriteHistogram()
        CALL Stop_All("WriteHistogram","WriteHistogram not currently coded for serial.")
    END SUBROUTINE WriteHistogram
#endif
! AJWT
! Bringing you a better FciMCPar.  A vision for the future...
!
!  This section contains parts of FciMCPar which are not dependent on MPI commands.
!  It's not yet complete, but at least compiles and runs

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

!This routine looks at the change in residual particle number over a number of cycles, and adjusts the 
!value of the diagonal shift in the hamiltonian in order to compensate for this
    SUBROUTINE UpdateDiagSftPar()
        USE CalcData , only : tGlobalSftCng
        INTEGER :: j,k,GrowthSteps,MaxCulls,error
        LOGICAL :: Changed

        Changed=.false.
        IF(tGlobalSftCng) THEN
!!            CALL MPI_AllReduce(NoCulls,MaxCulls,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,error)
            CALL MPIIReduce(NoCulls,1,MPI_MAX,MaxCulls)
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

            ELSEIF(tTruncInitiator.or.tDelayTruncInit) THEN
                WRITE(6,"(A2,A10,A16,A10,A16,A12,3A11,3A17,2A10,A13,A12,A13)") "#","Step","Shift","WalkerCng","GrowRate","TotWalkers","Annihil","NoDied","NoBorn","Proj.E","Av.Shift","Proj.E.ThisCyc",&
&               "NoatHF","NoatDoubs","AccRat","UniqueDets","IterTime"

                WRITE(15,"(A2,A10,A16,A10,A16,A12,3A13,3A17,2A10,A13,A12,A13,A17,A13,A13)") "#","Step","Shift","WalkerCng","GrowRate","TotWalkers","Annihil","NoDied","NoBorn","Proj.E","Av.Shift",&
&               "Proj.E.ThisCyc","NoatHF","NoatDoubs","AccRat","UniqueDets","IterTime","FracSpawnFrmSing","WalkDiffProc","TotImagTime"

                WRITE(16,"(A2,A10,2A15,2A16,2A20,2A18)") "# ","Step","No Aborted","NoAddedtoInit","FracDetsInit","FracWalksInit","NoDoubSpawns","NoExtraDoubs","InstAbortShift","AvAbortShift"

            ELSE
                WRITE(6,"(A)") "       Step     Shift      WalkerCng    GrowRate       TotWalkers    Annihil    NoDied    NoBorn    Proj.E          Av.Shift     Proj.E.ThisCyc   NoatHF NoatDoubs      AccRat     UniqueDets     IterTime"
                WRITE(15,"(A)") "#       Step     Shift      WalkerCng    GrowRate       TotWalkers    Annihil    NoDied    NoBorn    Proj.E          Av.Shift"&
&              // "Proj.E.ThisCyc   NoatHF NoatDoubs       AccRat     UniqueDets     IterTime    FracSpawnFromSing    WalkersDiffProc"
            
            ENDIF
        ENDIF

    END SUBROUTINE WriteFciMCStatsHeader

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

            ELSEIF(tTruncInitiator.or.tDelayTruncInit) THEN
                WRITE(15,"(I12,G16.7,I10,G16.7,I12,3I13,3G17.9,2I10,G13.5,I12,G13.5,G17.5,I13,G13.5)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,   &
 &                  INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AvDiagSft,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),IterTime,&
 &                  REAL(AllSpawnFromSing)/REAL(AllNoBorn),WalkersDiffProc,TotImagTime

                WRITE(16,"(I12,2I15,2G16.7,2I20,2G18.7)") Iter+PreviousCycles,AllNoAborted,AllNoAddedInitiators,(REAL(AllNoInitDets)/REAL(AllNoNonInitDets)),(REAL(AllNoInitWalk)/REAL(AllNoNonInitWalk)),&
 &                  AllNoDoubSpawns,AllNoExtraInitDoubs,DiagSftAbort,AvDiagSftAbort
                WRITE(6,"(I12,G16.7,I10,G16.7,I12,3I11,3G17.9,2I10,G13.5,I12,G13.5)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,    &
 &                  INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AvDiagSft,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),IterTime

            ELSE
!                WRITE(15,"(I12,G16.7,I9,G16.7,I12,3I11,3G17.9,2I10,2G13.5,2I6)") Iter+PreviousCycles,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,   &
! &                  AllTotWalkers,AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,ProjEIter,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
!                WRITE(6,"(I12,G16.7,I9,G16.7,I12,3I11,3G17.9,2I10,2G13.5,2I6)") Iter+PreviousCycles,DiagSft,AllTotWalkers-AllTotWalkersOld,AllGrowRate,    &
! &                  AllTotWalkers,AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,ProjEIter,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,AllMeanExcitLevel,AllMinExcitLevel,AllMaxExcitLevel
                WRITE(15,"(I12,G16.7,I10,G16.7,I12,3I13,3G17.9,2I10,G13.5,I12,G13.5,G13.5,I10,G13.5)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,   &
 &                  INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AvDiagSft,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),IterTime,REAL(AllSpawnFromSing)/REAL(AllNoBorn),WalkersDiffProc,TotImagTime
                WRITE(6,"(I12,G16.7,I10,G16.7,I12,3I11,3G17.9,2I10,G13.5,I12,G13.5)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,    &
 &                  INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AvDiagSft,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),IterTime
            ENDIF
            
            CALL FLUSH(6)
            CALL FLUSH(15)
            
        ENDIF

        RETURN

    END SUBROUTINE WriteFCIMCStats


    SUBROUTINE SetupParameters()
        use SystemData, only : tUseBrillouin,iRanLuxLev,tSpn,tHPHFInts,tRotateOrbs,tNoBrillouin,tROHF,tFindCINatOrbs,nOccBeta,nOccAlpha,tUHF
        use SystemData, only : tFixLz,LzTot,BasisFN
        USE mt95 , only : genrand_init
        use CalcData, only : EXCITFUNCS,tFCIMC
        use Calc, only : VirtCASorbs,OccCASorbs,FixShift,G_VMC_Seed
        use Determinants , only : GetH0Element3
        use SymData , only : nSymLabels,SymLabelList,SymLabelCounts
        use Logging , only : tTruncRODump
        use FciMCLoggingMOD , only : InitTriHElStats,InitSpinCoupHel
        use DetCalc, only : NMRKS,tagNMRKS,FCIDets
        use SymExcit3, only : CountExcitations3 
        use DetBitOps, only: CountBits
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet,HFDetTest(NEl),Seed,alpha,beta,symalpha,symbeta,endsymstate
        INTEGER :: DetLT,VecSlot,error,HFConn,MemoryAlloc,iMaxExcit,nStore(6),nJ(Nel),BRR2(nBasis),LargestOrb,nBits,HighEDet(NEl),iLutTemp(0:NIfTot)
        TYPE(HElement) :: rh,TempHii
        TYPE(BasisFn) HFSym
        REAL*8 :: TotDets,SymFactor,Choose
        CHARACTER(len=*), PARAMETER :: this_routine='SetupParameters'
        CHARACTER(len=12) :: abstr
        LOGICAL :: tSuccess,tFoundOrbs(nBasis),tTurnBackBrillouin
        REAL :: Gap
        INTEGER :: nSingles,nDoubles,HFLz

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
            IF(tTruncInitiator.or.tDelayTruncInit) OPEN(16,file='INITIATORStats',status='unknown')
        ENDIF

!Store information specifically for the HF determinant
        ALLOCATE(HFDet(NEl),stat=ierr)
        CALL LogMemAlloc('HFDet',NEl,4,this_routine,HFDetTag)
        do i=1,NEl
            HFDet(i)=FDet(i)
        enddo
        HFHash=CreateHash(HFDet)
        CALL GetSym(HFDet,NEl,G1,NBasisMax,HFSym)
        WRITE(6,"(A,I5)") "Symmetry of reference determinant is: ",INT(HFSym%Sym%S,4)
        IF(tFixLz) THEN
            CALL GetLz(HFDet,NEl,HFLz)
            WRITE(6,"(A,I5)") "Ml value of reference determinant is: ",HFLz
            IF(HFLz.ne.LzTot) THEN
                CALL Stop_All("SetupParameters","Chosen reference determinant does not have the same Lz value as indicated in the input.")
            ENDIF
        ENDIF
        
        IF(tFixLz.and.(.not.tNoBrillouin)) THEN
            WRITE(6,*) "Turning Brillouins theorem off since we are using non-canonical complex orbitals"
            tNoBrillouin=.true.
        ENDIF
        
!test the encoding of the HFdet to bit representation.
        ALLOCATE(iLutHF(0:NIfTot),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Cannot allocate memory for iLutHF")
        CALL EncodeBitDet(HFDet,iLutHF)
!Test that the bit operations are working correctly...
        CALL DecodeBitDet(HFDetTest,iLutHF)
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
        nBits = CountBits(iLutHF, NIfD, NEl)
        IF(nBits.ne.NEl) THEN
            CALL Stop_All(this_routine,"CountBits FAIL")
        ENDIF

!Check that the symmetry routines have set the symmetry up correctly...
        tSuccess=.true.
        tFoundOrbs(:)=.false.

        IF((.not.tHub).and.(.not.tUEG)) THEN
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
            WRITE(6,*) "Symmetry information of orbitals not the same in alpha and beta pairs."
            WRITE(6,*) "Symmetry now set up in terms of spin orbitals"
            WRITE(6,*) "I strongly suggest you check that the reference energy is correct."
            IF(.not.tNonUniRandExcits) CALL Stop_All(this_routine,'ERROR. Need to use the non-uniform random excitation generators when &
            & we have odd symmetries and the symmetry label lists are stored in spin orbitals.')
        ELSE
            WRITE(6,*) "Simply transferring this into a spin orbital representation."
        ENDIF
! From now on, the orbitals are contained in symlabellist2 and symlabelcounts2 rather than the original arrays.
! These are stored using spin orbitals.

        IF(tNonUniRandExcits) THEN
!Assume that if we want to use the non-uniform random excitation generator, we also want to use the NoSpinSym full excitation generators if they are needed. 
            tNoSpinSymExcitgens=.true.   
        ENDIF

!If using a CAS space truncation, write out this CAS space
        IF(tTruncCAS) THEN
            WRITE(6,*) "Truncated CAS space detected. Writing out CAS space..."
            DO I=NEl-OccCASorbs+1,NEl
                WRITE(6,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(6,G1(BRR(I))%SYM,.FALSE.)
                WRITE(6,'(I4)',advance='no') G1(BRR(I))%Ml
                WRITE(6,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
            ENDDO
            WRITE(6,*) "-----------------------------------------------------"
            DO I=NEl+1,NEl+VirtCASOrbs
                WRITE(6,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(6,G1(BRR(I))%SYM,.FALSE.)
                WRITE(6,'(I4)',advance='no') G1(BRR(I))%Ml
                WRITE(6,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
            ENDDO
        ENDIF

                                        
!Setup excitation generator for the HF determinant. If we are using assumed sized excitgens, this will also be assumed size.
        IF(.not.tNoSpinSymExcitgens) THEN
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

        ELSE
            IF(tUEG.or.tHub) THEN
                exflag=2
            ELSE
                exflag=3
            ENDIF
!Count all possible excitations - put into HFConn
            CALL CountExcitations3(HFDet,exflag,nSingles,nDoubles)
            HFConn=nSingles+nDoubles
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
            !IF(tLatticeGens) CALL Stop_All("SetupParameters","Cannot use HPHF with model systems currently.")
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
                CALL EncodeBitDet(HighEDet,iLutTemp)
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
            ELSEIF(tUHF) THEN
                WRITE(6,*) "High spin calculation - but single excitations will *NOT* be used to calculate energy as this is an unrestricted calculation."
            ELSE
                CALL Stop_All("SetupParameters","High-spin, restricted calculation detected, but single excitations are not being used to calculate the energy.  &
                              & Either use the UHF keyword, or turn off brillouins theorem using NOBRILLOUINS, ROHF or ROTATEDORBS.")
            ENDIF
!            tRotatedOrbs=.true.
!        ELSEIF(LMS.ne.0) THEN
!            CALL Stop_All(this_routine,"Ms not equal to zero, but tSpn is false. Error here")
        ENDIF

        IF((tHub.and.tReal).or.(tRotatedOrbs)) THEN
            tNoBrillouin=.true.
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
        SumDiagSftAbort=0.D0
        AvDiagSftAbort=0.D0
        NoAborted=0
        NoAddedInitiators=0
        NoInitDets=0
        NoNonInitDets=0
        NoInitWalk=0
        NoNonInitWalk=0
        NoDoubSpawns=0
        NoExtraInitDoubs=0
        TotImagTime=0.D0

!Also reinitialise the global variables - should not necessarily need to do this...
        AllSumENum=0.D0
        AllNoatHF=0
        AllNoatDoubs=0
        AllSumNoatHF=0.D0
        AllGrowRate=0.D0
        AllGrowRateAbort=0.D0
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
        AllNoAborted=0
        AllNoAddedInitiators=0
        AllNoInitDets=0
        AllNoNonInitDets=0
        AllNoInitWalk=0
        AllNoNonInitWalk=0
        AllNoDoubSpawns=0
        AllNoExtraInitDoubs=0
 

        IF(tHistSpawn.or.(tCalcFCIMCPsi.and.tFCIMC).or.tHistHamil) THEN
            ALLOCATE(HistMinInd(NEl))
            ALLOCATE(HistMinInd2(NEl))
            maxdet=0
            do i=1,nel
                maxdet=maxdet+2**(nbasis-i)
            enddo

            IF(.not.allocated(FCIDets)) THEN
                CALL Stop_All(this_routine,"A Full Diagonalization is required in the same calculation before histogramming can occur.")
            ENDIF

            IF(tHistHamil) THEN
                WRITE(6,*) "Histogramming total Hamiltonian, with Dets=", Det
                ALLOCATE(HistHamil(1:det,1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                HistHamil(:,:)=0.D0
                ALLOCATE(AvHistHamil(1:det,1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                AvHistHamil(:,:)=0.D0
                IF(iProcIndex.eq.0) THEN
                    ALLOCATE(AllHistHamil(1:det,1:det),stat=ierr)
                    IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                    AllHistHamil(:,:)=0.D0
                    ALLOCATE(AllAvHistHamil(1:det,1:det),stat=ierr)
                    IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                    AllAvHistHamil(:,:)=0.D0
                ENDIF
            ELSE
                WRITE(6,*) "Histogramming spawning wavevector, with Dets=", Det
                ALLOCATE(Histogram(1:det),stat=ierr)
                IF(ierr.ne.0) THEN
                    CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
                ENDIF
                Histogram(:)=0.D0
                ALLOCATE(AllHistogram(1:det),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All("SetupParameters","Error assigning memory for histogramming arrays (could deallocate NMRKS to save memory?)")
            ENDIF
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
            IF(iProcIndex.eq.Root) THEN
                ALLOCATE(AllHistogram(1:iNoBins))
                ALLOCATE(AllAttemptHist(1:iNoBins))
                ALLOCATE(AllSpawnHist(1:iNoBins))
                ALLOCATE(AllSinglesHist(1:iOffDiagNoBins))
                ALLOCATE(AllDoublesHist(1:iOffDiagNoBins))
                ALLOCATE(AllSinglesAttemptHist(1:iOffDiagNoBins))
                ALLOCATE(AllDoublesAttemptHist(1:iOffDiagNoBins))
                ALLOCATE(AllSinglesHistOccOcc(1:iOffDiagNoBins))
                ALLOCATE(AllSinglesHistOccVirt(1:iOffDiagNoBins))
                ALLOCATE(AllSinglesHistVirtOcc(1:iOffDiagNoBins))
                ALLOCATE(AllSinglesHistVirtVirt(1:iOffDiagNoBins))
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
     &      tPrintTriConnections.or.tHistTriConHEls.or.tCalcFCIMCPsi.or.tTruncCAS.or.tListDets.or.tPartFreezeCore.or.tPartFreezeVirt.or.tUEG.or.tHistHamil.or.TReadPops) THEN
            WRITE(6,*) "It is not possible to use to clean spawning routine..."
        ELSE
            WRITE(6,*) "Clean spawning routine in use..."
            tCleanRun=.true.
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

        IF(TPopsFile) THEN
            IF(mod(iWritePopsEvery,StepsSft).ne.0) CALL Warning(this_routine,"POPSFILE writeout should be a multiple of the update cycle length.")
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
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
            ValidSpawnedList(:)=InitialSpawnedSlots(:)
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
        IF(tPartFreezeVirt) THEN
            WRITE(6,'(A,I4,A,I5)') 'Partially freezing the highest ',NVirtPartFrozen,' virtual spin orbitals so that no more than ',NElVirtFrozen,' electrons occupy these orbitals.'
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
            IF((.not.tNonUniRandExcits).or.tStarOrbs.or.tTruncSpace.or.tTruncCAS.or.tListDets.or.tPartFreezeCore.or.tPartFreezeVirt) THEN
                CALL Stop_All("SetupParameters","Cannot use HighExcitsSing without Nonuniformrandexcits, or with starorbs or truncated spaces...")
            ENDIF
        ENDIF

        IF(tMultipleDetsSpawn) THEN
!We need to store a list of all double excitations of HF.
            CALL StoreDoubs()
        ENDIF

    END SUBROUTINE SetupParameters
    LOGICAL FUNCTION TestifDETinCAS(CASDet)
        INTEGER :: k,z,CASDet(NEl)
        LOGICAL :: tElecInVirt

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
        tElecInVirt=.false.
        do k=1,NEl      ! running over all electrons
            if (SpinInvBRR(CASDet(k)).gt.CASmax) THEN
                tElecInVirt=.true.
                EXIT            
! if at any stage an electron has an energy greater than the CASmax value, the determinant can be ruled out
! of the active space.  Upon identifying this, it is not necessary to check the remaining electrons.
            else
                if (SpinInvBRR(CASDet(k)).le.CASmin) THEN
                    z=z+1
                endif
! while running over all electrons, the number that occupy orbitals equal to or below the CASmin cutoff are
! counted.
            endif
        enddo

        if(tElecInVirt.or.(z.ne.CASmin)) THEN
! if an electron is in an orbital above the active space, or the inactive orbitals are not full, the determinant is automatically ruled out.            
            TestifDETinCAS=.false.
        else
! no orbital in virtual and the inactive orbitals are completely full - det in active space.        
            TestifDETinCAS=.true.
        endif
        
        RETURN

    END FUNCTION TestifDETinCAS

!This function will tell us whether we should allow attempted spawning at an excitation when we are truncating the space.
!We pass in the excitation level of the original particle, the two representations of the excitation (we only need the bit-representation of the excitation
!for HPHF) and the magnitude of the excitation (for determinant representation).
    LOGICAL FUNCTION CheckAllowedTruncSpawn(WalkExcitLevel,nJ,iLutnJ,IC)
        INTEGER :: nJ(NEl),WalkExcitLevel,iLutnJ(0:NIfTot),ExcitLevel,IC,iGetExcitLevel_2,i,NoInFrozenCore,TotalLz,MinVirt
        INTEGER :: kx,ky,kz ! For UEG

        CheckAllowedTruncSpawn=.true.

        IF(tTruncSpace.and.(.not.tTruncInitiator)) THEN
!We are truncating the space by excitation level
            IF(tHPHF) THEN
!With HPHF, we can't rely on this, since one excitation could be a single, and one a double. Also, IC is not returned.
                ExcitLevel = FindBitExcitLevel(iLutHF, iLutnJ, ICILevel)
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

        IF((tTruncCAS.and.(.not.tTruncInitiator)).and.CheckAllowedTruncSpawn) THEN
!This flag determines if the FCI space is restricted by whether the determinants are in the predescribed CAS.
            IF(.not.TestIfDetinCAS(nJ)) THEN
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
                CALL EncodeBitDet(nJ,iLutnJ)
            ENDIF
            do i=1,NAllowedDetList
!                WRITE(6,*) ILutnJ,AllowedDetList(0:NIfTot,i)
                IF(DetBitEQ(iLutnJ,AllowedDetList(0:NIfTot,i),NIfDBO)) THEN
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
!            ELSE
!Either the 'partially frozen core' is completely full, or it has the allowed number of holes or less.                
!Allowed to spawn, CheckAllowedTruncSpawn=.true.
!                CheckAllowedTruncSpawn=.true.
            ENDIF

        ENDIF

        IF(tPartFreezeVirt) THEN
!Want to check if the determinant we're about to spawn on has more than the restricted number of electrons in the partially frozen virtual orbitals.

!Run through the electrons in nJ, count the number in the partially frozen virtual orbitals - ie those occupying orbitals with energy (from BRR) 
!greater than that of the minimum unfrozen virtual.
!If this is greater than NElVirtFrozen then spawning is forbidden.
            NoInFrozenCore=0
            MinVirt=nBasis-NVirtPartFrozen
!BRR(i)=j: orbital i is the j-th lowest in energy  
            do i=1,NEl
                IF(SpinInvBRR(nJ(i)).gt.MinVirt) NoInFrozenCore=NoInFrozenCore+1
                IF(NoInFrozenCore.gt.NElVirtFrozen) THEN
!There are more electrons in the partially frozen virtual orbitals than has been specified as allowed.
                    CheckAllowedTruncSpawn=.false.
                    EXIT   ! Can exit out of the loop if this is satisfied, since excitation will definitely be accepted.
                ENDIF
            enddo
!Either the 'partially frozen virtual orbitals' are completely empty, or have the allowed number of electrons or less.                
!Allowed to spawn, CheckAllowedTruncSpawn=.true.
!                CheckAllowedTruncSpawn=.true.
!Want to be able to make it unable to spawn, but not able to spawn again.
        ENDIF

!The excitation generators should now only generate Lz allowed excitations if tFixLz is true.
!        IF(tFixLz) THEN
!            CALL GetLz(nJ,NEl,TotalLz)      !This could be improved by just checking that the change in momentum from the excitation was zero.
!            IF(TotalLz.ne.LzTot) THEN
!                CheckAllowedTruncSpawn=.false.
!                WRITE(6,*) "FALSE ",TotalLz
!                WRITE(6,*) nJ(:)
!                CALL Stop_All("CheckAllowedTruncSpawn","Should not get here with new excitation generators.")
!            ELSE
!                CheckAllowedTruncSpawn=.true.
!!                WRITE(6,*) "TRUE ",TotalLz
!            ENDIF
!        ENDIF

        IF(tUEG) THEN
!Check to see if this is an allowed excitation
!by summing kx, ky and kz to zero over all the electrons.
            kx=0
            ky=0
            kz=0
            do i=1,NEl
                kx=kx+G1(nJ(i))%k(1)
                ky=ky+G1(nJ(i))%k(2)
                kz=kz+G1(nJ(i))%k(3)
            enddo
            if( .not.((kx.eq.0) .and. (ky.eq.0) .and. (kz.eq.0)) ) then
                CheckAllowedTruncSpawn=.false.
            endif
        ENDIF


    END FUNCTION CheckAllowedTruncSpawn



!This is the same as BinSearchParts1, but this time, the list to search is passed in as an argument. The list goes from 1 to Length, but only between MinInd and MaxInd is actually searched.
    SUBROUTINE BinSearchParts3(iLut,List,Length,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: iLut(0:NIfTot),MinInd,MaxInd,PartInd
        INTEGER :: List(0:NIfTot,Length),Length
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(List(:,MaxInd),iLut(:),NIfDBO)
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
            Comp=DetBitLT(List(:,N),iLut(:),NIfDBO)

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
                    Comp=DetBitLT(List(:,i+1),iLut(:),NIfDBO)
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
    
    SUBROUTINE CalcApproxpDoubles(HFConn)
        use SystemData , only : tAssumeSizeExcitgen
        use SymData , only : SymClassSize
        use SymExcit3 , only : CountExcitations3
        INTEGER :: HFConn,PosExcittypes,iTotal,i
        integer :: nSing, nDoub, ncsf, ExcitInd

        ! TODO: A better approximation for ncsf.
        if (tCSF) then
            ncsf = 10
        else
            ncsf = 0
        endif

        IF(tHub.or.tUEG) THEN
            IF(tReal) THEN
                WRITE(6,*) "Since we are using a real-space hubbard model, only single excitations are connected."
                WRITE(6,*) "Setting pDoub to 0.D0"
                pDoubles=0.D0
                RETURN
            ELSE
                WRITE(6,*) "Since we are using a momentum-space hubbard model/UEG, only double excitaitons are connected."
                WRITE(6,*) "Setting pDoub to 1.D0"
                pDoubles=1.D0
                RETURN
            ENDIF
        ENDIF

        IF(tNoSpinSymExcitgens) THEN
!NSing=Number singles from HF, nDoub=No Doubles from HF

            WRITE(6,"(A)") "Calculating approximate pDoubles for use with excitation generator by looking a excitations from HF."
            exflag=3
            CALL CountExcitations3(HFDet,exflag,nSing,nDoub)
            iTotal=nSing + nDoub + ncsf

        ELSE

            WRITE(6,"(A)") "Calculating approximate pDoubles for use with excitation generator by looking a excitations from HF."
            IF(tAssumeSizeExcitgen) THEN
                PosExcittypes=SymClassSize*NEL+NIfD+4
                iTotal=HFExcit%ExcitData(1)
            ELSE
                PosExcittypes=HFExcit%ExcitData(2)
                iTotal=HFExcit%ExcitData(23)
            ENDIF
            IF(iTotal.ne.HFConn) THEN
                CALL Stop_All("CalcApproxpDoubles","Number of excitations from HF determinant has been confused somewhere...")
            ENDIF
!            WRITE(6,*) "**********"
!            do i=0,100
!                WRITE(6,*) i,HFExcit%ExcitData(PosExcittypes+i)
!            enddo
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

        ENDIF

        WRITE(6,"(I7,A,I7,A)") NDoub, " double excitations, and ",NSing," single excitations found from HF. This will be used to calculate pDoubles."

        IF(SinglesBias.ne.1.D0) THEN
            WRITE(6,*) "Singles Bias detected. Multiplying single excitation connectivity of HF determinant by ",SinglesBias," to determine pDoubles."
        ENDIF

        IF((NSing+nDoub+ncsf).ne.iTotal) THEN
            CALL Stop_All("CalcApproxpDoubles","Sum of number of singles and number of doubles does not equal total number of excitations")
        ENDIF
        IF((NSing.eq.0).or.(NDoub.eq.0)) THEN
            WRITE(6,*) "Number of singles or doubles found equals zero. pDoubles will be set to 0.95. Is this correct?"
            pDoubles = 0.95
            pSingles = 0.05
            RETURN
        elseif ((NSing < 0) .or. (NDoub < 0) .or. (ncsf < 0)) then
            call stop_all("CalcApproxpDoubles", &
                          "Number of singles, doubles or Yamanouchi symbols &
                          &found to be a negative number. Error here.")
        endif

        ! Set pDoubles to be the fraction of double excitations.
        ! If using CSFs, also consider only changing Yamanouchi Symbol
        if (tCSF) then
            pDoubles = real(nDoub,r2) / &
                   ((real(nSing,r2)*SinglesBias)+real(nDoub,r2)+real(ncsf,r2))
            pSingles = real(nSing,r2) / &
                   ((real(nSing,r2)*SinglesBias)+real(nDoub,r2)+real(ncsf,r2))

        else
            pDoubles = real(nDoub,r2) * SinglesBias / &
                   ((real(NSing,r2)*SinglesBias) + real(NDoub,r2))
            pSingles = real(nSing,r2) / &
                   ((real(nSing,r2)*SinglesBias) + real(nDoub,r2))
        endif

        IF(SinglesBias.ne.1.D0) THEN
            write (6, '("pDoubles set to ", f14.6, &
                       &" rather than (without bias): ", f14.6)') &
                       pDoubles, real(nDoub,r2) / real(iTotal,r2)
            write (6, '("pSingles set to ", f14.6, &
                       &" rather than (without bias): ", f14.6)') &
                       pSingles, real(nSing,r2) / real(iTotal,r2)

            WRITE(6,"(A,F14.6,A,F14.6)") "pDoubles set to: ",pDoubles, " rather than (without bias): ",real(nDoub,r2)/real(iTotal,r2)
        ELSE
            write (6, '("pDoubles set to: ", f14.6)') pDoubles
            write (6, '("pSingles set to: ", f14.6)') pSingles
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

        IF(ALLOCATED(SpinInvBRR)) RETURN
            
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
        INTEGER :: DetCurr(0:NIfTot),i,Elecs,j
        INTEGER(KIND=i2) :: CreateHashBit

        CreateHashBit=0
        Elecs=0
        do i=0,NIfTot
            do j=0,31
                IF(BTEST(DetCurr(i),j)) THEN
                    CreateHashBit=(1099511628211_8*CreateHashBit)+((i*32)+(j+1))*j
                    Elecs=Elecs+1
                    IF(Elecs.eq.NEl) RETURN
                ENDIF
            enddo
        enddo

    END FUNCTION CreateHashBit

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
            IF(IsNullDet(nJ)) EXIT
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
        use SymExcit3 , only : CountExcitations3,GenExcitations3
        INTEGER :: iMaxExcit,nStore(6),ExcitLength,nJ(NEl),ierr,iExcit,VecSlot,nSingles,ExcitMat3(2,2)
        INTEGER , ALLOCATABLE :: ExcitGenTemp(:)
        LOGICAL :: tAllExcitFound,tParity


        IF(tAssumeSizeExcitgen) THEN
            CALL Stop_All("StoreDoubs","Cannot have assumed sized excitation generators for full enumeration of determinants")
        ENDIF
        IF(tUseBrillouin) THEN
            CALL Stop_All("StoreDoubs","Cannot have Brillouin theorem as now storing singles too...")
        ENDIF

        
!NoDoubs here is actually the singles + doubles of HF
        IF(tNoSpinSymExcitgens) THEN
            exflag=3
            CALL CountExcitations3(HFDet,exflag,nSingles,NoDoubs)
            NoDoubs=nSingles+NoDoubs
        ELSE
            CALL GetSymExcitCount(HFExcit%ExcitData,NoDoubs)
        ENDIF

        ALLOCATE(DoublesDets(NEl,NoDoubs),stat=ierr)
        CALL LogMemAlloc('DoublesDets',NoDoubs*NEl,4,"StoreDoubs",DoublesDetsTag,ierr)
        DoublesDets(1:NEl,1:NoDoubs)=0
        
        VecSlot=1           !This is the next free slot in the DoublesDets array

        IF(tNoSpinSymExcitgens) THEN

            tAllExcitFound=.false.
            ExcitMat3(:,:)=0
!An exflag of anything but 1 or 2 indicates both the single and double excitations should be found.            
            exflag=3

            do while (.not.tAllExcitFound)
                CALL GenExcitations3(HFDet,iLutHF,nJ,exflag,ExcitMat3,tParity,tAllExcitFound)
                IF(tAllExcitFound) EXIT
                DoublesDets(1:NEl,VecSlot)=nJ(:)
                VecSlot=VecSlot+1
            enddo

        ELSE
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
                IF(IsNullDet(nJ)) EXIT
!                IF(iExcit.ne.2) THEN
!                    CALL Stop_All("StoreDoubles","Error - excitations other than doubles being generated in DoublesDets wavevector code")
!                ENDIF

                DoublesDets(1:NEl,VecSlot)=nJ(:)
                VecSlot=VecSlot+1

            enddo
        ENDIF

!This means that now NoDoubs is double excitations AND singles
!        NoDoubs=VecSlot-1

        IF(VecSlot.ne.(NoDoubs+1)) THEN
            WRITE(6,*) VecSlot,NoDoubs
            CALL Stop_All("StoreDoubs","Problem enumerating all double excitations")
        ENDIF

        DEALLOCATE(ExcitGenTemp)

    END SUBROUTINE StoreDoubs

!Initialize the Histogramming searching arrays if necessary
    SUBROUTINE InitHistMin()
        IF(tHistSpawn.or.tCalcFCIMCPsi.and.(Iter.ge.NHistEquilSteps)) THEN
            IF(Iter.eq.NHistEquilSteps) THEN
                IF(iProcIndex.eq.Root) WRITE(6,*) 'The iteration is equal to HISTEQUILSTEPS.  Beginning to histogram.'
            ENDIF
            HistMinInd(1:NEl)=FCIDetIndex(1:NEl)    !This is for the binary search when histogramming
        ENDIF
    END SUBROUTINE InitHistMin

!This routine sums in the energy contribution from a given walker and updates stats such as mean excit level
!AJWT added optional argument dProb which is a probability that whatever gave this contribution as generated.
!  It defaults to 1, and weights the contribution of this det. (Only in the projected energy)
    SUBROUTINE SumEContrib(DetCurr,ExcitLevel,WSign,iLutCurr,HDiagCurr,dProbFin)
        use SystemData, only : tNoBrillouin
        use CalcData, only: tFCIMC
        INTEGER :: DetCurr(NEl),ExcitLevel,i,HighIndex,LowIndex,iLutCurr(0:NIfTot),WSign,Bin
        INTEGER :: PartInd,iLutSym(0:NIfTot),OpenOrbs
        LOGICAL :: CompiPath,tSuccess,iLut2(0:NIfTot)
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
!            WRITE(6,*) 2,SumENum,(REAL(HOffDiag%v,r2)*WSign/dProbFin)     !This is simply the Hij*sign summed over the course of the update cycle

            
            
!        ELSE
!            AvSign=AvSign+REAL(WSign,r2)

        ELSEIF(ExcitLevel.eq.1) THEN
          if(tNoBrillouin) then
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
!            WRITE(6,*) 1,SumENum,(REAL(HOffDiag%v,r2)*WSign/dProbFin)     !This is simply the Hij*sign summed over the course of the update cycle
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
        IF((tHistSpawn.or.(tCalcFCIMCPsi.and.tFCIMC)).and.(Iter.ge.NHistEquilSteps)) THEN
            IF(ExcitLevel.eq.NEl) THEN
                CALL BinSearchParts2(iLutCurr,HistMinInd(ExcitLevel),Det,PartInd,tSuccess)
                if(tFCIMC) HistMinInd(ExcitLevel)=PartInd  !CCMC doesn't sum particle contributions in order, so we must search the whole space again
            ELSEIF(ExcitLevel.eq.0) THEN
                PartInd=1
                tSuccess=.true.
            ELSE
                CALL BinSearchParts2(iLutCurr,HistMinInd(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartInd,tSuccess)
                if(tFCIMC) HistMinInd(ExcitLevel)=PartInd  !CCMC doesn't sum particle contributions in order, so we must search the whole space again
            ENDIF
            IF(tSuccess) THEN
                IF(tHPHF) THEN
                    CALL FindExcitBitDetSym(iLutCurr,iLutSym)
                    IF(.not.DetBitEQ(iLutCurr,iLutSym)) THEN
                        IF(tFlippedSign) THEN
                            Histogram(PartInd)=Histogram(PartInd)-(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                        ELSE
                            Histogram(PartInd)=Histogram(PartInd)+(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                        ENDIF
                    ELSE
                        IF(tFlippedSign) THEN
                            Histogram(PartInd)=Histogram(PartInd)-REAL(WSign,r2)/dProbFin
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-REAL(WSign,r2)/dProbFin
                        ELSE
                            Histogram(PartInd)=Histogram(PartInd)+REAL(WSign,r2)/dProbFin
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+REAL(WSign,r2)/dProbFin
                        ENDIF
                    ENDIF
                ELSE
                    IF(tFlippedSign) THEN
                        Histogram(PartInd)=Histogram(PartInd)-REAL(WSign,r2)/dProbFin
                        IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-REAL(WSign,r2)/dProbFin
                    ELSE
                        Histogram(PartInd)=Histogram(PartInd)+REAL(WSign,r2)/dProbFin
                        IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+REAL(WSign,r2)/dProbFin
                    ENDIF
                ENDIF
                IF(tHPHF) THEN
!With HPHF space, we need to also include the spin-coupled determinant, which will have the same amplitude as the original determinant, unless it is antisymmetric.
                    IF(.not.DetBitEQ(iLutCurr,iLutSym)) THEN
                        IF(ExcitLevel.eq.NEl) THEN
                            CALL BinSearchParts2(iLutSym,FCIDetIndex(ExcitLevel),Det,PartInd,tSuccess)
                        ELSEIF(ExcitLevel.eq.0) THEN
                            PartInd=1
                            tSuccess=.true.
                        ELSE
                            CALL BinSearchParts2(iLutSym,FCIDetIndex(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartInd,tSuccess)
                        ENDIF
                        IF(tSuccess) THEN
                            CALL CalcOpenOrbs(iLutSym,OpenOrbs)
                            IF(tFlippedSign) THEN
                                IF(mod(OpenOrbs,2).eq.1) THEN
                                    Histogram(PartInd)=Histogram(PartInd)+(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                                ELSE
                                    Histogram(PartInd)=Histogram(PartInd)-(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                                ENDIF
                            ELSE
                                IF(mod(OpenOrbs,2).eq.1) THEN
                                    Histogram(PartInd)=Histogram(PartInd)-(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                                ELSE
                                    Histogram(PartInd)=Histogram(PartInd)+(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+(REAL(WSign,r2)/SQRT(2.0))/dProbFin
                                ENDIF
                            ENDIF
                        ELSE
                            WRITE(6,*) DetCurr(:)
                            WRITE(6,*) "***",iLutSym(0:NIfTot)
                            WRITE(6,*) "***",ExcitLevel,Det
                            CALL Stop_All("SumEContrib","Cannot find corresponding spin-coupled FCI determinant when histogramming")
                        ENDIF
                    ENDIF
                ENDIF
            ELSE
                WRITE(6,*) DetCurr(:)
                WRITE(6,*) "***",iLutCurr(0:NIfTot)
                WRITE(6,*) "***",ExcitLevel,HistMinInd(ExcitLevel),Det
                Call WriteBitDet(6,iLutCurr(0:NIfTot),.true.)
                CALL Stop_All("SumEContrib","Cannot find corresponding FCI determinant when histogramming")
            ENDIF
        ELSEIF(tHistEnergies) THEN
!This wil histogramm the energies of the particles, rather than the determinants themselves.
            Bin=INT(HDiagCurr/BinRange)+1
            IF(Bin.gt.iNoBins) THEN
                CALL Stop_All("SumEContrib","Histogramming energies higher than the arrays can cope with. Increase iNoBins or BinRange")
            ENDIF
            Histogram(Bin)=Histogram(Bin)+real(abs(WSign),r2)
        ENDIF

        IF(tPrintOrbOcc.and.(Iter.ge.StartPrintOrbOcc)) THEN
            do i=1,NEl
                OrbOccs(DetCurr(i))=OrbOccs(DetCurr(i))+ABS(WSign)
            enddo
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

END MODULE FciMCParMod

!The Exitgen manipulation routines are outside the module to allow them to be used in AnnihilationMod
!This routine copies an excitation generator from origExcit to NewExit, if the original claims that it is for the correct determinant
SUBROUTINE CopyExitgenPar(OrigExit,NewExit,DelOldCopy)
    use FCIMCParMod
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
    use FCIMCParMod
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
    use FCIMCParMod
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


!This is the same as BinSearchParts1, but this time, it searches though the full list of determinants created by the full diagonalizer when the histogramming option is on.
!This is outside the module so it is accessible to AnnihilateMod
SUBROUTINE BinSearchParts2(iLut,MinInd,MaxInd,PartInd,tSuccess)
    use SystemData , only : NIfTot,NIfDBO
    use DetCalc , only : FCIDets
    use DetBitOps, only: DetBitLT
    INTEGER :: iLut(0:NIfTot),MinInd,MaxInd,PartInd
    INTEGER :: i,j,N,Comp
    LOGICAL :: tSuccess

!    WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!    CALL FLUSH(6)
    i=MinInd
    j=MaxInd
    IF(i-j.eq.0) THEN
        Comp=DetBitLT(FCIDets(:,MaxInd),iLut(:),NIfDBO)
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
!        WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
        Comp=DetBitLT(FCIDets(:,N),iLut(:),NIfDBO)

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
                Comp=DetBitLT(FCIDets(:,i+1),iLut(:),NIfDBO)
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

    
