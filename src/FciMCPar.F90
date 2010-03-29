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
    use SystemData , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,nMsh,Arr,LMS,tHPHF,NIfD,NIfTot,NIfDBO,NIfY
    use SystemData , only : tHub,tReal,tRotatedOrbs,tFindCINatOrbs,tFixLz,LzTot,tUEG, tLatticeGens,tCSF
    use CalcData , only : InitWalkers,NMCyc,DiagSft,Tau,SftDamp,StepsSft,OccCASorbs,VirtCASorbs,tFindGroundDet
    use CalcData , only : NEquilSteps,TReadPops,tRegenDiagHEls
    use CalcData , only : GrowMaxFactor,CullFactor,TStartSinglePart,ScaleWalkers,MaxNoatHF,HFPopThresh
    use CalcData , only : tCCMC,tTruncCAS,tTruncInitiator,tDelayTruncInit,IterTruncInit,NShiftEquilSteps,tWalkContGrow,tMCExcits,NoMCExcits
    use HPHFRandExcitMod , only : FindExcitBitDetSym,GenRandHPHFExcit,GenRandHPHFExcit2Scratch
    use Determinants, only: FDet, get_helement, write_det
    USE DetCalc , only : ICILevel,nDet,Det,FCIDetIndex
    use GenRandSymExcitNUMod , only : GenRandSymExcitScratchNU,GenRandSymExcitNU,ScratchSize
    use IntegralsData , only : fck,NMax,UMat,tPartFreezeCore,NPartFrozen,NHolesFrozen,tPartFreezeVirt,NVirtPartFrozen,NElVirtFrozen
    USE UMatCache , only : GTID
    USE Logging , only : iWritePopsEvery,TPopsFile,iPopsPartEvery,tBinPops,tHistSpawn,iWriteHistEvery,tHistEnergies,IterShiftBlock,AllHistInitPops
    USE Logging , only : BinRange,iNoBins,OffDiagBinRange,OffDiagMax,AllHistInitPopsTag
    USE Logging , only : tPrintFCIMCPsi,tCalcFCIMCPsi,NHistEquilSteps,tPrintOrbOcc,StartPrintOrbOcc
    USE Logging , only : tHFPopStartBlock,tIterStartBlock,IterStartBlocking,HFPopStartBlocking,tInitShiftBlocking,tHistHamil,iWriteHamilEvery,HistInitPopsTag
    USE Logging , only : OrbOccs,OrbOccsTag,tPrintPopsDefault,iWriteBlockingEvery,tBlockEveryIteration,tHistInitPops,HistInitPopsIter,HistInitPops
    USE SymData , only : nSymLabels
    USE dSFMT_interface , only : genrand_real2_dSFMT
    USE Parallel
    USE FciMCData
    USE AnnihilationMod
    use DetBitops, only: EncodeBitDet, DecodeBitDet, DetBitEQ, DetBitLT
    use DetBitOps, only: FindExcitBitDet, FindBitExcitLevel
    use csf, only: get_csf_bit_yama, iscsf, csf_orbital_mask
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use util_mod, only: choose
    use constants, only: dp
    IMPLICIT NONE
    SAVE

    contains

#ifdef PARALLEL

    SUBROUTINE FciMCPar(Weight,Energyxw)
        use soft_exit, only : ChangeVars 
        use CalcData, only : iFullSpaceIter
        use UMatCache, only : UMatInd
        use FciMCLoggingMOD , only : FinaliseBlocking,FinaliseShiftBlocking,PrintShiftBlocking,PrintBlocking
        USE FciMCLoggingMOD , only : SumInErrorContrib,WriteInitPops
        use RotateOrbsMod , only : RotateOrbs
        use NatOrbsMod , only : PrintOrbOccs
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: i,j,error,HFConn
        CHARACTER(len=*), PARAMETER :: this_routine='FciMCPar'
        TYPE(HElement) :: Hamii
        LOGICAL :: TIncrement,tWritePopsFound,tSoftExitFound,tSingBiasChange
        REAL(4) :: s,etime,tstart(2),tend(2)
        INTEGER :: MaxWalkers,MinWalkers
        real*8 :: AllTotWalkers,MeanWalkers,Inpair(2),Outpair(2)

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

            s=etime(tstart)
            IF(tCCMC) THEN
                CALL PerformCCMCCycPar()
            ELSE
                CALL PerformFCIMCycPar()
            ENDIF
            s=etime(tend)
            IterTime=IterTime+(tend(1)-tstart(1))

            IF(tBlockEveryIteration) THEN
                Inpair(1)=REAL(HFIter,dp)
                Inpair(2)=ENumIter
                CALL MPIDSumRootArr(Inpair,2,Outpair,Root)
                IterEnergy=Outpair(2)/Outpair(1)
                IF(tErrorBlocking.and.(iProcIndex.eq.Root)) CALL SumInErrorContrib(Iter,Outpair(2),Outpair(1))
                ENumIter=0.D0
                HFIter=0
            ENDIF

            IF(mod(Iter,StepsSft).eq.0) THEN
!This will communicate between all nodes, find the new shift (and other parameters) and broadcast them to the other nodes.
                CALL CalcNewShift()

                IF((tTruncCAS.or.tTruncSpace.or.tHFInitiator).and.(Iter.gt.iFullSpaceIter).and.(iFullSpaceIter.ne.0)) THEN
!Test if we want to expand to the full space if an EXPANDSPACE variable has been set
                    IF(tHistSpawn.or.tCalcFCIMCPsi) THEN
                        IF(iProcIndex.eq.0) WRITE(6,*) "Unable to expand space since histgramming the wavefunction..."
                    ELSE
                        ICILevel=0
                        tTruncSpace=.false.
                        tTruncCAS=.false.
                        tHFInitiator=.false.
                        IF(tTruncInitiator) tTruncInitiator=.false.
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
                    CALL WriteToPopsfileParOneArr()
                ENDIF
                IF(tSingBiasChange) THEN
                    CALL CalcApproxpDoubles(HFConn)
                ENDIF
            
            ENDIF

            IF(mod(Iter,iWriteBlockingEvery).eq.0) THEN
                !Every 100 update cycles, write out a new blocking file.
                IF(tErrorBlocking.and.(Iter.gt.IterStartBlocking)) CALL PrintBlocking(Iter) 
                IF(tShiftBlocking.and.(Iter.gt.(VaryShiftIter+IterShiftBlock))) CALL PrintShiftBlocking(Iter)
            ENDIF

            IF(TPopsFile.and.(.not.tPrintPopsDefault).and.(mod(Iter,iWritePopsEvery).eq.0)) THEN
!This will write out the POPSFILE if wanted
                CALL WriteToPopsfileParOneArr()
            ENDIF
!            IF(TAutoCorr) CALL WriteHistogrammedDets()

            IF(tHistSpawn.and.(mod(Iter,iWriteHistEvery).eq.0)) THEN
                CALL WriteHistogram()
            ENDIF
            IF(tHistHamil.and.(mod(Iter,iWriteHamilEvery).eq.0)) THEN
                CALL WriteHamilHistogram()
            ENDIF

            IF(tHistInitPops.and.(MOD(Iter,HistInitPopsIter).eq.0)) THEN
                IF(iProcIndex.eq.0) WRITE(6,*) 'Writing out the spread of the initiator determinant populations.'
                CALL WriteInitPops(Iter+PreviousCycles)
            ENDIF

            Iter=Iter+1
!End of MC cycle
        enddo

        IF(TIncrement) Iter=Iter-1     !Reduce the iteration count for the POPSFILE since it is incremented upon leaving the loop (if done naturally)
        IF(TPopsFile) THEN
            CALL WriteToPopsfileParOneArr()
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

        IF(tHistEnergies) CALL WriteHistogramEnergies()

        IF(tPrintOrbOcc) THEN
            CALL PrintOrbOccs(OrbOccs)
            DEALLOCATE(OrbOccs)
            CALL LogMemDeAlloc(this_routine,OrbOccsTag)
        ENDIF

        IF(tHistInitPops) THEN
            DEALLOCATE(HistInitPops)
            CALL LogMemDeAlloc(this_routine,HistInitPopsTag)
            IF(iProcIndex.eq.0) THEN
                DEALLOCATE(AllHistInitPops)
                CALL LogMemDeAlloc(this_routine,AllHistInitPopsTag)
            ENDIF
        ENDIF

! Print out some load balancing stats nicely to end.
        CALL MPI_Reduce(TotWalkers,MaxWalkers,1,MPI_INTEGER,MPI_MAX,root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(TotWalkers,MinWalkers,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,error)
        CALL MPI_AllReduce(Real(TotWalkers,dp),AllTotWalkers,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
        if (iProcIndex.eq.Root) then
            MeanWalkers=AllTotWalkers/nProcessors
            write (6,'(/,1X,a55)') 'Load balancing information based on the last iteration:'
            write (6,'(1X,a33,1X,f18.10)') 'Mean number of walkers/processor:',MeanWalkers
            write (6,'(1X,a32,1X,i18)') 'Min number of walkers/processor:',MinWalkers
            write (6,'(1X,a32,1X,i18,/)') 'Max number of walkers/processor:',MaxWalkers
        end if

!Deallocate memory
        CALL DeallocFCIMCMemPar()

        IF(iProcIndex.eq.Root) THEN
            CLOSE(15)
            IF(tTruncInitiator.or.tDelayTruncInit) CLOSE(16)
        ENDIF
        IF(TDebug) CLOSE(11)

        RETURN

    END SUBROUTINE FciMCPar


!This is the heart of FCIMC, where the MC Cycles are performed.
    SUBROUTINE PerformFCIMCycPar()
        USE FciMCLoggingMOD , only : TrackSpawnAttempts
        use GenRandSymExcitCSF, only: GenRandSymCSFExcit
        USE CalcData , only : tAddtoInitiator,InitiatorWalkNo,tInitIncDoubs
        use detbitops, only : countbits
        INTEGER :: VecSlot,i,j,k,l,CopySign,Loop,iPartBloom
        INTEGER :: nJ(NEl),ierr,IC,Child,iCount,DetCurr(NEl),iLutnJ(0:NIfTot)
        REAL*8 :: Prob,rat,HDiag,HDiagCurr
        INTEGER :: iDie,WalkExcitLevel,Proc             !Indicated whether a particle should self-destruct on DetCurr
        INTEGER :: ExcitLevel,TotWalkersNew,iGetExcitLevel_2,error,length,temp,Ex(2,2),WSign,p,Scratch1(ScratchSize),Scratch2(ScratchSize),Scratch3(Scratchsize),FDetSym,FDetSpin
        LOGICAL :: tParity,tMainArr,tFilled,TestClosedShellDet,tHFFound,tHFFoundTemp
        INTEGER(KIND=i2) :: HashTemp
        TYPE(HElement) :: HDiagTemp
        CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: message
        REAL :: Gap
        
        CALL set_timer(Walker_Time,30)

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



        tHFFound=.false.
        tHFFoundTemp=.false.
        IF(tHFInitiator) THEN
            Proc=DetermineDetProc(iLutHF)   !This wants to return a value between 0 -> nProcessors-1
            IF(iProcIndex.ne.Proc) THEN
!The processor with the HF determinant on it will have to check through each determinant until its found.            
                tHFFound=.true.
            ENDIF
        ENDIF
        IF(tFlipHighPopFound) CALL FlipHighPopDet()
        
!VecSlot indicates the next free position in NewDets
        VecSlot=1
!Reset number at HF and doubles
        NoatHF=0
        NoatDoubs=0
        iPartBloom=0
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
        ValidSpawnedList(:)=InitialSpawnedSlots(:)
        tHFFound=.false.
        tHFFoundTemp=.false.
        IF(iProcIndex.ne.iHFProc) tHFFound=.true.
!The processor with the HF determinant on it will have to check through each determinant until its found. Once found, tHFFound is true and it no longer needs to be checked.           
        CALL InitHistMin()


        do j=1,TotWalkers
!j runs through all current walkers
!If we are rotoannihilating/direct annihilating, the sign indicates the sum of the signs on the determinant, and hence j loops over determinants, not particles.
            !WRITE(6,*) Iter,j,TotWalkers
!            CALL FLUSH(6)

!First, decode the bit-string representation of the determinant the walker is on, into a string of naturally-ordered integers
!            WRITE(6,*) 'CurrentDet (bit)',CurrentDets(:,j)
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j))

!Also, we want to find out the excitation level - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.

            IF(tTruncSpace.or.tHistSpawn.or.tCalcFCIMCPsi.or.tHistHamil) THEN
!We need to know the exact excitation level for truncated calculations.
                WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j),nel)
            ELSE
                WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j),2)
            ENDIF

            IF(tRegenDiagHEls) THEN
!We are not storing the diagonal hamiltonian elements for each particle. Therefore, we need to regenerate them.
!Need to find H-element!
                IF(DetBitEQ(CurrentDets(0:NIfTot,j),iLutHF,NIfDBO).and.(.not.(tHub.and.tReal))) THEN
!We know we are at HF - HDiag=0
                    HDiagCurr=0.D0
                ELSE
                    if (tHPHF) then
                        HDiagtemp = hphf_diag_helement (DetCurr,CurrentDets(:,j))
                    else
                        HDiagTemp = get_helement (DetCurr, DetCurr, 0)
                    endif
                    HDiagCurr=(REAL(HDiagTemp%v,dp))-Hii
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
                    EXIT
                ENDIF
            ENDIF

            ! Sum in any energy contribution from the determinant, including 
            ! other parameters, such as excitlevel info.
            ! This is where the projected energy is calculated.
            call SumEContrib (DetCurr, WalkExcitLevel, CurrentSign(j),CurrentDets(:,j), HDiagCurr, 1.d0)

            IF(tTruncInitiator) CALL CalcParentFlag(j,tHFFound,tHFFoundTemp)

            tFilled=.false.     !This is for regenerating excitations from the same determinant multiple times. There will be a time saving if we can store the excitation generators temporarily.
            IF(tMCExcits) THEN
!Multiple spawning attempt per walker.
                Loop=abs(CurrentSign(j))*NoMCExcits
            ELSE
!Here, we spawn each particle on the determinant in a seperate attempt.
                Loop=abs(CurrentSign(j))
            ENDIF

            do p=1,Loop
!we are simply looping over all the particles on the determinant

                IF(tHPHF) THEN
!                    CALL GenRandHPHFExcit(DetCurr,CurrentDets(:,j),nJ,iLutnJ,pDoubles,exFlag,Prob)
                    CALL GenRandHPHFExcit2Scratch(DetCurr,CurrentDets(:,j),nJ,iLutnJ,pDoubles,exFlag,Prob,Scratch1,Scratch2,tFilled,tGenMatHEl)
                elseif (tCSF) then
                    ! This is a bit of a hack based on the fact
                    ! that we mean something different by exFlag
                    ! than the normal determinential code.
                    exFlag = 7
                    call GenRandSymCSFExcit (DetCurr, CurrentDets(:,j), nJ, pSingles, pDoubles, IC, Ex, exFlag, Prob, Scratch1, Scratch2, Scratch3, tFilled, tParity)
                else
                    CALL GenRandSymExcitScratchNU(DetCurr,CurrentDets(:,j),nJ,pDoubles,IC,Ex,tParity,exFlag,Prob,Scratch1,Scratch2,tFilled)
!                    WRITE(6,'(A,8I3)') 'determinant generated for spawning',nJ
                ENDIF

!Calculate number of children to spawn
                IF(IsNullDet(nJ)) THEN
                    Child=0
                ELSEIF(((TTruncSpace.or.tTruncCAS).and.(.not.tTruncInitiator)).or.tPartFreezeCore.or.tPartFreezeVirt.or.tFixLz.or.tUEG) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.
!If we are using the CASStar - all spawns are allowed so no need to check.
!                    WRITE(6,*) 'cheking if a spawn is allowed'
!                    WRITE(6,*) 'tTruncSpace',tTruncSpace
!                    WRITE(6,*) 'tTruncCAS',tTruncCAS
!                    WRITE(6,*) 'tTruncInitiator',tTruncInitiator

                    IF(CheckAllowedTruncSpawn(WalkExcitLevel,nJ,iLutnJ,IC)) THEN
                        Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity)
                    ELSE
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                        Child=0
                    ENDIF
                ELSE
!SD Space is not truncated - allow attempted spawn as usual

!                    WRITE(6,'(A,3I20,A,8I3)') 'Child to be created from:',CurrentDets(:,j),' which is ',DetCurr(:)
!                    WRITE(6,'(A,3I20,A,8I3)') 'Child to be created on:',iLutnJ(:),' which is ',nJ(:)
                    !print*, 'attempt create'
                    Child=AttemptCreatePar(DetCurr,CurrentDets(:,j),CurrentSign(j),nJ,iLutnJ,Prob,IC,Ex,tParity)
                    !print*, 'attempted'
                    !if (child /= 0) then
                        !WRITE(6,'(A,3I20)') 'Child to be created on:',iLutnJ(:)
                        !WRITE(6,*) 'Child',Child
                    !endif

                ENDIF

                IF(Child.ne.0) THEN


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
                    ENDIF

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

                    Acceptances=Acceptances+ABS(Child)      !Sum the number of created children to use in acceptance ratio
                
                ENDIF   !End if child created

            enddo   !End of cycling over mulitple particles on same determinant.

            IF(tHFFoundTemp) tHFFound=.true.
            tHFFoundTemp=.false.

!We now have to decide whether the parent particle (j) wants to self-destruct or not...
!For rotoannihilation, we can have multiple particles on the same determinant - these can be stochastically killed at the same time.

            iDie=AttemptDiePar(DetCurr,HDiagCurr,WalkExcitLevel,CurrentSign(j))
            NoDied=NoDied+iDie          !Update death counter

!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births
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
                    SpawnedParts(:,ValidSpawnedList(iProcIndex))=CurrentDets(:,j)
                    SpawnedSign(ValidSpawnedList(iProcIndex))=CopySign
                    ValidSpawnedList(iProcIndex)=ValidSpawnedList(iProcIndex)+1
                ENDIF
            ENDIF

!Finish cycling over walkers
        enddo

!SumWalkersCyc calculates the total number of walkers over an update cycle on each process
        SumWalkersCyc=SumWalkersCyc+(INT(TotParts,i2))
!        WRITE(6,*) "Born, Die: ",NoBorn, NoDied

!Since VecSlot holds the next vacant slot in the array, TotWalkers will be one less than this.
        TotWalkersNew=VecSlot-1

!Output if there has been a particle bloom this iteration. A negative number indicates that particles were created from a single excitation.
        IF((iPartBloom.ne.0).and.(iProcIndex.eq.0)) THEN
            IF(tAddtoInitiator.and.(iPartBloom.gt.0)) THEN
                WRITE(6,"(A,I14,A,I8,A)") "Particle Blooms of more than 'n_add' in iteration ",Iter," :  A max of ",abs(iPartBloom)," particles created in one attempt from double excit."
            ELSEIF(tAddtoInitiator) THEN
                WRITE(6,"(A,I14,A,I8,A)") "Particle Blooms of more than 'n_add' in iteration ",Iter," :  A max of ",abs(iPartBloom)," particles created in one attempt from single excit."
            ELSEIF(iPartBloom.gt.0) THEN
                WRITE(6,"(A,I14,A,I8,A)") "LARGE Particle Blooms in iteration ",Iter," :  A max of ",abs(iPartBloom)," particles created in one attempt from double excit."
            ELSE
                WRITE(6,"(A,I14,A,I8,A)") "LARGE Particle Blooms in iteration ",Iter," :  A max of ",abs(iPartBloom)," particles created in one attempt from single excit."
            ENDIF
        ENDIF

        rat=REAL(TotWalkersNew,dp)/REAL(MaxWalkersPart,dp)
        IF(rat.gt.0.95) THEN
            WRITE(6,'(A)') "*WARNING* - Number of particles/determinants has increased to over 95% of MaxWalkersPart"
            CALL FLUSH(6)
        ENDIF

!Need to test whether any of the sublists are getting to the end of their allotted space.
        IF(nProcessors.gt.1) THEN
            do i=0,nProcessors-1
                rat=(ValidSpawnedList(i)-InitialSpawnedSlots(i))/(InitialSpawnedSlots(1)+0.D0)
!                    WRITE(6,*) rat,(ValidSpawnedList(i)-InitialSpawnedSlots(i)),InitialSpawnedSlots(1)
                IF(rat.gt.0.95) THEN
                    WRITE(6,'(A)') "*WARNING* - Highest processor spawned particles has reached over 95% of MaxSpawned"
                    CALL FLUSH(6)
                ENDIF
            enddo
        ELSE
            rat=(ValidSpawnedList(0)+0.D0)/(MaxSpawned+0.D0)
            IF(rat.gt.0.9) THEN
                WRITE(6,'(A)') "*WARNING* - Number of spawned particles has reached over 90% of MaxSpawned"
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
        
    END SUBROUTINE PerformFCIMCycPar
                        

    SUBROUTINE CalcParentFlag(j,tHFFound,tHFFoundTemp)
        USE CalcData , only : tAddtoInitiator,InitiatorWalkNo,tInitIncDoubs
        USE FciMCLoggingMOD, only : InitBinMin,InitBinIter
        INTEGER :: j,WalkExcitLevel,InitBinNo
        LOGICAL :: tParentInCAS,tHFFound,tHFFoundTemp

        IF(tTruncSpace) THEN
            IF(WalkExcitLevel.le.ICILevel) THEN
                ParentInitiator=0
                NoInitDets=NoInitDets+1.D0
                NoInitWalk=NoInitWalk+(ABS(REAL(CurrentSign(j))))
!Parent in allowed space.                        
            ELSEIF(tInitIncDoubs.and.(WalkExcitLevel.eq.2)) THEN
                ParentInitiator=0
                NoInitDets=NoInitDets+1.D0
                NoInitWalk=NoInitWalk+(ABS(REAL(CurrentSign(j))))
                NoExtraInitDoubs=NoExtraInitDoubs+1.D0
            ELSEIF(tAddtoInitiator.and.(ABS(CurrentSign(j)).gt.InitiatorWalkNo)) THEN
                ParentInitiator=0
                IF(mod(Iter,StepsSft).eq.0) NoAddedInitiators=NoAddedInitiators+1.D0
                NoInitDets=NoInitDets+1.D0
                NoInitWalk=NoInitWalk+(ABS(REAL(CurrentSign(j))))
            ELSE
                ParentInitiator=1
                NoNonInitDets=NoNonInitDets+1.D0
                NoNonInitWalk=NoNonInitWalk+(ABS(REAL(CurrentSign(j))))
!Parent outside allowed space.                        
            ENDIF
        ELSE
!This it the case where the fixed initiator space is defined using the CAS notation, or where it is limited to only the HF determinant.                           
!I.e. when tTruncCAS or tHFInitiator are true.
!If neither of these are true, the expandspace option must have been used and so the parent will always be in the initiator space. 
!                            tParentInCAS=TestIfDetInCAS(DetCurr)
            IF(tTruncCAS) THEN
                tParentInCAS=TestIfDetInCASBit(CurrentDets(:,j))
            ELSEIF(tHFFound) THEN
!The HF determinant has been found, the rest will therefore all be out of the active space.                            
                tParentInCAS=.false.
            ELSEIF(tHFFoundTemp) THEN
!The HF determinant has been found and we are still running over the particles on the same determinant.                              
                tParentInCAS=.true.
            ELSEIF(tHFInitiator) THEN
                tParentInCAS=DetBitEQ(CurrentDets(:,j),iLutHF,NIfDBO)
                IF(tParentInCAS) tHFFoundTemp=.true.
!The temp parameter means we have found the HF determinant but we are still in the loop of running over the walkers on this determinant.
!At the end of this loop, tHFFound becomes true becuase the HF has been and gone and the rest of the determinants need not be checked.
            ENDIF
            IF(tParentInCAS) THEN
                ParentInitiator=0
                NoInitDets=NoInitDets+1.D0
                NoInitWalk=NoInitWalk+(ABS(REAL(CurrentSign(j))))
!The parent walker from which we are attempting to spawn is in the active space - all children will carry this flag, and these spawn like usual.
            ELSEIF(tInitIncDoubs.and.(WalkExcitLevel.eq.2)) THEN
                ParentInitiator=0
                NoInitDets=NoInitDets+1.D0
                NoInitWalk=NoInitWalk+(ABS(REAL(CurrentSign(j))))
                NoExtraInitDoubs=NoExtraInitDoubs+1.D0
            ELSEIF(tAddtoInitiator.and.(ABS(CurrentSign(j)).gt.InitiatorWalkNo)) THEN
                ParentInitiator=0
                IF(mod(Iter,StepsSft).eq.0) NoAddedInitiators=NoAddedInitiators+1.D0
                NoInitDets=NoInitDets+1.D0
                NoInitWalk=NoInitWalk+(ABS(REAL(CurrentSign(j))))
            ELSE
                ParentInitiator=1
                NoNonInitDets=NoNonInitDets+1.D0
                NoNonInitWalk=NoNonInitWalk+(ABS(REAL(CurrentSign(j))))
!The parent from which we are attempting to spawn is outside the active space - children spawned on unoccupied determinants with this flag will be killed.
            ENDIF
        ENDIF

        IF(tHistInitPops.and.(MOD(Iter,HistInitPopsIter).eq.0)) THEN
            IF((ParentInitiator.eq.0).and.((ABS(CurrentSign(j)).gt.InitiatorWalkNo))) THEN
!Just summing in those determinants which are initiators. 

!Need to figure out which bin to put them in though.
                IF(CurrentSign(j).lt.0) THEN
                    InitBinNo=(FLOOR(((log(REAL(ABS(CurrentSign(j)))))-InitBinMin)/InitBinIter))+1
                    IF((InitBinNo.ge.1).and.(InitBinNo.le.25000)) THEN
                        HistInitPops(1,InitBinNo)=HistInitPops(1,InitBinNo)+1
                    ELSE
                        CALL Stop_All('CalcParentFlag','Trying to histogram outside the range of the bins.')
                    ENDIF

                ELSE
                    InitBinNo=(FLOOR(((log(REAL(CurrentSign(j))))-InitBinMin)/InitBinIter))+1
 
                    IF((InitBinNo.ge.1).and.(InitBinNo.le.25000)) THEN
                        HistInitPops(2,InitBinNo)=HistInitPops(2,InitBinNo)+1
                    ELSE
                        CALL Stop_All('CalcParentFlag','Trying to histogram outside the range of the bins.')
                    ENDIF
                ENDIF
            ENDIF
        ENDIF

        IF(tFlipHighPopSign.and.(ParentInitiator.eq.0)) THEN
            IF(CurrentSign(j).lt.MaxInitPop) THEN
                MaxInitPop=CurrentSign(j)
                HighPopFlip=j
            ENDIF

            tFlipHighPopFound=.true.     
        ENDIF





    END SUBROUTINE CalcParentFlag                    
    

    SUBROUTINE FlipHighPopDet()
!Found the highest population on each processor, need to find out which of these has the highest of all.
        INTEGER :: MaxPops(nProcessors),error,i,MaxPop,MaxPopProc,InitDetCurr(NEl)


        WRITE(6,*) 'Highest populated determinant has pop',CurrentSign(HighPopflip)

!Need some sort of mpi to send these all into an array, hopefully so that processor 0 in position 0, then 1, etc etc.

        CALL MPI_Gather(CurrentSign(HighPopFlip),1,MPI_INTEGER,MaxPops,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)

        IF(iProcIndex.eq.Root) THEN
            WRITE(6,*) 'MaxPops',MaxPops

            MaxPop=MaxPops(1)
            do i=2,nProcessors
                IF(MaxPops(i).lt.MaxPop) THEN
                    MaxPop=MaxPops(i)
                    MaxPopProc=i-1
                ENDIF
            enddo
        ENDIF

        CALL MPI_Bcast(MaxPopProc,1,MPI_INTEGER,Root,MPI_COMM_WORLD,error)
        IF(error.ne.0) CALL Stop_All('FlipHighPopDet','Error in MPI_Bcast.')


        IF(iProcIndex.eq.MaxPopProc) THEN

            WRITE(6,*) 'The most highly populated determinant with the opposite sign to the HF has ',CurrentSign(HighPopFlip),' walkers.'
            WRITE(6,*) 'Flipping the sign of this determinant.'
            WRITE(6,*) 'In bit representation this is: ',CurrentDets(:,HighPopFlip)
            CALL DecodeBitDet(InitDetCurr,CurrentDets(:,HighPopFlip))
            WRITE(6,*) 'which has orbitals: ',InitDetCurr

            CurrentSign(HighPopFlip)=(-1)*CurrentSign(HighPopFlip)

        ENDIF

        tFlipHighPopFound=.false.
        tFlipHighPopSign=.false.


    END SUBROUTINE FlipHighPopDet


    
!Every StepsSft steps, update the diagonal shift value (the running value for the correlation energy)
!We don't want to do this too often, since we want the population levels to acclimatise between changing the shifts
    SUBROUTINE CalcNewShift()
        USE FciMCLoggingMOD , only : PrintSpawnAttemptStats,InitErrorBlocking,SumInErrorContrib
        USE FciMCLoggingMOD , only : InitShiftErrorBlocking,SumInShiftErrorContrib
        INTEGER :: error,rc,MaxAllowedWalkers,MaxWalkersProc,MinWalkersProc
        INTEGER :: inpair(6),outpair(6)
        REAL*8 :: TempTotWalkers,TempTotParts
        REAL*8 :: TempSumNoatHF,MeanWalkers,TempSumWalkersCyc,TempAllSumWalkersCyc
        REAL*8 :: inpairreal(3),outpairreal(3),inpairInit(9),outpairInit(9)
        LOGICAL :: tReZeroShift

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
            AllNoatHF=-AllNoatHF
            NoatHF=-NoatHF
        ENDIF
        

!This first call will calculate the GrowRate for each processor, taking culling into account
!        WRITE(6,*) "Get Here"
!        CALL FLUSH(6)
        CALL UpdateDiagSftPar()

!Put a barrier here so all processes synchronise
        CALL MPIBarrier(error)

!We need to collate the information from the different processors
!Inpair and outpair are used to package variables to save on latency time
!        inpair(1)=TotWalkers
        inpair(1)=Annihilated
        inpair(2)=NoatDoubs
        inpair(3)=NoBorn
        inpair(4)=NoDied
        inpair(5)=HFCyc         !SumWalkersCyc is now an integer*8
        inpair(6)=SpawnFromSing
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
        CALL MPI_Reduce(inpair,outpair,6,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get Here 2"
!        CALL FLUSH(6)
!        AllTotWalkers=outpair(1)
        AllAnnihilated=outpair(1)
        AllNoatDoubs=outpair(2)
        AllNoBorn=outpair(3)
        AllNoDied=outpair(4)
        AllHFCyc=REAL(outpair(5),dp)
        AllSpawnFromSing=outpair(6)
!        AllTotParts=outpair(7)
!        AlliUniqueDets=REAL(outpair(9),dp)
        TempTotWalkers=REAL(TotWalkers,dp)
        TempTotParts=REAL(TotParts,dp)

        IF(tTruncInitiator) THEN
            inpairInit(1)=NoAborted
            inpairInit(2)=NoAddedInitiators
            inpairInit(3)=NoInitDets
            inpairInit(4)=NoNonInitDets
            inpairInit(5)=NoInitWalk
            inpairInit(6)=NoNonInitWalk
            inpairInit(7)=NoDoubSpawns
            inpairInit(8)=NoExtraInitDoubs
            inpairInit(9)=InitRemoved
 
!!            CALL MPI_Reduce(inpairInit,outpairInit,8,MPI_INTEGER,MPI_SUM,Root,MPI_COMM_WORLD,error)
            Call MPIDSumArr(inpairInit,9,outpairInit)
!            CALL MPI_Reduce(inpairinit,outpairinit,8,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)

            AllNoAborted=outpairInit(1)
            AllNoAddedInitiators=outpairInit(2)
            AllNoInitDets=outpairInit(3)
            AllNoNonInitDets=outpairInit(4)
            AllNoInitWalk=outpairInit(5)
            AllNoNonInitWalk=outpairInit(6)
            AllNoDoubSpawns=outpairInit(7)
            AllNoExtraInitDoubs=outpairInit(8)
            AllInitRemoved=outpairInit(9)
        ENDIF

        CALL MPI_Reduce(TempTotWalkers,AllTotWalkers,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(TempTotParts,AllTotParts,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)

        IF(iProcIndex.eq.0) THEN
            IF(AllTotWalkers.le.0.2) THEN
                WRITE(6,*) AllTotWalkers,TotWalkers
                CALL Stop_All("CalcNewShift","All particles have died. Consider choosing new seed, or raising shift value.")
            ENDIF
        ENDIF

!SumWalkersCyc is now an int*8, therefore is needs to be reduced as a real*8
        TempSumWalkersCyc=REAL(SumWalkersCyc,dp)
        TempAllSumWalkersCyc=0.D0
        CALL MPI_Reduce(TempSumWalkersCyc,TempAllSumWalkersCyc,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)

!        WRITE(6,*) "Get Here 3"
!        CALL FLUSH(6)
        


        IF(TSinglePartPhase) THEN
!Exit the single particle phase if the number of walkers exceeds the value in the input file.
!            CALL MPI_Barrier(MPI_COMM_WORLD,error)
            IF((AllTotParts.gt.(InitWalkers*nProcessors)).or.(AllNoatHF.gt.MaxNoatHF)) THEN
                WRITE(6,*) "Exiting the single particle growth phase - shift can now change"
                VaryShiftIter=Iter
                TSinglePartPhase=.false.
            ENDIF
!!Broadcast the fact that TSinglePartPhase may have changed to all processors - unfortunatly, have to do this broadcast every iteration.
!            CALL MPILBcast(TSinglePartPhase,1,root)
        ELSE
!Exit the single particle phase if the number of walkers exceeds the value in the input file.
!            CALL MPI_Barrier(MPI_COMM_WORLD,error)
!            IF(iProcIndex.eq.root) THEN     !Only exit phase if particle number is sufficient on head node.
                IF(AllNoatHF.lt.(MaxNoatHF-HFPopThresh)) THEN
                    WRITE(6,'(A)') "No at HF has fallen too low - reentering the single particle growth phase - particle number may grow again."
!                    VaryShiftIter=Iter
                    TSinglePartPhase=.true.
                    tReZeroShift=.true.
                ENDIF
!            ENDIF
!!Broadcast the fact that TSinglePartPhase may have changed to all processors - unfortunatly, have to do this broadcast every iteration.
!            CALL MPILBcast(TSinglePartPhase,1,root)
        ENDIF


!Cannot load-balance with direct annihilation, but still want max & min
        CALL MPI_Reduce(TotWalkers,MaxWalkersProc,1,MPI_INTEGER,MPI_MAX,root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(TotWalkers,MinWalkersProc,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,error)
        IF(iProcIndex.eq.Root) THEN
            WalkersDiffProc=MaxWalkersProc-MinWalkersProc

            MeanWalkers=AllTotWalkers/REAL(nProcessors,dp)
            IF((WalkersDiffProc.gt.NINT(MeanWalkers/10.D0).and.(AllTotParts.gt.(REAL(nProcessors*500,dp))))) THEN
                WRITE(6,"(A62,F20.10,2i12)") "Number of determinants assigned to each processor unbalanced. ", (WalkersDiffProc*10.D0)/REAL(MeanWalkers),MinWalkersProc,MaxWalkersProc
            ENDIF
        
        ENDIF


!AlliUniqueDets corresponds to the total number of unique determinants, summed over all iterations in the last update cycle, and over all processors.
!Divide by StepsSft to get the average number of unique determinants visited over a single iteration.
!        AlliUniqueDets=AlliUniqueDets/(REAL(StepsSft,dp))
        
        IF(GrowRate.eq.-1.D0) THEN
!tGlobalSftCng is on, and we want to calculate the change in the shift as a global parameter, rather than as a weighted average.
!This will only be a sensible value on the root.
            AllGrowRate=AllTotParts/AllTotPartsOld
            IF(tTruncInitiator) AllGrowRateAbort=(AllTotParts+AllNoAborted)/(AllTotPartsOld+AllNoAbortedOld)
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
!            AllMeanExcitLevel=AllMeanExcitLevel/real(nProcessors,dp)
!        ENDIF

!AvSign no longer calculated (but would be easy to put back in) - ACF much better bet...
!        AvSign=AvSign/real(SumWalkersCyc,dp)
!        AvSignHFD=AvSignHFD/real(SumWalkersCyc,dp)
!        CALL MPI_Reduce(AvSign,AllAvSign,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(AvSignHFD,AllAvSignHFD,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        IF(iProcIndex.eq.Root) THEN
!            AllAvSign=AllAvSign/real(nProcessors,dp)
!            AllAvSignHFD=AllAvSignHFD/real(nProcessors,dp)
!        ENDIF

!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        TempSumNoatHF=real(SumNoatHF,dp)
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
                    AvDiagSft=SumDiagSft/REAL(VaryShiftCycles,dp)
                ENDIF

                IF(tTruncInitiator) THEN
                    DiagSftAbort=DiagSftAbort-(log(AllGrowRateAbort)*SftDamp)/(Tau*(StepsSft+0.D0))
                    IF((Iter-VaryShiftIter).ge.NShiftEquilSteps) THEN
                        SumDiagSftAbort=SumDiagSftAbort+DiagSftAbort
                        AvDiagSftAbort=SumDiagSftAbort/REAL(VaryShiftCycles,dp)
                    ENDIF
                ENDIF
            ENDIF

!Calculate the instantaneous value of the 'shift' from the HF population
            HFShift=-1.D0/REAL(AllNoatHF,dp)*(REAL(AllNoatHF-OldAllNoatHF,dp)/(Tau*REAL(StepsSft,dp)))
            InstShift=-1.D0/AllTotParts*((AllTotParts-AllTotPartsOld)/(Tau*REAL(StepsSft,dp)))

            IF(AllSumNoatHF.ne.0.D0) THEN
!AllSumNoatHF can actually be 0 if we have equilsteps on.
                ProjectionE=AllSumENum/AllSumNoatHF
            ENDIF

!Calculate the projected energy where each update cycle contributes the same weight to the average for its estimator for the energy
            IF(AllHFCyc.ne.0.D0) THEN
                ProjEIterSum=ProjEIterSum+(AllENumCyc/AllHFCyc)
                HFPopCyc=HFPopCyc+1   !This is the number of iterations where we have a non-zero contribution from HF particles
                ProjEIter=ProjEIterSum/REAL(HFPopCyc,dp)
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

        AccRat=(REAL(Acceptances,dp))/TempSumWalkersCyc      !The acceptance ratio which is printed is only for the current node - not summed over all nodes

        CALL WriteFCIMCStats()


!This first bit checks if it is time to set up the blocking analysis.  This is obviously only done once, so these logicals become false once it is done. 
        IF(iProcIndex.eq.Root) THEN
            IF(tIterStartBlock) THEN
!If IterStartBlocking is positive, then start blocking when we are at that iteration. Otherwise, wait until out of fixed shift.
                IF(IterStartBlocking.gt.0) THEN
                    IF(Iter.ge.IterStartBlocking) THEN 
                        CALL InitErrorBlocking(Iter)
                        tIterStartBlock=.false.
                        tErrorBlocking=.true.
                    ENDIF
                ELSE
                    IF(.not.TSinglePartPhase) THEN
                        CALL InitErrorBlocking(Iter)
                        tIterStartBlock=.false.
                        tErrorBlocking=.true.
                    ENDIF
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
            IF(tErrorBlocking.and.(.not.tBlockEveryIteration)) CALL SumInErrorContrib(Iter,AllENumCyc,AllHFCyc)
            IF(tShiftBlocking.and.(Iter.ge.(VaryShiftIter+IterShiftBlock))) CALL SumInShiftErrorContrib(Iter,DiagSft)
        ENDIF

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
        Acceptances=0
        NoBorn=0
        SpawnFromSing=0
        NoDied=0
        ENumCyc=0.D0
!        ProjEIter=0.D0     Do not want to rezero, since otherwise, if there are no particles at HF in the next update cycle, it will print out zero.
        HFCyc=0
        NoAborted=0.D0
        NoAddedInitiators=0.D0
        NoInitDets=0.D0
        NoNonInitDets=0.D0
        NoInitWalk=0.D0
        NoNonInitWalk=0.D0
        NoDoubSpawns=0.D0
        NoExtraInitDoubs=0.D0
        InitRemoved=0.D0

!Reset TotWalkersOld so that it is the number of walkers now
        TotWalkersOld=TotWalkers
        TotPartsOld=TotParts
!Save the number at HF to use in the HFShift
        OldAllNoatHF=AllNoatHF

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
        AllNoatHF=0
        AllNoatDoubs=0
        AllNoBorn=0
        AllSpawnFromSing=0
        AllNoDied=0
        AllNoAborted=0.D0
        AllNoAddedInitiators=0.D0
        AllNoInitDets=0.D0
        AllNoNonInitDets=0.D0
        AllNoInitWalk=0.D0
        AllNoNonInitWalk=0.D0
        AllNoDoubSpawns=0.D0
        AllNoExtraInitDoubs=0.D0
        AllInitRemoved=0.D0



        RETURN
    END SUBROUTINE CalcNewShift

    
!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalcPar()
        use FciMCLoggingMOD , only : InitHistInitPops
        use SystemData , only : tRotateOrbs
        use CalcData , only : InitialPart
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet
        INTEGER :: DetLT,VecSlot,error,MemoryAlloc,Proc
        TYPE(HElement) :: rh,TempHii
        LOGICAL :: exists
        REAL*8 :: TotDets
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMCPar'
            
        IF(TReadPops) THEN
!Read in particles from multiple POPSFILES for each processor
            WRITE(6,*) "Reading in initial particle configuration from POPSFILES..."

            CALL ReadFromPopsFilePar()

        ELSE
!initialise the particle positions - start at HF with positive sign

!Set the maximum number of walkers allowed
            MaxWalkersPart=NINT(MemoryFacPart*InitWalkers)
!            WRITE(6,"(A,F14.5)") "Memory Factor for walkers is: ",MemoryFacPart
            WRITE(6,"(A,I14)") " Memory allocated for a maximum particle number per node of: ",MaxWalkersPart
            MaxSpawned=NINT(MemoryFacSpawn*InitWalkers)
!            WRITE(6,"(A,F14.5)") "Memory Factor for arrays used for spawning is: ",MemoryFacSpawn
!            WRITE(6,"(A,I14)") "Memory allocated for a maximum particle number per node for spawning of: ",MaxSpawned

!Put a barrier here so all processes synchronise
            CALL MPI_Barrier(MPI_COMM_WORLD,error)
!Allocate memory to hold walkers
            ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),4,this_routine,WalkVecDetsTag,ierr)
            WalkVecDets(0:NIfTot,1:MaxWalkersPart)=0

            IF(.not.tRegenDiagHEls) THEN
                ALLOCATE(WalkVecH(MaxWalkersPart),stat=ierr)
                CALL LogMemAlloc('WalkVecH',MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
                WalkVecH(:)=0.d0
            ELSE
                WRITE(6,"(A,F14.6,A)") " Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*8,dp)/1048576.D0," Mb/Processor"
            ENDIF
            
            ALLOCATE(WalkVecSign(MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecSign',MaxWalkersPart,4,this_routine,WalkVecSignTag,ierr)
!            ALLOCATE(WalkVec2Sign(MaxWalkersPart),stat=ierr)
!            CALL LogMemAlloc('WalkVec2Sign',MaxWalkersPart,4,this_routine,WalkVec2SignTag,ierr)
            MemoryAlloc=(NIfTot+1+3)*MaxWalkersPart*4    !Memory Allocated in bytes

            IF(tRegenDiagHEls) THEN
                MemoryAlloc=MemoryAlloc-(MaxWalkersPart*8)
            ENDIF
            
            WalkVecSign(:)=0

            WRITE(6,"(A,I12,A)") " Spawning vectors allowing for a total of ",MaxSpawned," particles to be spawned in any one iteration."
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

!Allocate pointers to the correct walker arrays
            CurrentDets=>WalkVecDets
            CurrentSign=>WalkVecSign
            IF(.not.tRegenDiagHEls) THEN
                CurrentH=>WalkVecH
            ENDIF
        
            iHFProc=DetermineDetProc(iLutHF)   !This wants to return a value between 0 -> nProcessors-1
            WRITE(6,*) "HF processor is: ",iHFProc

!Setup initial walker local variables
            IF(iProcIndex.eq.iHFProc) THEN
                CurrentDets(:,1)=iLutHF(:)
                IF(.not.tRegenDiagHEls) CurrentH(1)=0.D0
                IF(TStartSinglePart) THEN
                    CurrentSign(1)=InitialPart
                    TotWalkers=1
                    TotWalkersOld=1
                    TotParts=InitialPart
                    TotPartsOld=InitialPart
                    NoatHF=InitialPart
                ELSE
                    CurrentSign(1)=InitWalkers
                    TotWalkers=1
                    TotWalkersOld=1
                    TotParts=InitWalkers
                    TotPartsOld=InitWalkers
                ENDIF

            ELSE
                IF(tStartSinglePart) THEN
                    NoatHF=0
                    TotWalkers=0
                    TotWalkersOld=0
                    TotParts=0
                    TotPartsOld=0
                ELSE
                    TotWalkers=0
                    TotWalkersOld=0
                    TotParts=0
                    TotPartsOld=0
                ENDIF
            ENDIF

        
            IF(TStartSinglePart) THEN
!Initialise global variables for calculation on the root node
                IF(iProcIndex.eq.root) THEN
                    OldAllNoatHF=InitialPart
                    AllNoatHF=InitialPart
                    AllTotWalkers=1.D0
                    AllTotWalkersOld=1.D0
                    AllTotParts=REAL(InitialPart,dp)
                    AllTotPartsOld=REAL(InitialPart,dp)
                    AllNoAbortedOld=0.D0
                ENDIF
            ELSE
!In this, only one processor has initial particles.
                IF(iProcIndex.eq.Root) THEN
                    AllTotWalkers=1.D0
                    AllTotWalkersOld=1.D0
                    AllTotParts=REAL(InitWalkers,dp)
                    AllTotPartsOld=REAL(InitWalkers,dp)
                    AllNoAbortedOld=0.D0
                ENDIF
            ENDIF
        
            WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens + temp arrays) consists of : ",REAL(MemoryAlloc,dp)/1048576.D0," Mb/Processor"
            WRITE(6,*) "Only one array of memory to store main particle list allocated..."
            WRITE(6,*) "Initial memory allocation sucessful..."
            CALL FLUSH(6)

        ENDIF   !End if initial walkers method
            
!Put a barrier here so all processes synchronise
        CALL MPI_Barrier(MPI_COMM_WORLD,error)

        IF(tTruncInitiator.or.tDelayTruncInit) THEN
            IF(tDelayTruncInit) tTruncInitiator=.false.
        ENDIF


        IF(tPrintOrbOcc) THEN
            ALLOCATE(OrbOccs(nBasis),stat=ierr)
            CALL LogMemAlloc('OrbOccs',nBasis,8,this_routine,OrbOccsTag,ierr)
            OrbOccs(:)=0.D0
        ENDIF

        IF(tHistInitPops) THEN
            CALL InitHistInitPops()
        ENDIF
        tFlipHighPopFound=.false.
        tFlipHighPopSign=.false.
        MaxInitPop=0

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
        
    

!This routine is the same as WriteToPopsfilePar, but does not require two main arrays to hold the data.
!The root processors data will be stored in a temporary array while it recieves the data from the other processors.
!This routine will write out to a popsfile. It transfers all walkers to the head node sequentially, so does not want to be called too often
    SUBROUTINE WriteToPopsfileParOneArr()
        use util_mod, only: get_unique_filename
        use CalcData, only: iPopsFileNoWrite
        use Logging, only: tIncrementPops
        REAL*8 :: TempSumNoatHF
        INTEGER :: error,WalkersonNodes(0:nProcessors-1)
        INTEGER :: Stat(MPI_STATUS_SIZE),Tag,Total,i,j,k
        INTEGER , ALLOCATABLE :: OrigSign(:), OrigParts(:,:)
        INTEGER :: OrigSignTag=0,OrigPartsTag=0
        CHARACTER(len=*) , PARAMETER :: this_routine='WriteToPopsfileParOneArr'
        character(255) :: popsfile

        CALL MPI_Barrier(MPI_COMM_WORLD,error)  !sync

!First, make sure we have up-to-date information - again collect AllTotWalkers,AllSumNoatHF and AllSumENum...
!        CALL MPI_Reduce(TotWalkers,AllTotWalkers,1,MPI_INTEGER,MPI_Sum,root,MPI_COMM_WORLD,error)    
!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        TempSumNoatHF=real(SumNoatHF,dp)
        CALL MPIDSumRoot(TempSumNoatHF,1,AllSumNoatHF,Root)
        CALL MPIDSumRoot(SumENum,1,AllSumENum,Root)

!We also need to tell the root processor how many particles to expect from each node - these are gathered into WalkersonNodes
        CALL MPI_AllGather(TotWalkers,1,MPI_INTEGER,WalkersonNodes,1,MPI_INTEGER,MPI_COMM_WORLD,error)

        Tag=125
!        WRITE(6,*) "Get Here"
!        CALL FLUSH(6)

        IF(iProcIndex.eq.root) THEN
!First, check that we are going to receive the correct number of particles...
            Total=0
            do i=0,nProcessors-1
                Total=Total+INT(WalkersonNodes(i)/iPopsPartEvery)
            enddo
            AllTotWalkers=REAL(Total,dp)
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
                call get_unique_filename('POPSFILEHEAD',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
            ELSE
                call get_unique_filename('POPSFILE',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
            ENDIF
            OPEN(17,FILE=popsfile,Status='replace')
            WRITE(17,*) AllTotWalkers,"   TOTWALKERS (all nodes)"
            WRITE(17,*) DiagSft,"   DIAG SHIFT"
            WRITE(17,*) NINT(AllSumNoatHF,i2),"   SUMNOATHF (all nodes)"
            WRITE(17,*) AllSumENum,"   SUMENUM ( \sum<D0|H|Psi> - all nodes)"
            WRITE(17,*) Iter+PreviousCycles,"   PREVIOUS CYCLES"
            IF(tBinPops) THEN
                CLOSE(17)
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
                OPEN(17,FILE=popsfile,Status='replace',form='unformatted')
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


!This routine reads in particle configurations from a POPSFILE.
    SUBROUTINE ReadFromPopsfilePar()
        use util_mod, only: get_unique_filename
        use CalcData, only: iPopsFileNoRead
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn
        use Logging, only: tIncrementPops,tZeroProjE
        LOGICAL :: exists,First,tBinRead
        INTEGER :: AvWalkers,WalkerstoReceive(nProcessors)
        INTEGER*8 :: NodeSumNoatHF(nProcessors),TempAllSumNoatHF
        REAL*8 :: TempTotParts
        INTEGER :: TempInitWalkers,error,i,j,k,l,total,ierr,MemoryAlloc,Tag,iLutTemp(0:NIfTot),TempSign,Proc,CurrWalkers
        INTEGER :: Stat(MPI_STATUS_SIZE),AvSumNoatHF,VecSlot,IntegerPart,HFPointer,TempnI(NEl),ExcitLevel,VecInd,DetsMerged
        REAL*8 :: r,FracPart,TempTotWalkers,Gap
        TYPE(HElement) :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        character(255) :: popsfile
        
        PreviousCycles=0    !Zero previous cycles
        SumENum=0.D0
        TotParts=0
        SumNoatHF=0
        DiagSft=0.D0
        Tag=124             !Set Tag
        
        call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
        INQUIRE(FILE=popsfile,EXIST=exists)
        IF(exists) THEN
            OPEN(17,FILE=popsfile,Status='old')
            tBinRead=.false.
        ELSE
            tBinRead=.true.
            call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            INQUIRE(FILE=popsfile,EXIST=exists)
            IF(.not.exists) THEN
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"No POPSFILEs of any kind found.")
                ELSE
                    CALL Stop_All(this_routine,"POPSFILEBIN(.x) found, but POPSFILEHEAD(.x) also needed for header information")
                ENDIF
            ELSE
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"POPSFILEHEAD(.x) found, but no POPSFILEBIN(.x) for particle information - this is also needed")
                ELSE
                    call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                    OPEN(17,FILE=popsfile,Status='old')
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
            call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            OPEN(17,FILE=popsfile,Status='old',form='unformatted')
        ENDIF

        IF(iProcIndex.eq.Root) THEN

            AllSumNoatHF=REAL(TempAllSumNoatHF,dp)
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
            AvWalkers=NINT(AllTotWalkers/real(nProcessors,dp))

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
!Now, let the total space allocated for storing walkers which have been read in to be equal to the initwalkers from the input file.
!                InitWalkers=AvWalkers
            ELSE
                TSinglePartPhase=.true.
            ENDIF
            SumENum=AllSumENum/REAL(nProcessors,dp)     !Divide up the SumENum over all processors
            AvSumNoatHF=NINT(AllSumNoatHF/real(nProcessors,dp)) !This is the average Sumnoathf
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

        IF(MemoryFacPart.le.1.D0) THEN
            WRITE(6,*) 'MemoryFacPart must be larger than 1.0 when reading in a POPSFILE - increasing it to 1.50.'
            MemoryFacPart=1.50
        ENDIF
        
!Now we want to allocate memory on all nodes.
        MaxWalkersPart=NINT(MemoryFacPart*(NINT(InitWalkers*ScaleWalkers)))   !InitWalkers here is simply the average number of walkers per node, not actual
        MaxSpawned=NINT(MemoryFacSpawn*(NINT(InitWalkers*ScaleWalkers)))

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
        ALLOCATE(WalkVecSign(MaxWalkersPart),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating WalkVecSign array.')
        CALL LogMemAlloc('WalkVecSign',MaxWalkersPart,4,this_routine,WalkVecSignTag,ierr)
        WalkVecSign(:)=0

!Just allocating this here, so that the SpawnParts arrays can be used for sorting the determinants when using direct annihilation.

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

!Allocate pointers to the correct walker arrays...
        CurrentDets=>WalkVecDets
        CurrentSign=>WalkVecSign

!Need to now allocate other arrays
        IF(.not.tRegenDiagHEls) THEN
            ALLOCATE(WalkVecH(MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecH',MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
            WalkVecH(:)=0.d0
        ELSE
            WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*4*2,dp)/1048576.D0," Mb/Processor"
        ENDIF

        IF(.not.tRegenDiagHEls) THEN
            CurrentH=>WalkVecH
            NewH=>WalkVec2H
        ENDIF
        NewDets=>WalkVec2Dets
        NewSign=>WalkVec2Sign

! The hashing will be different in the new calculation from the one where the POPSFILE was produced, this means we must recalculate the processor each determinant wants to go to.                
! This is done by reading in all walkers to the root and then distributing them in the same way as the spawning steps are done - by finding the determinant and sending it there.

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

        CALL MPI_Barrier(MPI_COMM_WORLD,error)  !Sync

        IF(iProcIndex.eq.root) WRITE(6,'(I10,A)') INT(AllTotWalkers,i2)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers by a factor of: ",ScaleWalkers

! CurrWalkers is the number of determinants on a particular node, AllTotWalkers is the total over all nodes.
            IntegerPart=INT(ScaleWalkers)
            FracPart=ScaleWalkers-REAL(IntegerPart)

            do l=1,CurrWalkers
                CurrentSign(l)=CurrentSign(l)*IntegerPart
                r = genrand_real2_dSFMT() 
                IF(r.lt.FracPart) THEN
!Stochastically create another particle
                    IF(CurrentSign(l).lt.0) THEN
                        CurrentSign(l)=CurrentSign(l)-1
                    ELSE
                        CurrentSign(l)=CurrentSign(l)+1
                    ENDIF
                ENDIF
            enddo

            InitWalkers=NINT(InitWalkers*ScaleWalkers)  !New (average) number of initial particles for culling criteria
!Other parameters don't change (I think) because the number of determinants isn't changing.                
            TotWalkers=CurrWalkers
            TotWalkersOld=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                    AllTotWalkers=TotWalkers
                AllTotWalkersOld=AllTotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,i2)
            ENDIF

        ELSE
!We are not scaling the number of walkers...

            TotWalkers=CurrWalkers
            TotWalkersOld=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                AllTotWalkersOld=AllTotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,i2)
            ENDIF

        ENDIF
            
        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft

        MemoryAlloc=((NIfTot+1+3)*MaxWalkersPart*4)
        IF(tRegenDiagHEls) THEN
            MemoryAlloc=MemoryAlloc-(MaxWalkersPart*4*2)
        ENDIF

        WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,dp)/1048576.D0," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        CALL FLUSH(6)

        WRITE(6,*) "Excitgens will be regenerated when they are needed..."
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
                    HFPointer=j
                    First=.false.
                ENDIF

            ELSE
                IF(.not.tRegenDiagHEls) THEN
                    if (tHPHF) then
                        HElemTemp = hphf_diag_helement (TempnI, &
                                                        CurrentDets(:,j))
                    else
                        HElemTemp = get_helement (TempnI, TempnI, 0)
                    endif
                    CurrentH(j)=REAL(HElemTemp%v,dp)-Hii
                ENDIF

            ENDIF
            TotParts=TotParts+abs(CurrentSign(j))

        enddo

        TempTotParts=REAL(TotParts,dp)

        CALL MPI_Barrier(MPI_COMM_WORLD,error)  !Sync
        CALL MPI_Reduce(TempTotParts,AllTotParts,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)

        IF(iProcIndex.eq.root) AllTotPartsOld=AllTotParts

        WRITE(6,'(A,F20.1)') ' The total number of particles read from the POPSFILE is: ',AllTotParts

        RETURN

    END SUBROUTINE ReadFromPopsfilePar


!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
    INTEGER FUNCTION AttemptCreatePar(DetCurr,iLutCurr,WSign,nJ,iLutnJ,Prob,IC,Ex,tParity)
        use GenRandSymExcitNUMod , only : GenRandSymExcitBiased
        use Logging, only : CCMCDebug
        INTEGER :: DetCurr(NEl),nJ(NEl),IC,StoreNumTo,StoreNumFrom,DetLT,i,ExtraCreate,Ex(2,2),WSign
        INTEGER :: iLutCurr(0:NIfTot),Bin,iLutnJ(0:NIfTot),PartInd,ExcitLev,iLut(0:NIfTot),iLut2(0:NIfTot)
        LOGICAL :: tParity,SymAllowed,tSuccess
        integer :: yama(NIfY)
        REAL*8 :: Prob,r,rat
        TYPE(HElement) :: rh,rhcheck

        IF(tMCExcits) THEN
!If we are generating multiple excitations, then the probability of spawning on them must be reduced by the number of excitations generated.
!This is equivalent to saying that the excitation is likely to arise a factor of NoMCExcits more often.
            Prob=Prob*REAL(NoMCExcits,dp)
        ENDIF
            

!Calculate off diagonal hamiltonian matrix element between determinants
!        rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
        IF(tHPHF) THEN
            IF(tGenMatHEl) THEN
!The prob is actually prob/HEl, since the matrix element was generated at the same time as the excitation

                rat=Tau/abs(Prob)

                rh%v=Prob ! to get the signs right for later on.
!                WRITE(6,*) Prob, DetCurr(:),"***",nJ(:)
!                WRITE(6,*) "******"
!                CALL HPHFGetOffDiagHElement(DetCurr,nJ,iLutCurr,iLutnJ,rh)

            ELSE
!The IC given doesn't really matter. It just needs to know whether it is a diagonal or off-diagonal matrix element.
!However, the excitation generator can generate the same HPHF again. If this is done, the routine will send the matrix element back as zero.
                rh = hphf_off_diag_helement (DetCurr, nJ, iLutCurr, iLutnJ)
!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
                rat=Tau*abs(rh%v)/Prob
!                WRITE(6,*) Prob/rh%v, DetCurr(:),"***",nJ(:)
!                WRITE(6,*) "******"

            ENDIF
        ELSE
!Normal determinant spawn

            rh = get_helement (DetCurr, nJ, IC, Ex, tParity)
            !WRITE(6,*) rh%v

!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
            rat=Tau*abs(rh%v)/Prob
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
        r = genrand_real2_dSFMT() 
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
!            WRITE(6,"(A,F15.5,I5,G25.15,I8,G25.15)") "Outwards ", rat,ExtraCreate,real(rh%v),Prob
!        ENDIF

        IF(tHistEnergies) THEN
!First histogram off-diagonal matrix elements.
            Bin=INT((real(rh%v,dp)+OffDiagMax)/OffDiagBinRange)+1
            IF(Bin.le.0.or.Bin.gt.iOffDiagNoBins) THEN
                CALL Stop_All("AttemptCreatePar","Trying to histogram off-diagonal matrix elements, but outside histogram array bounds.")
            ENDIF
            IF(IC.eq.1) THEN
                SinglesAttemptHist(Bin)=SinglesAttemptHist(Bin)+(Tau/Prob)
                IF(AttemptCreatePar.ne.0) THEN
                    SinglesHist(Bin)=SinglesHist(Bin)+real(abs(AttemptCreatePar),dp)
                    IF(BRR(Ex(1,1)).le.NEl) THEN
                        IF(BRR(Ex(2,1)).le.NEl) THEN
                            SinglesHistOccOcc(Bin)=SinglesHistOccOcc(Bin)+real(abs(AttemptCreatePar),dp)
                        ELSE
                            SinglesHistOccVirt(Bin)=SinglesHistOccVirt(Bin)+real(abs(AttemptCreatePar),dp)
                        ENDIF
                    ELSE
                        IF(BRR(Ex(2,1)).le.NEl) THEN
                            SinglesHistVirtOcc(Bin)=SinglesHistVirtOcc(Bin)+real(abs(AttemptCreatePar),dp)
                        ELSE
                            SinglesHistVirtVirt(Bin)=SinglesHistVirtVirt(Bin)+real(abs(AttemptCreatePar),dp)
                        ENDIF
                    ENDIF
                ENDIF
            ELSEIF(IC.eq.2) THEN
                DoublesAttemptHist(Bin)=DoublesAttemptHist(Bin)+(Tau/Prob)
                IF(AttemptCreatePar.ne.0) THEN
                    DoublesHist(Bin)=DoublesHist(Bin)+real(abs(AttemptCreatePar),dp)
                ENDIF
            ENDIF

            IF(tHPHF) THEN
                rh = hphf_diag_helement (nJ, iLutnJ)
            ELSE
                rh = get_helement (nJ, nJ, 0)
            ENDIF
            Bin=INT((real(rh%v,dp)-Hii)/BinRange)+1
            IF(Bin.gt.iNoBins) THEN
                CALL Stop_All("AttemptCreatePar","Histogramming energies higher than the arrays can cope with. Increase iNoBins or BinRange")
            ENDIF
            IF(AttemptCreatePar.ne.0) THEN
                SpawnHist(Bin)=SpawnHist(Bin)+real(abs(AttemptCreatePar),dp)
!                WRITE(6,*) "Get Here!", real(abs(AttemptCreatePar),dp),Bin
            ENDIF
            AttemptHist(Bin)=AttemptHist(Bin)+(Tau/Prob)
        ENDIF

    END FUNCTION AttemptCreatePar



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
!If there are multiple particles, decide how many to kill in total...
        rat=Tau*(Kii-DiagSft)*abs(WSign)

        iKill=INT(rat)
        rat=rat-REAL(iKill)

!Stochastically choose whether to die or not
        r = genrand_real2_dSFMT() 
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
   

!This routine will change the reference determinant to DetCurr. It will also re-zero all the energy estimators, since they now correspond to
!projection onto a different determinant.
    SUBROUTINE ChangeRefDet(HDiagCurr,DetCurr,iLutCurr)
        use Determinants , only : GetH0Element3
        INTEGER :: iLutCurr(0:NIfTot),DetCurr(NEl),i,nStore(6),ierr,iMaxExcit
        INTEGER :: iLutTemp(0:NIfTot)
        INTEGER :: nJ(NEl)
        TYPE(HElement) :: TempHii
        REAL*8 :: HDiagCurr

!        CALL Stop_All("ChangeRefDet","This option does not currently work. Bug ghb24 if its needed")
!Problem is that we need to rerun the simulation from scratch, and particles currently in the simulation will keep on
!changing the reference since their diagonal K element will remain negative.

        do i=1,NEl
            FDet(i)=DetCurr(i)
        enddo

        WRITE(6,"(A)") "*** Changing the reference determinant ***"
        WRITE(6,"(A)") "Switching reference and zeroing energy counters - restarting simulation"
!        
!Initialise variables for calculation on each node
        Iter=1
        
        CALL DeallocFCIMCMemPar()
        IF(iProcIndex.eq.Root) THEN
            CLOSE(15)
            IF(tTruncInitiator.or.tDelayTruncInit) CLOSE(16)
!            IF(TAutoCorr) CLOSE(44)
        ENDIF
        IF(TDebug) CLOSE(11)
        CALL SetupParameters()
        CALL InitFCIMCCalcPar()


    END SUBROUTINE ChangeRefDet
    

!This routine will find the largest weighted MP1 determinants, from which we can construct energy level splitting dependant on the sign.
    

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


!        CALL CheckOrdering(PartList(:,1:Length),CurrentSign(1:Length),Length,.true.)

    END SUBROUTINE SortCompressListswH


    SUBROUTINE DeallocFCIMCMemPar()
        INTEGER :: i,error,length,temp
        CHARACTER(len=*), PARAMETER :: this_routine='DeallocFciMCMemPar'
        CHARACTER(LEN=MPI_MAX_ERROR_STRING) :: message

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
        IF(.not.tRegenDiagHEls) THEN
            DEALLOCATE(WalkVecH)
            CALL LogMemDealloc(this_routine,WalkVecHTag)
        ENDIF
        DEALLOCATE(SpawnVec)
        CALL LogMemDealloc(this_routine,SpawnVecTag)
        DEALLOCATE(SpawnVec2)
        CALL LogMemDealloc(this_routine,SpawnVec2Tag)
        DEALLOCATE(SpawnSignVec)
        CALL LogMemDealloc(this_routine,SpawnSignVecTag)
        DEALLOCATE(SpawnSignVec2)
        CALL LogMemDealloc(this_routine,SpawnSignVec2Tag)
        
        DEALLOCATE(HFDet)
        CALL LogMemDealloc(this_routine,HFDetTag)
        DEALLOCATE(iLutHF)

        IF(ALLOCATED(SpinInvBrr)) THEN
            CALL LogMemDealloc(this_routine,SpinInvBRRTag)
            DEALLOCATE(SpinInvBRR)
        ENDIF
        IF(ALLOCATED(CoreMask)) THEN
            DEALLOCATE(CoreMask)
            DEALLOCATE(ExtMask)
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
        USE FciMCLoggingMOD , only : PrintSpawnAttemptStats,InitErrorBlocking,SumInErrorContrib
        USE FciMCLoggingMOD , only : InitShiftErrorBlocking,SumInShiftErrorContrib
        INTEGER :: error,rc,MaxAllowedWalkers,MaxWalkersProc,MinWalkersProc
        INTEGER :: inpair(6),outpair(6)
        REAL*8 :: TempTotWalkers,TempTotParts
        REAL*8 :: TempSumNoatHF,MeanWalkers,TempSumWalkersCyc,TempAllSumWalkersCyc
        REAL*8 :: inpairreal(3),outpairreal(3),inpairInit(8),outpairInit(8)

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
                
!We need to collate the information from the different processors
!Inpair and outpair are used to package variables to save on latency time
!        inpair(1)=TotWalkers
        inpair(1)=Annihilated
        inpair(2)=NoatDoubs
        inpair(3)=NoBorn
        inpair(4)=NoDied
        inpair(5)=HFCyc         !SumWalkersCyc is now an integer*8
        inpair(6)=SpawnFromSing
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
         call MPIISumArr(inpair,6,outpair)
!        WRITE(6,*) "Get Here 2"
!        CALL FLUSH(6)
!        AllTotWalkers=outpair(1)
        AllAnnihilated=outpair(1)
        AllNoatDoubs=outpair(2)
        AllNoBorn=outpair(3)
        AllNoDied=outpair(4)
        AllHFCyc=REAL(outpair(5),dp)
        AllSpawnFromSing=outpair(6)
!        AllTotParts=outpair(7)
!        AlliUniqueDets=REAL(outpair(9),dp)
        TempTotWalkers=REAL(TotWalkers,dp)
        TempTotParts=REAL(TotParts,dp)

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
            Call MPIDSumArr(inpairInit,8,outpairInit)
!            CALL MPI_Reduce(inpairinit,outpairinit,8,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)

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
        Call MPIDSum(TempTotWalkers,1,AllTotWalkers)
        Call MPIDSum(TempTotParts,1,AllTotParts)
      

        IF(iProcIndex.eq.0) THEN
            IF(AllTotWalkers.le.0.2) THEN
                WRITE(6,*) AllTotWalkers,TotWalkers
                CALL Stop_All("CalcNewShift","All particles have died. Consider choosing new seed, or raising shift value.")
            ENDIF
        ENDIF

!SumWalkersCyc is now an int*8, therefore is needs to be reduced as a real*8
        TempSumWalkersCyc=REAL(SumWalkersCyc,dp)
        TempAllSumWalkersCyc=0.D0
!!        CALL MPI_Reduce(TempSumWalkersCyc,TempAllSumWalkersCyc,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        Call MPIDSum(TempSumWalkersCyc,1,TempAllSumWalkersCyc)
!        WRITE(6,*) "Get Here 3"
!        CALL FLUSH(6)

!Cannot load-balance with direct annihilation, but still want max & min
!!            CALL MPI_Reduce(TotWalkers,MaxWalkersProc,1,MPI_INTEGER,MPI_MAX,root,MPI_COMM_WORLD,error)
!!            CALL MPI_Reduce(TotWalkers,MinWalkersProc,1,MPI_INTEGER,MPI_MIN,root,MPI_COMM_WORLD,error)
        CALL MPIIReduce(TotWalkers,1,MPI_MAX,MaxWalkersProc)
        CALL MPIIReduce(TotWalkers,1,MPI_MIN,MinWalkersProc)
        
        IF(iProcIndex.eq.Root) THEN
            WalkersDiffProc=MaxWalkersProc-MinWalkersProc
        ENDIF

!AlliUniqueDets corresponds to the total number of unique determinants, summed over all iterations in the last update cycle, and over all processors.
!Divide by StepsSft to get the average number of unique determinants visited over a single iteration.
!        AlliUniqueDets=AlliUniqueDets/(REAL(StepsSft,dp))
        
        IF(GrowRate.eq.-1.D0) THEN
!tGlobalSftCng is on, and we want to calculate the change in the shift as a global parameter, rather than as a weighted average.
!This will only be a sensible value on the root.
            AllGrowRate=AllTotParts/AllTotPartsOld
            IF(tTruncInitiator) AllGrowRateAbort=(AllTotParts+AllNoAborted)/(AllTotPartsOld+AllNoAbortedOld)
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
!            AllMeanExcitLevel=AllMeanExcitLevel/real(nProcessors,dp)
!        ENDIF

!AvSign no longer calculated (but would be easy to put back in) - ACF much better bet...
!        AvSign=AvSign/real(SumWalkersCyc,dp)
!        AvSignHFD=AvSignHFD/real(SumWalkersCyc,dp)
!        CALL MPI_Reduce(AvSign,AllAvSign,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        CALL MPI_Reduce(AvSignHFD,AllAvSignHFD,1,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
!        IF(iProcIndex.eq.Root) THEN
!            AllAvSign=AllAvSign/real(nProcessors,dp)
!            AllAvSignHFD=AllAvSignHFD/real(nProcessors,dp)
!        ENDIF

!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        TempSumNoatHF=real(SumNoatHF,dp)
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
                    AvDiagSft=SumDiagSft/REAL(VaryShiftCycles,dp)
                ENDIF

                IF(tTruncInitiator) THEN
                    DiagSftAbort=DiagSftAbort-(log(AllGrowRateAbort)*SftDamp)/(Tau*(StepsSft+0.D0))
                    IF((Iter-VaryShiftIter).ge.NShiftEquilSteps) THEN
                        SumDiagSftAbort=SumDiagSftAbort+DiagSftAbort
                        AvDiagSftAbort=SumDiagSftAbort/REAL(VaryShiftCycles,dp)
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
                ProjEIter=ProjEIterSum/REAL(HFPopCyc,dp)
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

        AccRat=(REAL(Acceptances,dp))/TempSumWalkersCyc      !The acceptance ratio which is printed is only for the current node - not summed over all nodes

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
        Acceptances=0
        NoBorn=0
        SpawnFromSing=0
        NoDied=0
        ENumCyc=0.D0
!        ProjEIter=0.D0     Do not want to rezero, since otherwise, if there are no particles at HF in the next update cycle, it will print out zero.
        HFCyc=0
        NoAborted=0.D0
        NoAddedInitiators=0.D0
        NoInitDets=0.D0
        NoNonInitDets=0.D0
        NoInitWalk=0.D0
        NoNonInitWalk=0.D0
        NoDoubSpawns=0.D0
        NoExtraInitDoubs=0.D0

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
        AllNoatHF=0
        AllNoatDoubs=0
        AllNoBorn=0
        AllSpawnFromSing=0
        AllNoDied=0
        AllNoAborted=0.D0
        AllNoAddedInitiators=0.D0
        AllNoInitDets=0.D0
        AllNoNonInitDets=0.D0
        AllNoInitWalk=0.D0
        AllNoNonInitWalk=0.D0
        AllNoDoubSpawns=0.D0
        AllNoExtraInitDoubs=0.D0



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
            IF(tTruncInitiator.or.tDelayTruncInit) THEN
                WRITE(6,"(A2,A10,A16,A10,A16,A12,3A11,3A17,2A10,A13,A12,A13)") "#","Step","Shift","WalkerCng","GrowRate","TotWalkers","Annihil","NoDied","NoBorn","Proj.E","Av.Shift","Proj.E.ThisCyc",&
&               "NoatHF","NoatDoubs","AccRat","UniqueDets","IterTime"

                WRITE(15,"(A2,A10,A16,A10,A16,A12,3A13,3A17,2A10,A13,A12,A13,A17,A13,A13)") "#","Step","Shift","WalkerCng","GrowRate","TotWalkers","Annihil","NoDied","NoBorn","Proj.E","Av.Shift",&
&               "Proj.E.ThisCyc","NoatHF","NoatDoubs","AccRat","UniqueDets","IterTime","FracSpawnFrmSing","WalkDiffProc","TotImagTime"

                WRITE(16,"(A2,A10,2A15,2A16,2A20,5A18)") "# ","Step","No Aborted","NoAddedtoInit","FracDetsInit","FracWalksInit","NoDoubSpawns","NoExtraDoubs","InstAbortShift","AvAbortShift",&
&               "NoInitDets","NoNonInitDets","InitRemoved"

            ELSE
                WRITE(6,"(A)") "       Step     Shift      WalkerCng    GrowRate       TotWalkers    Annihil    NoDied    NoBorn    Proj.E          Av.Shift     Proj.E.ThisCyc   NoatHF NoatDoubs      AccRat     UniqueDets     IterTime"
                WRITE(15,"(A)") "#     1.Step   2.Shift    3.WalkerCng  4.GrowRate     5.TotWalkers  6.Annihil  7.NoDied  8.NoBorn  9.Proj.E       10.Av.Shift"&
&              // " 11.Proj.E.ThisCyc  12.NoatHF 13.NoatDoubs  14.AccRat  15.UniqueDets  16.IterTime 17.FracSpawnFromSing  18.WalkersDiffProc  19.TotImagTime  20.ProjE.ThisIter "&
&              // " 21.HFInstShift  22.TotInstShift"
            
            ENDIF
        ENDIF

    END SUBROUTINE WriteFciMCStatsHeader

    SUBROUTINE WriteFCIMCStats()

        IF(iProcIndex.eq.root) THEN

            WRITE(15,"(I12,G16.7,I10,G16.7,I12,3I13,3G17.9,2I10,G13.5,I12,G13.5,G17.5,I13,G13.5,3G17.9)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,   &
  &                INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AvDiagSft,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),IterTime,   &
  &                REAL(AllSpawnFromSing)/REAL(AllNoBorn),WalkersDiffProc,TotImagTime,IterEnergy,HFShift,InstShift
            WRITE(6,"(I12,G16.7,I10,G16.7,I12,3I11,3G17.9,2I10,G13.5,I12,G13.5)") Iter+PreviousCycles,DiagSft,INT(AllTotParts-AllTotPartsOld,i2),AllGrowRate,    &
  &                INT(AllTotParts,i2),AllAnnihilated,AllNoDied,AllNoBorn,ProjectionE,AvDiagSft,AllENumCyc/AllHFCyc,AllNoatHF,AllNoatDoubs,AccRat,INT(AllTotWalkers,i2),IterTime

            IF(tTruncInitiator.or.tDelayTruncInit) THEN
                WRITE(16,"(I12,2G15.1,2G16.7,2G20.1,5G18.7)") Iter+PreviousCycles,AllNoAborted,AllNoAddedInitiators,(REAL(AllNoInitDets)/REAL(AllNoNonInitDets)),&
 &              (REAL(AllNoInitWalk)/REAL(AllNoNonInitWalk)),AllNoDoubSpawns,AllNoExtraInitDoubs,DiagSftAbort,AvDiagSftAbort,AllNoInitDets,AllNoNonInitDets,AllInitRemoved
            ENDIF

            
            CALL FLUSH(6)
            CALL FLUSH(15)
            
        ENDIF

        RETURN

    END SUBROUTINE WriteFCIMCStats


    SUBROUTINE SetupParameters()
        use SystemData, only : tUseBrillouin,iRanLuxLev,tSpn,tHPHFInts,tRotateOrbs,tNoBrillouin,tROHF,tFindCINatOrbs,nOccBeta,nOccAlpha,tUHF
        use SystemData, only : tFixLz,LzTot,BasisFN,tBrillouinsDefault
        USE dSFMT_interface , only : dSFMT_init
        use CalcData, only : tFCIMC
        use CalcData , only : tRandomiseHashOrbs
        use Calc, only : VirtCASorbs,OccCASorbs,G_VMC_Seed
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn,TauFactor,StepsSftImag
        use Determinants , only : GetH0Element3
        use SymData , only : nSymLabels,SymLabelList,SymLabelCounts
        use Logging , only : tTruncRODump
        use DetCalc, only : NMRKS,tagNMRKS,FCIDets
        use SymExcit3, only : CountExcitations3 
        use DetBitOps, only: CountBits
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet,HFDetTest(NEl),Seed,alpha,beta,symalpha,symbeta,endsymstate
        INTEGER :: DetLT,VecSlot,error,HFConn,MemoryAlloc,iMaxExcit,nStore(6),nJ(Nel),BRR2(nBasis),LargestOrb,nBits,HighEDet(NEl),iLutTemp(0:NIfTot)
        TYPE(HElement) :: rh,TempHii
        TYPE(BasisFn) HFSym
        REAL*8 :: TotDets,SymFactor,r
        CHARACTER(len=*), PARAMETER :: this_routine='SetupParameters'
        CHARACTER(len=12) :: abstr
        LOGICAL :: tSuccess,tFoundOrbs(nBasis),tTurnBackBrillouin,FoundPair
        REAL :: Gap
        INTEGER :: nSingles,nDoubles,HFLz,ChosenOrb

!        CALL MPIInit(.false.)       !Initialises MPI - now have variables iProcIndex and nProcessors
        WRITE(6,*) ""
        WRITE(6,*) "Performing Parallel FCIMC...."
        WRITE(6,*) ""
        
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
            if (tReadPops) then
                ! Restart calculation.  Append to stats file (if it exists).
                OPEN(15,file='FCIMCStats',status='unknown',access='append')
            else
                OPEN(15,file='FCIMCStats',status='unknown')
            end if
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
        IF(tDebug) WRITE(6,"(A,I5)") "Symmetry of reference determinant is: ",INT(HFSym%Sym%S,4)
        IF(tFixLz) THEN
            CALL GetLz(HFDet,NEl,HFLz)
            WRITE(6,"(A,I5)") "Ml value of reference determinant is: ",HFLz
            IF(HFLz.ne.LzTot) THEN
                CALL Stop_All("SetupParameters","Chosen reference determinant does not have the same Lz value as indicated in the input.")
            ENDIF
        ENDIF

!Do a whole lot of tests to see if we can use Brillouins theorem or not.
        IF(tBrillouinsDefault) CALL CheckforBrillouins() 
        

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
        ELSE
            IF(tDebug) WRITE(6,*) "Simply transferring this into a spin orbital representation."
        ENDIF
! From now on, the orbitals are contained in symlabellist2 and symlabelcounts2 rather than the original arrays.
! These are stored using spin orbitals.
!Assume that if we want to use the non-uniform random excitation generator, we also want to use the NoSpinSym full excitation generators if they are needed. 

        IF(tTruncInitiator.and.(.not.tTruncCAS).and.(.not.tTruncSpace)) THEN
!Have not defined an initiatial initiator space using the EXCITE keyword or TRUNCATECAS.  The initial initiator space then defaults to only the HF determinant.                
            tHFInitiator=.true.
        ELSE
            tHFInitiator=.false.
        ENDIF

!If using a CAS space truncation, write out this CAS space
        IF(tTruncCAS) THEN
            IF(tTruncInitiator) THEN
                WRITE(6,'(A)') " *********** INITIATOR METHOD IN USE ***********"
                WRITE(6,'(A)') " Fixed initiator space defined using the CAS method."
            ELSE
                WRITE(6,*) "Truncated CAS space detected. Writing out CAS space..."
            ENDIF
            WRITE(6,'(A,I2,A,I2,A)') " In CAS notation, (spatial orbitals, electrons), this has been chosen as: (",(OccCASOrbs+VirtCASOrbs)/2,",",OccCASOrbs,")"
            DO I=NEl-OccCASorbs+1,NEl
                WRITE(6,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(6,G1(BRR(I))%SYM,.FALSE.)
                WRITE(6,'(I4)',advance='no') G1(BRR(I))%Ml
                WRITE(6,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
            ENDDO
            WRITE(6,'(A)') " -------------------------------------------------------------------------------------------------"
            DO I=NEl+1,NEl+VirtCASOrbs
                WRITE(6,'(6I7)',advance='no') I,BRR(I),G1(BRR(I))%K(1), G1(BRR(I))%K(2),G1(BRR(I))%K(3), G1(BRR(I))%MS
                CALL WRITESYM(6,G1(BRR(I))%SYM,.FALSE.)
                WRITE(6,'(I4)',advance='no') G1(BRR(I))%Ml
                WRITE(6,'(2F19.9)')  ARR(I,1),ARR(BRR(I),2)
            ENDDO
        ELSEIF(tTruncInitiator.and.tHFInitiator) THEN
            WRITE(6,'(A)') " *********** INITIATOR METHOD IN USE ***********"
            WRITE(6,'(A /)') " Starting with only the HF determinant in the fixed initiator space."
        ENDIF
 
!Setup excitation generator for the HF determinant. If we are using assumed sized excitgens, this will also be assumed size.
        IF(tUEG.or.tHub) THEN
            exflag=2
        ELSE
            exflag=3
        ENDIF
!Count all possible excitations - put into HFConn
        CALL CountExcitations3(HFDet,exflag,nSingles,nDoubles)
        HFConn=nSingles+nDoubles


!Initialise random number seed - since the seeds need to be different on different processors, subract processor rank from random number
        Seed=abs(G_VMC_Seed-iProcIndex)
        WRITE(6,*) "Value for seed is: ",Seed
!Initialise...
        CALL dSFMT_init(Seed)
        
        IF(tRandomiseHashOrbs) THEN
            ALLOCATE(RandomHash(nBasis),stat=ierr)
            IF(ierr.ne.0) THEN
                CALL Stop_All(this_routine,"Error in allocating RandomHash")
            ENDIF
            RandomHash(:)=0
            IF(iProcIndex.eq.root) THEN
                do i=1,nBasis
                    FoundPair=.false.
                    do while(.not.FoundPair)
                        r = genrand_real2_dSFMT()
                        ChosenOrb=INT(nBasis*r*1000)+1
                        do j=1,nBasis
                            IF(RandomHash(j).eq.ChosenOrb) EXIT
                        enddo
                        IF(j.eq.nBasis+1) THEN
                            RandomHash(i)=ChosenOrb
                            FoundPair=.true.
                        ELSE
                            FoundPair=.false.
                        ENDIF
                    enddo
                enddo

!                WRITE(6,*) "Random Orbital Indexing for hash:"
!                WRITE(6,*) RandomHash(:)
                do i=1,nBasis
                    IF((RandomHash(i).eq.0).or.(RandomHash(i).gt.nBasis*1000)) THEN
                        CALL Stop_All(this_routine,"Random Hash incorrectly calculated")
                    ENDIF
                    do j=i+1,nBasis
                        IF(RandomHash(i).eq.RandomHash(j)) THEN
                            CALL Stop_All(this_routine,"Random Hash incorrectly calculated")
                        ENDIF
                    enddo
                enddo
            ENDIF
            !Now broadcast to all processors
            CALL MPIIBCast(RandomHash,nBasis,Root)
        ENDIF

        IF(tHPHF) THEN
            !IF(tLatticeGens) CALL Stop_All("SetupParameters","Cannot use HPHF with model systems currently.")
            IF(tROHF.or.(LMS.ne.0)) CALL Stop_All("SetupParameters","Cannot use HPHF with high-spin systems.")
            tHPHFInts=.true.
        ENDIF

!Calculate Hii
        IF(tHPHF) THEN
            TempHii = hphf_diag_helement (HFDet, iLutHF)
        ELSE
            TempHii = get_helement (HFDet, HFDet, 0)
        ENDIF
        Hii=REAL(TempHii%v,dp)
        WRITE(6,*) "Reference Energy set to: ",Hii
        TempHii=GetH0Element3(HFDet)
        Fii=REAL(TempHii%v,dp)

!Find the highest energy determinant...
        IF(.not.tSpn) THEN
            do i=1,NEl
                HighEDet(i)=Brr(nBasis-(i-1))
            enddo
            IF(tHPHF) THEN
                call EncodeBitDet (HighEDet, iLutTemp)
                TempHii = hphf_diag_helement (HighEDet, iLutTemp)
            ELSE
                TempHii = get_helement (HighEDet, HighEDet, 0)
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

!Initialise variables for calculation on each node
        IterTime=0.0
        ProjectionE=0.D0
        AvSign=0.D0
        AvSignHFD=0.D0
        SumENum=0.D0
        SumNoatHF=0
        NoatHF=0
        NoatDoubs=0
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
        VaryShiftCycles=0
        AvDiagSft=0.D0
        SumDiagSft=0.D0
        SumDiagSftAbort=0.D0
        AvDiagSftAbort=0.D0
        NoAborted=0.D0
        NoAddedInitiators=0.D0
        NoInitDets=0.D0
        NoNonInitDets=0.D0
        NoInitWalk=0.D0
        NoNonInitWalk=0.D0
        NoDoubSpawns=0.D0
        NoExtraInitDoubs=0.D0
        InitRemoved=0.D0
        TotImagTime=0.D0
        HFIter=0
        ENumIter=0.D0
        IterEnergy=0.D0

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
        AllAnnihilated=0
        AllENumCyc=0.D0
        AllHFCyc=0.D0
!        AllDetsNorm=0.D0
        CullInfo(1:10,1:3)=0
        NoCulls=0
        AllNoAborted=0.D0
        AllNoAddedInitiators=0.D0
        AllNoInitDets=0.D0
        AllNoNonInitDets=0.D0
        AllNoInitWalk=0.D0
        AllNoNonInitWalk=0.D0
        AllNoDoubSpawns=0.D0
        AllNoExtraInitDoubs=0.D0
        AllInitRemoved=0.D0
 

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
        IF(tUseBrillouin) THEN
            WRITE(6,*) "Brillouin theorem specified, but this will not be in use with the non-uniform excitation generators."
        ENDIF
        WRITE(6,*) "Non-uniform excitation generators in use."
        CALL CalcApproxpDoubles(HFConn)
        IF(TauFactor.ne.0.D0) THEN
            WRITE(6,*) "TauFactor detected. Resetting Tau."
            Tau=TauFactor/REAL(HFConn,dp)
            WRITE(6,*) "Tau set to: ",Tau
        ENDIF
        IF(StepsSftImag.ne.0.D0) THEN
            WRITE(6,*) "StepsShiftImag detected. Resetting StepsShift."
            StepsSft=NINT(StepsSftImag/Tau)
            IF(StepsSft.eq.0) StepsSft=1
            WRITE(6,*) "StepsShift set to: ",StepsSft
        ENDIF
!        IF(tConstructNOs) THEN
!! This is the option for constructing the natural orbitals actually during a NECI calculation.  This is different (and probably a lot more complicated and doesn't 
!! currently work) from the FINDCINATORBS option which finds the natural orbitals given a final wavefunction.
!            ALLOCATE(OneRDM(nBasis,nBasis),stat=ierr)
!            CALL LogMemAlloc('OneRDM',nBasis*nBasis,8,this_routine,OneRDMTag,ierr)
!            OneRDM(:,:)=0.D0
!        ENDIF

        IF(TPopsFile) THEN
            IF(mod(iWritePopsEvery,StepsSft).ne.0) CALL Warning(this_routine,"POPSFILE writeout should be a multiple of the update cycle length.")
        ENDIF

        WRITE(6,*) "*Direct Annihilation* in use...Explicit load-balancing disabled."
        WRITE(6,*) ""
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

        IF(TReadPops) THEN
            IF(TStartSinglePart) CALL Stop_All("SetupParameters","ReadPOPS cannot work with StartSinglePart")
        ENDIF

        IF(.not.TReadPops) THEN
            WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft
        ENDIF
        WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp
        WRITE(6,*) "Maximum connectivity of HF determinant is: ",HFConn
        IF(TStartSinglePart) THEN
            TSinglePartPhase=.true.
            IF(TReadPops) THEN
                CALL Stop_All("SetupParameters","Cannot read in POPSFILE as well as starting with a single particle")
            ENDIF
        ELSE
            TSinglePartPhase=.false.
        ENDIF
        
        IF(ICILevel.ne.0) THEN
!We are truncating the excitations at a certain value
            TTruncSpace=.true.
            WRITE(6,'(A,I4)') "Truncating the S.D. space at determinants will an excitation level w.r.t. HF of: ",ICILevel
        ENDIF
        IF(tTruncCAS) THEN
!We are truncating the FCI space by only allowing excitations in a predetermined CAS space.
!The following line has already been written out if we are doing a CAS calculation.
!            WRITE(6,'(A,I4,A,I5)') "Truncating the S.D. space as determinants must be within a CAS of ",OccCASOrbs," , ",VirtCASOrbs
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

!Create the bit masks for the bit calculation of these properties.
            ALLOCATE(ExtMask(0:NIfD))
            ALLOCATE(CoreMask(0:NIfD))
            ExtMask(:)=0
            CoreMask(:)=0
            do i=1,nBasis
                IF(SpinInvBRR(i).gt.CASmax) THEN
                    !Orbital is in external space
                    ExtMask((SpinInvBRR(i)-1)/32) = ibset(ExtMask((i-1)/32),MOD((i-1),32))
!                    CoreMask = ibclr(CoreMask((i-1)/32),MOD((i-1),32))
                ELSEIF(SpinInvBRR(i).le.CASmin) THEN
                    !Orbital is in core space
!                    ExtMask = ibclr(ExtMask((i-1)/32),MOD((i-1),32))
                    CoreMask((SpinInvBRR(i)-1)/32) = ibset(CoreMask((i-1)/32),MOD((i-1),32))
!                    WRITE(6,*) "Setting orbital: ", SpinInvBRR(i)
                ELSE
                    !Orbital is in CAS space - these bits should already be cleared
!                    ExtMask = ibclr(ExtMask((i-1)/32),MOD((i-1),32))
!                    CoreMask = ibclr(CoreMask((i-1)/32),MOD((i-1),32))
                ENDIF
            enddo

!            do i=1,nBasis
!                IF(BTEST(CoreMask((i-1)/32),MOD((i-1),32))) THEN
!                    WRITE(6,"(I4)",advance='no') 1
!                ELSE
!                    WRITE(6,"(I4)",advance='no') 0
!                ENDIF
!            enddo
!            WRITE(6,*) ""

        ENDIF
        IF(tPartFreezeCore) THEN
            WRITE(6,'(A,I4,A,I5)') 'Partially freezing the lowest ',NPartFrozen,' spin orbitals so that no more than ',NHolesFrozen,' holes exist within this core.'
            CALL CreateSpinInvBRR()
        ENDIF
        IF(tPartFreezeVirt) THEN
            WRITE(6,'(A,I4,A,I5)') 'Partially freezing the highest ',NVirtPartFrozen,' virtual spin orbitals so that no more than ',NElVirtFrozen,' electrons occupy these orbitals.'
            CALL CreateSpinInvBRR()
        ENDIF

        SymFactor=(Choose(NEl,2)*Choose(nBasis-NEl,2))/(HFConn+0.D0)
        TotDets=1.D0
        do i=1,NEl
            WRITE(6,*) "Approximate excitation level population: ",i,NINT((Choose(NEl,i)*Choose(nBasis-NEl,i))/SymFactor)
            TotDets=TotDets+(Choose(NEl,i)*Choose(nBasis-NEl,i))/SymFactor
        enddo
        WRITE(6,*) "Approximate size of determinant space is: ",NINT(TotDets)

        IF(TStartSinglePart) THEN
            WRITE(6,"(A,F9.3,A,I9)") " Initial number of particles set to 1, and shift will be held at ",DiagSft," until particle number on root node gets to ",InitWalkers
        ELSE
            WRITE(6,*) " Initial number of walkers per processor chosen to be: ", InitWalkers
        ENDIF
 

    END SUBROUTINE SetupParameters

    SUBROUTINE CheckforBrillouins()
        use SystemData, only : tUseBrillouin,tNoBrillouin,tUHF
        use Determinants, only : tDefineDet
        INTEGER :: i,j
        LOGICAL :: tSpinPair
        
!Standard cases.
        IF((tHub.and.tReal).or.(tRotatedOrbs).or.((LMS.ne.0).and.(.not.tUHF))) THEN
!Open shell, restricted.            
            tNoBrillouin=.true.
        ELSE
!Closed shell restricted, or open shell unrestricted are o.k.            
            tNoBrillouin=.false.
            tUseBrillouin=.true.
        ENDIF

!Special case of complex orbitals.        
        IF(tFixLz.and.(.not.tNoBrillouin)) THEN
            WRITE(6,*) "Turning Brillouins theorem off since we are using non-canonical complex orbitals"
            tNoBrillouin=.true.
        ENDIF

!Special case of defining a det with LMS=0, but which is open shell. No Brillouins if it's a restricted HF calc.
        IF(tDefineDet.and.(LMS.eq.0).and.(.not.tUHF)) THEN
!If we are defining our own reference determinant, we want to find out if it is open shell or closed to know whether or not brillouins theorem holds.            
!If LMS/=0, then it is easy and must be open shell, otherwise we need to consider the occupied orbitals.
            do i=1,(NEl-1),2
!Assuming things will probably go alpha beta alpha beta, run through each alpha and see if there's a corresponding beta.
                tSpinPair=.false.
                IF(MOD(BRR(FDet(i)),2).ne.0) THEN
!Odd energy, alpha orbital.                    
                    IF(BRR(FDet(i+1)).ne.(BRR(FDet(i))+1)) THEN
!Check the next orbital to see if it's the beta (will be alpha+1 when ordered by energy). 
!If not, check the other orbitals for the beta, as it's possible the orbitals are ordered weird (?).
                        do j=1,NEl
                            IF(BRR(FDet(j)).eq.(BRR(FDet(i))+1)) tSpinPair=.true. 
                        enddo
                    ELSE
                        tSpinPair=.true.
                    ENDIF
                ELSE
!Even energy, beta orbital. The corresponding alpha will be beta-1.                    
                    IF(BRR(FDet(i+1)).ne.(BRR(FDet(i))-1)) THEN
                        do j=1,NEl
                            IF(BRR(FDet(j)).eq.(BRR(FDet(i))-1)) tSpinPair=.true. 
                        enddo
                    ELSE
                        tSpinPair=.true.
                    ENDIF
                ENDIF
                IF(.not.tSpinPair) EXIT
            enddo
            IF(.not.tSpinPair) THEN
!Open shell LMS=0 determinant.
!If restricted HF orbitals are being used, brillouins theorem does not hold.
                tNoBrillouin=.true.
                tUseBrillouin=.false.
                WRITE(6,'(A)') " Using an open shell reference determinant in a basis of restricted HF orbitals; Brillouins theorem is being turned off. "
            ENDIF
        ENDIF

    ENDSUBROUTINE CheckforBrillouins

    LOGICAL FUNCTION TestifDETinCASBit(iLutnI)
        INTEGER :: iLutnI(0:NIfD),i,Ext(0:NIfD),Core(0:NIfD)

        Ext(:) = iand(iLutnI(:),ExtMask(:))

        do i=0,NIfD
            IF(Ext(i).ne.0) THEN
                TestifDETinCASBit=.false.
                RETURN
            ENDIF
        enddo

        Core(:) = iand(iLutnI(:),CoreMask(:))

        do i=0,NIfD
            IF(Core(i).ne.CoreMask(i)) THEN
                TestifDETinCASBit=.false.
                RETURN
            ENDIF
        enddo

        TestifDETinCASBit=.true.

        RETURN
    END FUNCTION TestifDETinCASBit

    LOGICAL FUNCTION TestifDETinCAS(CASDet)
        INTEGER :: k,z,CASDet(NEl), orb
        LOGICAL :: tElecInVirt, bIsCsf

!        CASmax=NEl+VirtCASorbs
! CASmax is the max spin orbital number (when ordered energetically) within the chosen active space.
! Spin orbitals with energies larger than this maximum value must be unoccupied for the determinant
! to be in the active space.
!        CASmin=NEl-OccCASorbs   (These have been moved to the InitCalc subroutine so they're not calculated
! each time.
! CASmin is the max spin orbital number below the active space.  As well as the above criteria, spin 
! orbitals with energies equal to, or below that of the CASmin orbital must be completely occupied for 
! the determinant to be in the active space.

        bIsCsf = iscsf(CASDet)

        z=0
        tElecInVirt=.false.
        do k=1,NEl      ! running over all electrons
            ! TODO: is it reasonable to just apply the orbital mask anyway?
            !       it is probably faster than running iscsf...
            if (bIsCsf) then
                orb = iand(CASDet(k), csf_orbital_mask)
            else
                orb = CASDet(k)
            endif

            if (SpinInvBRR(orb).gt.CASmax) THEN
                tElecInVirt=.true.
                EXIT            
! if at any stage an electron has an energy greater than the CASmax value, the determinant can be ruled out
! of the active space.  Upon identifying this, it is not necessary to check the remaining electrons.
            else
                if (SpinInvBRR(orb).le.CASmin) THEN
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
!            IF(.not.TestIfDetinCAS(nJ)) THEN
            IF(.not.TestIfDetinCASBit(iLutnJ)) THEN
!Excitation not in allowed CAS space.
                CheckAllowedTruncSpawn=.false.
            ENDIF

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

        IF(tUEG.and.(.not.tLatticeGens)) THEN
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
        use CalcData , only : SinglesBias
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

!NSing=Number singles from HF, nDoub=No Doubles from HF

        WRITE(6,"(A)") " Calculating approximate pDoubles for use with excitation generator by looking a excitations from HF."
        exflag=3
        CALL CountExcitations3(HFDet,exflag,nSing,nDoub)
        iTotal=nSing + nDoub + ncsf

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
            pDoubles = real(nDoub,dp) / &
                   ((real(nSing,dp)*SinglesBias)+real(nDoub,dp)+real(ncsf,dp))
            pSingles = real(nSing,dp) / &
                   ((real(nSing,dp)*SinglesBias)+real(nDoub,dp)+real(ncsf,dp))

        else
            pDoubles = real(nDoub,dp) / &
                   ((real(NSing,dp)*SinglesBias) + real(NDoub,dp))
            pSingles = real(nSing,dp) * SinglesBias/ &
                   ((real(nSing,dp)*SinglesBias) + real(nDoub,dp))
        endif

        IF(SinglesBias.ne.1.D0) THEN
            write (6, '("pDoubles set to ", f14.6, &
                       &" rather than (without bias): ", f14.6)') &
                       pDoubles, real(nDoub,dp) / real(iTotal,dp)
            write (6, '("pSingles set to ", f14.6, &
                       &" rather than (without bias): ", f14.6)') &
                       pSingles, real(nSing,dp) / real(iTotal,dp)

            WRITE(6,"(A,F14.6,A,F14.6)") "pDoubles set to: ",pDoubles, " rather than (without bias): ",real(nDoub,dp)/real(iTotal,dp)
        ELSE
            write (6,'(A,F14.6)') " pDoubles set to: ", pDoubles
            write (6,'(A,F14.6)') " pSingles set to: ", pSingles
        ENDIF

        WRITE(6,'(A,F15.10)') " Assuming an average K_ij magnitude of approx 0.01, an appropriate tau is predicted to be around: ",(0.02*(1.D0/(REAL(NSing)+REAL(NDoub))))/0.01
!This is a rough guesstimate of what tau might like to be, assuming K_ij is approx 0.01 on average, and we want a probability of spawning to be about 0.02.        
!These are just stats taken from one system... will investigate further...

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

!This will store all the double excitations.
    SUBROUTINE StoreDoubs()
        use SystemData , only : tUseBrillouin
        use SymExcit3 , only : CountExcitations3,GenExcitations3
        INTEGER :: iMaxExcit,nStore(6),ExcitLength,nJ(NEl),ierr,iExcit,VecSlot,nSingles,ExcitMat3(2,2)
        LOGICAL :: tAllExcitFound,tParity

        IF(tUseBrillouin) THEN
            CALL Stop_All("StoreDoubs","Cannot have Brillouin theorem as now storing singles too...")
        ENDIF
        
!NoDoubs here is actually the singles + doubles of HF
        exflag=3
        CALL CountExcitations3(HFDet,exflag,nSingles,NoDoubs)
        NoDoubs=nSingles+NoDoubs

        ALLOCATE(DoublesDets(NEl,NoDoubs),stat=ierr)
        CALL LogMemAlloc('DoublesDets',NoDoubs*NEl,4,"StoreDoubs",DoublesDetsTag,ierr)
        DoublesDets(1:NEl,1:NoDoubs)=0
        
        VecSlot=1           !This is the next free slot in the DoublesDets array

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

!This means that now NoDoubs is double excitations AND singles
!        NoDoubs=VecSlot-1

        IF(VecSlot.ne.(NoDoubs+1)) THEN
            WRITE(6,*) VecSlot,NoDoubs
            CALL Stop_All("StoreDoubs","Problem enumerating all double excitations")
        ENDIF

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

!        MeanExcitLevel=MeanExcitLevel+real(ExcitLevel,dp)
!        IF(MinExcitLevel.gt.ExcitLevel) MinExcitLevel=ExcitLevel
!        IF(MaxExcitLevel.lt.ExcitLevel) MaxExcitLevel=ExcitLevel
!        DetsNorm=DetsNorm+REAL((WSign**2),dp)
        IF(ExcitLevel.eq.0) THEN
            IF(Iter.gt.NEquilSteps) SumNoatHF=SumNoatHF+WSign
            NoatHF=NoatHF+WSign
            HFCyc=HFCyc+WSign      !This is simply the number at HF*sign over the course of the update cycle 
            HFIter=HFIter+WSign
!            AvSign=AvSign+REAL(WSign,dp)
!            AvSignHFD=AvSignHFD+REAL(WSign,dp)
            
        ELSEIF(ExcitLevel.eq.2) THEN
            NoatDoubs=NoatDoubs+abs(WSign)
!At double excit - find and sum in energy
            IF(tHPHF) THEN
                HOffDiag = hphf_off_diag_helement (HFDet, DetCurr, iLutHF, &
                                                   iLutCurr)
            ELSE
                HOffDiag = get_helement (HFDet, DetCurr, ExcitLevel, iLutHF, &
                                         iLutCurr)
            ENDIF
            IF(Iter.gt.NEquilSteps) SumENum=SumENum+(REAL(HOffDiag%v,dp)*WSign/dProbFin)
!            AvSign=AvSign+REAL(WSign,dp)
!            AvSignHFD=AvSignHFD+REAL(WSign,dp)
            ENumCyc=ENumCyc+(REAL(HOffDiag%v,dp)*WSign/dProbFin)     !This is simply the Hij*sign summed over the course of the update cycle
            ENumIter=ENumIter+(REAL(HOffDiag%v,dp)*WSign/dProbFin)
!            WRITE(6,*) 2,SumENum,(REAL(HOffDiag%v,dp)*WSign/dProbFin)     !This is simply the Hij*sign summed over the course of the update cycle

            
            
!        ELSE
!            AvSign=AvSign+REAL(WSign,dp)

        ELSEIF(ExcitLevel.eq.1) THEN
          if(tNoBrillouin) then
!For the real-space hubbard model, determinants are only connected to excitations one level away, and brillouins theorem can not hold.
!For Rotated orbitals, brillouins theorem also cannot hold, and energy contributions from walkers on singly excited determinants must
!be included in the energy values along with the doubles.
            IF(tHPHF) THEN
                HOffDiag = hphf_off_diag_helement (HFDet, DetCurr, iLutHF, &
                                                   iLutCurr)
            ELSE
                HOffDiag = get_helement (HFDet, DetCurr, ExcitLevel, ilutHF, &
                                         iLutCurr)
            ENDIF
            IF(Iter.gt.NEquilSteps) SumENum=SumENum+(REAL(HOffDiag%v,dp)*WSign/dProbFin)
!            AvSign=AvSign+REAL(WSign,dp)
!            AvSignHFD=AvSignHFD+REAL(WSign,dp)
            ENumCyc=ENumCyc+(REAL(HOffDiag%v,dp)*WSign/dProbFin)     !This is simply the Hij*sign summed over the course of the update cycle
!            WRITE(6,*) 1,SumENum,(REAL(HOffDiag%v,dp)*WSign/dProbFin)     !This is simply the Hij*sign summed over the course of the update cycle
          endif 

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
                            Histogram(PartInd)=Histogram(PartInd)-(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                        ELSE
                            Histogram(PartInd)=Histogram(PartInd)+(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                        ENDIF
                    ELSE
                        IF(tFlippedSign) THEN
                            Histogram(PartInd)=Histogram(PartInd)-REAL(WSign,dp)/dProbFin
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-REAL(WSign,dp)/dProbFin
                        ELSE
                            Histogram(PartInd)=Histogram(PartInd)+REAL(WSign,dp)/dProbFin
                            IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+REAL(WSign,dp)/dProbFin
                        ENDIF
                    ENDIF
                ELSE
                    IF(tFlippedSign) THEN
                        Histogram(PartInd)=Histogram(PartInd)-REAL(WSign,dp)/dProbFin
                        IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-REAL(WSign,dp)/dProbFin
                    ELSE
                        Histogram(PartInd)=Histogram(PartInd)+REAL(WSign,dp)/dProbFin
                        IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+REAL(WSign,dp)/dProbFin
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
                                    Histogram(PartInd)=Histogram(PartInd)+(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                                ELSE
                                    Histogram(PartInd)=Histogram(PartInd)-(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                                ENDIF
                            ELSE
                                IF(mod(OpenOrbs,2).eq.1) THEN
                                    Histogram(PartInd)=Histogram(PartInd)-(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)-(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                                ELSE
                                    Histogram(PartInd)=Histogram(PartInd)+(REAL(WSign,dp)/SQRT(2.0))/dProbFin
                                    IF(tHistSpawn) InstHist(PartInd)=InstHist(PartInd)+(REAL(WSign,dp)/SQRT(2.0))/dProbFin
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
            Histogram(Bin)=Histogram(Bin)+real(abs(WSign),dp)
        ENDIF

        IF(tPrintOrbOcc.and.(Iter.ge.StartPrintOrbOcc)) THEN
            do i=1,NEl
                OrbOccs(DetCurr(i))=OrbOccs(DetCurr(i))+(WSign*WSign)
            enddo
        ENDIF

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
    
