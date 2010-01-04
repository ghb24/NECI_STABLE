!This module is to be used for various types of walker MC annihilation in serial and parallel.
MODULE AnnihilationMod
    use SystemData , only : NEl,tMerTwist,tHPHF,NIfTot,NIfDBO
    use CalcData , only : TRegenExcitgens,tAnnihilatebyRange,tUseGuide,tRegenDiagHEls,iInitGuideParts,iGuideDets,tKeepDoubleSpawns
    USE DetCalc , only : Det,FCIDetIndex
    USE Logging , only : tHistSpawn
    USE Parallel
    USE mt95 , only : genrand_real2
    USE FciMCData
    use DetBitOps, only: DetBitEQ, DetBitLT, FindBitExcitLevel
    use CalcData , only : tTruncInitiator
    IMPLICIT NONE

    contains

!This is a new annihilation algorithm. In this, determinants are kept on predefined processors, and newlyspawned particles are sent here so that all the annihilations are
!done on a predetermined processor, and not rotated around all of them.
    SUBROUTINE DirectAnnihilation(TotWalkersNew)
        integer :: i
        INTEGER :: MaxIndex,TotWalkersNew
!        WRITE(6,*) "Direct annihilation"
!        CALL FLUSH(6)

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
!        WRITE(6,*) "Transferred"
!        CALL FLUSH(6)
        CALL CompressSpawnedList(MaxIndex)

!        WRITE(6,*) "List compressed",MaxIndex,TotWalkersNew
!        CALL FLUSH(6)

!Binary search the main list and copy accross/annihilate determinants which are found.
!This will also remove the found determinants from the spawnedparts lists.

        CALL AnnihilateSpawnedParts(MaxIndex,TotWalkersNew)

!        WRITE(6,*) "Annihilation finished",MaxIndex,TotWalkersNew
!        CALL FLUSH(6)

!Put the surviving particles in the main list, maintaining order of the main list.
!Now we insert the remaining newly-spawned particles back into the original list (keeping it sorted), and remove the annihilated particles from the main list.
        CALL set_timer(Sort_Time,30)
        CALL InsertRemoveParts(MaxIndex,TotWalkersNew)

!       WRITE(6,*) "Surviving particles merged"
!       CALL FLUSH(6)

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
        
        CALL MPIAlltoAllI(sendcounts,1,recvcounts,1,error)

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
        CALL MPIAlltoAllvI(SpawnedSign(1:MaxSendIndex),sendcounts,disps,SpawnedSign2(1:MaxIndex),recvcounts,recvdisps,error)
        
!        WRITE(6,*) MaxIndex, "Recieved signs: "
!        do i=1,MaxIndex
!            WRITE(6,*) SpawnedSign2(i)
!        enddo

!Update the number of integers we need to send.
        do i=1,nProcessors
            sendcounts(i)=sendcounts(i)*(NIfTot+1)
            disps(i)=disps(i)*(NIfTot+1)
            recvcounts(i)=recvcounts(i)*(NIfTot+1)
            recvdisps(i)=recvdisps(i)*(NIfTot+1)
        enddo

!        WRITE(6,*) "Sent Particles: ", NINT(Gap),sendcounts(2)
!        do i=NINT(Gap)+1,NINT(Gap)+sendcounts(2)
!            write(6,*) i, '***', CountBits(spawnedparts(:,i), nifd)
!            WRITE(6,*) i,"***",SpawnedParts(:,i)
!        enddo
#ifdef PARALLEL
        CALL MPI_AlltoAllv(SpawnedParts(:,1:MaxSendIndex),sendcounts,disps,MPI_INTEGER,SpawnedParts2,recvcounts,recvdisps,MPI_INTEGER,MPI_COMM_WORLD,error)
#else
        SpawnedParts2(0:NIfTot,1:MaxIndex)=SpawnedParts(0:NIfTot,1:MaxSendIndex)
#endif

!        WRITE(6,*) MaxIndex, "Recieved particles: "
!        do i=1,MaxSpawned
!            IF(SpawnedParts2(1,i).ne.0) THEN
!            write(6,*) i, '***', CountBits(spawnedparts(:,i), nifd)
!                WRITE(6,*) SpawnedParts2(:,i)
!            ENDIF
!        enddo
        
        CALL halt_timer(Comms_Time)

    END SUBROUTINE SendProcNewParts

!This sorts and compresses the spawned list to make it easier for the rest of the annihilation process.
!This is not essential, but should proove worthwhile
    SUBROUTINE CompressSpawnedList(ValidSpawned)
        INTEGER :: VecInd,ValidSpawned,DetsMerged,ToRemove,i,SignProd,PartIndex,ExcitLevel
        LOGICAL :: tSuc

!We want to sort the list of newly spawned particles, in order for quicker binary searching later on. (this is not essential, but should proove faster)
!They should remain sorted after annihilation between spawned
        CALL SortBitDets(ValidSpawned,SpawnedParts(:,1:ValidSpawned), &
                         SpawnedSign(1:ValidSpawned))
        IF(tHistSpawn) HistMinInd2(1:NEl)=FCIDetIndex(1:NEl)

!First, we compress the list of spawned particles, so that they are only specified at most once in each processors list.
!During this, we transfer the particles over to SpawnedParts2
        IF(ValidSpawned.gt.0) THEN
            SpawnedParts2(0:NIfTot,1)=SpawnedParts(0:NIfTot,1)
            SpawnedSign2(1)=SpawnedSign(1)
        ENDIF
        VecInd=1
        DetsMerged=0
        ToRemove=0
        do i=2,ValidSpawned
            IF(.not.DetBitEQ(SpawnedParts(0:NIfTot,i),SpawnedParts2(0:NIfTot,VecInd),NIfDBO)) THEN
                IF(SpawnedSign2(VecInd).eq.0) ToRemove=ToRemove+1
                VecInd=VecInd+1
                SpawnedParts2(:,VecInd)=SpawnedParts(:,i)
                SpawnedSign2(VecInd)=SpawnedSign(i)
            ELSE
!The next determinant is equal to the current - want to look at the relative signs.                
                SignProd=SpawnedSign(i)*SpawnedSign2(VecInd)
                IF(SignProd.lt.0) THEN
!We are actually unwittingly annihilating, but just in serial... we therefore need to count it anyway.
                    Annihilated=Annihilated+2*(MIN(abs(SpawnedSign2(VecInd)),abs(SpawnedSign(i))))

                    IF(tTruncInitiator) THEN
!If we are doing a CAS star calculation, we also want to keep track of which parent the remaining walkers came from - those inside the active space or out.                
!This is only an issue if the two determinants we are merging have different Parent flags - otherwise they just keep whichever.
!As it is, the SpawnedParts2 determinant will have the parent flag that remains - just need to change this over if the number of walkers on SpawnedParts ends up dominating.
                        IF(SpawnedParts(NIfTot,i).ne.SpawnedParts2(NIfTot,VecInd)) THEN     ! Parent flags are not equal
                            IF(ABS(SpawnedSign(i)).gt.ABS(SpawnedSign2(VecInd))) SpawnedParts2(NIfTot,VecInd)=SpawnedParts(NIfTot,i)
                        ENDIF
                    ENDIF

                    IF(tHistSpawn) THEN
!We want to histogram where the particle annihilations are taking place.
                        ExcitLevel = FindBitExcitLevel(SpawnedParts(:,i), &
                                                       iLutHF, nel)
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
!                            WRITE(6,*) "***",SpawnedParts(0:NIfTot,i)
!                            CALL DecodeBitDet(TempDet,SpawnedParts(0:NIfTot,i))
!                            WRITE(6,*) "Full Det is: ",TempDet(:)
!                            IF(tHPHF) THEN
!                                CALL FindExcitBitDetSym(SpawnedParts(0:NIfD,i),iLutSym(:))
!                                WRITE(6,*) "*** Sym: ",iLutSym(:)
!                                CALL DecodeBitDet(TempDet,iLutSym(0:NIfTot))
!                                WRITE(6,*) "Full Sym Det is: ",TempDet(:)
!                            ENDIF
                            CALL Stop_All("CompressSpawnedList","Cannot find corresponding FCI determinant when histogramming")
                        ENDIF
                    ENDIF
                ELSEIF(tTruncInitiator) THEN
!This is the case where the determinants are the same but also have the same sign - so this usually doesn't matter except when we are doing CASStar calculations and 
!the parents are different.
!In this case we assume the determinants inside the CAS have spawned a second earlier - so the ones from outside the active space are spawning onto an occupied determinant
!and will therefore live - we can just make these equiv by treating them as they've all come from inside the active space.
                    IF(SpawnedParts(NIfTot,i).ne.SpawnedParts2(NIfTot,VecInd)) THEN     ! Parent flags are not equal
                        SpawnedParts2(NIfTot,VecInd)=0      ! Take all the walkers to have come from inside the CAS space.
                        IF(SpawnedSign2(VecInd).eq.0) SpawnedParts2(NIfTot,VecInd)=SpawnedParts(NIfTot,i) 
                        ! Think there might still be a case where SpawnedSign2 can be 0 - this means that the parent will be determined by SpawnedParts.
                        ! If its SpawnedParts that is 0 that's fine because the SpawnedParts2 flag is already carried across.
                    ELSEIF(tKeepDoubleSpawns.and.(SpawnedParts2(NIfTot,VecInd).eq.1)) THEN
!This is the option where if two determinants spawn onto another at the same time with the same sign, they are kept whether they've come from 
!inside or outside the active space.  This is different from before where two children spawned on the same determinant with the same sign, but both from outside the active
!space will be killed.
                        SpawnedParts2(NIfTot,VecInd)=0
                        NoDoubSpawns=NoDoubSpawns+1
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
                SpawnedParts(0:NIfTot,i-DetsMerged)=SpawnedParts2(0:NIfTot,i)
                SpawnedSign(i-DetsMerged)=SpawnedSign2(i)
            ENDIF
        enddo
        IF(DetsMerged.ne.ToRemove) THEN
            CALL Stop_All("CompressSpawnedList","Wrong number of entries removed from spawned list")
        ENDIF
        ValidSpawned=ValidSpawned-DetsMerged
        
    END SUBROUTINE CompressSpawnedList
    
    INTEGER FUNCTION DetermineDetProc(iLut)
        use systemdata , only: NIfDBO
        INTEGER :: iLut(0:NIfTot),i,j,Elecs!,TempDet(NEl),MurmurHash2Wrapper
        INTEGER(KIND=i2) :: Summ!,RangeofBins,NextBin

!        CALL DecodeBitDet(TempDet,iLut)
!        i=MurmurHash2Wrapper(TempDet,NEl,13)
!        write(6,*) i
        
        Summ=0
        Elecs=0
        lp2: do i=0,NIfDBO
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
        CALL MPIBarrier(error)
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

        CALL MPIBarrier(error)

!Allocate a buffer here to hold particles when using a buffered send...
!The buffer wants to be able to hold (MaxSpawned+1)x(NIfD+2) integers (*4 for in bytes). If we could work out the maximum ValidSpawned accross the determinants,
!it could get reduced to this... 
        IF(nProcessors.ne.1) THEN
            ALLOCATE(mpibuffer(8*(MaxSpawned+1)*(NIfTot+2)),stat=ierr)
            IF(ierr.ne.0) THEN
                CALL Stop_All("RotoAnnihilation","Error allocating memory for transfer buffers...")
            ENDIF
#ifdef PARALLEL
            CALL MPI_Buffer_attach(mpibuffer,8*(MaxSpawned+1)*(NIfTot+2),error)
#endif
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

#ifdef PARALLEL
!Detach buffers
            CALL MPI_Buffer_detach(mpibuffer,8*(MaxSpawned+1)*(NIfTot+2),error)
#endif
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


    
    SUBROUTINE AnnihilateBetweenSpawnedOneProc(ValidSpawned)
        use DetBitOps, only: DecodeBitDet
        INTEGER :: ValidSpawned,DetCurr(0:NIfTot),i,j,k,LowBound,HighBound,WSign
        INTEGER :: VecSlot,TotSign

        CALL SortBitDets(ValidSpawned,SpawnedParts(:,1:ValidSpawned), &
                         SpawnedSign(1:ValidSpawned))

        VecSlot=1
        i=1
        do while(i.le.ValidSpawned)
            LowBound=i
            DetCurr(0:NIfTot)=SpawnedParts(0:NIfTot,i)
            i=i+1
            do while(DetBitEQ(DetCurr(0:NIfTot),SpawnedParts(0:NIfTot,i),NIfDBO).and.(i.le.ValidSpawned))
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
                    SpawnedParts2(0:NIfTot,VecSlot)=DetCurr(0:NIfTot)
                    SpawnedSign2(VecSlot)=WSign
                    VecSlot=VecSlot+1
                enddo
            ENDIF

        enddo

        ValidSpawned=VecSlot-1

        do i=1,ValidSpawned
            SpawnedParts(0:NIfTot,i)=SpawnedParts2(0:NIfTot,i)
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
        use DetBitOps, only: DecodeBitDet
        INTEGER :: TotWalkersNew,ValidSpawned
        INTEGER :: i,DetsMerged,nJ(NEl)
        REAL*8 :: HDiag
        TYPE(HElement) :: HDiagTemp

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
                    CurrentDets(0:NIfTot,i-DetsMerged)=CurrentDets(0:NIfTot,i)
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
                CALL MergeLists(TotWalkersNew,MaxWalkersPart,ValidSpawned,SpawnedParts(0:NIfTot,1:ValidSpawned),SpawnedSign(1:ValidSpawned))
            ENDIF
        ELSE
            IF(TotWalkersNew.eq.0) THEN
!Merging algorithm will not work with no determinants in the main list.
                TotWalkersNew=ValidSpawned
                do i=1,ValidSpawned
                    CurrentDets(:,i)=SpawnedParts(:,i)
                    CurrentSign(i)=SpawnedSign(i)
!We want to calculate the diagonal hamiltonian matrix element for the new particle to be merged.
                    IF(DetBitEQ(CurrentDets(0:NIfTot,i),iLutHF,NIfDBO)) THEN
!We know we are at HF - HDiag=0
                        HDiag=0.D0
                    ELSE
                        CALL DecodeBitDet(nJ,CurrentDets(0:NIfTot,i))
                        IF(tHPHF) THEN
                            CALL HPHFGetDiagHElement(nJ,CurrentDets(0:NIfTot,i),HDiagTemp)
                        ELSE
                            HDiagTemp=GetHElement3(nJ,nJ,0)
                        ENDIF
                        HDiag=(REAL(HDiagTemp%v,8))-Hii
                    ENDIF
                    CurrentH(i)=HDiag
                enddo
            ELSE
                CALL MergeListswH(TotWalkersNew,MaxWalkersPart,ValidSpawned,SpawnedParts(0:NIfTot,1:ValidSpawned),SpawnedSign(1:ValidSpawned))
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
!                CompParts=DetBitLT(NewDets(0:NIfTot,IndParts),SpawnedParts(0:NIfD,IndSpawned))
!                IF(CompParts.eq.1) THEN
!!Want to move in the particle from NewDets (unless it wants to be annihilated)
!                    IF(NewSign(IndParts).ne.0) THEN
!!We want to keep this particle
!                        CurrentDets(0:NIfTot,VecInd)=NewDets(0:NIfTot,IndParts)
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
!                    CurrentDets(0:NIfTot,VecInd)=NewDets(0:NIfTot,IndParts)
!                    IF(.not.tRegenDiagHEls) CurrentH(VecInd)=NewH(IndParts)
!                    CurrentSign(VecInd)=NewSign(IndParts)+SpawnedSign(IndSpawned)
!                    IndParts=IndParts+1
!                    IndSpawned=IndSpawned+1
!
!                    do while(DetBitEQ(SpawnedParts(:,IndSpawned-1),SpawnedParts(:,IndSpawned)).and.(IndSpawned.le.ValidSpawned))
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
!                    CurrentDets(0:NIfTot,VecInd)=SpawnedParts(0:NIfTot,IndSpawned)
!                    CurrentSign(VecInd)=SpawnedSign(IndSpawned)
!                    IndSpawned=IndSpawned+1
!                    
!                    do while(DetBitEQ(SpawnedParts(:,IndSpawned-1),SpawnedParts(:,IndSpawned)).and.(IndSpawned.le.ValidSpawned))
!                        CurrentSign(VecInd)=CurrentSign(VecInd)+SpawnedSign(IndSpawned)
!                        IndSpawned=IndSpawned+1
!                    enddo
!                    
!                    TotParts=TotParts+abs(CurrentSign(VecInd))
!
!                    IF(.not.tRegenDiagHEls) THEN
!!Need to find H-element!
!                        IF(DetBitEQ(CurrentDets(0:NIfTot,VecInd),iLutHF)) THEN
!!We know we are at HF - HDiag=0
!                            HDiag=0.D0
!                            IF(tHub.and.tReal) THEN
!!Reference determinant is not HF
!                                CALL DecodeBitDet(nJ,CurrentDets(0:NIfTot,VecInd))
!                                HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                                HDiag=(REAL(HDiagTemp%v,r2))
!                            ENDIF
!                        ELSE
!                            CALL DecodeBitDet(nJ,CurrentDets(0:NIfTot,VecInd))
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
!                    CurrentDets(0:NIfTot,VecInd)=NewDets(0:NIfTot,i)
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
!                CurrentDets(0:NIfTot,VecInd)=SpawnedParts(0:NIfTot,IndSpawned)
!                CurrentSign(VecInd)=SpawnedSign(IndSpawned)
!                IndSpawned=IndSpawned+1
!
!                do while(DetBitEQ(SpawnedParts(:,IndSpawned-1),SpawnedParts(:,IndSpawned)).and.(IndSpawned.le.ValidSpawned))
!                    CurrentSign(VecInd)=CurrentSign(VecInd)+SpawnedSign(IndSpawned)
!                    IndSpawned=IndSpawned+1
!                enddo
!                    
!                TotParts=TotParts+abs(CurrentSign(VecInd))
!
!                IF(.not.tRegenDiagHEls) THEN
!!Need to find H-element!
!                    IF(DetBitEQ(CurrentDets(0:NIfTot,VecInd),iLutHF,NIfTot)) THEN
!!We know we are at HF - HDiag=0
!                        HDiag=0.D0
!                        IF(tHub.and.tReal) THEN
!!Reference determinant is not HF
!                            CALL DecodeBitDet(nJ,CurrentDets(0:NIfTot,VecInd))
!                            HDiagTemp=GetHElement2(nJ,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!                            HDiag=(REAL(HDiagTemp%v,r2))
!                        ENDIF
!                    ELSE
!                        CALL DecodeBitDet(nJ,CurrentDets(0:NIfTot,VecInd))
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
        use DetBitOps, only: DecodeBitDet
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
            CALL DecodeBitDet(nJ,SpawnedParts(0:NIfTot,i))
            HashArray1(i)=CreateHash(nJ)
!            IF(Iter.eq.1346.and.(HashArray1(i).eq.2905380077198165348)) THEN
!                WRITE(6,*) "Hash found, ",i,SpawnedSign(i),HashArray1(i),SpawnedParts(0:NIfTot,i)
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

        CALL MPIAlltoAllI(sendcounts,1,recvcounts,1,error)

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
        CALL MPIAlltoAllvI8(HashArray1(1:ValidSpawned),sendcounts,disps,HashArray2(1:MaxIndex),recvcounts,recvdisps,error)

!The signs of the hashes, index and CPU also need to be taken with them.
        CALL MPIAlltoAllvI(SpawnedSign(1:ValidSpawned),sendcounts,disps,SpawnedSign2(1:MaxIndex),recvcounts,recvdisps,error)
        CALL MPIAlltoAllvI(IndexTable1(1:ValidSpawned),sendcounts,disps,IndexTable2,recvcounts,recvdisps,error)
        CALL MPIAlltoAllvI(ProcessVec1(1:ValidSpawned),sendcounts,disps,ProcessVec2,recvcounts,recvdisps,error)

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

        CALL MPIAlltoAllI(sendcounts,1,recvcounts,1,error)

!We can now get recvdisps from recvcounts in the same way we obtained disps from sendcounts
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo

        ToAnnihilateonProc=recvdisps(nProcessors)+recvcounts(nProcessors)

        CALL MPIAlltoAllvI(IndexTable1(1:ToAnnihilateonProc),sendcounts,disps,IndexTable2,recvcounts,recvdisps,error)
        CALL MPIAlltoAllvI(SpawnedSign(1:ToAnnihilateonProc),sendcounts,disps,SpawnedSign2(1:),recvcounts,recvdisps,error)

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
        INTEGER :: iLut(0:NIfTot),MinInd,MaxInd,PartInd,DetArray(0:NIfTot,1:MaxInd)
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

        N=MinInd
        do while(N.le.MaxInd)
            Comp=DetBitLT(DetArray(:,N),iLut(:),NIfDBO)
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
        INTEGER :: iLut(0:NIfTot),MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(CurrentDets(:,MaxInd),iLut(:),NIfDBO)
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
            Comp=DetBitLT(CurrentDets(:,N),iLut(:),NIfDBO)

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
                    Comp=DetBitLT(CurrentDets(:,i+1),iLut(:),NIfDBO)
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


!In this routine, we want to search through the list of spawned particles. For each spawned particle, we binary search the list of particles on the processor
!to see if an annihilation event can occur. The annihilated particles are then removed from the spawned list
!to the whole list of spawned particles at the end of the routine.
!In the main list, we change the 'sign' element of the array to zero. These will be deleted at the end of the total annihilation step.
    SUBROUTINE AnnihilateSpawnedParts(ValidSpawned,TotWalkersNew)
        INTEGER :: ValidSpawned,MinInd,TotWalkersNew,PartInd,i,j,k,ToRemove,VecInd,SignProd,DetsMerged,PartIndex!,SearchInd,AnnihilateInd
        INTEGER :: ExcitLevel
        LOGICAL :: tSuccess,tSuc!,tSkipSearch

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
!                WRITE(6,'(3I20,A,3I20)') SpawnedParts(:,i),' equals ',CurrentDets(:,PartInd)
                
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
                            ExcitLevel = FindBitExcitLevel(SpawnedParts(:,i),&
                                                           iLutHF, nel)
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
                                WRITE(6,*) "***",SpawnedParts(0:NIftot,i)
                                Call WriteBitDet(6,SpawnedParts(0:NIfTot,i),.true.)
                                CALL Stop_All("AnnihilateSpawnedParts","Cannot find corresponding FCI determinant when histogramming")
                            ENDIF
                        ENDIF

                        CurrentSign(PartInd)=0
                        IF(SpawnedSign(i).eq.0) THEN
!The number of particles were equal and opposite. We want to remove this entry from the spawned list.
                            ToRemove=ToRemove+1
                        ELSEIF(tTruncInitiator) THEN
!If we are doing a CAS star calculation - then if the walkers that are left after annihilation have been spawned from determinants outside the active space,
!then it is like these have been spawned on an unoccupied determinant and they are killed.
                            IF(SpawnedParts(NIfTot,i).eq.1) THEN
                                NoAborted=NoAborted+ABS(SpawnedSign(i))
!                                WRITE(6,'(I20,A,3I20)') SpawnedSign(i),'walkers aborted from determinant:',SpawnedParts(:,i)
                                SpawnedSign(i)=0
                                ToRemove=ToRemove+1
                            ENDIF
                        ENDIF

                    ELSE
!There are more particles in the main list, than the spawned list. We want to annihilate all particles from the spawned list, but only some from main list.
                        CurrentSign(PartInd)=CurrentSign(PartInd)+SpawnedSign(i)
                        Annihilated=Annihilated+2*(abs(SpawnedSign(i)))
                        
                        IF(tHistSpawn) THEN
!We want to histogram where the particle annihilations are taking place.
                            ExcitLevel = FindBitExcitLevel(SpawnedParts(:,i),&
                                                           iLutHF, nel)
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
                                WRITE(6,*) "***",SpawnedParts(0:NIfTot,i)
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
                    ELSEIF(tTruncInitiator) THEN
!If doing a CAS star calculation - then if the signs on the current list is 0, and the walkers in the spawned list came from outside the cas space, these need to be killed.                        
                        IF(SpawnedParts(NIfTot,i).eq.1) THEN
                            NoAborted=NoAborted+ABS(SpawnedSign(i))
!                            WRITE(6,'(I20,A,3I20)') SpawnedSign(i),'walkers aborted from determinant:',SpawnedParts(:,i)
                            SpawnedSign(i)=0
                            ToRemove=ToRemove+1
                        ENDIF
                    ENDIF
                ENDIF


!                do while((DetBitEQ(SpawnedParts(:,i),CurrentDets(:,SearchInd))).and.(SearchInd.ge.1))
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
!                    do while((DetBitEQ(SpawnedParts(0:NIfTot,i),CurrentDets(0:NIfTot,SearchInd))).and.(SearchInd.le.TotWalkersNew))
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

            ELSEIF(tTruncInitiator) THEN
!Determinant in newly spawned list is not found in currentdets - usually this would mean the walkers just stay in this list and get merged later - but in this case we            
!want to check where the walkers came from - because if the newly spawned walkers are from a parent outside the active space they should be killed - as they have been
!spawned on an unoccupied determinant.
                IF(SpawnedParts(NIfTot,i).eq.1) THEN    !Walkers came from outside cas space.
                    NoAborted=NoAborted+ABS(SpawnedSign(i))
!                    WRITE(6,'(I20,A,3I20)') SpawnedSign(i),'walkers aborted from determinant:',SpawnedParts(:,i)
                    SpawnedSign(i)=0
                    ToRemove=ToRemove+1
                ENDIF
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
                    SpawnedParts2(0:NIfTot,i-DetsMerged)=SpawnedParts(0:NIfTot,i)
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
#ifdef PARALLEL
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
        CALL MPI_BSend(SpawnedParts(0:NIfTot,1:ValidSpawned),ValidSpawned*(NIfTot+1),MPI_INTEGER,MOD(iProcIndex+1,nProcessors),456,MPI_COMM_WORLD,error)
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

        CALL MPI_Recv(SpawnedParts2(0:NIfTot,1:ValidSpawned),ValidSpawned*(NIfTot+1),MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),456,MPI_COMM_WORLD,Stat,error)
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

#endif

    END SUBROUTINE RotateParticles

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
!        CALL MPIBarrier(error)
        CALL set_timer(Comms_Time,30)

        CALL MPIAlltoAllI(sendcounts,1,recvcounts,1,error)

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
        CALL MPIAlltoAllvI8(Hash2Array(1:TotWalkersNew),sendcounts,disps,HashArray(1:MaxIndex),recvcounts,recvdisps,error)        

!The signs of the hashes, index and CPU also need to be taken with them.
        CALL MPIAlltoAllvI(NewSign(1:TotWalkersNew),sendcounts,disps,CurrentSign,recvcounts,recvdisps,error)
        CALL MPIAlltoAllvI(IndexTable(1:TotWalkersNew),sendcounts,disps,Index2Table,recvcounts,recvdisps,error)
        CALL MPIAlltoAllvI(ProcessVec(1:TotWalkersNew),sendcounts,disps,Process2Vec,recvcounts,recvdisps,error)
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
!        CALL MPIBarrier(error)
        CALL set_timer(Comms_Time,30)

        CALL MPIAlltoAllI(sendcounts,1,recvcounts,1,error)

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
        CALL MPIAlltoAllvI(IndexTable(1:ToAnnihilateonProc),sendcounts,disps,Index2Table,recvcounts,recvdisps,error)
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

!********  The rest of the routines in this module, are for the annihilation of "minor" and "dominant" walkers as well as the guiding function  ********   
    
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
        CALL Sort2BitDetsPlus3(MinorValidSpawned,MinorSpawnDets(0:NIfTot,1:MinorValidSpawned),MinorSpawnParent(0:NIfTot,1:MinorValidSpawned),&
        &MinorSpawnSign(1:MinorValidSpawned))
        
!        IF(Iter.gt.1220) THEN
!            WRITE(6,*) 'sort'
!            CALL FLUSH(6)
!        ENDIF

! Make sure all processors have done this before carrying on.
        CALL MPIBarrier(error)


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
            ALLOCATE(mpibuffer(8*(MaxSpawned+1)*(NIfTot+3)),stat=ierr)
            IF(ierr.ne.0) THEN
                CALL Stop_All("RotoAnnihilateMinor","Error allocating memory for transfer buffers...")
            ENDIF
#ifdef PARALLEL
            CALL MPI_Buffer_attach(mpibuffer,8*(MaxSpawned+1)*(NIfTot+3),error)
#endif
            IF(error.ne.0) THEN
                CALL Stop_All("RotoAnnihilateMinor","Error allocating memory for transfer buffers...")
            ENDIF
        ENDIF
        CALL MPIBarrier(error)

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

            CALL MPIBarrier(error)
            CALL RotateMinorParticles(MinorValidSpawned)
 
!            IF(Iter.gt.1220) THEN
!                WRITE(6,*) 'last rotate'
!                CALL FLUSH(6)
!            ENDIF
!Detach buffers
#ifdef PARALLEL
            CALL MPI_Buffer_detach(mpibuffer,8*(MaxSpawned+1)*(NIfTot+3),error)
#endif
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
        CALL MPIBarrier(error)

        CALL MPIAlltoAllI(sendcounts(1:nProcessors),1,recvcounts(1:nProcessors),1,error)

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
        CALL MPIBarrier(error)

!Now send the chunks of hashes to the corresponding processors
!All the '2' arrays are like the 'All' arrays.
!TempMinorSpawnSign is the Signs from each processor, when just MinorSpawnSign is the 'All' array.
        CALL MPIAlltoAllvI8(HashArray(1:MinorValidSpawned),sendcounts,disps,Hash2Array(1:MaxIndex),recvcounts,recvdisps,error)        

!        tWrite=.false.
!        IF(MinorValidSpawned.gt.3) THEN
!            WRITE(6,*) 'TempMinorSpawnSign'
!            do i=1,MinorValidSpawned
!                WRITE(6,*) TempMinorSpawnSign(i)
!            enddo
!            tWrite=.true.
!        ENDIF

!The signs of the hashes, index and CPU also need to be taken with them.
        CALL MPIAlltoAllvI(TempMinorSpawnSign(1:MinorValidSpawned),sendcounts,disps,MinorSpawnSign(1:MaxIndex),recvcounts,recvdisps,error)
        CALL MPIAlltoAllvI(IndexTable(1:MinorValidSpawned),sendcounts,disps,Index2Table(1:MaxIndex),recvcounts,recvdisps,error)
        CALL MPIAlltoAllvI(ProcessVec(1:MinorValidSpawned),sendcounts,disps,Process2Vec(1:MaxIndex),recvcounts,recvdisps,error)
        
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
        CALL MPIBarrier(error)

        CALL MPIAlltoAllI(sendcounts,1,recvcounts,1,error)

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

        CALL MPIBarrier(error)

!Perform another matrix transpose of the annihilation data using MPI_AlltoAllv, to send the data back to its correct Processor
!The signs of the hashes, index and CPU also need to be taken with them. (CPU does not need to be taken - every element of CPU should be equal to the rank of the processor+1)
!Hash also does not need to be taken, but will be taken as a precaution
        CALL MPIAlltoAllvI(IndexTable(1:ToAnnihilateonProc),sendcounts,disps,Index2Table,recvcounts,recvdisps,error)


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
        INTEGER :: i,MinorValidSpawned,error
#ifdef PARALLEL
        INTEGER :: Stat(MPI_STATUS_SIZE)

! This is the number of particles spawned (and still alive).  Must be sent with the arrays so the next processor knows the size.        
        MinorSpawnSign(0)=MinorValidSpawned

!Send the signs of the particles (number sent is in the first element)
        CALL MPI_BSend(MinorSpawnSign(0:MinorValidSpawned),MinorValidSpawned+1,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),123,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in sending signs")
        ENDIF

!...then send the particles themselves...
        CALL MPI_BSend(MinorSpawnDets(0:NIfTot,1:MinorValidSpawned),MinorValidSpawned*(NIfTot+1),MPI_INTEGER,MOD(iProcIndex+1,nProcessors),456,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in sending particles")
        ENDIF


!...and then send the parents of the walkers...
        CALL MPI_BSend(MinorSpawnParent(0:NIfTot,1:MinorValidSpawned),MinorValidSpawned*(NIfTot+1),MPI_INTEGER,MOD(iProcIndex+1,nProcessors),789,MPI_COMM_WORLD,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in sending particle parents")
        ENDIF

        CALL MPIBarrier(error)

!Receive signs (let it receive the maximum possible (only the first ValidSpawned will be updated.))
        CALL MPI_Recv(MinorSpawnSign2(0:MaxSpawned),MaxSpawned+1,MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),123,MPI_COMM_WORLD,Stat,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in receiving signs")
        ENDIF

!Update the ValidSpawned variable for this new set of data we are about to receive...
        MinorValidSpawned=MinorSpawnSign2(0)

        CALL MPI_Recv(MinorSpawnDets2(0:NIfTot,1:MinorValidSpawned),MinorValidSpawned*(NIfTot+1),MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),456,MPI_COMM_WORLD,Stat,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in receiving particles")
        ENDIF

        CALL MPI_Recv(MinorSpawnParent2(0:NIfTot,1:MinorValidSpawned),MinorValidSpawned*(NIfTot+1),MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),789,MPI_COMM_WORLD,Stat,error)
        IF(error.ne.MPI_SUCCESS) THEN
            CALL Stop_All("RotateParticles","Error in receiving particle parents")
        ENDIF

        do i=1,MinorValidSpawned
            MinorSpawnDets(0:NIfTot,i)=MinorSpawnDets2(0:NIfTot,i)
            MinorSpawnParent(0:NIfTot,i)=MinorSpawnParent2(0:NIfTot,i)
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

#endif

    END SUBROUTINE RotateMinorParticles



    SUBROUTINE AnnihilateMinorSpawnedParts(MinorValidSpawned,NoMinorWalkersNew)
        INTEGER :: i,j,k,MinorValidSpawned,NoMinorWalkersNew,MinInd,MaxDetInd,MinDetInd,SumDetPop,SumMinorDetPop
        INTEGER :: ToRemove,DetsMerged,SignProd,PartInd,EqDetPopsTag,ierr,FinalMinorDet
        LOGICAL :: DetsEq,tSuccess,tAnnihilateOne
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
            IF((i+j).le.MinorValidSpawned) DetsEq=DetBitEQ(MinorSpawnDets(0:NIfTot,i),MinorSpawnDets(0:NIfTot,i+j),NIfDBO)
            do while (DetsEq)
                SumMinorDetPop=SumMinorDetPop+MinorSpawnSign(i+j)
                j=j+1
                IF((i+j).gt.MinorValidSpawned) EXIT
                DetsEq=DetBitEQ(MinorSpawnDets(0:NIfTot,i),MinorSpawnDets(0:NIfTot,i+j),NIfDBO)
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
                DetsEq=DetBitEQ(MinorSpawnDets(0:NIfTot,i),MinorStarDets(0:NIfTot,PartInd-j),NIfDBO)
                do while (DetsEq)
                    ! If the determinant is still equal, add the walkers on it to SumDetPop, and this index becomes the minimum.
                    SumDetPop=SumDetPop+MinorStarSign(PartInd-j)
                    MinDetInd=PartInd-j
                    j=j+1
                    DetsEq=DetBitEQ(MinorSpawnDets(0:NIfTot,i),MinorStarDets(0:NIfTot,PartInd-j),NIfDBO)
                    ! If this is true, the walkers on the next determinant will be added.
                enddo

                ! Now check those above the determinant found.
                j=1
                DetsEq=DetBitEQ(MinorSpawnDets(0:NIfTot,i),MinorStarDets(0:NIfTot,PartInd+j),NIfDBO)
                do while (DetsEq)
                    ! If the determinant is still equal, add the walkers on it to SumDetPop, and this index becomes the minimum.
                    SumDetPop=SumDetPop+MinorStarSign(PartInd+j)
                    MaxDetInd=PartInd+j
                    j=j+1
                    DetsEq=DetBitEQ(MinorSpawnDets(0:NIfTot,i),MinorStarDets(0:NIfTot,PartInd+j),NIfDBO)
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
                        DetsEq=DetBitEQ(MinorSpawnParent(0:NIfTot,i),MinorStarParent(0:NIfTot,j),NIfDBO)
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
!                        WRITE(6,*) MinorStarDets(0:NIfTot,i),'*',MinorStarParent(0:NIfTot,i)
!                        do j=MinDetInd,MaxDetInd
!                            WRITE(6,*) MinorSpawnDets(0:NIfTot,j),'*',MinorSpawnParent(0:NIfTot,j)
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
                    MinorSpawnDets2(0:NIfTot,i-DetsMerged)=MinorSpawnDets(0:NIfTot,i)
                    MinorSpawnSign2(i-DetsMerged)=MinorSpawnSign(i)
                    MinorSpawnParent2(0:NIfTot,i-DetsMerged)=MinorSpawnParent(0:NIfTot,i)
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
                MinorSpawnDets(0:NIfTot,i)=MinorSpawnDets2(0:NIfTot,i)
                MinorSpawnParent(0:NIfTot,i)=MinorSpawnParent2(0:NIfTot,i)
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
                    MinorStarDets(0:NIfTot,i-DetsMerged)=MinorStarDets(0:NIfTot,i)
                    MinorStarSign(i-DetsMerged)=MinorStarSign(i)
                    MinorStarParent(0:NIfTot,i-DetsMerged)=MinorStarParent(0:NIfTot,i)
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

            CALL MergeListswH2(NoMinorWalkersNew,MaxWalkersPart,MinorValidSpawned,MinorSpawnDets(0:NIfTot,1:MinorValidSpawned),&
            &MinorSpawnParent(0:NIfTot,1:MinorValidSpawned),MinorSpawnSign(1:MinorValidSpawned))
            
        ENDIF

        NoMinorWalkers=NoMinorWalkersNew    


    END SUBROUTINE InsertRemoveMinorParts

    
!Do a binary search the guiding function dets, between the indices of MinInd and MaxInd. If successful, tSuccess will be true and 
!PartInd will be a coincident determinant. If there are multiple values, the chosen one may be any of them...
!If failure, then the index will be one less than the index that the particle would be in if it was present in the list.
!(or close enough!)
    SUBROUTINE BinSearchGuideParts(iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: iLut(0:NIfTot),MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
        i=MinInd
        j=MaxInd
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(GuideFuncDets(:,N),iLut(:),NIfDBO)

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
                    Comp=DetBitLT(GuideFuncDets(:,i+1),iLut(:),NIfDBO)
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
        INTEGER :: iLut(0:NIfTot),MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
        i=MinInd
        j=MaxInd
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(MinorStarDets(:,N),iLut(:),NIfDBO)

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
                    Comp=DetBitLT(MinorStarDets(:,i+1),iLut(:),NIfDBO)
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
        INTEGER :: iLut(0:NIfTot),MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp,AllExcDets(0:NIfTot,MinInd:MaxInd)
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
        i=MinInd
        j=MaxInd
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(AllExcDets(:,N),iLut(:),NIfDBO)

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
                    Comp=DetBitLT(AllExcDets(:,i+1),iLut(:),NIfDBO)
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
    
    SUBROUTINE RotoAnnihilGuidingFunc(ValidSpawned)
! This routine takes the spawned particles (that have already been annihilated with themselves and with the main wavefunction), and tries to 
! annihilate them with the guiding function.
! The guiding function itself may be annihilated also, but does not spawn or die by itself.
! However, if the guiding function is completely annihilated on one processor, the spawned particles must be rotated around the other processors
! to check for possible annihilations there.
        use SystemData , only : G1,nBasis,Brr,NMsh,NMax,Alat,ECore,nBasis,nBasisMax
        use IntegralsData , only : UMat,fck
        use Determinants , only : GetHElement2
        use DetBitOps, only: DecodeBitDet
        INTEGER :: i,j,n,ValidSpawned,InitNoDetstoRotate,NoDetstoRotate,CombSign,error
        INTEGER :: ExcitLevel,DoubDet(NEl)
        TYPE(HElement) :: HDoubTemp
        REAL*8 :: Hdoub
        LOGICAL :: tRotateSpawnedTemp,tRotateSpawned,tDetinSpawnList,DetsEq
#ifdef PARALLEL
        INTEGER :: Stat(MPI_STATUS_SIZE)
#endif


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
                    DetsEq=DetBitEQ(SpawnedParts(0:NIfTot,i),GuideFuncDets(0:NIfTot,j),NIfDBO)
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
                                DetstoRotate(0:NIfTot,NoDetstoRotate)=SpawnedParts(0:NIfTot,i)
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
                            DetstoRotate(0:NIfTot,NoDetstoRotate)=SpawnedParts(0:NIfTot,i)
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
        CALL MPIAllReduceLORScal(tRotateSpawnedTemp,tRotateSpawned,error)


!The allocated DetstoRotate arrays are as big as iGuideDets (the number of determinants in the guiding function), but will only need to rotate 
!a portion of these (those determinants for which the guiding function has had all its walkers annihilated).

        IF(tRotateSpawned) THEN !If NoDetstoRotate is greater than 0 on any processor, need to rotate all arrays, otherwise will overwrite each other etc.

            InitNoDetstoRotate=NoDetstoRotate
            !For now, rotate this same sized array each time, even if not completely necessary.
            !We are currently just making the determinant (and its sign) 0 if we no longer want to rotate it, but could in the future remove it from the
            !array.  Probably doesn't make all that much difference because will never find a 0 determinant in the guiding function.

            do n=1,nProcessors-1
            !Rotate the DetstoRotate, SigntoRotate and NoDetstoRotate values.

#ifdef PARALLEL

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
                CALL MPI_BSend(DetstoRotate(0:NIfTot,1:InitNoDetstoRotate),(NIfTot+1)*InitNoDetstoRotate,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),456,MPI_COMM_WORLD,error)
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
                CALL MPI_Recv(DetstoRotate2(0:NIfTot,1:InitNoDetstoRotate),InitNoDetstoRotate*(NIfTot+1),MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),456,MPI_COMM_WORLD,Stat,error)
                IF(error.ne.MPI_SUCCESS) THEN
                    CALL Stop_All("RotoAnnihilGuidingFunc","Error in receiving particles")
                ENDIF

                do i=1,InitNoDetstoRotate
                    SigntoRotate(i)=SigntoRotate2(i)
                    DetstoRotate(0:NIfTot,i)=DetstoRotate2(0:NIfTot,i)
                enddo
#endif

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
                            DetsEq=DetBitEQ(DetstoRotate(0:NIfTot,i),GuideFuncDets(0:NIfTot,j),NIfDBO)

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

#ifdef PARALLEL

            !Send the sign of those we want to rotate to the next processor.
            CALL MPI_BSend(SigntoRotate(0:InitNoDetstoRotate),InitNoDetstoRotate+1,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),123,MPI_COMM_WORLD,error)
            IF(error.ne.MPI_SUCCESS) THEN
                CALL Stop_All("RotoAnnihilGuidingFunc","Error in sending signs")
            ENDIF

            !Then send the determinants
            CALL MPI_BSend(DetstoRotate(0:NIfTot,1:InitNoDetstoRotate),(NIfTot+1)*InitNoDetstoRotate,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),456,MPI_COMM_WORLD,error)
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
            CALL MPI_Recv(DetstoRotate2(0:NIfTot,1:InitNoDetstoRotate),InitNoDetstoRotate*(NIfTot+1),MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),456,MPI_COMM_WORLD,Stat,error)
            IF(error.ne.MPI_SUCCESS) THEN
                CALL Stop_All("RotoAnnihilGuidingFunc","Error in receiving particles")
            ENDIF

            do i=1,InitNoDetstoRotate
                SigntoRotate(i)=SigntoRotate2(i)
                DetstoRotate(0:NIfTot,i)=DetstoRotate2(0:NIfTot,i)
            enddo

#endif

            !Now add all the remaining DetstoRotate2 and their signs to the SpawnedPart and SpawnedSign lists.
            !Since I just copied the determinants from SpawnedPart to DetstoRotate, any in DetstoRotate will already been in SpawnedPart.
            !Need to just search for it and overwrite its spin value.
            do j=1,InitNoDetstoRotate
                tDetinSpawnList=.false.
                do i=1,ValidSpawned
                    DetsEq=.false.
                    DetsEq=DetBitEQ(SpawnedParts(0:NIfTot,i),DetstoRotate(0:NIfTot,j),NIfDBO)
                    IF(DetsEq) THEN
                        SpawnedSign(i)=SigntoRotate(j)
                        tDetinSpawnList=.true.
                        EXIT
                    ENDIF
                enddo
                IF(.not.tDetinSpawnList) THEN
                    WRITE(6,*) 'Determinant from rotate list : ',DetstoRotate(0:NIfTot,j)
!                    do i=1,ValidSpawned
!                        WRITE(6,*) SpawnedParts(0:NIfTot,i)
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
            ExcitLevel = FindBitExcitLevel(GuideFuncDets(:,i), iLutHF, 2)
            IF(ExcitLevel.eq.2) THEN
                DoubDet(:)=0
                CALL DecodeBitDet(DoubDet,GuideFuncDets(0:NIfTot,i))
                HdoubTemp=GetHElement2(HFDet,DoubDet,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,ExcitLevel,ECore)
                HDoub=REAL(HDoubTemp%v,r2)
                GuideFuncDoub=GuideFuncDoub+(GuideFuncSign(i)*Hdoub)
            ENDIF
        enddo


    ENDSUBROUTINE RotoAnnihilGuidingFunc



END MODULE AnnihilationMod
