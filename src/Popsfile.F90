MODULE PopsfileMod

    use SystemData, only: nel, tHPHF, tFixLz, tCSF, nBasis, tNoBrillouin
    use CalcData, only: tTruncInitiator,DiagSft,tWalkContGrow,nEquilSteps,ScaleWalkers, &
                        tReadPopsRestart, tRegenDiagHEls,InitWalkers, tReadPopsChangeRef, &
                        nShiftEquilSteps,iWeightPopRead
    use DetBitOps, only: DetBitLT,FindBitExcitLevel,DetBitEQ
    use Determinants, only : get_helement,write_det
    use hphf_integrals, only: hphf_diag_helement
    USE dSFMT_interface , only : genrand_real2_dSFMT
    use FciMCData
    use bit_reps
    use Parallel
    use AnnihilationMod, only: DetermineDetProc
    USE Logging , only : iWritePopsEvery,tPopsFile,iPopsPartEvery,tBinPops
    USE Logging , only : tPrintPopsDefault
    use sort_mod

    implicit none

    contains

    !   V.3 POPSFILE ROUTINES   !
!This routine reads in particle configurations from a POPSFILE v.3.
!EndPopsList is the number of entries in the POPSFILE to read, and ReadBatch is the number of determinants
!which can be read in in a single batch.
    SUBROUTINE ReadFromPopsfilev3(EndPopsList,ReadBatch,CurrWalkers64,CurrParts,CurrHF)
        integer(8) , intent(in) :: EndPopsList  !Number of entries in the POPSFILE.
        integer , intent(in) :: ReadBatch       !Size of the batch of determinants to read in in one go.
        integer(int64) , intent(out) :: CurrWalkers64    !Number of determinants which end up on a given processor.
        integer(int64) , dimension(lenof_sign) , intent(out) :: CurrParts
        integer , dimension(lenof_sign) , intent(out) :: CurrHF
        integer :: CurrWalkers
        integer :: iunit,i,j,BatchReadTag,ierr,PopsInitialSlots(0:nProcessors-1)
        real(8) :: BatchSize
        integer :: PopsSendList(0:nProcessors-1),proc,sendcounts(nProcessors),disps(nProcessors)
        integer :: MaxSendIndex,recvcount,err
        integer(n_int) , allocatable :: BatchRead(:,:)
        integer(n_int) :: WalkerTemp(0:NIfTot)
        integer(8) :: Det,AllCurrWalkers,TempCurrWalkers
        logical :: FormPops,BinPops,tReadAllPops,tStoreDet
        integer , dimension(lenof_sign) :: SignTemp
        integer :: TempNI(NEl) 
        character(len=*) , parameter :: this_routine='ReadFromPopsfilev3'
        HElement_t :: HElemTemp
        !variables from header file
        logical :: tPop64Bit,tPopHPHF,tPopLz
        integer :: iPopLenof_sign,iPopNEl,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot
        integer(8) :: iPopAllTotWalkers
        real(8) :: PopDiagSft
        integer(8) , dimension(lenof_sign) :: PopSumNoatHF
        HElement_t :: PopAllSumENum

        call open_pops_head(iunit,formpops,binpops)
        IF(FormPops) THEN
            call ReadPopsHeadv3(iunit,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot)

                if(EndPopsList.ne.iPopAllTotWalkers) then
                    call stop_all(this_routine,"Error in assessing number of entries in POPSFILE")
                endif

        ELSEIF(BinPops) THEN
            if(iProcIndex.eq.root) then
                close(iunit)    !iunit here refers to the header file.
                OPEN(iunit,FILE='POPSFILEBIN',Status='old',form='unformatted')
            endif
        ENDIF
        
        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,I4,A)") "Although ",EndPopsList," configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            ENDIF
            if(ScaleWalkers.ne.1) call warning(this_routine,"ScaleWalkers parameter found, but not implemented in POPSFILE v3 - ignoring.")

        ENDIF

        BatchSize=REAL(ReadBatch,dp)/REAL(nProcessors,dp)
        if(iProcIndex.eq.Root) then
            !Create PopsInitialSlots
            do i=0,nProcessors-1
                PopsInitialSlots(i)=NINT(BatchSize*i)+1
            enddo
            !Allocate array to store particle to distribute
            allocate(BatchRead(0:NIfTot,1:ReadBatch),stat=ierr)
            CALL LogMemAlloc('BatchRead',ReadBatch*(NIfTot+1),size_n_int,this_routine,BatchReadTag,ierr)
        endif

        CurrHF=0        !Number of HF walkers on each node.
        CurrParts=0     !Number of walkers on each node.
        CurrWalkers=0   !Number of determinants on each node.
        Det=1
        tReadAllPops=.false.
        do while(.not.tReadAllPops)

            if(iProcIndex.eq.Root) then

                !Get ready for reading in the next batch of walkers
                BatchRead(:,:)=0
                PopsSendList(:)=PopsInitialSlots(:)

                do while(Det.le.EndPopsList)

                    tStoreDet=.false.
                    do while(.not.tStoreDet)
                        if(BinPops) then
                            read(iunit) WalkerTemp(:)
                        else
                            read(iunit,*) WalkerTemp(:)
                        endif
                        Det=Det+1

                        if(iWeightPopRead.ne.0) then
                            call extract_sign(WalkerTemp(:),SignTemp)
                            do i=1,lenof_sign
                                if(SignTemp(i).ge.iWeightPopRead) then
                                    tStoreDet=.true.
                                    exit
                                endif
                            enddo
                        else
                            tStoreDet=.true.
                        endif
                    enddo

                    proc = DetermineDetProc(WalkerTemp)
                    BatchRead(:,PopsSendList(proc)) = WalkerTemp(:)
                    PopsSendList(proc) = PopsSendList(proc) + 1
                    if(proc.ne.(nProcessors-1)) then
                        if(PopsInitialSlots(proc+1)-PopsSendList(proc).lt.2) then
                            exit  !Now distribute the particles
                        endif
                    else
                        if(ReadBatch-PopsSendList(proc).lt.2) then
                            exit  !Now distribute the particles
                        endif
                    endif

                enddo

                if(Det.gt.EndPopsList) tReadAllPops=.true.

                do j=0,nProcessors-1
                    sendcounts(j+1)=(PopsSendList(j)-(NINT(BatchSize*j)+1))*(NIfTot+1)
                    disps(j+1)=(NINT(BatchSize*j))*(NIfTot+1)
                enddo
                MaxSendIndex=disps(nProcessors)+sendcounts(nProcessors)

            endif

            !Now scatter the particles read in to their correct processors.
            call MPIScatter(sendcounts,recvcount,err)
            if(err.ne.0) call stop_all(this_routine,"MPI scatter error")
            call MPIScatterV(BatchRead(:,1:MaxSendIndex),sendcounts,disps,CurrentDets(:,CurrWalkers+1:MaxWalkersPart),recvcount,err)
            if(err.ne.0) call stop_all(this_routine,"MPI error")
            CurrWalkers=CurrWalkers+recvcount/(NIfTot+1)
            call MPIBCast(tReadAllPops)

        enddo

        close(iunit)

        !Test we have still got all determinants
        TempCurrWalkers=int(CurrWalkers,8)
        call MPISum(TempCurrWalkers,1,AllCurrWalkers)
        if(iProcIndex.eq.Root) then
            if((iWeightPopRead.eq.0).and.(AllCurrWalkers.ne.EndPopsList)) then
                call Stop_All(this_routine,"Not all walkers accounted for when reading in")
            endif
        endif

        if(iProcIndex.eq.Root) then
            deallocate(BatchRead)
            CALL LogMemDealloc(this_routine,BatchReadTag)
        endif

        !Order the determinants on all the lists.
        call sort (currentdets(:,1:CurrWalkers))

        !Run through all determinants on each node, and calculate the total number of walkers, and noathf
        do i=1,CurrWalkers
            call extract_sign(CurrentDets(:,i),SignTemp)
            CurrParts=CurrParts+abs(SignTemp)
            if(DetBitEQ(CurrentDets(:,i),iLutRef,NIfDBO)) then
                if(CurrHF(1).ne.0) then
                    call stop_all(this_routine,"HF already found, but shouldn't have")
                endif
                CurrHF=CurrHF+SignTemp 
                IF(.not.tRegenDiagHEls) CurrentH(i)=0.D0
            else
                if(.not.tRegenDiagHEls) THEN
                !Calculate diagonal matrix element
                    call decode_bit_det (TempnI, currentDets(:,i))
                    if (tHPHF) then
                        HElemTemp = hphf_diag_helement (TempnI,CurrentDets(:,i))
                    else
                        HElemTemp = get_helement (TempnI, TempnI, 0)
                    endif
                    CurrentH(i)=REAL(HElemTemp,dp)-Hii
                endif
            endif
        enddo

        CurrWalkers64=int(CurrWalkers,int64)    !Since this variable is eventually going to be
                                                !Totwalkers, it wants to be a 64 bit int.
    
    end subroutine ReadFromPopsfilev3

    subroutine CheckPopsParams(tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                    iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                    PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot,WalkerListSize)
        logical , intent(in) :: tPop64Bit,tPopHPHF,tPopLz
        integer , intent(in) :: iPopLenof_sign,iPopNel,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot
        integer(8) , intent(in) :: iPopAllTotWalkers
        real(8) , intent(in) :: PopDiagSft
        integer(8) , dimension(lenof_sign) , intent(in) :: PopSumNoatHF
        HElement_t , intent(in) :: PopAllSumENum
        integer , intent(out) :: WalkerListSize
        character(len=*) , parameter :: this_routine='CheckPopsParams'

        !Ensure all NIF and symmetry options the same as when popsfile was written out.
#ifdef __INT64
        if(.not.tPop64Bit) call stop_all(this_routine,"Popsfile created with 32 bit walkers, but now using 64 bit.")
#else
        if(tPop64Bit) call stop_all(this_routine,"Popsfile created with 64 bit walkers, but now using 32 bit.")
#endif
        if(tPopHPHF.neqv.tHPHF) call stop_all(this_routine,"Popsfile HPHF and input HPHF not same")
        if(tPopLz.neqv.tFixLz) call stop_all(this_routine,"Popsfile Lz and input Lz not same")
        if(iPopLenof_sign.ne.lenof_sign) call stop_all(this_routine,"Popsfile lenof_sign and input lenof_sign not same")
        if(iPopNEl.ne.NEl) call stop_all(this_routine,"Popsfile NEl and input NEl not same")
        if(PopNIfD.ne.NIfD) call stop_all(this_routine,"Popsfile NIfD and calculated NIfD not same")
        if(PopNIfY.ne.NIfY) call stop_all(this_routine,"Popsfile NIfY and calculated NIfY not same")
        if(PopNIfSgn.ne.NIfSgn) call stop_all(this_routine,"Popsfile NIfSgn and calculated NIfSgn not same")
        if(PopNIfFlag.ne.NIfFlag) call stop_all(this_routine,"Popsfile NIfFlag and calculated NIfFlag not same")
        if(PopNIfTot.ne.NIfTot) call stop_all(this_routine,"Popsfile NIfTot and calculated NIfTot not same")


        IF(.not.tWalkContGrow) THEN
!If we want the walker number to be stable, take the shift from the POPSFILE, otherwise, keep the input value.
            DiagSft=PopDiagSft
        ENDIF

        if(PopDiagSft.eq.0.D0) then
            !If the popsfile has a shift of zero, continue letting the population grow
            tWalkContGrow=.true.
            DiagSft=PopDiagSft
        endif

        if(tWalkContGrow) then
            !If continuing to grow, ensure we can allocate enough memory for what we hope to get the walker population to,
            !rather than the average number of determinants in the popsfile.
            WalkerListSize=max(initwalkers,NINT(real(iPopAllTotWalkers,8)/real(nProcessors,8)))
        else
            WalkerListSize=NINT(real(iPopAllTotWalkers,8)/real(nProcessors,8))
        endif

        AllSumNoatHF=PopSumNoatHF
        AllSumENum=PopAllSumENum
        PreviousCycles=iPopIter
    
    end subroutine CheckPopsParams


!Routine for reading in from iunit the header information from a popsile v3 file.
    subroutine ReadPopsHeadv3(iunithead,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot)
        integer , intent(in) :: iunithead
        logical , intent(out) :: tPop64Bit,tPopHPHF,tPopLz
        integer , intent(out) :: iPopLenof_sign,iPopNel,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot
        integer(8) , intent(out) :: iPopAllTotWalkers
        real(8) , intent(out) :: PopDiagSft
        integer(8) , dimension(lenof_sign) , intent(out) :: PopSumNoatHF
        HElement_t , intent(out) :: PopAllSumENum
        character(len=24) :: junk,junk2,junk3,junk4,junk5
        character(255) :: FirstLine
        integer :: PopsVersion

        PopsVersion=FindPopsfileVersion(iunithead)
        if(PopsVersion.ne.3) call stop_all("ReadPopsfileHeadv3","Wrong popsfile version for this routine.")
            
        if(iProcIndex.eq.root) then
            read(iunithead,'(A,L,A,L,A,L,A,I5,A,I7)') junk,tPop64Bit,junk2,tPopHPHF,junk3,tPopLz,junk4,iPopLenof_sign,junk5,iPopNEl
            read(iunithead,*) iPopAllTotWalkers
            read(iunithead,*) PopDiagSft 
            read(iunithead,*) PopSumNoatHF 
            read(iunithead,*) PopAllSumENum 
            read(iunithead,*) iPopIter 
            read(iunithead,*) PopNIfD 
            read(iunithead,*) PopNIfY 
            read(iunithead,*) PopNIfSgn 
            read(iunithead,*) PopNIfFlag 
            read(iunithead,*) PopNIfTot
        endif
        !Broadcast the read in values from the header to all nodes.
        call MPIBCast(tPop64Bit)
        call MPIBCast(tPopHPHF)
        call MPIBCast(tPopLz)
        call MPIBCast(iPopLenof_sign)
        call MPIBCast(iPopNEl)
        call MPIBCast(iPopAllTotWalkers)
        call MPIBCast(PopDiagSft)
        call MPIBCast(PopSumNoatHF)
        call MPIBCast(PopAllSumENum)
        call MPIBCast(iPopIter)
        call MPIBCast(PopNIfD)
        call MPIBCast(PopNIfY)
        call MPIBCast(PopNIfSgn)
        call MPIBCast(PopNIfFlag)
        call MPIBCast(PopNIfTot)

    end subroutine ReadPopsHeadv3
    
    !NOTE: This should only be used for the v3 POPSFILEs, since we only open the POPSFILE on the head node.
    subroutine open_pops_head(iunithead,formpops,binpops)
        use util_mod, only: get_free_unit
        integer , intent(out) :: iunithead
        logical , intent(out) :: formpops,binpops

        if(iProcIndex.eq.root) then
            iunithead=get_free_unit()
            inquire(file='POPSFILE',exist=formpops)
            if(formpops) then
                open(iunithead,file='POPSFILE',status='old')
                binpops=.false.
            else
                inquire(file='POPSFILEBIN',exist=binpops)
                if(binpops) then
                    open(iunithead,file='POPSFILEHEAD',status='old')
                else 
                    call stop_all("open_pops_head","No POPSFILEs detected...")
                endif
            endif
            rewind(iunithead)
        endif
        call MPIBCast(binpops) 
        call MPIBCast(formpops) 

    end subroutine open_pops_head


!Return the version number of the popsfile
    integer function FindPopsfileVersion(iunithead)
        integer, intent(in) :: iunithead
        logical :: formpops,binpops
        character(255) :: FirstLine

        if(iProcIndex.eq.root) then
            rewind(iunithead)
            read(iunithead,'(a255)') FirstLine

            if(index(FirstLine,'VERSION').eq.0) then
                FindPopsfileVersion=1
            else
                rewind(iunithead)
                read(iunithead,*) FirstLine,FirstLine,FirstLine,FindPopsfileVersion
            endif
        endif
        call MPIBCast(FindPopsfileVersion)

    end function FindPopsfileVersion


!This routine is the same as WriteToPopsfilePar, but does not require two main arrays to hold the data.
!The root processors data will be stored in a temporary array while it recieves the data from the other processors.
!This routine will write out to a popsfile. It transfers all walkers to the head node sequentially, so does not want to be called too often
    SUBROUTINE WriteToPopsfileParOneArr(Dets,nDets)
        use util_mod, only: get_unique_filename, get_free_unit
        use CalcData, only: iPopsFileNoWrite
        use Logging, only: tIncrementPops
        use constants, only: size_n_int,MpiDetInt,n_int
        integer(int64),intent(in) :: nDets !The number of occupied entries in Dets
        integer(kind=n_int),intent(in) :: Dets(0:nIfTot,1:nDets)
        INTEGER :: error
        integer(int64) :: WalkersonNodes(0:nProcessors-1)
        INTEGER :: Tag,Total,i,j,k
        INTEGER(KIND=n_int), ALLOCATABLE :: Parts(:,:)
        INTEGER :: PartsTag=0
        INTEGER :: nMaxDets
        integer :: iunit
        CHARACTER(len=*) , PARAMETER :: this_routine='WriteToPopsfileParOneArr'
        character(255) :: popsfile
        INTEGER, DIMENSION(lenof_sign) :: TempSign


        IF(lenof_sign.ne.1) THEN
            WRITE(6,*) "Cannot write complex walkers out to POPSFILE yet..."
        ENDIF

        CALL MPIBarrier(error)  !sync

!First, make sure we have up-to-date information - again collect AllTotWalkers,AllSumNoatHF and AllSumENum...
!        CALL MPI_Reduce(TotWalkers,AllTotWalkers,1,MPI_INTEGER,MPI_Sum,root,MPI_COMM_WORLD,error)    
!Calculate the energy by summing all on HF and doubles - convert number at HF to a real since no int*8 MPI data type
        CALL MPISum(SumNoatHF,1,AllSumNoatHF)
        CALL MPISum(SumENum,1,AllSumENum)

!We also need to tell the root processor how many particles to expect from each node - these are gathered into WalkersonNodes
        CALL MPIAllGather(nDets,WalkersonNodes,error)

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
#ifdef __INT64
            IF(iPopsPartEvery.ne.1) THEN
                IF(tBinPops) THEN
                    WRITE(6,"(A,I12,A)") "Writing a 64-bit binary reduced POPSFILEBIN, printing a total of ",INT(AllTotWalkers,int64), " particles."
                ELSE
                    WRITE(6,"(A,I12,A)") "Writing a 64-bit reduced POPSFILE, printing a total of ",INT(AllTotWalkers,int64), " particles."
                ENDIF
            ELSE
                IF(tBinPops) THEN
                    WRITE(6,*) "Writing to 64-bit binary POPSFILEBIN..."
                ELSE
                    WRITE(6,*) "Writing to 64-bit POPSFILE..."
                ENDIF
            ENDIF
#else
            IF(iPopsPartEvery.ne.1) THEN
                IF(tBinPops) THEN
                    WRITE(6,"(A,I12,A)") "Writing a binary reduced POPSFILEBIN, printing a total of ",INT(AllTotWalkers,int64), " particles."
                ELSE
                    WRITE(6,"(A,I12,A)") "Writing a reduced POPSFILE, printing a total of ",INT(AllTotWalkers,int64), " particles."
                ENDIF
            ELSE
                IF(tBinPops) THEN
                    WRITE(6,*) "Writing to binary POPSFILEBIN..."
                ELSE
                    WRITE(6,*) "Writing to POPSFILE..."
                ENDIF
            ENDIF
#endif
            IF(tBinPops) THEN
                call get_unique_filename('POPSFILEHEAD',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
            ELSE
                call get_unique_filename('POPSFILE',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
            ENDIF
            iunit = get_free_unit()
            OPEN(iunit,FILE=popsfile,Status='replace')
            WRITE(iunit,"(A)") "# POPSFILE VERSION 2"
#ifdef __INT64
            WRITE(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') '64BitDets=',.TRUE.,'HPHF=',tHPHF,'Lz=',tFixLz,'Initiator=',tTruncInitiator
#else
            WRITE(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') '64BitDets=',.FALSE.,'HPHF=',tHPHF,'Lz=',tFixLz,'Initiator=',tTruncInitiator
#endif
            WRITE(iunit,*) AllTotWalkers,"   TOTWALKERS (all nodes)"
            WRITE(iunit,*) DiagSft,"   DIAG SHIFT"
            WRITE(iunit,*) AllSumNoatHF,"   SUMNOATHF (all nodes)"
            WRITE(iunit,*) AllSumENum,"   SUMENUM ( \sum<D0|H|Psi> - all nodes)"
            WRITE(iunit,*) Iter+PreviousCycles,"   PREVIOUS CYCLES"
            IF(tBinPops) THEN
                CLOSE(iunit)
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.true.,iPopsFileNoWrite,popsfile)
                OPEN(iunit,FILE=popsfile,Status='replace',form='unformatted')
            ENDIF

            IF(tBinPops) THEN
                do j=1,nDets
!First write out walkers on head node
                    IF(mod(j,iPopsPartEvery).eq.0) THEN
                        call extract_sign(Dets(:,j),TempSign)
                        WRITE(iunit) Dets(0:NIfDBO,j),TempSign(:)
                    ENDIF
                enddo
            ELSE
                do j=1,nDets
!First write out walkers on head node
                    IF(mod(j,iPopsPartEvery).eq.0) THEN
                        do k=0,NIfDBO
                            WRITE(iunit,"(I24)",advance='no') Dets(k,j)
                        enddo
                        call extract_sign(Dets(:,j),TempSign)
                        WRITE(iunit,*) TempSign(:)
                    ENDIF
                enddo
            ENDIF
!            WRITE(6,*) "Written out own walkers..."
!            CALL FLUSH(6)

!Now, we copy the head nodes data to a new array...
            nMaxDets=maxval(WalkersOnNodes)
            ALLOCATE(Parts(0:NIfTot,nMaxDets),stat=error)
            CALL LogMemAlloc('Parts',int(nMaxDets,int32)*(NIfTot+1),size_n_int,this_routine,PartsTag,error)

!Now we need to receive the data from each other processor sequentially
!We can overwrite the head nodes information, since we have now stored it elsewhere.
            do i=1,nProcessors-1
!Run through all other processors...receive the data...
                j=WalkersonNodes(i)*(NIfTot+1)
                CALL MPIRecv(Parts(0:NIfTot,1:WalkersonNodes(i)),j,i,Tag,error)
!                WRITE(6,*) "Recieved walkers for processor ",i
!                CALL FLUSH(6)
                
!Then write it out...
                IF(tBinPops) THEN
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            call extract_sign(Parts(:,j),TempSign)
                            WRITE(iunit) Parts(0:NIfDBO,j),TempSign(:)
                        ENDIF
                    enddo
                ELSE
                    do j=1,WalkersonNodes(i)
                        IF(mod(j,iPopsPartEvery).eq.0) THEN
                            do k=0,NIfDBO
                                WRITE(iunit,"(I24)",advance='no') Parts(k,j)
                            enddo
                            call extract_sign(Parts(:,j),TempSign)
                            WRITE(iunit,*) TempSign(:)
                        ENDIF
                    enddo
                ENDIF
!                WRITE(6,*) "Writted out walkers for processor ",i
!                CALL FLUSH(6)

            enddo

            CLOSE(iunit)

!Deallocate memory for temporary storage of information.
            DEALLOCATE(Parts)
            CALL LogMemDealloc(this_routine,PartsTag)

        ELSE
!All other processors need to send their data to root...
            j=nDets*(NIfTot+1)
            CALL MPISend(Dets(0:NIfTot,1:nDets),j,root,Tag,error)
!            WRITE(6,*) "Have sent info to head node..."
!            CALL FLUSH(6)
        ENDIF

!Reset the values of the global variables
        AllSumNoatHF = 0
        AllSumENum=0.D0
        AllTotWalkers=0.D0

        RETURN

    END SUBROUTINE WriteToPopsfileParOneArr

!This routine reads in particle configurations from a POPSFILE.
    SUBROUTINE ReadFromPopsfilePar()
        use util_mod, only: get_unique_filename, get_free_unit
        use CalcData, only: iPopsFileNoRead
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn,iWeightPopRead
        use Logging, only: tIncrementPops,tZeroProjE
        use constants, only: size_n_int,bits_n_int
        LOGICAL :: exists,tBinRead
        INTEGER :: AvWalkers,WalkerstoReceive(nProcessors)
        INTEGER*8 :: NodeSumNoatHF(nProcessors)
        integer(int64) :: TempTotParts(lenof_sign),TempCurrWalkers
        INTEGER :: TempInitWalkers,error,i,j,l,total,ierr,MemoryAlloc,Tag,Proc,CurrWalkers,ii
        INTEGER , DIMENSION(lenof_sign) :: TempSign
        INTEGER*8 :: iLutTemp64(0:nBasis/64+1)
        INTEGER :: iLutTemp32(0:nBasis/32+1)
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        INTEGER :: AvSumNoatHF,IntegerPart,TempnI(NEl),ExcitLevel
        INTEGER :: NIfWriteOut,pos,orb,PopsVersion, iunit
        REAL*8 :: r,FracPart,Gap,DiagSftTemp
        HElement_t :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        character(255) :: popsfile,FirstLine
        character(len=24) :: junk,junk2,junk3,junk4
        LOGICAL :: tPop64BitDets,tPopHPHF,tPopLz,tPopInitiator
        integer(n_int) :: ilut_largest(0:NIfTot)
        integer :: sign_largest

        IF(lenof_sign.ne.1) CALL Stop_All("ReadFromPopsfilePar","Popsfile does not work with complex walkers")
        
        PreviousCycles=0    !Zero previous cycles
        SumENum=0.D0
        TotParts=0
        SumNoatHF=0
        DiagSft=0.D0
        Tag=124             !Set Tag
        
        call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
        iunit = get_free_unit()
        INQUIRE(FILE=popsfile,EXIST=exists)
        IF(exists) THEN
            OPEN(iunit,FILE=popsfile,Status='old')
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
                    OPEN(iunit,FILE=popsfile,Status='old')
                ENDIF
            ENDIF
        ENDIF

        READ(iunit,'(a255)') FirstLine

        IF(INDEX(FirstLine,'VERSION').eq.0) THEN
!No version number to be found
            PopsVersion=1
            REWIND(iunit)
        ELSE
            !Found version - which number is it?
            REWIND(iunit)
            READ(iunit,*) FirstLine,FirstLine,FirstLine,PopsVersion
        ENDIF
        WRITE(6,"(A,I5,A)") "Version",PopsVersion," POPSFILE detected"

!Read in initial data on processors which have a popsfile
        IF(PopsVersion.eq.2) THEN
            READ(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') junk,tPop64BitDets,junk2,tPopHPHF,junk3,tPopLz,junk4,tPopInitiator
        ELSE
            WRITE(6,'(A)') "Reading in from depreciated POPSFILE - assuming that parameters are the same as when POPSFILE was written"
        ENDIF
        READ(iunit,*) AllTotWalkers
        READ(iunit,*) DiagSftTemp
        READ(iunit,*) AllSumNoatHF
        READ(iunit,*) AllSumENum
        READ(iunit,*) PreviousCycles

        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,I4,A)") "Although ",AllTotWalkers," configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            ENDIF
        ENDIF

        IF(.not.tWalkContGrow) THEN
!If we want the walker number to continue growing, then take the diagonal shift from the input, rather than the POPSFILE.
            DiagSft=DiagSftTemp
        ENDIF

        IF(DiagSftTemp.eq.0.D0) THEN
            tWalkContGrow=.true.
            DiagSft=DiagSftTemp
        ENDIF

        IF(tBinRead) THEN
!Test for the end of the file.
!If this is not the end of the file, there is one more keyword that tells us the calculation had not entered variable shift mode yet.
!Want to put this test at the end of the non-binary file too.
            CLOSE(iunit)
            call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            OPEN(iunit,FILE=popsfile,Status='old',form='unformatted')
        ENDIF

        IF(iProcIndex.eq.Root) THEN

            WRITE(6,*) "Number of cycles in previous simulation: ",PreviousCycles
            IF(NEquilSteps.gt.0) THEN
                WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
                NEquilSteps=0
            ENDIF
            IF(TZeroProjE) THEN
!Reset energy estimator
                WRITE(6,*) "Resetting projected energy counters to zero..."
                AllSumENum=0.D0
                AllSumNoatHF = 0
            ENDIF

!Need to calculate the number of walkers each node will receive...
            AvWalkers=NINT(AllTotWalkers/real(nProcessors,dp))

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
            if (total /= AllTotWalkers) then
                CALL Stop_All("ReadFromPopsfilePar","All Walkers not accounted for when reading in from POPSFILE")
            endif
            
!InitWalkers needs to be reset for the culling criteria
            IF(.not.tWalkContGrow) THEN
!Now, let the total space allocated for storing walkers which have been read in to be equal to the initwalkers from the input file.
!                InitWalkers=AvWalkers
            ELSE
                TSinglePartPhase=.true.
            ENDIF
            SumENum=AllSumENum/REAL(nProcessors,dp)     !Divide up the SumENum over all processors
            AvSumNoatHF = AllSumNoatHF(1)/nProcessors !This is the average Sumnoathf
            do i=1,nProcessors-1
                NodeSumNoatHF(i)=INT(AvSumNoatHF,int64)
            enddo
            NodeSumNoatHF(nProcessors)=AllSumNoatHF(1)-INT((AvSumNoatHF*(nProcessors-1)),int64)

            ProjectionE=AllSumENum/real(AllSumNoatHF(1),dp)
                
!Reset the global variables
            AllSumENum=0.D0
            AllSumNoatHF = 0

        ENDIF

        CALL MPIBarrier(error)  !Sync

!Now we need to scatter the WalkerstoReceive to each node, and allocate the desired memory to each node...
!Broadcast info which needs to go to all processors
        CALL MPIBCast(DiagSft)
        CALL MPIBCast(SumENum)
        CALL MPIBCast(InitWalkers)
        CALL MPIBCast(NEquilSteps)
        CALL MPIBCast(NShiftEquilSteps)
        CALL MPIBCast(TSinglePartPhase)
!        CALL MPI_BCast(tChangenProcessors,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
!Scatter the number of walkers each node will receive to TempInitWalkers, and the SumNoatHF for each node which is distributed approximatly equally
        CALL MPIScatter(WalkerstoReceive,TempInitWalkers,error)
        CALL MPIScatter(NodeSumNoatHF,SumNoatHF(1),error)

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

        CALL MPIBarrier(error)  !Sync

!Allocate memory to hold walkers at least temporarily
        ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating WalkVecDets array.')
        CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),size_n_int,this_routine,WalkVecDetsTag,ierr)
        WalkVecDets(:,:)=0
        MemoryAlloc=(NIfTot+1)*MaxWalkersPart*size_n_int    !Memory Allocated in bytes

!Just allocating this here, so that the SpawnParts arrays can be used for sorting the determinants when using direct annihilation.
        WRITE(6,"(A,I12,A)") " Spawning vectors allowing for a total of ",MaxSpawned," particles to be spawned in any one iteration."
        ALLOCATE(SpawnVec(0:NIfTot,MaxSpawned),stat=ierr)
        CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVecTag,ierr)
        SpawnVec(:,:)=0
        ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
        CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVec2Tag,ierr)
        SpawnVec2(:,:)=0
        MemoryAlloc=MemoryAlloc+MaxSpawned*(1+NIfTot)*2*size_n_int

!Point at correct spawning arrays
        SpawnedParts=>SpawnVec
        SpawnedParts2=>SpawnVec2
!Allocate pointer to the correct walker array...
        CurrentDets=>WalkVecDets

!Need to now allocate other arrays
        IF(.not.tRegenDiagHEls) THEN
            ALLOCATE(WalkVecH(MaxWalkersPart),stat=ierr)
            CALL LogMemAlloc('WalkVecH',MaxWalkersPart,8,this_routine,WalkVecHTag,ierr)
            WalkVecH(:)=0.d0
            MemoryAlloc=MemoryAlloc+8*MaxWalkersPart
        ELSE
            WRITE(6,"(A,F14.6,A)") "Diagonal H-Elements will not be stored. This will *save* ",REAL(MaxWalkersPart*8,dp)/1048576.D0," Mb/Processor"
        ENDIF

        IF(.not.tRegenDiagHEls) THEN
            CurrentH=>WalkVecH
        ENDIF

! The hashing will be different in the new calculation from the one where the POPSFILE was produced, this means we must recalculate the processor each determinant wants to go to.                
! This is done by reading in all walkers to the root and then distributing them in the same way as the spawning steps are done - by finding the determinant and sending it there.
        IF((PopsVersion.ne.1).and.tHPHF.and.(.not.tPopHPHF)) THEN
            CALL Stop_All(this_routine,"HPHF on, but HPHF was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.tFixLz.and.(.not.tPopLz)) THEN
            CALL Stop_All(this_routine,"Lz on, but Lz was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tHPHF).and.tPopHPHF) THEN
            CALL Stop_All(this_routine,"HPHF off, but HPHF was used for creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tFixLz).and.tPopLz) THEN
            CALL Stop_All(this_routine,"Lz off, but Lz was used for creation of the POPSFILE")
        ENDIF
        ! TODO: Add tests for CSFs here.
        IF(PopsVersion.eq.1) THEN
            tPop64BitDets=.false.
            NIfWriteOut=nBasis/32
            IF(tCSF) NIfWriteOut=NIfWriteOut+1
        ELSE
            IF(.not.tPop64BitDets) THEN
                NIfWriteOut=nBasis/32
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ELSE
                NIfWriteOut=nBasis/64
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ENDIF
        ENDIF

        CurrWalkers=0
        sign_largest = 0
        ilut_largest = 0
        do i=1,AllTotWalkers
            iLutTemp(:)=0
            IF(PopsVersion.ne.1) THEN
                IF(tBinRead) THEN
                    IF(tPop64BitDets) THEN
                        READ(iunit) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ELSE
                    IF(tPop64BitDets) THEN
                        READ(iunit,*) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ENDIF
            ELSE
                !POPSFILE v. 1 only printed out 32 bit determinant strings.
                IF(tBinRead) THEN
                    READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                ELSE
                    READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                ENDIF
            ENDIF

#ifdef __INT64
            if (.not.tPop64BitDets) then
                ! If we are using 64 bit integers, but have read in 32 bit 
                ! integers, then we need to convert them.
                do ii=0,nBasis/32
                    do j=0,31
                        if(btest(iLutTemp32(ii),j)) then
                            orb=(ii*32)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp64(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp64(0:NIfDBO)
            endif

#else
            ! If we are using 32 bit integers, but have read in 64 bit 
            ! integers, then we need to convert them.
            if (tPop64BitDets) then
                do ii=0,nBasis/64
                    do j=0,63
                        if(btest(iLutTemp64(ii),j)) then
                            orb=(ii*64)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp32(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp32(0:NIfDBO)
            endif
        
#endif
            call decode_bit_det (TempnI, iLutTemp)
            Proc=DetermineDetProc(iLutTemp)   !This wants to return a value between 0 -> nProcessors-1
            IF((Proc.eq.iProcIndex).and.(abs(TempSign(1)).ge.iWeightPopRead)) THEN
                CurrWalkers=CurrWalkers+1
                call encode_bit_rep(CurrentDets(:,CurrWalkers),iLutTemp(0:NIfDBO),TempSign,0)   !Do not need to send a flag here...
                                                                                                !TODO: Add flag for complex walkers to read in both
            ENDIF

            ! Keep track of what the most highly weighted determinant is
            if (abs(TempSign(1)) > sign_largest) then
                sign_largest = abs(TempSign(1))
                ilut_largest = iLutTemp
            endif
        enddo
        CLOSE(iunit)
        TempCurrWalkers=REAL(CurrWalkers,dp)

        ! Sort the lists so that they are in order if we change the number
        ! of processors.
        call sort (currentdets(:,1:CurrWalkers))

        ! Check that the bit-det comparisons agree that it is in order.
        do i=2,currwalkers
            if(DetBitLT(CurrentDets(:,i),CurrentDets(:,i-1),NIfDBO) == 1) then
                print*, 'Walkers: ', i-1, i
                print*, 'bit reps: '
                print*, currentdets(:, i-1)
                print*, currentdets(:, i)
                call stop_all (this_routine, 'Main list out of order')
            endif
        enddo

        CALL MPIBarrier(error)  !Sync
        CALL MPIAllReduce(TempCurrWalkers,MPI_SUM,AllTotWalkers)

        IF(iProcIndex.eq.root) WRITE(6,'(I10,A)') INT(AllTotWalkers,int64)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers by a factor of: ",ScaleWalkers

! CurrWalkers is the number of determinants on a particular node, AllTotWalkers is the total over all nodes.
            IntegerPart=INT(ScaleWalkers)
            FracPart=ScaleWalkers-REAL(IntegerPart)

            do l=1,CurrWalkers
                call extract_sign(CurrentDets(:,l),TempSign)
                TempSign=TempSign*IntegerPart
                r = genrand_real2_dSFMT() 
                IF(r.lt.FracPart) THEN
!Stochastically create another particle
                    IF(TempSign(1).lt.0) THEN
                        TempSign(1)=TempSign(1)-1
                    ELSE
                        TempSign(1)=TempSign(1)+1
                    ENDIF
                ENDIF
                call encode_sign(CurrentDets(:,l),TempSign)
            enddo

            InitWalkers=NINT(InitWalkers*ScaleWalkers)  !New (average) number of initial particles for culling criteria
!Other parameters don't change (I think) because the number of determinants isn't changing.                
            TotWalkers=CurrWalkers
            TotWalkersOld=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                aLlTotWalkersOld=AllTotWalkers
                iter_data_fciqmc%tot_parts_old = AllTotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ELSE
!We are not scaling the number of walkers...

            TotWalkers=CurrWalkers
            TotWalkersOld=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                AllTotWalkersOld=AllTotWalkers
                iter_data_fciqmc%tot_parts_old = AllTotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ENDIF
            
        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft
        WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,dp)/1048576.D0," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        WRITE(6,*) "Excitgens will be regenerated when they are needed..."
        CALL FLUSH(6)

        ! If we are changing the reference determinant to the largest
        ! weighted one in the file, do it here
        if (tReadPopsChangeRef .or. tReadPopsRestart) then
            if (.not. DetBitEq(ilut_largest, iLutRef, NIfDBO)) then
                
                ! Set new reference
                iLutRef = ilut_largest
                call decode_bit_det (ProjEDet, iLutRef)
                tNoBrillouin = .true.

                ! Recalculate the reference E
                if (tHPHF) then
                    HElemTemp = hphf_diag_helement (ProjEDet, iLutRef)
                else
                    HElemTemp = get_helement (ProjEDet, ProjEDet, 0)
                endif
                Hii = real(HElemTemp, dp)

                ! Output info on root node.
                if (iProcIndex == root) then
                    write(6, '(a)', advance='no') &
                        "Changing projected energy reference determinant to: "
                    call write_det (6, ProjEDet, .true.)
                    write (6, '(a)') &
                        "Ensuring that Brillouin's theorem is no longer used."
                    write (6, '(a,g25.15)') &
                        "Reference energy now set to: ", Hii
                endif
            endif
        endif

!Now find out the data needed for the particles which have been read in...
        TotParts=0
        do j=1,TotWalkers
            call decode_bit_det (TempnI, currentDets(:,j))
            Excitlevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j), 2)
            IF(Excitlevel.eq.0) THEN
                IF(.not.tRegenDiagHEls) CurrentH(j)=0.D0
            ELSE
                IF(.not.tRegenDiagHEls) THEN
                    if (tHPHF) then
                        HElemTemp = hphf_diag_helement (TempnI, &
                                                        CurrentDets(:,j))
                    else
                        HElemTemp = get_helement (TempnI, TempnI, 0)
                    endif
                    CurrentH(j)=REAL(HElemTemp,dp)-Hii
                ENDIF

            ENDIF
            call extract_sign(CurrentDets(:,j),TempSign)
            TotParts=TotParts+abs(TempSign(1))

        enddo

        TempTotParts=REAL(TotParts,dp)

        CALL MPIBarrier(error)  !Sync
        CALL MPIReduce(TempTotParts,MPI_SUM,AllTotParts)

        IF(iProcIndex.eq.root) AllTotPartsOld=AllTotParts
        write(6,'(A,i20)') ' The total number of particles read from the POPSFILE is: ',AllTotParts(1)

        if (tReadPopsRestart) then
            tPopsAlreadyRead = .true.
            call ChangeRefDet (ProjEDet)
            tPopsAlreadyRead = .false.
        endif

    END SUBROUTINE ReadFromPopsfilePar

!This routine reads in particle configurations from a POPSFILE.
! It's a bastardisation of ReadFromPopsfilePar but only reads a popsfile into an array and leaves the rest until later.
! nDets comes in with the max number of amplitudes possible in Dets.
! Dets is then filled with all the amplitudes in the POPSFILE and nDets is returned with that number.

! This will very likely still be tweaking the inner heart-strings of FciMCPar.  Caveatis stulti!  AJWT
! GHB says he will incorporate this functionality into a rewrite of ReadFromPopsfilePar. 19/8/2010
    SUBROUTINE ReadFromPopsfileOnly(Dets,nDets)
        use util_mod, only: get_unique_filename, get_free_unit
        use CalcData, only: iPopsFileNoRead
        use CalcData , only : MemoryFacPart,MemoryFacAnnihil,MemoryFacSpawn,iWeightPopRead
        use Logging, only: tIncrementPops,tZeroProjE
        use constants, only: size_n_int,bits_n_int
        integer(int64),intent(inout) :: nDets !The number of occupied entries in Dets
        integer(kind=n_int),intent(out) :: Dets(0:nIfTot,1:nDets)
        LOGICAL :: exists,tBinRead
        integer(int64) :: TempTotParts(lenof_sign),TempCurrWalkers
        INTEGER :: TempInitWalkers,error,i,j,l,total,ierr,MemoryAlloc,Tag,Proc,CurrWalkers,ii
        INTEGER , DIMENSION(lenof_sign) :: TempSign
        INTEGER*8 :: iLutTemp64(0:nBasis/64+1)
        INTEGER :: iLutTemp32(0:nBasis/32+1)
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        INTEGER :: AvSumNoatHF,IntegerPart,TempnI(NEl),ExcitLevel
        INTEGER :: NIfWriteOut,pos,orb,PopsVersion, iunit
        REAL*8 :: r,FracPart,Gap,DiagSftTemp
        HElement_t :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        character(255) :: popsfile,FirstLine
        character(len=24) :: junk,junk2,junk3,junk4
        LOGICAL :: tPop64BitDets,tPopHPHF,tPopLz,tPopInitiator
        integer(n_int) :: ilut_largest(0:NIfTot)
        integer :: sign_largest

        IF(lenof_sign.ne.1) CALL Stop_All("ReadFromPopsfilePar","Popsfile does not work with complex walkers")
        
        PreviousCycles=0    !Zero previous cycles
        SumENum=0.D0
        TotParts=0
        SumNoatHF=0
        DiagSft=0.D0
        Tag=124             !Set Tag
        
        call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
        iunit = get_free_unit()
        INQUIRE(FILE=popsfile,EXIST=exists)
        IF(exists) THEN
            OPEN(iunit,FILE=popsfile,Status='old')
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
                    OPEN(iunit,FILE=popsfile,Status='old')
                ENDIF
            ENDIF
        ENDIF

        READ(iunit,'(a255)') FirstLine

        IF(INDEX(FirstLine,'VERSION').eq.0) THEN
!No version number to be found
            PopsVersion=1
            REWIND(iunit)
        ELSE
            !Found version - which number is it?
            REWIND(iunit)
            READ(iunit,*) FirstLine,FirstLine,FirstLine,PopsVersion
        ENDIF
        WRITE(6,"(A,I5,A)") "Version",PopsVersion," POPSFILE detected"



!Read in initial data on processors which have a popsfile
        IF(PopsVersion.eq.2) THEN
            READ(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') junk,tPop64BitDets,junk2,tPopHPHF,junk3,tPopLz,junk4,tPopInitiator
        ELSE
            WRITE(6,'(A)') "Reading in from depreciated POPSFILE - assuming that parameters are the same as when POPSFILE was written"
        ENDIF
        READ(iunit,*) AllTotWalkers
        READ(iunit,*) DiagSftTemp
        READ(iunit,*) AllSumNoatHF
        READ(iunit,*) AllSumENum
        READ(iunit,*) PreviousCycles

        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,I4,A)") "Although ",AllTotWalkers," configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            ENDIF
        ENDIF

        IF(.not.tWalkContGrow) THEN
!If we want the walker number to continue growing, then take the diagonal shift from the input, rather than the POPSFILE.
            DiagSft=DiagSftTemp
        ENDIF

        IF(DiagSftTemp.eq.0.D0) THEN
            tWalkContGrow=.true.
            DiagSft=DiagSftTemp
        ENDIF

        IF(tBinRead) THEN
!Test for the end of the file.
!If this is not the end of the file, there is one more keyword that tells us the calculation had not entered variable shift mode yet.
!Want to put this test at the end of the non-binary file too.
            CLOSE(iunit)
            call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            OPEN(iunit,FILE=popsfile,Status='old',form='unformatted')
        ENDIF

        IF(iProcIndex.eq.Root) THEN

            WRITE(6,*) "Number of cycles in previous simulation: ",PreviousCycles
            IF(NEquilSteps.gt.0) THEN
                WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
                NEquilSteps=0
            ENDIF
            IF(TZeroProjE) THEN
!Reset energy estimator
                WRITE(6,*) "Resetting projected energy counters to zero..."
                AllSumENum=0.D0
                AllSumNoatHF = 0
            ENDIF

!InitWalkers needs to be reset for the culling criteria
            IF(.not.tWalkContGrow) THEN
!Now, let the total space allocated for storing walkers which have been read in to be equal to the initwalkers from the input file.
!                InitWalkers=AvWalkers
            ELSE
                TSinglePartPhase=.true.
            ENDIF
            SumENum=AllSumENum/REAL(nProcessors,dp)     !Divide up the SumENum over all processors
            AvSumNoatHF = AllSumNoatHF(1)/nProcessors !This is the average Sumnoathf

            ProjectionE=AllSumENum/real(AllSumNoatHF(1),dp)
                
!Reset the global variables
            AllSumENum=0.D0
            AllSumNoatHF = 0

        ENDIF

        CALL MPIBarrier(error)  !Sync

!Now we need to scatter the WalkerstoReceive to each node, and allocate the desired memory to each node...
!Broadcast info which needs to go to all processors
        CALL MPIBCast(DiagSft)
        CALL MPIBCast(SumENum)
        CALL MPIBCast(InitWalkers)
        CALL MPIBCast(NEquilSteps)
        CALL MPIBCast(NShiftEquilSteps)
        CALL MPIBCast(TSinglePartPhase)
!        CALL MPI_BCast(tChangenProcessors,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
!Scatter the number of walkers each node will receive to TempInitWalkers, and the SumNoatHF for each node which is distributed approximatly equally

        IF(MemoryFacPart.le.1.D0) THEN
            WRITE(6,*) 'MemoryFacPart must be larger than 1.0 when reading in a POPSFILE - increasing it to 1.50.'
            MemoryFacPart=1.50
        ENDIF
        
        CALL MPIBarrier(error)  !Sync

        if(AllTotWalkers>nDets) CALL Stop_All(this_routine,'Not enough memory to read in POPSFILE.')

! The hashing will be different in the new calculation from the one where the POPSFILE was produced, this means we must recalculate the processor each determinant wants to go to.                
! This is done by reading in all walkers to the root and then distributing them in the same way as the spawning steps are done - by finding the determinant and sending it there.
        IF((PopsVersion.ne.1).and.tHPHF.and.(.not.tPopHPHF)) THEN
            CALL Stop_All(this_routine,"HPHF on, but HPHF was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.tFixLz.and.(.not.tPopLz)) THEN
            CALL Stop_All(this_routine,"Lz on, but Lz was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tHPHF).and.tPopHPHF) THEN
            CALL Stop_All(this_routine,"HPHF off, but HPHF was used for creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tFixLz).and.tPopLz) THEN
            CALL Stop_All(this_routine,"Lz off, but Lz was used for creation of the POPSFILE")
        ENDIF
        ! TODO: Add tests for CSFs here.
        IF(PopsVersion.eq.1) THEN
            tPop64BitDets=.false.
            NIfWriteOut=nBasis/32
            IF(tCSF) NIfWriteOut=NIfWriteOut+1
        ELSE
            IF(.not.tPop64BitDets) THEN
                NIfWriteOut=nBasis/32
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ELSE
                NIfWriteOut=nBasis/64
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ENDIF
        ENDIF

        CurrWalkers=0
        sign_largest = 0
        ilut_largest = 0
        write(6,*) "Reading in ", AllTotWalkers, " walkers"
        do i=1,AllTotWalkers
            iLutTemp(:)=0
            IF(PopsVersion.ne.1) THEN
                IF(tBinRead) THEN
                    IF(tPop64BitDets) THEN
                        READ(iunit) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ELSE
                    IF(tPop64BitDets) THEN
                        READ(iunit,*) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ENDIF
            ELSE
                !POPSFILE v. 1 only printed out 32 bit determinant strings.
                IF(tBinRead) THEN
                    READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                ELSE
                    READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                ENDIF
            ENDIF

#ifdef __INT64
            if (.not.tPop64BitDets) then
                ! If we are using 64 bit integers, but have read in 32 bit 
                ! integers, then we need to convert them.
                do ii=0,nBasis/32
                    do j=0,31
                        if(btest(iLutTemp32(ii),j)) then
                            orb=(ii*32)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp64(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp64(0:NIfDBO)
            endif

#else
            ! If we are using 32 bit integers, but have read in 64 bit 
            ! integers, then we need to convert them.
            if (tPop64BitDets) then
                do ii=0,nBasis/64
                    do j=0,63
                        if(btest(iLutTemp64(ii),j)) then
                            orb=(ii*64)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp32(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp32(0:NIfDBO)
            endif
        
#endif
            call decode_bit_det (TempnI, iLutTemp)
            Proc=0  !DetermineDetProc(iLutTemp)   !This wants to return a value between 0 -> nProcessors-1
            IF((Proc.eq.iProcIndex).and.(abs(TempSign(1)).ge.iWeightPopRead)) THEN
                CurrWalkers=CurrWalkers+1
                call encode_bit_rep(Dets(:,CurrWalkers),iLutTemp(0:NIfDBO),TempSign,0)   !Do not need to send a flag here...
                                                                                                !TODO: Add flag for complex walkers to read in both
            ENDIF

            ! Keep track of what the most highly weighted determinant is
            if (abs(TempSign(1)) > sign_largest) then
                sign_largest = abs(TempSign(1))
                ilut_largest = iLutTemp
            endif
        enddo
        CLOSE(iunit)
        TempCurrWalkers=REAL(CurrWalkers,dp)

        ! Sort the lists so that they are in order if we change the number
        ! of processors.
        call sort(Dets(:,1:CurrWalkers))

        ! Check that the bit-det comparisons agree that it is in order.
        do i=2,currwalkers
            if(DetBitLT(Dets(:,i),Dets(:,i-1),NIfDBO) == 1) then
                print*, 'Walkers: ', i-1, i
                print*, 'bit reps: '
                print*, dets(:, i-1)
                print*, dets(:, i)
                call stop_all (this_routine, 'Main list out of order')
            endif
        enddo

        CALL MPIBarrier(error)  !Sync
        CALL MPIAllReduce(TempCurrWalkers,MPI_SUM,AllTotWalkers)

        IF(iProcIndex.eq.root) WRITE(6,'(I10,A)') INT(AllTotWalkers,int64)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers by a factor of: ",ScaleWalkers

! CurrWalkers is the number of determinants on a particular node, AllTotWalkers is the total over all nodes.
            IntegerPart=INT(ScaleWalkers)
            FracPart=ScaleWalkers-REAL(IntegerPart)

            do l=1,CurrWalkers
                call extract_sign(Dets(:,l),TempSign)
                TempSign=TempSign*IntegerPart
                r = genrand_real2_dSFMT() 
                IF(r.lt.FracPart) THEN
!Stochastically create another particle
                    IF(TempSign(1).lt.0) THEN
                        TempSign(1)=TempSign(1)-1
                    ELSE
                        TempSign(1)=TempSign(1)+1
                    ENDIF
                ENDIF
                call encode_sign(Dets(:,l),TempSign)
            enddo

            InitWalkers=NINT(InitWalkers*ScaleWalkers)  !New (average) number of initial particles for culling criteria
!Other parameters don't change (I think) because the number of determinants isn't changing.                
            nDets=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ELSE
!We are not scaling the number of walkers...

            nDets=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ENDIF
            
        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft
        WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,dp)/1048576.D0," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        WRITE(6,*) "Excitgens will be regenerated when they are needed..."
        CALL FLUSH(6)

!Now find out the data needed for the particles which have been read in...
        TotParts=0
        do j=1,nDets
            call extract_sign(Dets(:,j),TempSign)
            TotParts=TotParts+abs(TempSign(1))
        enddo

        TempTotParts=REAL(TotParts,dp)

        CALL MPIBarrier(error)  !Sync
        CALL MPIReduce(TempTotParts,MPI_SUM,AllTotParts)

        IF(iProcIndex.eq.root) AllTotPartsOld=AllTotParts
        write(6,'(A,i20)') ' The total number of particles read from the POPSFILE is: ',AllTotParts(1)

        if (tReadPopsRestart) then
            tPopsAlreadyRead = .true.
            call ChangeRefDet (ProjEDet)
            tPopsAlreadyRead = .false.
        endif

    END SUBROUTINE ReadFromPopsfileOnly

END MODULE PopsfileMod
