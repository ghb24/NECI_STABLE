! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
!This is a parallel MPI version of the FciMC code.
!All variables refer to values per processor

MODULE FciMCLoggingMod

    USE Global_utilities
    USE Parallel_neci
    USE Logging , only : tSaveBlocking,tBlockEveryIteration,HistInitPops,HistInitPopsTag,AllHistInitPops,AllHistInitPopsTag
    use SystemData, only: NEl
    use bit_reps, only: NIfTot, NIfDBO
    USE SymData , only : nSymLabels
    USE Determinants , only : get_helement, get_helement_excit
    USE CalcData , only : NMCyc,StepsSft
    use DetBitOps, only: DetBitEQ, FindExcitBitDet, FindBitExcitLevel
    use constants, only: dp,n_int
    use MemoryManager, only: TagIntType

    IMPLICIT NONE
    save

    real(dp) :: NoNotAccept,NoAccept,TotHElNotAccept,TotHElAccept,MaxHElNotAccept,MinHElAccept
    real(dp) :: NoPosSpinCoup,NoNegSpinCoup,SumPosSpinCoup,SumNegSpinCoup,SumHFCon,SumSpinCon,InitBinMin,InitBinIter,InitBinMax

    real(dp) , ALLOCATABLE :: CurrBlockSum(:),BlockSum(:),BlockSqrdSum(:)
    real(dp) , ALLOCATABLE :: CurrShiftBlockSum(:),ShiftBlockSum(:),ShiftBlockSqrdSum(:)
    INTEGER(TagIntType) :: CurrBlockSumTag,BlockSumTag,BlockSqrdSumTag
    INTEGER :: TotNoBlockSizes,StartBlockIter
    INTEGER(TagIntType) :: CurrShiftBlockSumTag,ShiftBlockSumTag,ShiftBlockSqrdSumTag
    INTEGER :: TotNoShiftBlockSizes,StartShiftBlockIter




    contains

    SUBROUTINE InitErrorBlocking(Iter)
        CHARACTER(len=*), PARAMETER :: this_routine='InitErrorBlocking'
        INTEGER :: ierr,Iter

! First want to find out how many different block sizes we will get if the calculation goes to completion. 
! The number of iterations we are doing the blocking for is NMCYC - the current iteration.
! The current iteration is the iteration when this routine is called - i.e. when the blocking analysis starts - the default will
! be when the HF population reaches 100.  This prevents problems faced when the HF walker dies.

! The block sizes increase in powers of 2.  The number of blocks is therefore log base 2 of the number of iterations.

! log_2(Iteration) = log_10(Iteration) / log_10(2)

        StartBlockIter=Iter 
        ! This is the iteration the blocking was started.

        IF(tBlockEveryIteration) THEN
            TotNoBlockSizes=FLOOR( (LOG10(REAL(NMCyc-StartBlockIter))) / (LOG10(2.0_dp)) )
        ELSE
            TotNoBlockSizes=FLOOR( (LOG10((REAL(NMCyc-StartBlockIter))/(REAL(StepsSft)))) / (LOG10(2.0_dp)) )
        ENDIF
        WRITE(6,*) 'Beginning blocking analysis of the errors in the projected energies.'
        WRITE(6,"(A,I6)") "The total number of different block sizes possible is: ",TotNoBlockSizes
        ! The blocks will have size 1,2,4,8,....,2**TotNoBlockSizes
        ! In the below arrays, the element i will correspond to block size 2**(i), 
        !but the arrays go from 0 -> TotNoBlockSizes.

! Then need to allocate three arrays of this size - one for the current block, 
!one for the sum of the blocks over the elapsed cycles,
! and one for the sum of the squares of the blocks over the elapsed iterations.

        ALLOCATE(CurrBlockSum(0:TotNoBlockSizes),stat=ierr)
        CALL LogMemAlloc('CurrBlockSum',TotNoBlockSizes,8,this_routine,CurrBlockSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating CurrBlockSum.')
        CurrBlockSum(:)=0.0_dp
 
        ALLOCATE(BlockSum(0:TotNoBlockSizes),stat=ierr)
        CALL LogMemAlloc('BlockSum',TotNoBlockSizes,8,this_routine,BlockSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating BlockSum.')
        BlockSum(:)=0.0_dp
 
        ALLOCATE(BlockSqrdSum(0:TotNoBlockSizes),stat=ierr)
        CALL LogMemAlloc('BlockSqrdSum',TotNoBlockSizes,8,this_routine,BlockSqrdSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating BlockSqrdSum.')
        BlockSqrdSum(:)=0.0_dp

    END SUBROUTINE InitErrorBlocking




    SUBROUTINE InitShiftErrorBlocking(Iter)
        CHARACTER(len=*), PARAMETER :: this_routine='InitShiftErrorBlocking'
        INTEGER :: ierr,Iter


        StartShiftBlockIter=Iter 
        ! This is the iteration the blocking was started.

        TotNoShiftBlockSizes=FLOOR( (LOG10((REAL(NMCyc-StartShiftBlockIter))/(REAL(StepsSft)))) / (LOG10(2.0_dp)) )
        WRITE(6,*) 'Beginning blocking analysis of the errors in the shift.'
        WRITE(6,"(A,I6)") "The total number of different block sizes possible is: ",TotNoShiftBlockSizes
        ! The blocks will have size 1,2,4,8,....,2**TotNoBlockSizes
        ! In the below arrays, the element i will correspond to block size 2**(i), 
        !but the arrays go from 0 -> TotNoBlockSizes.

! Then need to allocate three arrays of this size - one for the current block, 
!one for the sum of the blocks over the elapsed cycles,
! and one for the sum of the squares of the blocks over the elapsed iterations.

        ALLOCATE(CurrShiftBlockSum(0:TotNoShiftBlockSizes),stat=ierr)
        CALL LogMemAlloc('CurrShiftBlockSum',TotNoShiftBlockSizes,8,this_routine,CurrShiftBlockSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating CurrShiftBlockSum.')
        CurrShiftBlockSum(:)=0.0_dp
 
        ALLOCATE(ShiftBlockSum(0:TotNoShiftBlockSizes),stat=ierr)
        CALL LogMemAlloc('ShiftBlockSum',TotNoShiftBlockSizes,8,this_routine,ShiftBlockSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating ShiftBlockSum.')
        ShiftBlockSum(:)=0.0_dp
 
        ALLOCATE(ShiftBlockSqrdSum(0:TotNoShiftBlockSizes),stat=ierr)
        CALL LogMemAlloc('ShiftBlockSqrdSum',TotNoShiftBlockSizes,8,this_routine,ShiftBlockSqrdSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating ShiftBlockSqrdSum.')
        ShiftBlockSqrdSum(:)=0.0_dp



    END SUBROUTINE InitShiftErrorBlocking

 

    SUBROUTINE RestartBlocking(Iter)
        INTEGER :: Iter

        StartBlockIter=0
        StartBlockIter=Iter+StepsSft
        ! ChangeVars gets called at the end of the run, wont actually start until the next iteration.

        TotNoBlockSizes=0
        IF(tBlockEveryIteration) THEN
            TotNoBlockSizes=FLOOR( (LOG10(REAL(NMCyc-StartBlockIter))) / (LOG10(2.0_dp)) )
        ELSE
            TotNoBlockSizes=FLOOR( (LOG10((REAL(NMCyc-StartBlockIter))/(REAL(StepsSft)))) / (LOG10(2.0_dp)) )
        ENDIF

        if (allocated(CurrBlockSum)) CurrBlockSum(:)=0.0_dp
        if (allocated(BlockSum)) BlockSum(:)=0.0_dp
        if (allocated(BlockSqrdSum)) BlockSqrdSum(:)=0.0_dp

    END SUBROUTINE RestartBlocking
 


    SUBROUTINE RestartShiftBlocking(Iter)
        INTEGER :: Iter

        StartShiftBlockIter=0
        StartShiftBlockIter=Iter+StepsSft
        ! ChangeVars gets called at the end of the run, wont actually start until the next iteration.

        TotNoShiftBlockSizes=0
        TotNoShiftBlockSizes=FLOOR( (LOG10((REAL(NMCyc-StartShiftBlockIter))/(REAL(StepsSft)))) / (LOG10(2.0_dp)) )

        if (allocated(CurrShiftBlockSum)) CurrShiftBlockSum(:)=0.0_dp
        if (allocated(ShiftBlockSum)) ShiftBlockSum(:)=0.0_dp
        if (allocated(ShiftBlockSqrdSum)) ShiftBlockSqrdSum(:)=0.0_dp

    END SUBROUTINE RestartShiftBlocking



    SUBROUTINE SumInErrorContrib(Iter,AllENumCyc,AllHFCyc)
        INTEGER :: i,Iter,NoContrib
        real(dp) :: AllENumCyc,AllHFCyc

! First we need to find out what number contribution to the blocking this is.
        IF(tBlockEveryIteration) THEN
            NoContrib=(Iter-StartBlockIter)+1
        ELSE
            NoContrib=((Iter-StartBlockIter)/StepsSft)+1
        ENDIF

! We then need to add the contribution from this cycle to all the relevant blocks.        
        
        do i=0,TotNoBlockSizes

! First of all sum the energy contributions into each of the blocks.
            CurrBlockSum(i)=CurrBlockSum(i)+(AllENumCyc/AllHFCyc)
            ! The values contained in this array is the sum of the energy contributions to that block.
            ! They have not been divided by the number of contributions yet.

            ! Find out if we have completed a block.
            ! If we have, then add the average value from that block to the overall sum array (BlockSum),
            ! and zero the CurrBlockSum.
            IF(MOD(NoContrib,(2**i)).eq.0) THEN
                BlockSum(i)=BlockSum(i)+(CurrBlockSum(i)/(2**i))
                BlockSqrdSum(i)=BlockSqrdSum(i)+((CurrBlockSum(i)/(2**i))**2)
                CurrBlockSum(i)=0.0_dp
            ENDIF

        enddo

    END SUBROUTINE SumInErrorContrib



    SUBROUTINE SumInShiftErrorContrib(Iter,DiagSft)
        INTEGER :: i,Iter,NoContrib
        real(dp) :: DiagSft

! First we need to find out what number contribution to the blocking this is.
        NoContrib=((Iter-StartShiftBlockIter)/StepsSft)+1

! We then need to add the contribution from this cycle to all the relevant blocks.        
        
        do i=0,TotNoShiftBlockSizes

! First of all sum the energy contributions into each of the blocks.
            CurrShiftBlockSum(i)=CurrShiftBlockSum(i)+DiagSft
            ! The values contained in this array is the sum of the energy contributions to that block.
            ! They have not been divided by the number of contributions yet.

            ! Find out if we have completed a block.
            ! If we have, then add the average value from that block to the overall sum array (BlockSum),
            ! and zero the CurrBlockSum.
            IF(MOD(NoContrib,(2**i)).eq.0) THEN
                ShiftBlockSum(i)=ShiftBlockSum(i)+(CurrShiftBlockSum(i)/(2**i))
                ShiftBlockSqrdSum(i)=ShiftBlockSqrdSum(i)+((CurrShiftBlockSum(i)/(2**i))**2)
                CurrShiftBlockSum(i)=0.0_dp
            ENDIF

        enddo

    END SUBROUTINE SumInShiftErrorContrib



    SUBROUTINE PrintBlocking(Iter)
        use util_mod, only: get_free_unit
        INTEGER :: i,NoBlocks,Iter,NoBlockSizes,NoContrib, iunit
        real(dp) :: MeanEn,MeanEnSqrd,StandardDev,Error,ErrorinError
        CHARACTER(len=30) :: abstr

!First find out how many blocks would have been formed with the number of iterations actually performed. 

        NoContrib=0
        IF(tBlockEveryIteration) THEN
            NoContrib=(Iter-StartBlockIter)+1
        ELSE
            NoContrib=((Iter-StartBlockIter)/StepsSft)+1
        ENDIF
        IF(NoContrib.eq.1) RETURN

        NoBlockSizes=0
        NoBlockSizes=FLOOR( (LOG10(REAL(NoContrib-1)))/ (LOG10(2.0_dp)))

        iunit = get_free_unit()
        IF(tSaveBlocking) THEN
!We are saving the blocking files - write to new file.
            abstr=''
            write(abstr,'(I12)') Iter
            abstr='BLOCKINGANALYSIS-'//adjustl(abstr)
            OPEN(iunit,file=abstr,status='unknown')
        ELSE
            OPEN(iunit,file='BLOCKINGANALYSIS',status='unknown')
        ENDIF

        WRITE(iunit,'(I4,A,I4)') NoBlockSizes,' blocks were formed with sizes from 1 to ',(2**(NoBlockSizes))
        WRITE(iunit,'(3A16,5A20)') '1.Block No.','2.Block Size  ','3.No. of Blocks','4.Mean E', &
            '5.Mean E^2','6.SD','7.Error','8.ErrorinError'

        do i=0,NoBlockSizes-1
            
            ! First need to find out how many blocks of this particular size contributed to the final sum in BlockSum.
            ! NoContrib is the total number of contributions to the blocking throughout the simulation.
            NoBlocks=0
            NoBlocks=FLOOR(REAL(NoContrib/(2**i)))
            ! This finds the lowest integer multiple of 2**i (the block size). 

            MeanEn=BlockSum(i)/REAL(NoBlocks)
            MeanEnSqrd=BlockSqrdSum(i)/REAL(NoBlocks)

            StandardDev=SQRT(MeanEnSqrd-(MeanEn**2))
            IF(StandardDev.eq.0.0_dp) THEN
                Error=0.0_dp
                ErrorinError=0.0_dp
            ELSE
                Error=StandardDev/SQRT(REAL(NoBlocks-1))
                ErrorinError=Error/SQRT(2.0_dp*(REAL(NoBlocks-1)))
            ENDIF

!This is from the blocking paper, and indicates the error in the blocking error, due to the limited number of blocks available.

            WRITE(iunit,'(3I16,5F20.10)') i,(2**i),NoBlocks,MeanEn,MeanEnSqrd,StandardDev,Error,ErrorinError
            
        enddo

        CLOSE(iunit)


    END SUBROUTINE PrintBlocking


    SUBROUTINE PrintShiftBlocking(Iter)
        use util_mod, only: get_free_unit
        INTEGER :: i,NoBlocks,Iter,NoBlockSizes,NoContrib,iunit
        real(dp) :: MeanShift,MeanShiftSqrd,StandardDev,Error,ErrorinError

!First find out how many blocks would have been formed with the number of iterations actually performed. 

        NoContrib=0
        NoContrib=((Iter-StartShiftBlockIter)/StepsSft)+1

        NoBlockSizes=0
        NoBlockSizes=FLOOR( (LOG10(REAL(NoContrib-1)))/ (LOG10(2.0_dp)))

        iunit = get_free_unit()
        OPEN(iunit,file='SHIFTBLOCKINGANALYSIS',status='unknown')
        WRITE(iunit,'(I4,A,I4)') NoBlockSizes,' blocks were formed with sizes from 1 to ',(2**(NoBlockSizes))
        WRITE(iunit,'(3A16,5A20)') '1.Block No.','2.Block Size  ','3.No. of Blocks','4.Mean Shift', &
            '5.Mean Shift^2','6.SD','7.Error','8.ErrorinError'

        do i=0,NoBlockSizes-1
            
            ! First need to find out how many blocks of this particular size contributed to the final sum in BlockSum.
            ! NoContrib is the total number of contributions to the blocking throughout the simulation.

            NoBlocks=0
            NoBlocks=FLOOR(REAL(NoContrib/(2**i)))
            ! This finds the lowest integer multiple of 2**i (the block size). 

            MeanShift=ShiftBlockSum(i)/REAL(NoBlocks)
            MeanShiftSqrd=ShiftBlockSqrdSum(i)/REAL(NoBlocks)

            StandardDev=SQRT(MeanShiftSqrd-(MeanShift**2))
            IF(StandardDev.eq.0.0_dp) THEN
                Error=0.0_dp
                ErrorinError=0.0_dp
            ELSE
                Error=StandardDev/SQRT(REAL(NoBlocks-1))
                ErrorinError=Error/SQRT(2.0_dp*(REAL(NoBlocks-1)))
            ENDIF

!This is from the blocking paper, and indicates the error in the blocking error, due to the limited number of blocks available.

            WRITE(iunit,'(3I16,5F20.10)') i,(2**i),NoBlocks,MeanShift,MeanShiftSqrd,StandardDev,Error,ErrorinError
            
        enddo

        CLOSE(iunit)

    END SUBROUTINE PrintShiftBlocking



    SUBROUTINE FinaliseBlocking(Iter)
        INTEGER :: Iter
        CHARACTER(len=*), PARAMETER :: this_routine='FinaliseBlocking'


        CALL PrintBlocking(Iter)

        DEALLOCATE(CurrBlockSum)
        CALL LogMemDeAlloc(this_routine,CurrBlockSumTag)
 
        DEALLOCATE(BlockSum)
        CALL LogMemDeAlloc(this_routine,BlockSumTag)
 
        DEALLOCATE(BlockSqrdSum)
        CALL LogMemDeAlloc(this_routine,BlockSqrdSumTag)
 

    END SUBROUTINE FinaliseBlocking



    SUBROUTINE FinaliseShiftBlocking(Iter)
        INTEGER :: Iter
        CHARACTER(len=*), PARAMETER :: this_routine='FinaliseShiftBlocking'


        CALL PrintShiftBlocking(Iter)

        DEALLOCATE(CurrShiftBlockSum)
        CALL LogMemDeAlloc(this_routine,CurrShiftBlockSumTag)
 
        DEALLOCATE(ShiftBlockSum)
        CALL LogMemDeAlloc(this_routine,ShiftBlockSumTag)
 
        DEALLOCATE(ShiftBlockSqrdSum)
        CALL LogMemDeAlloc(this_routine,ShiftBlockSqrdSumTag)
 

    END SUBROUTINE FinaliseShiftBlocking


    SUBROUTINE InitHistInitPops()
        USE CalcData , only : InitiatorWalkNo
        INTEGER :: ierr
        character(*), parameter :: this_routine = 'InitHistInitPops'
    

        if (allocated(HistInitPops)) then
            deallocate (HistInitPops, stat=ierr)
            call LogMemDealloc (this_routine, HistInitPopsTag, ierr)
        endif

        allocate (HistInitPops(2,25000), stat=ierr)
        call LogMemAlloc('HistInitPops', 50000, 4, this_routine, &
                         HistInitPopsTag, ierr)
        HistInitPops = 0

        if (iProcIndex.eq.0) then
            if (allocated(AllHistInitPops)) then
                deallocate (AllHistInitPops, stat=ierr)
                call LogMemDealloc (this_routine, AllHistInitPopsTag, ierr)
            endif

            allocate (AllHistInitPops(2,25000), stat=ierr)
            CALL LogMemAlloc('AllHistInitPops', 50000, 4, this_routine, &
                             AllHistInitPopsTag, ierr)
            AllHistInitPops = 0
        endif

        InitBinMin=log(REAL(InitiatorWalkNo+1))
        InitBinMax=log(1000000.0_dp)
        InitBinIter=ABS(InitBinMax-InitBinMin)/25000.0


    END SUBROUTINE InitHistInitPops


    SUBROUTINE WriteInitPops(Iter)
        use util_mod, only: get_free_unit
        CHARACTER(len=21) :: abstr
        INTEGER :: i,Iter,error, iunit
        real(dp) :: InitBinCurr

!This will open a file called InitPops-"Iter" on unit number 17.
        abstr=''
        write(abstr,'(I12)') Iter
        abstr='InitPops-'//adjustl(abstr)

        call MPIReduce(HistInitPops,MPI_SUM,AllHistInitPops)

        IF(iProcIndex.eq.0) THEN
            iunit = get_free_unit()
            OPEN(iunit,FILE=abstr,STATUS='unknown')

            InitBinCurr=(-1)*InitBinMax            
            do i=25000,1,-1
                IF(AllHistInitPops(1,i).ne.0) WRITE(iunit,'(F20.10,2I20)') &
                InitBinCurr,(-1)*(NINT(EXP(ABS(InitBinCurr)))),AllHistInitPops(1,i)
                InitBinCurr=InitBinCurr+InitBinIter
            enddo

            InitBinCurr=InitBinMin
            do i=1,25000
                IF(AllHistInitPops(2,i).ne.0) WRITE(iunit,'(F20.10,2I20)') &
                InitBinCurr,NINT(EXP(InitBinCurr)),AllHistInitPops(2,i)
                InitBinCurr=InitBinCurr+InitBinIter
            enddo
 
            CLOSE(iunit)
            AllHistInitPops(:,:)=0
        ENDIF
        HistInitPops(:,:)=0

    END SUBROUTINE WriteInitPops


    SUBROUTINE TrackSpawnAttempts(Child,DetCurr,nJ,IC,Ex,tParity)
        INTEGER :: Child,DetCurr(NEl),nJ(NEl),IC,Ex(2,2)
        LOGICAL :: tParity
        HElement_t :: HEl

!        WRITE(6,*) 'Child',Child
!        WRITE(6,*) 'DetCurr',DetCurr
!        WRITE(6,*) 'nJ',nJ
!        WRITE(6,*) 'iLutnJ',iLutnJ
!        CALL neci_flush(6)
!        stop

        ! Need to find the H element between the current determinant and that which we're trying to spawn on.
        HEl = get_helement_excit (DetCurr, nJ, IC, Ex, tParity)
            
        IF(Child.eq.0) THEN
            ! Spawn not accepted.
            NoNotAccept=NoNotAccept+1.0_dp
            TotHElNotAccept=TotHElNotAccept+ABS(REAL(HEl,dp))
            IF(ABS(REAL(HEl,dp)).gt.ABS(MaxHElNotAccept)) MaxHElNotAccept=ABS(REAL(HEl,dp))
        ELSE
            ! Spawn accepted.
            NoAccept=NoAccept+1.0_dp
            TotHElAccept=TotHElAccept+ABS(REAL(HEl,dp))
            IF((MinHElAccept.eq.0.0_dp).or.(ABS(REAL(HEl,dp)).lt.ABS(MinHElAccept))) MinHElAccept=ABS(REAL(HEl,dp))
        ENDIF

    ENDSUBROUTINE TrackSpawnAttempts
        

    SUBROUTINE PrintSpawnAttemptStats(Iteration)
        use util_mod, only: get_free_unit
        real(dp) :: AllStats(4),AcceptStats(4),AllMaxHElNotAccept(1:nProcessors),AllMinHElAccept(1:nProcessors)
        INTEGER :: i,error,Iteration, iunit

        ! Need to distribute the max and min values to all processors, but only if it has changed.        
        AcceptStats(1)=TotHElNotAccept       ! Total not accepted
        AcceptStats(2)=TotHElAccept          ! Total accepted
        AcceptStats(3)=NoNotAccept
        AcceptStats(4)=NoAccept
        AllStats(:)=0.0_dp
        AllMaxHElNotAccept(:)=0.0_dp
        AllMinHElAccept(:)=0.0_dp

        call MPIReduce(AcceptStats,MPI_Sum,AllStats)
        call MPIGather(MaxHElNotAccept,AllMaxHElNotAccept(1:nProcessors),error)
        call MPIGather(MinHElAccept,AllMinHElAccept(1:nProcessors),error)

        IF(iProcIndex.eq.Root) THEN 
!            WRITE(6,*) 'AllMinHElAccept',AllMinHElAccept
            CALL neci_flush(6)
            MaxHElNotAccept=ABS(AllMaxHElNotAccept(1))
            do i=2,nProcessors
                IF(ABS(AllMaxHElNotAccept(i)).gt.ABS(MaxHElNotAccept)) MaxHElNotAccept=ABS(AllMaxHElNotAccept(i))
            enddo

            MinHElAccept=0.0_dp
            IF(AllStats(4).gt.0.0_dp) THEN
                do i=1,nProcessors
                    IF(AllMinHElAccept(i).ne.0.0_dp) THEN
                        MinHElAccept=ABS(AllMinHElAccept(i))
                        EXIT
                    ENDIF
                enddo
                do i=1,nProcessors
                    IF((AllMinHElAccept(i).ne.0.0_dp).and.(ABS(AllMinHElAccept(i)).lt.ABS(MinHElAccept))) &
                    MinHElAccept=ABS(AllMinHElAccept(i))
                enddo
            ENDIF

            iunit = get_free_unit()
            open(iunit,file="SpawnAttemptStats",status="unknown")
            WRITE(iunit,'(I10,2F20.1,5F20.6)') Iteration,AllStats(3),AllStats(4), &
            AllStats(3)/AllStats(4),AllStats(1)/AllStats(3),AllStats(2)/AllStats(4),MaxHElNotAccept,MinHElAccept
            close(iunit)
        ENDIF

    ENDSUBROUTINE PrintSpawnAttemptStats

ENDMODULE FciMCLoggingMod

