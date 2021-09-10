!This is a parallel MPI version of the FciMC code.
!All variables refer to values per processor

MODULE FciMCLoggingMod

    USE Global_utilities
    USE Parallel_neci
    USE LoggingData, only: tSaveBlocking, tBlockEveryIteration, HistInitPops, HistInitPopsTag, AllHistInitPops, AllHistInitPopsTag
    use SystemData, only: NEl, tGUGA
    use bit_reps, only: NIfTot
    USE SymData, only: nSymLabels
    USE Determinants, only: get_helement, get_helement_excit
    use CalcData, only: NMCyc, StepsSft, InitiatorWalkNo
    use DetBitOps, only: DetBitEQ, FindExcitBitDet, FindBitExcitLevel
    use constants, only: dp, n_int, maxExcit
    use MemoryManager, only: TagIntType
    use FciMCData, only: HighPopNeg, HighPopPos, MaxInitPopNeg, MaxInitPopPos
    use fortran_strings, only: str

    IMPLICIT NONE
    save

    real(dp) :: NoNotAccept, NoAccept, TotHElNotAccept, TotHElAccept, MaxHElNotAccept, MinHElAccept
    real(dp) :: NoPosSpinCoup, NoNegSpinCoup, SumPosSpinCoup, SumNegSpinCoup, SumHFCon, SumSpinCon, InitBinMin, InitBinIter, InitBinMax

    real(dp), ALLOCATABLE :: CurrBlockSum(:), BlockSum(:), BlockSqrdSum(:)
    real(dp), ALLOCATABLE :: CurrShiftBlockSum(:), ShiftBlockSum(:), ShiftBlockSqrdSum(:)
    INTEGER(TagIntType) :: CurrBlockSumTag, BlockSumTag, BlockSqrdSumTag
    INTEGER :: TotNoBlockSizes, StartBlockIter
    INTEGER(TagIntType) :: CurrShiftBlockSumTag, ShiftBlockSumTag, ShiftBlockSqrdSumTag
    INTEGER :: TotNoShiftBlockSizes, StartShiftBlockIter

contains

    SUBROUTINE InitErrorBlocking(Iter)
        CHARACTER(len=*), PARAMETER :: this_routine = 'InitErrorBlocking'
        INTEGER :: ierr, Iter

! First want to find out how many different block sizes we will get if the calculation goes to completion.
! The number of iterations we are doing the blocking for is NMCYC - the current iteration.
! The current iteration is the iteration when this routine is called - i.e. when the blocking analysis starts - the default will
! be when the HF population reaches 100.  This prevents problems faced when the HF walker dies.

! The block sizes increase in powers of 2.  The number of blocks is therefore log base 2 of the number of iterations.

! log_2(Iteration) = log_10(Iteration) / log_10(2)

        StartBlockIter = Iter
        ! This is the iteration the blocking was started.

        IF (tBlockEveryIteration) THEN
            TotNoBlockSizes = FLOOR((LOG10(REAL(NMCyc - StartBlockIter, dp))) / (LOG10(2.0_dp)))
        ELSE
            TotNoBlockSizes = FLOOR((LOG10((REAL(NMCyc - StartBlockIter, dp)) / (REAL(StepsSft, dp)))) / (LOG10(2.0_dp)))
        end if
        write(stdout, *) 'Beginning blocking analysis of the errors in the projected energies.'
        write(stdout, "(A,I6)") "The total number of different block sizes possible is: ", TotNoBlockSizes
        ! The blocks will have size 1,2,4,8,....,2**TotNoBlockSizes
        ! In the below arrays, the element i will correspond to block size 2**(i),
        !but the arrays go from 0 -> TotNoBlockSizes.

! Then need to allocate three arrays of this size - one for the current block,
!one for the sum of the blocks over the elapsed cycles,
! and one for the sum of the squares of the blocks over the elapsed iterations.

        allocate(CurrBlockSum(0:TotNoBlockSizes), stat=ierr)
        CALL LogMemAlloc('CurrBlockSum', TotNoBlockSizes, 8, this_routine, CurrBlockSumTag, ierr)
        IF (ierr /= 0) CALL Stop_All(this_routine, 'Problem allocating CurrBlockSum.')
        CurrBlockSum(:) = 0.0_dp

        allocate(BlockSum(0:TotNoBlockSizes), stat=ierr)
        CALL LogMemAlloc('BlockSum', TotNoBlockSizes, 8, this_routine, BlockSumTag, ierr)
        IF (ierr /= 0) CALL Stop_All(this_routine, 'Problem allocating BlockSum.')
        BlockSum(:) = 0.0_dp

        allocate(BlockSqrdSum(0:TotNoBlockSizes), stat=ierr)
        CALL LogMemAlloc('BlockSqrdSum', TotNoBlockSizes, 8, this_routine, BlockSqrdSumTag, ierr)
        IF (ierr /= 0) CALL Stop_All(this_routine, 'Problem allocating BlockSqrdSum.')
        BlockSqrdSum(:) = 0.0_dp

    END SUBROUTINE InitErrorBlocking

    SUBROUTINE InitShiftErrorBlocking(Iter)
        CHARACTER(len=*), PARAMETER :: this_routine = 'InitShiftErrorBlocking'
        INTEGER :: ierr, Iter

        StartShiftBlockIter = Iter
        ! This is the iteration the blocking was started.

        TotNoShiftBlockSizes = FLOOR((LOG10((REAL(NMCyc - StartShiftBlockIter, dp)) / (REAL(StepsSft, dp)))) / (LOG10(2.0_dp)))
        write(stdout, *) 'Beginning blocking analysis of the errors in the shift.'
        write(stdout, "(A,I6)") "The total number of different block sizes possible is: ", TotNoShiftBlockSizes
        ! The blocks will have size 1,2,4,8,....,2**TotNoBlockSizes
        ! In the below arrays, the element i will correspond to block size 2**(i),
        !but the arrays go from 0 -> TotNoBlockSizes.

! Then need to allocate three arrays of this size - one for the current block,
!one for the sum of the blocks over the elapsed cycles,
! and one for the sum of the squares of the blocks over the elapsed iterations.

        allocate(CurrShiftBlockSum(0:TotNoShiftBlockSizes), stat=ierr)
        CALL LogMemAlloc('CurrShiftBlockSum', TotNoShiftBlockSizes, 8, this_routine, CurrShiftBlockSumTag, ierr)
        IF (ierr /= 0) CALL Stop_All(this_routine, 'Problem allocating CurrShiftBlockSum.')
        CurrShiftBlockSum(:) = 0.0_dp

        allocate(ShiftBlockSum(0:TotNoShiftBlockSizes), stat=ierr)
        CALL LogMemAlloc('ShiftBlockSum', TotNoShiftBlockSizes, 8, this_routine, ShiftBlockSumTag, ierr)
        IF (ierr /= 0) CALL Stop_All(this_routine, 'Problem allocating ShiftBlockSum.')
        ShiftBlockSum(:) = 0.0_dp

        allocate(ShiftBlockSqrdSum(0:TotNoShiftBlockSizes), stat=ierr)
        CALL LogMemAlloc('ShiftBlockSqrdSum', TotNoShiftBlockSizes, 8, this_routine, ShiftBlockSqrdSumTag, ierr)
        IF (ierr /= 0) CALL Stop_All(this_routine, 'Problem allocating ShiftBlockSqrdSum.')
        ShiftBlockSqrdSum(:) = 0.0_dp

    END SUBROUTINE InitShiftErrorBlocking

    SUBROUTINE RestartBlocking(Iter)
        INTEGER :: Iter

        StartBlockIter = 0
        StartBlockIter = Iter + StepsSft
        ! ChangeVars gets called at the end of the run, wont actually start until the next iteration.

        TotNoBlockSizes = 0
        IF (tBlockEveryIteration) THEN
            TotNoBlockSizes = FLOOR((LOG10(REAL(NMCyc - StartBlockIter, dp))) / (LOG10(2.0_dp)))
        ELSE
            TotNoBlockSizes = FLOOR((LOG10((REAL(NMCyc - StartBlockIter, dp)) / (REAL(StepsSft, dp)))) / (LOG10(2.0_dp)))
        end if

        if (allocated(CurrBlockSum)) CurrBlockSum(:) = 0.0_dp
        if (allocated(BlockSum)) BlockSum(:) = 0.0_dp
        if (allocated(BlockSqrdSum)) BlockSqrdSum(:) = 0.0_dp

    END SUBROUTINE RestartBlocking

    SUBROUTINE RestartShiftBlocking(Iter)
        INTEGER :: Iter

        StartShiftBlockIter = 0
        StartShiftBlockIter = Iter + StepsSft
        ! ChangeVars gets called at the end of the run, wont actually start until the next iteration.

        TotNoShiftBlockSizes = 0
        TotNoShiftBlockSizes = FLOOR((LOG10((REAL(NMCyc - StartShiftBlockIter, dp)) / (REAL(StepsSft, dp)))) / (LOG10(2.0_dp)))

        if (allocated(CurrShiftBlockSum)) CurrShiftBlockSum(:) = 0.0_dp
        if (allocated(ShiftBlockSum)) ShiftBlockSum(:) = 0.0_dp
        if (allocated(ShiftBlockSqrdSum)) ShiftBlockSqrdSum(:) = 0.0_dp

    END SUBROUTINE RestartShiftBlocking

    SUBROUTINE SumInErrorContrib(Iter, AllENumCyc, AllHFCyc)
        INTEGER :: i, Iter, NoContrib
        real(dp) :: AllENumCyc, AllHFCyc

! First we need to find out what number contribution to the blocking this is.
        IF (tBlockEveryIteration) THEN
            NoContrib = (Iter - StartBlockIter) + 1
        ELSE
            NoContrib = ((Iter - StartBlockIter) / StepsSft) + 1
        end if

! We then need to add the contribution from this cycle to all the relevant blocks.

        do i = 0, TotNoBlockSizes

! First of all sum the energy contributions into each of the blocks.
            CurrBlockSum(i) = CurrBlockSum(i) + (AllENumCyc / AllHFCyc)
            ! The values contained in this array is the sum of the energy contributions to that block.
            ! They have not been divided by the number of contributions yet.

            ! Find out if we have completed a block.
            ! If we have, then add the average value from that block to the overall sum array (BlockSum),
            ! and zero the CurrBlockSum.
            IF (MOD(NoContrib, (2**i)) == 0) THEN
                BlockSum(i) = BlockSum(i) + (CurrBlockSum(i) / (2**i))
                BlockSqrdSum(i) = BlockSqrdSum(i) + ((CurrBlockSum(i) / (2**i))**2)
                CurrBlockSum(i) = 0.0_dp
            end if

        end do

    END SUBROUTINE SumInErrorContrib

    SUBROUTINE SumInShiftErrorContrib(Iter, DiagSft)
        INTEGER :: i, Iter, NoContrib
        real(dp) :: DiagSft

! First we need to find out what number contribution to the blocking this is.
        NoContrib = ((Iter - StartShiftBlockIter) / StepsSft) + 1

! We then need to add the contribution from this cycle to all the relevant blocks.

        do i = 0, TotNoShiftBlockSizes

! First of all sum the energy contributions into each of the blocks.
            CurrShiftBlockSum(i) = CurrShiftBlockSum(i) + DiagSft
            ! The values contained in this array is the sum of the energy contributions to that block.
            ! They have not been divided by the number of contributions yet.

            ! Find out if we have completed a block.
            ! If we have, then add the average value from that block to the overall sum array (BlockSum),
            ! and zero the CurrBlockSum.
            IF (MOD(NoContrib, (2**i)) == 0) THEN
                ShiftBlockSum(i) = ShiftBlockSum(i) + (CurrShiftBlockSum(i) / (2**i))
                ShiftBlockSqrdSum(i) = ShiftBlockSqrdSum(i) + ((CurrShiftBlockSum(i) / (2**i))**2)
                CurrShiftBlockSum(i) = 0.0_dp
            end if

        end do

    END SUBROUTINE SumInShiftErrorContrib

    SUBROUTINE PrintBlocking(Iter)
        use util_mod, only: get_free_unit
        INTEGER :: i, NoBlocks, Iter, NoBlockSizes, NoContrib, iunit
        real(dp) :: MeanEn, MeanEnSqrd, StandardDev, Error, ErrorinError
        CHARACTER(len=30) :: abstr

!First find out how many blocks would have been formed with the number of iterations actually performed.

        NoContrib = 0
        IF (tBlockEveryIteration) THEN
            NoContrib = (Iter - StartBlockIter) + 1
        ELSE
            NoContrib = ((Iter - StartBlockIter) / StepsSft) + 1
        end if
        IF (NoContrib == 1) RETURN

        NoBlockSizes = 0
        NoBlockSizes = FLOOR((LOG10(REAL(NoContrib - 1, dp))) / (LOG10(2.0_dp)))

        iunit = get_free_unit()
        IF (tSaveBlocking) THEN
!We are saving the blocking files - write to new file.
            abstr = 'BLOCKINGANALYSIS-'//str(Iter)
            open(iunit, file=abstr, status='unknown')
        ELSE
            open(iunit, file='BLOCKINGANALYSIS', status='unknown')
        end if

        write(iunit, '(I4,A,I4)') NoBlockSizes, ' blocks were formed with sizes from 1 to ', (2**(NoBlockSizes))
        write(iunit, '(3A16,5A20)') '1.Block No.', '2.Block Size  ', '3.No. of Blocks', '4.Mean E', &
            '5.Mean E^2', '6.SD', '7.Error', '8.ErrorinError'

        do i = 0, NoBlockSizes - 1

            ! First need to find out how many blocks of this particular size contributed to the final sum in BlockSum.
            ! NoContrib is the total number of contributions to the blocking throughout the simulation.
            NoBlocks = 0
            NoBlocks = FLOOR(REAL(NoContrib / (2**i), dp))
            ! This finds the lowest integer multiple of 2**i (the block size).

            MeanEn = BlockSum(i) / REAL(NoBlocks, dp)
            MeanEnSqrd = BlockSqrdSum(i) / REAL(NoBlocks, dp)

            StandardDev = SQRT(MeanEnSqrd - (MeanEn**2))
            IF (abs(StandardDev) < 1.0e-12_dp) THEN
                Error = 0.0_dp
                ErrorinError = 0.0_dp
            ELSE
                Error = StandardDev / SQRT(REAL(NoBlocks - 1, dp))
                ErrorinError = Error / SQRT(2.0_dp * (REAL(NoBlocks - 1, dp)))
            end if

!This is from the blocking paper, and indicates the error in the blocking error, due to the limited number of blocks available.

            write(iunit, '(3I16,5F20.10)') i, (2**i), NoBlocks, MeanEn, MeanEnSqrd, StandardDev, Error, ErrorinError

        end do

        close(iunit)

    END SUBROUTINE PrintBlocking

    SUBROUTINE PrintShiftBlocking(Iter)
        use util_mod, only: get_free_unit
        INTEGER :: i, NoBlocks, Iter, NoBlockSizes, NoContrib, iunit
        real(dp) :: MeanShift, MeanShiftSqrd, StandardDev, Error, ErrorinError

!First find out how many blocks would have been formed with the number of iterations actually performed.

        NoContrib = 0
        NoContrib = ((Iter - StartShiftBlockIter) / StepsSft) + 1

        NoBlockSizes = 0
        NoBlockSizes = FLOOR((LOG10(REAL(NoContrib - 1, dp))) / (LOG10(2.0_dp)))

        iunit = get_free_unit()
        open(iunit, file='SHIFTBLOCKINGANALYSIS', status='unknown')
        write(iunit, '(I4,A,I4)') NoBlockSizes, ' blocks were formed with sizes from 1 to ', (2**(NoBlockSizes))
        write(iunit, '(3A16,5A20)') '1.Block No.', '2.Block Size  ', '3.No. of Blocks', '4.Mean Shift', &
            '5.Mean Shift^2', '6.SD', '7.Error', '8.ErrorinError'

        do i = 0, NoBlockSizes - 1

            ! First need to find out how many blocks of this particular size contributed to the final sum in BlockSum.
            ! NoContrib is the total number of contributions to the blocking throughout the simulation.

            NoBlocks = 0
            NoBlocks = FLOOR(REAL(NoContrib / (2**i), dp))
            ! This finds the lowest integer multiple of 2**i (the block size).

            MeanShift = ShiftBlockSum(i) / REAL(NoBlocks, dp)
            MeanShiftSqrd = ShiftBlockSqrdSum(i) / REAL(NoBlocks, dp)

            StandardDev = SQRT(MeanShiftSqrd - (MeanShift**2))
            IF (abs(StandardDev) < 1.0e-12_dp) THEN
                Error = 0.0_dp
                ErrorinError = 0.0_dp
            ELSE
                Error = StandardDev / SQRT(REAL(NoBlocks - 1, dp))
                ErrorinError = Error / SQRT(2.0_dp * (REAL(NoBlocks - 1, dp)))
            end if

!This is from the blocking paper, and indicates the error in the blocking error, due to the limited number of blocks available.

            write(iunit, '(3I16,5F20.10)') i, (2**i), NoBlocks, MeanShift, MeanShiftSqrd, StandardDev, Error, ErrorinError

        end do

        close(iunit)

    END SUBROUTINE PrintShiftBlocking

    SUBROUTINE FinaliseBlocking(Iter)
        INTEGER :: Iter
        CHARACTER(len=*), PARAMETER :: this_routine = 'FinaliseBlocking'

        CALL PrintBlocking(Iter)

        DEallocate(CurrBlockSum)
        CALL LogMemDeAlloc(this_routine, CurrBlockSumTag)

        DEallocate(BlockSum)
        CALL LogMemDeAlloc(this_routine, BlockSumTag)

        DEallocate(BlockSqrdSum)
        CALL LogMemDeAlloc(this_routine, BlockSqrdSumTag)

    END SUBROUTINE FinaliseBlocking

    SUBROUTINE FinaliseShiftBlocking(Iter)
        INTEGER :: Iter
        CHARACTER(len=*), PARAMETER :: this_routine = 'FinaliseShiftBlocking'

        CALL PrintShiftBlocking(Iter)

        DEallocate(CurrShiftBlockSum)
        CALL LogMemDeAlloc(this_routine, CurrShiftBlockSumTag)

        DEallocate(ShiftBlockSum)
        CALL LogMemDeAlloc(this_routine, ShiftBlockSumTag)

        DEallocate(ShiftBlockSqrdSum)
        CALL LogMemDeAlloc(this_routine, ShiftBlockSqrdSumTag)

    END SUBROUTINE FinaliseShiftBlocking

    SUBROUTINE InitHistInitPops()
        USE CalcData, only: InitiatorWalkNo
        INTEGER :: ierr
        character(*), parameter :: this_routine = 'InitHistInitPops'

        if (allocated(HistInitPops)) then
            deallocate(HistInitPops, stat=ierr)
            call LogMemDealloc(this_routine, HistInitPopsTag, ierr)
        end if

        allocate(HistInitPops(2, 25000), stat=ierr)
        call LogMemAlloc('HistInitPops', 50000, 4, this_routine, &
                         HistInitPopsTag, ierr)
        HistInitPops = 0

        if (iProcIndex == 0) then
            if (allocated(AllHistInitPops)) then
                deallocate(AllHistInitPops, stat=ierr)
                call LogMemDealloc(this_routine, AllHistInitPopsTag, ierr)
            end if

            allocate(AllHistInitPops(2, 25000), stat=ierr)
            CALL LogMemAlloc('AllHistInitPops', 50000, 4, this_routine, &
                             AllHistInitPopsTag, ierr)
            AllHistInitPops = 0
#ifdef DEBUG_
        else
            ! in debug mode, we have to allocate this on all procs
            if (allocated(AllHistInitPops)) then
                deallocate(AllHistInitPops, stat=ierr)
            end if

            allocate(AllHistInitPops(2, 1), stat=ierr)
            AllHistInitPops = 0
#endif
        end if

        InitBinMin = log(REAL(InitiatorWalkNo + 1, dp))
        InitBinMax = log(1000000.0_dp)
        InitBinIter = ABS(InitBinMax - InitBinMin) / 25000.0

    END SUBROUTINE InitHistInitPops

    SUBROUTINE WriteInitPops(Iter)
        use util_mod, only: get_free_unit
        CHARACTER(len=21) :: abstr
        INTEGER :: i, Iter, iunit
        real(dp) :: InitBinCurr

!This will open a file called InitPops-"Iter" on unit number 17.
        abstr = 'InitPops-'//str(Iter)

        call MPIReduce(HistInitPops, MPI_SUM, AllHistInitPops)

        IF (iProcIndex == 0) THEN
            iunit = get_free_unit()
            open(iunit, FILE=abstr, STATUS='unknown')

            InitBinCurr = (-1) * InitBinMax
            do i = 25000, 1, -1
                IF (AllHistInitPops(1, i) /= 0) write(iunit, '(F20.10,2I20)') &
                    InitBinCurr, (-1) * (NINT(EXP(ABS(InitBinCurr)))), AllHistInitPops(1, i)
                InitBinCurr = InitBinCurr + InitBinIter
            end do

            InitBinCurr = InitBinMin
            do i = 1, 25000
                IF (AllHistInitPops(2, i) /= 0) write(iunit, '(F20.10,2I20)') &
                    InitBinCurr, NINT(EXP(InitBinCurr)), AllHistInitPops(2, i)
                InitBinCurr = InitBinCurr + InitBinIter
            end do

            close(iunit)
            AllHistInitPops(:, :) = 0
        end if
        HistInitPops(:, :) = 0

    END SUBROUTINE WriteInitPops

    SUBROUTINE TrackSpawnAttempts(Child, DetCurr, nJ, IC, Ex, tParity)
        INTEGER :: Child, DetCurr(NEl), nJ(NEl), IC, Ex(2, maxExcit)
        LOGICAL :: tParity
        HElement_t(dp) :: HEl

!        write(stdout,*) 'Child',Child
!        write(stdout,*) 'DetCurr',DetCurr
!        write(stdout,*) 'nJ',nJ
!        write(stdout,*) 'iLutnJ',iLutnJ
!        CALL neci_flush(6)
!        stop

        ! Need to find the H element between the current determinant and that which we're trying to spawn on.
        if (tGUGA) then
            call stop_all("TrackSpawnAttempts", "modify for GUGA!")
        end if
        HEl = get_helement_excit(DetCurr, nJ, IC, Ex, tParity)

        IF (Child == 0) THEN
            ! Spawn not accepted.
            NoNotAccept = NoNotAccept + 1.0_dp
            TotHElNotAccept = TotHElNotAccept + ABS(REAL(HEl, dp))
            IF (ABS(REAL(HEl, dp)) > ABS(MaxHElNotAccept)) MaxHElNotAccept = ABS(REAL(HEl, dp))
        ELSE
            ! Spawn accepted.
            NoAccept = NoAccept + 1.0_dp
            TotHElAccept = TotHElAccept + ABS(REAL(HEl, dp))
            IF (abs(MinHElAccept) < 1.0e-12_dp .or. (ABS(REAL(HEl, dp)) < ABS(MinHElAccept))) MinHElAccept = ABS(REAL(HEl, dp))
        end if

    ENDSUBROUTINE TrackSpawnAttempts

    SUBROUTINE PrintSpawnAttemptStats(Iteration)
        use util_mod, only: get_free_unit
        real(dp) :: AllStats(4), AcceptStats(4), AllMaxHElNotAccept(1:nProcessors), AllMinHElAccept(1:nProcessors)
        INTEGER :: i, error, Iteration, iunit

        ! Need to distribute the max and min values to all processors, but only if it has changed.
        AcceptStats(1) = TotHElNotAccept       ! Total not accepted
        AcceptStats(2) = TotHElAccept          ! Total accepted
        AcceptStats(3) = NoNotAccept
        AcceptStats(4) = NoAccept
        AllStats(:) = 0.0_dp
        AllMaxHElNotAccept(:) = 0.0_dp
        AllMinHElAccept(:) = 0.0_dp

        call MPIReduce(AcceptStats, MPI_Sum, AllStats)
        call MPIGather(MaxHElNotAccept, AllMaxHElNotAccept(1:nProcessors), error)
        call MPIGather(MinHElAccept, AllMinHElAccept(1:nProcessors), error)

        IF (iProcIndex == Root) THEN
!            write(stdout,*) 'AllMinHElAccept',AllMinHElAccept
            CALL neci_flush(6)
            MaxHElNotAccept = ABS(AllMaxHElNotAccept(1))
            do i = 2, nProcessors
                IF (ABS(AllMaxHElNotAccept(i)) > ABS(MaxHElNotAccept)) MaxHElNotAccept = ABS(AllMaxHElNotAccept(i))
            end do

            MinHElAccept = 0.0_dp
            IF (AllStats(4) > 0.0_dp) THEN
                do i = 1, nProcessors
                    IF (abs(AllMinHElAccept(i)) > 1.0e-12_dp) THEN
                        MinHElAccept = ABS(AllMinHElAccept(i))
                        EXIT
                    end if
                end do
                do i = 1, nProcessors
                    IF (abs(AllMinHElAccept(i)) > 1.0e-12_dp .and. (ABS(AllMinHElAccept(i)) < ABS(MinHElAccept))) &
                        MinHElAccept = ABS(AllMinHElAccept(i))
                end do
            end if

            iunit = get_free_unit()
            open(iunit, file="SpawnAttemptStats", status="unknown")
            write(iunit, '(I10,2F20.1,5F20.6)') Iteration, AllStats(3), AllStats(4), &
                AllStats(3) / AllStats(4), AllStats(1) / AllStats(3), AllStats(2) / AllStats(4), MaxHElNotAccept, MinHElAccept
            close(iunit)
        end if

    ENDSUBROUTINE PrintSpawnAttemptStats

    subroutine HistInitPopulations(SignCurr, VecSlot)

        integer, intent(in) :: VecSlot
        real(dp), intent(in) :: SignCurr
        integer :: InitBinNo

        if (abs(SignCurr) > InitiatorWalkNo) then
            ! Just summing in those determinants which are initiators.
            ! Need to figure out which bin to put them in though.
            InitBinNo = floor((log(abs(SignCurr)) - InitBinMin) / &
                              InitBinIter) + 1
            if (InitBinNo >= 1 .and. InitBinNo <= 25000) then
                if (SignCurr < 0) then
                    HistInitPops(1, InitBinNo) = HistInitPops(1, InitBinNo) + 1
                else
                    HistInitPops(2, InitBinNo) = HistInitPops(2, InitBinNo) + 1
                end if
!            else
!                call stop_all (this_routine, 'Trying to histogram outside&
!                              & the range of the bins.')
            end if
        end if

        if (SignCurr < MaxInitPopNeg) then
            MaxInitPopNeg = SignCurr
            HighPopNeg = VecSlot
        end if
        if (SignCurr > MaxInitPopPos) then
            MaxInitPopPos = SignCurr
            HighPopPos = VecSlot
        end if

    end subroutine HistInitPopulations

ENDMODULE FciMCLoggingMod

