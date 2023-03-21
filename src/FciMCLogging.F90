!This is a parallel MPI version of the FciMC code.
!All variables refer to values per processor

MODULE FciMCLoggingMod

    USE Global_utilities
    USE Parallel_neci
    USE LoggingData, only: tSaveBlocking, tBlockEveryIteration, HistInitPops, HistInitPopsTag, AllHistInitPops, AllHistInitPopsTag
    use SystemData, only: NEl, tGUGA
    use bit_rep_data, only: NIfTot
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

