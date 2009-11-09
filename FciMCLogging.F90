!This is a parallel MPI version of the FciMC code.
!All variables refer to values per processor

MODULE FciMCLoggingMod

    USE Global_utilities
    USE Parallel
    USE HElem , only : HElement
    USE Logging , only : tPrintTriConnections,TriConMax,NoTriConBins,tHistTriConHEls,NoTriConHElBins,TriConHElSingMax,TriConHElDoubMax
    USE Logging , only : tPrintHElAccept
    USE SystemData , only : NEl,NIfD,NIfTot
    USE SymData , only : nSymLabels
    USE Determinants , only : GetHElement3,GetHElement4
    use GenRandSymExcitNUMod , only : GenRandSymExcitScratchNU,ScratchSize
    USE CalcData , only : NMCyc,StepsSft
    use DetBitOps, only: DetBitEQ

    IMPLICIT NONE
    save

    INTEGER , PARAMETER :: Root=0   !This is the rank of the root processor
    INTEGER , PARAMETER :: r2=kind(0.d0)
    REAL*8 , ALLOCATABLE :: SignCohTriHist(:,:),SignIncohTriHist(:,:),SignCohHFTriHist(:,:),SignIncohHFTriHist(:,:),TriConnHElHistSing(:,:),TriConnHElHistDoub(:,:)
    REAL*8 , ALLOCATABLE :: AllSignCohTriHist(:,:),AllSignIncohTriHist(:,:),AllSignCohHFTriHist(:,:),AllSignIncohHFTriHist(:,:)
    REAL*8 , ALLOCATABLE :: AllTriConnHElHistSing(:,:),AllTriConnHElHistDoub(:,:),TriHjkHistSing(:,:),TriHjkHistDoub(:,:),AllTriHjkHistSing(:,:),AllTriHjkHistDoub(:,:)
    REAL*8 :: NoSignCohTri,NoSignInCohTri,SignCohTri,SignInCohTri,TriConHEls(3,2) 
    INTEGER :: SignCohTriHistTag,SignIncohTriHistTag,SignCohHFTriHistTag,SignIncohHFTriHistTag,TriConnHElHistSingTag,TriConnHElHistDoubTag
    INTEGER :: AllSignCohTriHistTag,AllSignIncohTriHistTag,AllSignCohHFTriHistTag,AllSignIncohHFTriHistTag
    INTEGER :: AllTriConnHElHistSingTag,AllTriConnHElHistDoubTag,TriHjkHistSingTag,TriHjkHistDoubTag,AllTriHjkHistSingTag,AllTriHjkHistDoubTag
    REAL*8 :: NoNotAccept,NoAccept,TotHElNotAccept,TotHElAccept,MaxHElNotAccept,MinHElAccept
    REAL*8 :: NoPosSpinCoup,NoNegSpinCoup,SumPosSpinCoup,SumNegSpinCoup,SumHFCon,SumSpinCon

    REAL*8 , ALLOCATABLE :: CurrBlockSum(:),BlockSum(:),BlockSqrdSum(:)
    REAL*8 , ALLOCATABLE :: CurrShiftBlockSum(:),ShiftBlockSum(:),ShiftBlockSqrdSum(:)
    INTEGER :: CurrBlockSumTag,BlockSumTag,BlockSqrdSumTag,TotNoBlockSizes,StartBlockIter
    INTEGER :: CurrShiftBlockSumTag,ShiftBlockSumTag,ShiftBlockSqrdSumTag,TotNoShiftBlockSizes,StartShiftBlockIter



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

        TotNoBlockSizes=FLOOR( (LOG10((REAL(NMCyc-StartBlockIter))/(REAL(StepsSft)))) / (LOG10(2.D0)) )
        WRITE(6,*) 'Beginning blocking analysis of the errors in the projected energies.'
        WRITE(6,"(A,I6)") "The total number of different block sizes possible is: ",TotNoBlockSizes
        ! The blocks will have size 1,2,4,8,....,2**TotNoBlockSizes
        ! In the below arrays, the element i will correspond to block size 2**(i), but the arrays go from 0 -> TotNoBlockSizes.

! Then need to allocate three arrays of this size - one for the current block, one for the sum of the blocks over the elapsed cycles,
! and one for the sum of the squares of the blocks over the elapsed iterations.

        ALLOCATE(CurrBlockSum(0:TotNoBlockSizes),stat=ierr)
        CALL LogMemAlloc('CurrBlockSum',TotNoBlockSizes,8,this_routine,CurrBlockSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating CurrBlockSum.')
        CurrBlockSum(:)=0.D0
 
        ALLOCATE(BlockSum(0:TotNoBlockSizes),stat=ierr)
        CALL LogMemAlloc('BlockSum',TotNoBlockSizes,8,this_routine,BlockSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating BlockSum.')
        BlockSum(:)=0.D0
 
        ALLOCATE(BlockSqrdSum(0:TotNoBlockSizes),stat=ierr)
        CALL LogMemAlloc('BlockSqrdSum',TotNoBlockSizes,8,this_routine,BlockSqrdSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating BlockSqrdSum.')
        BlockSqrdSum(:)=0.D0

    END SUBROUTINE InitErrorBlocking




    SUBROUTINE InitShiftErrorBlocking(Iter)
        CHARACTER(len=*), PARAMETER :: this_routine='InitShiftErrorBlocking'
        INTEGER :: ierr,Iter


        StartShiftBlockIter=Iter 
        ! This is the iteration the blocking was started.

        TotNoShiftBlockSizes=FLOOR( (LOG10((REAL(NMCyc-StartShiftBlockIter))/(REAL(StepsSft)))) / (LOG10(2.D0)) )
        WRITE(6,*) 'Beginning blocking analysis of the errors in the shift.'
        WRITE(6,"(A,I6)") "The total number of different block sizes possible is: ",TotNoShiftBlockSizes
        ! The blocks will have size 1,2,4,8,....,2**TotNoBlockSizes
        ! In the below arrays, the element i will correspond to block size 2**(i), but the arrays go from 0 -> TotNoBlockSizes.

! Then need to allocate three arrays of this size - one for the current block, one for the sum of the blocks over the elapsed cycles,
! and one for the sum of the squares of the blocks over the elapsed iterations.

        ALLOCATE(CurrShiftBlockSum(0:TotNoShiftBlockSizes),stat=ierr)
        CALL LogMemAlloc('CurrShiftBlockSum',TotNoShiftBlockSizes,8,this_routine,CurrShiftBlockSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating CurrShiftBlockSum.')
        CurrShiftBlockSum(:)=0.D0
 
        ALLOCATE(ShiftBlockSum(0:TotNoShiftBlockSizes),stat=ierr)
        CALL LogMemAlloc('ShiftBlockSum',TotNoShiftBlockSizes,8,this_routine,ShiftBlockSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating ShiftBlockSum.')
        ShiftBlockSum(:)=0.D0
 
        ALLOCATE(ShiftBlockSqrdSum(0:TotNoShiftBlockSizes),stat=ierr)
        CALL LogMemAlloc('ShiftBlockSqrdSum',TotNoShiftBlockSizes,8,this_routine,ShiftBlockSqrdSumTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating ShiftBlockSqrdSum.')
        ShiftBlockSqrdSum(:)=0.D0



    END SUBROUTINE InitShiftErrorBlocking

 

    SUBROUTINE RestartBlocking(Iter)
        INTEGER :: Iter

        StartBlockIter=0
        StartBlockIter=Iter+StepsSft
        ! ChangeVars gets called at the end of the run, wont actually start until the next iteration.

        TotNoBlockSizes=0
        TotNoBlockSizes=FLOOR( (LOG10((REAL(NMCyc-StartBlockIter))/(REAL(StepsSft)))) / (LOG10(2.D0)) )

        CurrBlockSum(:)=0.D0
        BlockSum(:)=0.D0
        BlockSqrdSum(:)=0.D0

    END SUBROUTINE RestartBlocking
 


    SUBROUTINE RestartShiftBlocking(Iter)
        INTEGER :: Iter

        StartShiftBlockIter=0
        StartShiftBlockIter=Iter+StepsSft
        ! ChangeVars gets called at the end of the run, wont actually start until the next iteration.

        TotNoShiftBlockSizes=0
        TotNoShiftBlockSizes=FLOOR( (LOG10((REAL(NMCyc-StartShiftBlockIter))/(REAL(StepsSft)))) / (LOG10(2.D0)) )

        CurrShiftBlockSum(:)=0.D0
        ShiftBlockSum(:)=0.D0
        ShiftBlockSqrdSum(:)=0.D0

    END SUBROUTINE RestartShiftBlocking



    SUBROUTINE SumInErrorContrib(Iter,AllENumCyc,AllHFCyc)
        INTEGER :: i,Iter,NoContrib
        REAL*8 :: AllENumCyc,AllHFCyc

! First we need to find out what number contribution to the blocking this is.
        NoContrib=((Iter-StartBlockIter)/StepsSft)+1

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
                CurrBlockSum(i)=0.D0
            ENDIF

        enddo

    END SUBROUTINE SumInErrorContrib



    SUBROUTINE SumInShiftErrorContrib(Iter,DiagSft)
        INTEGER :: i,Iter,NoContrib
        REAL*8 :: DiagSft

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
                CurrShiftBlockSum(i)=0.D0
            ENDIF

        enddo

    END SUBROUTINE SumInShiftErrorContrib



    SUBROUTINE PrintBlocking(Iter)
        INTEGER :: i,NoBlocks,Iter,NoBlockSizes,NoContrib
        REAL*8 :: MeanEn,MeanEnSqrd,StandardDev,Error,ErrorinError

!First find out how many blocks would have been formed with the number of iterations actually performed. 

        NoContrib=0
        NoContrib=((Iter-StartBlockIter)/StepsSft)+1

        NoBlockSizes=0
        NoBlockSizes=FLOOR( (LOG10(REAL(NoContrib-1)))/ (LOG10(2.D0)))

        OPEN(62,file='BLOCKINGANALYSIS',status='unknown')
        WRITE(62,'(I4,A,I4)') NoBlockSizes,' blocks were formed with sizes from 1 to ',(2**(NoBlockSizes))
        WRITE(62,'(3A16,5A20)') '1.Block No.','2.Block Size  ','3.No. of Blocks','4.Mean E','5.Mean E^2','6.SD','7.Error','8.ErrorinError'

        do i=0,NoBlockSizes-1
            
            ! First need to find out how many blocks of this particular size contributed to the final sum in BlockSum.
            ! NoContrib is the total number of contributions to the blocking throughout the simulation.

            NoBlocks=0
            NoBlocks=FLOOR(REAL(NoContrib/(2**i)))
            ! This finds the lowest integer multiple of 2**i (the block size). 

            MeanEn=BlockSum(i)/REAL(NoBlocks)
            MeanEnSqrd=BlockSqrdSum(i)/REAL(NoBlocks)

            StandardDev=SQRT(MeanEnSqrd-(MeanEn**2))
            IF(StandardDev.eq.0.D0) THEN
                Error=0.D0
                ErrorinError=0.D0
            ELSE
                Error=StandardDev/SQRT(REAL(NoBlocks-1))
                ErrorinError=Error/SQRT(2.D0*(REAL(NoBlocks-1)))
            ENDIF

!This is from the blocking paper, and indicates the error in the blocking error, due to the limited number of blocks available.

            WRITE(62,'(3I16,5F20.10)') i,(2**i),NoBlocks,MeanEn,MeanEnSqrd,StandardDev,Error,ErrorinError
            
        enddo

        CLOSE(62)


    END SUBROUTINE PrintBlocking


    SUBROUTINE PrintShiftBlocking(Iter)
        INTEGER :: i,NoBlocks,Iter,NoBlockSizes,NoContrib
        REAL*8 :: MeanShift,MeanShiftSqrd,StandardDev,Error,ErrorinError

!First find out how many blocks would have been formed with the number of iterations actually performed. 

        NoContrib=0
        NoContrib=((Iter-StartShiftBlockIter)/StepsSft)+1

        NoBlockSizes=0
        NoBlockSizes=FLOOR( (LOG10(REAL(NoContrib-1)))/ (LOG10(2.D0)))

        OPEN(62,file='SHIFTBLOCKINGANALYSIS',status='unknown')
        WRITE(62,'(I4,A,I4)') NoBlockSizes,' blocks were formed with sizes from 1 to ',(2**(NoBlockSizes))
        WRITE(62,'(3A16,5A20)') '1.Block No.','2.Block Size  ','3.No. of Blocks','4.Mean Shift','5.Mean Shift^2','6.SD','7.Error','8.ErrorinError'

        do i=0,NoBlockSizes-1
            
            ! First need to find out how many blocks of this particular size contributed to the final sum in BlockSum.
            ! NoContrib is the total number of contributions to the blocking throughout the simulation.

            NoBlocks=0
            NoBlocks=FLOOR(REAL(NoContrib/(2**i)))
            ! This finds the lowest integer multiple of 2**i (the block size). 

            MeanShift=ShiftBlockSum(i)/REAL(NoBlocks)
            MeanShiftSqrd=ShiftBlockSqrdSum(i)/REAL(NoBlocks)

            StandardDev=SQRT(MeanShiftSqrd-(MeanShift**2))
            IF(StandardDev.eq.0.D0) THEN
                Error=0.D0
                ErrorinError=0.D0
            ELSE
                Error=StandardDev/SQRT(REAL(NoBlocks-1))
                ErrorinError=Error/SQRT(2.D0*(REAL(NoBlocks-1)))
            ENDIF

!This is from the blocking paper, and indicates the error in the blocking error, due to the limited number of blocks available.

            WRITE(62,'(3I16,5F20.10)') i,(2**i),NoBlocks,MeanShift,MeanShiftSqrd,StandardDev,Error,ErrorinError
            
        enddo

        CLOSE(62)

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



    SUBROUTINE InitTriHElStats()
        INTEGER :: i,ierr
        REAL*8 :: BinVal,BinIter
        CHARACTER(len=*), PARAMETER :: this_routine='InitTriHElStats'

        IF(tPrintTriConnections) THEN
            NoSignCohTri=0.D0
            NoSignInCohTri=0.D0
            SignCohTri=0.D0
            SignInCohTri=0.D0
            IF(iProcIndex.eq.root) THEN
                OPEN(77,file='TriConnTotals',status='unknown')
                WRITE(77,"(A12,2A24,8A20)") "1.Iteration","2.No. Sign Coh Loops","3.No. Sign Incoh Loops","4.Sign Coh Tot","5.Sign Incoh Tot","6.Sign Coh/Iter",&
                                                &"7.Sign Incoh/Iter","8.Sign Coh/Loop","9.Sign Incoh/Loop","10.Ratio No.","11.Ratio Val."
            ENDIF

! Set up histogramms.
            ALLOCATE(SignCohTriHist(2,NoTriConBins),stat=ierr)
            CALL LogMemAlloc('SignCohTriHist',2*NoTriConBins,8,this_routine,SignCohTriHistTag,ierr)
            ALLOCATE(SignIncohTriHist(2,NoTriConBins),stat=ierr)
            CALL LogMemAlloc('SignIncohTriHist',2*NoTriConBins,8,this_routine,SignIncohTriHistTag,ierr)
 
            ALLOCATE(SignCohHFTriHist(2,NoTriConBins),stat=ierr)
            CALL LogMemAlloc('SignCohHFTriHist',2*NoTriConBins,8,this_routine,SignCohHFTriHistTag,ierr)
            ALLOCATE(SignIncohHFTriHist(2,NoTriConBins),stat=ierr)
            CALL LogMemAlloc('SignIncohHFTriHist',2*NoTriConBins,8,this_routine,SignIncohHFTriHistTag,ierr)
 
            SignCohTriHist(:,:)=0.D0
            SignIncohTriHist(:,:)=0.D0
            SignCohHFTriHist(:,:)=0.D0
            SignIncohHFTriHist(:,:)=0.D0

            IF(iProcIndex.eq.Root) THEN
                ALLOCATE(AllSignCohTriHist(2,NoTriConBins),stat=ierr)
                CALL LogMemAlloc('AllSignCohTriHist',2*NoTriConBins,8,this_routine,AllSignCohTriHistTag,ierr)
                ALLOCATE(AllSignIncohTriHist(2,NoTriConBins),stat=ierr)
                CALL LogMemAlloc('AllSignIncohTriHist',2*NoTriConBins,8,this_routine,AllSignIncohTriHistTag,ierr)
     
                ALLOCATE(AllSignCohHFTriHist(2,NoTriConBins),stat=ierr)
                CALL LogMemAlloc('AllSignCohHFTriHist',2*NoTriConBins,8,this_routine,AllSignCohHFTriHistTag,ierr)
                ALLOCATE(AllSignIncohHFTriHist(2,NoTriConBins),stat=ierr)
                CALL LogMemAlloc('AllSignIncohHFTriHist',2*NoTriConBins,8,this_routine,AllSignIncohHFTriHistTag,ierr)
     
                AllSignCohTriHist(:,:)=0.D0
                AllSignIncohTriHist(:,:)=0.D0
                AllSignCohHFTriHist(:,:)=0.D0
                AllSignIncohHFTriHist(:,:)=0.D0
            ENDIF
 
            BinIter=ABS(TriConMax)/REAL(NoTriConBins)

            BinVal=0.D0
            do i=1,NoTriConBins
                SignCohTriHist(1,i)=BinVal
                SignIncohTriHist(1,i)=(-1)*BinVal
                SignCohHFTriHist(1,i)=BinVal
                SignIncohHFTriHist(1,i)=(-1)*BinVal
                BinVal=BinVal+BinIter
            enddo
        ENDIF

        IF(tHistTriConHEls) THEN
            TriConHEls(:,:)=0.D0
            ! TriConHEls(1,1) - number of singles
            ! TriConHEls(1,2) - sum of single elements
            ! TriConHEls(2,1) - number of doubles
            ! TriConHEls(2,2) - sum of double elements
            ! TriConHEls(3,1) - number of Hjk elements
            ! TriConHEls(3,2) - sum of Hjk elements
            ALLOCATE(TriConnHElHistSing(2,NoTriConHElBins),stat=ierr)
            CALL LogMemAlloc('TriConnHElHistSing',2*NoTriConHElBins,8,this_routine,TriConnHElHistSingTag,ierr)
            ALLOCATE(TriConnHElHistDoub(2,NoTriConHElBins),stat=ierr)
            CALL LogMemAlloc('TriConnHElHistDoub',2*NoTriConHElBins,8,this_routine,TriConnHElHistDoubTag,ierr)
            ALLOCATE(TriHjkHistSing(2,NoTriConHElBins),stat=ierr)
            CALL LogMemAlloc('TriHjkHistSing',2*NoTriConHElBins,8,this_routine,TriHjkHistSingTag,ierr)
            ALLOCATE(TriHjkHistDoub(2,NoTriConHElBins),stat=ierr)
            CALL LogMemAlloc('TriHjkHistDoub',2*NoTriConHElBins,8,this_routine,TriHjkHistDoubTag,ierr)

            TriConnHElHistSing(:,:)=0.D0
            TriConnHElHistDoub(:,:)=0.D0
            TriHjkHistSing(:,:)=0.D0
            TriHjkHistDoub(:,:)=0.D0
 
            IF(iProcIndex.eq.Root) THEN
                ALLOCATE(AllTriConnHElHistSing(2,NoTriConHElBins),stat=ierr)
                CALL LogMemAlloc('AllTriConnHElHistSing',2*NoTriConHElBins,8,this_routine,AllTriConnHElHistSingTag,ierr)
                ALLOCATE(AllTriConnHElHistDoub(2,NoTriConHElBins),stat=ierr)
                CALL LogMemAlloc('AllTriConnHElHistDoub',2*NoTriConHElBins,8,this_routine,AllTriConnHElHistDoubTag,ierr)
                ALLOCATE(AllTriHjkHistSing(2,NoTriConHElBins),stat=ierr)
                CALL LogMemAlloc('AllTriHjkHistSing',2*NoTriConHElBins,8,this_routine,AllTriHjkHistSingTag,ierr)
                ALLOCATE(AllTriHjkHistDoub(2,NoTriConHElBins),stat=ierr)
                CALL LogMemAlloc('AllTriHjkHistDoub',2*NoTriConHElBins,8,this_routine,AllTriHjkHistDoubTag,ierr)

                AllTriConnHElHistSing(:,:)=0.D0
                AllTriConnHElHistDoub(:,:)=0.D0
                AllTriHjkHistSing(:,:)=0.D0
                AllTriHjkHistDoub(:,:)=0.D0
            ENDIF
     
            BinIter=ABS(2*TriConHElSingMax)/REAL(NoTriConHElBins)
            BinVal=(-1)*TriConHElSingMax
            do i=1,NoTriConHElBins
                TriConnHElHistSing(1,i)=BinVal
                TriHjkHistSing(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo
 
            BinIter=ABS(2*TriConHElDoubMax)/REAL(NoTriConHElBins)
            BinVal=(-1)*TriConHElDoubMax
            do i=1,NoTriConHElBins
                TriConnHElHistDoub(1,i)=BinVal
                TriHjkHistDoub(1,i)=BinVal
                BinVal=BinVal+BinIter
            enddo

        ENDIF


        IF(tPrintHElAccept) THEN
            NoNotAccept=0.D0
            NoAccept=0.D0
            TotHElNotAccept=0.D0
            TotHElAccept=0.D0
            MaxHElNotAccept=0.D0
            MinHElAccept=0.D0
            IF(iProcIndex.eq.root) THEN
                OPEN(84,file='HElsAcceptance',status='unknown')
                WRITE(84,'(A10,7A20)') "Iteration","No. Not Accepted","No. Accepted","Ratio NotAcc/Acc","Av.HEl Not Accept","Av.HEl Accept","Max HEl Not Accept","Min HEl Accept"
            ENDIF
        ENDIF

    ENDSUBROUTINE InitTriHElStats




    SUBROUTINE FindSpinCoupHEl(iLutHF,iLutCurr)
! Fed into here is a doubly excited occupied determinant - want to take the two excited orbitals and flip their spins.
! Then find the coupling H element between the original and spin-flipped determinants and add it to the stats.
        USE HPHFRandExcitMod , only : FindExcitBitDetSym 
        use DetBitOps, only: DecodeBitdet
        use SystemData, only: NIfTot, nel, NIfD
        INTEGER :: iLutCurr(0:NIfTot),iLutHF(0:NIfTot),i
        INTEGER :: iLutSym(0:NIfTot),nI(NEl),nJ(NEl),nHF(NEl),Ex(2,2)
        TYPE(HElement) :: SpinCoupHEl,HElHFI,HElHFJ
        LOGICAL :: DetsEqSpinCoup,tSign

! First get the spin flipped determinant.
! Can do this in two ways.  Either flip the spin of all electrons - this means that doubly occupied spat orbs will be unchanged or

        CALL FindExcitBitDetSym(iLutCurr(0:NIfD),iLutSym(0:NIfD))


! - just flip the sign of the two excited electrons.
!        CALL GetBitExcitation(iLutHF,iLutCurr,Ex,tSign)
        ! Electrons Ex(1,1) and Ex(1,2) are excited to Ex(2,1) and Ex(2,2)


! Now find the HElement between these two determinants.        
        
        CALL DecodeBitDet(nI,iLutCurr(0:NIfTot))
        CALL DecodeBitDet(nJ,iLutSym(0:NIfTot))

        CALL DecodeBitDet(nHF,iLutHF(0:NIfTot))

! Want to replace the excited electrons in nI with the spin flipped versions.
!        nJ(:)=nI(:)
!        do i=1,NEl
!            IF(nJ(i).eq.Ex(2,1)) THEN
!                IF(MOD(nJ(i),2).eq.0) THEN
!                    nJ(i)=Ex(2,1)-1
!                ELSE
!                    nJ(i)=Ex(2,1)+1
!                ENDIF
!            ELSEIF(nJ(i).eq.Ex(2,2)) THEN
!                IF(MOD(nJ(i),2).eq.0) THEN
!                    nJ(i)=Ex(2,2)-1
!                ELSE
!                    nJ(i)=Ex(2,2)+1
!                ENDIF
!            ENDIF
!        enddo


        DetsEqSpinCoup=.false.
        DetsEqSpinCoup=DetBitEQ(iLutCurr(0:NIfTot),iLutSym(0:NIfTot))

        HElHFI=GetHElement3(nHF,nI,-1)
        HElHFJ=GetHElement3(nHF,nJ,-1)

        IF(.not.DetsEqSpinCoup) THEN

            SpinCoupHEl=GetHElement3(nI,nJ,-1)

            IF((((REAL(HElHFI%v,r2)).lt.0.D0).and.((REAL(HElHFJ%v,r2)).gt.0.D0)).or.(((REAL(HElHFI%v,r2)).gt.0.D0).and.((REAL(HElHFJ%v,r2)).lt.0.D0))) THEN
!                WRITE(6,*) '*'
!                WRITE(6,'(A30,F15.6,A30,F15.6)') 'HEl between HF and one det : ',REAL(HElHFI%v,r2),' and to the spin coupled : ',REAL(HElHFJ%v,r2)
!                WRITE(6,*) 'HFDet',nHF(:)
!                WRITE(6,*) 'First Det',nI(:)
!                WRITE(6,*) 'Second Det',nJ(:)

                IF((REAL(SpinCoupHEl%v,r2)).lt.0.D0) THEN
                    NoNegSpinCoup=NoNegSpinCoup+1.D0
                    SumNegSpinCoup=SumNegSpinCoup+REAL(SpinCoupHEl%v,r2)
                ELSEIF((REAL(SpinCoupHEl%v,r2)).gt.0.D0) THEN
                    NoPosSpinCoup=NoPosSpinCoup+1.D0
                    SumPosSpinCoup=SumPosSpinCoup+REAL(SpinCoupHEl%v,r2)
                ENDIF
!                WRITE(6,*) 'Spin coupled HEl',REAL(SpinCoupHEl%v,r2)            
                SumHFCon=SumHFCon+ABS(REAL(HElHFI%v,r2))
                SumSpinCon=SumSpinCon+ABS(REAL(SpinCoupHEl%v,r2))

            ENDIF

            IF(((((REAL(HElHFI%v,r2)).lt.0.D0).and.((REAL(HElHFJ%v,r2)).lt.0.D0)).or.(((REAL(HElHFI%v,r2)).gt.0.D0).and.((REAL(HElHFJ%v,r2)).gt.0.D0)))&
            &.and.(REAL(SpinCoupHEl%v,r2).ne.0.D0)) THEN
                WRITE(6,*) '*'
                WRITE(6,'(A30,F15.6,A30,F15.6)') 'HEl between HF and one det : ',REAL(HElHFI%v,r2),' and to the spin coupled : ',REAL(HElHFJ%v,r2)
                WRITE(6,*) 'HFDet',nHF(:)
                WRITE(6,*) 'First Det',nI(:)
                WRITE(6,*) 'Second Det',nJ(:)

                WRITE(6,*) 'Spin coupled HEl',REAL(SpinCoupHEl%v,r2)            
                CALL FLUSH(6)
                WRITE(6,*) '******* Determinants have the same sign with HF, but non-zero connecting element.'
!                CALL Stop_All("FindSpinCoupHEl","Determinants have the same sign with HF, but non-zero connecting element.")
            ENDIF
        ENDIF


    ENDSUBROUTINE FindSpinCoupHEl



    SUBROUTINE PrintSpinCoupHEl(Iteration)
        REAL*8 :: SpinCoupHElStats(4),AllSpinCoupHElStats(4)
        INTEGER :: error,Iteration

!Write to files the sum of the sign coherent and incoherent triangles. 
        SpinCoupHElStats(1)=NoPosSpinCoup
        SpinCoupHElStats(2)=NoNegSpinCoup
        SpinCoupHElStats(3)=SumPosSpinCoup
        SpinCoupHElStats(4)=SumNegSpinCoup
        AllSpinCoupHElStats(:)=0.D0

#ifdef PARALLEL
        CALL MPI_Reduce(SpinCoupHElStats,AllSpinCoupHElStats,4,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
#else        
        AllSpinCoupHElStats=SpinCoupHElStats
#endif

        IF(iProcIndex.eq.Root) THEN

!            WRITE(87,'(I8,10F19.6)') Iteration,AllSpinCoupHElStats(1),AllSpinCoupHElStats(2),AllSpinCoupHElStats(3),AllSpinCoupHElStats(4),(AllSpinCoupHElStats(3)+AllSpinCoupHElStats(4)),&
            WRITE(87,'(I8,11F18.6)') Iteration,AllSpinCoupHElStats(1),AllSpinCoupHElStats(2),AllSpinCoupHElStats(3),AllSpinCoupHElStats(4),&
                                             &(AllSpinCoupHElStats(1)/REAL(Iteration)),(AllSpinCoupHElStats(2)/REAL(Iteration)),(AllSpinCoupHElStats(3)/REAL(Iteration)),&
!                                             &(AllSpinCoupHElStats(4)/REAL(Iteration)),((AllSpinCoupHElStats(3)+AllSpinCoupHElStats(4))/REAL(Iteration))
                                             &(AllSpinCoupHElStats(4)/REAL(Iteration)),SumHFCon,SumSpinCon,SumSpinCon/SumHFCon
        ENDIF

    ENDSUBROUTINE PrintSpinCoupHEl
 



    SUBROUTINE FindTriConnections(DetCurr,iLutnJ,iLutHF,nJ,IC,Ex,pDoubles,tFilled,tParity,Scratch1,Scratch2,exflag)
        TYPE(HElement) :: Hjk,Hij,Hik,HEl
        INTEGER :: iLutnJ(0:NIfTot),k,DetCurr(NEl),nJ(NEl),IC,Ex(2,2),Scratch1(ScratchSize),Scratch2(ScratchSize)
        INTEGER :: nK(NEl),IC2,IC3,Ex2(2,2),iLutnJ2(0:NIfTot),iLutnK(0:NIfTot),BinNo,NoPos,NoNeg,ICgen,iLutHF(0:NIfTot),exflag
        LOGICAL :: tParity2,DetsEqTri,tHF,tFilled,tParity
        REAL*8 :: Prob2,pDoubles
        

        CALL GenRandSymExcitScratchNU(DetCurr,iLutnJ,nK,pDoubles,IC2,Ex2,tParity2,exFlag,Prob2,Scratch1,Scratch2,tFilled)

        ! Need to check that the determinant we just generated is not the same as nJ.
        DetsEqTri=.false.

        ! These routines find the bit representation of nJ and nK given the excitation matrices Ex and Ex2 respectively.
        CALL FindExcitBitDet(iLutnJ,iLutnJ2,IC,Ex,NIfD)
        CALL FindExcitBitDet(iLutnJ,iLutnK,IC2,Ex2,NIfD)

        DetsEqTri=DetBitEQ(iLutnJ2(0:NIfTot),iLutnK(0:NIfTot))

        IF(.not.DetsEqTri) THEN
            ! Add the connecting elements to the relevant sum.

            ! First quickly test if any of the determinants are the HF.
            tHF=.false.
            tHF=DetBitEQ(iLutHF(0:NIfTot),iLutnJ(0:NIfTot))
            IF(.not.tHF) tHF=DetBitEQ(iLutHF(0:NIfTot),iLutnJ2(0:NIfTot))
            IF(.not.tHF) tHF=DetBitEQ(iLutHF(0:NIfTot),iLutnK(0:NIfTot))

            ! Calculate Hjk first (connecting element between two excitations), because if this is 0, no need to go further.
            CALL FindBitExcitLevel(iLutnJ2,iLutnK,IC3,NEl)
            Hjk=GetHElement3(nJ,nK,IC3)

            ! Histogram and add in the Hjk elements - regardless of whether or not this is 0.
            ! If the connection is not via a double or a single, the element will not be histogrammed, but it will always be 0,
            ! and this will be added into the sum.

            TriConHEls(3,1)=TriConHEls(3,1)+1.D0
            TriConHEls(3,2)=TriConHEls(3,2)+ABS(REAL(Hjk%v,r2))
            IF(IC3.eq.1) THEN
                BinNo=CEILING((REAL(Hjk%v,r2)+TriConHElSingMax)*NoTriConHElBins)/(2*TriConHElSingMax)
                TriHjkHistSing(2,BinNo)=TriHjkHistSing(2,BinNo)+1.D0
            ELSEIF(IC3.eq.2) THEN
                BinNo=CEILING((REAL(Hjk%v,r2)+TriConHElDoubMax)*NoTriConHElBins)/(2*TriConHElDoubMax)
                TriHjkHistDoub(2,BinNo)=TriHjkHistDoub(2,BinNo)+1.D0
            ENDIF 

            ! Now histogram all the stats from the whole loops.
            IF((REAL(Hjk%v,r2)).ne.0.D0) THEN
                NoPos=0
                NoNeg=0
                IF((REAL(Hjk%v,r2)).gt.0.D0) NoPos=NoPos+1
                IF((REAL(Hjk%v,r2)).lt.0.D0) NoNeg=NoNeg+1

                Hij=GetHElement4(DetCurr,nJ,IC,Ex,tParity)
                IF((REAL(Hij%v,r2)).gt.0.D0) NoPos=NoPos+1
                IF((REAL(Hij%v,r2)).lt.0.D0) NoNeg=NoNeg+1

                Hik=GetHElement4(DetCurr,nK,IC2,Ex2,tParity2)
                IF((REAL(Hik%v,r2)).gt.0.D0) NoPos=NoPos+1
                IF((REAL(Hik%v,r2)).lt.0.D0) NoNeg=NoNeg+1

                ! If there are 1 or 3 positive elements, the triangular connection is 'sign coherent'.
                ! i.e. if a walker starts with a positive sign at i, it would return to i with a positive sign after completing the 
                ! three cycle loop.
                ! If there are 0 or 2 positive elements, the walker would return with the opposite sign from its starting point, 
                ! and the loop is considered 'sign incoherent'.

                IF((NoPos.eq.1).or.(NoPos.eq.3)) THEN
                    SignCohTri=SignCohTri+(REAL(Hjk%v,r2)*REAL(Hij%v,r2)*REAL(Hik%v,r2))
                    NoSignCohTri=NoSignCohTri+1.D0
                    BinNo=CEILING(((REAL(Hjk%v,r2)*REAL(Hij%v,r2)*REAL(Hik%v,r2))*NoTriConBins)/TriConMax)
                    IF((BinNo.gt.NoTriConBins)) THEN
                        WRITE(6,*) 'The value about to be histogrammed is :',(REAL(Hjk%v,r2)*REAL(Hij%v,r2)*REAL(Hik%v,r2))
                        CALL FLUSH(6)
                        CALL Stop_All('PerformFCIMCCycle','Trying to histogram the sign coherent triangles of determinants, &
                                                                                & but a value is outside the chosen range.')
                    ELSEIF((BinNo.le.0).and.((REAL(Hjk%v,r2)*REAL(Hij%v,r2)*REAL(Hik%v,r2)).ne.0.D0)) THEN
                        WRITE(6,*) 'The value about to be histogrammed is :',(REAL(Hjk%v,r2)*REAL(Hij%v,r2)*REAL(Hik%v,r2))
                        CALL FLUSH(6)
                        CALL Stop_All('PerformFCIMCCycle','Trying to histogram the sign coherent triangles of determinants, &
                                                                                & but a value is below 0.')
                    ELSEIF(BinNo.gt.0) THEN
                        SignCohTriHist(2,BinNo)=SignCohTriHist(2,BinNo)+1.D0
                        IF(tHF) SignCohHFTriHist(2,BinNo)=SignCohHFTriHist(2,BinNo)+1.D0
                    ENDIF
                ELSEIF((NoNeg.eq.1).or.(NoNeg.eq.3)) THEN
                    SignIncohTri=SignIncohTri+(REAL(Hjk%v,r2)*REAL(Hij%v,r2)*REAL(Hik%v,r2))
                    NoSignIncohTri=NoSignIncohTri+1.D0
                    BinNo=CEILING((ABS((REAL(Hjk%v,r2)*REAL(Hij%v,r2)*REAL(Hik%v,r2)))*NoTriConBins)/TriConMax)
                    IF((BinNo.gt.NoTriConBins)) THEN
                        WRITE(6,*) 'The value about to be histogrammed is :',(REAL(Hjk%v,r2)*REAL(Hij%v,r2)*REAL(Hik%v,r2))
                        CALL FLUSH(6)
                        CALL Stop_All('PerformFCIMCCycle','Trying to histogram the sign coherent triangles of determinants, &
                                                                                & but a value is outside the chosen range.')
                    ELSEIF((BinNo.le.0).and.((REAL(Hjk%v,r2)*REAL(Hij%v,r2)*REAL(Hik%v,r2)).ne.0.D0)) THEN
                        WRITE(6,*) 'The value about to be histogrammed is :',(REAL(Hjk%v,r2)*REAL(Hij%v,r2)*REAL(Hik%v,r2))
                        CALL FLUSH(6)
                        CALL Stop_All('PerformFCIMCCycle','Trying to histogram the sign coherent triangles of determinants, &
                                                                                & but a value is below 0.')
                    ELSEIF(BinNo.gt.0) THEN 
                        SignIncohTriHist(2,BinNo)=SignIncohTriHist(2,BinNo)+1.D0
                        IF(tHF) SignIncohHFTriHist(2,BinNo)=SignIncohHFTriHist(2,BinNo)+1.D0
                    ENDIF
                ENDIF

                ! TriConHEls(1,1) - number of singles
                ! TriConHEls(1,2) - sum of single elements
                ! TriConHEls(2,1) - number of doubles
                ! TriConHEls(2,2) - sum of double elements
         
                k=1
                do while (k.le.3)
                    ! consider each of the 3 H elements, whose excitation levels have been calculated.
                    IF(k.eq.1) THEN
                        ICgen=IC3
                        HEl=Hjk
                    ELSEIF(k.eq.2) THEN
                        ICgen=IC2
                        HEl=Hik
                    ELSEIF(k.eq.3) THEN
                        ICgen=IC
                        HEl=Hij
                    ELSE
                        WRITE(6,*) 'error in k'
                        CALL FLUSH(6)
                        stop
                    ENDIF

                    ! add the H elements to the appropriate histogram, depending on their excitation level.
                    IF(ICgen.eq.1) THEN
                        TriConHEls(1,1)=TriConHEls(1,1)+1.D0
                        TriConHEls(1,2)=TriConHEls(1,2)+ABS(REAL(HEl%v,r2))
                        BinNo=CEILING((REAL(HEl%v,r2)+TriConHElSingMax)*NoTriConHElBins)/(2*TriConHElSingMax)
                        TriConnHElHistSing(2,BinNo)=TriConnHElHistSing(2,BinNo)+1.D0
                    ELSEIF(ICgen.eq.2) THEN
                        TriConHEls(2,1)=TriConHEls(2,1)+1.D0
                        TriConHEls(2,2)=TriConHEls(2,2)+ABS(REAL(HEl%v,r2))
                        BinNo=CEILING((REAL(HEl%v,r2)+TriConHElDoubMax)*NoTriConHElBins)/(2*TriConHElDoubMax)
                        TriConnHElHistDoub(2,BinNo)=TriConnHElHistDoub(2,BinNo)+1.D0
                    ELSE
                        WRITE(6,*) 'H element value : ',REAL(HEl%v,r2)
                        WRITE(6,*) 'IC (excitation level) : ',ICgen
                        CALL Stop_All('PerformFCIMCCycle','An H element is neither a single nor double, but it is supposedly &
                                                           & connected.')
                    ENDIF

                    IF(BinNo.gt.NoTriConHElBins) THEN
                        WRITE(6,*) 'The value about to be histogrammed is :',(REAL(HEl%v,r2))
                        WRITE(6,*) 'With excitation level : ',ICgen
                        CALL FLUSH(6)
                        CALL Stop_All('PerformFCIMCCycle','Trying to histogram an H element in a triangle of determinants, &
                                                                                & but the value is outside the chosen range.')
                    ENDIF
                    IF((BinNo.le.0).and.(REAL(HEl%v,r2).ne.0.D0)) THEN
                        WRITE(6,*) 'The value about to be histogrammed is :',(REAL(HEl%v,r2))
                        WRITE(6,*) 'With excitation level : ',ICgen
                        WRITE(6,*) 'Bin number : ',BinNo
                        CALL FLUSH(6)
                        CALL Stop_All('PerformFCIMCCycle','Trying to histogram an H element in a triangle of determinants, &
                                                                                & but the value is below 0.')
                    ENDIF

                    k=k+1
                enddo

            ENDIF
        ENDIF

    ENDSUBROUTINE FindTriConnections


    SUBROUTINE TrackSpawnAttempts(Child,DetCurr,j,nJ,iLutnJ,IC,Ex,tParity)
        INTEGER :: Child,DetCurr(NEl),j,nJ(NEl),iLutnJ(0:NIfTot),IC,Ex(2,2)
        LOGICAL :: tParity
        TYPE(HElement) :: HEl

!        WRITE(6,*) 'Child',Child
!        WRITE(6,*) 'DetCurr',DetCurr
!        WRITE(6,*) 'nJ',nJ
!        WRITE(6,*) 'iLutnJ',iLutnJ
!        CALL FLUSH(6)
!        stop

        ! Need to find the H element between the current determinant and that which we're trying to spawn on.
        HEl=GetHElement4(DetCurr,nJ,IC,Ex,tParity)
            
        IF(Child.eq.0) THEN
            ! Spawn not accepted.
            NoNotAccept=NoNotAccept+1.D0
            TotHElNotAccept=TotHElNotAccept+ABS(REAL(HEl%v,r2))
            IF(ABS(REAL(HEl%v,r2)).gt.ABS(MaxHElNotAccept)) MaxHElNotAccept=ABS(REAL(HEl%v,r2))
        ELSE
            ! Spawn accepted.
            NoAccept=NoAccept+1.D0
            TotHElAccept=TotHElAccept+ABS(REAL(HEl%v,r2))
            IF((MinHElAccept.eq.0.D0).or.(ABS(REAL(HEl%v,r2)).lt.ABS(MinHElAccept))) MinHElAccept=ABS(REAL(HEl%v,r2))
        ENDIF

    ENDSUBROUTINE TrackSpawnAttempts
        

    SUBROUTINE PrintSpawnAttemptStats(Iteration)
        REAL*8 :: AllStats(4),AcceptStats(4),AllMaxHElNotAccept(1:nProcessors),AllMinHElAccept(1:nProcessors)
        INTEGER :: i,error,Iteration

        ! Need to distribute the max and min values to all processors, but only if it has changed.        
        AcceptStats(1)=TotHElNotAccept       ! Total not accepted
        AcceptStats(2)=TotHElAccept          ! Total accepted
        AcceptStats(3)=NoNotAccept
        AcceptStats(4)=NoAccept
        AllStats(:)=0.D0
        AllMaxHElNotAccept(:)=0.D0
        AllMinHElAccept(:)=0.D0

!        WRITE(6,*) 'MinHElAccept',MinHElAccept
        CALL FLUSH(6)
!        CALL MPI_Barrier(MPI_COMM_WORLD,error)

#ifdef PARALLEL
        CALL MPI_Reduce(AcceptStats,AllStats,4,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Gather(MaxHElNotAccept,1,MPI_DOUBLE_PRECISION,AllMaxHElNotAccept(1:nProcessors),1,MPI_DOUBLE_PRECISION,Root,MPI_COMM_WORLD,error)
        CALL MPI_Gather(MinHElAccept,1,MPI_DOUBLE_PRECISION,AllMinHElAccept(1:nProcessors),1,MPI_DOUBLE_PRECISION,Root,MPI_COMM_WORLD,error)
#else
        AllStats=AcceptStats
        AllMaxHElNotAccept=MaxHElNotAccept
        AllMinHElAccept=MinHElAccept
#endif        


        IF(iProcIndex.eq.Root) THEN 
!            WRITE(6,*) 'AllMinHElAccept',AllMinHElAccept
            CALL FLUSH(6)
            MaxHElNotAccept=ABS(AllMaxHElNotAccept(1))
            do i=2,nProcessors
                IF(ABS(AllMaxHElNotAccept(i)).gt.ABS(MaxHElNotAccept)) MaxHElNotAccept=ABS(AllMaxHElNotAccept(i))
            enddo

            MinHElAccept=0.D0
            IF(AllStats(4).gt.0.D0) THEN
                do i=1,nProcessors
                    IF(AllMinHElAccept(i).ne.0.D0) THEN
                        MinHElAccept=ABS(AllMinHElAccept(i))
                        EXIT
                    ENDIF
                enddo
                do i=1,nProcessors
                    IF((AllMinHElAccept(i).ne.0.D0).and.(ABS(AllMinHElAccept(i)).lt.ABS(MinHElAccept))) MinHElAccept=ABS(AllMinHElAccept(i))
                enddo
            ENDIF

            WRITE(84,'(I10,2F20.1,5F20.6)') Iteration,AllStats(3),AllStats(4),AllStats(3)/AllStats(4),AllStats(1)/AllStats(3),AllStats(2)/AllStats(4),MaxHElNotAccept,MinHElAccept
        ENDIF

    ENDSUBROUTINE PrintSpawnAttemptStats

    SUBROUTINE PrintTriConnStats(Iteration)
        REAL*8 :: TriConStats(4),AllTriConStats(4)
        INTEGER :: error,Iteration

!Write to files the sum of the sign coherent and incoherent triangles. 
        TriConStats(1)=NoSignCohTri
        TriConStats(2)=NoSignIncohTri
        TriConStats(3)=SignCohTri
        TriConStats(4)=SignIncohTri
        AllTriConStats(:)=0.D0

#ifdef PARALLEL
        CALL MPI_Reduce(TriConStats,AllTriConStats,4,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
#else
        AllTriConStats=TriConStats
#endif

        IF(iProcIndex.eq.Root) THEN
            WRITE(77,"(I12,2F24.2,8F20.10)") Iteration,AllTriConStats(1),AllTriConStats(2),AllTriConStats(3),AllTriConStats(4),(AllTriConStats(3)/(Iteration)),&
                                             &(AllTriConStats(4)/(Iteration)),(AllTriConStats(3)/AllTriConStats(1)),(AllTriConStats(4)/AllTriConStats(2)),&
                                             &(AllTriConStats(1)/AllTriConStats(2)),(ABS(AllTriConStats(3)/AllTriConStats(4)))
        ENDIF

    ENDSUBROUTINE PrintTriConnStats
 
 
    SUBROUTINE PrintTriConnHist()
        INTEGER :: error,i
        CHARACTER(len=*), PARAMETER :: this_routine='PrintTriConnHist'

#ifdef PARALLEL
        CALL MPI_Reduce(SignCohTriHist,AllSignCohTriHist,2*NoTriConBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SignIncohTriHist,AllSignIncohTriHist,2*NoTriConBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)

        CALL MPI_Reduce(SignCohHFTriHist,AllSignCohHFTriHist,2*NoTriConBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(SignIncohHFTriHist,AllSignIncohHFTriHist,2*NoTriConBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
#else
        AllSignCohTriHist=SignCohTriHist
        AllSignIncohTriHist=SignIncohTriHist
        AllSignCohHFTriHist=SignCohHFTriHist
        AllSignIncohHFTriHist=SignIncohHFTriHist
#endif

        IF(iProcIndex.eq.Root) THEN
            OPEN(78,file='TriConnHistograms',status='unknown')
            WRITE(78,"(4A25)") "Sign Coh Bin Value","No. in bin","SignIncoh Bin Value","No. in bin"
            OPEN(79,file='TriConnHFHistograms',status='unknown')
            WRITE(79,"(4A25)") "Sign Coh Bin Value","No. in bin","SignIncoh Bin Value","No. in bin"
 
            do i=1,NoTriConBins
                IF((AllSignCohTriHist(2,i).ne.0.D0).or.(AllSignIncohTriHist(2,i).ne.0.D0)) THEN
                    WRITE(78,"(4F25.10)") SignCohTriHist(1,i),AllSignCohTriHist(2,i),SignIncohTriHist(1,i),AllSignIncohTriHist(2,i)
                ENDIF
                IF((AllSignCohHFTriHist(2,i).ne.0.D0).or.(AllSignIncohHFTriHist(2,i).ne.0.D0)) THEN
                    WRITE(79,"(4F25.10)") SignCohHFTriHist(1,i),AllSignCohHFTriHist(2,i),SignIncohHFTriHist(1,i),AllSignIncohHFTriHist(2,i)
                ENDIF
            enddo 
            CLOSE(78)
            CLOSE(79)
        ENDIF

#ifdef PARALLEL
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
#endif

        DEALLOCATE(SignIncohTriHist)
        CALL LogMemDealloc(this_routine,SignIncohTriHistTag)
        DEALLOCATE(SignCohTriHist)
        CALL LogMemDealloc(this_routine,SignCohTriHistTag)
        DEALLOCATE(SignCohHFTriHist)
        CALL LogMemDealloc(this_routine,SignCohHFTriHistTag)
        DEALLOCATE(SignIncohHFTriHist)
        CALL LogMemDealloc(this_routine,SignIncohHFTriHistTag)

        IF(iProcIndex.eq.Root) THEN
            DEALLOCATE(AllSignCohTriHist)
            CALL LogMemDealloc(this_routine,AllSignCohTriHistTag)
            DEALLOCATE(AllSignIncohTriHist)
            CALL LogMemDealloc(this_routine,AllSignIncohTriHistTag)
            DEALLOCATE(AllSignCohHFTriHist)
            CALL LogMemDealloc(this_routine,AllSignCohHFTriHistTag)
            DEALLOCATE(AllSignIncohHFTriHist)
            CALL LogMemDealloc(this_routine,AllSignIncohHFTriHistTag)
        ENDIF


    ENDSUBROUTINE PrintTriConnHist


    SUBROUTINE PrintTriConnHElHist()
        INTEGER :: error,i
        CHARACTER(len=*), PARAMETER :: this_routine='PrintTriConnHElHist'
        REAL*8 :: AllTriConHEls(3,2)


        AllTriConHEls(:,:)=0.D0
#ifdef PARALLEL
        CALL MPI_Reduce(TriConHEls,AllTriConHEls,6,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
#else
        AllTriConHEls=TriConHEls
#endif
        ! TriConHEls(1,1) - number of singles
        ! TriConHEls(1,2) - sum of single elements
        ! TriConHEls(2,1) - number of doubles
        ! TriConHEls(2,2) - sum of double elements
        ! TriConHEls(3,1) - number of Hjk elements
        ! TriConHEls(3,2) - sum of Hjk elements
        IF(iProcIndex.eq.Root) THEN
            WRITE(6,*) "***"
            WRITE(6,*) "*** Stats for determinants connected in triangular forms. ***"
            WRITE(6,*) "Number of single H elements included in the histograms : ",AllTriConHEls(1,1)
            WRITE(6,*) "These elements sum to : ",AllTriConHEls(1,2)
            WRITE(6,*) "which amounts to an average SINGLE H element size of : ",(AllTriConHEls(1,2)/AllTriConHEls(1,1)) 
            WRITE(6,*) "***"
            WRITE(6,*) "Number of double H elements included in the histograms : ",AllTriConHEls(2,1)
            WRITE(6,*) "These elements sum to : ",AllTriConHEls(2,2)
            WRITE(6,*) "which amounts to an average DOUBLE H element size of : ",(AllTriConHEls(2,2)/AllTriConHEls(2,1)) 
            WRITE(6,*) "***"
            WRITE(6,*) "The average size of all H elements is then : ",((AllTriConHEls(2,2)+AllTriConHEls(1,2))/(AllTriConHEls(1,1)+AllTriConHEls(2,1)))
            WRITE(6,*) "***"
            WRITE(6,*) "***"
            WRITE(6,*) "Number of Hjk elements histogrammed : ",AllTriConHEls(3,1)
            WRITE(6,*) "These elements sum to : ",AllTriConHEls(3,2)
            WRITE(6,*) "which amounts to an average Hjk elements size of : ",(AllTriConHEls(3,2)/AllTriConHEls(3,1))
            WRITE(6,*) "***"
        ENDIF

#ifdef PARALLEL
        CALL MPI_Reduce(TriConnHElHistSing,AllTriConnHElHistSing,2*NoTriConHElBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(TriConnHElHistDoub,AllTriConnHElHistDoub,2*NoTriConHElBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)

        CALL MPI_Reduce(TriHjkHistSing,AllTriHjkHistSing,2*NoTriConHElBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
        CALL MPI_Reduce(TriHjkHistDoub,AllTriHjkHistDoub,2*NoTriConHElBins,MPI_DOUBLE_PRECISION,MPI_SUM,Root,MPI_COMM_WORLD,error)
#else
        AllTriConnHElHistSing=TriConnHElHistSing
        AllTriConnHElHistDoub=TriConnHElHistDoub
        AllTriHjkHistSing=TriHjkHistSing
        AllTriHjkHistDoub=TriHjkHistDoub
#endif

        IF(iProcIndex.eq.Root) THEN
            OPEN(80,file='TriConHElHistSing',status='unknown')
            WRITE(80,"(2A25)") "Bin Value","No. in bin"
            OPEN(81,file='TriConHElHistDoub',status='unknown')
            WRITE(81,"(2A25)") "Bin Value","No. in bin"
 
            OPEN(82,file='TriHjkHistSing',status='unknown')
            WRITE(82,"(2A25)") "Bin Value","No. in bin"
            OPEN(83,file='TriHjkHistDoub',status='unknown')
            WRITE(83,"(2A25)") "Bin Value","No. in bin"

            do i=1,NoTriConHElBins
                IF(AllTriConnHElHistSing(2,i).ne.0.D0) THEN
                    WRITE(80,"(2F25.10)") TriConnHElHistSing(1,i),AllTriConnHElHistSing(2,i)
                ENDIF
                IF(AllTriConnHElHistDoub(2,i).ne.0.D0) THEN
                    WRITE(81,"(2F25.10)") TriConnHElHistDoub(1,i),AllTriConnHElHistDoub(2,i)
                ENDIF
                IF(AllTriHjkHistSing(2,i).ne.0.D0) THEN
                    WRITE(82,"(2F25.10)") TriHjkHistSing(1,i),AllTriHjkHistSing(2,i)
                ENDIF
                IF(AllTriHjkHistDoub(2,i).ne.0.D0) THEN
                    WRITE(83,"(2F25.10)") TriHjkHistDoub(1,i),AllTriHjkHistDoub(2,i)
                ENDIF
            enddo 
            CLOSE(80)
            CLOSE(81)
            CLOSE(82)
            CLOSE(83)
        ENDIF

#ifdef PARALLEL
        CALL MPI_Barrier(MPI_COMM_WORLD,error)
#endif

        DEALLOCATE(TriConnHElHistSing)
        CALL LogMemDealloc(this_routine,TriConnHElHistSingTag)
        DEALLOCATE(TriConnHElHistDoub)
        CALL LogMemDealloc(this_routine,TriConnHElHistDoubTag)
        DEALLOCATE(TriHjkHistSing)
        CALL LogMemDealloc(this_routine,TriHjkHistSingTag)
        DEALLOCATE(TriHjkHistDoub)
        CALL LogMemDealloc(this_routine,TriHjkHistDoubTag)

        IF(iProcIndex.eq.Root) THEN
            DEALLOCATE(AllTriConnHElHistSing)
            CALL LogMemDealloc(this_routine,AllTriConnHElHistSingTag)
            DEALLOCATE(AllTriConnHElHistDoub)
            CALL LogMemDealloc(this_routine,AllTriConnHElHistDoubTag)
            DEALLOCATE(AllTriHjkHistSing)
            CALL LogMemDealloc(this_routine,AllTriHjkHistSingTag)
            DEALLOCATE(AllTriHjkHistDoub)
            CALL LogMemDealloc(this_routine,AllTriHjkHistDoubTag)
        ENDIF
        
    ENDSUBROUTINE PrintTriConnHElHist


!These are available to both serial and parallel
    SUBROUTINE InitSpinCoupHEl()

        NoNegSpinCoup=0.D0
        NoPosSpinCoup=0.D0
        SumNegSpinCoup=0.D0
        SumPosSpinCoup=0.D0
        SumHFCon=0.D0
        SumSpinCon=0.D0
        
        IF(iProcIndex.eq.root) THEN
            OPEN(87,file='SpinCoupHEl',status='unknown')
!            WRITE(87,'(A8,10A19)') "1.Iter","2.No.Pos HEls","3.No.Neg HEls","4.Sum Pos HEl","5.Sum Neg HEl","6.Net Sum HEl","7.No.Pos/Iter","8.No.Neg/Iter","9.Sum Pos/Iter","10.Sum Neg/Iter","11.Net Sum/Iter"
            WRITE(87,'(A8,11A18)') "1.Iter","2.No.Pos HEls","3.No.Neg HEls","4.Sum Pos HEl","5.Sum Neg HEl","6.No.Pos/Iter","7.No.Neg/Iter","8.Sum Pos/Iter","9.Sum Neg/Iter",&
            &"10.Sum HF HEls","11.Sum SpinCoup","12.HF HEl/SpinHEl"
        ENDIF

    ENDSUBROUTINE InitSpinCoupHEl

ENDMODULE FciMCLoggingMod

