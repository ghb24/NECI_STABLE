MODULE MCStat
      use constants, only: dp,int64
      IMPLICIT NONE
      TYPE BlockStats
         real(dp), POINTER :: wCurBlock(:)  !(0:)
         real(dp), POINTER :: wBlockSum(:)  !(0:)
         real(dp), POINTER :: wBlockSumSq(:)  !(0:)
         INTEGER                       iBlocks
         INTEGER                       iBMax
         !Histogram average value of block:
         Integer                     bucketCount
         integer(int64), Pointer       :: buckets(:,:)  !(0:,0:)
         real(dp), Pointer          :: bucketMin(:)  !(0:)
         real(dp), Pointer          :: bucketMax(:)  !(0:)
      END TYPE
      TYPE BlockStatsCov
         real(dp), POINTER :: wBlockSum(:)  !(0:)
         INTEGER                       iBlocks
         INTEGER                       iBMax
      END TYPE
      TYPE MCStats
         TYPE(BlockStats)           BlockDeltaSign
         TYPE(BlockStats)           BlockSign
         TYPE(BlockStats)           BlockRatio
         TYPE(BlockStatsCov)        BlockSignDeltaSign
         !  Magically, F90 will know the relevant numbers of rows and columns in this once it has been created.
         integer(int64), POINTER         :: nGen(:,:)  !(0:,0:)
         integer(int64), POINTER         :: nAcc(:,:)  !(0:,0:)
         real(dp), POINTER :: wWeighting(:)  !(0:)
         real(dp), POINTER :: wWeightingSq(:)  !(0:)
         real(dp), POINTER :: wGraphWeight(:)  !(0:)
         real(dp), POINTER :: wGraphWeightSq(:)  !(0:)
         real(dp), POINTER :: wDelta(:)  !(0:)
         real(dp), POINTER :: wDeltaSq(:)  !(0:)
         real(dp), POINTER :: wWeightedDelta(:)  !(0:)
         real(dp), POINTER :: wWeightedDeltaSq(:)  !(0:)
         real(dp), POINTER :: wTrees(:)  !(0:)
         real(dp), POINTER :: wNonTreesPos(:)  !(0:)
         real(dp), POINTER :: wNonTreesNeg(:)  !(0:)
         integer(int64), POINTER         :: nGraphs(:)  !(0:)
         integer(int64), POINTER         :: nNonTreesNeg(:)  !(0:)
         integer(int64), POINTER         :: nNonTreesPos(:)  !(0:)
         integer(int64), POINTER         :: nTrees(:)  !(0:)
         integer(int64)                     iAccTot
         INTEGER                       iVMax
         integer(int64)                     iSeqLen
         INTEGER                       nSeqs
         real(dp)                        fSeqLenSq
         real(dp)               woWeight
         real(dp)               woDelta
         INTEGER                       ioClass
         real(dp)               wRefValue
         real(dp)               wRefWeight
         real(dp)                        foProb
      END TYPE
      CONTAINS
         SUBROUTINE GetStats(MCS,iV,wAvgWeighting,wAvgWeightedValue,wAvgDelta)
            TYPE(MCStats) MCS
            INTEGER iV
            real(dp) wAvgWeighting,wAvgWeightedValue,wAvgDelta
            wAvgWeighting=MCS%wWeighting(iV)/(0.0_dp+MCS%nGraphs(iV))
            wAvgWeightedValue=MCS%wRefValue*wAvgWeighting+MCS%wWeightedDelta(iV)/(0.0_dp+MCS%nGraphs(iV))
            wAvgDelta=MCS%wDelta(iV)/(0.0_dp+MCS%nGraphs(iV))
         END subroutine
!  The constructor
         SUBROUTINE CreateBlockStats(BS,iBlocks)
            TYPE(BlockStats) BS
            INTEGER iBlocks
            real(dp) bucketWidth
            real(dp) bucketCentre
            Integer i
            ALLOCATE(BS%wCurBlock(0:iBlocks))
            ALLOCATE(BS%wBlockSum(0:iBlocks))
            ALLOCATE(BS%wBlockSumSq(0:iBlocks))
            BS%wCurBlock=(0.0_dp)
            BS%wBlockSum=(0.0_dp)
            BS%wBlockSumSq=(0.0_dp)
            BS%iBlocks=iBlocks
            !Histogram setup
            BS%bucketCount = 250
            bucketCentre = -0.07603822
            Allocate(BS%buckets(0:iBlocks,0:BS%bucketCount))
            Allocate(BS%bucketMin(0:iBlocks))
            Allocate(BS%bucketMax(0:iBlocks))
            Do i=0,iBlocks
               bucketWidth=2*15/(2**i)**0.5
               !bucketWidth=3*3*3*0.452*(2**i)**-0.23
               BS%bucketMin(i)=bucketCentre - 0.5*bucketWidth
               BS%bucketMax(i)=bucketCentre + 0.5*bucketWidth
            EndDo
         END subroutine
         SUBROUTINE CreateBlockStatsCov(BSC,iBlocks)
            TYPE(BlockStatsCov) BSC
            INTEGER iBlocks
            ALLOCATE(BSC%wBlockSum(0:iBlocks))
            BSC%wBlockSum=(0.0_dp)
            BSC%iBlocks=iBlocks
         END subroutine
         SUBROUTINE DestroyBlockStatsCov(BSC)
            TYPE(BlockStatsCov) BSC
            DEALLOCATE(BSC%wBlockSum)
         END subroutine
         SUBROUTINE DestroyBlockStats(BS)
            TYPE(BlockStats) BS
            DEALLOCATE(BS%wCurBlock)
            DEALLOCATE(BS%wBlockSum)
            DEALLOCATE(BS%wBlockSumSq)
         END subroutine
         SUBROUTINE Create(MCS,iV,iMaxCycles,wRefValue,wRefWeight)
            TYPE(MCStats) MCS
            INTEGER iV
            INTEGER iMaxCycles
            INTEGER i
            INTEGER iBlocks
            real(dp) wRefValue, wRefWeight
            iBlocks=iMaxCycles
!LOG(iMaxCycles+0.0_dp)/LOG(2.0_dp)+1
            MCS%iVMax=iV
            i=1
            CALL CreateBlockStats(MCS%BlockDeltaSign,iBlocks)
            CALL CreateBlockStats(MCS%BlockSign,iBlocks)
            CALL CreateBlockStats(MCS%BlockRatio,iBlocks)
            CALL CreateBlockStatsCov(MCS%BlockSignDeltaSign,iBlocks)
            ALLOCATE(MCS%nGen(0:iV,0:iV))
            ALLOCATE(MCS%nAcc(0:iV,0:iV))
            MCS%nGen=0
            MCS%nAcc=0
            ALLOCATE(MCS%wWeighting(0:iV))
            ALLOCATE(MCS%wWeightingSq(0:iV))
            ALLOCATE(MCS%wGraphWeight(0:iV))
            ALLOCATE(MCS%wGraphWeightSq(0:iV))
            ALLOCATE(MCS%wDelta(0:iV))
            ALLOCATE(MCS%wDeltaSq(0:iV))
            ALLOCATE(MCS%wWeightedDelta(0:iV))
            ALLOCATE(MCS%wWeightedDeltaSq(0:iV))
            MCS%wWeighting=(0.0_dp)
            MCS%wWeightingSq=(0.0_dp)
            MCS%wGraphWeight=(0.0_dp)
            MCS%wGraphWeightSq=(0.0_dp)
            MCS%wDelta=(0.0_dp)
            MCS%wWeightedDelta=(0.0_dp)
            MCS%wDeltaSq=(0.0_dp)
            MCS%wWeightedDeltaSq=(0.0_dp)
            ALLOCATE(MCS%wTrees(0:iV))
            ALLOCATE(MCS%wNonTreesPos(0:iV))
            ALLOCATE(MCS%wNonTreesNeg(0:iV))
            MCS%wTrees=(0.0_dp)
            MCS%wNonTreesPos=(0.0_dp)
            MCS%wNonTreesneg=(0.0_dp)
            ALLOCATE(MCS%nGraphs(0:iV))
            ALLOCATE(MCS%nNonTreesNeg(0:iV))
            ALLOCATE(MCS%nNonTreesPos(0:iV))
            ALLOCATE(MCS%nTrees(0:iV))
            MCS%nGraphs=0
            MCS%nNonTreesNeg=0
            MCS%nNonTreesPos=0
            MCS%nTrees=0
            MCS%iSeqLen=0
            MCS%nSeqs=0
            MCS%fSeqLenSq=0.0_dp
            MCS%wRefValue=wRefValue
            MCS%wRefWeight=wRefWeight
         END subroutine
         SUBROUTINE Delete(MCS)
            TYPE(MCStats) MCS
            DEALLOCATE(MCS%nGen)
            DEALLOCATE(MCS%nAcc)
            DEALLOCATE(MCS%wWeighting)
            DEALLOCATE(MCS%wWeightingSq)
            DEALLOCATE(MCS%wGraphWeight)
            DEALLOCATE(MCS%wGraphWeightSq)
            DEALLOCATE(MCS%wDelta)
            DEALLOCATE(MCS%wDeltaSq)
            DEALLOCATE(MCS%wWeightedDelta)
            DEALLOCATE(MCS%wWeightedDeltaSq)
            DEALLOCATE(MCS%wTrees)
            DEALLOCATE(MCS%wNonTreesPos)
            DEALLOCATE(MCS%wNonTreesNeg)
            DEALLOCATE(MCS%nGraphs)
            DEALLOCATE(MCS%nNonTreesNeg)
            DEALLOCATE(MCS%nNonTreesPos)
            DEALLOCATE(MCS%nTrees)
            CALL DestroyBlockStats(MCS%BlockDeltaSign)
            CALL DestroyBlockStats(MCS%BlockSign)
            CALL DestroyBlockStats(MCS%BlockRatio)
            CALL DestroyBlockStatsCov(MCS%BlockSignDeltaSign)
            
         END subroutine
!  What is passed in as wSETilde is the sign of the weight times Etilde.  wSign is the sign of the weight
!  We store Delta=ETilde-ETReference
         SUBROUTINE AddGraph(M,nTimes, iV, wWeighting,wValue,wGraphWeight,iClass,iTree,iAcc,ioV,igV,tLog,fProb,TBLOCKING)
            IMPLICIT NONE
            TYPE(MCStats) M
            integer(int64) nTimes
            real(dp) fProb
            INTEGER iV,iTree,iAcc,ioV,igV,iClass
            real(dp) wWeighting,wValue,wWeightedValue,wGraphWeight,wDelta
            real(dp) cc,ave1,ave2,hh,top,bot,calc
            LOGICAL tLog,tNewSeq,tNewPower,TBLOCKING
            tNewPower=.true.    !Really, this should be only true occasionally when we want to write stats out,
                                !But since this is decreciated code, this will remove compile warnings.
            wWeightedValue=wWeighting*wValue
            IF(M%nGraphs(0).EQ.0.OR.iAcc.EQ.0.OR.(iV.EQ.ioV.AND.iV.EQ.1)) THEN
               M%iSeqLen=M%iSeqLen+nTimes
               tNewSeq=.FALSE.
            ELSE
               M%nSeqs=M%nSeqs+1
               M%fSeqLenSq=M%fSeqLenSq+(M%iSeqLen+0.0_dp)**2
               tNewSeq=.TRUE.
            ENDIF
            IF(tNewSeq) THEN
               ave1=M%wWeightedDelta(0)/M%nGraphs(0)
               top=ave1
!               std1=sqrt(M%wDeltaSq(0)/M%nGraphs(0)-ave1*ave1)
               !M%wWeighting(0) gives the running sum of the weightings for the vertex level begin sampled
                ave2=M%wWeighting(0)/M%nGraphs(0)
               !bot is Wref + running average of the weighting - wRefWeight is the running sum of weights over all previous vertex levels
               bot=M%wRefWeight+ave2
!               std2=sqrt(1.0_dp-ave2*ave2)
!               cc=std1/ave1+std2/ave2
               CALL CalcStDev(M,cc)
               hh=ave1/ave2
               !calc gives us the running deltas for the MC run
               calc=(top/bot)+M%wRefValue
               !Gives convergence of energy estimator and denominator
!               write (22,"(2G25.16)") calc, bot
!               write(12,*) M%nGraphs(0), calc !a running count of how delta changes
               
                !To create movie
!                WRITE(6,*) M%wRefValue  
!                IF((ioV.eq.2).and.(MOD(M%nGraphs(0),500).eq.0)) THEN    
!              !      WRITE(102,*) "*********"
!                    WRITE(102,"(I12,I4,2G25.16)") M%nGraphs(0),ioV,calc,cc
!                ENDIF
!                IF((ioV.eq.3).and.(MOD(M%nGraphs(0),20000).eq.0)) THEN
!                    WRITE(103,"(I12,I4,2G25.16)") M%nGraphs(0),ioV,calc,cc
!                ENDIF
!                IF((ioV.eq.4).and.(MOD(M%nGraphs(0),1000).eq.0)) THEN
!                    WRITE(104,"(I12,I4,2G25.16)") M%nGraphs(0),iov,calc,cc
!                ENDIF
                
!VMC file - No. graphs, Sequence length, vertex level, Class, 
               IF(tLog) WRITE(22,"(I20,I15,2I3,8G25.16)") M%nGraphs(0),M%iSeqLen,ioV,M%ioClass,M%woWeight, &
      &             M%woDelta,ave2,ave1,cc,M%foProb,hh,calc
               M%iSeqLen=1
            ENDIF
            IF(iAcc.GT.0) THEN
               M%nAcc(ioV,igV)=M%nAcc(ioV,igV)+nTimes
               M%iAccTot=M%iAccTot+nTimes
            ENDIF 
            wDelta=wValue-M%wRefValue
            CALL AddN(M%nGraphs,iV,nTimes) 
            CALL AddWS(M%wGraphWeight,M%wGraphWeightSq,iV,nTimes,wGraphWeight)
            CALL AddWS(M%wWeighting,M%wWeightingSq,iV,nTimes,wWeighting)
            CALL AddWS(M%wDelta,M%wDeltaSq,iV,nTimes,wDelta)
            CALL AddWS(M%wWeightedDelta,M%wWeightedDeltaSq,iV,nTimes,wDelta*wWeighting)
            IF(iTree.EQ.1) THEN
               CALL AddN(M%nTrees,iV,nTimes) 
               CALL AddW(M%wTrees,iV,nTimes,wGraphWeight)
            ELSE
               IF(wGraphWeight.GT.0.0_dp) THEN
                  CALL AddN(M%nNonTreesPos,iV,nTimes)
                  CALL AddW(M%wNonTreesPos,iV,nTimes,wGraphWeight)
               ELSE
                  CALL AddN(M%nNonTreesNeg,iV,nTimes)
                  CALL AddW(M%wNonTreesNeg,iV,nTimes,wGraphWeight)
               ENDIF
            ENDIF

            !Only do blocking analysis if blocking tag is on
            IF(TBLOCKING) THEN
            
!               CALL AddToBlockStats(M%BlockDeltaSign,wDelta*wSign,nTimes,M%wSDelta(0),M%nGraphs(0),tNewPower)
                CALL AddToBlockStatsII(M%BlockSignDeltaSign,M%BlockSign,M%BlockDeltaSign,M%BlockRatio, &
      &                 wDelta*wWeighting,wWeighting,nTimes,M%wWeightedDelta(0),M%wWeighting(0),M%nGraphs(0),M)
               
                IF(tNewPower.or.(iV.EQ.0)) THEN
!.. Write out the blocking file every time we go past another power of 2

!                   OPEN(23,FILE="MCBLOCKS",STATUS="UNKNOWN")
!                   Call WriteBlockStats(23,M%BlockDeltaSign,M%nGraphs(0)) 
!                   CLOSE(23)
               
!                   OPEN(23,FILE="MCBLOCKSHIST",STATUS="UNKNOWN")
!                   Call WriteHistogram(23,M%BlockDeltaSign)
!                   CLOSE(23)
             
!                   OPEN(23,FILE="MCBLOCKS2",STATUS="UNKNOWN")
!                   Call WriteBlockStats(23,M%BlockSign,M%nGraphs(0)) 
!                   CLOSE(23)

!                   OPEN(23,FILE="MCBLOCKSRatio",STATUS="UNKNOWN")
!                   Call WriteBlockStats(23,M%BlockRatio,M%nGraphs(0))
!                   CLOSE(23)
               
                    OPEN(23,FILE="MCBLOCKS",STATUS="UNKNOWN")
                    CALL WriteBlockStatsII(23,M)
                    CLOSE(23)
                ENDIF
            ENDIF
             
            M%ioClass=iClass
            M%woWeight=wGraphWeight
            M%woDelta=wDelta
            M%foProb=fProb
            M%nGen(ioV,igV)=M%nGen(ioV,igV)+nTimes
         END subroutine
         
      SUBROUTINE WriteBlockStats(iUnit,M,nGraphs)
         TYPE(BlockStats) M
         INTEGER iUnit, i
         integer(int64) nBlocks,nGraphs
         real(dp) blockVar, blockAvg, blockError, blockErrorError
!         WRITE(iUnit,*) "it is WriteBlockStats which is printed" 
         WRITE(iUnit,*) "#MCBLOCKS for ",nGraphs," steps."
         DO i=0,M%iBMax
            nBlocks=Int(nGraphs/2**i)
            blockAvg=M%wBlockSum(i)/nBlocks
            blockVar=M%wBlockSumSq(i)/nBlocks-blockAvg**2
!This /(nBlocks-1) comes from Flyvbjerg's paper, and is because we are working out the best estimator of the stdev,not the actual stdev.
            blockError=SQRT(ABS(blockVar/(nBlocks-1.0_dp)))
            blockErrorError=blockError/SQRT(ABS(2.0_dp*(nBlocks-1.0_dp)))
            WRITE(iUnit,"(I7,4G25.16)") i,blockError,blockErrorError,blockAvg,blockVar
         ENDDO
      END subroutine

      Subroutine CalcStDev(MCStat, estimatedError)
         TYPE(MCStats) MCStat
         real(dp) estimatedError
         Call EstimateError(MCStat, estimatedError, 0)
      End subroutine

      Subroutine EstimateError(MCStat, estimatedError, blockIndex)
         TYPE(MCStats) MCStat
         TYPE(BlockStats) BlockWeightDelta, BlockWeight
         TYPE(BlockStatscov) BlockProduct
         INTEGER blockIndex
         integer(int64) nBlocks
         real(dp) weightAvg, weightDeltaAvg
         real(dp) weightVar, weightDeltaVar, coVar, estimatedVar, estimatedError, jj, kk

         BlockProduct = MCStat%BlockSignDeltaSign
         BlockWeightDelta = MCStat%BlockDeltaSign
         BlockWeight =  MCStat%BlockSign

         nBlocks=Int(MCStat%nGraphs(0)/2**blockIndex)
         !Calculate averages of numerator & denominator
         weightDeltaAvg=(BlockWeightDelta%wBlockSum(blockIndex))/nBlocks
         weightAvg=(BlockWeight%wBlockSum(blockIndex))/nBlocks
         !Calculate variances of n&d and the covariance between them
         weightDeltaVar=(BlockWeightDelta%wBlockSumSq(blockIndex))/nBlocks-weightDeltaAvg**2
         weightVar=(BlockWeight%wBlockSumSq(blockIndex))/nBlocks-weightAvg**2
         coVar=(BlockProduct%wBlockSum(blockIndex))/nBlocks-weightAvg*weightDeltaAvg
         !Estimate the variance & error in the ratio (Data Reduction and Error Analysis for Physical Sceinces (8C.13) p58-64)
         jj=(weightDeltaAvg/weightAvg)**2
         kk=weightDeltaVar/weightDeltaAvg**2+weightVar/weightAvg**2-2*coVar/(weightAvg*weightDeltaAvg)
         estimatedVar=jj*kk
         estimatedError=Sqrt(Abs(estimatedVar/(nBlocks-1.0_dp)))
      End subroutine

      SUBROUTINE WriteBlockStatsII(iUnit,MCStat)
         TYPE(MCStats) MCStat
         TYPE(BlockStats) BlockRatio
         INTEGER iUnit, i
         integer(int64) nBlocks
         real(dp) ratioAvg, ratioVar, estimatedError, ratioError
         BlockRatio = MCStat%BlockRatio
         
!         WRITE(iUnit,*) "It is writeblockstatsII which is printed"
         Write(iUnit,*) "#MCBLOCKS for ",MCStat%nGraphs(0)," steps."
         DO i=0,BlockRatio%iBMax
            nBlocks=Int(MCStat%nGraphs(0)/2**i)
            !Calculate an estimate for the error using blocks of this size (2**i)
            Call EstimateError(MCStat, estimatedError, i)
            !Find average, variance & error in the blocking of the actaul ratio 
            ratioAvg=(BlockRatio%wBlockSum(i))/nBlocks
            ratioVar=(BlockRatio%wBlockSumSq(i))/nBlocks-ratioAvg**2
            ratioError=Sqrt(ratioVar/(nBlocks-1))
            !Write to 'blocking' file:
            Write(iUnit, "(I3, 3G25.16)") i, ratioAvg, ratioError, estimatedError
            !ee=/(SQRT(ABS(2.0_dp*(nBlocks-1.0_dp))))
         ENDDO
      END subroutine

      
!.. Deal with blocking in a new way.
!.. wCurBlock(i) contains the remainder of the sum of values after having removed blocks of length 2**i
!
!.. i.e. For a sequence with indicidual values 1 2 0 2 1 3 2
!..  we would have wCurBlock being
!
!i  2**i   wCurBlock
!4    16        0
!3     8       11=1+2+0+2+1+3+2 
!2     4        6=1+3+2
!1     2        2
!0     1        0
!
!  Here we have the i=4 and subsequent elements =0 (although their actual value is 11) to speed the algorithm.
!  We only deal with blocks up to just larger than our sequence length, but as our sequence grows, we can extend
!  the blocking by copying the value from the largest non-zero block length

!  wBlockSum contains the sum of the average values of all blocks of length up 2**i in the sequence:
!i  2**i   wBlockSum
!4    16        0
!3     8        0 
!2     4        5/4=[(1+2+0+2)]/4
!1     2        9/2=[(1+2)+(0+2)+(1+3)]/2
!0     1        11 =[1+2+0+2+1+3+2]/1
!

!  We now need to be able to add a chain of repeated wSigns of length nTimes
!
!  To do this we consider each level, i, of wCurBlock.  If non-zero, we take enough samples from nTimes to fill to 
!   the length of the block, add these to what's present, divide by the block length, and add to the wBlockSum(i)
!   (and also to wBlockSumSq)
!  We then zero the wCurBlock(i).
!  Now we may find how many complete blocks of length 2**i we can form from the remainder of nTimes, and add 
!  their average (wSign) into wBlockSum that many times.
!  Finally, any remainder of nTimes*wSign is added to wCurBlock(i)

!  An example of adding three samples of value 5 to this system is given below

!i  2**i   wCurBlock(i)   n in wCurBlock(i)    wBlockSum(i)
!4    16        0          0                       0
!3     8       11          7                       0
!2     4        6          3                       5/4
!1     2        2          1                       9/2
!0     1        0          0                       11

!First i=0.  No remainders to deal with and all is added to wBlockSum
!i  2**i   wCurBlock(i)   n in wCurBlock(i)    wBlockSum(i)
!4    16        0          0                       0
!3     8       11          7                       0
!2     4        6          3                       5/4
!1     2        2          1                       9/2
!0     1        0          0                       11+5+5+5

!i=1.  One of the 5s is added to wCurBlock, and then averaged.  The other two form a block of two and are also averaged.
!This leaves no remainder
!i  2**i   wCurBlock(i)   n in wCurBlock(i)    wBlockSum(i)
!4    16        0          0                       0
!3     8       11          7                       0
!2     4        6          3                       5/4
!1     2        0          0                       9/2+(2+5)/2 + (5+5)/2
!0     1        0          0                       26

!i=2.  One of the 5s is added to wCurBlock, and then averaged.  The other two form a remainder.
!i  2**i   wCurBlock(i)   n in wCurBlock(i)    wBlockSum(i)
!4    16        0          0                       0
!3     8       11          7                       0
!2     4        5+5        2                       5/4+(6+5)/4
!1     2        0          0                       26/2
!0     1        0          0                       26

!i=3.  As we're about to go over length 8, we extend the 8 values into the 16s before we lose them. One of the 5s is added to wCurBlock, and then averaged.  The other two form a remainder.
!i  2**i   wCurBlock(i)   n in wCurBlock(i)    wBlockSum(i)
!4    16       11          7                       0
!3     8       5+5         2                       (11+5)/8
!2     4       10          2                       16/4
!1     2        0          0                       26/2
!0     1        0          0                       26

!i=4.  All of the 5s is added to wCurBlock as a remainder.
!i  2**i   wCurBlock(i)   n in wCurBlock(i)    wBlockSum(i)
!4    16       11+5+5+5    10                       0
!3     8       5+5         2                       (11+5)/8
!2     4       10          2                       16/4
!1     2        0          0                       26/2
!0     1        0          0                       26


    SUBROUTINE AddToBlockStatsII(BSDS,BS,BDS,BlockRatio,wDS,wS,nTimes,wValTotDS,wValTotS,nGraphs,M)
         TYPE(MCStats) M
         TYPE(BlockStats) BS, BDS, BlockRatio
         TYPE(BlockStatsCov) BSDS
         real(dp) wDS,wS,wValTotDS,wValTotS
         integer(int64) no,nc,nt,iBlockSize,nGraphs,nTimes
         real(dp) bbS,bbDS
         INTEGER i,iobmax

            iBlockSize=1
            i=0
            ioBMax=BS%iBMax
            DO WHILE(iBlockSize.LE.nGraphs)
               ! iBlockSize is length of the cur block
               nt=nTimes
               no=nGraphs-nTimes !number of graphs before this cycle
               nc=MOD(no,iBlockSize) !number of samples currently in wCurBlock
! If we don't already have a sum for this block size, get it from the stats
               IF(.NOT.abs(BS%wCurBlock(i)).gt.0.0_dp.AND.no.LT.iBlockSize) BS%wCurBlock(i)=(wValTotS)-(nTimes)*wS
               IF(.NOT.abs(BDS%wCurBlock(i)).gt.0.0_dp.AND.no.LT.iBlockSize) BDS%wCurBlock(i)=(wValTotDS)-(nTimes)*wDS
               IF(nc+nt.GE.iBlockSize.AND.nc.NE.0) THEN
!  Add enough from the new set to fill the old to nn, and send to BlockSum
                  bbS=(wS*(iBlockSize-nc)+BS%wCurBlock(i))/(iBlockSize)
                  bbS=bbS+M%wRefWeight
                  bbDS=(wDS*(iBlockSize-nc)+BDS%wCurBlock(i))/(iBlockSize)
                  BSDS%wBlockSum(i)=BSDS%wBlockSum(i)+(bbS*bbDS)
                  BS%wBlockSum(i)=BS%wBlockSum(i)+bbS
                  BS%wBlockSumSq(i)=BS%wBlockSumSq(i)+(bbS*bbS)
                  BDS%wBlockSum(i)=BDS%wBlockSum(i)+bbDS
                  BDS%wBlockSumSq(i)=BDS%wBlockSumSq(i)+(bbDS*bbDS)
                  If (bbS.ne.0) Then
                      BlockRatio%wBlockSum(i)=BlockRatio%wBlockSum(i)+(bbDS/bbs)
                      BlockRatio%wBlockSumSq(i)=BlockRatio%wBlockSumSq(i)+(bbDS/bbs)*(bbDS/bbs)
                  End If
                  If (bbDS.NE.0) Then
                      Call AddToHistogram(BDS, bbDS, i, 1)
                  EndIf
                  nt=nt-(iBlockSize-nc)
                  BS%wCurBlock(i)=0.0_dp
                  BDS%wCurBlock(i)=0.0_dp
              ENDIF
!  Now add in all blocks of length nn from the remainder
               bbS=(Int(nt/iBlockSize))*wS !there are nn lots of wVal - but we want the average
               If (bbS.ne.0) bbS=bbS+M%wRefWeight
               bbDS=(Int(nt/iBlockSize))*wDS
               BSDS%wBlockSum(i)=BSDS%wBlockSum(i)+(bbS*bbDS)
               BS%wBlockSum(i)=BS%wBlockSum(i)+bbS
               BS%wBlockSumSq(i)=BS%wBlockSumSq(i)+(bbS*bbS)
               BDS%wBlockSum(i)=BDS%wBlockSum(i)+bbDS
               BDS%wBlockSumSq(i)=BDS%wBlockSumSq(i)+(bbDS*bbDS)
               If (bbS.ne.0) Then
                  BlockRatio%wBlockSum(i)=BlockRatio%wBlockSum(i)+(bbDS/bbs)
                  BlockRatio%wBlockSumSq(i)=BlockRatio%wBlockSumSq(i)+(bbDS/bbs)*(bbDS/bbs)
               End If
               If (bbDS.NE.0) Then 
                   Call AddToHistogram(BDS, wDS, i, Int(nt/iBlockSize)) 
               EndIf
!  Add in the remainder
               BS%wCurBlock(i)=BS%wCurBlock(i)+(MOD(nt,iBlockSize))*wS
               BDS%wCurBlock(i)=BDS%wCurBlock(i)+(MOD(nt,iBlockSize))*wDS
               iBlockSize=iBlockSize*2
               i=i+1
            ENDDO
            BS%iBMax=i-2
            BDS%iBMax=i-2
            BlockRatio%iBMax=i-2
    END subroutine

    Subroutine AddToHistogram(BS, newBlock, blockSizeIndex, blockCount)
         Type(BlockStats) BS
         real(dp) newBlock
         Integer blockSizeIndex, bucketIndex, blockCount
         real(dp) bucketWidth, newBlockV

         newBlockV = newBlock         

         If (newBlockV.GT.BS%bucketMax(blockSizeIndex)) Then 
             BS%buckets(blockSizeIndex, BS%bucketCount)=BS%buckets(blockSizeIndex, BS%bucketCount)+blockCount
             !Print *, blockSizeIndex, bucketIndex, BS%buckets(blockSizeIndex, BS%bucketCount)
         ElseIf (newBlockV.LT.BS%bucketMin(blockSizeIndex)) Then
             BS%buckets(blockSizeIndex, 1)=BS%buckets(blockSizeIndex, 1)+blockCount
             !Print *, blockSizeIndex, bucketIndex, BS%buckets(blockSizeIndex, 1)
         Else
             bucketWidth = BS%bucketMax(blockSizeIndex)-BS%bucketMin(blockSizeIndex)
             newBlockV = newBlockV - BS%bucketMin(blockSizeIndex)
             bucketIndex = Int(newBlockV/bucketWidth*Real(BS%bucketCount))
             BS%buckets(blockSizeIndex, bucketIndex)=BS%buckets(blockSizeIndex, bucketIndex)+blockCount
             !Print *, blockSizeIndex, bucketIndex, BS%buckets(blockSizeIndex, bucketIndex)
         EndIf
    End subroutine

    Subroutine WriteHistogram(fileNumber, BS)
         Type(BlockStats) BS
         Integer fileNumber, iBuckets, iBlockSize
         real(dp) bucketWidth         
         
         Do iBuckets=1, BS%bucketCount
             Do iBlockSize=0, BS%iBMax
                 bucketWidth = BS%bucketMax(iBlockSize)-BS%bucketMin(iBlockSize)
                 Write(fileNumber, "(F25.16)", ADVANCE='NO') BS%bucketMin(iBlockSize)+bucketWidth/BS%bucketCount*(iBuckets-0.5)
                 Write(fileNumber, "(I7)", ADVANCE='NO') BS%buckets(iBlockSize, iBuckets)
             EndDo
             Write(fileNumber, "(A)") ""
         EndDo
    End   subroutine

      SUBROUTINE WriteStats(M,iUnit)
         TYPE(MCStats) M
         INTEGER iUnit
         real(dp) rStDev
         real(dp) iC
         iC=M%nGraphs(0)+0.0_dp
         CALL CalcStDev(M,rStDev)
         WRITE(iUnit,"(I3,I10,11G25.16,I10,G25.16)") 0,M%nGraphs(0),                         &
     &                   M%wWeighting(0)/iC,                                             &
     &              SQRT((M%wWeightingSq(0)/iC-(M%wWeighting(0)*M%wWeighting(0)/(iC*iC)))), &
     &                   M%wDelta(0)/iC,                                                       &
     &                   M%wWeightedDelta(0)/M%wWeighting(0),                                   &
     &                   M%wDeltaSq(0)/iC,                                          &
     &                   (M%iAccTot+0.0_dp)/M%nGraphs(0),                              &
     &                   (M%nTrees(0)+0.0_dp)/M%nGraphs(0),                            &
     &                   (M%nNonTreesPos(0)+0.0_dp)/M%nGraphs(0),                      &
     &                   (M%nNonTreesNeg(0)+0.0_dp)/M%nGraphs(0),                      &
     &                   (M%wWeightedDelta(0)/iC),                                         &
     &                   M%wRefValue,                                            &
     &                   M%nSeqs,                                                   &
     &                   rStDev
!SUMDLWDB/ICOUNT
      END subroutine
!  Calculate the Standard deviation of the result stored in MCStats.
!  The result is given by  wETReference + <sign*Delta>/<sign>.
!  The Standard deviation regards sign as V and Delta as U, and uses the co-variance of U and V, and has X=<UV>/<U>

!  (s(X)/X)^2   = (s(UV)/UV)^2 +(s(V)/V)^2 - 2 c(UV,V)^2 /(UV V)

!    c(UV,V)^2    = <((uv)_i -<uv>) (v_i-<v>)> = ... = <(uv)_i v_i> - <uv> <v> = <u> - <uv> <v>

! => (s(X)/X)^2 = (s(UV)/UV)^2 + (s(V)/V)^2 + 2 ( 1 - U/(UV V) )



!OLD
!  (s(UV)/UV)^2 = (s(U)/U)^2 + (s(V)/V)^2 + 2 (c(UV,V)/UV)^2
!         c(UV) = {$nind=$n*(0.5**$corrt);}
         SUBROUTINE CalcStDevOLD(M,rStDev)
            TYPE(MCStats) M
            real(dp) rStDev
            real(dp) x,sxbx2,sx,uv,suv2,u,su2,v,sv2,n,nind
            n=M%nGraphs(0)
            nind=M%nseqs
!            WRITE(6,*) "N,NIND:",N,NIND

            uv=M%wWeightedDelta(0)/n
            suv2=M%wDeltaSq(0)/n-uv*uv
            suv2=suv2/nind
!            WRITE(6,*) "UV, SUV2,%",uv,suv2, suv2/(uv*uv)
         
            v=M%wWeighting(0)/n
            sv2=1.0_dp-v*v
            sv2=sv2/nind
!            WRITE(6,*) "V, SV2,%", v, sv2, sv2/(v*v)

            u=M%wDelta(0)/n
            su2=M%wDeltaSq(0)/n-u*u
            su2=su2/nind
!            WRITE(6,*) "U,SU2", u,su2,su2/(u*u)
!            WRITE(6,*) "", u
!            WRITE(6,*) "C(UV,V)i^2",u-uv*v
            x=uv/v
!            WRITE(6,*) "X", x
!            WRITE(6,*) "2(1-u/(uv v))",2*(1-u/(uv*v))
            sxbx2=suv2/(uv*uv)+sv2/(v*v) !+ 2*(1-u/(uv*v))
!            WRITE(6,*) "SXBX2",SXBX2
            sx=sqrt(abs(sxbx2))*abs(x)
!            WRITE(6,*) "SX",SX
            if(x.eq.0) sx=0
            rStDev=sx
         END  subroutine
            


!.. rejig the sums so the result is Sum w Delta / Sum w
         SUBROUTINE WriteLongStats(M,iUnit,OW,OE,Time)
            INTEGER iUnit
            TYPE(MCStats) M
            CHARACTER(len=20) STR2
            INTEGER I,J
            real(dp) OW,OE,iC
            real(dp) Time,fAveSeqLen
            iC=(M%nGraphs(0))
            WRITE(iUnit,"(I12,2G25.16,F19.7,2I12,G25.12)") M%iVMax,M%wWeighting(0)/iC-OW,M%wWeighting(0)/iC, &
      &             Time,M%nGraphs(0),M%nGraphs(0)-M%nGraphs(1),M%wDelta(0)/iC-OE
            WRITE(STR2,"(A,I5,A)") "(A,",M%iVMax+1,"I)"
            WRITE(iUnit,STR2) "GRAPHS(V)",(M%nGraphs(I),I=0,M%iVMax)
            WRITE(iUnit,STR2) "TREES(V)",(M%nTrees(I),I=0,M%iVMax)
            WRITE(iUnit,STR2) "NON-TR+(V)",(M%nNonTreesPos(I),I=0,M%iVMax)
            WRITE(iUnit,STR2) "NON-TR-(V)",(M%nNonTreesNeg(I),I=0,M%iVMax)
            WRITE(STR2,"(A,I5,A)") "(A,",M%iVMax+1,"G)"
            WRITE(iUnit,STR2) "WGHTT(V)",(M%wTrees(I),I=0,M%iVMax)
            WRITE(iUnit,STR2) "WGHT+(V)",(M%wNonTreesPos(I),I=0,M%iVMax)
            WRITE(iUnit,STR2) "WGHT-(V)",(M%wNonTreesNeg(I),I=0,M%iVMax)
            WRITE(STR2,"(A,I5,A)") "(A,I2,",M%iVMax,"I)"
            DO J=1,M%iVMax
               WRITE(iUnit,STR2) "GEN->",J,(M%nGen(I,J),I=1,M%iVMax)
            ENDDO
            DO J=1,M%iVMax
               WRITE(iUnit,STR2) "ACC->",J,(M%nAcc(I,J),I=1,M%iVMax)
            ENDDO
            WRITE(iUnit,*) "Sequences: ",M%nSeqs
            fAveSeqLen=(M%nGraphs(0)+0.0_dp)/M%nSeqs
            WRITE(iUnit,*) "Seq Len: ",fAveSeqLen,"+-",SQRT((M%fSeqLenSq/M%nSeqs)-fAveSeqLen**2)
         END subroutine

!.. just give the additional components for this vertex level
         SUBROUTINE WriteLongStats2(M,iUnit,OW,Time)
            INTEGER iUnit
            TYPE(MCStats) M
            CHARACTER(len=20) STR2
            INTEGER I,J
            real(dp) OW,iC,wAvgWeighting,wAvgWeightedValue,wAvgDelta
            real(dp) Time,fAveSeqLen
            iC=(M%nGraphs(0))
            Call GetStats(M,0,wAvgWeighting,wAvgWeightedValue,wAvgDelta)
            WRITE(iUnit,"(I12,2G25.16,F19.7,2I12,G25.12)") M%iVMax,wAvgWeighting,wAvgWeighting+OW,Time,M%nGraphs(0), &
      &         M%nGraphs(0)-M%nGraphs(1),wAvgWeightedValue
            WRITE(STR2,"(A,I5,A)") "(A,",M%iVMax+1,"I10)"
            WRITE(iUnit,STR2) "GRAPHS(V)",(M%nGraphs(I),I=0,M%iVMax)
            WRITE(iUnit,STR2) "TREES(V)",(M%nTrees(I),I=0,M%iVMax)
            WRITE(iUnit,STR2) "NON-TR+(V)",(M%nNonTreesPos(I),I=0,M%iVMax)
            WRITE(iUnit,STR2) "NON-TR-(V)",(M%nNonTreesNeg(I),I=0,M%iVMax)
            WRITE(STR2,"(A,I5,A)") "(A,",M%iVMax+1,"G25.16)"
            WRITE(iUnit,STR2) "WGHTT(V)",(M%wTrees(I),I=0,M%iVMax)
            WRITE(iUnit,STR2) "WGHT+(V)",(M%wNonTreesPos(I),I=0,M%iVMax)
            WRITE(iUnit,STR2) "WGHT-(V)",(M%wNonTreesNeg(I),I=0,M%iVMax)
            WRITE(STR2,"(A,I5,A)") "(A,I2,",M%iVMax,"I10)"
            DO J=1,M%iVMax
               WRITE(iUnit,STR2) "GEN->",J,(M%nGen(I,J),I=1,M%iVMax)
            ENDDO
            DO J=1,M%iVMax
               WRITE(iUnit,STR2) "ACC->",J,(M%nAcc(I,J),I=1,M%iVMax)
            ENDDO
            WRITE(iUnit,*) "Sequences: ",M%nSeqs
            fAveSeqLen=(M%nGraphs(0)+0.0_dp)/M%nSeqs
            WRITE(iUnit,*) "Seq Len: ",fAveSeqLen,"+-",SQRT((M%fSeqLenSq/M%nSeqs)-fAveSeqLen**2)
         END subroutine
         SUBROUTINE AddWS(w,wSq,iV,nTimes,wV)
            INTEGER iV
            real(dp) :: w(0:iV),wSq(0:iV)
            real(dp) wV,t
            integer(int64) nTimes
            t=nTimes+0.0_dp
            t=t*wV
            w(0)=w(0)+t
            w(iV)=w(iV)+t
            t=t*wV
            wSq(0)=wSq(0)+t
            wSq(iV)=wSq(iV)+t
         END subroutine
         SUBROUTINE AddW(w,iV,nTimes,wV)
            INTEGER iV
            real(dp) :: w(0:iV)
            real(dp) wV,t
            integer(int64) nTimes
            t=nTimes+0.0_dp
            t=t*wV
            w(0)=w(0)+t
            w(iV)=w(iV)+t
         END subroutine
         SUBROUTINE AddN(n,iV,nV)
            INTEGER iV
            integer(int64) :: n(0:iV),nV
            n(0)=n(0)+nV
            n(iV)=n(iV)+nV
         END subroutine
            
END MODULE MCStat         
         
