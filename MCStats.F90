MODULE MCStats
      USE HElement
      IMPLICIT NONE
      TYPE BlockStats
         TYPE(HDElement), POINTER :: wCurBlock(0:)
         TYPE(HDElement), POINTER :: wBlockSum(0:)
         TYPE(HDElement), POINTER :: wBlockSumSq(0:)
         INTEGER                       iBlocks
         INTEGER                       iBMax
         !Histogram average value of block:
         Integer                     bucketCount
         Integer*8, Pointer       :: buckets(0:,0:)
         Real*8, Pointer          :: bucketMin(0:)
         Real*8, Pointer          :: bucketMax(0:)
      END TYPE
      TYPE BlockStatsCov
         TYPE(HDElement), POINTER :: wBlockSum(0:)
         INTEGER                       iBlocks
         INTEGER                       iBMax
      END TYPE
      TYPE MCStats
         TYPE(BlockStats)           BlockDeltaSign
         TYPE(BlockStats)           BlockSign
         TYPE(BlockStatsCov)           BlockSignDeltaSign
         !  Magically, F90 will know the relevant numbers of rows and columns in this once it has been created.
         INTEGER*8, POINTER         :: nGen(0:,0:)
         INTEGER*8, POINTER         :: nAcc(0:,0:)
         TYPE(HDElement), POINTER :: wSign(0:)
         TYPE(HDElement), POINTER :: wSignSq(0:)
         TYPE(HDElement), POINTER :: wWeight(0:)
         TYPE(HDElement), POINTER :: wWeightSq(0:)
         TYPE(HDElement), POINTER :: wDelta(0:)
         TYPE(HDElement), POINTER :: wSDelta(0:)
         TYPE(HDElement), POINTER :: wDeltaSq(0:)
         TYPE(HDElement), POINTER :: wTrees(0:)
         TYPE(HDElement), POINTER :: wNonTreesPos(0:)
         TYPE(HDElement), POINTER :: wNonTreesNeg(0:)
         INTEGER*8, POINTER         :: nGraphs(0:)
         INTEGER*8, POINTER         :: nNonTreesNeg(0:)
         INTEGER*8, POINTER         :: nNonTreesPos(0:)
         INTEGER*8, POINTER         :: nTrees(0:)
         INTEGER*8                     iAccTot
         INTEGER                       iVMax
         INTEGER*8                     iSeqLen
         INTEGER                       nSeqs
         REAL*8                        fSeqLenSq
         TYPE(HDElement)               woWeight
         TYPE(HDElement)               woDelta
         INTEGER                       ioClass
         TYPE(HDElement)               wETReference
         REAL*8                        foProb

      END TYPE
      CONTAINS
         SUBROUTINE GetStats(MCS,iV,wSign,wETilde,wsETilde)
            TYPE(MCStats) MCS
            INTEGER iV
            TYPE(HDElement) wSign,wETilde,wsETilde
            wSign=MCS%wSign(iV)/HDElement(0.D0+MCS%nGraphs(iV))
            wETilde=MCS%wETReference*wSign+MCS%wSDelta(iV)/HDElement(0.D0+MCS%nGraphs(iV))
            wsETilde=MCS%wDelta(iV)/HDElement(0.D0+MCS%nGraphs(iV))
         END
!  The constructor
         SUBROUTINE CreateBlockStats(BS,iBlocks)
            TYPE(BlockStats) BS
            INTEGER iBlocks
            Real*8 bucketWidth
            Real*8 bucketCentre
            Integer i
            ALLOCATE(BS%wCurBlock(0:iBlocks))
            ALLOCATE(BS%wBlockSum(0:iBlocks))
            ALLOCATE(BS%wBlockSumSq(0:iBlocks))
            BS%wCurBlock=HDElement(0.D0)
            BS%wBlockSum=HDElement(0.D0)
            BS%wBlockSumSq=HDElement(0.D0)
            BS%iBlocks=iBlocks
            !Histogram setup
            BS%bucketCount = 50
            bucketCentre = -0.7514263584275787E-01
            Allocate(BS%buckets(0:iBlocks,0:BS%bucketCount))
            Allocate(BS%bucketMin(0:iBlocks))
            Allocate(BS%bucketMax(0:iBlocks))
            Do i=1,iBlocks
               bucketWidth=3*500/(2**iBlocks)**0.5
               BS%bucketMin(i)=bucketCentre - 0.5*bucketWidth
               BS%bucketMax(i)=bucketCentre + 0.5*bucketWidth
            EndDo
         END
         SUBROUTINE CreateBlockStatsCov(BSC,iBlocks)
            TYPE(BlockStatsCov) BSC
            INTEGER iBlocks
            ALLOCATE(BSC%wBlockSum(0:iBlocks))
            BSC%wBlockSum=HDElement(0.D0)
            BSC%iBlocks=iBlocks
         END
         SUBROUTINE DestroyBlockStatsCov(BSC)
            TYPE(BlockStatsCov) BSC
            DEALLOCATE(BSC%wBlockSum)
         END
         SUBROUTINE DestroyBlockStats(BS)
            TYPE(BlockStats) BS
            DEALLOCATE(BS%wCurBlock)
            DEALLOCATE(BS%wBlockSum)
            DEALLOCATE(BS%wBlockSumSq)
         END
         SUBROUTINE Create(MCS,iV,iMaxCycles,wETReference)
            TYPE(MCStats) MCS
            INTEGER iV
            INTEGER iMaxCycles
            INTEGER i,j
            INTEGER iBlocks
            TYPE(HDElement) wETReference
            iBlocks=iMaxCycles
!LOG(iMaxCycles+0.D0)/LOG(2.D0)+1
            MCS%iVMax=iV
            i=HDElementSize
            CALL CreateBlockStats(MCS%BlockDeltaSign,iBlocks)
            CALL CreateBlockStats(MCS%BlockSign,iBlocks)
            CALL CreateBlockStatsCov(MCS%BlockSignDeltaSign,iBlocks)
            ALLOCATE(MCS%nGen(0:iV,0:iV))
            ALLOCATE(MCS%nAcc(0:iV,0:iV))
            MCS%nGen=0
            MCS%nAcc=0
            ALLOCATE(MCS%wSign(0:iV))
            ALLOCATE(MCS%wSignSq(0:iV))
            ALLOCATE(MCS%wWeight(0:iV))
            ALLOCATE(MCS%wWeightSq(0:iV))
            ALLOCATE(MCS%wDelta(0:iV))
            ALLOCATE(MCS%wDeltaSq(0:iV))
            ALLOCATE(MCS%wSDelta(0:iV))
            MCS%wSign=HDElement(0.D0)
            MCS%wSignSq=HDElement(0.D0)
            MCS%wWeight=HDElement(0.D0)
            MCS%wWeightSq=HDElement(0.D0)
            MCS%wDelta=HDElement(0.D0)
            MCS%wSDelta=HDElement(0.D0)
            MCS%wDeltaSq=HDElement(0.D0)
            ALLOCATE(MCS%wTrees(0:iV))
            ALLOCATE(MCS%wNonTreesPos(0:iV))
            ALLOCATE(MCS%wNonTreesNeg(0:iV))
            MCS%wTrees=HDElement(0.D0)
            MCS%wNonTreesPos=HDElement(0.D0)
            MCS%wNonTreesneg=HDElement(0.D0)
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
            MCS%fSeqLenSq=0.D0
            MCS%wETReference=wETReference
         END
         SUBROUTINE Delete(MCS)
            TYPE(MCStats) MCS
            DEALLOCATE(MCS%nGen)
            DEALLOCATE(MCS%nAcc)
            DEALLOCATE(MCS%wSign)
            DEALLOCATE(MCS%wSignSq)
            DEALLOCATE(MCS%wWeight)
            DEALLOCATE(MCS%wWeightSq)
            DEALLOCATE(MCS%wDelta)
            DEALLOCATE(MCS%wDeltaSq)
            DEALLOCATE(MCS%wSDelta)
            DEALLOCATE(MCS%wTrees)
            DEALLOCATE(MCS%wNonTreesPos)
            DEALLOCATE(MCS%wNonTreesNeg)
            DEALLOCATE(MCS%nGraphs)
            DEALLOCATE(MCS%nNonTreesNeg)
            DEALLOCATE(MCS%nNonTreesPos)
            DEALLOCATE(MCS%nTrees)
            CALL DestroyBlockStats(MCS%BlockDeltaSign)
            CALL DestroyBlockStats(MCS%BlockSign)
            CALL DestroyBlockStatsCov(MCS%BlockSignDeltaSign)
            
         END
!  What is passed in as wSETilde is the sign of the weight times Etilde.  wSign is the sign of the weight
!  We store Delta=ETilde-ETReference
         SUBROUTINE AddGraph(M,nTimes, iV, wSign,wSETilde,wWeight,iClass,iTree,iAcc,ioV,igV,tLog,fProb)
            IMPLICIT NONE
            TYPE(MCStats) M
            INTEGER*8 nTimes
            REAL*8 fProb
            INTEGER iV,iTree,iAcc,ioV,igV,iClass
            TYPE(HDElement) wSign,wETilde,wSETilde,wWeight,bb,ss,mm,ee,wVal,wDelta
            INTEGER i,ioBMax,j
            REAL*8 cc,ave1,ave2,std1,std2
            INTEGER*8 no,nc,nt,nn,nnn
            LOGICAL tLog,tNewSeq,tNewPower
            SAVE nnn
            wETilde=wSign*wSETilde
            IF(M%nGraphs(0).EQ.0.OR.iAcc.EQ.0.OR.(iV.EQ.ioV.AND.iV.EQ.1)) THEN
               M%iSeqLen=M%iSeqLen+nTimes
               tNewSeq=.FALSE.
            ELSE
               M%nSeqs=M%nSeqs+1
               M%fSeqLenSq=M%fSeqLenSq+(M%iSeqLen+0.D0)**2
               tNewSeq=.TRUE.
            ENDIF
            IF(tNewSeq) THEN
               ave1=M%wSDelta(0)%v/M%nGraphs(0)
!               std1=sqrt(M%wDeltaSq(0)%v/M%nGraphs(0)-ave1*ave1)
               ave2=M%wSign(0)%v/M%nGraphs(0)
!               std2=sqrt(1.D0-ave2*ave2)
!               cc=std1/ave1+std2/ave2
               CALL CalcStDev(M,cc)
               IF(tLog) WRITE(22,"(I20,I15,2I3,6G25.16)") M%nGraphs(0),M%iSeqLen,ioV,M%ioClass,M%woWeight%v,M%woDelta%v,ave2,ave1,cc,M%foProb
               M%iSeqLen=1
            ENDIF
            IF(iAcc.GT.0) THEN
               M%nAcc(ioV,igV)=M%nAcc(ioV,igV)+nTimes
               M%iAccTot=M%iAccTot+nTimes
            ENDIF 
            wDelta=wETilde-M%wETReference
            CALL AddN(M%nGraphs,iV,nTimes) 
            CALL AddWS(M%wWeight,M%wWeightSq,iV,nTimes,wWeight)
            CALL AddWS(M%wSign,M%wSignSq,iV,nTimes,wSign)
            CALL AddWS(M%wDelta,M%wDeltaSq,iV,nTimes,wDelta)
            CALL AddW(M%wSDelta,iV,nTimes,wDelta*wSign)
            IF(iTree.EQ.1) THEN
               CALL AddN(M%nTrees,iV,nTimes) 
               CALL AddW(M%wTrees,iV,nTimes,wWeight)
            ELSE
               IF(wWeight%v.GT.0.D0) THEN
                  CALL AddN(M%nNonTreesPos,iV,nTimes)
                  CALL AddW(M%wNonTreesPos,iV,nTimes,wWeight)
               ELSE
                  CALL AddN(M%nNonTreesNeg,iV,nTimes)
                  CALL AddW(M%wNonTreesNeg,iV,nTimes,wWeight)
               ENDIF
            ENDIF

!!            CALL AddToBlockStats(M%BlockDeltaSign,wDelta*wSign,nTimes,M%wSDelta(0),M%nGraphs(0),tNewPower)

               CALL AddToBlockStatsII(M%BlockSignDeltaSign,M%BlockSign,M%BlockDeltaSign,wDelta*wSign,wSign,nTimes,M%wSDelta(0),M%wSign(0),M%nGraphs(0),tNewPower)
               
               IF(tNewPower.or.iV.EQ.0) THEN
               OPEN(23,FILE="MCBLOCKS",STATUS="UNKNOWN")
               Call WriteBlockStats(23,M%BlockDeltaSign,M%nGraphs(0)) 
               CLOSE(23)
            ENDIF

            
            !!            CALL AddToBlockStats(M%BlockSign,wSign,nTimes,M%wSign(0),M%nGraphs(0),tNewPower)
!.. Write out the blocking file every time we go past another power of 2
             
            IF(tNewPower.or.iV.EQ.0) THEN
               OPEN(23,FILE="MCBLOCKS2",STATUS="UNKNOWN")
               Call WriteBlockStats(23,M%BlockSign,M%nGraphs(0)) 
               CLOSE(23)
             ENDIF

            IF(tNewPower.or.iV.EQ.0) THEN
               OPEN(23,FILE="MCBLOCKS3",STATUS="UNKNOWN")
               CALL WriteBlockStatsII(23,M%BlockSignDeltaSign,M%BlockDeltaSign,M%BlockSign,M%nGraphs(0))
               CLOSE(23)
            ENDIF
             
            M%ioClass=iClass
            M%woWeight=wWeight
            M%woDelta=wDelta
            M%foProb=fProb
            M%nGen(ioV,igV)=M%nGen(ioV,igV)+nTimes
         END
         
      SUBROUTINE WriteBlockStats(iUnit,M,nGraphs)
         TYPE(BlockStats) M
         INTEGER iUnit
         INTEGER*8 nn,nGraphs
         TYPE(HDElement) mm,ee,ss
         REAL*8 cc
         INTEGER i
               WRITE(iUnit,*) "#MCBLOCKS for ",nGraphs," steps."
               DO i=0,M%iBMax
                  nn=Int(nGraphs/2**i)
                  mm=M%wBlockSum(i)/HDElement(nn+0.D0)
                  cc=DREAL(M%wBlockSumSq(i)/HDElement(nn+0.D0)-mm*mm)
!This /(nn-1) comes from Flyvbjerg's paper, and is because we are working out the best estimator of the stdev,not the actual stdev.
                  ss=SQRT(ABS(cc/(nn-1.D0)))
!                  ss=SQRT(cc/nn)
!ee is the error in the estimator ss.  
                  ee=ss/HDElement(SQRT(ABS(2.D0*(nn-1.D0))))
                  WRITE(iUnit,"(I,4G25.16)") i,ss,ee,mm
                  !WRITE(23,"(2I8,4G25.16)") i,nn,ss,ee,mm
                  !nn=nn/2
               ENDDO
      END

      SUBROUTINE WriteBlockStatsII(iUnit,BSDS,BDS,BS,nGraphs)
      TYPE(BlockStats) BDS,BS
      TYPE(BlockStatsCov) BSDS
      INTEGER iUnit
      INTEGER*8 nn,nGraphs
      TYPE(HDElement) ss,jj,ee,kk,mm,meanS,meanDS,varDS,varS,covar
      REAL*8 cc
      INTEGER i

          WRITE(iUnit,*) "#MCBLOCKS for ",nGraphs," steps."
          DO i=0,BS%iBMax
              nn=Int(nGraphs/2**i)
              mm=BSDS%wBlockSum(i)/HDElement(nn+0.D0)
              meanS=BS%wBlockSum(i)/HDElement(nn+0.D0)
              meanDS=BDS%wBlockSum(i)/HDElement(nn+0.D0)
              varDS=DREAL(BDS%wBlockSumSq(i)/HDElement(nn+0.D0)-(meanDS*meanDS))
              varS=DREAL(BS%wBlockSumSq(i)/HDElement(nn+0.D0)-(meanS*meanS))
              covar=(BSDS%wBlockSum(i)/HDElement(nn+0.D0))-(meanS*meanDS)
              jj=(meanDS/(meanS*HDElement(SQRT(nn-1.D0))))
              kk=(varDS/(meanDS*meanDS))+(varS/(meanS*meanS))-((HDElement(2.D0)*covar)/(meanDS*meanS))              
              ss=ABS((DREAL(jj))*SQRT(DREAL(kk)))
              ee=ss/HDElement(SQRT(ABS(2.D0*(nn-1.D0))))
              WRITE(iUnit,"(I,4G25.16)") i,ss,ee,mm
          ENDDO
      END

      
!.. Deal with blocking in a new way.
!.. wCurBlock(i) contains the remainder of the sum of values after having removed blocks of length 2**i
!
!.. i.e. For a sequence  1 2 0 2 1 3 2
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


    SUBROUTINE AddToBlockStatsII(BSDS,BS,BDS,wDS,wS,nTimes,wValTotDS,wValTotS,nGraphs,tNewPower)
         TYPE(BlockStats) BS, BDS
         TYPE(BlockStatsCov) BSDS
         LOGICAL tNewPower
         TYPE(HDElement) wDS,wS,wValTotDS,wValTotS
         INTEGER*8 no,nc,nt,nn,nGraphs,nTimes
         TYPE(HDElement) bbS,bbDS
         INTEGER i,iobmax

            nn=1
            i=0
            ioBMax=BS%iBMax
            DO WHILE(nn.LE.nGraphs)
               nt=nTimes
! nc is the number of samples in the wCurBlock.
! nn is the length of the cur block
               no=nGraphs-nTimes
               nc=MOD(no,nn)
! If we don't already have a sum for this block size, get it from the stats
               IF(.NOT.(BS%wCurBlock(i).AGT.0.D0).AND.no.LT.nn) BS%wCurBlock(i)=(wValTotS)-HDElement(nTimes)*wS
               IF(.NOT.(BDS%wCurBlock(i).AGT.0.D0).AND.no.LT.nn) BDS%wCurBlock(i)=(wValTotDS)-HDElement(nTimes)*wDS
               IF(nc+nt.GE.nn.AND.nc.NE.0) THEN
!  Add enough from the new set to fill the old to nn, and send to BlockSum
                  bbS=(wS*HDElement(nn-nc)+BS%wCurBlock(i))/HDElement(nn)
                  bbDS=(wDS*HDElement(nn-nc)+BDS%wCurBlock(i))/HDElement(nn)
                  BSDS%wBlockSum(i)=BSDS%wBlockSum(i)+(bbS*bbDS)
                  BS%wBlockSum(i)=BS%wBlockSum(i)+bbS
                  BS%wBlockSumSq(i)=BS%wBlockSumSq(i)+(bbS*bbS)
                  BDS%wBlockSum(i)=BDS%wBlockSum(i)+bbDS
                  BDS%wBlockSumSq(i)=BDS%wBlockSumSq(i)+(bbDS*bbDS)
                  nt=nt-(nn-nc)
                  BS%wCurBlock(i)=0.D0
                  BDS%wCurBlock(i)=0.D0
              ENDIF
!  Now add in all blocks of length nn from the remainder
               bbS=HDElement(Int(nt/nn))*wS !there are nn lots of wVal - but we want the average
               bbDS=HDElement(Int(nt/nn))*wDS
               BSDS%wBlockSum(i)=BSDS%wBlockSum(i)+(bbS*bbDS)
               BS%wBlockSum(i)=BS%wBlockSum(i)+bbS
               BS%wBlockSumSq(i)=BS%wBlockSumSq(i)+(bbS*bbS)
               BDS%wBlockSum(i)=BDS%wBlockSum(i)+bbDS
               BDS%wBlockSumSq(i)=BDS%wBlockSumSq(i)+(bbDS*bbDS)
!  Add in the remainder
               BS%wCurBlock(i)=BS%wCurBlock(i)+HDElement(MOD(nt,nn))*wS
               BDS%wCurBlock(i)=BDS%wCurBlock(i)+HDElement(MOD(nt,nn))*wDS
               nn=nn*2
               i=i+1
            ENDDO
            BS%iBMax=i-2
            BDS%iBMax=i-2
            
    END

    Subroutine AddToHistogram(BS, newBlock, blockSizeIndex)
         Type(BlockStats) BS
         Type(HDElement) newBlock
         Integer blockSizeIndex, bucketIndex
         Real*8 bucketWidth, newBlockV

         newBlockV = newBlock%v         

         If (newBlockV.GT.BS%bucketMax(blockSizeIndex)) Then 
             BS%buckets(blockSizeIndex, BS%bucketCount)=BS%buckets(blockSizeIndex, BS%bucketCount)+1
         ElseIf (newBlockV.LT.BS%bucketMin(blockSizeIndex)) Then
             BS%buckets(blockSizeIndex, 1)=BS%buckets(blockSizeIndex, 1)+1
         Else
             bucketWidth = BS%bucketMax(blockSizeIndex)-BS%bucketMin(blockSizeIndex)
             bucketIndex = Int(newBlockV/bucketWidth*Real(BS%bucketCount))
             BS%buckets(blockSizeIndex, bucketIndex)=BS%buckets(blockSizeIndex, bucketIndex)+1
         EndIf
    End

    Subroutine WriteHistogram(fileNumber, BS)
         Type(BlockStats) BS
         Integer fileNumber, iBuckets, iBlockSize
         
         Do iBuckets=1, BS%bucketCount
             Do iBlockSize=0, BS%iBMax
                 Write(fileNumber, "") BS%buckets(iBlockSize, iBuckets)
             EndDo
             !Write(fileNumber, "A") "\n"
         EndDo
    End  

           !  This routine is now never called 
            SUBROUTINE AddToBlockStats(M,wVal,nTimes,wValTot,nGraphs,tNewPower)
         TYPE(BlockStats) M
         LOGICAL tNewPower
         TYPE(HDElement) wVal,wvalTot
         INTEGER*8 no,nc,nt,nn,nnn,nGraphs,nTimes
         TYPE(HDElement)bb,ss,mm,ee
         REAL*8 cc
         INTEGER i,j,iobmax

!.. We enable the new blocking counting method if TRUE, and the old one if FALSE
         IF(.TRUE.) THEN
            nn=1
            i=0
            ioBMax=M%iBMax
            DO WHILE(nn.LE.nGraphs)
               nt=nTimes
! nc is the number of samples in the wCurBlock.
! nn is the length of the cur block
               no=nGraphs-nTimes
               nc=MOD(no,nn)
! If we don't already have a sum for this block size, get it from the stats
               !IF(.NOT.(M%wCurBlock(i).AGT.0.D0).AND.no.LT.nn) M%wCurBlock(i)=(M%wDelta(0))-HDElement(nTimes)*wVal
               IF(.NOT.(M%wCurBlock(i).AGT.0.D0).AND.no.LT.nn) M%wCurBlock(i)=(wValTot)-HDElement(nTimes)*wVal
               IF(nc+nt.GE.nn.AND.nc.NE.0) THEN
!  Add enough from the new set to fill the old to nn, and send to BlockSum
                  bb=(wVal*HDElement(nn-nc)+M%wCurBlock(i))/HDElement(nn)
                  M%wBlockSum(i)=M%wBlockSum(i)+bb
                  M%wBlockSumSq(i)=M%wBlockSumSq(i)+bb*bb
                  nt=nt-(nn-nc)
                  M%wCurBlock(i)=0.D0
               ENDIF
!  Now add in all blocks of length nn from the remainder
               bb=HDElement(Int(nt/nn))*wVal !there are nn lots of wVal - but we want the average
               M%wBlockSum(i)=M%wBlockSum(i)+bb
               M%wBlockSumSq(i)=M%wBlockSumSq(i)+bb*bb
!  Add in the remainder
               M%wCurBlock(i)=M%wCurBlock(i)+HDElement(MOD(nt,nn))*wVal
               nn=nn*2
               i=i+1
            ENDDO
            M%iBMax=i-2
         ELSE
!.. Old Deal with blocking
            IF(nGraphs.GT.nTimes) THEN
               nnn=nnn+1
               ioBMax=M%iBMax
               DO j=1,nTimes
                  i=0
!            nn=(M%nGraphs(0)-nTimes+j)*2
                  nn=nnn*2
                  M%wCurBlock(0)=M%wCurBlock(0)+wVal
                  DO WHILE(.NOT.BTEST(nn,0))
                     M%wCurBlock(i+1)=M%wCurBlock(i+1)+M%wCurBlock(i)
                     bb=M%wCurBlock(i)/HDElement(0.D0+2**i)
                     M%wBlockSum(i)=M%wBlockSum(i)+bb
                     bb=bb*bb
                     M%wBlockSumSq(i)=M%wBlockSumSq(i)+bb
                     M%wCurBlock(i)=0.D0
                     i=i+1
                     IF(i-2.GT.ioBMax) M%iBMax=i-2
                     IF(M%iBMax.GT.M%iBlocks) STOP "TOO many cycles"
                     nn=nn/2
                  ENDDO
               ENDDO
               DO i=0,M%iBMax+1
                  cc=M%wBlockSum(i)%v/(nnn/2**i)
                  cc=M%wBlockSumSq(i)%v/(nnn/2**i)-cc*cc
                  cc=sqrt(abs(cc/(nnn/2**i-1.D0)))
               ENDDO
            ENDIF
         ENDIF
         tNewPower=M%iBMax.NE.ioBMax
      END
      SUBROUTINE WriteStats(M,iUnit)
         TYPE(MCStats) M
         INTEGER iUnit
         REAL*8 Time,rStDev
         TYPE(HDElement) iC
         iC=M%nGraphs(0)+0.D0
         CALL CalcStDev(M,rStDev)
         WRITE(iUnit,"(I3,I10,11G25.16,I10,G25.16)") 0,M%nGraphs(0),                         &
     &                   M%wSign(0)/iC,                                             &
     &              SQRT(DREAL(M%wSignSq(0)/iC-(M%wSign(0)*M%wSign(0)/(iC*iC)))), &
     &                   M%wDelta(0)/iC,                                                       &
     &                   M%wSDelta(0)/M%wSign(0),                                   &
     &                   M%wDeltaSq(0)/iC,                                          &
     &                   (M%iAccTot+0.D0)/M%nGraphs(0),                              &
     &                   (M%nTrees(0)+0.D0)/M%nGraphs(0),                            &
     &                   (M%nNonTreesPos(0)+0.D0)/M%nGraphs(0),                      &
     &                   (M%nNonTreesNeg(0)+0.D0)/M%nGraphs(0),                      &
     &                   (M%wSDelta(0)/iC),                                         &
     &                   M%wETReference,                                            &
     &                   M%nSeqs,                                                   &
     &                   rStDev
!SUMDLWDB/ICOUNT
      END
!  Calculate the Standard deviation of the result stored in MCStats.
!  The result is given by  wETReference + <sign*Delta>/<sign>.
!  The Standard deviation regards sign as V and Delta as U, and uses the co-variance of U and V, and has X=<UV>/<U>

!  (s(X)/X)^2   = (s(UV)/UV)^2 +(s(V)/V)^2 - 2 c(UV,V)^2 /(UV V)

!    c(UV,V)^2    = <((uv)_i -<uv>) (v_i-<v>)> = ... = <(uv)_i v_i> - <uv> <v> = <u> - <uv> <v>

! => (s(X)/X)^2 = (s(UV)/UV)^2 + (s(V)/V)^2 + 2 ( 1 - U/(UV V) )



!OLD
!  (s(UV)/UV)^2 = (s(U)/U)^2 + (s(V)/V)^2 + 2 (c(UV,V)/UV)^2
!         c(UV) = {$nind=$n*(0.5**$corrt);}
         SUBROUTINE CalcStDev(M,rStDev)
            TYPE(MCStats) M
            REAL*8 rStDev
            REAL*8 x,sxbx2,sx,uv,suv2,u,su2,v,sv2,n,nind
            n=M%nGraphs(0)
            nind=M%nseqs/2
!            WRITE(6,*) "N,NIND:",N,NIND

            uv=M%wSDelta(0)%v/n
            suv2=M%wDeltaSq(0)%v/n-uv*uv
            suv2=suv2/nind
!            WRITE(6,*) "UV, SUV2,%",uv,suv2, suv2/(uv*uv)
         
            v=M%wSign(0)%v/n
            sv2=1.D0-v*v
            sv2=sv2/nind
!            WRITE(6,*) "V, SV2,%", v, sv2, sv2/(v*v)

            u=M%wDelta(0)%v/n
            su2=M%wDeltaSq(0)%v/n-u*u
            su2=su2/nind
!            WRITE(6,*) "U,SU2", u,su2,su2/(u*u)
!            WRITE(6,*) "", u
!            WRITE(6,*) "C(UV,V)i^2",u-uv*v
            x=uv/v
!            WRITE(6,*) "X", x
!            WRITE(6,*) "2(1-u/(uv v))",2*(1-u/(uv*v))
            sxbx2=suv2/(uv*uv)+sv2/(v*v) !+ 2*(1-u/(uv*v))
!            WRITE(6,*) "SXBX2",SXBX2
            sx=sqrt(abs(sxbx2))*x
!            WRITE(6,*) "SX",SX
            if(x.eq.0) sx=0
            rStDev=sx
         END 
            


         SUBROUTINE WriteLongStats(M,iUnit,OW,OE,Time)
            INTEGER iUnit
            TYPE(MCStats) M
            CHARACTER*20 STR2
            INTEGER I,J
            TYPE(HDElement) OW,OE,iC
            REAL*8 Time,fAveSeqLen
            iC=HDElement(M%nGraphs(0))
            WRITE(iUnit,"(I12,2G25.16,F19.7,2I12,G25.12)") ,M%iVMax,M%wSign(0)/iC-OW,M%wSign(0)/iC,Time,M%nGraphs(0),M%nGraphs(0)-M%nGraphs(1),M%wDelta(0)/iC-OE
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
            fAveSeqLen=(M%nGraphs(0)+0.D0)/M%nSeqs
            WRITE(iUnit,*) "Seq Len: ",fAveSeqLen,"+-",SQRT((M%fSeqLenSq/M%nSeqs)-fAveSeqLen**2)
         END
         SUBROUTINE AddWS(w,wSq,iV,nTimes,wV)
            TYPE(HDElement) :: w(0:iV),wSq(0:iV)
            TYPE(HDElement) wV,t
            INTEGER iV
            INTEGER*8 nTimes
            t=nTimes+0.D0
            t=t*wV
            w(0)=w(0)+t
            w(iV)=w(iV)+t
            t=t*wV
            wSq(0)=wSq(0)+t
            wSq(iV)=wSq(iV)+t
         END
         SUBROUTINE AddW(w,iV,nTimes,wV)
            TYPE(HDElement) :: w(0:iV)
            TYPE(HDElement) wV,t
            INTEGER iV
            INTEGER*8 nTimes
            t=nTimes+0.D0
            t=t*wV
            w(0)=w(0)+t
            w(iV)=w(iV)+t
         END
         SUBROUTINE AddN(n,iV,nV)
            INTEGER*8 :: n(0:iV),nV
            INTEGER iV
            n(0)=n(0)+nV
            n(iV)=n(iV)+nV
         END
            
END MODULE MCStats         
         
