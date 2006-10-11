MODULE MCStats
      USE HElement
      IMPLICIT NONE
      TYPE MCStats
         TYPE(HDElement), POINTER :: wCurBlock(0:)
         TYPE(HDElement), POINTER :: wBlockSum(0:)
         TYPE(HDElement), POINTER :: wBlockSumSq(0:)
!  Magically, F90 will know the relevant numbers of rows and columns in this once it has been created.
         INTEGER*8, POINTER         :: nGen(0:,0:)
         INTEGER*8, POINTER         :: nAcc(0:,0:)
         TYPE(HDElement), POINTER :: wValue(0:)
         TYPE(HDElement), POINTER :: wValueSq(0:)
         TYPE(HDElement), POINTER :: wWeight(0:)
         TYPE(HDElement), POINTER :: wWeightSq(0:)
         TYPE(HDElement), POINTER :: wETilde(0:)
         TYPE(HDElement), POINTER :: wETildeSq(0:)
         TYPE(HDElement), POINTER :: wTrees(0:)
         TYPE(HDElement), POINTER :: wNonTreesPos(0:)
         TYPE(HDElement), POINTER :: wNonTreesNeg(0:)
         INTEGER*8, POINTER         :: nGraphs(0:)
         INTEGER*8, POINTER         :: nNonTreesNeg(0:)
         INTEGER*8, POINTER         :: nNonTreesPos(0:)
         INTEGER*8, POINTER         :: nTrees(0:)
         INTEGER                       iBlocks
         INTEGER                       iBMax
         INTEGER*8                     iAccTot
         INTEGER                       iVMax
         INTEGER*8                     iSeqLen
         INTEGER                       nSeqs
         REAL*8                        fSeqLenSq

      END TYPE
      CONTAINS
         SUBROUTINE GetStats(MCS,iV,wValue,wETilde)
            TYPE(MCStats) MCS
            INTEGER iV
            TYPE(HDElement) wValue,wETilde
            wValue=MCS%wValue(iV)/HDElement(0.D0+MCS%nGraphs(iV))
            wETilde=MCS%wETilde(iV)/HDElement(0.D0+MCS%nGraphs(iV))
         END
!  The constructor
         SUBROUTINE Create(MCS,iV,iMaxCycles)
            TYPE(MCStats) MCS
            INTEGER iV
            INTEGER iMaxCycles
            INTEGER i,j
            INTEGER iBlocks
            iBlocks=iMaxCycles
            MCS%iBlocks=iBlocks
!LOG(iMaxCycles+0.D0)/LOG(2.D0)+1
            MCS%iVMax=iV
            i=HDElementSize
            ALLOCATE(MCS%wCurBlock(0:iBlocks))
            ALLOCATE(MCS%wBlockSum(0:iBlocks))
            ALLOCATE(MCS%wBlockSumSq(0:iBlocks))
            MCS%wCurBlock=HDElement(0.D0)
            MCS%wBlockSum=HDElement(0.D0)
            MCS%wBlockSumSq=HDElement(0.D0)
            ALLOCATE(MCS%nGen(0:iV,0:iV))
            ALLOCATE(MCS%nAcc(0:iV,0:iV))
            MCS%nGen=0
            MCS%nAcc=0
            ALLOCATE(MCS%wValue(0:iV))
            ALLOCATE(MCS%wValueSq(0:iV))
            ALLOCATE(MCS%wWeight(0:iV))
            ALLOCATE(MCS%wWeightSq(0:iV))
            ALLOCATE(MCS%wETilde(0:iV))
            ALLOCATE(MCS%wETildeSq(0:iV))
            MCS%wValue=HDElement(0.D0)
            MCS%wValueSq=HDElement(0.D0)
            MCS%wWeight=HDElement(0.D0)
            MCS%wWeightSq=HDElement(0.D0)
            MCS%wETilde=HDElement(0.D0)
            MCS%wETildeSq=HDElement(0.D0)
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
         END
         SUBROUTINE Delete(MCS)
            TYPE(MCStats) MCS
            DEALLOCATE(MCS%wCurBlock)
            DEALLOCATE(MCS%wBlockSum)
            DEALLOCATE(MCS%wBlockSumSq)
            DEALLOCATE(MCS%nGen)
            DEALLOCATE(MCS%nAcc)
            DEALLOCATE(MCS%wValue)
            DEALLOCATE(MCS%wValueSq)
            DEALLOCATE(MCS%wWeight)
            DEALLOCATE(MCS%wWeightSq)
            DEALLOCATE(MCS%wETilde)
            DEALLOCATE(MCS%wETildeSq)
            DEALLOCATE(MCS%wTrees)
            DEALLOCATE(MCS%wNonTreesPos)
            DEALLOCATE(MCS%wNonTreesNeg)
            DEALLOCATE(MCS%nGraphs)
            DEALLOCATE(MCS%nNonTreesNeg)
            DEALLOCATE(MCS%nNonTreesPos)
            DEALLOCATE(MCS%nTrees)
           
         END
         SUBROUTINE AddGraph(M,nTimes, iV, wValue,wETilde,wWeight,iTree,iAcc,ioV,igV)
            TYPE(MCStats) M
            INTEGER*8 nTimes
            INTEGER iV,iTree,iAcc,ioV,igV
            TYPE(HDElement) wValue,wETilde,wWeight,bb,ss,mm,ee,wVal
            INTEGER i,ioBMax,j
            REAL*8 cc
            INTEGER*8 no,nc,nt,nn
            IF(iAcc.EQ.0.OR.(iV.EQ.ioV.AND.iV.EQ.1)) THEN
               M%iSeqLen=M%iSeqLen+nTimes
            ELSE
               M%nSeqs=M%nSeqs+1
               M%fSeqLenSq=M%fSeqLenSq+(M%iSeqLen+0.D0)**2
!               WRITE(6,*) M%fSeqLenSq
               M%iSeqLen=1
            ENDIF
            M%nGen(ioV,igV)=M%nGen(ioV,igV)+nTimes
            IF(iAcc.GT.0) THEN
               M%nAcc(ioV,igV)=M%nAcc(ioV,igV)+nTimes
               M%iAccTot=M%iAccTot+nTimes
            ENDIF 
            CALL AddN(M%nGraphs,iV,nTimes) 
            CALL AddWS(M%wWeight,M%wWeightSq,iV,nTimes,wWeight)
            CALL AddWS(M%wValue,M%wValueSq,iV,nTimes,wValue)
            CALL AddWS(M%wETilde,M%wETildeSq,iV,nTimes,wETilde)
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

!  We now need to be able to add a chain of repeated wValues of length nTimes
!
!  To do this we consider each level, i, of wCurBlock.  If non-zero, we take enough samples from nTimes to fill to 
!   the length of the block, add these to what's present, divide by the block length, and add to the wBlockSum(i)
!   (and also to wBlockSumSq)
!  We then zero the wCurBlock(i).
!  Now we may find how many complete blocks of length 2**i we can form from the remainder of nTimes, and add 
!  their average (wValue) into wBlockSum that many times.
!  Finally, any remainder of nTimes*wValue is added to wCurBlock(i)

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


            wVal=wETilde
!            WRITE(6,*) nTimes,wETilde
            IF(.TRUE.) THEN
            nn=1
            i=0
            ioBMax=M%iBMax
!            WRITE(6,*) M%wETilde
            DO WHILE(nn.LE.M%nGraphs(0))
               nt=nTimes
! nc is the number of samples in the wCurBlock.
! nn is the length of the cur block
               no=M%nGraphs(0)-nTimes
               nc=MOD(no,nn)
! If we don't already have a sum for this block size, get it from the stats
               IF(.NOT.(M%wCurBlock(i).AGT.0.D0).AND.no.LT.nn) M%wCurBlock(i)=M%wETilde(0)-HDElement(nTimes)*wVal
               IF(nc+nt.GE.nn.AND.nc.NE.0) THEN
!  Add enough from the new set to fill the old to nn, and send to BlockSum
                  bb=(wVal*HDElement(nn-nc)+M%wCurBlock(i))/HDElement(nn)
                  M%wBlockSum(i)=M%wBlockSum(i)+bb
                  M%wBlockSumSq(i)=M%wBlockSumSq(i)+bb*bb
                  nt=nt-(nn-nc)
                  M%wCurBlock(i)=0.D0
               ENDIF
!  Now add in all blocks of length nn from the remainder
               bb=HDElement((nt/nn))
               M%wBlockSum(i)=M%wBlockSum(i)+bb*wVal
               M%wBlockSumSq(i)=M%wBlockSumSq(i)+bb*wVal*wVal
!  Add in the remainder
               M%wCurBlock(i)=M%wCurBlock(i)+HDElement(MOD(nt,nn))*wVal
               nn=nn*2
               i=i+1
            ENDDO
            M%iBMax=i-2
            ELSE
!.. Old Deal with blocking
            ioBMax=M%iBMax
            DO j=1,nTimes
            i=0
            nn=(M%nGraphs(0)-nTimes+j)*2
            M%wCurBlock(0)=M%wCurBlock(0)+wVal
            DO WHILE(.NOT.BTEST(nn,0))
               M%wCurBlock(i+1)=M%wCurBlock(i+1)+M%wCurBlock(i)
               bb=M%wCurBlock(i)/HDElement(0.D0+2**i)
               M%wBlockSum(i)=M%wBlockSum(i)+bb
               bb=bb*bb
               M%wBlockSumSq(i)=M%wBlockSumSq(i)+bb
               M%wCurBlock(i)=0.D0
               i=i+1
               IF(i-1.GT.ioBMax) M%iBMax=i-2
               IF(M%iBMax.GT.M%iBlocks) STOP "TOO many cycles"
               nn=nn/2
            ENDDO
            ENDDO
            ENDIF
!            DO i=M%iBMax,0,-1
!               WRITE(6,"(I,3F,I)") i,m%wCurBlock(i),m%wBlockSum(i),m%wBlockSumSq(i),MOD(M%nGraphs(0),2**i)
!            ENDDO
!.. Write out the blocking file every time we go past another power of 2
            IF(M%iBMax.NE.ioBMax) THEN
               OPEN(23,FILE="MCBLOCKS",STATUS="UNKNOWN")
               nn=M%nGraphs(0)
               DO i=0,M%iBMax
                  mm=M%wBlockSum(i)/HDElement(nn+0.D0)
                  cc=DREAL(M%wBlockSumSq(i)/HDElement(nn+0.D0)-mm*mm)
                  ss=SQRT(ABS(cc)/(nn-1.D0))
                  ee=ss/HDElement(SQRT(2.D0*(nn-1.D0)))
                  WRITE(23,"(I,3G25.16)") i,ss,ee,mm
                  nn=nn/2
               ENDDO
               CLOSE(23)
             ENDIF
         END
         SUBROUTINE WriteStats(M,iUnit)
            TYPE(MCStats) M
            INTEGER iUnit
            REAL*8 Time
            TYPE(HDElement) iC
            iC=M%nGraphs(0)+0.D0
            WRITE(iUnit,"(I3,I10,10G25.16)") 0,M%nGraphs(0)/1000,                    &
     &                   M%wValue(0)/iC,                                             &
     &              SQRT(DREAL(M%wValueSq(0)/iC-(M%wValue(0)*M%wValue(0)/(iC*iC)))), &
     &                   0.D0,                                                       &
     &                   M%wETilde(0)/M%wValue(0),                                   &
     &                   M%wETildeSq(0)/iC,                                          &
     &                   (M%iAccTot+0.D0)/M%nGraphs(0),                              &
     &                   (M%nTrees(0)+0.D0)/M%nGraphs(0),                            &
     &                   (M%nNonTreesPos(0)+0.D0)/M%nGraphs(0),                      &
     &                   (M%nNonTreesNeg(0)+0.D0)/M%nGraphs(0),                      &
     &                   0.D0
!SUMDLWDB/ICOUNT
         END
         SUBROUTINE WriteLongStats(M,iUnit,OW,OE,Time)
            INTEGER iUnit
            TYPE(MCStats) M
            CHARACTER*20 STR2
            INTEGER I,J
            TYPE(HDElement) OW,OE,iC
            REAL*8 Time,fAveSeqLen
            iC=HDElement(M%nGraphs(0))
            WRITE(iUnit,"(I12,2G25.16,F19.7,2I12,F19.7)") ,M%iVMax,M%wValue(0)/iC-OW,M%wValue(0)/iC,Time,M%nGraphs(0),M%nGraphs(0)-M%nGraphs(1),M%wETilde(0)/iC-OE
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
         
