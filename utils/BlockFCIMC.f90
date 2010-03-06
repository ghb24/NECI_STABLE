PROGRAM BlkFCIMC

      IMPLICIT NONE
      CHARACTER(Len=1000) :: FirstLineRubbish
      INTEGER :: i,Iter,WalkCng,TotWalkers,Annihil,Died,Born,StartIter,ierr,TotPoints
      INTEGER :: TotBlkSize,NoBlocks,BlockSize,VecInd,Blk,j
      REAL*8 :: Sft,GrowRate,ProjE,AvBlock,SumBlock,Mean,MeanSq,SD,Error,ErrorinError,AvShift,ProjEInst
      REAL*8 :: AvBlockEn,SumBlockEn,MeanEn,MeanSqEn,SDEn,ErrorEn,ErrorinErrorEn
      REAL*8 , ALLOCATABLE :: Shifts(:),Energies(:),BlkAv(:),BlkAvEn(:)

      OPEN(UNIT=11,FILE='FCIMCStats',STATUS='OLD',ACTION='READ',    &
        POSITION='REWIND')
      OPEN(UNIT=12,FILE='BlockingInfo',STATUS='UNKNOWN')
      OPEN(UNIT=13,FILE='BlockingInfoEnergies',STATUS='UNKNOWN')
      WRITE(12,"(A)") '#  Blocklength   Log_2 BL     MeanBlocks              SD-Blocks &
                       &      Err-Blocks           Err in Err-Blocks'
      WRITE(13,"(A)") '#  Blocklength   Log_2 BL     MeanBlocks              SD-Blocks &
                       &      Err-Blocks           Err in Err-Blocks'

      READ(11,*) FirstLineRubbish
    
      WRITE(6,*) "Enter the iteration number from which you want to start the blocking analysis:"
      READ(*,*) StartIter

      i=0
      do while(.true.)
          READ(11,'(I12,G16.7,I10,G16.7,I12,3I13,3G17.9)',END=99) Iter,Sft,WalkCng, &
                GrowRate,TotWalkers,Annihil,Died,Born,ProjE,AvShift,ProjEInst
          IF(Iter.ge.StartIter) THEN
              i=i+1
          ENDIF
      enddo

99    CONTINUE

      TotPoints=i
      WRITE(6,*) "Number of data points to be used in blocking = ", TotPoints

!Now allocate memory for the blocking data
      ALLOCATE(Shifts(TotPoints),stat=ierr)
      IF(ierr.ne.0) THEN
          STOP 'Error in allocation'
      ENDIF
      ALLOCATE(Energies(TotPoints),stat=ierr)
      IF(ierr.ne.0) THEN
          STOP 'Error in allocation'
      ENDIF
      Shifts(:)=0.D0
      Energies(:)=0.D0

!Reread in data to fill arrays
      REWIND 11
      READ(11,*) FirstLineRubbish
      i=0
      do while(.true.)
          READ(11,'(I12,G16.7,I10,G16.7,I12,3I13,3G17.9)',END=98) Iter,Sft,WalkCng, &
                 GrowRate,TotWalkers,Annihil,Died,Born,ProjE,AvShift,ProjEInst
          IF(Iter.ge.StartIter) THEN
              i=i+1
              Shifts(i)=Sft
              Energies(i)=ProjEInst
          ENDIF
      enddo

98    CONTINUE
      CLOSE(11)

!Allocate data to hold averages of each block
      ALLOCATE(BlkAv(TotPoints),stat=ierr)
      ALLOCATE(BlkAvEn(TotPoints),stat=ierr)
      IF(ierr.ne.0) THEN
          STOP 'Error in allocation'
      ENDIF
      BlkAv(:)=0.D0
      BlkAvEn(:)=0.D0

!      do i=1,TotPoints
!          WRITE(7,'(I12,G16.7)') 500,Shifts(i)
!      enddo

!Find information about the maximum block size, and the number of different block sizes
      do i=1,1000
          IF(2**i.gt.TotPoints) EXIT
      enddo

!We subtract one here, since the loop will exit when Blocklength > number of points
!The second one is removed, since one block will not be enough for an estimate of the error.
      TotBlkSize=i-2
      WRITE(6,*) "Largest Blocksize is: ",2**TotBlkSize

      do i=0,TotBlkSize     !Loop over all blocks

          BlockSize=2**i
          NoBlocks=0            !This is the counter for the number of blocks with a given blocksize
          VecInd=1              !Start at the beginning of the list of data to create the blocks
          BlkAv(:)=0.D0         !Rezero the block averages
          BlkAvEn(:)=0.D0         !Rezero the block averages

          do while(VecInd.le.TotPoints)  !Carry on looping over the data until we have 

              Blk=0             !This is the size of the block we are growing
              SumBlock=0.D0     !This is the sum of the elements in the block
              SumBlockEn=0.D0     !This is the sum of the elements in the block
              do while(Blk.lt.BlockSize)
                  SumBlock=SumBlock+Shifts(VecInd)
                  SumBlockEn=SumBlockEn+Energies(VecInd)
                  Blk=Blk+1
                  VecInd=VecInd+1
                  IF(VecInd.gt.TotPoints) EXIT      !Break out of the loop when we reach the end of the list
              enddo

!Check if we have exited because we have reached the end of a block, or if we have reached the end of the datafile
              IF(Blk.eq.Blocksize) THEN
!Block is complete - find the average, and move onto the next block of the data
                  AvBlock=SumBlock/REAL(BlockSize)
                  AvBlockEn=SumBlockEn/REAL(BlockSize)
                  NoBlocks=NoBlocks+1
                  BlkAv(NoBlocks)=AvBlock      !Store the average of each block
                  BlkAvEn(NoBlocks)=AvBlockEn      !Store the average of each block
              ELSE
!We have reached the end of the file - chuck the last few data points for this blocksize
                  IF(VecInd.ne.(TotPoints+1)) THEN
                      STOP 'Should not have exited loop - more data can be read to fill Blk'
                  ENDIF
                  WRITE(6,*) "Number of wasted datapoints for Blocksize ",BlockSize, ", is: ",Blk
                  EXIT      !Exit out of loop over blocks for the given blocksize
              ENDIF

          enddo

!Now we need to calculate relevant information for this blocksize and write it out before moving onto the next blocksize.
          WRITE(6,*) "For a Blocksize of ",BlockSize," , number of blocks created = ",NoBlocks

!First, find the mean & sd of all blocks...
          Mean=0.D0     
          MeanEn=0.D0     
          do j=1,NoBlocks
              Mean=Mean+BlkAv(j)
              MeanSq=MeanSq+(BlkAv(j)**2)
              MeanEn=MeanEn+BlkAvEn(j)
              MeanSqEn=MeanSqEn+(BlkAvEn(j)**2)
          enddo
          Mean=Mean/REAL(NoBlocks)
          MeanSq=MeanSq/REAL(NoBlocks)
          MeanEn=MeanEn/REAL(NoBlocks)
          MeanSqEn=MeanSqEn/REAL(NoBlocks)

          SD=SQRT(MeanSq-(Mean**2))
          Error=SD/SQRT(REAL(NoBlocks-1))

          SDEn=SQRT(MeanSqEn-(MeanEn**2))
          ErrorEn=SDEn/SQRT(REAL(NoBlocks-1))

!This is from the blocking paper, and indicates the error in the blocking error, due to the limited number of blocks available.
          ErrorinError=Error/SQRT(2.D0*(REAL(NoBlocks-1)))

          ErrorinErrorEn=ErrorEn/SQRT(2.D0*(REAL(NoBlocks-1)))

!Write out info
          WRITE(12,'(2I12,4G25.16)') BlockSize,i,Mean,SD,Error,ErrorinError
          
          WRITE(13,'(2I12,4G25.16)') BlockSize,i,MeanEn,SDEn,ErrorEn,ErrorinErrorEn

      enddo
      CLOSE(12)
      CLOSE(13)


END PROGRAM BlkFCIMC
