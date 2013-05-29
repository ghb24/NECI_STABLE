PROGRAM BlkFCIMC

      IMPLICIT NONE
      INTEGER, PARAMETER :: MaxFields=100
      CHARACTER(Len=1000) :: FirstLineRubbish,Line,data1(MaxFields)
      INTEGER :: i,Iter,WalkCng,TotWalkers,Annihil,Died,Born,StartIter,ierr,TotPoints
      INTEGER :: TotBlkSize,NoBlocks,BlockSize,VecInd,Blk,j,k,fields,err
      REAL*8 :: Sft,GrowRate,ProjE,AvBlock,SumBlock,Mean,MeanSq,SD,Error,ErrorinError,AvShift,ProjEInst
      REAL*8 :: AvBlockEn,SumBlockEn,MeanEn,MeanSqEn,SDEn,ErrorEn,ErrorinErrorEn
      REAL*8 :: AvBlockHF,SumBlockHF,MeanHF,MeanSqHF,SDHF,ErrorHF,ErrorinErrorHF
      REAL*8 :: AvBlockNum,SumBlockNum,MeanNum,MeanSqNum,SDNum,ErrorNum,ErrorinErrorNum
      REAL*8 :: MaxDeviation,MaxDeviationEn,CorrCoeff
      REAL*8 :: MaxDeviationNum,MaxDeviationHF,FinalVal,TrueMean,TrueMeanEn,TrueMeanHF,TrueMeanNum
      REAL*8 :: CrossCorr,Covar,FinalErr,LargestErr,LargestErrHF,LargestErrNum,LargestErrEn
      REAL*8 , ALLOCATABLE :: Shifts(:),Energies(:),BlkAv(:),BlkAvEn(:),Num(:),HF(:),BlkAvNum(:),BlkAvHF(:)
      LOGICAL :: tSeperateBlock

      INTEGER :: AllNoatHF,AllNoatDoubs,WalkersDiffProc
      REAL*8 :: AccRat,IterTime,FracSing,TotImagTime,IterEnergy,HFShift,InstShift,TotEnergy,EnergyNum,EnergyHF
      INTEGER*8 :: AllTotWalkers

      OPEN(UNIT=11,FILE='FCIMCStats',STATUS='OLD',ACTION='READ',    &
        POSITION='REWIND',FORM='FORMATTED')
      OPEN(UNIT=12,FILE='BlockingInfo',STATUS='UNKNOWN')
      OPEN(UNIT=13,FILE='BlockingInfoEnergies',STATUS='UNKNOWN')
      WRITE(12,"(A)") '#  Blocklength   Log_2 BL     MeanBlocks              &
      &SD-Blocks               Err-Blocks           Err in Err-Blocks'
      WRITE(13,"(A)") '#  Blocklength   Log_2 BL     MeanBlocks              &
      &SD-Blocks               Err-Blocks           Err in Err-Blocks'

      READ(11,*) FirstLineRubbish
      READ(11,*) FirstLineRubbish
      READ(11,*) FirstLineRubbish

!Cunning trick to calculate the number of fields in the FCIMCStats file, to tell us how up-to-date it is
      READ(11,"(A)") Line
!      WRITE(6,*) Line
      do i=MaxFields,1,-1
          READ(Line,*,iostat=err) (data1(k),k=1,i)
          IF(err.eq.0) THEN
              fields=i
              EXIT
          ENDIF
      enddo
      IF(i.le.1) STOP 'Error in counting number of fields'

      WRITE(6,*) "Number of fields = ",fields

      IF(fields.lt.25) THEN
          WRITE(6,*) "OLD FCIMCStats file detected. Cannot do seperate blocking. Energy may represent a biased estimate"
          tSeperateBlock=.false.
      ELSE
          WRITE(6,*) "Seperate numerator and denominator contributions to the projected energy will be blocked."
          WRITE(6,*) "BlockingInfoNum and BlockingInfoHF should be seperately analysed to calculate true energy and error."
          WRITE(6,*) "BlockingInfoEnergies may represent a biased estimator."
          tSeperateBlock=.true.
      ENDIF

      IF(tSeperateBlock) THEN
          OPEN(UNIT=14,FILE='BlockingInfoNum',STATUS='UNKNOWN')
          OPEN(UNIT=15,FILE='BlockingInfoHF',STATUS='UNKNOWN')
      ENDIF

      REWIND(11)
      READ(11,*) FirstLineRubbish
      READ(11,*) FirstLineRubbish
      READ(11,*) FirstLineRubbish
    
      WRITE(6,*) "Enter the iteration number from which you want to start the blocking analysis:"
      READ(*,*) StartIter

      i=0
      do while(.true.)
          READ(11,'(I12,G16.7,I10,G16.7,I12,3I13,3G17.9)',END=99) Iter,Sft,WalkCng,GrowRate,TotWalkers,&
          Annihil,Died,Born,ProjE,AvShift,ProjEInst
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
      IF(tSeperateBlock) THEN
          ALLOCATE(Num(TotPoints),stat=ierr)
          IF(ierr.ne.0) THEN
              STOP 'Error in allocation'
          ENDIF
          ALLOCATE(HF(TotPoints),stat=ierr)
          IF(ierr.ne.0) THEN
              STOP 'Error in allocation'
          ENDIF
          Num(:)=0.D0
          HF(:)=0.D0
      ENDIF
      Shifts(:)=0.D0
      Energies(:)=0.D0

!Reread in data to fill arrays
      REWIND(11)
      READ(11,*) FirstLineRubbish
      READ(11,*) FirstLineRubbish
      READ(11,*) FirstLineRubbish
      i=0
      do while(.true.)

          IF(tSeperateBlock) THEN
              READ(11,"(I12,G16.7,I10,G16.7,I12,3I13,3G17.9,2I10,G13.5,I12,G13.5,G17.5,I13,G13.5,6G17.9)",END=98) Iter,Sft,WalkCng,GrowRate,   &
                   TotWalkers,Annihil,Died,Born,ProjE,AvShift,ProjEInst,AllNoatHF,AllNoatDoubs,AccRat,AllTotWalkers,IterTime,   &
                   FracSing,WalkersDiffProc,TotImagTime,IterEnergy,HFShift,InstShift,TotEnergy,EnergyHF,EnergyNum
          ELSE
              READ(11,'(I12,G16.7,I10,G16.7,I12,3I13,3G17.9)',END=98) Iter,Sft,WalkCng,GrowRate,TotWalkers,&
              Annihil,Died,Born,ProjE,AvShift,ProjEInst
              EnergyNum=0.D0
              EnergyHF=0.D0
          ENDIF
          IF(Iter.ge.StartIter) THEN
              i=i+1
              Shifts(i)=Sft
              Energies(i)=ProjEInst
              IF(tSeperateBlock) THEN
                  Num(i)=EnergyNum
                  HF(i)=EnergyHF
              ENDIF
          ENDIF
      enddo

98    CONTINUE
      CLOSE(11)

      Mean=0.D0
      MeanEn=0.D0
      MeanNum=0.D0
      MeanHF=0.D0
      do i=1,TotPoints
          Mean=Mean+Shifts(i)
          MeanEn=MeanEn+Energies(i)
		  IF(tSeperateBlock) THEN
			  MeanNum=MeanNum+Num(i)
			  MeanHF=MeanHF+HF(i)
		  ENDIF
      enddo
      Mean=Mean/REAL(TotPoints,8)
      MeanEn=MeanEn/REAL(TotPoints,8)
      MeanNum=MeanNum/REAL(TotPoints,8)
      MeanHF=MeanHF/REAL(TotPoints,8)

      MaxDeviation=0.D0
      MaxDeviationEn=0.D0
      MaxDeviationNum=0.D0
      MaxDeviationHF=0.D0
      do i=1,TotPoints
          IF(ABS(Shifts(i)-Mean).gt.MaxDeviation) MaxDeviation=ABS(Shifts(i)-Mean)
          IF(ABS(Energies(i)-MeanEn).gt.MaxDeviationEn) MaxDeviationEn=ABS(Energies(i)-MeanEn)
		  IF(tSeperateBlock) THEN
			  IF(ABS(Num(i)-MeanNum).gt.MaxDeviationNum) MaxDeviationNum=ABS(Num(i)-MeanNum)
			  IF(ABS(HF(i)-MeanHF).gt.MaxDeviationHF) MaxDeviationHF=ABS(HF(i)-MeanHF)
		  ENDIF
      enddo

      WRITE(6,"(A,4G25.16)") "Mean values for the shift, projected energy, numerator for projected E and HF are: ",Mean,MeanEn,MeanNum,MeanHF
      WRITE(6,"(A,4G25.16)") "Maximum deviations from the mean shift, projected energy, numerator for proj E and HF are: ",MaxDeviation,MaxDeviationEn,MaxDeviationNum,MaxDeviationHF


!Allocate data to hold averages of each block
      ALLOCATE(BlkAv(TotPoints),stat=ierr)
      ALLOCATE(BlkAvEn(TotPoints),stat=ierr)
      ALLOCATE(BlkAvNum(TotPoints),stat=ierr)
      ALLOCATE(BlkAvHF(TotPoints),stat=ierr)
      IF(ierr.ne.0) THEN
          STOP 'Error in allocation'
      ENDIF
      BlkAv(:)=0.D0
      BlkAvEn(:)=0.D0
      BlkAvNum(:)=0.D0
      BlkAvHF(:)=0.D0

      LargestErr=0.D0
      LargestErrEn=0.D0
      LargestErrNum=0.D0
      LargestErrHF=0.D0

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
          NoBlocks=0              !This is the counter for the number of blocks with a given blocksize
          VecInd=1                !Start at the beginning of the list of data to create the blocks
          BlkAv(:)=0.D0           !Rezero the block averages
          BlkAvEn(:)=0.D0         !Rezero the block averages
          BlkAvNum(:)=0.D0         !Rezero the block averages
          BlkAvHF(:)=0.D0         !Rezero the block averages

          do while(VecInd.le.TotPoints)  !Carry on looping over the data until we have 

              Blk=0             !This is the size of the block we are growing
              SumBlock=0.D0     !This is the sum of the elements in the block
              SumBlockEn=0.D0     !This is the sum of the elements in the block
              SumBlockNum=0.D0
              SumBlockHF=0.D0
              do while(Blk.lt.BlockSize)
                  SumBlock=SumBlock+Shifts(VecInd)
                  SumBlockEn=SumBlockEn+Energies(VecInd)
				  IF(tSeperateBlock) THEN
					  SumBlockNum=SumBlockNum+Num(VecInd)
					  SumBlockHF=SumBlockHF+HF(VecInd)
				  ENDIF
                  Blk=Blk+1
                  VecInd=VecInd+1
                  IF(VecInd.gt.TotPoints) EXIT      !Break out of the loop when we reach the end of the list
              enddo

!Check if we have exited because we have reached the end of a block, or if we have reached the end of the datafile
              IF(Blk.eq.Blocksize) THEN
!Block is complete - find the average, and move onto the next block of the data
                  AvBlock=SumBlock/REAL(BlockSize)
                  AvBlockEn=SumBlockEn/REAL(BlockSize)
                  AvBlockNum=SumBlockNum/REAL(BlockSize)
                  AvBlockHF=SumBlockHF/REAL(BlockSize)
                  NoBlocks=NoBlocks+1
                  BlkAv(NoBlocks)=AvBlock      !Store the average of each block
                  BlkAvEn(NoBlocks)=AvBlockEn      !Store the average of each block
                  BlkAvNum(NoBlocks)=AvBlockNum
                  BlkAvHF(NoBlocks)=AvBlockHF
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
          MeanSq=0.D0     
          MeanEn=0.D0     
          MeanSqEn=0.D0
          MeanNum=0.D0
          MeanSqNum=0.D0
          MeanHF=0.D0
          MeanSqHF=0.D0
          do j=1,NoBlocks
!              WRITE(6,*) BlkAvEn(j)
              Mean=Mean+BlkAv(j)
              MeanSq=MeanSq+(BlkAv(j)**2)
              MeanEn=MeanEn+BlkAvEn(j)
              MeanSqEn=MeanSqEn+(BlkAvEn(j)**2)
              MeanNum=MeanNum+BlkAvNum(j)
              MeanSqNum=MeanSqNum+(BlkAvNum(j)**2)
              MeanHF=MeanHF+BlkAvHF(j)
              MeanSqHF=MeanSqHF+(BlkAvHF(j)**2)
          enddo
          Mean=Mean/REAL(NoBlocks)
          MeanSq=MeanSq/REAL(NoBlocks)
          MeanEn=MeanEn/REAL(NoBlocks)
          MeanSqEn=MeanSqEn/REAL(NoBlocks)
          MeanNum=MeanNum/REAL(NoBlocks)
          MeanSqNum=MeanSqNum/REAL(NoBlocks)
          MeanHF=MeanHF/REAL(NoBlocks)
          MeanSqHF=MeanSqHF/REAL(NoBlocks)

          SD=SQRT(MeanSq-(Mean**2))
          Error=SD/SQRT(REAL(NoBlocks-1))

          SDEn=SQRT(MeanSqEn-(MeanEn**2))
          ErrorEn=SDEn/SQRT(REAL(NoBlocks-1))

          SDNum=SQRT(MeanSqNum-(MeanNum**2))
          ErrorNum=SDNum/SQRT(REAL(NoBlocks-1))
          
          SDHF=SQRT(MeanSqHF-(MeanHF**2))
          ErrorHF=SDHF/SQRT(REAL(NoBlocks-1))

!This is from the blocking paper, and indicates the error in the blocking error, due to the limited number of blocks available.
          ErrorinError=Error/SQRT(2.D0*(REAL(NoBlocks-1)))

          ErrorinErrorEn=ErrorEn/SQRT(2.D0*(REAL(NoBlocks-1)))
          ErrorinErrorNum=ErrorNum/SQRT(2.D0*(REAL(NoBlocks-1)))
          ErrorinErrorHF=ErrorHF/SQRT(2.D0*(REAL(NoBlocks-1)))

          !Find the largest errors, neglecting the three largest block sizes
          IF(i.lt.TotBlkSize-2) THEN
              IF(Error.gt.LargestErr) LargestErr=Error
              IF(ErrorEn.gt.LargestErrEn) LargestErrEn=ErrorEn
              IF(ErrorNum.gt.LargestErrNum) LargestErrNum=ErrorNum
              IF(ErrorHF.gt.LargestErrHF) LargestErrHF=ErrorHF
          ENDIF

          !Save the true mean
          IF(i.eq.0) THEN
              TrueMean=Mean
              TrueMeanEn=MeanEn
              TrueMeanNum=MeanNum
              TrueMeanHF=MeanHF
          ENDIF
              
!Write out info
          WRITE(12,'(2I12,4G25.16)') BlockSize,i,Mean,SD,Error,ErrorinError
          WRITE(13,'(2I12,4G25.16)') BlockSize,i,MeanEn,SDEn,ErrorEn,ErrorinErrorEn
          IF(tSeperateBlock) THEN
              WRITE(14,'(2I12,4G25.16)') BlockSize,i,MeanNum,SDNum,ErrorNum,ErrorinErrorNum
              WRITE(15,'(2I12,4G25.16)') BlockSize,i,MeanHF,SDHF,ErrorHF,ErrorinErrorHF
          ENDIF

      enddo
      CLOSE(12)
      CLOSE(13)
      IF(tSeperateBlock) THEN
          CLOSE(14)
          CLOSE(15)

          Covar=0.D0
          MeanSqHF=0.D0
          MeanSqNum=0.D0
          do i=1,TotPoints
              Covar=Covar+((Num(i)-TrueMeanNum)*(HF(i)-TrueMeanHF))
              MeanSqHF=MeanSqHF+(HF(i)**2)
              MeanSqNum=MeanSqNum+(Num(i)**2)
          enddo
          Covar=Covar/REAL(TotPoints,8)
          MeanSqHF=MeanSqHF/REAL(TotPoints,8)
          MeanSqNum=MeanSqNum/REAL(TotPoints,8)

          SDHF=SQRT(MeanSqHF-(TrueMeanHF**2))
          SDNum=SQRT(MeanSqNum-(TrueMeanNum**2))

          CorrCoeff=Covar/(SDHF*SDNum)

          WRITE(6,*) "Correlation Coefficient is: ",CorrCoeff
          
          FinalVal=TrueMeanNum/TrueMeanHF
          CrossCorr=(2.D0*CorrCoeff*LargestErrNum*LargestErrHF)/(TrueMeanHF*TrueMeanNum)
          FinalErr=ABS(FinalVal)*SQRT(((LargestErrHF/TrueMeanHF)**2)+((LargestErrNum/TrueMeanNum)**2)-CrossCorr)

!          WRITE(6,*) "Relative err num: ",((LargestErrNum/TrueMeanNum)**2)
!          WRITE(6,*) "Relative err HF: ",((LargestErrHF/TrueMeanHF)**2)
!          WRITE(6,*) "Cross Corr: ",CrossCorr

          WRITE(6,"(A,G20.10,A,G20.10,A,G20.10,A)") "Assuming that error in the HF values is: ",LargestErrHF, " and in the energy numerator: ",LargestErrNum, " and in shift: ",LargestErr," we can estimate the error to be:"

          WRITE(6,"(A,G25.15,A,G25.15)") "***ENERGY-PROJ*** ",FinalVal, " +/- ", FinalErr
          WRITE(6,"(A,G25.15,A,G25.15)") "***ENERGY-SHIFT***",TrueMean, " +/- ",LargestErr
          CALL FLUSH(6)

      ENDIF


END PROGRAM BlkFCIMC
