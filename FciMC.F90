!Choose Tau - find HiiMax
!How to calculate wavevector - add this facility

MODULE FciMCMod
    USE System , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax,
    USE Calc , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
    USE Calc , only : TReadPops,ScaleWalkers
    USE Determinants , only : FDet,GetHElement2
    USE Integrals , only : fck,NMax,nMsh,UMat
    USE Logging , only : TPopsFile,TCalcWavevector,WavevectorPrint
    USE MemoryManager , only : LogMemAlloc,LogMemDealloc
    USE HElem
    IMPLICIT NONE
    SAVE
    TYPE Walker
!Det indicates the determinant that the walker is currently at. WSign: +ve = .true. , -ve = .false.
!Initially, dets are HF and sign is positive
        INTEGER :: Det(NEl)=FDet(1:NEl)
        LOGICAL :: WSign=.true.
    END TYPE Walker

    TYPE(Walker) , POINTER :: WalkVec(:),WalkVec2(:)
    INTEGER :: WalkVecTag=0,WalkVec2Tag=0

!This is the memory in bytes of a single walker
    INTEGER :: WalkerElSize=(NEl+1)*4

    INTEGER Seed,MaxWalkers,TotWalkers,TotWalkersOld,PreviousNMCyc

    TYPE(HElement) :: Hii

    contains

    SUBROUTINE FciMC(Weight,Energyxw)
        IMPLICIT NONE
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: i,iSub

        CALL TISET('FCIMC',iSub)

        IF(HElementSize.gt.1) THEN
            CALL STOPGM("StarDiagMC","StarDiagMC cannot function with complex orbitals.")
        ENDIF

        WRITE(6,*) ""
        WRITE(6,*) "Performing FCIMC...."

        OPEN(15,file='StarMCStats',status='unknown')

        IF(TReadPops)
            WRITE(6,*) "Reading in POPSFILE to restart calculation..."
            OPEN(17,FILE='POPSFILE',Status='old')
            READ(17,*) InitWalkers
            READ(17,*) DiagSft
            READ(17,*) PreviousNMCyc
            WRITE(6,*) "Initial number of walkers read to be: ", InitWalkers
        ELSE
            WRITE(6,*) "Initial number of walkers chosen to be: ", InitWalkers
        ENDIF

        WRITE(6,*) "Damping parameter for Diag Shift set to: ", SftDamp

!Initialise random number seed
        Seed=G_VMC_Seed

        IF(DiagSft.gt.0.D0) THEN
            CALL StopGM("StarDiagMC","Intial value of DiagSft should be negative.")
        ELSE
            WRITE(6,*) "Initial Diagonal Shift (Ecorr guess) is: ", DiagSft
        ENDIF

        CALL InitFCIMCCalc()

!Print out initial starting configurations
        WRITE(6,*) ""
        WRITE(6,*) "       Step  Shift  WalkerChange  GrowRate  TotWalkers"
        WRITE(15,*) "#       Step  Shift  WalkerChange  GrowRate  TotWalkers"
        IF(TReadPops) THEN
            WRITE(6,"(I9,G16.7,I9,G16.7,I9)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
            WRITE(15,"(I9,G16.7,I9,G16.7,I9)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
        ELSE
            WRITE(6,"(I9,G16.7,I9,G16.7,I9)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
            WRITE(15,"(I9,G16.7,I9,G16.7,I9)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
        ENDIF

!Start MC simulation...
        do Iter=1,NMCyc
            CALL PerformFCIMCyc()
        enddo

!This is the heart of FCIMC, where the MC Cycles are performed
    SUBROUTINE PerformFCIMCyc()
        IMPLICIT NONE
        INTEGER :: VecSlot,i,j,k,l 
        
!VecSlot indicates the next free position in WalkVec2
            VecSlot=1

            do j=1,TotWalkers
!j runs through all current walkers

                CALL SetupExcitGen(WalkVec(j)%Det,...)


                IF((WalkVec(j)%Det).eq.1) THEN
!We are at HF - treat this walker slightly differently, since it is attached to all excits
!Run through all double excits and determine whether to create
!Change, so that now we are only selecting a single excitation at random - 
!if we do this, then we should increase prob by No.Excits
!                    r=Ran2(Seed)
!                    r=r*(NList-2)+2
                    do k=2,NList
!                    do k=int(r),int(r)
!Prob of creating a new walker on the excit is given by tau*abs(Hij)
                        rat=tau*abs(List(k,2))

                        IF(rat.gt.Ran2(Seed)) THEN
!Determine sign, and create new walker at k
                            WalkVec2(VecSlot)%Det=k
!If product of signs is +ve, create -ve walker, else create +ve walker
                            IF(WalkVec(j)%WSign) THEN
!Walker is positive
                                IF(List(k,2).gt.0.D0) THEN
                                    WalkVec2(VecSlot)%WSign=.false. !-ve walker
                                ELSE
                                    WalkVec2(VecSlot)%WSign=.true. !+ve walker
                                ENDIF
                            ELSE
!Walker is negative
                                IF(List(k,2).gt.0.D0) THEN
                                    WalkVec2(VecSlot)%WSign=.true. !+ve walker
                                ELSE
                                    WalkVec2(VecSlot)%WSign=.false. !-ve walker
                                ENDIF
                            ENDIF
!Increase walker number, and the next available slot
                            VecSlot=VecSlot+1
                        ENDIF

                    enddo

                ELSE
!We are at an excitation - only possibility is to create walker at HF
                    DetCurr=WalkVec(j)%Det
                    rat=tau*abs(List(DetCurr,2))

                    IF(rat.gt.Ran2(Seed)) THEN
!Create new walker at HF
                        WalkVec2(VecSlot)%Det=1
!If product of signs is +ve, create -ve walker, else create +ve walker
                        IF(WalkVec(j)%WSign) THEN
!Walker is positive
                            IF(List(DetCurr,2).gt.0.D0) THEN
                                WalkVec2(VecSlot)%WSign=.false. !-ve walker
                            ELSE
                                WalkVec2(VecSlot)%WSign=.true. !+ve walker
                            ENDIF
                        ELSE
!Walker is negative
                            IF(List(DetCurr,2).gt.0.D0) THEN
                                WalkVec2(VecSlot)%WSign=.true. !+ve walker
                            ELSE
                                WalkVec2(VecSlot)%WSign=.false. !-ve walker
                            ENDIF
                        ENDIF
!Increase walker number, and the next available slot
                        VecSlot=VecSlot+1
                    ENDIF

!We have finished looking at excitations of walker determinant
                ENDIF

!Next we have to decide if the walker wants to self-destruct or not
!Kill with prob (Hii-DiagSft)*tau - this should ALWAYS be positive
                rat=tau*(List((WalkVec(j)%Det),0)-DiagSft)
                IF(Ran2(Seed).gt.rat) THEN
!This walker is spared - copy him across to the new WalkVec - in the same position
!If the walker isn't copied across, it has self-destructed
                    WalkVec2(VecSlot)%Det=WalkVec(j)%Det
                    WalkVec2(VecSlot)%WSign=WalkVec(j)%WSign
                    VecSlot=VecSlot+1
                ENDIF

!Finish cycling over walkers
            enddo

!Since VecSlot holds the next vacant slot in the array, TotWalkersNew will be one less than this.
            TotWalkersNew=VecSlot-1
            rat=(TotWalkersNew+0.D0)/(MaxWalkers+0.D0)
            IF(rat.gt.0.9) THEN
                WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
            ENDIF

!In this next section, we annihilate pairs of walkers on the same determinant with opposing signs
!This does not have to be performed every cycle - however, we have to at the moment,since it is the only way we transfer the determinants between arrays
!              IF(mod(i,1).eq.0) THEN
            IF(.true.) THEN

                IF(.NOT.TBinCancel) THEN
!This method of canceling down pairs of opposite sign walkers is based on ordering the list of determinants
!Once ordered, each block of walkers on similar determinants can be analysed, and the residual walker concentration moved to WalkVec
!First need to order the vector of new walkers according to determinant the walker is on - this is an NlogN scaling operation
                    CALL SORTIW(TotWalkersNew,WalkVec2(1:TotWalkersNew))

                    j=1
!j is the counter over all TotWalkersNew - it indicates when we have reached the end of the entire list
                    DetCurr=WalkVec2(j)%Det
!DetCurr is the current determinant
                    VecSlot=1
                    do while(j.le.TotWalkersNew)
!Loop over all walkers
                        TotWalkersDet=0
                        do while ((WalkVec2(j)%Det.eq.DetCurr).and.(j.le.TotWalkersNew))
!Loop over all walkers on DetCurr and count residual number after cancelling
                            IF(WalkVec2(j)%WSign) THEN
                                TotWalkersDet=TotWalkersDet+1
                            ELSE
                                TotWalkersDet=TotWalkersDet-1
                            ENDIF
                            j=j+1
                        enddo
!Transfer residual population into VecSlot, along with residual sign
                        IF(TotWalkersDet.gt.0) THEN
                            do l=1,abs(TotWalkersDet)
                                WalkVec(VecSlot)%Det=DetCurr
                                WalkVec(VecSlot)%WSign=.true.
                                VecSlot=VecSlot+1
                            enddo
                        ELSE
                            do l=1,abs(TotWalkersDet)
                                WalkVec(VecSlot)%Det=DetCurr
                                WalkVec(VecSlot)%WSign=.false.
                                VecSlot=VecSlot+1
                            enddo
                        ENDIF
!Update new current determinant
                        DetCurr=WalkVec2(j)%Det
                    enddo
!The new number of cancelled determinants is given by one less than VecSlot again.
                    TotWalkers=VecSlot-1

                    IF(TCalcWavevector.and.(mod(i,WavevectorPrint).eq.0)) THEN
!We want to compare the evolution of wavefunction to exact wavefunction, so calculate it every WavevectorPrint
!Use the redundant List(i,1) vector to store wavefunctions, depicted by the residual concentration of walkers - zero it
                        do j=1,NList
                            List(j,1)=0.D0
                        enddo
                        do j=1,TotWalkers
!Run through all walkers, and bin walkers in the correct wavevector component
                            IF(WalkVec(j)%WSign) THEN
                                List(WalkVec(j)%Det,1)=List(WalkVec(j)%Det,1)+1.D0
                            ELSE
                                List(WalkVec(j)%Det,1)=List(WalkVec(j)%Det,1)-1.D0
                            ENDIF
                        enddo
!Now, normalise wavevector
                        Norm=0.D0
                        do j=1,NList
                            Norm=Norm+(List(j,1)**2)
                        enddo
                        Norm=SQRT(Norm)
                        do j=1,NList
                            List(j,1)=List(j,1)/Norm
                        enddo
                        do j=1,toprint
                            WRITE(14,*) j,List(j,1)
                        enddo
                        WRITE(14,*) ""

                    ENDIF

                ELSE
!In this method of cancelling, we simply histogram the walkers according to the determinant they're on - scales with size of system - i.e. badly
!Use the redundant List(i,1) vector to store wavefunctions, depicted by the residual concentration of walkers - zero it
                    do j=1,NList
                        List(j,1)=0.D0
                    enddo

                    do j=1,TotWalkersNew
!Run through all walkers

                        IF((WalkVec2(j)%Det).gt.NList) THEN
                            WRITE(6,*) "Serious problem here..."
                            STOP 'Serious problem here...'
                        ENDIF

                        IF(WalkVec2(j)%WSign) THEN
!Walker is positive - add to determinant contribution
                            List((WalkVec2(j)%Det),1)=List((WalkVec2(j)%Det),1)+1.D0
                        ELSE
!Walker is negative - subtract from determinant contribution
                            List((WalkVec2(j)%Det),1)=List((WalkVec2(j)%Det),1)-1.D0
                        ENDIF

                    enddo

                    VecSlot=1
                    Norm=0.D0
                    do j=1,NList
!Run through all determinants - normalise, find new number of walkers, and find new WalkVec
                        NWalk=nint(List(j,1))
                        IF(NWalk.ne.0) THEN
!Norm will be used to normalise the eigenvector
                            Norm=Norm+(List(j,1)**2)
                            IF((List(j,1)/abs(List(j,1))).lt.0.D0) THEN
!Component has overall negative sign
                                DetSign=.false.
                            ELSE
                                DetSign=.true.
                            ENDIF
                            do k=1,abs(NWalk)
                                WalkVec(VecSlot)%Det=j
                                WalkVec(VecSlot)%WSign=DetSign
                                VecSlot=VecSlot+1
                            enddo

                        ENDIF

                    enddo

                    IF(TCalcWavevector.and.(mod(i,Wavevectorprint).eq.0)) THEN
!If we want to, we can renormalise the trial vector - only do this every 50 cycles
                        Norm=SQRT(Norm)
                        do j=1,NList
                            List(j,1)=List(j,1)/Norm
                        enddo
                        do j=1,toprint
                            WRITE(14,*) j,List(j,1)
                        enddo
                        WRITE(14,*) ""
                    ENDIF

!The new number of cancelled determinants is given by one less than VecSlot again.
                    TotWalkers=VecSlot-1

                ENDIF

!End of walker cancellation
            ENDIF
!             WRITE(6,"(I9,G16.7,I9,G16.7,I9)") i,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers

!Now we want to calculate the change in the shift of the diagonal elements...
!We only want to change every so often, so that walker numbers have a chance to acclimatise between shift change
!             StepsSft=100
            IF(mod(i,StepsSft).eq.0) THEN
                GrowRate=(TotWalkers+0.D0)/(TotWalkersOld+0.D0)
                DiagSft=DiagSft-(log(GrowRate))/(SftDamp*Tau*(StepsSft+0.D0))
                IF((DiagSft).gt.0.D0) THEN
                    WRITE(6,*) "***WARNING*** - DiagSft trying to become positive..."
                    STOP
                ENDIF
!Write out MC cycle number, Shift, Change in Walker no, Growthrate, New Total Walkers
                IF(TReadPops) THEN
                    WRITE(15,"(I9,G16.7,I9,G16.7,I9)") i+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                    WRITE(6,"(I9,G16.7,I9,G16.7,I9)") i+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                ELSE
                    WRITE(15,"(I9,G16.7,I9,G16.7,I9)") i,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                    WRITE(6,"(I9,G16.7,I9,G16.7,I9)") i,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                ENDIF
                CALL FLUSH(15)
                CALL FLUSH(6)
                TotWalkersOld=TotWalkers
            ENDIF

!Finish MC Cycle
        enddo

        RETURN

    END SUBROUTINE PerformFCIMCyc

!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalc()
        IMPLICIT NONE
        INTEGER :: ierr,i,j,k,l

!Set the maximum number of walkers allowed
        MaxWalkers=100*InitWalkers

!Allocate memory to hold walkers
        ALLOCATE(WalkVec(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec',MaxWalkers,WalkerElSize,this_routine,WalkVecTag,ierr)
        ALLOCATE(WalkVec2(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec2',MaxWalkers,WalkerElSize,this_routine,WalkVec2Tag,ierr)

        IF(TReadPops) THEN
            IF((ABS(ScaleWalkers-1.D0)).lt.1.D-8) THEN
!Read in walker positions
                do i=1,InitWalkers
                    READ(17,*) WalkVec(i)%Det,WalkVec(i)%WSign
                enddo
            ELSE
!Read in walker positions - we will scale these later...
                do i=1,InitWalkers
                    READ(17,*) WalkVec2(i)%Det,WalkVec2(i)%WSign
                enddo
                WRITE(6,*) "Scaling number of walkers by: ",ScaleWalkers
                ReadWalkers=InitWalkers
                InitWalkers=0
!First, count the total number of initial walkers on each determinant - sort into list
                CALL SORTIW(ReadWalkers,WalkVec2(1:ReadWalkers))

                j=1
!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
                DetCurr=WalkVec2(j)%Det
!DetCurr is the current determinant
                do while(j.le.ReadWalkers)
!Loop over all walkers
                    TotWalkersDet=0
                    do while ((WalkVec2(j)%Det.eq.DetCurr).and.(j.le.ReadWalkers))
!Loop over all walkers on DetCurr and count residual number after cancelling
                        IF(WalkVec2(j)%WSign) THEN
                            TotWalkersDet=TotWalkersDet+1
                        ELSE
                            TotWalkersDet=TotWalkersDet-1
                        ENDIF
                        j=j+1
                    enddo
                    DetCurr=WalkVec2(j)%Det
!Count total number of initial walkers
                    InitWalkers=InitWalkers+abs(nint((TotWalkersDet+0.D0)*ScaleWalkers))
                enddo
                WRITE(6,*) "Total number of walkers is now: ",InitWalkers
!Set the new maximum number of walkers allowed
                MaxWalkers=100*InitWalkers

!Deallocate old memory block for WalkVec
                DEALLOCATE(WalkVec)
                CALL LogMemDealloc(this_routine,WalkVecTag)

!Allocate memory to hold new maximum number of walkers
                ALLOCATE(WalkVec(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('WalkVec',MaxWalkers,5,this_routine,WalkVecTag,ierr)

!Now multiply them up...
                j=1
                VecSlot=1
!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
                DetCurr=WalkVec2(j)%Det
!DetCurr is the current determinant
                do while(j.le.ReadWalkers)
!Loop over all walkers
                    TotWalkersDet=0
                    do while ((WalkVec2(j)%Det.eq.DetCurr).and.(j.le.ReadWalkers))
!Loop over all walkers on DetCurr and count residual number after cancelling
                        IF(WalkVec2(j)%WSign) THEN
                            TotWalkersDet=TotWalkersDet+1
                        ELSE
                            TotWalkersDet=TotWalkersDet-1
                        ENDIF
                        j=j+1
                    enddo
!Now multiply up the number of walkers, and insert into WalkVec
                    TotWalkersDet=nint((TotWalkersDet+0.D0)*ScaleWalkers)
                    IF(TotWalkersDet.gt.0) THEN
                        do l=1,abs(TotWalkersDet)
                            WalkVec(VecSlot)%Det=DetCurr
                            WalkVec(VecSlot)%WSign=.true.
                            VecSlot=VecSlot+1
                        enddo
                    ELSE
                        do l=1,abs(TotWalkersDet)
                            WalkVec(VecSlot)%Det=DetCurr
                            WalkVec(VecSlot)%WSign=.false.
                            VecSlot=VecSlot+1
                        enddo
                    ENDIF
                    DetCurr=WalkVec2(j)%Det
                enddo
                IF((VecSlot-1).ne.InitWalkers) THEN
                    WRITE(6,*) "Problem scaling up walker number - exiting..."
                    STOP 'Problem scaling up walker number - exiting...'
                ENDIF

!Now deallocate and reallocate WalkVec2 with correct number of total walkers
                DEALLOCATE(WalkVec2)
                CALL LogMemDealloc(this_routine,WalkVec2Tag)
                ALLOCATE(WalkVec2(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('WalkVec2',MaxWalkers,5,this_routine,WalkVec2Tag,ierr)

            ENDIF

!End of reading in POPSFILE
            CLOSE(17)
        ENDIF
!Setup initial trial walker position - already they are all set to be at HF to start (with positive sign)
!TotWalkers contains the number of current walkers at each step
        TotWalkers=InitWalkers
        TotWalkersOld=InitWalkers

        RETURN

    END SUBROUTINE InitFCIMCCalc






