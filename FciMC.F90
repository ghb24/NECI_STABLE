MODULE FciMCMod
    USE System , only : NEl,Alat,Brr,ECore,G1,nBasis,nBasisMax
    USE Calc , only : InitWalkers,NMCyc,G_VMC_Seed,DiagSft,Tau,SftDamp,StepsSft
    USE Calc , only : TReadPops,ScaleWalkers,TMCExcitSpace,NoMCExcits
    USE Determinants , only : FDet,GetHElement2
    USE Integrals , only : fck,NMax,nMsh,UMat
    USE Logging , only : TPopsFile,TCalcWavevector,WavevectorPrint
    USE MemoryManager , only : LogMemAlloc,LogMemDealloc
    USE HElem
    IMPLICIT NONE
    SAVE

    INTEGER , POINTER :: WalkVecDets(:,:),WalkVec2Dets(:,:)
    LOGICAL , POINTER :: WalkVecSign(:),WalkVec2Sign(:)
    INTEGER :: WalkVecDetsTag=0,WalkVec2DetsTag=0,WalkVecSignTag=0,WalkVec2SignTag=0

    INTEGER :: Seed,MaxWalkers,TotWalkers,TotWalkersOld,PreviousNMCyc,Iter
    INTEGER :: exFlag=3

    REAL*8 :: GrowRate

    TYPE(HElement) :: Hii

    contains

    SUBROUTINE FciMC(Weight,Energyxw)
        IMPLICIT NONE
        TYPE(HDElement) :: Weight,Energyxw
        INTEGER :: i,iSub
        CHARACTER(len=*), PARAMETER :: this_routine='FCIMC'

        CALL TISET('FCIMC',iSub)

        IF(HElementSize.gt.1) THEN
            CALL STOPGM("StarDiagMC","StarDiagMC cannot function with complex orbitals.")
        ENDIF

        WRITE(6,*) ""
        WRITE(6,*) "Performing FCIMC...."

        IF(TCalcWaveVector) THEN
            WRITE(6,*) "Wavevector calculation is only available in star MC..."
        ENDIF

        OPEN(15,file='FCIMCStats',status='unknown')

        IF(TReadPops) THEN
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

!Calculate Hii
        Hii=GetHElement2(FDet,FDet,NEl,nBasisMax,G1,nBasis,Brr,nMsh,fck,NMax,ALat,UMat,0,ECore)

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
!TotWalkersOld is the number of walkers last time the shift was changed
        IF(TReadPops) THEN
            WRITE(6,"(I9,G16.7,I9,G16.7,I9)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
            WRITE(15,"(I9,G16.7,I9,G16.7,I9)") PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
        ELSE
            WRITE(6,"(I9,G16.7,I9,G16.7,I9)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
            WRITE(15,"(I9,G16.7,I9,G16.7,I9)") 0,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
        ENDIF
        
!Reset TotWalkersOld so that it is the number of walkers now
        TotWalkersOld=TotWalkers

!Start MC simulation...
        do Iter=1,NMCyc
            
            CALL PerformFCIMCyc()
            
            IF(mod(Iter,StepsSft).eq.0) THEN
!Every StepsSft steps, update the diagonal shift value (the running value for the correlation energy)
!We don't want to do this too often, since we want the population levels to acclimatise between changing the shifts
                CALL UpdateDiagSft()

!Write out MC cycle number, Shift, Change in Walker no, Growthrate, New Total Walkers
                IF(TReadPops) THEN
                    WRITE(15,"(I9,G16.7,I9,G16.7,I9)") Iter+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                    WRITE(6,"(I9,G16.7,I9,G16.7,I9)") Iter+PreviousNMCyc,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                ELSE
                    WRITE(15,"(I9,G16.7,I9,G16.7,I9)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                    WRITE(6,"(I9,G16.7,I9,G16.7,I9)") Iter,DiagSft,TotWalkers-TotWalkersOld,GrowRate,TotWalkers
                ENDIF
                CALL FLUSH(15)
                CALL FLUSH(6)

!Reset TotWalkersOld so that it is the number of walkers now
                TotWalkersOld=TotWalkers
            ENDIF

!End of MC cycle
        enddo

        IF(TPopsFile) THEN
!Print out current state of simulation, so it can be restarted if desired...
            OPEN(16,file='POPSFILE',status='unknown')
            WRITE(16,*) TotWalkers, "   TOTWALKERS"
            WRITE(16,*) DiagSft, "   DIAG SHIFT"
            IF(TReadPops) THEN
                WRITE(16,*) NMCyc+PreviousNMCyc, "   NO. CYCLES"
            ELSE
                WRITE(16,*) NMCyc, "   MC CYCLES"
            ENDIF
            do i=1,TotWalkers
                WRITE(16,*) WalkVecDets(:,i),WalkVecSign(i)
            enddo
            CLOSE(16)
        ENDIF

        Weight=HDElement(1.D0)
        Energyxw=HDElement(2.D0*DiagSft)

!Deallocate memory
        DEALLOCATE(WalkVecDets)
        CALL LogMemDealloc(this_routine,WalkVecDetsTag)
        DEALLOCATE(WalkVec2Dets)
        CALL LogMemDealloc(this_routine,WalkVec2DetsTag)
        DEALLOCATE(WalkVecSign)
        CALL LogMemDealloc(this_routine,WalkVecSignTag)
        DEALLOCATE(WalkVec2Sign)
        CALL LogMemDealloc(this_routine,WalkVec2SignTag)

        CLOSE(15)

        CALL TIHALT('FCIMC',iSub)

        RETURN

    END SUBROUTINE FciMC

!This is the heart of FCIMC, where the MC Cycles are performed
    SUBROUTINE PerformFCIMCyc()
        USE System , only : Arr
        IMPLICIT NONE
        INTEGER :: VecSlot,i,j,k,l,DetCurr(NEl),iMaxExcit,nExcitMemLen,nStore(6)
        INTEGER :: nJ(NEl),ierr,nExcitTag=0,IC,Child,iSubCyc,TotWalkersNew
        REAL*8 :: Prob,rat
        INTEGER , ALLOCATABLE :: nExcit(:)
        CHARACTER(len=*), PARAMETER :: this_routine='PerformFCIMCyc'
        
        CALL TISET('MCyc',iSubCyc)
        
!VecSlot indicates the next free position in WalkVec2
        VecSlot=1

        do j=1,TotWalkers
!j runs through all current walkers
            do k=1,NEl
                DetCurr(k)=WalkVecDets(k,j)
            enddo

!Setup excit generators for this determinant (This can be reduced to an order N routine later for abelian symmetry.
            iMaxExcit=0
            CALL IAZZERO(nStore,6)
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,exFlag)
            ALLOCATE(nExcit(nExcitMemLen),stat=ierr)
            CALL LogMemAlloc('nExcit',nExcitMemLen,4,this_routine,nExcitTag,ierr)
            nExcit(1)=0
            CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcit,nJ,iMaxExcit,0,nStore,exFlag)

            IF(TMCExcitSpace) THEN
!Excitations are picked stochastically

                do i=1,NoMCExcits

                    CALL GenRandSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,nExcit,nJ,Seed,IC,0,UMat,Arr,Prob)

                    Child=AttemptCreate(DetCurr,WalkVecSign(j),nJ,Prob,IC)
                    IF(Child.eq.1) THEN
!We have successfully created a positive child at nJ
                        do k=1,NEl
                            WalkVec2Dets(k,VecSlot)=nJ(k)
                        enddo
                        WalkVec2Sign(VecSlot)=.true.
                        VecSlot=VecSlot+1
                    ELSEIF(Child.eq.-1) THEN
!We have successfully created a negative child at nJ
                        do k=1,NEl
                            WalkVec2Dets(k,VecSlot)=nJ(k)
                        enddo
                        WalkVec2Sign(VecSlot)=.false.
                        VecSlot=VecSlot+1
                    ENDIF

                enddo

            ELSE
!Run through all possible excitations of each walker

                do while(.true.)
                    CALL GenSymExcitIt2(DetCurr,NEl,G1,nBasis,nBasisMax,.FALSE.,nExcit,nJ,IC,0,nStore,exFlag)
                    IF(nJ(1).eq.0) EXIT

                    Child=AttemptCreate(DetCurr,WalkVecSign(j),nJ,1.D0,IC)
                    IF(Child.eq.1) THEN
!We have successfully created a positive child at nJ
                        do k=1,NEl
                            WalkVec2Dets(k,VecSlot)=nJ(k)
                        enddo
                        WalkVec2Sign(VecSlot)=.true.
                        VecSlot=VecSlot+1
                    ELSEIF(Child.eq.-1) THEN
!We have successfully created a negative child at nJ
                        do k=1,NEl
                            WalkVec2Dets(k,VecSlot)=nJ(k)
                        enddo
                        WalkVec2Sign(VecSlot)=.false.
                        VecSlot=VecSlot+1
                    ENDIF

                enddo

            ENDIF

!We now have to decide whether the parent particle (j) wants to self-destruct or not...
            IF(.NOT.AttemptDie(DetCurr)) THEN
!This indicates that the particle is spared...copy him across to WalkVec2
                do k=1,NEl
                    WalkVec2Dets(k,VecSlot)=DetCurr(k)
                enddo
                WalkVec2Sign(VecSlot)=WalkVecSign(j)
                VecSlot=VecSlot+1
            ENDIF

!Destroy excitation generators for current walker
            DEALLOCATE(nExcit)
            CALL LogMemDealloc(this_routine,nExcitTag)

!Finish cycling over walkers
        enddo
                
!Since VecSlot holds the next vacant slot in the array, TotWalkersNew will be one less than this.
        TotWalkersNew=VecSlot-1
        rat=(TotWalkersNew+0.D0)/(MaxWalkers+0.D0)
        IF(rat.gt.0.9) THEN
            WRITE(6,*) "*WARNING* - Number of walkers has increased to over 90% of MaxWalkers"
        ENDIF

!This routine now cancels down the particles with opposing sign on each determinant
!This routine does not necessarily need to be called every Iter, but it does at the moment, since it is the only way to 
!transfer the residual particles back onto WalkVec
        CALL AnnihilatePairs(TotWalkersNew)
        
        CALL TIHALT('MCyc',iSubCyc)

        RETURN

    END SUBROUTINE PerformFCIMCyc

!This routine looks at the change in residual particle number over a number of cycles, and adjusts the 
!value of the diagonal shift in the hamiltonian in order to compensate for this
    SUBROUTINE UpdateDiagSft()
        IMPLICIT NONE

        GrowRate=(TotWalkers+0.D0)/(TotWalkersOld+0.D0)
        DiagSft=DiagSft-(log(GrowRate))/(SftDamp*Tau*(StepsSft+0.D0))
        IF((DiagSft).gt.0.D0) THEN
            WRITE(6,*) "***WARNING*** - DiagSft trying to become positive..."
            STOP
        ENDIF

    END SUBROUTINE UpdateDiagSft


!This routine cancels out particles of opposing sign on the same determinant.
    SUBROUTINE AnnihilatePairs(TotWalkersNew)
        IMPLICIT NONE
        INTEGER :: TotWalkersNew,j,k,l,DetCurr(NEl),VecSlot,TotWalkersDet
        INTEGER :: DetLT

!First, it is necessary to sort the list of determinants
        CALL SortDets(TotWalkersNew,WalkVec2Dets(:,1:TotWalkersNew),NEl,WalkVec2Sign(1:TotWalkersNew),1)

!Once ordered, each block of walkers on similar determinants can be analysed, and the residual walker concentration moved to WalkVec
        j=1
!j is the counter over all uncancelled walkers - it indicates when we have reached the end of the list of total walkers
        do k=1,NEl
!DetCurr is the current determinant
            DetCurr(k)=WalkVec2Dets(k,j)
        enddo
        VecSlot=1

        do while(j.le.TotWalkersNew)
!Loop over all walkers
            TotWalkersDet=0
            do while ((DetLT(WalkVec2Dets(:,j),DetCurr,NEl).eq.0).and.(j.le.TotWalkersNew))
!Loop over all walkers on DetCurr and count residual number after cancelling
                IF(WalkVec2Sign(j)) THEN
                    TotWalkersDet=TotWalkersDet+1
                ELSE
                    TotWalkersDet=TotWalkersDet-1
                ENDIF
                j=j+1
            enddo
!Transfer residual population into VecSlot, along with residual sign
            IF(TotWalkersDet.gt.0) THEN
!Positive sign particles want to populate this determinant
                do l=1,abs(TotWalkersDet)
                    do k=1,NEl
                        WalkVecDets(k,VecSlot)=DetCurr(k)
                    enddo
                    WalkVecSign(VecSlot)=.true.
                    VecSlot=VecSlot+1
                enddo
            ELSE
!Negative sign particles want to populate this determinant
                do l=1,abs(TotWalkersDet)
                    do k=1,NEl
                        WalkVecDets(k,VecSlot)=DetCurr(k)
                    enddo
                    WalkVecSign(VecSlot)=.false.
                    VecSlot=VecSlot+1
                enddo
            ENDIF
!Now update the current determinant
            do k=1,NEl
                DetCurr(k)=WalkVec2Dets(k,j)
            enddo
        enddo
!The new number of residual cancelled walkers is given by one less that VecSlot again.
        TotWalkers=VecSlot-1

        RETURN

    END SUBROUTINE AnnihilatePairs

!This initialises the calculation, by allocating memory, setting up the initial walkers, and reading from a file if needed
    SUBROUTINE InitFCIMCCalc()
        IMPLICIT NONE
        INTEGER :: ierr,i,j,k,l,DetCurr(NEl),ReadWalkers,TotWalkersDet
        INTEGER :: DetLT,VecSlot
        CHARACTER(len=*), PARAMETER :: this_routine='InitFCIMC'

!Set the maximum number of walkers allowed
        MaxWalkers=100*InitWalkers

!Allocate memory to hold walkers
        ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
        ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
        ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)
        ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
        CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)

        IF(TReadPops) THEN
            IF((ABS(ScaleWalkers-1.D0)).lt.1.D-8) THEN
!Read in walker positions
                do i=1,InitWalkers
                    READ(17,*) WalkVecDets(:,i),WalkVecSign(i)
                enddo
            ELSE
!Read in walker positions - we will scale these later...
                do i=1,InitWalkers
                    READ(17,*) WalkVec2Dets(:,i),WalkVec2Sign(i)
                enddo
                WRITE(6,*) "Scaling number of walkers by: ",ScaleWalkers
                ReadWalkers=InitWalkers
                InitWalkers=0
!First, count the total number of initial walkers on each determinant - sort into list
                CALL SortDets(ReadWalkers,WalkVec2Dets(:,1:ReadWalkers),NEl,WalkVec2Sign(1:ReadWalkers),1)

                j=1
!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
                do k=1,NEl
!DetCurr is the current determinant
                    DetCurr(k)=WalkVec2Dets(k,j)
                enddo

                do while(j.le.ReadWalkers)
!Loop over all walkers
                    TotWalkersDet=0
                    do while ((DetLT(WalkVec2Dets(:,j),DetCurr,NEl).eq.0).and.(j.le.ReadWalkers))
!Loop over all walkers on DetCurr and count residual number after cancelling
                        IF(WalkVec2Sign(j)) THEN
                            TotWalkersDet=TotWalkersDet+1
                        ELSE
                            TotWalkersDet=TotWalkersDet-1
                        ENDIF
                        j=j+1
                    enddo
!Now update the current determinant
                    do k=1,NEl
                        DetCurr(k)=WalkVec2Dets(k,j)
                    enddo
!Count total number of initial walkers
                    InitWalkers=InitWalkers+abs(nint((TotWalkersDet+0.D0)*ScaleWalkers))
                enddo
                WRITE(6,*) "Total number of walkers is now: ",InitWalkers
!Set the new maximum number of walkers allowed
                MaxWalkers=100*InitWalkers

!Deallocate old memory block for WalkVec
                DEALLOCATE(WalkVecDets)
                CALL LogMemDealloc(this_routine,WalkVecDetsTag)
                DEALLOCATE(WalkVecSign)
                CALL LogMemDealloc(this_routine,WalkVecSignTag)

!Allocate memory to hold new maximum number of walkers
                ALLOCATE(WalkVecDets(NEl,MaxWalkers),stat=ierr)
                CALL LogMemAlloc('WalkVecDets',MaxWalkers*NEl,4,this_routine,WalkVecDetsTag,ierr)
                ALLOCATE(WalkVecSign(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('WalkVecSign',MaxWalkers,4,this_routine,WalkVecSignTag,ierr)

!Now multiply them up...
                j=1
                VecSlot=1
!j is the counter over all read in walkers - it indicates when we have reached the end of the entire list
                do k=1,NEl
                    DetCurr(k)=WalkVec2Dets(k,j)
                enddo
!DetCurr is the current determinant
                do while(j.le.ReadWalkers)
!Loop over all walkers
                    TotWalkersDet=0
                    do while ((DetLT(WalkVec2Dets(:,j),DetCurr,NEl).eq.0).and.(j.le.ReadWalkers))
!Loop over all walkers on DetCurr and count residual number after cancelling
                        IF(WalkVec2Sign(j)) THEN
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
                            do k=1,NEl
                                WalkVecDets(k,VecSlot)=DetCurr(k)
                            enddo
                            WalkVecSign(VecSlot)=.true.
                            VecSlot=VecSlot+1
                        enddo
                    ELSE
                        do l=1,abs(TotWalkersDet)
                            do k=1,NEl
                                WalkVecDets(k,VecSlot)=DetCurr(k)
                            enddo
                            WalkVecSign(VecSlot)=.false.
                            VecSlot=VecSlot+1
                        enddo
                    ENDIF
                    do k=1,NEl
                        DetCurr(k)=WalkVec2Dets(k,j)
                    enddo
                enddo
                IF((VecSlot-1).ne.InitWalkers) THEN
                    WRITE(6,*) "Problem scaling up walker number - exiting..."
                    STOP 'Problem scaling up walker number - exiting...'
                ENDIF

!Now deallocate and reallocate WalkVec2 with correct number of total walkers
                DEALLOCATE(WalkVec2Dets)
                CALL LogMemDealloc(this_routine,WalkVec2DetsTag)
                DEALLOCATE(WalkVec2Sign)
                CALL LogMemDealloc(this_routine,WalkVec2SignTag)
                ALLOCATE(WalkVec2Dets(NEl,MaxWalkers),stat=ierr)
                CALL LogMemAlloc('WalkVec2Dets',MaxWalkers*NEl,4,this_routine,WalkVec2DetsTag,ierr)
                ALLOCATE(WalkVec2Sign(MaxWalkers),stat=ierr)
                CALL LogMemAlloc('WalkVec2Sign',MaxWalkers,4,this_routine,WalkVec2SignTag,ierr)

            ENDIF

!End of reading in POPSFILE
            CLOSE(17)

        ELSE
!If not reading in from POPSFILE, then we need to initialise the particle positions - start at HF with positive sign

            do j=1,InitWalkers
                do k=1,NEl
                    WalkVecDets(k,j)=FDet(k)
                enddo
                WalkVecSign(j)=.true.
            enddo

        ENDIF

!TotWalkers contains the number of current walkers at each step
        TotWalkers=InitWalkers
        TotWalkersOld=InitWalkers

        RETURN

    END SUBROUTINE InitFCIMCCalc

!This function tells us whether we should create a child particle on nJ, from a parent particle on DetCurr with sign WSign, created with probability Prob
!It returns zero if we are not going to create a child, or -1/+1 if we are to create a child, giving the sign of the new particle
    INTEGER FUNCTION AttemptCreate(DetCurr,WSign,nJ,Prob,IC)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),nJ(NEl),IC
        LOGICAL :: WSign
        REAL*8 :: Prob,Ran2,rat
        TYPE(HElement) :: rh

!Calculate off diagonal hamiltonian matrix element between determinants
        rh=GetHElement2(DetCurr,nJ,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,IC,ECore)
!Divide by the probability of creating the excitation to negate the fact that we are only creating a few determinants
        rat=Tau*abs(rh%v)/Prob
!Stochastically choose whether to create or not according to Ran2
        IF(rat.gt.Ran2(Seed)) THEN
!Child is created - what sign is it?
            IF(WSign) THEN
!Parent particle is positive
                IF((rh%v).gt.0.D0) THEN
                    AttemptCreate=-1     !-ve walker created
                ELSE
                    AttemptCreate=1      !+ve walker created
                ENDIF

            ELSE
!Parent particle is negative
                IF((rh%v).gt.0.D0) THEN
                    AttemptCreate=1      !+ve walker created
                ELSE
                    AttemptCreate=-1     !-ve walker created
                ENDIF
            ENDIF

        ELSE
!No child particle created
            AttemptCreate=0
        ENDIF

        RETURN

    END FUNCTION AttemptCreate

!This function tells us whether we should kill the particle at determinant DetCurr
    LOGICAL FUNCTION AttemptDie(DetCurr)
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),DetLT
        TYPE(HElement) :: rh
        REAL*8 :: Ran2

!Test if determinant is FDet - in a strongly single-configuration problem, this will save time
        IF(DetLT(DetCurr,FDet,NEl).eq.0) THEN
            rh=0.D0
        ELSE
!Calculate the diagonal hamiltonian matrix element for the determinant
            rh=GetHElement2(DetCurr,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!Subtract from the diagonal the value of the lowest hamiltonian matrix element
            rh=rh-Hii
        ENDIF

!Subtract the current value of the shift and multiply by tau
        rh=HElement(Tau)*(rh-(HElement(DiagSft)))
        IF((rh%v).lt.0.D0) THEN
            WRITE(6,*) "Serious problem with -ve death probabilities..."
            STOP "Serious problem with -ve death probabilities..."
        ENDIF
    
!Stochastically choose whether to die or not
        IF(rh.agt.Ran2(Seed)) THEN
!Kill particle
            AttemptDie=.true.
        ELSE
!Particle survives to fight another day...
            AttemptDie=.false.
        ENDIF

        RETURN

    END FUNCTION AttemptDie

END MODULE FciMCMod 

!This is a determinant comparison function
!DetA and DetB are two determinants. Function will return 0 if they are the same, -1 if A.lt.B and +1 if A.gt.B when ordered.
INTEGER FUNCTION DetLT(DetA,DetB,NEl)
    IMPLICIT NONE
    INTEGER :: DetA(NEl),DetB(NEl),NEl,i

    do i=1,NEl
        IF(DetA(i).lt.DetB(i)) THEN
            DetLT=-1
            RETURN
        ELSEIF(DetA(i).gt.DetB(i)) THEN
            DetLT=+1
            RETURN
        ENDIF
    enddo
    DetLT=0
    RETURN

END FUNCTION DetLT


