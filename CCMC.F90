#ifdef PARALLEL
! based on GHB's FciMCPar, this is a version in Coupled Cluster space.


!Note that the main routine called by FciMCPar is PerformCCMCCycPar and is outside the module, so
!FciMCPar doesn't need to use the module

!Based on AttemptDiePar
!This function tells us whether we should kill the particle at determinant DetCurr
!If also diffusing, then we need to know the probability with which we have spawned. This will reduce the death probability.
!The function allows multiple births(if +ve shift) or deaths from the same particle.
!The returned number is the number of deaths if positive, and the number of births if negative.
!Multiple particles can be attempted to die at the same time - here, |WSign| > 1 and the probability of a death will be multiplied by |WSign|
!dProb is an extra probability factor the death probability is multiplied by 
    INTEGER FUNCTION AttemptDieProbPar(DetCurr,Kii,IC,WSign,dProb)
        USE FciMCParMod
        IMPLICIT NONE
        INTEGER :: DetCurr(NEl),iKill,IC,WSign
!        TYPE(HElement) :: rh,rhij
        REAL*8 :: r,rat,Kii
        REAL*8 dProb
        LOGICAL :: tDETinCAS


!Calculate the diagonal hamiltonian matrix element for the determinant
!        rh=GetHElement2(DetCurr,DetCurr,NEl,nBasisMax,G1,nBasis,Brr,NMsh,fck,NMax,ALat,UMat,0,ECore)
!Subtract from the diagonal the value of the lowest hamiltonian matrix element
!        rh=rh-Hii

!Subtract the current value of the shift and multiply by tau
        IF(TFixShiftShell.and.(.not.TSinglePartPhase)) THEN
!            IF((IC.eq.0).or.(IC.eq.2)) THEN
            IF(IC.le.ShellFix) THEN
                rat=Tau*(Kii-FixShift)
            ELSE
                rat=Tau*(Kii-DiagSft)
            ENDIF
        ELSEIF(tFixShiftKii.and.(.not.TSinglePartPhase)) THEN
            IF(Kii.le.FixedKiiCutoff) THEN
                rat=Tau*(Kii-FixShift)
            ELSE
                rat=Tau*(Kii-DiagSft)
            ENDIF
        ELSEIF(tFixCASShift.and.(.not.TSinglePartPhase)) THEN
! The 'TestifDETinCAS' function finds out if the determinant is in the complete active space or not.  If it is, the shift is fixed at 
! the chosen fixed value (FixShift), if not the shift remains as the changing DiagSft value.
           
            tDETinCAS=TestifDETinCAS(DetCurr)
            
            IF(tDETinCAS) THEN
                rat=Tau*(Kii-FixShift)
            ELSE
                rat=Tau*(Kii-DiagSft)
            ENDIF
        ELSE
            rat=Tau*(Kii-DiagSft)
        ENDIF

!If there are multiple particles, decide how many to kill in total...
        rat=rat*abs(WSign)
        rat=rat*dProb

        iKill=INT(rat)
        rat=rat-REAL(iKill)

!Stochastically choose whether to die or not
        IF(tMerTwist) THEN
            CALL genrand_real2(r) 
        ELSE
            CALL RANLUX(r,1)
        ENDIF
        IF(abs(rat).gt.r) THEN
            IF(rat.ge.0.D0) THEN
!Depends whether we are trying to kill or give birth to particles.
                iKill=iKill+1
            ELSE
                iKill=iKill-1
            ENDIF
        ENDIF

        AttemptDieProbPar=iKill

        RETURN

    END FUNCTION AttemptDieProbPar

!Add the excitation in iLutnJ to iLutnI and return it in iLutnI.  iSgn is
!updated with the relevant permutation or set to zero if the excitation is
!disallowed.
SUBROUTINE AddBitExcitor(iLutnI,iLutnJ,iLutRef,nIfD,iSgn)
   use SystemData, only : nEl
   INTEGER nIfD
   INTEGER iLutnI(0:nIfD), iLutnJ(0:nIfD),iLutRef(0:nIfD)
   INTEGER iLutTmp(0:nIfD)
   INTEGER T1,T2,T3
   INTEGER iSgn
! We need to run through the bits of J and I concurrently, setting bits of I
   INTEGER i,j
!NB This is wrong as yet, but just to test compile - it doesn't set the sign correctly.
!First the occ
   do i=0,nIfD
! First mask so we only have 'occupied' orbitals.  The occupieds which are zero
! mean they have been excited from.
      T1=IAND(iLutnI(i),iLutRef(i))
      T2=IAND(iLutnJ(i),iLutRef(i))
      IF(IAND(NOT(T1),NOT(T2)).NE.0) THEN
!The two excitors are exciting from at least one of the same orbitals.  This
!gives a sign of zero.
         iSgn=0
      ENDIF
      iLutTmp(i)=IAND(T1,T2)  !Combine the excitors
   enddo
   do i=0,nIfD
! Now mask so we only have 'virtual' orbitals.  The virtuals which are set mean they have been excited to.
      T3=NOT(iLutRef(i))
      T1=IAND(iLutnI(i),T3)
      T2=IAND(iLutnJ(i),T3)
      IF(IAND(T1,T2).NE.0) THEN
!The two excitors are exciting to at least one of the same orbitals.  This
!gives a sign of zero.
         iSgn=0
      ENDIF
      iLutnI(i)=IOR(iLutTmp(i),IOR(T1,T2)) 
!Combine the excitors and the excitor from and put them back into where they should go
   enddo
   return
END SUBROUTINE AddBitExcitor

!Note that the main routine called by FciMCPar is PerformCCMCCycPar and is outside the module, so
!FciMCPar doesn't need to use the module
    
!Based on PerformCleanFCIMCycPar()

    SUBROUTINE PerformCCMCCycPar()
      USe FCIMCParMod
      Use Determinants, only: GetHElement3
        INTEGER :: VecSlot,i,j,k,l,CopySign,iPartBloom
        INTEGER :: nJ(NEl),ierr,IC,Child,DetCurr(NEl),iLutnJ(0:NIfD)
        REAL*8 :: Prob,rat,HDiagCurr,r
        INTEGER :: iDie,WalkExcitLevel,Proc
        INTEGER :: ExcitLevel,TotWalkersNew,iGetExcitLevel_2,Ex(2,2),WSign,p,Scratch1(2,nSymLabels),Scratch2(2,nSymLabels)
        LOGICAL :: tParity,tFilled
        
! We select up to nEl excitors at a time and store them here
        INTEGER :: SelectedExcitors(0:NIfD,nEl)     
        INTEGER :: SelectedExcitorIndices(nEl)     

! Temporary Storage
        INTEGER iLutnI(0:nIfD)

!The sign of the resultant composite
        INTEGER iSgn


! The probability that we select another excitor (in the cumulative selection
! phase)
        REAL*8 dProbSelNewExcitor
! The prob that we choose the number of Excitors we have done
        REAL*8 dProbNumExcit
! The probl that this (composite) excitor decomposed into the birth/death excitor
        REAL*8 dProbDecompose
! The index of the excitor we've decomposed into.
        INTEGER iPartDie

! The maximum excitation level
        INTEGER nMaxSelExcitors

! The index of the reference det
        INTEGER iHFDet

! The number of excitors we select to make the composite.
        INTEGER iCompositeSize
        TYPE(HElement) Htmp

        iHFDet=1
        dProbSelNewExcitor=0.5

        IF (tTruncSpace) THEN
            nMaxSelExcitors=ICILevel
        ELSE
            nMaxSelExcitors=nEl
        ENDIF
        write(6,*) "Max excitors can be selected.", nMaxSelExcitors

        IF(TDebug.and.(mod(Iter,10).eq.0)) THEN
            WRITE(11,*) Iter,TotWalkers,NoatHF,NoatDoubs,MaxIndex,TotParts
            CALL FLUSH(11)
        ENDIF
        
        CALL set_timer(Walker_Time,30)
        WRITE(6,*) "Number of particles, excitors:",TotParts, TotWalkers
        
!VecSlot indicates the next free position in NewDets
        VecSlot=1
!Reset number at HF and doubles
        NoatHF=0
        NoatDoubs=0
        DetsNorm=0.D0
        iPartBloom=0
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
        ValidSpawnedList(:)=InitialSpawnedSlots(:)

!We loop over a number of samples of cluster space
! Each sample selects a number of excitors together.

! As the number of walkers in the HF reference det is the normalization, we loop
! over each walker there and use it a number of times
! We take the number of walkers as the number of samples to begin with.
        iHFDet=1
        NoatHF=abs(CurrentSign(iHFDet))
        write(6,*) "HF det"
        call WriteBitDet(6,iLutHF,.true.)
        do j=1,TotWalkers
         call WriteBitDet(6,CurrentDets(:,j),.false.)
         call WriteBitEx(6,iLutHF,CurrentDets(:,j),.true.)
        enddo
        do j=1,NoatHF
! Now select a sample of walkers.  We need up to as many walkers as electrons.
            dProbNumExcit=1   !prob of this number of excitors will go here
            dProb=1           !prob of these excitors given this number of excitors goes here
!NB it's possible to select no excitors with this loop.
            do i=1,nMaxSelExcitors
! Calculate the probability that we've reached this far in the loop
               dProbNumExcit=dProbNumExcit*dProbSelNewExcitor

! decide not to choose another walker with this prob.
               call genrand_real2(r)  !On GHB's advice
               if(r.lt.dProbSelNewExcitor) exit
! Select a new random walker
               call genrand_real2(r)  !On GHB's advice
               k=2+floor(r*(TotWalkers-1))    !This selects a unique excitor (which may be multiply populated)
               Write(6,*) "Selected excitor",k
               SelectedExcitorIndices(i)=k
               SelectedExcitors(:,i)=CurrentDets(:,k)
               dProb=(dProb*abs(CurrentSign(k)))/TotParts
               write(6,*) "Prob ",i,": ",(abs(CurrentSign(k))+0.d0)/TotParts," Cuml:", dProb
            enddo
            WRITE(6,*) 'prob out of sel routine.',dProbNumExcit
            if(i.gt.nMaxSelExcitors) THEN !We've been limited by the max number of excitations
               ! Let s be dProbSelNewExcitor, and X be nMaxSelExcitors
               !  The sum of all levels up to X-1 is
               !  (s - s^X) / (1 - s).  We take 1-this to be the prob of
               !  choosing this level
               dProbNumExcit= 1- (dProbSelNewExcitor - dProbNumExcit) / (1-dProbSelNewExcitor)
            ENDIF
            iCompositeSize=i-1  !Save the number of excitors we've selected
            WRITE(6,*) "Iteration ",Iter," Excitors in composite:", iCompositeSize
            do i=1,iCompositeSize
               call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.true.)
            enddo
            Write(6,*) "Level chosen Prob      : ",dProbNumExcit
            Write(6,*) "Select Prob given level: ",dProb
            dProb=dProbNumExcit*dProb
            iLutnI(:)=iLutHF(:)
            iSgn=sign(CurrentSign(iHFDet),1)
            do i=1,iCompositeSize 
               call AddBitExcitor(iLutnI,SelectedExcitors(:,i),iLutHF,nIfD,iSgn)
               Write(6,*) "Results of addition ",i, "Sign ",iSgn
               call WriteBitEx(6,iLutHF,iLutnI,.true.)
            enddo
            CALL FLUSH(6)

            WRITE(6,*) "Chosen det/excitor is:"
            call WriteBitDet(6,iLutnI,.true.)
            CALL FLUSH(6)

!First, decode the bit-string representation of the determinant the walker is on, into a string of naturally-ordered integers
            CALL DecodeBitDet(DetCurr,iLutnI(:),NEl,NIfD)

!Also, we want to find out the excitation level of the determinant - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.
            IF(tTruncSpace) THEN
!We need to know the exact excitation level for truncated calculations.
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,NEl)
            ELSE
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,2)
            ENDIF
!HDiags are stored.
            HDiagCurr=CurrentH(j)

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
            CALL SumEContrib(DetCurr,WalkExcitLevel,CurrentSign(j),CurrentDets(:,j),HDiagCurr)

            tFilled=.false.     !This is for regenerating excitations from the same determinant multiple times. There will be a time saving if we can store the excitation generators temporarily.

!Here, we spawn each particle on the determinant in a seperate attempt.
!we are simply looping over all the particles on the determinant
!This will only be a help if most determinants are multiply occupied.

            CALL GenRandSymExcitScratchNU(DetCurr,iLutnI,nJ,pDoubles,IC,Ex,tParity,exFlag,Prob,Scratch1,Scratch2,tFilled)

!Calculate number of children to spawn
            IF(TTruncSpace) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.

                IF(WalkExcitLevel.eq.(ICILevel-1)) THEN
!The current walker is one below the excitation cutoff - if IC is a double, then could go over - we need to check

                    IF(IC.eq.2) THEN
!Need to check excitation level of excitation
                        ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,ICILevel)
                        IF(ExcitLevel.gt.ICILevel) THEN
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                            Child=0
                        ELSE
                            Child=AttemptCreatePar(DetCurr,iLutnI,iSgn,nJ,iLutnJ,Prob,IC,Ex,tParity,1,.false.)

                        ENDIF
                    ELSE
                        Child=AttemptCreatePar(DetCurr,iLutnI,iSgn,nJ,iLutnJ,Prob,IC,Ex,tParity,1,.false.)
                    ENDIF

                ELSEIF(WalkExcitLevel.eq.ICILevel) THEN
!Walker is at the excitation cutoff level - all possible excitations could be disallowed - check the actual excitation level
                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,ICILevel)
                    IF(ExcitLevel.gt.ICILevel) THEN
!Attempted excitation is above the excitation level cutoff - do not allow the creation of children
                        Child=0
                    ELSE
                        Child=AttemptCreatePar(DetCurr,iLutnI,iSgn,nJ,iLutnJ,Prob,IC,Ex,tParity,1,.false.)
                    ENDIF
                ELSE
!Excitation cannot be in a dissallowed excitation level - allow it as normal
                    Child=AttemptCreatePar(DetCurr,iLutnI,iSgn,nJ,iLutnJ,Prob,IC,Ex,tParity,1,.false.)
                ENDIF 
             
            ELSE
!SD Space is not truncated - allow attempted spawn as usual

                Child=AttemptCreatePar(DetCurr,iLutnI,iSgn,nJ,iLutnJ,Prob,IC,Ex,tParity,1,.false.)

            ENDIF
             
            IF(Child.ne.0) THEN
!We want to spawn a child - find its information to store

!                    WRITE(6,*) "Spawning particle to:",nJ(:)
!                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,NEl)
!                    WRITE(6,*) "Excitlevel:", ExcitLevel
                 NoBorn=NoBorn+abs(Child)     !Update counter about particle birth
                 IF(IC.eq.1) THEN
                     SpawnFromSing=SpawnFromSing+abs(Child)
                 ENDIF

                 IF(abs(Child).gt.25) THEN
!If more than 25 particles are created in one go, then log this fact and print out later that this has happened.
                     IF(abs(Child).gt.abs(iPartBloom)) THEN
                         IF(IC.eq.1) THEN
                             iPartBloom=-abs(Child)
                         ELSE
                             iPartBloom=abs(Child)
                         ENDIF
                     ENDIF
!                        WRITE(6,"(A,I10,A)") "LARGE PARTICLE BLOOM - ",Child," particles created in one attempt."
!                        WRITE(6,"(A,I5)") "Excitation: ",IC
!                        WRITE(6,"(A,G25.10)") "PROB IS: ",Prob
!                        CALL FLUSH(6)
                 ENDIF

!We need to calculate the bit-representation of this new child. This can be done easily since the ExcitMat is known.
                 IF(.not.tHPHF) CALL FindExcitBitDet(iLutnI,iLutnJ,IC,Ex,NIfD)

!In direct annihilation, we spawn particles into a seperate array, but we do not store them contiguously in the SpawnedParts/SpawnedSign arrays.
!The processor that the newly-spawned particle is going to be sent to has to be determined, and then it will get put into the the appropriate element determined by ValidSpawnedList.

                 Proc=DetermineDetProc(iLutnJ)   !This wants to return a value between 0 -> nProcessors-1
!                    WRITE(6,*) iLutnJ(:),Proc,ValidSpawnedList(Proc),Child,TotWalkers
!                    CALL FLUSH(6)
                 SpawnedParts(:,ValidSpawnedList(Proc))=iLutnJ(:)
                 SpawnedSign(ValidSpawnedList(Proc))=Child
                 ValidSpawnedList(Proc)=ValidSpawnedList(Proc)+1

                 Acceptances=Acceptances+ABS(Child)      !Sum the number of created children to use in acceptance ratio
             
            ENDIF   !End if child created


!We now have to decide whether the parent particle (j) wants to die or not...
!we can have multiple particles on the same determinant - these can be stochastically killed at the same time.

! We have to decompose our composite excitor into one of its parts.  
            IF(iCompositeSize.GT.0) THEN
               dProbDecompose=1.0/iCompositeSize
               call genrand_real2(r)  !On GHB's advice
               k=1+floor(r*iCompositeSize)
               iPartDie=SelectedExcitorIndices(k)
            ELSE
               dProbDecompose=1
               iPartDie=iHFDet
            ENDIF 
            Write(6,'(A,I)',advance='no') "Killing at excitor: ",iPartDie
            call WriteDet(6,CurrentDets(:,iPartDie),nEl,.true.)
            call WriteBitEx(6,iLutHF,CurrentDets(:,iPartDie),.true.)
!Now get the full representation of the dying excitor
            CALL DecodeBitDet(DetCurr,CurrentDets(:,iPartDie),NEl,NIfD)
            Htmp=GetHElement3(DetCurr,DetCurr,0)
            HDiagCurr=REAL(Htmp%v,r2)
            HDiagCurr=HDiagCurr-Hii

!Die with a probability as normal, but biased because we only selected this set
!of excitors with dProb.
            iDie=AttemptDieProbPar(DetCurr,HDiagCurr,WalkExcitLevel,iSgn,1/dProb)
            NoDied=NoDied+iDie          !Update death counter


!iDie can be positive to indicate the number of deaths, or negative to indicate the number of births
!We slot the particles back into the same array and position VecSlot if the particle survives. If it dies, then j increases, moving onto the next
!entry, but VecSlot remains where it is, meaning that j should never be less that VecSlot

            IF(iSgn.le.0) THEN
                CopySign=iSgn+iDie    !Copy sign is the total number of particles x sign that we want to copy accross.
                IF(CopySign.gt.0) THEN
!If we are copying to the main array, we have to ensure that we maintain sign-coherence in the array. Therefore, if we are spawning anti-particles,
!it wants to go in the spawning array, rather than the main array, so it has a chance to annihilate. However, since anti-particles should not be created
!in normal circumstances, we will remove this possibility.
                    CALL Stop_All("PerformFCIMCyc","Creating anti-particles")
                ENDIF
            ELSE
                CopySign=iSgn-iDie
                IF(CopySign.lt.0) THEN
                    CALL Stop_All("PerformFCIMCyc","Creating anti-particles")
                ENDIF
            ENDIF

!We only have a single array, therefore surviving particles are simply transferred back into the original array.
!We should not get to the case where we want to overwrite particles that we haven't even got to yet.
            IF(CopySign.ne.0) THEN

                CurrentDets(:,VecSlot)=CurrentDets(:,iPartDie)
                CurrentSign(VecSlot)=CopySign
                CurrentH(VecSlot)=HDiagCurr
                VecSlot=VecSlot+1

            ENDIF   !To kill if

!Finish cycling over determinants
        enddo


!***Birth/death processes finished. Tidy up and then annihilate.

!SumWalkersCyc calculates the total number of walkers over an update cycle on each process
        SumWalkersCyc=SumWalkersCyc+(INT(TotParts,i2))
!        WRITE(6,*) "Born, Die: ",NoBorn, NoDied

!Since VecSlot holds the next vacant slot in the main determinant array, TotWalkers will be one less than this.
        TotWalkersNew=VecSlot-1

!Output if there has been a particle bloom this iteration. A negative number indicates that particles were created from a single excitation.
        IF(iPartBloom.ne.0) THEN
            WRITE(6,"(A,I10,A)") "LARGE PARTICLE BLOOMS in iteration ",Iter
            IF(iPartBloom.gt.0) THEN
                WRITE(6,"(A,I10,A)") "A max of ",abs(iPartBloom)," particles created in one attempt from double excit."
            ELSE
                WRITE(6,"(A,I10,A)") "A max of ",abs(iPartBloom)," particles created in one attempt from single excit."
            ENDIF
        ENDIF


        rat=(TotWalkersNew+0.D0)/(MaxWalkersPart+0.D0)
        IF(rat.gt.0.95) THEN
            WRITE(6,*) "*WARNING* - Number of particles/determinants has increased to over 95% of MaxWalkersPart"
            CALL FLUSH(6)
        ENDIF

!Need to test whether any of the sublists in the spawning array are getting to the end of their allotted space.
        IF(nProcessors.gt.1) THEN
            do i=0,nProcessors-1
                rat=(ValidSpawnedList(i)-InitialSpawnedSlots(i))/(InitialSpawnedSlots(1)+0.D0)
                IF(rat.gt.0.95) THEN
                    WRITE(6,*) "*WARNING* - Highest processor spawned particles has reached over 95% of MaxSpawned"
                    CALL FLUSH(6)
                ENDIF
            enddo
        ELSE
            rat=(ValidSpawnedList(0)+0.D0)/(MaxSpawned+0.D0)
            IF(rat.gt.0.9) THEN
                WRITE(6,*) "*WARNING* - Number of spawned particles has reached over 90% of MaxSpawned"
                CALL FLUSH(6)
            ENDIF
        ENDIF
        
        CALL halt_timer(Walker_Time)
        CALL set_timer(Annihil_Time,30)
!        CALL MPI_Barrier(MPI_COMM_WORLD,error)
!        WRITE(6,*) "Get into annihilation"
!        CALL FLUSH(6)

!This is the direct annihilation algorithm. The newly spawned walkers should be in a seperate array (SpawnedParts) and the other list should be ordered.
        CALL DirectAnnihilation(TotWalkersNew)

        CALL halt_timer(Annihil_Time)
        
    END SUBROUTINE PerformCCMCCycPar


#endif
