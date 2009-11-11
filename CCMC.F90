#include "macros.h"
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
        use SystemData, only: nEl, tMerTwist
        use CalcData , only : TFixShiftShell,ShellFix,FixShift,DiagSft,FixedKiiCutoff,Tau,tFixShiftKii
        use CalcData , only : tFixCASShift
        use FciMCData
        use FciMCParMod, only: TestifDETinCAS
        use mt95
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



!Note that the main routine called by FciMCPar is PerformCCMCCycPar and is outside the module, so
!FciMCPar doesn't need to use the module
    
!Based on PerformCleanFCIMCycPar()

!First we loop over particles at the HF det, and for each, we select a group of excitors to make a composite excitor.
! From this we make an excitation, and accept or reject a particle created there.

! Next we decide whether or not that composite excitor should die.  Its death is caused by the death of one of its composites.  We save a list of the excitors which need to die.

!Then we go through the list of excitors which should die, and reduce their weight on the list, and remove if necessary.
!Next we annihilate.
#ifdef PARALLEL
    SUBROUTINE PerformCCMCCycPar()
      USe FCIMCParMod
      use CCMCData
      Use Determinants, only: GetHElement3
      Use Logging, only: CCMCDebug
      use DetBitOps, only: DecodeBitDet
        IMPLICIT NONE
        INTEGER :: VecSlot,i,j,k,l,CopySign,iPartBloom
        INTEGER :: nJ(NEl),ierr,IC,Child,DetCurr(NEl),iLutnJ(0:NIfTot)
        REAL*8 :: Prob,rat,HDiagCurr,r
        INTEGER :: iDie,WalkExcitLevel,Proc
        INTEGER :: ExcitLevel,TotWalkersNew,iGetExcitLevel_2,Ex(2,2),WSign,p,Scratch1(ScratchSize),Scratch2(ScratchSize)
        LOGICAL :: tParity,tFilled
        
! We select up to nEl excitors at a time and store them here
        INTEGER :: SelectedExcitors(0:NIfTot,nEl)     
        INTEGER :: SelectedExcitorIndices(nEl)     

! Temporary Storage
        INTEGER iLutnI(0:nIfTot)

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
!and the number at it
        INTEGER HFCount

! The number of excitors we select to make the composite.
        INTEGER iCompositeSize
        TYPE(HElement) Htmp

        INTEGER, allocatable :: iKillDetIndices(:,:)
        INTEGER iDeaths
!Debug print level higher is more
        INTEGER iDebug

        INTEGER TotRealWalkers
        REAL*8 dProbNorm,dClusterProb,dProb
        INTEGER iExcitor
        INTEGER iMinEx,iMaxEx,iMaxExTemp
!iMaxExcitorSelections is the number of times we decide to loop for each particle.
!  it includes the first (exact) loop where we consider the particle by itself.
        INTEGER iMaxExcitorSelections
        LOGICAL tSuccess

        REAL*8 dNGenComposite  ! The number of ways the composite could've been generated.

        REAL*8 dT1Sq
        INTEGER AttemptDieProbPar
!        REAL*8 AJWTProjE
        INTEGER iCumlExcits,iLeftHere(nEl)
        INTEGER iCurrentCompositeSize
        INTEGER,save :: nClusterBirths
        INTEGER,save :: nClusterDeaths   
        INTEGER,save :: nClusterChildren
        INTEGER nClusters
        REAL*8 dClusterProbs
        nClusters=0
        dClusterProbs=0

        IF(mod(Iter,StepsSft).eq.0) THEN
    
           nClusterBirths=0 
           nClusterDeaths=0 
           nClusterChildren=0 
        endif
        if(tCCMCFCI) then
         iMaxExcitorSelections=1
        else
         iMaxExcitorSelections=2
        endif

        if(Iter.eq.1) dT1SqCuml=0 
        iDebug=CCMCDebug
!Number of excitors dying
        iDeaths=0

!        AJWTProjE=0
        dT1Sq=0
        IF(iDebug.gt.0) WRITE(6,*) "Entering CCMC Cycle"
        iHFDet=1
        dProbSelNewExcitor=0.5


        IF(TDebug.and.(mod(Iter,10).eq.0)) THEN
            WRITE(11,*) Iter,TotWalkers,NoatHF,NoatDoubs,MaxIndex,TotParts
            CALL FLUSH(11)
        ENDIF


        IF(tHistSpawn.or.tCalcFCIMCPsi) HistMinInd(1:NEl)=FCIDetIndex(1:NEl)    !This is for the binary search when histogramming
         !This info is destroyed by SumEContrib and needs to be reset each cycle
        
        CALL set_timer(Walker_Time,30)
        IF(iDebug.gt.0) WRITE(6,*) "Number of particles, excitors:",TotParts, TotWalkers
        
!Reset number at HF and doubles
        NoatHF=0
        NoatDoubs=0
        iPartBloom=0
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
        ValidSpawnedList(:)=InitialSpawnedSlots(:)

!We loop over a number of samples of cluster space
! Each sample selects a number of excitors together.

! As the number of walkers in the HF reference det is the normalization, we loop
! over each walker there and use it a number of times
! We take the number of walkers as the number of samples to begin with.
        CALL BinSearchParts(iLutHF,1,TotWalkers,iHFDet,tSuccess)
        if(.not.tSuccess) then
            WRITE(6,*) "WARNING: Cannot find HF det in particle list"
            HFcount=1
            iHFDet=-1
        else
            HFcount=abs(CurrentSign(iHFDet))
        endif

        allocate(iKillDetIndices(2,TotParts*2))

        IF(iDebug.gt.1) THEN
         write(6,*) "HF det"
         call WriteBitDet(6,iLutHF,.true.)
         write(6,*) "Particle list"
         do j=1,TotWalkers
            write(6,'(i7)',advance='no') CurrentSign(j)
            call WriteBitEx(6,iLutHF,CurrentDets(:,j),.true.)
!            do l=0,nIfD
!               Write(6,'(i)', advance='no') CurrentDets(l,j)
!            enddo
!            write(6,*)
         enddo
        endif
!TotRealWalkers gets updated with the number of non-zero walkers at each stage.
        TotRealWalkers=TotWalkers
        iCumlExcits=0
        do j=1,TotWalkers
!#if 0
          IF(iDebug.gt.4) WRITE(6,*) "Iteration ",Iter,':',j
! Deal with T_1^2
          CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),WalkExcitLevel,2)
          if(WalkExcitLevel.eq.1) then
            do l=j,TotWalkers
               CALL FindBitExcitLevel(iLutHF,CurrentDets(:,l),WalkExcitLevel,2)
               if(WalkExcitLevel.eq.1) then
                  iSgn=CurrentSign(j)
                  iSgn=iSgn*CurrentSign(l)
                  iLutnI(:)=CurrentDets(:,j)
                  call AddBitExcitor(iLutnI,CurrentDets(:,l),iLutHF,iSgn)
                  if(iSgn.ne.0) then
                     CALL DecodeBitDet(DetCurr,iLutnI(:))
                     Htmp=GetHElement3(HFDet, DetCurr,2)
                     dT1Sq=dT1Sq+(Real(Htmp%v,r2)*iSgn)
                     !WRITE(6,'(A,I,2G)', advance='no') 'T1',iSgn,real(Htmp%v,r2),dT1Sq
                     !call WriteBitEx(6,iLutHF,CurrentDets(:,j),.false.)
                     !call WriteBitEx(6,iLutHF,CurrentDets(:,l),.false.)
                     !call WriteBitEx(6,iLutHF,iLutnI,.true.)
                  endif
               endif
            enddo
          endif 
!          if(WalkExcitLevel.eq.1.or.WalkExcitLevel.eq.2) then
!            CALL DecodeBitDet(DetCurr,CurrentDets(:,j))
!            Htmp=GetHElement3(HFDet, DetCurr,WalkExcitLevel)
!            AJWTProjE=AJWTProjE+(Real(Htmp%v,r2)*CurrentSign(j))
!          endif
!#endif
!Go through all particles on the current walker det
          do l=1,abs(CurrentSign(j))
            iCumlExcits=iCumlExcits+1
            if(j.eq.iHFDet) then
               iMaxEx=1
            else
               if(tExactCluster) then  !Go through all possible clusters explicitly
                  IF (tTruncSpace) THEN
                     nMaxSelExcitors=ICILevel+2  !We need to be able to couple (say) 4 singles to make a quad and then spawn back to the sdoubles space
                  ELSE
                     nMaxSelExcitors=nEl
                  ENDIF
                  iMaxExTemp=1
                  iMaxEx=0
                  k=TotParts-iCumlExcits  !iCumlExcits includes this det.
                  if(j.lt.iHFDet) k=k-HFcount
!Count the number of allowed composites - this allows for all numbers of composites
                  if(iDebug.gt.5) WRITE(6,*) "Counting Excitations:  Level,#, Cuml"
                  do i=k,max(k-(nMaxSelExcitors-1)+1,1),-1  !-1 because we've already chosen this excitor.  +1 to exclude the end (off by 1 problem)
                     iMaxExTemp=(iMaxExTemp*i)/(k+1-i)
                     iMaxEx=iMaxEx+iMaxExTemp
                     if(iDebug.gt.5) WRITE(6,*) k+1-i,iMaxExTemp,iMaxEx
                  enddo
                  if(k.eq.0) iMaxEx=0
                  iMaxEx=iMaxEx+1  !Account for the non-cluster
                  if(iDebug.gt.4) write(6,*) iMaxEx-1, ' combined excitations available from det ',j
!iMaxEx is now the number of possible clusters which can be made from this and walkers after it.  we go through them in turn.
               else
                  iMaxEx=iMaxExcitorSelections
               endif
            endif
! tNewExcitor will be set at the end of this loop if we decide to come 
            do iExcitor=1,iMaxEx


               if(iExcitor.eq.1) then  !Deal with all excitors singly
                  if(iDebug.gt.4) write(6,*) 'Plain old excitor.'
!just select the single excitor
                  iLutnI(:)=CurrentDets(:,j)
                  iSgn=sign(1,CurrentSign(j))
                  dClusterProb=1 
                  dProbNorm=1
                  iCompositeSize=0
                  CALL DecodeBitDet(DetCurr,iLutnI(:))
!Also take into account the contributions from the dets in the list
                  HDiagCurr=CurrentH(j)
                  if(tHistSpawn) then
                     CALL FindBitExcitLevel(iLutHF,iLutnI,WalkExcitLevel,nEl)
                  else
                     CALL FindBitExcitLevel(iLutHF,iLutnI,WalkExcitLevel,2)
                  endif
                  CALL SumEContrib(DetCurr,WalkExcitLevel,iSgn,iLutnI,HDiagCurr,1.d0)
               else
                  if(iDebug.gt.4) write(6,*) 'Excitor composite number ',iExcitor
! Now select a sample of walkers.  We need up to as many walkers as electrons.
                  dProbNumExcit=1
!dProbSelNewExcitor   !prob of this number of excitors will go here
                  IF (tTruncSpace) THEN
                     nMaxSelExcitors=ICILevel+2  !We need to be able to couple (say) 4 singles to make a quad and then spawn back to the sdoubles space
                  ELSE
                     nMaxSelExcitors=nEl
                  ENDIF
            
                  if(nMaxSelExcitors.gt.(TotWalkers-1)) nMaxSelExcitors=TotWalkers-1
                  IF(iDebug.gt.4) write(6,*) "Max excitors can be selected.", nMaxSelExcitors
                  if(nMaxSelExcitors.lt.2) exit  !If we can't choose a new excit, we leave the loop
                  dClusterProb=1           !prob of these excitors given this number of excitors goes here
   !NB it's possible to select no excitors with this loop.
                  dProbNorm=1
!We force the first excitor to be the current det.
                  SelectedExcitorIndices(1)=j
                  SelectedExcitors(:,1)=CurrentDets(:,j)

!We keep track of the number of ways we could've generated this composite
                  dNGenComposite=abs(CurrentSign(j))


                  if(tExactCluster) then  !Each time we're here, generate another cluster, up to all of them.
                     if(iExcitor.eq.2) then
!Setup the exact cluster if we need to 
!Pretend we've just finished a single excitor
                        iCurrentCompositeSize=1
                        iLeftHere(1)=abs(CurrentSign(j))-1
                     endif
                     if(iDebug.gt.5) then
                       WRITE(6,*) "EXACT CLUSTER IN"
                        do k=2,iCurrentCompositeSize 
                           if(iDebug.gt.5) WRITE(6,'(A1,I5,A1,I5,A1)',advance='no') '[',SelectedExcitorIndices(k),',',iLeftHere(k),'] '
                        enddo
                        WRITE(6,*)
                     endif
!For a given CompositeSize we generate all clusters, then increase composite size.
!To generate a new cluster, we start from the end and try to point to the next excitor.  This is complicated by the indexing
!SelectedExcitors(i) is the excitor, and iLeftHere(i) is the index of the current excitor in that list of excitors
! (counting down from abs(CurrentSigns(SelectedExcitors(i)))-1 to 0.  
!If there are no more of the selectedexcitor, then we try the next excitor.
!If there are no more excitors to choose from, we step back in the selectedexcitor list to the previous index and try to increment that.
! After having incremented somewhere in the list, we fill in up to the composite size from then on.
! If this isn't possible (because we've enumerated all composite excitors of that composite size), we move to a larger composite size, and fill that in afresh.

!Go from the end and try to increment
                     do i=iCurrentCompositeSize,2,-1
                        if(iLeftHere(i).gt.0) then  !We've got at least one  more at this level, so just use that
                           iLeftHere(i)=iLeftHere(i)-1
                           exit
                        else  !try to get the next excitor
                           SelectedExcitorIndices(i)=SelectedExcitorIndices(i)+1
                           !we're not allowed to select the HF det.
                           if (SelectedExcitorIndices(i).eq.iHFDet) SelectedExcitorIndices(i)=SelectedExcitorIndices(i)+1
                           !if we've reached the end, then we need to go back up and try to increment that one.
                           if (SelectedExcitorIndices(i).gt.TotWalkers) cycle
                           !Hooray - we've got an allowed excitor
                           iLeftHere(i)=abs(CurrentSign(SelectedExcitorIndices(i)))-1
                           SelectedExcitors(:,i)=CurrentDets(:,SelectedExcitorIndices(i))
                           exit
                        endif
                     enddo

                     k=i

                     do while (k.le.iCurrentCompositeSize)
                        !If we fell off the end of the previous loop, then we need to increase the composite size.
                        if (i.eq.1) iCurrentCompositeSize=iCurrentCompositeSize+1
                        if (iCurrentCompositeSize.gt.nMaxSelExcitors) exit
                        if (iDebug.gt.5) WRITE(6,*) "Filling from ",i+1," to ",iCurrentCompositeSize
                        ! now go through from the one after the one we incrememented and fill the rest in sequentially
                        do k=i+1,iCurrentCompositeSize
                           if(iLeftHere(k-1).gt.0) then  !We've got at least one  more at this level, so just use that
                              iLeftHere(k)=iLeftHere(k-1)-1
                              SelectedExcitorIndices(k)=SelectedExcitorIndices(k-1)
                              SelectedExcitors(:,k)=CurrentDets(:,SelectedExcitorIndices(k))
                              cycle
                           else  !try to get the next excitor
                              SelectedExcitorIndices(k)=SelectedExcitorIndices(k-1)+1
                              !we're not allowed to select the HF det.
                              if (SelectedExcitorIndices(k).eq.iHFDet) SelectedExcitorIndices(k)=SelectedExcitorIndices(k)+1
                              !if we've reached the end, then we need to increase increase the composite size.
                              ! we do this by exiting the k loop, which will loop us round again, incrementing the composite size if i=1
                              if (SelectedExcitorIndices(k).gt.TotWalkers) then
                                 if(iDebug.gt.5) WRITE(6,*) "Out of excitors.  Increasing loop size."
                                 i=1
                                 exit
                              endif
                              !Hooray - we've got an allowed excitor
                              iLeftHere(k)=abs(CurrentSign(SelectedExcitorIndices(k)))-1
                              SelectedExcitors(:,k)=CurrentDets(:,SelectedExcitorIndices(k))
                           endif
                        enddo !for k
                     enddo !for having successfully filled out the excitor list
   
                     if (iCurrentCompositeSize.gt.nMaxSelExcitors) then
                        if(iDebug.gt.5) WRITE(6,*) "Not allowed >",nMaxSelExcitors," so exiting this excitor." 
                        exit  !we've reached the end prematurely.  This is allowed
                     endif
                     dClusterProb=1
!Account for intermediate Normalization

                     iCompositeSize=iCurrentCompositeSize
                     if(iDebug.gt.5) WRITE(6,*) "EXACT CLUSTER OUT"
                     do k=2,iCompositeSize 
                        if(iDebug.gt.5) WRITE(6,'(A1,I5,A1,I5,A1)',advance='no') '[',SelectedExcitorIndices(k),',',iLeftHere(k),'] '
                        dProbNorm=dProbNorm*HFCount
                        dNGenComposite=dNGenComposite*abs(CurrentSign(SelectedExcitorIndices(k)))  !For each new excit added to the composite, we multiply up to count the number of ways we could've generated it.
                     enddo
                     if(iDebug.gt.5) WRITE(6,*)
                  else
! We wish to sum:
! sum_A [ (1/2) sum_{B\ne A} [ (1/3) sum_{C\ne A,B} ... ]]
!  since the sums over B and C are too large a space to do fully, we stochastically sum them:
! sum_A [ (1/2)  1/X sum_i^X 1/p(B_i)  [  (1/3)  1/Y sum_j^Y 1/p(C_j) ... ]]
!   where X and Y are the number of samples taken in the sums.

! If B_i and C_j were selected uniformly randomly out of N-1 particles,  (remember we've already chosen A), 
!  then p(B_i) would be 1/(N-1)  and  p(C_j) would be 1/(N-2)
!  The total prefactor would then be
!  [1/(XY)]  [(N-1) (N-2) / ( 2 . 3 )]     somwehat akin to the combinatorial coefficient.

!  We wish to sample X and Y = 1, but p(B_i) and p(C_j) are not uniform.
!  All in all, we then wish to use the total generation probability of this composite stochastically.
!
!  In this case it would be  [  p(B_i)  p(C_j) / ( 2 . 3 )  ]  p(Choosing this level)

!  In fact, we need the contribution of the composite excitor to vary as (HFcount)^-(iCompositeSize-1)
!    to account for the intermediate normalization.  We make this happen by multiplying the probability of 
!    having generated the excitor by (HFcount)^(iCompositeSize-1)  (i.e. one HFcount for each level)
                     do i=2,nMaxSelExcitors
   ! Calculate the probability that we've reached this far in the loop
                        !We must have at least one excitor, so we cannot exit here.
                        if(i.gt.2) then
                           call genrand_real2(r)  !On GHB's advice
                           if(r.lt.dProbSelNewExcitor) exit
                        endif
            
      ! decide not to choose another walker with this prob.
                        dProbNumExcit=dProbNumExcit*dProbSelNewExcitor
      ! Select a new random walker
                        k=0
                        do while (k.eq.0)
                           call genrand_real2(r)  !On GHB's advice
                           k=1+floor(r*(TotWalkers-1))    !This selects a unique excitor (which may be multiply populated)
                           if(k.ge.iHFDet) k=k+1          !We're not allowed to select the HF det
                           if(CurrentSign(k).eq.0) k=0
                        enddo
                        dNGenComposite=dNGenComposite*abs(CurrentSign(k))  !For each new excit added to the composite, we multiply up to count the number of ways we could've generated it.
                        IF(iDebug.gt.5) Write(6,*) "Selected excitor",k
                        SelectedExcitorIndices(i)=k
                        SelectedExcitors(:,i)=CurrentDets(:,k)
                        dClusterProb=(dClusterProb*abs(CurrentSign(k)))/((TotParts-HFcount))

!Account for intermediate Normalization
                        dProbNorm=dProbNorm*HFCount

!Account for possible orderings of selection
                        dProbNorm=dProbNorm*i
                        IF(iDebug.gt.5) WRITE(6,*) "TotParts,HFCount:",TotParts,HFcount
                        IF(iDebug.gt.5) write(6,*) "Prob ",i,": ",(abs(CurrentSign(k))+0.d0)/(TotParts-HFcount)," Cuml:", dClusterProb
                     enddo
                     IF(iDebug.gt.5) WRITE(6,*) 'prob out of sel routine.',dProbNumExcit
                     if(i.gt.nMaxSelExcitors) THEN !We've been limited by the max number of excitations
                        ! Let s be dProbSelNewExcitor, and X be nMaxSelExcitors
                        !  The sum of all levels up to X-1 is
                        !  (s - s^X) / (1 - s).  We take 1-this to be the prob of
                        !  choosing this level
                        dProbNumExcit= 1- (dProbSelNewExcitor - dProbNumExcit) / (1-dProbSelNewExcitor)
                     ENDIF
                     iCompositeSize=i-1  !Save the number of excitors we've selected   
                     dClusterProb=dClusterProb*(iMaxEx-1)
                  endif
                  IF(iDebug.gt.3) then
                     WRITE(6,*) " Excitors in composite:", iCompositeSize
                     do i=1,iCompositeSize
                        write(6,'(I2)',advance='no') sign(1,CurrentSign(SelectedExcitorIndices(i)))
                        call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.true.)
                     enddo
                     Write(6,*) "Level chosen Prob      : ",dProbNumExcit
                     Write(6,*) "Select Prob given level: ",dClusterProb
                     Write(6,*) "Prob norm              : ",dProbNorm
                  endif
                  dClusterProb=dProbNumExcit*dClusterProb
! Lets check the probs
                  nClusters=nClusters+1
                  dClusterProbs=dClusterProbs+1/dClusterProb
                  iLutnI(:)=CurrentDets(:,j)
                  iSgn=sign(1,CurrentSign(j)) !The sign of the first excitor
                  do i=2,iCompositeSize 
                     iSgn=iSgn*sign(1,CurrentSign(SelectedExcitorIndices(i)))
                     call AddBitExcitor(iLutnI,SelectedExcitors(:,i),iLutHF,iSgn)
                     IF(iDebug.gt.3) Write(6,*) "Results of addition ",i, "Sign ",iSgn,':'
                     if(iSgn.eq.0) exit
                     IF(iDebug.gt.3) call WriteBitEx(6,iLutHF,iLutnI,.true.)
                  enddo
                  IF(iDebug.gt.0) CALL FLUSH(6)
                  if(iSgn.eq.0) cycle
               endif

               IF(iDebug.gt.4) then
                  WRITE(6,*) "Chosen det/excitor is:"
                  call WriteBitDet(6,iLutnI,.true.)
                  CALL FLUSH(6)
               endif

!First, decode the bit-string representation of the determinant the walker is on, into a string of naturally-ordered integers
               CALL DecodeBitDet(DetCurr,iLutnI(:))
               CALL FindBitExcitLevel(iLutHF,iLutnI(:),WalkExcitLevel,nEl)
               if(iDebug.gt.4) WRITE(6,*) "Excitation Level ", WalkExcitLevel



               tFilled=.false.     !This is for regenerating excitations from the same determinant multiple times. There will be a time saving if we can store the excitation generators temporarily.

!Here, we spawn each particle on the determinant in a seperate attempt.
!we are simply looping over all the particles on the determinant
!This will only be a help if most determinants are multiply occupied.

               CALL GenRandSymExcitScratchNU(DetCurr,iLutnI,nJ,pDoubles,IC,Ex,tParity,exFlag,Prob,Scratch1,Scratch2,tFilled)
!We need to calculate the bit-representation of this new child. This can be done easily since the ExcitMat is known.
               IF(.not.tHPHF) CALL FindExcitBitDet(iLutnI,iLutnJ,IC,Ex,NIfD)
               IF(iDebug.gt.4) then
                   WRITE(6,*) "Random excited det level ",iC
                   call WriteDet(6,nJ,nEl,.true.)
                   Write(6,*) "Prob ex|from",Prob
               endif
!Prob is the Prob of choosing nJ from nI
!dClusterProb is the Probability of having chosen this cluster excitor (normalized such that <1/dClusterProb> = 1)
!dProbNorm is the renormalization factor for this level of excitors - decreasing by HFcount for each extra level of excitors
               Prob=Prob*dClusterProb*dProbNorm  !Now include the prob of choosing the det we spawned from
               if(iDebug.gt.4) Write(6,*) "Prob ex tot",Prob
               if(iCompositeSize.gt.1) dProb=dProb
!Calculate number of children to spawn
               IF(TTruncSpace) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.
                   IF(CheckAllowedTruncSpawn(WalkExcitLevel,nJ,iLutnJ,IC)) THEN
                       Child=AttemptCreatePar(DetCurr,iLutnI,iSgn,nJ,iLutnJ,Prob,IC,Ex,tParity,1,.false.)
                   ELSE
                       Child=0
                   ENDIF
               ELSE
!SD Space is not truncated - allow attempted spawn as usual
                   Child=AttemptCreatePar(DetCurr,iLutnI,iSgn,nJ,iLutnJ,Prob,IC,Ex,tParity,1,.false.)
               ENDIF

               IF(iDebug.gt.4.or.((iDebug.eq.3.or.iDebug.eq.4).and.Child.ne.0)) THEN
   !We've not printed this out before
                  call WriteBitEx(6,iLutHF,iLutnI,.false.)
                  write(6,'(A)',advance='no') ,' ==> '
                  call WriteBitEx(6,iLutHF,iLutnJ,.false.)
                  WRITE(6,'(A,I7)',advance='no') "Children:",Child             
                  if(iDebug.eq.3.and.iCompositeSize.gt.1) THEN
                     write(6,'(A)',advance='no') ' from '
                     do i=1,iCompositeSize
                        call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.false.)
                     enddo
                  endif
                  write(6,*)
               endif

               IF(Child.ne.0) THEN
                  if(iCompositeSize.gt.1) nClusterChildren=nClusterChildren+1 
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



   ! We have to decompose our composite excitor into one of its parts.  
               IF(iCompositeSize.GT.0) THEN
                  dProbDecompose=dNGenComposite
!                  dProb=dProb/HFcount
                  call genrand_real2(r)  !On GHB's advice
                  k=1+floor(r*iCompositeSize)
                  iPartDie=SelectedExcitorIndices(k)
   !Now get the full representation of the dying excitor
                  CALL DecodeBitDet(DetCurr,iLutnI)
                  Htmp=GetHElement3(DetCurr,DetCurr,0)
                  HDiagCurr=REAL(Htmp%v,r2)
                  HDiagCurr=HDiagCurr-Hii
                  iSgn=sign(1,CurrentSign(iPartDie))
               ELSE
                  dProbDecompose=1
                  iPartDie=j
                  HDiagCurr=CurrentH(j)
               ENDIF 
               dProb=dClusterProb*dProbDecompose

!Also, we want to find out the excitation level of the determinant - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.
               IF(tTruncSpace) THEN
!We need to know the exact excitation level for truncated calculations.
                   CALL FindBitExcitLevel(iLutHF,iLutnI,WalkExcitLevel,NEl)
               ELSE
                   CALL FindBitExcitLevel(iLutHF,iLutnI,WalkExcitLevel,2)
               ENDIF
!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
!               WRITE(6,'(A,3G,2I)') "PPP",dProb,dProbNorm,dProb*dProbNorm,WalkExcitLevel,iCompositeSize
!               CALL SumEContrib(DetCurr,WalkExcitLevel,iSgn,iLutnI,HDiagCurr,(dProb*dProbNorm))
!HDiags are stored.
!               if(iExcitor.eq.1) THEN
!                  HDiagCurr=CurrentH(j)

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
!                  CALL SumEContrib(DetCurr,WalkExcitLevel,CurrentSign(j),CurrentDets(:,j),HDiagCurr,1.D0)
!               endif


               IF(iDebug.gt.4) then
                  Write(6,'(A,I7)',advance='no') "Killing at excitor: ",iPartDie
                  call WriteBitEx(6,iLutHF,CurrentDets(:,iPartDie),.false.)
                  Write(6,'(A,G25.16)') "Prob: ",dProb*dProbNorm
               endif

   !Die with a probability as normal, but biased because we only selected this set
   !of excitors with dProb.
               if(iCompositeSize.gt.1) dProb=dProb
               iDie=AttemptDieProbPar(DetCurr,HDiagCurr,WalkExcitLevel,iSgn,1/(dProb*dProbNorm))
               IF(iDebug.gt.4.or.((iDebug.eq.3.or.iDebug.eq.4).and.iDie.ne.0)) then
                  Write(6,'(A,I7)',advance='no') " Killing at excitor: ",iPartDie
                  Write(6,'(A)',advance='no') " chosen "
                  call WriteBitEx(6,iLutHF,CurrentDets(:,iPartDie),.false.)
                  WRITE(6,'(A,I7)',advance='no') "Number died ",iDie
                  if(iCompositeSize.gt.1) then
                     WRITE(6,'(A)',advance='no') " from "
                     do i=1,iCompositeSize
                        call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.false.)
                     enddo
                  endif
                  WRITE(6,*)
                  if(iDebug.eq.4) Write(6,'(A,G25.16)') "Prob: ",dProb*dProbNorm
               endif
               NoDied=NoDied+iDie          !Update death counter
               if(iCompositeSize.gt.1.and.iDie.gt.0) nClusterDeaths=nClusterDeaths+1 
               if(iCompositeSize.gt.1.and.iDie.lt.0) nClusterBirths=nClusterBirths+1 
               
               iDie=iDie*iSgn
               if(iDie.ne.0) then
                  iDeaths=iDeaths+1  !The number of deaths we need to modify in the particle list, not the sum of the number that died
                  iKillDetIndices(1,iDeaths)=iDie
                  iKillDetIndices(2,iDeaths)=iPartDie
               endif
!               TotParts=TotParts-abs(CurrentSign(iPartDie))
!               CurrentSign(iPartDie)=CurrentSign(iPartDie)-iDie
!               if(CurrentSign(iPartDie).eq.0) TotRealWalkers=TotRealWalkers-1
   !            TotParts=TotParts+abs(CurrentSign(iPartDie))
!Finish cycling over composites from a particle
            enddo 
!Finish cycling over particles in a det
          enddo     
!Finish cycling over determinants
        enddo
!#if 0
        if(iDebug.gt.2) then
           write(6,*) "nClusters:",nClusters
           write(6,*) "<1/dClusterProb>",dClusterProbs/nClusters
        endif
        dT1Sq=(dT1Sq/HFCount)/HFCount
        if(.not.tCCMCFCI) ENumCyc=ENumCyc+dT1Sq
        dT1SqCuml=dT1SqCuml+dT1Sq
        if (iMaxExcitorSelections.gt.1) then
!            ENumCyc=ENumCyc+dT1Sq
!            IF(Iter.gt.NEquilSteps) SumENum=SumENum+dT1Sq
        endif
!        write(78,'(I,4G)') Iter, dT1Sq,AJWTProjE/HFCount,AJWTProjE,dT1SqCuml/Iter
!#endif

!Now we know what needs to die, we kill it.
        do j=1,iDeaths
            IF(iDebug.gt.4) Write(6,*) "Killing ",iPartDie,iPartDie
            iDie=iKillDetIndices(1,j)
            iPartDie=iKillDetIndices(2,j)
! Actually do this earlier
            CurrentSign(iPartDie)=CurrentSign(iPartDie)-iDie
       enddo
!End the loop over killing dets.
!Deallocate saved loop.
       deallocate(iKillDetIndices)
      IF(iDebug.gt.0) write(6,*) "Total Born: ",NoBorn
      IF(iDebug.gt.0) write(6,*) "Total Died: ",NoDied
      IF(iDebug.gt.0) write(6,*) "Compressing walkers."
!Now traverse the list of walkers, removing walkers which have nothing remaining on them
!VecSlot indicates the next free position in NewDets
       VecSlot=1
       do j=1,TotWalkers
!Do a look over all walkers to accumulate stats, and kill if necessary
!Also, we want to find out the excitation level of the determinant - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.
            IF(tTruncSpace) THEN
!We need to know the exact excitation level for truncated calculations.
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),WalkExcitLevel,NEl)
            ELSE
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),WalkExcitLevel,2)
            ENDIF
!HDiags are stored.
            HDiagCurr=CurrentH(j)
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j))

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info

            iSgn=sign(1,CurrentSign(j))
            iSgn=CurrentSign(j)
!            CALL SumEContrib(DetCurr,WalkExcitLevel,iSgn,CurrentDets(:,j),HDiagCurr,1.d0)
            CopySign=CurrentSign(j)
            IF(CopySign.ne.0.or.WalkExcitLevel.eq.0) THEN
                CurrentDets(:,VecSlot)=CurrentDets(:,j)
                CurrentSign(VecSlot)=CopySign
                CurrentH(VecSlot)=CurrentH(j)
                VecSlot=VecSlot+1
            ENDIF   !To kill if
        enddo
!Since VecSlot holds the next vacant slot in the main determinant array, TotWalkers will be one less than this.
        TotWalkersNew=VecSlot-1
        IF(iDebug.gt.0) write(6,*) "TotNewWalkers ",TotWalkersNew

!***Birth/death processes finished. Tidy up and then annihilate.

!SumWalkersCyc calculates the total number of walkers over an update cycle on each process
        SumWalkersCyc=SumWalkersCyc+(INT(TotParts,i2))
!        WRITE(6,*) "Born, Die: ",NoBorn, NoDied


        write(79,'(5I7)') NoBorn,NoDied,nClusterChildren,nClusterDeaths,nClusterBirths

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
        IF(iDebug.gt.0) WRITE(6,*) "Beginning Annihilation:",TotWalkersNew
        CALL DirectAnnihilation(TotWalkersNew)

        CALL halt_timer(Annihil_Time)
        IF(iDebug.gt.0) WRITE(6,*) "Leaving CCMC Cycle"
        
    END SUBROUTINE PerformCCMCCycPar
#else 
    SUBROUTINE PerformCCMCCycPar()
    END SUBROUTINE PerformCCMCCycPar
#endif


   SUBROUTINE CCMCStandalone(Weight,Energyxw)
      Use global_utilities
      use SystemData, only: nEl,nIfD,nIfTot
      use CCMCData, only: tCCMCFCI,dInitAmplitude
      use DetCalc, only: Det       ! The number of Dets/Excitors in FCIDets
      use DetCalc, only: FCIDets   ! (NIfTot, Det).  Lists all allowed excitors in compressed form
      use DetCalc, only:FCIDetIndex! (nEl).  The index of the different excitation levels
      use CalcData, only: NMCyc    ! The number of MC Cycles
      use CalcData, only: StepsSft ! The number of steps between shift updates
      use FciMCData, only: pDoubles,exFlag,TTruncSpace,Hii,Iter
      use FciMCData, only: TotParts,TotWalkers,TotWalkersOld,TotPartsOld,AllTotPartsOld,AllTotWalkersOld
      use FciMCData, only : HFDet
      use FciMCParMod, only: iLutHF
      use FciMCParMod, only: CheckAllowedTruncSpawn, SetupParameters,BinSearchParts3
      use FciMCParMod, only: CalcNewShift,SumEContrib,InitHistMin
      use FciMCParMod, only: WriteHistogram
      Use Logging, only: CCMCDebug
      USE Logging , only : tHistSpawn,iWriteHistEvery
      USE SymData , only : nSymLabels
      USE Determinants , only : FDet,GetHElement2,GetHElement4,GetHElement3
      USE DetCalc , only : ICILevel,Det,FCIDetIndex
      use CalcData, only: Tau,DiagSft,InitWalkers
      use GenRandSymExcitNUMod , only : GenRandSymExcitScratchNU,ScratchSize
      USE HElem
      USE mt95 , only : genrand_real2
      use DetBitOps, only: DecodeBitDet
      IMPLICIT NONE
      REAL*8 Weight,EnergyxW
      INTEGER , PARAMETER :: r2=kind(0.d0)
      CHARACTER(len=*), PARAMETER :: this_routine='CCMCStandalone'
      INTEGER ierr
      REAL*8, ALLOCATABLE :: Amplitude(:,:) !(Det, 2)
      INTEGER tagAmplitude


      INTEGER nSelects              ! The number of selections of clusters to make
      REAL*8 dProbSelNewExcitor
! The prob that we choose the number of Excitors we have done
      INTEGER iNumExcitors          ! The number of non-zero excitors (excluding the ref det)
      REAL*8 dTotAbsAmpl            ! The total of the absolute amplitudes

      REAL*8 dAmp
      LOGICAL tSuc
      INTEGER PartIndex
      INTEGER iCurAmpList,iOldAmpList
      INTEGER iDebug
      REAL*8  dProb,dProbNorm,r,rat,dCurTot,dProbNumExcit,dClusterProb,dProbDecompose
      INTEGER i,j,k,l
      INTEGER iPartDie
      INTEGER :: SelectedExcitors(0:NIfTot,nEl)     
      INTEGER :: SelectedExcitorIndices(nEl)     
! Temporary Storage
      INTEGER iLutnI(0:nIfTot),iLutnJ(0:nIfTot)
      INTEGER nJ(nEl)
      INTEGER DetCurr(0:nIfTot)
      INTEGER :: ExcitLevel,iGetExcitLevel_2,Ex(2,2)
      INTEGER , ALLOCATABLE :: Scratch1(:),Scratch2(:) !Both (ScratchSize)
      INTEGER iCompositeSize,WalkExcitLevel,IC,iMaxExcit
      REAL*8 dNGenComposite,Prob
      LOGICAL :: tParity,tFilled
      INTEGER iSgn
      TYPE(HElement) rh,HTmp
      REAL*8 HDiagCurr
      REAL*8 dTmp,dPsiTotAbsAmpl,dClusterNorm
      INTEGER iShiftLeft
      REAL*8 WalkerScale
      REAL*8 dT1Sq
      REAL*8 dTolerance
      REAL*8 dAveProbSel(0:4)
      WRITE(6,*) "Entering CCMC Standalone..."
      WRITE(6,*) "pDoubles",pDoubles
      nSelects=InitWalkers
      iDebug=CCMCDebug
      dProbSelNewExcitor=0.7

      dTolerance=1e-16
      iCurAmpList=1  !Start with list 1

! Setup Memory
      Allocate(Amplitude(Det,2),stat=ierr)
      LogAlloc(ierr,'Amplitude',Det*2,8,tagAmplitude)

! Now setup the amplitude list.  Let's start with nothing initially, and
      Amplitude(:,:)=0
!  place ampl 1 in the HF det
      Amplitude(1,iCurAmpList)=dInitAmplitude
      iNumExcitors=0
      dTotAbsAmpl=Amplitude(1,iCurAmpList)

      iShiftLeft=StepsSft
      WalkerScale=10000000/dInitAmplitude
      TotWalkers=WalkerScale*dTotAbsAmpl
      TotParts=WalkerScale*dTotAbsAmpl
      TotWalkersOld=WalkerScale*dTotAbsAmpl
      TotPartsOld=WalkerScale*dTotAbsAmpl
      AllTotWalkersOld=WalkerScale*dTotAbsAmpl
      AllTotPartsOld=WalkerScale*dTotAbsAmpl
      Call SetupParameters()
      Call InitHistMin() !Setup Histogramming arrays if needed 
      ALLOCATE(Scratch1(ScratchSize))
      ALLOCATE(Scratch2(ScratchSize))

! Each cycle we select combinations of excitors randomly, and spawn and birth/die from them
      Iter=1
      do while (Iter.le.NMCyc)
! Copy the old Amp list to the new
         iOldAmpList=iCurAmpList
         iCurAmpList=3-iCurAmpList
         Amplitude(:,iCurAmpList)=Amplitude(:,iOldAmpList)
         IF(iDebug.gt.1) THEN
            write(6,*) "Cycle ",Iter
            write(6,*) "Current, Old Ampl list", iCurAmpList,iOldAmpList
            write(6,*) "HF det"
            call WriteBitDet(6,iLutHF,.true.)
            write(6,*) "Particle list"
            do j=1,Det
               if(iDebug.gt.4.or.abs(Amplitude(j,iCurAmpList)).gt.(dTolerance*dInitAmplitude)) THEN
                  write(6,'(G17.9)',advance='no') Amplitude(j,iCurAmpList)
                  call WriteBitEx(6,iLutHF,FciDets(:,j),.true.)
               ENDIF
   !            do l=0,nIfD
   !               Write(6,'(i)', advance='no') CurrentDets(l,j)
   !            enddo
   !            write(6,*)
            enddo
         endif
         iNumExcitors=0
         dTotAbsAmpl=0
         do j=2,Det
            if(abs(Amplitude(j,iCurAmpList)).gt.dTolerance*dInitAmplitude) then
               iNumExcitors=iNumExcitors+1
               dTotAbsAmpl=dTotAbsAmpl+abs(Amplitude(j,iCurAmpList))
            ENDIF
         enddo
         iC=0
         dT1Sq=0
         do j=1,Det
            do while(j.ge.FCIDetIndex(iC+1).or.FCIDetIndex(iC).eq.FCIDetIndex(iC+1))  !Need to take into account if (e.g.) singles are empty (FCIDI(0:3) = 1 2 2 3, we want j=2 to get to iC=2 not iC=1
               iC=iC+1
            enddo
            CALL DecodeBitDet(DetCurr,FCIDets(:,j))
            i=1
            if(Amplitude(j,iCurAmpList).lt.0) i=-1
!            call WriteBitEx(6,iLutHF,FCIDets(:,j),.false.)
!            WRITE(6,*) "J",j,iC
            dAmp=abs(Amplitude(j,iCurAmpList))/Amplitude(1,iCurAmpList)
            if(dAmp.ne.0.d0)  call SumEContrib(DetCurr,iC,i,FCIDets(:,j),dTmp,1/dAmp)
! Deal with T_1^2
             if(iC.eq.1.and..not.tCCMCFCI) then
                do l=j+1,FCIDetIndex(2)-1
                  iSgn=1
                  dAmp=Amplitude(j,iCurAmpList)*Amplitude(l,iCurAmpList)
                  iLutnI(:)=FCIDets(:,j)
                  call AddBitExcitor(iLutnI,FCIDets(:,l),iLutHF,iSgn)
                  if(iSgn.ne.0.and.dAmp.ne.0.d0) then
 !                    WRITE(6,*) "L",l
                     CALL DecodeBitDet(DetCurr,iLutnI(:))
                     Htmp=GetHElement3(HFDet, DetCurr,2)
                     dAmp=dAmp/(Amplitude(1,iCurAmpList)**2)
                     dT1Sq=dT1Sq+(Real(Htmp%v,r2)*iSgn)*dAmp
!                     WRITE(6,*) dAmp*iSgn
                    call SumEContrib(DetCurr,2,iSgn,iLutnI(:),dTmp,1/dAmp)
!                    WRITE(6,'(A,2G)') 'T1',iSgn*dAmp*real(Htmp%v,r2),dT1Sq
                     !call WriteBitEx(6,iLutHF,CurrentDets(:,j),.false.)
                     !call WriteBitEx(6,iLutHF,CurrentDets(:,l),.false.)
                     !call WriteBitEx(6,iLutHF,iLutnI,.true.)
                  endif
               enddo
            endif
         enddo
         if(iDebug.gt.1) WRITE(6,*) "dET1Sq: ",dT1Sq
         IF(iDebug.gt.1) THEN
            WRITE(6,*) "Total non-zero excitors: ",iNumExcitors
            WRITE(6,*) "Total abs Amplitudes: ",dTotAbsAmpl
         endif 
!  First calculate the energy
         Iter=Iter+1
!  The total cluster wavefunction has the following structure
!
!  Psi_CC = N_0 exp( T / N_0) = N_0 + T + T^2 / (2! N_0) + T^3 / (3! (N_0)^2) + ...
!
!  Instead of the full Psi_CC, we use |_ Psi_CC _| which is truncated at the excitation level corresponding to the number of excitation operators.
!
!  We pick clusters from each level of excitation separately.
!  Were we to use Psi_CC, each level n would have an sum absolute amplitude of N_0 (T/N_0)^n, so the total sum absolute amplitude is N_0 exp(T/N_0).
!  For |_ Psi_CC _| we must truncate this sum manually.
         dPsiTotAbsAmpl=0
         dTmp=Amplitude(1,iOldAmpList)
         do i=0,iNumExcitors
            if(abs(dTmp).lt.1d-8*abs(dPsiTotAbsAmpl)) exit
            dPsiTotAbsAmpl=dPsiTotAbsAmpl+dTmp
            dTmp=dTmp*dTotAbsAmpl/Amplitude(1,iOldAmpList)
            if(i.gt.0) dTmp=dTmp/i
         enddo
         if(iDebug.gt.2) WRITE(6,"(A,I5,A,I5)") "Finished calculating total cluster amplitude at level ",i-1,"/",iNumExcitors
         if(iDebug.gt.1) WRITE(6,*) "Total cluster amplitude: ",dPsiTotAbsAmpl
!         TotParts=dPsiTotAbsAmpl*WalkerScale
         TotParts=(dTotAbsAmpl+abs(Amplitude(1,iCurAmpList)))*WalkerScale
         TotWalkers=TotParts
!  Loop over cluster selections
         dAveProbSel(:)=0
         do j=1,nSelects 
            if(iDebug.gt.3) write(6,*) "Selection ",j
!  Now pick the cluster

! Old comments
! We wish to sum:
! sum_A [ (1/2) sum_{B\ne A} [ (1/3) sum_{C\ne A,B} ... ]]
!  since the sums over B and C are too large a space to do fully, we stochastically sum them:
! sum_A [ (1/2)  1/X sum_i^X 1/p(B_i)  [  (1/3)  1/Y sum_j^Y 1/p(C_j) ... ]]
!   where X and Y are the number of samples taken in the sums.

! If B_i and C_j were selected uniformly randomly out of N-1 particles,  (remember we've already chosen A), 
!  then p(B_i) would be 1/(N-1)  and  p(C_j) would be 1/(N-2)
!  The total prefactor would then be
!  [1/(XY)]  [(N-1) (N-2) / ( 2 . 3 )]     somwehat akin to the combinatorial coefficient.

!  We wish to sample X and Y = 1, but p(B_i) and p(C_j) are not uniform.
!  All in all, we then wish to use the total generation probability of this composite stochastically.
!
!  In this case it would be  [  p(B_i)  p(C_j) / ( 2 . 3 )  ]  p(Choosing this level)

!  In fact, we need the contribution of the composite excitor to vary as (HFcount)^-(iCompositeSize-1)
!    to account for the intermediate normalization.  We make this happen by multiplying the probability of 
!    having generated the excitor by (HFcount)^(iCompositeSize-1)  (i.e. one HFcount for each level)
            dProbNumExcit=dProbSelNewExcitor
            dProbNorm=nSelects
            dClusterProb=1
            dNGenComposite=1
            dClusterNorm=Amplitude(1,iOldAmpList)
            iMaxExcit=min(nEl,iNumExcitors)
!For FCI using amplitudes we only make a single selection
            if(tCCMCFCI) iMaxExcit=min(1,iMaxExcit)
            do i=1,iMaxExcit
   ! Calculate the probability that we've reached this far in the loop
                        !We must have at least one excitor, so we cannot exit here.
               call genrand_real2(r)  !On GHB's advice
               if(r.lt.dProbSelNewExcitor) exit
! Calculate the cluster noormalization for this level:  Used to /i
               dClusterNorm=(dClusterNorm/(i))*dTotAbsAmpl/Amplitude(1,iOldAmpList)
   
! decide not to choose another walker with this prob.
               dProbNumExcit=dProbNumExcit*(1-dProbSelNewExcitor)
! Select a new random walker
               call genrand_real2(r)  !On GHB's advice
               r=r*dTotAbsAmpl
               dCurTot=0
               do k=2,Det
                  dCurTot=dCurTot+abs(Amplitude(k,iOldAmpList))
                  if(dCurTot.ge.r) exit
               enddo
               if(k.gt.Det) THEN
                  write(6,*) "Position in weight list: ", r
                  write(6,*) "Current Total: ",dCurTot
                  write(6,*) "Max Total: ", dTotAbsAmpl
                  write(6,*) Amplitude
                  CALL Stop_All("CCMCStandalone","Invalid Excitor selected")
               endif
               IF(iDebug.gt.5) Write(6,*) "Selected excitor",k
               SelectedExcitorIndices(i)=k
               SelectedExcitors(:,i)=FCIDets(:,k)
               dNGenComposite=dNGenComposite*abs(Amplitude(k,iOldAmpList))  !For each new excit added to the composite, we multiply up to count the number of ways we could've generated it.
               dClusterProb=dClusterProb*abs(Amplitude(k,iOldAmpList))/dTotAbsAmpl  !For each new excit added to the composite, we multiply up to count the number of ways we could've generated it.
!               dProbNorm=dProbNorm/dTotAbsAmpl
            enddo
            IF(iDebug.gt.5) WRITE(6,*) 'prob out of sel routine.',dProbNumExcit
            if(i.gt.iMaxExcit) THEN !We've been limited by the max number of excitations
               ! Let s be dProbSelNewExcitor, and X be nMaxSelExcitors
               !  The sum of all levels from 0 up to X-1 is
               !  s ( 1 - (1-s)^X)/(1-(1-s)) = 1 - (1-s)^X
               !  We take 1-this to be the prob of choosing this level, i.e. (1-s)^X
               !  dProbNumExcit is currently s*(1-s)^X
               dProbNumExcit= dProbNumExcit / dProbSelNewExcitor
            ENDIF
!dProbNumExcit is the prob of having chosen this level
            iCompositeSize=i-1  !Save the number of excitors we've selected   
            if(iCompositeSize.eq.0) then !We haven't correctly calculated dNGenComposite for zero excitors.  Do so now.
               dNGenComposite=dNGenComposite*abs(Amplitude(1,iOldAmpList))  !For each new excit added to the composite, we multiply up to count the number of ways we could've generated it.
            endif
!            WRITE(115,"(30I5)") SelectedExcitorIndices(1:iCompositeSize)

!Up to here dProbNorm merely accounts for the number of selections made
            dProbNorm=dProbNorm*dProbNumExcit
!At This point dProbNorm is the number to divide any contribution from this cluster by to account for its selection.
!  It does not take into account the amplitude weighting of each cluster, so the
!  sum of all  1/dProbNorm = # levels
!  dClusterNorm is the probability that this cluster was chosen, given the level had already been selected.
!  This includes multiple selections of the same excitor as well as combinations of excitors which produce a 0 sign.
            dAveProbSel(iCompositeSize)=dAveProbSel(iCompositeSize)+1/dProbNorm

!Now form the cluster
            IF(iDebug.gt.3) then
               WRITE(6,*) " Excitors in composite:", iCompositeSize
               do i=1,iCompositeSize
                  write(6,'(I2)',advance='no') int(sign(1.d0,Amplitude(SelectedExcitorIndices(i),iOldAmpList)))
                  call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.true.)
               enddo
               Write(6,*) "Level chosen Prob      : ",dProbNumExcit
               Write(6,*) "Select Prob given level: ",dClusterProb
               Write(6,*) "Prob norm              : ",dProbNorm

!The Cluster Normalization for this level (i.e. what the actual magnitude of this level should contribute to the sum.  (TotAbsAmp**iLevel / (iLevel!) ) * N0 / (N_0^iLevel)  
               Write(6,*) "Cluster norm           : ",dClusterNorm
            endif
            iLutnI(:)=iLutHF(:)
            iSgn=1 !The sign of the first excitor
            do i=1,iCompositeSize 
               iSgn=iSgn*int(sign(1.d0,Amplitude(SelectedExcitorIndices(i),iOldAmpList)))
               call AddBitExcitor(iLutnI,SelectedExcitors(:,i),iLutHF,iSgn)
               IF(iDebug.gt.3) Write(6,*) "Results of addition ",i, "Sign ",iSgn,':'
               if(iSgn.eq.0) exit
               IF(iDebug.gt.3) call WriteBitEx(6,iLutHF,iLutnI,.true.)
            enddo
            iF(iDebug.gt.0) CALL FLUSH(6)
            if(iSgn.eq.0) cycle

            IF(iDebug.gt.4) then
               WRITE(6,*) "Chosen det/excitor is:"
               call WriteBitDet(6,iLutnI,.true.)
               CALL FLUSH(6)
            endif

!First, decode the bit-string representation of the determinant the walker is on, into a string of naturally-ordered integers
            CALL DecodeBitDet(DetCurr,iLutnI(:))
            CALL FindBitExcitLevel(iLutHF,iLutnI(:),WalkExcitLevel,nEl)
            if(iDebug.gt.4) WRITE(6,*) "Excitation Level ", WalkExcitLevel

            tFilled=.false.     !This is for regenerating excitations from the same determinant multiple times. There will be a time saving if we can store the excitation generators temporarily.

!Here, we spawn each particle on the determinant in a seperate attempt.
!we are simply looping over all the particles on the determinant
!This will only be a help if most determinants are multiply occupied.

!            call WriteDet(6,DetCurr,nEl,.true.)
            CALL GenRandSymExcitScratchNU(DetCurr,iLutnI,nJ,pDoubles,IC,Ex,tParity,exFlag,Prob,Scratch1,Scratch2,tFilled)
!We need to calculate the bit-representation of this new child. This can be done easily since the ExcitMat is known.
            CALL FindExcitBitDet(iLutnI,iLutnJ,IC,Ex,NIfD)
            IF(iDebug.gt.4) then
                WRITE(6,*) "Random excited det level ",iC
                call WriteDet(6,nJ,nEl,.true.)
                Write(6,*) "Prob ex|from",Prob
            endif
!Prob is the Prob of choosing nJ from nI
!dClusterProb is the Probability of having chosen this cluster excitor (normalized such that <1/dClusterProb> = 1)
!dProbNorm is the renormalization factor for this level of excitors - decreasing by HFcount for each extra level of excitors
!            Prob=Prob*dProbNorm  !/dClusterProb  !Now include the prob of choosing the det we spawned from
            if(iDebug.gt.4) Write(6,*) "Prob ex tot",Prob
            if(iCompositeSize.gt.1) dProb=dProb
!Calculate amplitude to spawn
            IF(TTruncSpace) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.
                IF(CheckAllowedTruncSpawn(WalkExcitLevel,nJ,iLutnJ,IC)) THEN
                  rh=GetHElement4(DetCurr,nJ,IC,Ex,tParity)
                ELSE
                  rh=HElement(0)
                ENDIF
            ELSE
!SD Space is not truncated - allow attempted spawn as usual
                  rh=GetHElement4(DetCurr,nJ,IC,Ex,tParity)
            ENDIF

            rat=-iSgn*Tau*rh%v*dClusterNorm/(Prob*dProbNorm)
            IF(iDebug.gt.4.or.((iDebug.eq.3.or.iDebug.eq.4))) THEN
!We've not printed this out before
               call WriteBitEx(6,iLutHF,iLutnI,.false.)
               write(6,'(A)',advance='no') ,' ==> '
               call WriteBitEx(6,iLutHF,iLutnJ,.false.)
               WRITE(6,'(A,G25.16)',advance='no') "Children:",rat
               if(iDebug.eq.3.and.iCompositeSize.gt.1) THEN
                  write(6,'(A)',advance='no') ' from '
                  do i=1,iCompositeSize
                     call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.false.)
                  enddo
               endif
               write(6,*)
            endif
            if(abs(rat).gt.1e-4*dTolerance*dInitAmplitude) then
!Now add in a contribution from the child
               CALL FindBitExcitLevel(iLutHF,iLutnJ(:),IC,nEl)
               IF(IC.eq.NEl) THEN
                   CALL BinSearchParts3(iLutnJ(:),FCIDets(:,:),Det,FCIDetIndex(IC),Det,PartIndex,tSuc)
               ELSEIF(IC.eq.0) THEN
                   PartIndex=1
                   tSuc=.true.
               ELSE
                   CALL BinSearchParts3(iLutnJ(:),FCIDets(:,:),Det,FCIDetIndex(IC),FCIDetIndex(IC+1)-1,PartIndex,tSuc)
               ENDIF
            
               if(.not.tSuc) THEN      
                  WRITE(6,*) "Cannot find excitor "
                  call WriteBitEx(6,iLutHF,iLutnJ,.true.)
                  call WriteBitDet(6,iLutnJ,.true.)
                  WRITE(6,*) "Excitation Level: ",IC
                  WRITE(6,*) "Dets ",FCIDetIndex(IC),' to ',FCIDetIndex(IC+1)-1
                  call Stop_All("CCMCStandalone","Cannot find excitor in list.")
               endif
               Amplitude(PartIndex,iCurAmpList)=Amplitude(PartIndex,iCurAmpList)+rat
            endif
! Now deal with birth/death.
   ! We have to decompose our composite excitor into one of its parts.  
            IF(iCompositeSize.GT.1) THEN
!  We modify the composite (t_a t_b t_c) -> (t_a t_b t_c - x) by changing just one of the parts
!  t_a -> t_a (1- x/(t_a t_b t_c)).   
               dProbDecompose=1.0/iCompositeSize
               dProbDecompose=.0713
               call genrand_real2(r)  !On GHB's advice
               k=1+floor(r*iCompositeSize)
               iPartDie=SelectedExcitorIndices(k)
            ELSEIF(iCompositeSize.EQ.1) THEN
               dProbDecompose=1
               iPartDie=SelectedExcitorIndices(1)
            ELSE
               dProbDecompose=1
               iPartDie=1
            ENDIF 
            Htmp=GetHElement3(DetCurr,DetCurr,0)
            HDiagCurr=REAL(Htmp%v,r2)
            HDiagCurr=HDiagCurr-Hii
            dProb=dProbDecompose !dClusterProb*dProbDecompose

!Also, we want to find out the excitation level of the determinant - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.
!               IF(tTruncSpace) THEN
!We need to know the exact excitation level for truncated calculations.
!                   CALL FindBitExcitLevel(iLutHF,iLutnI,WalkExcitLevel,NEl)
!               ELSE
!                   CALL FindBitExcitLevel(iLutHF,iLutnI,WalkExcitLevel,2)
!               ENDIF
!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
!               WRITE(6,'(A,3G,2I)') "PPP",dProb,dProbNorm,dProb*dProbNorm,WalkExcitLevel,iCompositeSize
!               CALL SumEContrib(DetCurr,WalkExcitLevel,iSgn,iLutnI,HDiagCurr,(dProb*dProbNorm))
!HDiags are stored.
!               if(iExcitor.eq.1) THEN
!                  HDiagCurr=CurrentH(j)

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
!                  CALL SumEContrib(DetCurr,WalkExcitLevel,CurrentSign(j),CurrentDets(:,j),HDiagCurr,1.D0)
!               endif


            IF(iDebug.gt.4) then
               Write(6,'(A,I7)',advance='no') "Killing at excitor: ",iPartDie
               call WriteBitEx(6,iLutHF,FCIDets(:,iPartDie),.false.)
               Write(6,'(A,G25.16)') "Prob: ",dProb*dProbNorm
            endif

!  This will be the amount we wish to subtract from t_x

!dProb = 1
!            rat=Tau*(HDiagCurr-DiagSft)*dClusterNorm*dProb/dProbNorm !(dProb*dProbNorm)
            rat=Tau*(HDiagCurr-DiagSft)*dProb/(dProbNorm*dClusterProb) !(dProb*dProbNorm)
            IF(iDebug.ge.4) then
               WRITE(6,*) "Product Contributions to Number Died:"
               WRITE(6,*) "Energy difference: ",HDiagCurr-DiagSft
               WRITE(6,*) "Tau              : ",Tau
               WRITE(6,*) "Sign             : ",iSgn
               WRITE(6,*) "Cluster Prob     : ",dClusterProb
               WRITE(6,*) "Cluster Norm     : ",dClusterNorm
               WRITE(6,*) "dProb            : ",dProb
               WRITE(6,*) "1/dProbNorm      : ",1/dProbNorm
               WRITE(6,*) "Delta t_x        : ",-rat
               WRITE(6,*) "dNGenComposite   : ",dNGenComposite
            endif 

!! rat is what we wish to modify t_a t_b t_c by. (but positive - we'll actually want to subtract it)
!!  To do this we modify the chosen part (e.g. t_a) by
!! t_a (1 - rat / (t_a t_b t_c) )
!
!! dNGenComposite = t_a t_b t_c
!            rat= rat/dNGenComposite
!! t_a(new) = t_a(new)+ t_a(old) * rat
            IF(iDebug.ge.4) then
               WRITE(6,*) "Death ratio      : ",rat
            endif
         
            r=Amplitude(iPartDie,iOldAmpList)*rat 
            Amplitude(iPartDie,iCurAmpList)=Amplitude(iPartDie,iCurAmpList)-r
! The two next lines are the ancien regime
!            r=iSgn*Tau*(HDiagCurr-DiagSft)*dClusterNorm/dProbNorm !(dProb*dProbNorm)
!            Amplitude(iPartDie,iCurAmpList)=Amplitude(iPartDie,iCurAmpList)-r
            IF(iDebug.gt.4.or.((iDebug.eq.3.or.iDebug.eq.4))) then
               Write(6,'(A,I7)',advance='no') " Killing at excitor: ",iPartDie
               Write(6,'(A)',advance='no') " chosen "
               call WriteBitEx(6,iLutHF,FCIDets(:,iPartDie),.false.)
               WRITE(6,'(A,G25.16)',advance='no') "Number died ",r
               if(iCompositeSize.gt.1) then
                  WRITE(6,'(A)',advance='no') " from "
                  do i=1,iCompositeSize
                     call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.false.)
                  enddo
               endif
               WRITE(6,*)
               if(iDebug.eq.4) Write(6,'(A,G25.16)') "Prob: ",dProb*dProbNorm
            endif
             

         enddo ! Cluster choices
         if(iDebug.ge.2) WRITE(6,*) "Total Selection Prob (should =1): ",sum(dAveProbSel)
         if(iDebug.ge.2) WRITE(6,*) dAveProbSel
         

! Collate stats
! Calc Shift
         iShiftLeft=iShiftLeft-1

!TotWalkers is used for this and is WalkerScale* total of all amplitudes
         if(iShiftLeft.eq.0)  Call CalcNewShift()
         if(iShiftLeft.le.0)  iShiftLeft=StepsSft

         IF(tHistSpawn.and.(mod(Iter,iWriteHistEvery).eq.0)) THEN
             CALL WriteHistogram()
         ENDIF

      enddo !MC Cycles
      
      LogDealloc(tagAmplitude)
      DeAllocate(Amplitude)
   END SUBROUTINE CCMCStandalone

!Add the excitation in iLutnJ to iLutnI and return it in iLutnI.  iSgn is
!updated with the relevant permutation or set to zero if the excitation is
!disallowed.
SUBROUTINE AddBitExcitor(iLutnI,iLutnJ,iLutRef,iSgn)
   use SystemData, only : nEl,nIfD, NIfTot
   IMPLICIT NONE
   INTEGER iLutnI(0:nIfTot), iLutnJ(0:nIfTot),iLutRef(0:nIfTot)
   INTEGER iLutTmp(0:nIfTot)
   INTEGER T1,T2,T3
   INTEGER iSgn
! We need to run through the bits of J and I concurrently, setting bits of I
   INTEGER i,j

   integer iIlevel,iJlevel,iTmpLevel,iI,iJ
   CALL FindBitExcitLevel(iLutRef,iLutnI,iIlevel,nEl)
   CALL FindBitExcitLevel(iLutRef,iLutnJ,iJlevel,nEl)
   iI=iIlevel
   iJ=iJlevel
   iTmpLevel=0
!   WRITE(6,'(A,I)',advance='no') 'x1,x2,s: ',iSgn
!   call WriteBitEx(6,iLutRef,iLutnI,.false.)
!   call WriteBitEx(6,iLutRef,iLutnJ,.true.)

!We specify the definitions of the excitors in terms of creation and annihilation operators:

!  e.g. t_ij^ab  means   a^+_a a^+_b a_j a_i  (where i<j and a<b).
!  our determinants are specified by an ordered list of occupied orbitals:
!    i,j,k  (i<j<k)
!  This correpsonds to a^+_i a^+_j a^+_k |0>

! Thus an excitor applied to the HF det will yield the postive determinant

!To compose two excitors, we must move the annihilation operators to the right,
! and the creation operators to the left.  Each switch of operators incurs a sign flip.

!First the occ
   do i=0,nIfD
! First mask so we only have 'occupied' orbitals.  The occupieds which are zero
! mean they have been excited from.
      T1=IAND(iLutnI(i),iLutRef(i))
      T2=IAND(iLutnJ(i),iLutRef(i))
      IF(IAND(IEOR(T1,iLutRef(i)),IEOR(T2,iLutRef(i))).NE.0) THEN
!The two excitors are exciting from at least one of the same orbitals.  This
!gives a sign of zero.
!         WRITE(6,'(A,B,B)') 'T1,T2: ',T1,T2
!         WRITE(6,'(A,B)') "Overlapping from",IAND(IEOR(T1,iLutRef(i)),IEOR(T2,iLutRef(i)))
         iSgn=0
         return
      ENDIF
      iLutTmp(i)=IAND(T1,T2)  !Combine the excitors
!We now (very naively, but I can't think of a better way to do it) work out the sign change for the occupied reorder.
!  Tmp lies  between I and J, and will receive the eventual complete excitor.
!  Tmp contains no creation operators as yet, so the right of I can move freely to the left of Tmp
!  J contains a string of iILevel creation operators which must also be taken into account.
! We step through both I and J simultaneously, and move any annihilation operators from either to Tmp,
! with the appropriate sign change.
! iI,iJ and iTmpLevel indicate the number of annihilators left in I J or Tmp.
      T1=IAND(NOT(iLutnI(i)),iLutRef(i))
      T2=IAND(NOT(iLutnJ(i)),iLutRef(i))
      T3=0
      do j=0,31
!      WRITE(6,'(A,B33.32,B33.32,B33.32)') 'T1,T3,T2: ',T1,T3,T2
         if(BTEST(T1,j)) THEN
            iTmpLevel=iTmpLevel+1
            iI=iI-1
            T1=IBCLR(T1,j)
            T3=IBSET(T3,j)
         else if(BTest(T2,j)) then
!Shift the relevant bit from J to tmp.
            if(btest(iJLevel+iJ+iTmpLevel,0)) iSgn=-iSgn
            iJ=iJ-1
            T2=IBCLR(T2,j)
            T3=IBSET(T3,j)
         endif
!         Write(6,*) j,iTmpLevel,iI,iJ,iSgn
      enddo
   enddo
   iI=iILevel
   iJ=iJLevel 
   iTmpLevel=0
!   WRITE(6,*) "Creation"
   do i=nIfD,0,-1  !Go backwards for the sign
! Now mask so we only have 'virtual' orbitals.  The virtuals which are set mean they have been excited to.
      T3=NOT(iLutRef(i))
      T1=IAND(iLutnI(i),T3)
      T2=IAND(iLutnJ(i),T3)
      IF(IAND(T1,T2).NE.0) THEN
!The two excitors are exciting to at least one of the same orbitals.  This
!gives a sign of zero.
!         WRITE(6,'(A,B)') "Overlapping to",IAND(T1,T2)
         iSgn=0
      ENDIF
      iLutnI(i)=IOR(iLutTmp(i),IOR(T1,T2)) 
!Combine the excitors and the excitor from and put them back into where they should go
!We now (very naively, but I can't think of a better way to do it) work out the sign change for the occupied reorder.
!  Tmp lies  between I and J, and will receive the eventual complete excitor.
!  The right of I can still move freely to the left of Tmp
!  J contains a string of iJ creation operators which must also be taken into account,
!   as well as the iJLevel+iILevel annihilation operators and iTmpLevel creation operators in Tmp
!  totalling  (iILevel+iJLevel)+iI+iTmpLeel
! We step through both I and J simultaneously, and move any annihilation operators from either to Tmp,
! with the appropriate sign change.
! iI,iJ and iTmpLevel indicate the number of annihilators left in I J or Tmp.
      T3=0
      do j=31,0,-1
!      WRITE(6,'(A,B33.32,B33.32,B33.32)') 'T1,T3,T2: ',T1,IOR(T3,ieor(iLutTmp(i),iLutRef(i))),T2
         if(BTEST(T1,j)) THEN
            iTmpLevel=iTmpLevel+1
            iI=iI-1
            T1=IBCLR(T1,j)
            T3=IBSET(T3,j)
         else if(BTest(T2,j)) then
!Shift the relevant bit from J to tmp. 
            if(btest((iILevel+iJLevel+iJ)+iTmpLevel,0)) iSgn=-iSgn
            iJ=iJ-1
            T2=IBCLR(T2,j)
            T3=IBSET(T3,j)
         endif
!         Write(6,*) j,iTmpLevel,iI,iJ,iSgn
      enddo
   enddo
!   write(6,*) 'x final sign:',iSgn
   return
END SUBROUTINE AddBitExcitor
