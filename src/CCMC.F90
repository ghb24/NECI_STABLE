MODULE CCMC
    use Determinants, only: get_helement, write_det, write_det_len
    use sort_mod
    use constants, only: dp, int64, n_int, end_n_int,lenof_sign, int32
    use CCMCData, only: ExcitToDetSign,AddBitExcitor
    use ClusterList
    use bit_rep_data, only: NIfDBO,NIfTot
    use bit_reps, only: encode_det
    use FciMCParMod, only: calculate_new_shift_wrapper, iter_data_ccmc
   IMPLICIT NONE
   CONTAINS

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
    INTEGER FUNCTION AttemptDieProbPar(Kii,WSign,dProb)
        use SystemData, only: nEl
        use CalcData , only : DiagSft,Tau
        use FciMCData
        use FciMCParMod, only: TestifDETinCAS
        use DetBitOps, only: FindExcitBitDet, FindBitExcitLevel
        use dSFMT_interface
        IMPLICIT NONE
        INTEGER :: iKill
!        HElement_t :: rh,rhij
        REAL*8 :: r,rat,Kii
        REAL*8 dProb
        integer, dimension(lenof_sign), intent(in) :: wSign

!If there are multiple particles, decide how many to kill in total...
        rat=Tau*(Kii-DiagSft)*abs(WSign(1))*dProb

        iKill=INT(rat)
        rat=rat-REAL(iKill)

!Stochastically choose whether to die or not
        r = genrand_real2_dSFMT() 
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
    SUBROUTINE PerformCCMCCycParInt()
      USe FCIMCParMod
      use CCMCData
      Use Logging, only: CCMCDebug
      use bit_reps, only: decode_bit_det
        IMPLICIT NONE
        INTEGER :: VecSlot,i,j,k,l,CopySign
        INTEGER :: nJ(NEl),IC,DetCurr(NEl)
        INTEGER, DIMENSION(lenof_sign) :: Child
        INTEGER(KIND=n_int) :: iLutnJ(0:NIfTot)
        REAL*8 :: Prob,rat,HDiagCurr,r
        INTEGER :: iDie,WalkExcitLevel,Proc
        INTEGER :: TotWalkersNew,Ex(2,2),Scratch1(ScratchSize),Scratch2(ScratchSize)
        LOGICAL :: tParity,tFilled
        
! We select up to nEl excitors at a time and store them here
        INTEGER(KIND=n_int) :: SelectedExcitors(0:NIfTot,nEl)     
        INTEGER :: SelectedExcitorIndices(nEl)     

! Temporary Storage
        INTEGER(KIND=n_int) :: iLutnI(0:nIfTot)

        ! Unused
        integer :: scratch3(scratchsize)

        ! The sign of the resultant composite
        integer, dimension(lenof_sign) :: iSgn
        


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
        HElement_t Htmp

        INTEGER, allocatable :: iKillDetIndices(:,:)
        INTEGER iDeaths
!Debug print level higher is more
        INTEGER iDebug

        INTEGER TotRealWalkers
        REAL*8 dProbNorm,dClusterProb,dProb
        INTEGER iExcitor
        INTEGER iMaxEx,iMaxExTemp
!iMaxExcitorSelections is the number of times we decide to loop for each particle.
!  it includes the first (exact) loop where we consider the particle by itself.
        INTEGER iMaxExcitorSelections
        LOGICAL tSuccess

        REAL*8 dNGenComposite  ! The number of ways the composite could've been generated.

        REAL*8 dT1Sq
!        REAL*8 AJWTProjE
        INTEGER iCumlExcits,iLeftHere(nEl)
        INTEGER iCurrentCompositeSize
        INTEGER,save :: nClusterBirths
        INTEGER,save :: nClusterDeaths   
        INTEGER,save :: nClusterChildren
        INTEGER nClusters
        REAL*8 dClusterProbs

        !Changes to allow compatibility with the new packaged walkers.
        INTEGER, DIMENSION(lenof_sign) :: TempSign,TempSign2,TempSign3
        nClusters=0
        dClusterProbs=0

        ! Reset counters
        iter_data_ccmc%nborn = 0
        iter_data_ccmc%ndied = 0
        iter_data_ccmc%nannihil = 0

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


        IF(TDebug.and.(mod(Iter,10).eq.0)) THEN
            WRITE(11,*) Iter,TotWalkers,NoatHF,NoatDoubs,MaxIndex,TotParts(1)
            CALL FLUSH(11)
        ENDIF


        IF(tHistSpawn.or.tCalcFCIMCPsi) HistMinInd(1:NEl)=FCIDetIndex(1:NEl)    !This is for the binary search when histogramming
         !This info is destroyed by SumEContrib and needs to be reset each cycle
        
        CALL set_timer(Walker_Time,30)
        IF(iDebug.gt.0) WRITE(6,*) "Number of particles, excitors:",TotParts(1), TotWalkers
        
        iPartBloom=0
!ValidSpawndList now holds the next free position in the newly-spawned list, but for each processor.
        ValidSpawnedList(:)=InitialSpawnedSlots(:)

!We loop over a number of samples of cluster space
! Each sample selects a number of excitors together.

! As the number of walkers in the HF reference det is the normalization, we loop
! over each walker there and use it a number of times
! We take the number of walkers as the number of samples to begin with.
        CALL BinSearchParts(iLutHF, 1, int(TotWalkers,int32), iHFDet,tSuccess)
        if(.not.tSuccess) then
            WRITE(6,*) "WARNING: Cannot find HF det in particle list"
            HFcount=1
            iHFDet=-1
        else
            call extract_sign(CurrentDets(:,iHFDet),TempSign)
            HFcount=abs(TempSign(1))
        endif

        allocate(iKillDetIndices(2,TotParts(1)*2))

        IF(iDebug.gt.1) THEN
         write(6,*) "HF det"
         call WriteBitDet(6,iLutHF,.true.)
         write(6,*) "Particle list"
         do j=1,TotWalkers
            call extract_sign(CurrentDets(:,j),TempSign)
            write(6,'(i7)',advance='no') TempSign(1)
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
          call extract_sign(CurrentDets(:,j),TempSign)
          iSgn(1)=TempSign(1)
          IF(iDebug.gt.4) WRITE(6,*) "Iteration ",Iter,':',j
! Deal with T_1^2
          WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j), 2)
          if(WalkExcitLevel.eq.1) then
            do l=j,TotWalkers
               WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,l), 2)
               if(WalkExcitLevel.eq.1) then
                  call extract_sign(CurrentDets(:,l),TempSign2)
                  iSgn(1)=iSgn(1)*TempSign2(1)
                  iLutnI(:)=CurrentDets(:,j)
                  call AddBitExcitor(iLutnI,CurrentDets(:,l),iLutHF,iSgn(1))
                  if(iSgn(1).ne.0) then
                     call decode_bit_det (DetCurr, iLutnI)
                     Htmp = get_helement (HFDet, DetCurr, iLutHF, iLutnI)
                     dT1Sq=dT1Sq+(Real(Htmp,dp)*iSgn(1))
                     !WRITE(6,'(A,I,2G)', advance='no') 'T1',iSgn,real(Htmp,dp),dT1Sq
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
!            AJWTProjE=AJWTProjE+(Real(Htmp,dp)*CurrentSign(j))
!          endif
!#endif
!Go through all particles on the current walker det
          do l=1,abs(TempSign(1))    !  abs(CurrentSign(j))
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
                  k=TotParts(1)-iCumlExcits  !iCumlExcits includes this det.
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


               call extract_sign(CurrentDets(:,j),TempSign2)
               if(iExcitor.eq.1) then  !Deal with all excitors singly
                  if(iDebug.gt.4) write(6,*) 'Plain old excitor.'
!just select the single excitor
                  iLutnI(:)=CurrentDets(:,j)
                  iSgn(1)=sign(1,TempSign2(1))  ! sign(1,CurrentSign(j))
                  dClusterProb=1 
                  dProbNorm=1
                  iCompositeSize=0
                  call decode_bit_det (DetCurr, iLutnI)
!Also take into account the contributions from the dets in the list
                  HDiagCurr=CurrentH(j)
                  if(tHistSpawn) then
                     WalkExcitLevel = FindBitExcitLevel(iLutHF, iLutnI, nel)
                  else
                     WalkExcitLevel = FindBitExcitLevel(iLutHF, iLutnI, 2)
                  endif
                  CALL SumEContrib(DetCurr,WalkExcitLevel,iSgn(1),iLutnI,HDiagCurr,1.d0)
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
                  dNGenComposite=abs(TempSign2(1))  ! abs(CurrentSign(j))


                  if(tExactCluster) then  !Each time we're here, generate another cluster, up to all of them.
                     if(iExcitor.eq.2) then
!Setup the exact cluster if we need to 
!Pretend we've just finished a single excitor
                        iCurrentCompositeSize=1
                        iLeftHere(1)=abs(TempSign2(1))-1    ! abs(CurrentSign(j))-1
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
                           call extract_sign(CurrentDets(:,SelectedExcitorIndices(i)),TempSign3)
                           iLeftHere(i)=abs(TempSign3(1))-1
!                           iLeftHere(i)=abs(CurrentSign(SelectedExcitorIndices(i)))-1
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
                              call extract_sign(CurrentDets(:,SelectedExcitorIndices(k)),TempSign3)
                              iLeftHere(k)=abs(TempSign3(1))-1
!                              iLeftHere(k)=abs(CurrentSign(SelectedExcitorIndices(k)))-1
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
                        call extract_sign(CurrentDets(:,SelectedExcitorIndices(k)),TempSign3)
                        dNGenComposite=dNGenComposite*abs(TempSign3(1))  !For each new excit added to the composite, we multiply up to count the number of ways we could've generated it.
!                        dNGenComposite=dNGenComposite*abs(CurrentSign(SelectedExcitorIndices(k)))  !For each new excit added to the composite, we multiply up to count the number of ways we could've generated it.
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
                           r = genrand_real2_dSFMT()  !On GHB's advice
                           if(r.lt.dProbSelNewExcitor) exit
                        endif
            
      ! decide not to choose another walker with this prob.
                        dProbNumExcit=dProbNumExcit*dProbSelNewExcitor
      ! Select a new random walker
                        k=0
                        do while (k.eq.0)
                           r = genrand_real2_dSFMT()  !On GHB's advice
                           k=1+floor(r*(TotWalkers-1))    !This selects a unique excitor (which may be multiply populated)
                           if(k.ge.iHFDet) k=k+1          !We're not allowed to select the HF det
                           call extract_sign(CurrentDets(:,k),TempSign3)
                           if(TempSign3(1).eq.0) k=0
                           ! if(CurrentSign(k).eq.0) k=0
                        enddo
                        call extract_sign(CurrentDets(:,k),TempSign3)
                        dNGenComposite=dNGenComposite*abs(TempSign3(1))  !For each new excit added to the composite, we multiply up to count the number of ways we could've generated it.
                        ! dNGenComposite=dNGenComposite*abs(CurrentSign(k))  !For each new excit added to the composite, we multiply up to count the number of ways we could've generated it.
                        IF(iDebug.gt.5) Write(6,*) "Selected excitor",k
                        SelectedExcitorIndices(i)=k
                        SelectedExcitors(:,i)=CurrentDets(:,k)
                        dClusterProb=(dClusterProb*abs(TempSign3(1)))/((TotParts(1)-HFcount))

!Account for intermediate Normalization
                        dProbNorm=dProbNorm*HFCount

!Account for possible orderings of selection
                        dProbNorm=dProbNorm*i
                        IF(iDebug.gt.5) WRITE(6,*) "TotParts,HFCount:",TotParts(1),HFcount
                        IF(iDebug.gt.5) write(6,*) "Prob ",i,": ",(abs(TempSign3(1))+0.d0)/(TotParts(1)-HFcount)," Cuml:", dClusterProb
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
                        call extract_sign(CurrentDets(:,SelectedExcitorIndices(i)),TempSign3)
                        write(6,'(I2)',advance='no') sign(1,TempSign3(1))
                        call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.true.)
                     enddo
                     Write(6,"(A,G25.17)") "   Level chosen Prob      : ",dProbNumExcit
                     Write(6,"(A,G25.17)") "   Select Prob given level: ",dClusterProb
                     Write(6,"(A,G25.17)") "   Prob norm              : ",dProbNorm
                  endif
                  dClusterProb=dProbNumExcit*dClusterProb
! Lets check the probs
                  nClusters=nClusters+1
                  dClusterProbs=dClusterProbs+1/dClusterProb
                  iLutnI(:)=CurrentDets(:,j)
                  call extract_sign(CurrentDets(:,j),TempSign3)
                  iSgn(1)=sign(1,TempSign3(1)) !The sign of the first excitor
                  ! iSgn(1)=sign(1,CurrentSign(j)) !The sign of the first excitor
                  do i=2,iCompositeSize 
                     call extract_sign(CurrentDets(:,SelectedExcitorIndices(i)),TempSign3)
                     iSgn(1)=iSgn(1)*sign(1,TempSign3(1))
                     ! iSgn(1)=iSgn(1)*sign(1,CurrentSign(SelectedExcitorIndices(i)))
                     call AddBitExcitor(iLutnI,SelectedExcitors(:,i),iLutHF,iSgn(1))
                     IF(iDebug.gt.3) Write(6,*) "Results of addition ",i, "Sign ",iSgn(1),':'
                     if(iSgn(1).eq.0) exit
                     IF(iDebug.gt.3) call WriteBitEx(6,iLutHF,iLutnI,.true.)
                  enddo
                  IF(iDebug.gt.0) CALL FLUSH(6)
                  if(iSgn(1).eq.0) cycle
               endif

               IF(iSgn(1)/=0.and.iDebug.gt.4) then
                  WRITE(6,*) "Chosen det/excitor is:"
                  WRITE(6,"(A)",advance="no") "  "
                  call WriteBitDet(6,iLutnI,.true.)
                  CALL FLUSH(6)
               endif

               ! First, decode the bit-string representation of the 
               ! determinant the walker is on, into a string of 
               ! naturally-ordered integers
               call decode_bit_det (DetCurr, iLutnI)
               WalkExcitLevel = FindBitExcitLevel(iLutHF, iLutnI, nel)
               if(iDebug.gt.4) WRITE(6,*) "Excitation Level ", WalkExcitLevel



               tFilled=.false.     !This is for regenerating excitations from the same determinant multiple times. There will be a time saving if we can store the excitation generators temporarily.

!Here, we spawn each particle on the determinant in a seperate attempt.
!we are simply looping over all the particles on the determinant
!This will only be a help if most determinants are multiply occupied.

               call gen_rand_excit (DetCurr, iLutnI, nJ, iLutnJ, exFlag, IC, &
                                    Ex, tParity, Prob, tFilled, Scratch1, &
                                    Scratch2, Scratch3)
               if(.not.IsNullDet(nJ)) then  !Check it hasn't given us a null determinant as it couldn't find one in a sensible time.
!We need to calculate the bit-representation of this new child. This can be done easily since the ExcitMat is known.
                  IF(.not.tHPHF) CALL FindExcitBitDet(iLutnI,iLutnJ,IC,Ex)
                  IF(iDebug.gt.4) then
                      WRITE(6,*) "  Random excited det level ",iC
                      WRITE(6,"(A)",advance="no") "   "
                      call write_det (6, nJ, .true.)
                      Write(6,*) "  Prob ex|from",Prob
                  endif
!Prob is the Prob of choosing nJ from nI
!dClusterProb is the Probability of having chosen this cluster excitor (normalized such that <1/dClusterProb> = 1)
!dProbNorm is the renormalization factor for this level of excitors - decreasing by HFcount for each extra level of excitors
                  Prob=Prob*dClusterProb*dProbNorm  !Now include the prob of choosing the det we spawned from
                  if(iDebug.gt.4) Write(6,*) "  Prob ex tot",Prob
                  if(iCompositeSize.gt.1) dProb=dProb
!Calculate number of children to spawn
                  IF(TTruncSpace) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.
                      IF(CheckAllowedTruncSpawn(WalkExcitLevel,nJ,iLutnJ,IC)) THEN
                          Child=AttemptCreatePar(DetCurr,iLutnI,iSgn,nJ,iLutnJ,Prob,IC,Ex,tParity)
                      ELSE
                          Child=0
                      ENDIF
                  ELSE
!SD Space is not truncated - allow attempted spawn as usual
                      Child=AttemptCreatePar(DetCurr,iLutnI,iSgn,nJ,iLutnJ,Prob,IC,Ex,tParity)
                  ENDIF

                  IF(iDebug.gt.4.or.((iDebug.eq.3.or.iDebug.eq.4).and.Child(1).ne.0)) THEN
   !We've not printed this out before
                     write(6,"(A)",advance="no") "    "
                     call WriteBitEx(6,iLutHF,iLutnI,.false.)
                     write(6,'(A)',advance='no') ' ==> '
                     call WriteBitEx(6,iLutHF,iLutnJ,.false.)
                     WRITE(6,'(A,I7)',advance='no') "Children:",Child(1)             
                     if(iDebug.eq.3.and.iCompositeSize.gt.1) THEN
                        write(6,'(A)',advance='no') ' from '
                        do i=1,iCompositeSize
                           call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.false.)
                        enddo
                     endif
                     write(6,*)
                  endif

                  IF(Child(1).ne.0) THEN
                     if(iCompositeSize.gt.1) nClusterChildren=nClusterChildren+1 
   !We want to spawn a child - find its information to store

   !                    WRITE(6,*) "Spawning particle to:",nJ(:)
   !                    ExcitLevel=iGetExcitLevel_2(HFDet,nJ,NEl,NEl)
   !                    WRITE(6,*) "Excitlevel:", ExcitLevel
                       NoBorn=NoBorn+abs(Child(1))     !Update counter about particle birth
                       iter_data_ccmc%nborn(1) = iter_data_ccmc%nborn(1) + abs(Child(1))
                       IF(IC.eq.1) THEN
                           SpawnFromSing=SpawnFromSing+abs(Child(1))
                       ENDIF

                       IF(abs(Child(1)).gt.25) THEN
   !If more than 25 particles are created in one go, then log this fact and print out later that this has happened.
                           IF(abs(Child(1)).gt.abs(iPartBloom)) THEN
                               IF(IC.eq.1) THEN
                                   iPartBloom=-abs(Child(1))
                               ELSE
                                   iPartBloom=abs(Child(1))
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
                       call encode_det(SpawnedParts(:,ValidSpawnedList(Proc)),iLutnJ)
                       call encode_sign(SpawnedParts(:,ValidSpawnedList(Proc)),Child)
                       ! SpawnedParts(:,ValidSpawnedList(Proc))=iLutnJ(:)
                       ! SpawnedSign(ValidSpawnedList(Proc))=Child
                       ValidSpawnedList(Proc)=ValidSpawnedList(Proc)+1

                       Acceptances=Acceptances+ABS(Child(1))      !Sum the number of created children to use in acceptance ratio
                
                  ENDIF   !End if child created

               ENDIF !.not.IsNullDet(nJ)

   ! We have to decompose our composite excitor into one of its parts.  
               IF(iCompositeSize.GT.0) THEN
                  dProbDecompose=dNGenComposite
!                  dProb=dProb/HFcount
                  r = genrand_real2_dSFMT()  !On GHB's advice
                  k=1+floor(r*iCompositeSize)
                  iPartDie=SelectedExcitorIndices(k)

                  ! Now get the full representation of the dying excitor
                  call decode_bit_det (DetCurr, iLutnI)
                  Htmp = get_helement (DetCurr, DetCurr, 0)
                  HDiagCurr=REAL(Htmp,dp)
                  HDiagCurr=HDiagCurr-Hii
                  call extract_sign(CurrentDets(:,iPartDie),TempSign3)
                  iSgn=sign(1,TempSign3(1))
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
                   WalkExcitLevel = FindBitExcitLevel(iLutHF, iLutnI, nel)
               ELSE
                   WalkExcitLevel = FindBitExcitLevel(iLutHF, iLutnI, 2)
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
               iDie=AttemptDieProbPar(HDiagCurr,iSgn,1/(dProb*dProbNorm))
               IF(iDebug.gt.4.or.(iDebug.eq.4.and.iDie.ne.0)) then
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
               iter_data_ccmc%ndied(1) = iter_data_ccmc%ndied(1) + iDie
               if(iCompositeSize.gt.1.and.iDie.gt.0) nClusterDeaths=nClusterDeaths+1 
               if(iCompositeSize.gt.1.and.iDie.lt.0) nClusterBirths=nClusterBirths+1 
               
               iDie=iDie*iSgn(1)
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

            call extract_sign(CurrentDets(:,iPartDie),TempSign3)
            TempSign3(1)=TempSign3(1)-iDie
            call encode_sign(CurrentDets(:,iPartDie),TempSign3)
            ! CurrentSign(iPartDie)=CurrentSign(iPartDie)-iDie
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
                WalkExcitLevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j),&
                                                   nel)
            ELSE
                WalkExcitLevel = FindBitExcitLevel(iLutHF,CurrentDets(:,j), 2)
            ENDIF

            ! HDiags are stored.
            HDiagCurr=CurrentH(j)
            call decode_bit_det (DetCurr, CurrentDets(:,j))

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info

            call extract_sign(CurrentDets(:,j),TempSign3)
            iSgn(1)=TempSign3(1)
            !   iSgn(1)=CurrentSign(j)
!            CALL SumEContrib(DetCurr,WalkExcitLevel,iSgn,CurrentDets(:,j),HDiagCurr,1.d0)
            ! CopySign=CurrentSign(j)
            CopySign=TempSign3(1)
            IF(CopySign.ne.0.or.WalkExcitLevel.eq.0) THEN
                call encode_det(CurrentDets(:,VecSlot),CurrentDets(:,j))
                call encode_sign(CurrentDets(:,VecSlot),TempSign3)
                ! CurrentDets(:,VecSlot)=CurrentDets(:,j)
                ! CurrentSign(VecSlot)=CopySign
                CurrentH(VecSlot)=CurrentH(j)
                VecSlot=VecSlot+1
            ENDIF   !To kill if
        enddo
!Since VecSlot holds the next vacant slot in the main determinant array, TotWalkers will be one less than this.
        TotWalkersNew=VecSlot-1
        IF(iDebug.gt.0) write(6,*) "TotNewWalkers ",TotWalkersNew

!***Birth/death processes finished. Tidy up and then annihilate.

!SumWalkersCyc calculates the total number of walkers over an update cycle on each process
        SumWalkersCyc=SumWalkersCyc+(INT(TotParts(1),int64))
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
        CALL DirectAnnihilation(TotWalkersNew, iter_data_ccmc)

        CALL halt_timer(Annihil_Time)
        IF(iDebug.gt.0) WRITE(6,*) "Leaving CCMC Cycle"
        
        ! Update counters
        iter_data_ccmc%update_growth = iter_data_ccmc%update_growth + iter_data_ccmc%nborn - iter_data_ccmc%ndied - iter_data_ccmc%nannihil
        iter_data_ccmc%update_iters = iter_data_ccmc%update_iters + 1
        
    END SUBROUTINE PerformCCMCCycParInt
#else 
    SUBROUTINE PerformCCMCCycParInt()
    END SUBROUTINE PerformCCMCCycParInt
#endif



! We create a transition matrix where each element correpsonds to a cluster.  We encode cluster X=(a_{X_i}) = (x,y,z)  as an 'bit' string in base Det-1 sum_{i=1}^{|X|} (X_i-1)*(Det-1)**(i-1)
   INTEGER FUNCTION GetClusterIndex(Clust,iSize,TL)
      use DetCalcData, only: Det       ! The number of Dets/Excitors in FCIDets
      Use CCMCData, only: CCTransitionLog
      IMPLICIT NONE
      INTEGER iSize,Clust(iSize)
      INTEGER i
      INTEGER ind
      INTEGER base
      Type(CCTransitionLog) TL
      if(.not.TL%tNonUniq) then
         GetClusterIndex=GetUniqClusterIndex(Clust,iSize,TL)
         return
      endif
      ind=0
      base=1
      do I=1,iSize
         ind=ind+base*(Clust(i)-1)
         base=base*Det
      enddo
      GetClusterIndex=ind
!      WRITE(6,"(I3,A,20I3)") ind," ClIn",Clust
   END FUNCTION GetClusterIndex

!For Uniq clusters - i.e. where the indices are ordered, there's a compact indexing scheme
! (a1,a2,a3,...,aD) 1<=a1<a2<a3<...<aD<=N
!
!Through some combinatorial identities (see AJWT's 28/1/10 notes) we can create a unique index with

! \sum_{i=1}^D Binomial(N-a_{i-1}, D-i+1) - Binomial(N-a_i +1, D-i+1) + 1
! This runs from 1... Binomial(N, D), where N is the total number of allowed excitors,
! and D is the number of excitors in the cluster (iSize)

! Since we are (also) allowed all sizes from i=0 ... iSize-1, we add in sum_i=0^iSize-1 Binomial(N,i) too

!Heaven help you if you give it an unordered cluster.

!NB Clust's indices are 2-based (i.e. the first 1-, 2-, and 3- clusters are : [2],[2,3],[2,3,4])
   INTEGER FUNCTION GetUniqClusterIndex(Clust,iSize,TL)
      use CCMCData, only: CCTransitionLog
      IMPLICIT NONE
      INTEGER iSize,Clust(iSize)
      INTEGER i
      INTEGER ind
      TYPE(CCTransitionLog) TL

      INTEGER Binomial
      ind=0
      do I=0,iSize-1
         ind=ind+Binomial(TL%nExcitors,i)
      enddo
      do I=1,iSize
         if(i>1) then
            ind=ind+Binomial(TL%nExcitors-(Clust(i-1)-1),iSize-i+1)
         else
            ind=ind+Binomial(TL%nExcitors,iSize-i+1)
         endif
         ind=ind-Binomial(TL%nExcitors-(Clust(i)-1)+1,iSize-i+1)
      enddo
      GetUniqClusterIndex=ind
   END FUNCTION GetUniqClusterIndex

! We create a transition matrix where each element correpsonds to a cluster.  We encode cluster X=(a_{X_i}) = (x,y,z)  as an 'bit' string in base Det-1 sum_{i=1}^{|X|} (X_i-1)*(Det-1)**(i-1)
   INTEGER FUNCTION GetClusterIndLevel(oind)
      use DetCalcData, only: Det       ! The number of Dets/Excitors in FCIDets
      IMPLICIT NONE
      INTEGER i
      INTEGER ind,oind
      ind=oind
      i=0
      do while(ind.ne.0)
         i=i+1
         ind=ind/(Det)
      enddo
      GetClusterIndLevel=i
!      WRITE(6,"(I3,A,20I3)") ind," ClIn",Clust
   END FUNCTION GetClusterIndLevel

!Write a Cluster index, decompressing it into its individual excitors.
! We create a transition matrix where each element correpsonds to a cluster.  We encode cluster X=(a_{X_i}) = (x,y,z)  as an 'bit' string in base Det sum_{i=1}^{|X|} (X_i)*(Det)**(i-1)
   SUBROUTINE WriteClusterInd(iUnit,ind0,lTerm,TL,tNonUniq)
      use DetCalcData, only: Det,FciDets       ! The number of Dets/Excitors in FCIDets
      use FCIMCParMod, only: iLutHF
      use CCMCData, only: Cluster,CCTransitionLog,WriteCluster
      IMPLICIT NONE
      INTEGER ind,ind0,i
      LOGICAL lTerm
      INTEGER iUnit
      TYPE(Cluster) C
      TYPE(CCTransitionLog) TL
      LOGICAL, optional :: tNonUniq
      LOGICAL tNonU
      if(present(tNonUniq)) then
         tNonU=tNonUniq
      else
         tNonU=TL%tNonUniq
      endif
         
      if(.not.tNonU) then
         call InitCluster(C)
         call UnpackUniqClusterInd(C,ind0,TL)
         do i=1,C%iSize
            C%SelectedExcitors(:,i)=FCIDets(:,C%SelectedExcitorIndices(i))
         enddo
         call WriteCluster(iUnit,C,lTerm)
         return
      endif
      ind=ind0
      WRITE(iUnit,'(A)',advance='no') '['
      write(iUnit,'(I4)',advance='no') ind0
      do while(ind.ne.0)
         i=mod(ind,Det)+1
         call WriteBitEx(iUnit,iLutHF,FCIDets(:,i),.false.)
         write(iUnit,'(A)',advance='no') ','
         ind=ind/(Det)
      enddo
      WRITE(iUnit,'(A)',advance='no') ']'
      if(lTerm) WRITE(iUnit,*)
   END SUBROUTINE WriteClusterInd

!Unpack a uniq cluster index.  The size of the cluster is indicated by the index, and clusters are indexed consecutively with
!  size 0, followed by all size 1, followed by all size 2...
   subroutine UnpackUniqClusterInd(C,ind0,TL)
      use CCMCData, only: Cluster,CCTransitionLog
      TYPE(Cluster) C
      TYPE(CCTransitionLog) TL
      INTEGER ind0,ind,start,i,loc,N

      INTEGER Binomial
      ind=ind0
      start=1
      C%iSize=0
      N=TL%nExcitors
!      write(6,*) "Ind0,N",ind0,N
      do while (ind.ge.start)
         ind=ind-start
         C%iSize=C%iSize+1
         start=(start*(N-C%iSize+1))/C%iSize  !Create Binomial(nExcitors, iSize)
!         write(6,*) "ind,start,size",ind,start,C%iSize
      enddo
!      write(6,*) "iSize",C%iSize
!ind is now the index within the set of size iSize.
! sequentially subtract off layers until we've unpacked.  This is horribly inefficient.
      loc=1 !loc is the index next to unpack
      do while (loc.le.C%iSize)
         if(loc.eq.1) then
            C%SelectedExcitorIndices(loc)=0
         else
            C%SelectedExcitorIndices(loc)=C%SelectedExcitorIndices(loc-1)
         endif
         start=-1  !Start will be the start of the set of indices for C%SelectedExcitorIndices(loc)+1
         do while (ind.ge.start)  !If our index is in that set (or beyond), we keep searching.
            if(start.gt.0) ind=ind-start  !Remove another layer
            C%SelectedExcitorIndices(loc)=C%SelectedExcitorIndices(loc)+1
            ! Number of ways of choosing the remaining number of excitors in the cluster
            ! from the remaining number of total excitors (in an ordered fashion), if we
            ! were to have the next value of C%SelectedExcitorIndices(loc)
            start=Binomial(N-C%SelectedExcitorIndices(loc),C%iSize-loc) 
!            write(6,*) "ind,start,loc",ind,start,loc,C%SelectedExcitorIndices(loc)
         enddo
!         WRITE(6,*) "Position ",loc,": ",C%SelectedExcitorIndices(loc)
         loc=loc+1
      enddo
      do i=1,C%iSize  !We've created indices of excitors, which need to be shifted by 1 to give indices in the Dets list
            C%SelectedExcitorIndices(i)=C%SelectedExcitorIndices(i)+1
      enddo
   end subroutine UnpackUniqClusterInd

!This function will return a new determinant on which to spawn if we haven't gone beyond the number which we're meant to spawn.
!False will be returned if we're done.
!True will be returned if an attempt was made to generate a spanwer.  It is poassible that with random spawning it was not poassible to 
! feasibly generate a spawnee.  a null determinant is returned in this case (which can be tested for by the function IsNullDet)

LOGICAL FUNCTION GetNextSpawner(S,iDebug)
   use CCMCData
   use FciMCData, only: pDoubles,tTruncSpace
   use FciMCParMod, only: CheckAllowedTruncSpawn
   use SystemData, only : nEl
   use GenRandSymExcitNUMod , only : gen_rand_excit, scratchsize
   use SymExcit3, only: GenExcitations3
   use DetBitOps, only: FindExcitBitDet
   IMPLICIT NONE
   TYPE(Spawner) S
   INTEGER iDebug

   LOGICAL tFilled
   LOGICAL tParity
   LOGICAL tDone

    ! unused
    integer(kind=n_int) :: iLutnJ(0:nifTot)
    integer :: scratch3(scratchsize)

   tDone=.false.
   S%iIndex=S%iIndex+1
   S%bValid=.false.
   if(.not.S%tFull) THEN
      if(S%iIndex.gt.S%nSpawnings) THEN
         GetNextSpawner=.false.
      else
         tFilled=S%iIndex.gt.1     !This is for regenerating excitations from the same determinant multiple times. There will be a time saving if we can store the excitation generators temporarily.
         call gen_rand_excit (S%C%DetCurr, S%C%iLutDetCurr, S%nJ, iLutnJ, &
                              S%exFlag, S%iExcitLevel, S%ExcitMat, tParity, &
                              S%dProbSpawn, tFilled, S%Scratch1, S%Scratch2, &
                              scratch3)
         GetNextSpawner=.true.
         S%dProbSpawn=S%dProbSpawn*S%nSpawnings
      endif
   ELSE
!      WRITE(6,*) tDone,S%dProbSpawn
!      write(6,*) S%ExcitMat,tParity
!      Write(6,*) "Getting Excitations"
      CALL GenExcitations3(S%C%DetCurr,S%C%iLutDetCurr,S%nJ,S%exFlag,S%ExcitMat,tParity,tDone)
!      call WriteDet(6,S%nJ,nEl,.false.)
!      WRITE(6,*) tDone
      if(S%ExcitMat(1,2).eq.0) then
         S%iExcitLevel=1
      else
         S%iExcitLevel=2
      endif 
!      WRITE(6,*) tDone,S%dProbSpawn
!      write(6,*) S%ExcitMat,tParity
      GetNextSpawner=.not.tDone
      if(tDone) S%nJ(1)=0
   ENDIF
   if(GetNextSpawner.and.IsNullDet(S%nJ)) then
      if(iDebug.gt.4) WRITE(6,*) " GetNextSpawner",S%iIndex
      if(iDebug.gt.4) WRITE(6,*) "  Excitation not found"
   else if(GetNextSpawner) then  !Check it hasn't given us a null determinant as it couldn't find one in a sensible time.
!We need to calculate the bit-representation of this new child. This can be done easily since the ExcitMat is known.
      S%bValid=.true.
      if(iDebug.gt.4) WRITE(6,*) " GetNextSpawner",S%iIndex
      CALL FindExcitBitDet(S%C%iLutDetCurr,S%iLutnJ,S%iExcitLevel,S%ExcitMat)
      IF(iDebug.gt.4) then
          WRITE(6,*) "  Random excited det level ",S%iExcitLevel
          write(6,"(A)",advance="no") "   "
          call write_det(6, S%nJ, .true.)
          Write(6,*) "  Prob ex|from",S%dProbSpawn
      endif
!Prob is the Prob of choosing nJ from nI
!dClusterProb is the Probability of having chosen this cluster excitor (normalized such that <1/dClusterProb> = 1)
      if(iDebug.gt.4) Write(6,*) "  Prob ex tot",S%dProbSpawn
!Calculate amplitude to spawn
      IF(TTruncSpace) THEN
!We have truncated the excitation space at a given excitation level. See if the spawn should be allowed.
          IF(CheckAllowedTruncSpawn(S%C%iExcitLevel,S%nJ,S%iLutnJ,S%iExcitLevel)) THEN
            S%HIJ = get_helement (S%C%DetCurr, S%nJ, S%iExcitLevel, &
                                  S%ExcitMat, tParity)
          ELSE
            S%HIJ=(0)
          ENDIF
      ELSE
!SD Space is not truncated - allow attempted spawn as usual
            S%HIJ = get_helement (S%C%DetCurr, S%nJ, S%iExcitLevel, &
                                        S%ExcitMat, tParity)
      ENDIF
   ENDIF
END FUNCTION GetNextSpawner


SUBROUTINE InitCluster(C)
   use CCMCData
   use SystemData, only: nEl
   IMPLICIT NONE
   TYPE(Cluster) C
   allocate(C%SelectedExcitors(0:NIfTot,nEl))
   allocate(C%SelectedExcitorIndices(nEl))
   allocate(C%iLutDetCurr(0:NIfTot))
   allocate(C%DetCurr(nEl))
END SUBROUTINE InitCluster


SUBROUTINE DelCluster(C)
   use CCMCData
   IMPLICIT NONE
   TYPE(Cluster) C
   deallocate(C%SelectedExcitors)
   deallocate(C%SelectedExcitorIndices)
   deallocate(C%iLutDetCurr)
   deallocate(C%DetCurr)
end subroutine DelCluster

!iMaxSize is the maximum number of excitors in a cluster.
SUBROUTINE InitClustSelectorFull(CS,iMaxSize)
   use CCMCData
   IMPLICIT NONE
   TYPE(ClustSelector) CS
   INTEGER iMaxSize
   CS%tFull=.true.
   CS%iMaxSize=iMaxSize
   CS%tDynamic=.false.
   Call InitCluster(CS%C)
END SUBROUTINE InitClustSelectorFull

SUBROUTINE InitClustSelectorRandom(CS,iMaxSize,nSelects,dProbSelNewEx)
   use CCMCData
   IMPLICIT NONE
   TYPE(ClustSelector) CS
   INTEGER iMaxSize,nSelects
   REAL*8 dProbSelNewEx
   if(nSelects<0) then
      CS%tDynamic=.true.
   else
      CS%tDynamic=.false.
   endif
   CS%tFull=.false.
   CS%iMaxSize=iMaxSize
   CS%nSelects=nSelects
   CS%dProbSelNewExcitor=dProbSelNewEx
   Call InitCluster(CS%C)
END SUBROUTINE InitClustSelectorRandom

SUBROUTINE ResetClustSelector(CS,iRefPos)
   use CCMCData
   IMPLICIT NONE
   TYPE(ClustSelector) CS
   INTEGER iRefPos
   CS%iIndex=0
   CS%C%iSize=0
   CS%iRefPos=iRefPos
END SUBROUTINE ResetClustSelector

!Takes an ordered tuple of length iSize, and gives the next one in sequence.
!iMin..iMax are the max extent of the values in the tuple.
!If Tuple(1) is 0 then it initializes the tuple.
!Afterwards, if tDone is set then it has run out of tuples.
SUBROUTINE IncrementOrderedTuple(Tuple,iSize,iMin,iMax,tDone)
   IMPLICIT NONE
   INTEGER Tuple(iSize),iSize,iMax,iMin
   LOGICAL tDone
   INTEGER i
   if(iSize.eq.0) then
      tDone=.true.
      return
   endif
   if(Tuple(1).lt.iMin) then
      i=1
   else
      i=iSize
   endif
! i is the index in the tuple we're currently trying to increment.
   do while (i.le.iSize)
      Tuple(i)=Tuple(i)+1
      if(Tuple(i).gt.(iMax-(iSize-i))) then
!If we've gone beyond what is the max possible value for this slot, then we move back, otherwise we move forward
         i=i-1
         if(i.eq.0) then
            tDone=.true.
            return
         endif
      else
         i=i+1
         if(i.le.iSize) Tuple(i)=Tuple(i-1)
      endif
   enddo
   tDone=.false.
   return
END SUBROUTINE IncrementOrderedTuple
   
   

!Reset a Spawner to use a new (collapsed) Cluster. This takes a pointer to the cluster for reference purposes, so the cluster and any structure it is derived from must have a TARGET attribute.
SUBROUTINE ResetSpawner(S,C,nSpawn)
   use SymExcit3, only: CountExcitations3
   use CCMCData
   use SystemData, only: nel
   use FciMCData, only: exFlag
   IMPLICIT NONE
   TYPE(Spawner) S
   INTEGER nSpawn
   Type(Cluster),target :: C
   S%C=>C
   S%iIndex=0
   S%ExcitMat(:,:)=0
   S%nSpawnings=nSpawn
   if(S%tFull) then
!      call CountExcitations3(C%DetCurr,exFlag,nS,nD)
      S%dProbSpawn=1.D0
!      S%dProbSpawn=1.D0/(nD+nS)
   endif
   S%exFlag=3
END SUBROUTINE ResetSpawner

!iMaxExcitLevel is the maximum excitation level in an excitor
SUBROUTINE InitSpawner(S,tFull,iMaxExcitLevel)
   use CCMCData, only: Spawner
   use GenRandSymExcitNUMod , only : ScratchSize
   use SystemData, only: nel
   IMPLICIT NONE
   TYPE(Spawner) S
   LOGICAL tFull
   INTEGER iMaxExcitLevel
   S%tFull=tFull
   S%iMaxExcitLevel=iMaxExcitLevel
   if(.not.tFull) then
      allocate(S%Scratch1(ScratchSize))
      allocate(S%Scratch2(ScratchSize))
   endif
   allocate(S%nJ(nEl))
   allocate(S%iLutnJ(0:nIfTot))
END SUBROUTINE InitSpawner



!This runs over all singles and doubles in the list of excitors and initiates according to the MP1 wavefunction

! tFCI                         If set, then use (1+T)|D_0> rather than (1+T+T^2 /2)|D0> as the wavefunction
! Amplitude(nExcit)            is the amplitude of each excitor
! nExcit                       is the number of excitors
! ExcitList(0:nIfTot,nExcit)   contains the bit-compressed list of excitors
! ExcitLevelIndex(0:nEl+1)     is the index of the first det of each excitation level in ExcitList
! HFDet(nEl)                   is the reference determinant on which the excitors are based
SUBROUTINE InitMP1Amplitude(tFCI,Amplitude,nExcit,ExcitList,ExcitLevelIndex,dInitAmp,dTotAbsAmpl)
   use CCMCData
   use SystemData, only: nel
   use FciMCData, only: HFDet
   use FciMCParMod, only: iLutHF,SumEContrib,BinSearchParts3
   use Determinants, only: GetH0Element3
   use bit_reps, only: decode_bit_det
   use constants, only: dp
   IMPLICIT NONE
   LOGICAL tFCI
   REAL*8 Amplitude(nExcit)
   INTEGER nExcit
   INTEGER(KIND=n_int) ExcitList(0:nIfTot,nExcit)
   INTEGER ExcitLevelIndex(0:nEl+1)
   REAL*8 dInitAmp,dTotAbsAmpl
   INTEGER iC,j,l,iSgn
   REAL*8 dAmp
   INTEGER DetCurr(nEl)
   HElement_t HTmp,H0Tmp,H0HF
   INTEGER(KIND=n_int) iLutnI(0:nIfTot)
   INTEGER PartIndex
   LOGICAL tSuc

   iC=0
   H0HF=GetH0Element3(HFDet)
   Amplitude(:)=0
   Amplitude(1)=dInitAmp
   dTotAbsAmpl=0
   do j=1,nExcit
      do while(j.ge.ExcitLevelIndex(iC+1).or.ExcitLevelIndex(iC).eq.ExcitLevelIndex(iC+1))  !Need to take into account if (e.g.) singles are empty (FCIDI(0:3) = 1 2 2 3, we want j=2 to get to iC=2 not iC=1
         iC=iC+1
      enddo
      call decode_bit_det (DetCurr, ExcitList(:,j))
      if(iC.ge.1) then
         Htmp = get_helement (HFDet,  DetCurr, iC, iLutHF, ExcitList(:,j))
         H0tmp=GetH0Element3(DetCurr)
         H0tmp=H0tmp-H0HF
         Amplitude(j)=Amplitude(j)-dInitAmp*(Htmp)/(H0tmp)
         if(iC.eq.1.and..not.tFCI) then
            do l=ExcitLevelIndex(1),j-1
               iSgn=1
               dAmp=Amplitude(j)*Amplitude(l)
               iLutnI(:)=ExcitList(:,j)
               call AddBitExcitor(iLutnI,ExcitList(:,l),iLutHF,iSgn)
               if(iSgn.ne.0.and.dAmp.ne.0.d0) then
                  CALL BinSearchParts3(iLutnI,ExcitList,nExcit,ExcitLevelIndex(2),ExcitLevelIndex(3)-1,PartIndex,tSuc)
                  dAmp=dAmp/Amplitude(1)
                  if(tSuc) then
                     Amplitude(PartIndex)=Amplitude(PartIndex)-dAmp*iSgn  !Remove this amount from the double excitation
                  else
                     Call WriteBitDet(6,iLutnI,.true.)
                     Call Stop_All("InitMP1Amplitude","Failed to find det.")
                  endif
               endif
            enddo
         endif
         Amplitude(j)=Amplitude(j)*ExcitToDetSign(iLutHF,ExcitList(:,j),IC)
      endif
      dTotAbsAmpl=dTotAbsAmpl+abs(Amplitude(j))
   enddo
END SUBROUTINE

!This runs over all amplitudes and initializsed them with small random values
! Amplitude(nExcit)            is the amplitude of each excitor
! nExcit                       is the number of excitors
SUBROUTINE InitRandAmplitude(Amplitude,nExcit,dInitAmp,dTotAbsAmpl)
   use CCMCData
   use SystemData, only: nEl
   use FciMCData, only: HFDet
   use FciMCParMod, only: iLutHF,SumEContrib,BinSearchParts3
   use Determinants, only: GetH0Element3
   use constants, only: dp
   use dSFMT_interface , only : genrand_real2_dSFMT
   IMPLICIT NONE
   REAL*8 Amplitude(nExcit)
   INTEGER nExcit
   REAL*8 dInitAmp,dTotAbsAmpl

   INTEGER j
   Amplitude(:)=0
   Amplitude(1)=dInitAmp
   dTotAbsAmpl=dInitAmp
   do j=2,nExcit
      Amplitude(j)=dInitAmp*(genrand_real2_dSFMT()-0.5)/1000
      dTotAbsAmpl=dTotAbsAmpl+abs(Amplitude(j))
   enddo
END SUBROUTINE

subroutine AttemptSpawn(S,C,Amplitude,dTol,TL,WalkerScale,iDebug)
   use SystemData, only: nEl
   use FciMCParMod, only: iLutHF
   use CCMCData, only: Spawner, Cluster,CCTransitionLog
   use FciMCParMod, only: BinSearchParts3
   Use CalcData, only: Tau
   use DetBitOps, only: FindBitExcitLevel
   use DetCalcData, only: FCIDets   ! (0:NIfDBO, Det).  Lists all allowed excitors in compressed form
   use DetCalcData, only:FCIDetIndex! (0:nEl+1).  The index of the different excitation levels
   use DetCalcData, only: Det       ! The number of Dets/Excitors in FCIDets
   Use Logging, only: tCCMCLogTransitions
   use FciMCData, only: Iter
   use CalcData, only: NEquilSteps
   implicit none
   type(Spawner) S
   type(Cluster) C
   real*8 Amplitude(:)
   real*8 dTol
   TYPE(CCTransitionLog) TL               ! Store data on transitions
   real*8 WalkerScale 
   integer iDebug

   real*8 rat
   integer i,j
   integer IC
   logical tSuc
   integer PartIndex
   IF(iDebug.gt.4) THEN
      WRITE(6,*) "  HIJ: ",S%HIJ
   ENDIF
   rat=-C%iSgn*Tau*S%HIJ*C%dAbsAmplitude/(S%dProbSpawn*C%dProbNorm*C%dClusterNorm)  

! C%dAbsAmplitude is there so that the change in the amp depends on the current amp.

   IF(iDebug.gt.3) THEN
!We've not printed this out before
      write(6,"(A)",advance="no") "    "
      call WriteBitEx(6,iLutHF,C%iLutDetCurr,.false.)
      write(6,'(A)',advance='no') ' ==> '
      call WriteBitEx(6,iLutHF,S%iLutnJ,.false.)
      WRITE(6,'(A,G25.16)',advance='no') "Children:",rat
      if(iDebug.eq.3.and.C%iSize.gt.1) THEN
         write(6,'(A)',advance='no') ' from '
         do i=1,C%iSize
            call WriteBitEx(6,iLutHF,C%SelectedExcitors(:,i),.false.)
         enddo
      endif
      write(6,*)
   endif
   if(abs(rat).gt.1e-4*dTol) then
!Now add in a contribution from the child
      IC = FindBitExcitLevel(iLutHF, S%iLutnJ(:), nEl)
      CALL BinSearchParts3(S%iLutnJ(:),FCIDets(:,:),Det,FCIDetIndex(IC),FCIDetIndex(IC+1)-1,PartIndex,tSuc)
      if(.not.tSuc) THEN      
         WRITE(6,*) "Cannot find excitor "
         call WriteBitEx(6,iLutHF,S%iLutnJ,.true.)
         call WriteBitDet(6,S%iLutnJ,.true.)
         WRITE(6,*) "Excitation Level: ",IC
         i=FCIDetIndex(IC)
         j=FCIDetIndex(IC+1)-1
         WRITE(6,*) "Dets ",i,' to ',j
         call WriteExcitorList(6,Amplitude(i:j),FciDets(:,i:j),i-1,j-i+1,0.d0,"Excitors in that level")
         call Stop_All("CCMCStandalone","Cannot find excitor in list.")
      endif
! We need to calculate the sign change from excitor to det:
!   Here we convert from a det back to an excitor.
      rat=rat*ExcitToDetSign(iLutHF,S%iLutnJ,IC)
      Amplitude(PartIndex)=Amplitude(PartIndex)+rat
      iter_data_ccmc%nborn=iter_data_ccmc%nborn+abs(rat*WalkerScale)
      if(tCCMCLogTransitions.and.Iter.gt.NEquilSteps) then
         call LogTransition(TL,C%SelectedExcitorIndices(:),C%iSize,PartIndex,rat,C%dProbNorm)
      endif
   endif
end subroutine AttemptSpawn


!Take cluster C and make an anti-excitor corresponding to its collapsed version to take into account its death.
subroutine AttemptDie(C,CurAmpl,OldAmpl,TL,WalkerScale,iDebug)
   use CCMCData, only: Cluster,CCTransitionLog
   use FciMCData, only: Hii
   Use CalcData, only: Tau,DiagSft
   use DetCalcData, only: FCIDets   ! (0:NIfDBO, Det).  Lists all allowed excitors in compressed form
   use DetCalcData, only:FCIDetIndex! (0:nEl+1).  The index of the different excitation levels
   use DetCalcData, only: Det       ! The number of Dets/Excitors in FCIDets
   use constants, only: dp
   use FciMCParMod, only: iLutHF
   Use Logging, only: lLogTransitions=>tCCMCLogTransitions
   use FciMCData, only: Iter
   use CalcData, only: NEquilSteps
   use FciMCParMod, only: BinSearchParts3
   use dSFMT_interface , only : genrand_real2_dSFMT
   
   implicit none
   Type(Cluster) C
   real*8 CurAmpl(:),OldAmpl(:)
   integer iDebug
   TYPE(CCTransitionLog) TL               ! Store data on transitions
   real*8 WalkerScale 

   INTEGER iC,iPartDie
   LOGICAL tSuc
   real*8 r,rat,HDiagCurr
   HElement_t Htmp
   integer i

! We have to decompose our composite excitor into one of its parts.  
   IF(C%iSize.GT.1) THEN
!This is an old version of death which may still work, but hasn't been tested
!!  We modify the composite (t_a t_b t_c) -> (t_a t_b t_c - x) by changing just one of the parts
!!  t_a -> t_a (1- x/(t_a t_b t_c)).   
! We need this in to keep the random number sequence.
! We need this in to keep the random number sequence.
      r = genrand_real2_dSFMT()  !On GHB's advice
!      k=1+floor(r*C%iSize)
!      iPartDie=C%SelectedExcitorIndices(k)
!We try an alternative death method - by creating an antiparticle in the excitor corresponding to this cluster
!
!      IC=C%iExcitLevel
!      CALL BinSearchParts3(C%iLutDetCurr(:),FCIDets(:,:),Det,FCIDetIndex(IC),FCIDetIndex(IC+1)-1,iPartDie,tSuc)
!To make things compatible from a cluster to the excitor,
! we divide by N_T ^|X| | t_x1 t_x2 ... | / |X|! and multiply by N_T t_X 
!NB N_T = 1/N_0

!We try an alternative death method - by creating an antiparticle in the excitor corresponding to this cluster

      IC=C%iExcitLevel
      CALL BinSearchParts3(C%iLutDetCurr(:),FCIDets(:,:),Det,FCIDetIndex(IC),FCIDetIndex(IC+1)-1,iPartDie,tSuc)
!To make things compatible from a cluster to the excitor,
! we divide by N_T ^|X| | t_x1 t_x2 ... | / |X|! and multiply by N_T t_X 
!NB N_T = 1/N_0
               
!dAbsAmplitude is | t_x1 t_x2 ... | / N_0 ^|X| 
   
!dAbsAmplitude is | t_x1 t_x2 ... | / N_0 ^|X| 
   ELSEIF(C%iSize.EQ.1) THEN
      iPartDie=C%SelectedExcitorIndices(1)
      IC=C%iExcitLevel
   ELSE
      iPartDie=1
      IC=0
   ENDIF 

   Htmp = get_helement (C%DetCurr, C%DetCurr, 0)
   HDiagCurr=Htmp
   HDiagCurr=HDiagCurr-Hii

   IF(iDebug.gt.4) then
      Write(6,'(A,I7)',advance='no') "Killing at excitor: ",iPartDie
      call WriteBitEx(6,iLutHF,FCIDets(:,iPartDie),.false.)
      Write(6,'(A,G25.16)') "Prob: ",C%dProbNorm
   endif

!  This will be the amount we wish to subtract from t_x

!dProb = 1
   rat=Tau*(HDiagCurr-DiagSft)/(C%dProbNorm*C%dClusterProb) !(dProb*dProbNorm)  !The old version
   rat=rat*(C%dAbsAmplitude*C%iSgn)

!   Here we convert from a det back to an excitor.
   rat=rat*ExcitToDetSign(iLutHF,FCIDets(:,iPartDie),IC)

! We've now calculated rat fully
   r=rat
   rat=rat/abs(OldAmpl(iPartDie))  !Take into account we're killing at a different place from the cluster

   IF(iDebug.ge.4) then
      WRITE(6,*) "   Product Contributions to Number Died:"
      WRITE(6,*) "    Energy difference: ",HDiagCurr-DiagSft
      WRITE(6,*) "    Tau              : ",Tau
      WRITE(6,*) "    Sign             : ",C%iSgn
      WRITE(6,*) "    Cluster Prob     : ",C%dClusterProb
      WRITE(6,*) "    Cluster Norm     : ",C%dClusterNorm
      WRITE(6,*) "    1/dProbNorm      : ",1/C%dProbNorm
      WRITE(6,*) "    dAbsAmplitude    : ",C%dAbsAmplitude
   endif 

!! rat is what we wish to modify t_a t_b t_c by. (but positive - we'll actually want to subtract it)
!!  To do this we modify the chosen part (e.g. t_a) by
!! t_a (1 - rat / (t_a t_b t_c) )
!
!! dAbsAmplitude = t_a t_b t_c
!            rat= rat/dAbsAmplitude
!! t_a(new) = t_a(new)+ t_a(old) * rat
   IF(iDebug.ge.4) then
      WRITE(6,*) "Death ratio      : ",rat
   endif

   CurAmpl(iPartDie)=CurAmpl(iPartDie)-r
   iter_data_ccmc%ndied=iter_data_ccmc%ndied+abs(r*WalkerScale)
   if(lLogTransitions.and.Iter.gt.NEquilSteps) then
      call LogTransition(TL,C%SelectedExcitorIndices(:),C%iSize,iPartDie,-r,C%dProbNorm)
   endif
   IF(iDebug.gt.3.) then
      Write(6,'(A,I7)',advance='no') " Killing at excitor: ",iPartDie
      Write(6,'(A)',advance='no') " chosen "
      call WriteBitEx(6,iLutHF,FCIDets(:,iPartDie),.false.)
      WRITE(6,'(A,G25.16)',advance='no') "Number died ",r
      if(C%iSize.gt.1) then
         WRITE(6,'(A)',advance='no') " from "
         do i=1,C%iSize
            call WriteBitEx(6,iLutHF,C%SelectedExcitors(:,i),.false.)
         enddo
      endif
      WRITE(6,*)
      if(iDebug.eq.4) Write(6,'(A,G25.16)') "Prob: ",C%dProbNorm
   endif
end subroutine AttemptDie

subroutine AttemptSpawnParticle(S,C,iDebug,SpawnList,SpawnAmps,nSpawned,nMaxSpawn)
   use SystemData, only: nEl
   use FciMCParMod, only: iLutHF
   use CCMCData, only: Spawner, Cluster
   Use CalcData, only: Tau
   use DetBitOps, only: FindBitExcitLevel
   USE dSFMT_interface , only : genrand_real2_dSFMT
   implicit none
   type(Spawner) S
   type(Cluster) C
   INTEGER(KIND=n_int) :: SpawnList(0:nIfTot,*)
   INTEGER SpawnAmps(:)
   integer nSpawned
   integer nMaxSpawn
   integer iSpawnAmp
   integer iDebug

   real*8 rat,r
   integer i
   integer IC
   IF(iDebug.gt.4) THEN
      WRITE(6,*) "  HIJ: ",S%HIJ
   ENDIF
   rat=-C%iSgn*Tau*S%HIJ*C%dAbsAmplitude/(S%dProbSpawn*C%dProbNorm*C%dClusterNorm)  

! C%dAbsAmplitude is there so that the change in the amp depends on the current amp.

   IF(iDebug.gt.3) THEN
!We've not printed this out before
      write(6,"(A)",advance="no") "    "
      call WriteBitEx(6,iLutHF,C%iLutDetCurr,.false.)
      write(6,'(A)',advance='no') ' ==> '
      call WriteBitEx(6,iLutHF,S%iLutnJ,.false.)
      WRITE(6,'(A,G25.16)',advance='no') " Children ratio:",rat
      if(iDebug.eq.3.and.C%iSize.gt.1) THEN
         write(6,'(A)',advance='no') ' from '
         do i=1,C%iSize
            call WriteBitEx(6,iLutHF,C%SelectedExcitors(:,i),.false.)
         enddo
      endif
      write(6,*)
   endif
!   Here we convert from a det back to an excitor.
   IC = FindBitExcitLevel(iLutHF, S%iLutnJ(:), nEl)
   rat=rat*ExcitToDetSign(iLutHF,S%iLutnJ,IC)
   r=abs(rat)
   iSpawnAmp=floor(r)
   if ((r-iSpawnAmp)>genrand_real2_dSFMT()) iSpawnAmp=iSpawnAmp+1
   if(iSpawnAmp>0) then

      nSpawned=nSpawned+1 !The index into the spawning list
      iter_data_ccmc%nborn=iter_data_ccmc%nborn+1
      if(nSpawned>nMaxSpawn) call Stop_All("AttemptSpawnParticle","Not enough space in spawning list.")
      if(rat>0) then
         SpawnAmps(nSpawned)=iSpawnAmp
      else
         SpawnAmps(nSpawned)=-iSpawnAmp
      endif
      SpawnList(:,nSpawned)=S%iLutnJ(:)
      IF(iDebug.gt.3) THEN
   !We've not printed this out before
         WRITE(6,*) "  Spawned ",iSpawnAmp
      endif
   else if(iDebug>3) then
      WRITE(6,*) "  No Spawning"
   endif
end subroutine AttemptSpawnParticle
!Take cluster C and make an anti-excitor corresponding to its collapsed version to take into account its death.
subroutine AttemptDieParticle(C,iDebug,SpawnList,SpawnAmps,nSpawned)
   use CCMCData, only: Cluster,CCTransitionLog
   use FciMCData, only: Hii
   Use CalcData, only: Tau,DiagSft
   use constants, only: dp
   use FciMCParMod, only: iLutHF
   use dSFMT_interface , only : genrand_real2_dSFMT
   use SystemData, only : nEl
   
   implicit none
   Type(Cluster) C
   INTEGER(KIND=n_int) :: SpawnList(0:nIfTot,*)
   INTEGER SpawnAmps(:)
   INTEGER nSpawned
   integer iDebug

   INTEGER iC,iPartDie
   real*8 r,rat,HDiagCurr
   HElement_t Htmp
   integer i
   integer iSpawnAmp

! We have to decompose our composite excitor into one of its parts.  
   IF(C%iSize.GT.1) THEN
!This is an old version of death which may still work, but hasn't been tested
!!  We modify the composite (t_a t_b t_c) -> (t_a t_b t_c - x) by changing just one of the parts
!!  t_a -> t_a (1- x/(t_a t_b t_c)).   
! We need this in to keep the random number sequence.
! We need this in to keep the random number sequence.
      r = genrand_real2_dSFMT()  !On GHB's advice
!      k=1+floor(r*C%iSize)
!      iPartDie=C%SelectedExcitorIndices(k)
!We try an alternative death method - by creating an antiparticle in the excitor corresponding to this cluster
!
!      IC=C%iExcitLevel
!      CALL BinSearchParts3(C%iLutDetCurr(:),FCIDets(:,:),Det,FCIDetIndex(IC),FCIDetIndex(IC+1)-1,iPartDie,tSuc)
!To make things compatible from a cluster to the excitor,
! we divide by N_T ^|X| | t_x1 t_x2 ... | / |X|! and multiply by N_T t_X 
!NB N_T = 1/N_0

!We try an alternative death method - by creating an antiparticle in the excitor corresponding to this cluster

      IC=C%iExcitLevel
!      CALL BinSearchParts3(C%iLutDetCurr(:),FCIDets(:,:),Det,FCIDetIndex(IC),FCIDetIndex(IC+1)-1,iPartDie,tSuc)
!To make things compatible from a cluster to the excitor,
! we divide by N_T ^|X| | t_x1 t_x2 ... | / |X|! and multiply by N_T t_X 
!NB N_T = 1/N_0
               
!dAbsAmplitude is | t_x1 t_x2 ... | / N_0 ^|X| 
   
!dAbsAmplitude is | t_x1 t_x2 ... | / N_0 ^|X| 
   ELSEIF(C%iSize.EQ.1) THEN
      iPartDie=C%SelectedExcitorIndices(1)
      IC=C%iExcitLevel
   ELSE
      iPartDie=1
      IC=0
   ENDIF 

   Htmp = get_helement (C%DetCurr, C%DetCurr, 0)
   HDiagCurr=Htmp
   HDiagCurr=HDiagCurr-Hii

!   IF(iDebug.gt.4) then
!      Write(6,'(A,I7)',advance='no') "Killing at excitor: ",iPartDie
!      call WriteBitEx(6,iLutHF,FCIDets(:,iPartDie),.false.)
!      Write(6,'(A,G25.16)') "Prob: ",C%dProbNorm
!   endif

!  This will be the amount we wish to subtract from t_x

!dProb = 1
   rat=Tau*(HDiagCurr-DiagSft)/(C%dProbNorm*C%dClusterProb) !(dProb*dProbNorm)  !The old version
   rat=rat*(C%dAbsAmplitude*C%iSgn)

!   Here we convert from a det back to an excitor.
   rat=rat*ExcitToDetSign(iLutHF,C%iLutDetCurr,IC)

! We've now calculated rat fully
   r=rat
!   rat=rat/abs(OldAmpl(iPartDie))  !Take into account we're killing at a different place from the cluster

   IF(iDebug.ge.4) then
      WRITE(6,*) "   Death ratio       : ",r
      WRITE(6,*) "   Product Contributions to Number Died:"
      WRITE(6,*) "    Energy difference: ",HDiagCurr-DiagSft
      WRITE(6,*) "    Tau              : ",Tau
      WRITE(6,*) "    Sign             : ",C%iSgn
      WRITE(6,*) "    Cluster Prob     : ",C%dClusterProb
      WRITE(6,*) "    Cluster Norm     : ",C%dClusterNorm
      WRITE(6,*) "    1/dProbNorm      : ",1/C%dProbNorm
      WRITE(6,*) "    dAbsAmplitude    : ",C%dAbsAmplitude
   endif 

!! rat is what we wish to modify t_a t_b t_c by. (but positive - we'll actually want to subtract it)
!!  To do this we modify the chosen part (e.g. t_a) by
!! t_a (1 - rat / (t_a t_b t_c) )
!
!! dAbsAmplitude = t_a t_b t_c
!            rat= rat/dAbsAmplitude
!! t_a(new) = t_a(new)+ t_a(old) * rat
!   IF(iDebug.ge.4) then
!      WRITE(6,*) "Death ratio      : ",rat
!   endif
   rat=-rat  !The we're creating death   
   r=abs(rat)
   iSpawnAmp=floor(r)
   if ((r-iSpawnAmp)>genrand_real2_dSFMT()) iSpawnAmp=iSpawnAmp+1
   if(iSpawnAmp>0) then
      nSpawned=nSpawned+1 !The index into the spawning list
      iter_data_ccmc%ndied=iter_data_ccmc%ndied+1
      IF(iDebug.gt.3.) then
         Write(6,'(A,I7)',advance='no') " Killing at excitor: ",iPartDie
         Write(6,'(A)',advance='no') " chosen "
         call WriteBitEx(6,iLutHF,SpawnList(:,nSpawned),.false.)
         WRITE(6,'(A,G25.16)',advance='no') "Number died ",r
         if(C%iSize.gt.1) then
            WRITE(6,'(A)',advance='no') " from "
            do i=1,C%iSize
               call WriteBitEx(6,iLutHF,C%SelectedExcitors(:,i),.false.)
            enddo
         endif
         WRITE(6,*)
         if(iDebug.eq.4) Write(6,'(A,G25.16)') "Prob: ",C%dProbNorm
      endif
      if(rat>0) then
         SpawnAmps(nSpawned)=iSpawnAmp
      else
         SpawnAmps(nSpawned)=-iSpawnAmp
      endif
      SpawnList(:,nSpawned)=C%iLutDetCurr(:)
   else if(iDebug>3) then
      write(6,*) "  No Death"
   endif
end subroutine AttemptDieParticle

SUBROUTINE CCMCStandalone(Weight,Energyxw)
   Use global_utilities
   use SystemData, only: nEl
   use Parallel, only: iProcIndex
   use FciMCData, only: root
   use CCMCData, only: tCCMCFCI,dInitAmplitude,dProbSelNewExcitor,tExactCluster,tExactSpawn,nSpawnings,tCCBuffer
   use CCMCData, only: ClustSelector,Spawner,CCTransitionLog,nClustSelections,tExactEnergy
   use DetCalcData, only: Det       ! The number of Dets/Excitors in FCIDets
   use DetCalcData, only: FCIDets   ! (0:NIfDBO, Det).  Lists all allowed excitors in compressed form
   use DetCalcData, only:FCIDetIndex! (0:nEl+1).  The index of the different excitation levels
   use CalcData, only: NMCyc    ! The number of MC Cycles
   use CalcData, only: StepsSft ! The number of steps between shift updates
   use CalcData, only: TStartMP1
   use FciMCData, only: Iter
   use FciMCData, only: TotParts,TotWalkers,TotWalkersOld,TotPartsOld,AllTotPartsOld,AllTotWalkersOld,AllTotParts
   use FciMCData, only: tTruncSpace
   use FciMCData, only: ProjectionE
   use FciMCParMod, only: iLutHF
   use FciMCParMod, only: CheckAllowedTruncSpawn, SetupParameters,BinSearchParts3
   use FciMCParMod, only: InitHistMin, calculate_new_shift_wrapper
   use FciMCData, only: NoatHF,NoatDoubs
   use FciMCParMod, only: WriteHistogram,SumEContrib
   Use Logging, only: CCMCDebug,tCCMCLogTransitions,tCCMCLogUniq
   USE Logging , only : tHistSpawn,iWriteHistEvery
   USE DetCalcData , only : ICILevel
   use CalcData, only: NEquilSteps
   use FciMCParMod, only: WriteFciMCStats, WriteFciMCStatsHeader
   use constants, only: dp
   use CCMCData, only: WriteCluster
   use ClusterList
   use CalcData, only: DiagSft
   IMPLICIT NONE
   real(dp) Weight,EnergyxW
   TYPE(AmplitudeList_double),target :: AL

   INTEGER iNumExcitors          ! The number of non-zero excitors (excluding the ref det)
   REAL*8 dTotAbsAmpl            ! The total of the absolute amplitudes

   INTEGER iCurAmpList,iOldAmpList  !the index of current and previous amplitude lists
   INTEGER iDebug
   INTEGER i
! Temporary Storage
   INTEGER PartIndex,IC          ! Used in buffering
   LOGICAL tSuc                  ! Also used in buffering

   INTEGER iOldTotWalkers        ! Info user for update to calculate shift
   INTEGER iShiftLeft            ! Number of steps left until we recalculate shift
   REAL*8 WalkerScale            ! Scale factor for turning floating point amplitudes into integer walkers.
   REAL*8 dProjE                 ! Stores the Projected Energy
   REAL*8 dTolerance             ! The tolerance for when to regard a value as zero
   REAL*8 dAveTotAbsAmp          ! Average of Total absolute amplitude over all post-equil cycles
   REAL*8 dAveNorm               ! Average of Normalization (ampl of Ref) over all post-equil cycles
   REAL*8 dAmpPrintTol           ! What size amplitudes do we bother printing 
   LOGICAL lLogTransitions       ! Do we log transitions

   TYPE(ClustSelector),target :: CSMain   ! A normal ClustSelector based on the current amplitudes
   TYPE(ClustSelector),target :: CSBuff   ! This is used when we're doing buffered CC
   TYPE(ClustSelector),pointer :: CS      ! This will point to the appropriate selector
   LOGICAL tPostBuffering                 ! Set after prebuffering
   TYPE(AmplitudeList_double), target :: ALBuffer !(Det)  used for buffered CC, storing intermediate amplitudes from cluster generation
   LOGICAL tMoreClusters                  ! Indicates we've not finished selecting clusters 
   TYPE(AmplitudeList_double), pointer :: OldAL          ! The previous cycle's amplitudes
   INTEGER oldALIndex                     !The index in OldAL of the list to use

   TYPE(Spawner) S                        ! A spawner used to generate spanees from a cluster

   INTEGER nAmpl                          ! Number of excitors on the Amplitude array
   INTEGER nBuffAmpl                      ! Number of determinants in the buffered Amplitude array
   INTEGER nCurAmpl                       ! Set to nAmpl or nBuffAmpl depending where OldAL points
   INTEGER iExcitLevelCluster             ! The maximum excitation level a cluster can be at
   INTEGER iMaxAmpLevel                   ! The maximum excitation level of a stored amplitude

   TYPE(CCTransitionLog) TL               ! Store data on transitions
   INTEGER iRefPos
   INTEGER, DIMENSION(lenof_sign) :: TempSign   !ghb24: For compatibility with new walker arrays & routines
   TYPE(timer) :: CCMC_time



   WRITE(6,*) "Entering CCMC Standalone..."
   CCMC_time%timer_name='CCMC Standalone'
   call set_timer(CCMC_time,20)

   iRefPos=1  !Always first element
   iDebug=CCMCDebug

   Call SetupParameters()

   ! Reset counters
   iter_data_ccmc%nborn = 0
   iter_data_ccmc%ndied = 0
   iter_data_ccmc%nannihil = 0

   lLogTransitions=tCCMCLogTransitions
   dTolerance=0 !1e-16
   iCurAmpList=1  !Start with list 1

   if(ICILevel/=0)  then 
      nAmpl=FciDetIndex(ICILevel+1)-1
      iMaxAmpLevel=ICILevel
   else
      nAmpl=Det
      iMaxAmpLevel=nEl
   endif
   write(6,*) "Number of stored amplitudes: ",nAmpl
! Setup Memory
   call AllocateAmplitudeList(AL,nAmpl,2)

   if(tTruncSpace) then
      if(tCCMCFCI) then
         iExcitLevelCluster=min(ICILevel,nEl)
      else
         iExcitLevelCluster=min(ICILevel+2,nEl)
      endif
      write(6,*) "Cluster excitation level cutoff: ",iExcitLevelCluster
   else
      iExcitLevelCluster=nEl
      write(6,*) "No cluster excitation level cutoff"
   endif
   if(tCCBuffer) then !CCBuffer uses a third array to store the sum of collapsed clusters' amplitudes, and spawns from that
      if(ICILevel/=0)  then 
         nBuffAmpl=FciDetIndex(iExcitLevelCluster+1)-1
      else
         nBuffAmpl=Det
      endif
      WRITE(6,*) "Buffered Amplitudes:",nBuffAmpl
      call AllocateAmplitudeList(ALBuffer,nBuffAmpl,1)
   endif

! Now setup the amplitude list.  Let's start with nothing initially, and
   AL%Amplitude(:,:)=0
!  place ampl 1 in the HF det
   if(tStartMP1) then
      write(6,*) "Initializing with MP1 amplitudes."
      CALL InitMP1Amplitude(tCCMCFCI,AL%Amplitude(:,iCurAmpList),nAmpl,FciDets,FCIDetIndex,dInitAmplitude,dTotAbsAmpl)
   elseif(ICILevel==1) then
      write(6,*) "Initializing with random amplitudes."
      CALL InitRandAmplitude(AL%Amplitude(:,iCurAmpList),nAmpl,dInitAmplitude,dTotAbsAmpl)
   else
      AL%Amplitude(1,iCurAmpList)=dInitAmplitude
      iNumExcitors=0
      dTotAbsAmpl=AL%Amplitude(1,iCurAmpList)
   endif
   dAmpPrintTol=(dTolerance*dInitAmplitude)
   if(iDebug.ge.4) dAmpPrintTol=0

   iShiftLeft=StepsSft-1  !So the first one comes at StepsSft
   if(iProcIndex.eq.root) then
      WalkerScale=100000/dInitAmplitude
   else
      WalkerScale=0
   endif
   TotWalkers=0
!WalkerScale*dTotAbsAmpl
   TotParts(1)=0
!WalkerScale*dTotAbsAmpl
   TotWalkersOld=0
!WalkerScale*dTotAbsAmpl
   TotPartsOld(1)=0
!WalkerScale*dTotAbsAmpl
   AllTotWalkersOld=1
   AllTotParts(1)=WalkerScale*dTotAbsAmpl
   AllTotPartsOld(1)=WalkerScale*dTotAbsAmpl
   iOldTotWalkers=WalkerScale*dTotAbsAmpl
   iter_data_ccmc%tot_parts_old = WalkerScale * dTotAbsAmpl
   dAveTotAbsAmp=0
   dAveNorm=0
   Iter=1
   IF(tHistSpawn) THEN
      Call InitHistMin() !Setup Histogramming arrays if needed 
   ENDIF
   CALL WriteFciMCStatsHeader()

   if(tCCMCFCI) THEN
      iNumExcitors=1  !Only one excitor allowed
   ELSE
      IF (tTruncSpace) THEN !we can go up to two beyond the max level of truncation
         iNumExcitors=ICILevel+2  !We need to be able to couple (say) 4 singles to make a quad and then spawn back to the sdoubles space
      ELSE
         iNumExcitors=nEl  !Otherwise just the number of electrons is the max number of excitors
      ENDIF
   ENDIF

   if(lLogTransitions) call InitTransitionLog(TL,nAmpl,iNumExcitors,.not.tCCMCLogUniq)

   if(tExactCluster) then
      CALL InitClustSelectorFull(CSMain,iNumExcitors)
   else
      CALL InitClustSelectorRandom(CSMain,iNumExcitors,nClustSelections,dProbSelNewExcitor)
   endif
   if(tCCBuffer) then
      CALL InitClustSelectorFull(CSBuff,1)
   endif

   CALL InitSpawner(S,tExactSpawn,ICILevel)

! Each cycle we select combinations of excitors randomly, and spawn and birth/die from them
   do while (Iter.le.NMCyc)
      if(iDebug>3)  write(79,*) "Cycle", Iter
      if(iDebug>3)  write(89,*) "Cycle", Iter
! Copy the old Amp list to the new
      iOldAmpList=iCurAmpList
      iCurAmpList=3-iCurAmpList
      AL%Amplitude(:,iCurAmpList)=AL%Amplitude(:,iOldAmpList)
      IF(iDebug.gt.1) THEN
         write(6,*) "Cycle ",Iter
         call WriteExcitorList(6,AL%Amplitude(:,iCurAmpList),FciDets,0,nAmpl,dAmpPrintTol,"Excitor list")
      endif
      call CalcTotals(iNumExcitors,dTotAbsAmpl,AL%Amplitude(:,iCurAmpList),nAmpl,dTolerance*dInitAmplitude,WalkerScale,iRefPos,iOldTotWalkers,iDebug)
      if(tExactEnergy) then
         CALL CalcClusterEnergy(tCCMCFCI,AL%Amplitude(:,iCurAmpList),nAmpl,FciDets,FCIDetIndex,iRefPos,iDebug,dProjE)
      else
         dProjE=ProjectionE
      endif
! Collate stats


      IF(iDebug.gt.1) THEN
         WRITE(6,*) "Total non-zero excitors: ",iNumExcitors
         WRITE(6,"(A,G30.22)") "Total abs Amplitudes: ",dTotAbsAmpl
         WRITE(6,*) "Projected Energy: ",dProjE
      endif 



!  Loop over cluster selections
!  Point to the main cluster selector, not the buffer
      CS=>CSMain
      call ResetClustSelector(CS,iRefPos)
      if(tCCBuffer) then
         call ResetClustSelector(CSBuff,iRefPos)
         ALBuffer%Amplitude(:,1)=0
      endif
      tMoreClusters=.true.
      tPostBuffering=.false.
      nCurAmpl=nAmpl
      OldAL=>AL
      OldALIndex=iOldAmpList
!A standard run has PostBuffering .false. and selects purely from the main cluster selector.
!A buffered run first has PostBuffering .false., selecting from the main CS, and storing the cumulative amplitudes in ALBuffer
!        second, it has   PostBuffering .true.,  selecting from the CSBuff and spawning and dying from there.
      call AccumulateAmplitudeList(OldAL,nCurAmpl,OldALIndex,iRefPos)

      do while (tMoreClusters)
         if(.not.tPostBuffering) then
            i=min(iNumExcitors,nEl)
            tMoreClusters=GetNextCluster(CS,FciDets,nAmpl,AL,iOldAmpList,dTotAbsAmpl,i,iDebug)
            if(tCCBuffer.and..not.tMoreClusters) then
               if(iDebug.gt.2) WRITE(6,*) "Buffering Complete.  Now Spawning."
               tPostBuffering=.true.
               nCurAmpl=nBuffAmpl
               CS=>CSBuff
               OldAL=>ALBuffer
               OldALIndex=1
               if(iDebug.gt.2) call WriteExcitorList(6,ALBuffer%Amplitude(:,OldALIndex),FciDets,0,nBuffAmpl,dAmpPrintTol,"Cluster expanded wavefunction")
               call AccumulateAmplitudeList(OldAL,nCurAmpl,OldALIndex,iRefPos)
            endif
         endif
         if(tPostBuffering) then  !If we've finished buffering, we read from the buffer (which is now pointed to by OldAL
            i=min(iNumExcitors,nEl)
            tMoreClusters=GetNextCluster(CSBuff,FciDets,nBuffAmpl,OldAL,OldALIndex,dTotAbsAmpl,i,iDebug)
         endif
         if(.not.tMoreClusters) exit
!Now form the cluster
         IF(iDebug.gt.3) then
            write(6,*) "Selection ",CS%iIndex
            WRITE(6,*) " Excitors in composite:", CS%C%iSize
            do i=1,CS%C%iSize
               call WriteBitEx(6,iLutHF,CS%C%SelectedExcitors(:,i),.true.)
            enddo
            Write(6,"(A,G25.17)") "   Select Prob given level: ",CS%C%dClusterProb
            Write(6,"(A,G25.17)") "   Prob norm              : ",CS%C%dProbNorm
            Write(6,"(A,G25.17)") "   Cluster norm           : ",CS%C%dClusterNorm
         endif

!The final logic tells it whether to convert from an excitor to a det.
         CALL CollapseCluster(CS%C,iLutHF,OldAL%Amplitude(:,OldALIndex),nCurAmpl,iDebug,.not.(tCCBuffer.and.tPostBuffering))
         IF(CS%C%iSgn/=0.and.iDebug.gt.4) then
            WRITE(6,*) "Chosen det/excitor is:"
            WRITE(6,"(A)",advance="no") "  "
            call WriteBitDet(6,CS%C%iLutDetCurr,.true.)
            CALL FLUSH(6)
         endif
         if(tTruncSpace.and.CS%C%iExcitLevel>iExcitLevelCluster) cycle !Don't try to die if we're truncated

         if(lLogTransitions.and.Iter.gt.NEquilSteps) call LogCluster(TL,CS%C)
         if(iDebug.gt.5) call WriteCluster(6,CS%C,.true.)
         if(CS%C%iSgn.eq.0) then
            if(iDebug.gt.3) write(6,*) "Sign Zero so moving to next cluster"
            cycle  !No point in doing anything
         endif

         !If we're in buffer mode, save in the buffer and carry on.
         if(tCCBuffer.and..not.tPostBuffering) then
            IC=CS%C%iExcitLevel
            CALL BinSearchParts3(CS%C%iLutDetCurr(:),FCIDets(:,:),Det,FCIDetIndex(IC),FCIDetIndex(IC+1)-1,PartIndex,tSuc)
            if(.not.tSuc) then
               call write_det (6, CS%C%DetCurr, .true.)
               Call Stop_All("CCMCStandalone","Failed to find det.")
            endif
            ALBuffer%Amplitude(PartIndex,1)=ALBuffer%Amplitude(PartIndex,1)+CS%C%dAbsAmplitude*CS%C%iSgn
            cycle
         endif

         if(iDebug.gt.3) WRITE(6,*) "Cluster Amplitude: ",CS%C%iSgn*CS%C%dAbsAmplitude 
!         if(iDebug.gt.3) WRITE(6,*) " Cluster Prob: ",CS%C%dSelectionProb
         if(.not.tExactEnergy.and.CS%C%iExcitLevel.le.2) then
            if(iDebug>3) then
               call WriteCluster(79,CS%C,.false.)
               call WriteCluster(89,CS%C,.false.)
               write(79,*)  (CS%C%dAbsAmplitude/AL%Amplitude(1,iOldAmpList))/CS%C%dSelectionProb,CS%C%iSgn !,CS%C%dSelectionProb,CS%C%dAbsAmplitude
               write(89,*)  CS%C%iSgn ,CS%C%dSelectionProb,CS%C%dProbNorm,CS%C%dAbsAmplitude,CS%C%dSelectionNorm
            endif
            TempSign(1)=CS%C%iSgn
            CALL SumEContrib(CS%C%DetCurr,CS%C%iExcitLevel,TempSign,CS%C%iLutDetCurr,0.d0,1/CS%C%dSelectionNorm)
         endif
!Now consider a number of possible spawning events
         CALL ResetSpawner(S,CS%C,nSpawnings)

!GetNextSpawner will generate either all possible spawners sequentially, or a single randomly chosen one (or none at all, if the randomly chosen one is disallowed)
         do while (GetNextSpawner(S,iDebug))
            if(.not.S%bValid) cycle
            call AttemptSpawn(S,CS%C,AL%Amplitude(:,iCurAmpList),dTolerance*dInitAmplitude,TL,WalkerScale,iDebug)
         enddo !GetNextSpawner
! Now deal with birth/death.
         if((.not.tTruncSpace).or.CS%C%iExcitLevel<=iMaxAmpLevel)          &
  &            call AttemptDie(CS%C,AL%Amplitude(:,iCurAmpList),AL%Amplitude(:,iOldAmpList),TL,WalkerScale,iDebug)
      enddo ! Cluster choices


      IF(tHistSpawn.and.(mod(Iter,iWriteHistEvery).eq.0)) THEN
          CALL WriteHistogram()
      ENDIF
      if(Iter.gt.NEquilSteps) then
         dAveTotAbsAmp=dAveTotAbsAmp+dTotAbsAmpl
         dAveNorm=dAveNorm+AL%Amplitude(1,iCurAmpList)
      endif
! Calc Shift
      iShiftLeft=iShiftLeft-1

!TotWalkers is used for this and is WalkerScale* total of all amplitudes
      NoAtHF=AL%Amplitude(iRefPos,iCurAmpList)
      if(iShiftLeft.le.0)  Call calculate_new_shift_wrapper(iter_data_ccmc, &
                                                            TotParts)
      if(iShiftLeft.le.0)  iShiftLeft=StepsSft
      Iter=Iter+1
!Reset number at HF and doubles
      NoatHF=0
      NoatDoubs=0
   enddo !MC Cycles

      if(iDebug.gt.1) call WriteExcitorList(6,AL%Amplitude(:,iCurAmpList),FciDets,0,nAmpl,dAmpPrintTol,"Final Excitor list")

! Find the largest 10 amplitudes in each level
      call WriteMaxExcitorList(6,AL%Amplitude(:,iCurAmpList),FciDets,FCIDetIndex,iMaxAmpLevel,10)
      nullify(OldAL)
      nullify(CS)
      if(lLogTransitions) then
         call WriteTransitionLog(6,TL,AL%Amplitude(:,iCurAmpList),NMCyc-NEquilSteps,dAveTotAbsAmp,dAveNorm)
         call DeleteTransitionLog(TL)
   endif
   call DeallocateAmplitudeList(AL)
   if(tCCBuffer) then
      call DeAllocateAmplitudeList(ALBuffer)
   endif
   Weight=0.D0
   Energyxw=ProjectionE
   call halt_timer(CCMC_time)
END SUBROUTINE CCMCStandalone

SUBROUTINE CCMCStandaloneParticle(Weight,Energyxw)
   Use global_utilities
   use SystemData, only: nEl
   use Parallel, only: iProcIndex
   use FciMCData, only: root
   use CCMCData, only: tCCMCFCI,dInitAmplitude,dProbSelNewExcitor,tExactCluster,tExactSpawn,nSpawnings,tCCBuffer
   use CCMCData, only: WriteCluster
   use CCMCData, only: ClustSelector,Spawner,nClustSelections
   use CalcData, only: NMCyc    ! The number of MC Cycles
   use CalcData, only: StepsSft ! The number of steps between shift updates
   use CalcData, only: TStartMP1
   use FciMCData, only: Iter
   use FciMCData, only: TotParts,TotWalkers,TotWalkersOld,TotPartsOld,AllTotPartsOld,AllTotWalkersOld,AllTotParts
   use FciMCData, only: NoatHF,NoatDoubs
   use FciMCData, only: tTruncSpace
   use FciMCData, only: ProjectionE
   use FciMCParMod, only: iLutHF
   use FciMCParMod, only: CheckAllowedTruncSpawn, SetupParameters,BinSearchParts3
   use FciMCParMod, only: InitHistMin, calculate_new_shift_wrapper
   use FciMCParMod, only: WriteHistogram,SumEContrib
   Use Logging, only: CCMCDebug,tCCMCLogTransitions,tCCMCLogUniq
   USE Logging , only : tHistSpawn,iWriteHistEvery
   USE DetCalcData , only : ICILevel
   use CalcData, only: InitWalkers,NEquilSteps
   use FciMCParMod, only: WriteFciMCStats, WriteFciMCStatsHeader
   USE dSFMT_interface , only : genrand_real2_dSFMT
   use CalcData, only: MemoryFacSpawn
   use AnnihilationMod, only: AnnihilationInterface
   use CalcData, only: DiagSft
   use CalcData, only: TStartSinglePart
   use timing, only: print_timing_report
   IMPLICIT NONE
   real(dp) Weight,EnergyxW
   CHARACTER(len=*), PARAMETER :: this_routine='CCMCStandaloneParticle'
   TYPE(AmplitudeList_int),target :: AL
   INTEGER(kind=n_int), allocatable :: DetList(:,:)
   INTEGER  tagDetList

   INTEGER iNumExcitors          ! The number of non-zero excitors (excluding the ref det)
   REAL*8 dTotAbsAmpl            ! The total of the absolute amplitudes

   INTEGER iDebug
   INTEGER i,iMin
! Temporary Storage

   INTEGER iOldTotWalkers        ! Info user for update to calculate shift
   INTEGER iShiftLeft            ! Number of steps left until we recalculate shift
   REAL*8 WalkerScale            ! Scale factor for turning floating point amplitudes into integer walkers.
   REAL*8 dTolerance             ! The tolerance for when to regard a value as zero
   REAL*8 dAveTotAbsAmp          ! Average of Total absolute amplitude over all post-equil cycles
   REAL*8 dAveNorm               ! Average of Normalization (ampl of Ref) over all post-equil cycles
   INTEGER dAmpPrintTol           ! What size amplitudes do we bother printing 

   TYPE(ClustSelector) :: CS   ! A normal ClustSelector based on the current amplitudes
   LOGICAL tMoreClusters                  ! Indicates we've not finished selecting clusters 

   TYPE(Spawner) S                        ! A spawner used to generate spanees from a cluster

   INTEGER nAmpl                          ! Number of excitors on the Amplitude array
   INTEGER iExcitLevelCluster             ! The maximum excitation level a cluster can be at
   INTEGER iMaxAmpLevel                   ! The maximum excitation level of a stored amplitude

   INTEGER, parameter :: iCurAmpList=1     !Just for futureproofness at the moment - there is only one
   
   INTEGER(kind=n_int), allocatable :: SpawnList(:,:)
   INTEGER, allocatable ::  SpawnAmps(:)
   INTEGER tagSpawnList,tagSpawnAmps
   INTEGER nSpawned,nMaxSpawn

   LOGICAL tS

   LOGICAL tShifting                      ! Are we in variable shift mode

   integer nMaxAmpl

   integer iRefPos      !The location of the reference det in the Amplitude array

   !ghb24: Changes to allow compatibility with the new packaged walkers.
   INTEGER, DIMENSION(lenof_sign) :: TempSign
   TYPE(timer) :: CCMC_time


   WRITE(6,*) "Entering CCMC Standalone Particle..."
   CCMC_time%timer_name='CCMC Standalone Particle'
   call set_timer(CCMC_time,20)

!   Spawntime%timer_name='SpawnTime'
!   Dietime%timer_name='DieTime'
!   Etime%timer_name='ETime'
   iRefPos=1
   iDebug=CCMCDebug

   Call SetupParameters()

   ! Reset counters
   iter_data_ccmc%nborn = 0
   iter_data_ccmc%ndied = 0
   iter_data_ccmc%nannihil = 0

   dTolerance=0 !1e-16

   if(ICILevel/=0)  then 
      iMaxAmpLevel=ICILevel
   else
      iMaxAmpLevel=nEl
   endif
   nMaxAmpl=InitWalkers
   WalkerScale=1
! Setup Memory
   write(6,*) "Max Amplitude List size: ", nMaxAmpl
   call AllocateAmplitudeList(AL,nMaxAmpl,1)
   Allocate(DetList(0:nIfTot,nMaxAmpl))
   LogAlloc(ierr,'DetList',(nIfTot+1)*nMaxAmpl,4,tagDetList)
   nMaxSpawn=MemoryFacSpawn*nMaxAmpl
   Allocate(SpawnList(0:nIfTot,nMaxSpawn))
   LogAlloc(ierr,'SpawnList',(nIfTot+1)*nMaxAmpl,4,tagSpawnList)
   Allocate(SpawnAmps(nMaxSpawn))
   LogAlloc(ierr,'SpawnAmps',nMaxAmpl,4,tagSpawnAmps)

   if(tTruncSpace) then
      if(tCCMCFCI) then
         iExcitLevelCluster=min(ICILevel,nEl)
      else
         iExcitLevelCluster=min(ICILevel+2,nEl)
      endif
      write(6,*) "Cluster excitation level cutoff: ",iExcitLevelCluster
   else
      iExcitLevelCluster=nEl
      write(6,*) "No cluster excitation level cutoff"
   endif

! Now setup the amplitude list.  Let's start with nothing initially, and
   AL%Amplitude(:,:)=0
!  !  place ampl 1 in the HF det
!   if(tStartMP1) then
!      write(6,*) "Initializing with MP1 amplitudes."
!      CALL InitMP1Amplitude(tCCMCFCI,Amplitude(:,iCurAmpList),nAmpl,FciDets,FCIDetIndex,dInitAmplitude,dTotAbsAmpl)
!   else
   if(TStartSinglePart) then
      AL%Amplitude(iRefPos,iCurAmpList)=dInitAmplitude
      tShifting=.false.
   else
      AL%Amplitude(iRefPos,iCurAmpList)=dInitAmplitude
      tShifting=.true.
   endif
   DetList(:,1)=iLutHF 
      nAmpl=1
      iNumExcitors=0
   dTotAbsAmpl=AL%Amplitude(iRefPos,iCurAmpList)
!   endif
   dAmpPrintTol=(dTolerance*dInitAmplitude)
   if(iDebug.ge.4) dAmpPrintTol=0

   iShiftLeft=StepsSft-1  !So the first one comes at StepsSft
   TotWalkers=0
!WalkerScale*dTotAbsAmpl
   TotParts(1)=0
!WalkerScale*dTotAbsAmpl
   TotWalkersOld=0
!WalkerScale*dTotAbsAmpl
   TotPartsOld(1)=0
!WalkerScale*dTotAbsAmpl
   AllTotWalkersOld=1
!WalkerScale*dTotAbsAmpl
   AllTotParts(1)=WalkerScale*dTotAbsAmpl
   AllTotPartsOld(1)=WalkerScale*dTotAbsAmpl
   iOldTotWalkers=WalkerScale*dTotAbsAmpl
   iter_data_ccmc%tot_parts_old = WalkerScale * dTotAbsAmpl
   dAveTotAbsAmp=0
   dAveNorm=0
   Iter=1
   IF(tHistSpawn) THEN
      Call InitHistMin() !Setup Histogramming arrays if needed 
   ENDIF
   CALL WriteFciMCStatsHeader()

   if(tCCMCFCI) THEN
      iNumExcitors=1  !Only one excitor allowed
   ELSE
      IF (tTruncSpace) THEN !we can go up to two beyond the max level of truncation
         iNumExcitors=ICILevel+2  !We need to be able to couple (say) 4 singles to make a quad and then spawn back to the sdoubles space
      ELSE
         iNumExcitors=nEl  !Otherwise just the number of electrons is the max number of excitors
      ENDIF
   ENDIF

   CALL InitClustSelectorRandom(CS,iNumExcitors,nClustSelections,dProbSelNewExcitor)

   CALL InitSpawner(S,tExactSpawn,ICILevel)

! Each cycle we select combinations of excitors randomly, and spawn and birth/die from them
   do while (Iter.le.NMCyc)
!Find teh HF det
      CALL BinSearchParts3(iLutHF,DetList,nAmpl,1,nAmpl,iRefPos,tS)
      if(.not.tS) call Stop_All("CCMCStandaloneParticle","Failed to find HF det.")
! Collate stats
      nSpawned=0
      IF(iDebug.gt.1) THEN
         write(6,*) "Cycle ",Iter
         call WriteExcitorList(6,AL%Amplitude(:,iCurAmpList),DetList,0,nAmpl,dAmpPrintTol,"Excitor list")
      endif


      call CalcTotals(iNumExcitors,dTotAbsAmpl,AL%Amplitude(:,iCurAmpList),nAmpl,dTolerance*dInitAmplitude,WalkerScale,iRefPos,iOldTotWalkers,iDebug)
!      if(.not.tShifting) then
!         if(iNumExcitors>dInitAmplitude) then
!            tShifting=.true.
!         else
!            iShiftLeft=2
!         endif
!      endif

      
      IF(iDebug.gt.1) THEN
         WRITE(6,*) "Total non-zero excitors: ",iNumExcitors
         WRITE(6,"(A,G30.22)") "Total abs Amplitudes: ",dTotAbsAmpl
      endif 


! Calc Shift
      iShiftLeft=iShiftLeft-1

      NoAtHF=AL%Amplitude(iRefPos,iCurAmpList)
!TotWalkers is used for this and is WalkerScale* total of all amplitudes
      if(iShiftLeft.le.0)  Call calculate_new_shift_wrapper(iter_data_ccmc, &
                                                            TotParts)
      if(iShiftLeft.le.0)  iShiftLeft=StepsSft
!Reset number at HF and doubles
      NoatHF=0
      NoatDoubs=0

!  Loop over cluster selections
!  Point to the main cluster selector, not the buffer
      call ResetClustSelector(CS,iRefPos)
      tMoreClusters=.true.
      iMin=min(iNumExcitors,nEl)
      call AccumulateAmplitudeList(AL,nAmpl,iCurAmpList,iRefPos)
      do while (tMoreClusters)
         tMoreClusters=GetNextCluster(CS,DetList,nAmpl,AL,iCurAmpList,dTotAbsAmpl,iMin,iDebug)
         if(.not.tMoreClusters) exit
!Now form the cluster
         IF(iDebug.gt.3) then
            write(6,*) "Selection ",CS%iIndex
            WRITE(6,*) " Excitors in composite:", CS%C%iSize
            do i=1,CS%C%iSize
               call WriteBitEx(6,iLutHF,CS%C%SelectedExcitors(:,i),.true.)
            enddo
            Write(6,"(A,G25.17)") "   Select Prob given level: ",CS%C%dClusterProb
            Write(6,"(A,G25.17)") "   Prob norm              : ",CS%C%dProbNorm
            Write(6,"(A,G25.17)") "   Cluster norm           : ",CS%C%dClusterNorm
         endif

!The final logic tells it whether to convert from an excitor to a det.
         CALL CollapseCluster(CS%C,iLutHF,AL%Amplitude(:,iCurAmpList),nAmpl,iDebug,.true.)
         IF(CS%C%iSgn/=0.and.iDebug.gt.4) then
            WRITE(6,*) " Chosen det/excitor is:"
            WRITE(6,"(A)",advance="no") "  "
            call WriteBitDet(6,CS%C%iLutDetCurr,.true.)
            CALL FLUSH(6)
         endif
         if(tTruncSpace.and.CS%C%iExcitLevel>iExcitLevelCluster) cycle !Don't try to die if we're truncated

         if(iDebug.gt.5) then
            write(6, "(A)", advance='no') "  "
            call WriteCluster(6,CS%C,.true.)
         endif
         if(CS%C%iSgn.eq.0) then
            if(iDebug.gt.3) write(6,*) "Sign Zero so moving to next cluster"
            cycle  !No point in doing anything
         endif

         if(iDebug.gt.3) WRITE(6,*) " Cluster Amplitude: ",CS%C%iSgn*CS%C%dAbsAmplitude 
!         if(iDebug.gt.3) WRITE(6,*) " Cluster Prob: ",CS%C%dSelectionProb
         TempSign(1)=CS%C%iSgn
!         call set_timer(Etime,20)
         CALL SumEContrib(CS%C%DetCurr,CS%C%iExcitLevel,TempSign,CS%C%iLutDetCurr,0.d0,1/CS%C%dSelectionNorm)
!         call halt_timer(Etime)
!         call set_timer(Spawntime,20)
!Now consider a number of possible spawning events
         CALL ResetSpawner(S,CS%C,nSpawnings)

!GetNextSpawner will generate either all possible spawners sequentially, or a single randomly chosen one (or none at all, if the randomly chosen one is disallowed)
         do while (GetNextSpawner(S,iDebug))
            if(.not.S%bValid) cycle
            call AttemptSpawnParticle(S,CS%C,iDebug,SpawnList(:,:),SpawnAmps(:),nSpawned,nMaxSpawn)
         enddo !GetNextSpawner
! Now deal with birth/death.
!         call halt_timer(Spawntime)
!         call set_timer(Dietime,20)
         if((.not.tTruncSpace).or.CS%C%iExcitLevel<=iMaxAmpLevel)          &
  &         call AttemptDieParticle(CS%C,iDebug,SpawnList,SpawnAmps,nSpawned)
!         call halt_timer(Dietime)
      enddo ! Cluster choices

! At this point SpawnList contains a set of newly spawned particles and SpawnAmps the amount spawned
      if(nSpawned>0) then
         if(iDebug>2) write(6,*) "Calling Annihilation with ", nSpawned, " spawned."
         if(iDebug>2) call WriteExcitorList(6,SpawnAmps,SpawnList,0,nSpawned,dAmpPrintTol,"Spawned list")
         call AnnihilationInterface(nAmpl,DetList,AL%Amplitude(:,iCurAmpList),nMaxAmpl,nSpawned,SpawnList,SpawnAmps,nMaxSpawn,iter_data_ccmc)
      else
         if(iDebug>2) write(6,*) "No spawnings in toto."
      endif 
      


      IF(tHistSpawn.and.(mod(Iter,iWriteHistEvery).eq.0)) THEN
          CALL WriteHistogram()
      ENDIF
      if(Iter.gt.NEquilSteps) then
         dAveTotAbsAmp=dAveTotAbsAmp+dTotAbsAmpl
         dAveNorm=dAveNorm+AL%Amplitude(iRefPos,iCurAmpList)
      endif
      Iter=Iter+1
   enddo !MC Cycles

   if(iDebug.gt.1) call WriteExcitorList(6,AL%Amplitude(:,iCurAmpList),DetList,0,nAmpl,dAmpPrintTol,"Final Excitor list")

! Find the largest 10 amplitudes in each level
!   call WriteMaxExcitorList(6,AL%Amplitude(:,iCurAmpList),DetList,FCIDetIndex,iMaxAmpLevel,10)
   LogDealloc(tagSpawnList)
   Deallocate(SpawnList)
   LogDealloc(tagSpawnAmps)
   Deallocate(SpawnAmps)
   call DeallocateAmplitudeList(AL)
   LogDealloc(tagDetList)
   Deallocate(DetList)
   Weight=0.D0
   Energyxw=ProjectionE
   call halt_timer(CCMC_time)
END SUBROUTINE CCMCStandaloneParticle


   SUBROUTINE InitTransitionLog(TL,nAmpl,nMaxClusterSize,tNonUniq)
      use CCMCData, only: CCTransitionLog
      IMPLICIT NONE
      TYPE(CCTransitionLog) TL
      INTEGER nAmpl,nMaxClusterSize
      LOGICAL tNonUniq

      INTEGER i

      INTEGER Binomial
      TL%tNonUniq=tNonUniq
      TL%nExcitors=nAmpl-1

      TL%nMaxSize=min(nMaxClusterSize,TL%nExcitors)
      ! We create a transition matrix where each element correpsonds to a cluster.  We encode cluster X=(a_{X_i}) = (x,y,z)  as an 'bit' string in base Det sum_{i=1}^{|X|} (X_i)*(Det)**(i-1)
      if(tNonUniq) then
         WRITE(6,*) "Using non-unique logging"
         TL%MaxIndex=(nAmpl)**TL%nMaxSize
      else
         WRITE(6,*) "Using unique logging"
         TL%MaxIndex=-1
         do i=0,TL%nMaxSize
            TL%MaxIndex=TL%MaxIndex+Binomial(TL%nExcitors,i)
         enddo
      endif
      WRITE(6,*) "Max Cluster Index:",TL%MaxIndex
      write(6,*) "Logging memory>",(8*2*2*(TL%MaxIndex+1)*(nAmpl+1))/1048576, " Mb"
      allocate(TL%dProbTransition(2,2,0:TL%MaxIndex,0:nAmpl))
      allocate(TL%dProbClust(2,0:TL%MaxIndex))
      allocate(TL%dProbUniqClust(2,-1:TL%MaxIndex))
      TL%dProbTransition(:,:,:,:)=0
      TL%dProbClust(:,:)=0
      TL%dProbUniqClust(:,:)=0
   end subroutine InitTransitionLog
   SUBROUTINE WriteTransitionLog(iUnit,TL,Amplitude,nCycles,dTotAbsAmp,dTotNorm)
      use CCMCData, only: CCTransitionLog
      IMPLICIT NONE
      TYPE(CCTransitionLog) TL
      INTEGER iUnit
      REAL*8 Amplitude(:)
      INTEGER nCycles
      REAL*8 dTotAbsAmp,dTotNorm
      REAL*8 dAveTotAbsAmp,dAveNorm
      INTEGER i,j,k
      REAL*8 r
      REAL*8 flin,flout,flint,floutt
      dAveTotAbsAmp=dTotAbsAmp/nCycles
      dAveNorm=dTotNorm/nCycles
      WRITE(6,*) "Transition Log for last ",nCycles," cycles."
   
      WRITE(iUnit,*) "Cluster Probabilities:  Debiased number per cycle;   Number per cycle"
      do i=0,TL%MaxIndex
         if(TL%dProbClust(2,i).ne.0) then
           Call WriteClusterInd(iUnit,i,.false.,TL)
           WRITE(iUnit,"(2G25.17)") TL%dProbClust(1,i)/nCycles,TL%dProbClust(2,i)/nCycles
         endif
      enddo
      WRITE(iUnit,*) "Unique Cluster Probabilities:  Debiased number per cycle;   Number per cycle"
      do i=0,TL%MaxIndex
         if(TL%dProbUniqClust(2,i).ne.0) then
           Call WriteClusterInd(iUnit,i,.false.,TL)
           WRITE(iUnit,"(2G25.17)") TL%dProbUniqClust(1,i)/nCycles,TL%dProbUniqClust(2,i)/nCycles
         endif
      enddo
      WRITE(iUnit,"(A,2G25.17)") "[Invalid] ",TL%dProbUniqClust(1,-1)/nCycles,TL%dProbUniqClust(2,-1)/nCycles

      WRITE(iUnit,*) "Transition Probabilities: Net value per cycle;   Net number per cycle"
      do i=0,TL%MaxIndex
         flin=0
         flout=0
         do j=0,TL%nExcitors+1 
            if(TL%dProbTransition(1,2,i,j).ne.0) then
               Call WriteClusterInd(iUnit,i,.false.,TL)
               WRITE(iUnit,"(A)",advance='no') "=>"
               Call WriteClusterInd(iUnit,j,.false.,TL)
               WRITE(iUnit,"(2G25.17)") TL%dProbTransition(1,1,i,j)/nCycles,TL%dProbTransition(1,2,i,j)/nCycles
            endif
            flout=flout+TL%dProbTransition(1,1,i,j)/nCycles
         enddo
         if(i<=TL%nExcitors+1) then
            do j=0,TL%MaxIndex
               flin=flin+TL%dProbTransition(1,1,j,i)/nCycles
            enddo
         endif
         flint=flint+flin
         floutt=floutt+flout
         if(flin.ne.0.or.flout.ne.0) then
            Call WriteClusterInd(iUnit,i,.false.,TL)
            write(iUnit,"(A,2G20.10)") "  Flux in, out",flin,flout
         endif
         
      enddo
      WRITE(iUnit,"(A)") "Transition Probabilities normalized by # FromCluster: value per cycle;   number per cycle"
      do i=0,TL%MaxIndex
         do j=0,TL%nExcitors+1 
            if(TL%dProbTransition(1,2,i,j).ne.0) then
               Call WriteClusterInd(iUnit,i,.false.,TL)
               WRITE(iUnit,"(A)",advance='no') "=>"
               Call WriteClusterInd(iUnit,j,.false.,TL)
               if(TL%tNonUniq) THEN
                  WRITE(iUnit,"(2G25.17)") TL%dProbTransition(1,1,i,j)/TL%dProbClust(2,i),TL%dProbTransition(1,2,i,j)/TL%dProbClust(2,i)
               else
                  WRITE(iUnit,"(2G25.17)") TL%dProbTransition(1,1,i,j)/TL%dProbUniqClust(2,i),TL%dProbTransition(1,2,i,j)/TL%dProbUniqClust(2,i)
               endif
            endif
         enddo
      enddo
      WRITE(iUnit,"(A)") "Renormalized Transition Probabilities"
      do i=0,TL%MaxIndex
         do j=0,TL%nExcitors+1 
            if(TL%dProbTransition(1,2,i,j).ne.0) then
               Call WriteClusterInd(iUnit,i,.false.,TL)
               WRITE(iUnit,"(A)",advance='no') "=>"
               Call WriteClusterInd(iUnit,j,.false.,TL)
               k=GetClusterIndLevel(i)
               r=Amplitude(1)*(dAveTotAbsAmp/dAveNorm)**k
               if(TL%tNonUniq) THEN
                  WRITE(iUnit,"(2G25.17)") TL%dProbTransition(2,1,i,j)/(TL%dProbClust(2,i)*r),TL%dProbTransition(2,2,i,j)/(TL%dProbClust(2,i)*r)
               else
                  WRITE(iUnit,"(2G25.17)") TL%dProbTransition(2,1,i,j)/(TL%dProbUniqClust(2,i)*r),TL%dProbTransition(2,2,i,j)/(TL%dProbUniqClust(2,i)*r)
               endif
!                  WRITE(iUnit,"(2G25.17)") TL%dProbTransition(2,1,i,j)/TL%dProbClust(2,i)
            endif
         enddo
      enddo
      WRITE(iUnit,*) "Biased Renormalized Transition Probabilities"
      do i=0,TL%MaxIndex
         do j=0,TL%nExcitors+1 
            if(TL%dProbTransition(1,2,i,j).ne.0) then
               Call WriteClusterInd(iUnit,i,.false.,TL)
               WRITE(iUnit,"(A)",advance='no') "=>"
               Call WriteClusterInd(iUnit,j,.false.,TL)
               k=GetClusterIndLevel(i)
               r=Amplitude(1)*(dAveTotAbsAmp/dAveNorm)**k
               if(TL%tNonUniq) THEN
                  WRITE(iUnit,"(2G25.17)") TL%dProbTransition(1,1,i,j)/(TL%dProbClust(2,i)*r),TL%dProbTransition(1,2,i,j)/(TL%dProbClust(2,i))
               else
                  WRITE(iUnit,"(2G25.17)") TL%dProbTransition(1,1,i,j)/(TL%dProbUniqClust(2,i)*r),TL%dProbTransition(1,2,i,j)/(TL%dProbUniqClust(2,i))
               endif
            endif
         enddo
      enddo
      WRITE(iUnit,"(A)") "Transition values/number transitions: Plain;   Debiased"
      do i=0,TL%MaxIndex
         do j=0,TL%nExcitors+1
            if(TL%dProbTransition(1,2,i,j).ne.0) then
               Call WriteClusterInd(iUnit,i,.false.,TL)
               WRITE(iUnit,"(A)",advance='no') "=>"
               Call WriteClusterInd(iUnit,j,.false.,TL)
               WRITE(iUnit,"(2G25.17)") TL%dProbTransition(1,1,i,j)/TL%dProbTransition(1,2,i,j),TL%dProbTransition(2,1,i,j)/TL%dProbTransition(1,2,i,j)
            endif
         enddo
      enddo
   end subroutine 
   subroutine DeleteTransitionLog(TL)
      use CCMCData, only: CCTransitionLog
      IMPLICIT NONE
      TYPE(CCTransitionLog) TL
      deallocate(TL%dProbClust)
      deallocate(TL%dProbTransition)
   end subroutine 
   subroutine LogTransition(TL,C1,C1size,C2,value,dProbNorm)
      use CCMCData, only: CCTransitionLog
      IMPLICIT NONE
      TYPE(CCTransitionLog) TL
      INTEGER C1size
      INTEGER C1(C1size),C2
      INTEGER i1,i2
      REAL*8 value,dProbNorm 
      i1=GetClusterIndex(C1,C1size,TL)
      i2=C2-1
   ! We create a transition matrix where each element correpsonds to a cluster.  We encode cluster X=(a_{X_i}) = (x,y,z)  as an 'bit' string in base Det sum_{i=1}^{|X|} (X_i)*(Det)**(i-1)
      TL%dProbTransition(1,1,i1,i2)=TL%dProbTransition(1,1,i1,i2)+value
      TL%dProbTransition(1,2,i1,i2)=TL%dProbTransition(1,2,i1,i2)+1
      TL%dProbTransition(2,1,i1,i2)=TL%dProbTransition(2,1,i1,i2)+value*dProbNorm
      TL%dProbTransition(2,2,i1,i2)=TL%dProbTransition(2,2,i1,i2)+1*dProbNorm
!      WRITE(6,"(A)",advance='no') "LT: "
!      Call WriteClusterInd(6,i1,.false.,TL)
!      WRITE(6,"(A)",advance='no') "=>"
!      Call WriteClusterInd(6,i2,.false.,TL)
!      WRITE(6,"(2I4,4G25.17)") i1,i2,value,dProbNorm,TL%dProbTransition(1,1,i1,i2),TL%dProbTransition(1,2,i1,i2)
   end subroutine LogTransition

!Return in Inds(1:iS) sorted excitor indices of a cluster C.  If iS<C%iSize, then the cluster is certainly invalid.
   subroutine GetUniqCluster(C,Inds,iS)
      use CCMCData, only: Cluster
      implicit none
      Type(Cluster) C
      integer Inds(C%iSize)
      integer iS
      integer i,j
!Now make a sorted unique version
      Inds(:)=C%SelectedExcitorIndices(:)
      iS=C%iSize
      call sort (Inds)
      i=2
      j=1
      do while (j.le.iS.and.i.le.C%iSize)
         if(Inds(i)==Inds(j)) then !remove this one (and shorten)
            iS=iS-1
            i=i+1
         else
            j=j+1
            Inds(j)=Inds(i)
            i=i+1
         endif
      enddo
   end subroutine GetUniqCluster
   
   subroutine LogCluster(TL,C)
      use CCMCData, only: CCTransitionLog,Cluster
      implicit None
      TYPE(CCTransitionLog) TL
      Type(Cluster) C
      INTEGER I,iS
      INTEGER Inds(C%iSize)
      if(TL%tNonUniq) then
         i=GetClusterIndex(C%SelectedExcitorIndices(:),C%iSize,TL)
   ! dProbNorm is the prob that a cluster in this level would've been chosen had they been equally weighted
   ! 
   !  dClusterNorm is the probability that this cluster was chosen, given the level had already been selected.
   !  This includes multiple selections of the same excitor as well as combinations of excitors which produce a 0 sign.
         TL%dProbClust(1,i)=TL%dProbClust(1,i)+1/(C%dProbNorm*C%dClusterNorm)
         TL%dProbClust(2,i)=TL%dProbClust(2,i)+1
      endif
      call GetUniqCluster(C,Inds,iS)
      if(iS==C%iSize) then
         i=GetClusterIndex(Inds(:),iS,TL)
      else
         i=-1
      endif
! dProbNorm is the prob that a cluster in this level would've been chosen had they been equally weighted
! 
!  dClusterNorm is the probability that this cluster was chosen, given the level had already been selected.
!  This includes multiple selections of the same excitor as well as combinations of excitors which produce a 0 sign.
      TL%dProbUniqClust(1,i)=TL%dProbUniqClust(1,i)+1/(C%dProbNorm*C%dClusterNorm)
      TL%dProbUniqClust(2,i)=TL%dProbUniqClust(2,i)+1
!      WRITE(6,"(A)",advance='no') "LC: "
!      Call WriteCluster(6,C,.false.)
!      Call WriteClusterInd(6,i,.false.,TL)
!      WRITE(6,"(I)") i
   end subroutine LogCluster

END MODULE CCMC

#ifdef PARALLEL
SUBROUTINE PerformCCMCCycPar
   use CCMC
   Call PerformCCMCCycParInt()
end subroutine PerformCCMCCycPar
#endif

