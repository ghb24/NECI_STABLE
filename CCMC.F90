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
SUBROUTINE AddBitExcitor(iLutnI,iLutnJ,iLutRef,iSgn)
   use SystemData, only : nEl,nIfD
   IMPLICIT NONE
   INTEGER iLutnI(0:nIfD), iLutnJ(0:nIfD),iLutRef(0:nIfD)
   INTEGER iLutTmp(0:nIfD)
   INTEGER T1,T2,T3
   INTEGER iSgn
! We need to run through the bits of J and I concurrently, setting bits of I
   INTEGER i,j

   integer iIlevel,iJlevel,iTmpLevel,iI,iJ
   CALL FindBitExcitLevel(iLutRef,iLutnI,nIfD,iIlevel,nEl)
   CALL FindBitExcitLevel(iLutRef,iLutnJ,nIfD,iJlevel,nEl)
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

!Note that the main routine called by FciMCPar is PerformCCMCCycPar and is outside the module, so
!FciMCPar doesn't need to use the module
    
!Based on PerformCleanFCIMCycPar()

!First we loop over particles at the HF det, and for each, we select a group of excitors to make a composite excitor.
! From this we make an excitation, and accept or reject a particle created there.

! Next we decide whether or not that composite excitor should die.  Its death is caused by the death of one of its composites.  We save a list of the excitors which need to die.

!Then we go through the list of excitors which should die, and reduce their weight on the list, and remove if necessary.
!Next we annihilate.

    SUBROUTINE PerformCCMCCycPar()
      USe FCIMCParMod
      use CCMCData
      Use Determinants, only: GetHElement3
        IMPLICIT NONE
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
        REAL*8 dProbNorm,dProb
        INTEGER iExcitor
        INTEGER iMinEx,iMaxEx
!iMaxExcitorSelections is the number of times we decide to loop for each particle.
!  it includes the first (exact) loop where we consider the particle by itself.
        INTEGER iMaxExcitorSelections
        LOGICAL tSuccess

        REAL*8 dT1Sq
        INTEGER AttemptDieProbPar
        REAL*8 AJWTProjE
  

        if(Iter.eq.1) dT1SqCuml=0 
        iMaxExcitorSelections=1
        iDebug=4
!Number of excitors dying
        iDeaths=0

        AJWTProjE=0
        dT1Sq=0
        IF(iDebug.gt.0) WRITE(6,*) "Entering CCMC Cycle"
        iHFDet=1
        dProbSelNewExcitor=0.5


        IF(TDebug.and.(mod(Iter,10).eq.0)) THEN
            WRITE(11,*) Iter,TotWalkers,NoatHF,NoatDoubs,MaxIndex,TotParts
            CALL FLUSH(11)
        ENDIF
        
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
        if(.not.tSuccess) STOP "Cannot find HF det in particle list"
        HFcount=abs(CurrentSign(iHFDet))

        allocate(iKillDetIndices(2,TotParts*2))

        IF(iDebug.gt.1) THEN
         write(6,*) "HF det"
         call WriteBitDet(6,iLutHF,.true.)
         write(6,*) "Particle list"
         do j=1,TotWalkers
          write(6,'(i)',advance='no') CurrentSign(j)
          call WriteBitEx(6,iLutHF,CurrentDets(:,j),.true.)
         enddo
        endif
!TotRealWalkers gets updated with the number of non-zero walkers at each stage.
        TotRealWalkers=TotWalkers
        do j=1,TotWalkers
!#if 0
          write(89,'(I)',advance='no') CurrentSign(j)
          IF(iDebug.gt.4) WRITE(6,*) "Iteration ",Iter,':',j
! Deal with T_1^2
          CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,2)
          if(WalkExcitLevel.eq.1) then
            do l=j,TotWalkers
               CALL FindBitExcitLevel(iLutHF,CurrentDets(:,l),NIfD,WalkExcitLevel,2)
               if(WalkExcitLevel.eq.1) then
                  iSgn=CurrentSign(j)
                  iSgn=iSgn*CurrentSign(l)
                  iLutnI(:)=CurrentDets(:,j)
                  call AddBitExcitor(iLutnI,CurrentDets(:,l),iLutHF,iSgn)
                  if(iSgn.ne.0) then
                     CALL DecodeBitDet(DetCurr,iLutnI(:),NEl,NIfD)
                     Htmp=GetHElement3(HFDet, DetCurr,2)
                     dT1Sq=dT1Sq+(Real(Htmp%v,r2)*iSgn)
                     WRITE(6,'(A,I,2G)', advance='no') 'T1',iSgn,real(Htmp%v,r2),dT1Sq
                     call WriteBitEx(6,iLutHF,CurrentDets(:,j),.false.)
                     call WriteBitEx(6,iLutHF,CurrentDets(:,l),.false.)
                     call WriteBitEx(6,iLutHF,iLutnI,.true.)
                  endif
               endif
            enddo
          endif 
          if(WalkExcitLevel.eq.1.or.WalkExcitLevel.eq.2) then
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j),NEl,NIfD)
            Htmp=GetHElement3(HFDet, DetCurr,WalkExcitLevel)
            AJWTProjE=AJWTProjE+(Real(Htmp%v,r2)*CurrentSign(j))
          endif
!#endif
!Go through all particles on the current walker det
          do l=1,abs(CurrentSign(j))
            if(j.eq.iHFDet) then
               iMaxEx=1
            else
               iMaxEx=iMaxExcitorSelections
            endif
! tNewExcitor will be set at the end of this loop if we decide to come 
            do iExcitor=1,iMaxEx
               if(iExcitor.eq.1) then
!just select the single excitor
                  iLutnI(:)=CurrentDets(:,j)
                  iSgn=sign(1,CurrentSign(j))
                  dProb=1 
                  dProbNorm=1
                  iCompositeSize=0
               else
                  stop
! Now select a sample of walkers.  We need up to as many walkers as electrons.
                  dProbNumExcit=1
!dProbSelNewExcitor   !prob of this number of excitors will go here
                  IF (tTruncSpace) THEN
                     nMaxSelExcitors=ICILevel
                  ELSE
                     nMaxSelExcitors=nEl
                  ENDIF
            
                  if(nMaxSelExcitors.gt.(TotWalkers-1)) nMaxSelExcitors=TotRealWalkers-1
                  IF(iDebug.gt.4) write(6,*) "Max excitors can be selected.", nMaxSelExcitors
                  if(nMaxSelExcitors.lt.2) exit  !If we can't choose a new excit, we leave the loop
                  dProb=1           !prob of these excitors given this number of excitors goes here
   !NB it's possible to select no excitors with this loop.
                  dProbNorm=1
!We force the first excitor to be the current det.
                  SelectedExcitorIndices(1)=j
                  SelectedExcitors(:,1)=CurrentDets(:,j)

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
                     IF(iDebug.gt.5) Write(6,*) "Selected excitor",k
                     SelectedExcitorIndices(i)=k
                     SelectedExcitors(:,i)=CurrentDets(:,k)
                     dProb=(dProb*abs(CurrentSign(k)))/((TotParts-HFcount))

!Account for intermediate Normalization
                     dProbNorm=dProbNorm*HFCount

!Account for possible orderings of selection
                     dProbNorm=dProbNorm/i

                     IF(iDebug.gt.5) WRITE(6,*) "TotParts,HFCount:",TotParts,HFcount
                     IF(iDebug.gt.5) write(6,*) "Prob ",i,": ",(abs(CurrentSign(k))+0.d0)/(TotParts-HFcount)," Cuml:", dProb
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
                  IF(iDebug.gt.3) then
                     WRITE(6,*) " Excitors in composite:", iCompositeSize
                     do i=1,iCompositeSize
                        write(6,'(I2)',advance='no') sign(1,CurrentSign(SelectedExcitorIndices(i)))
                        call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.true.)
                     enddo
                     Write(6,*) "Level chosen Prob      : ",dProbNumExcit
                     Write(6,*) "Select Prob given level: ",dProb
                     Write(6,*) "Prob norm              : ",dProbNorm
                  endif
                  dProb=dProbNumExcit*dProb
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
               CALL DecodeBitDet(DetCurr,iLutnI(:),NEl,NIfD)



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
               Prob=Prob*dProb*dProbNorm  !Now include the prob of choosing the det we spawned from
               if(iDebug.gt.4) Write(6,*) "Prob ex tot",Prob
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

               IF(iDebug.gt.4.or.((iDebug.eq.3.or.iDebug.eq.4).and.Child.ne.0)) THEN
   !We've not printed this out before
                  call WriteBitEx(6,iLutHF,iLutnI,.false.)
                  write(6,'(A)',advance='no') ,' ==> '
                  call WriteBitEx(6,iLutHF,iLutnJ,.false.)
                  WRITE(6,'(A,I)',advance='no') "Children:",Child             
                  if(iDebug.eq.3.and.iCompositeSize.gt.1) THEN
                     write(6,'(A)',advance='no') ' from '
                     do i=1,iCompositeSize
                        call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.false.)
                     enddo
                  endif
                  write(6,*)
               endif

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
                  dProbDecompose=1.0/iCompositeSize
                  call genrand_real2(r)  !On GHB's advice
                  k=1+floor(r*iCompositeSize)
                  iPartDie=SelectedExcitorIndices(k)
   !Now get the full representation of the dying excitor
                  CALL DecodeBitDet(DetCurr,iLutnI,NEl,NIfD)
                  Htmp=GetHElement3(DetCurr,DetCurr,0)
                  HDiagCurr=REAL(Htmp%v,r2)
                  HDiagCurr=HDiagCurr-Hii
                  iSgn=sign(1,CurrentSign(iPartDie))
               ELSE
                  dProbDecompose=1
                  iPartDie=j
                  HDiagCurr=CurrentH(j)
!                  CALL SumEContrib(DetCurr,WalkExcitLevel,iSgn,iLutnI,HDiagCurr,1.d0)
               ENDIF 

!Also, we want to find out the excitation level of the determinant - we only need to find out if its connected or not (so excitation level of 3 or more is ignored.
!This can be changed easily by increasing the final argument.
               IF(tTruncSpace) THEN
!We need to know the exact excitation level for truncated calculations.
                   CALL FindBitExcitLevel(iLutHF,iLutnI,NIfD,WalkExcitLevel,NEl)
               ELSE
                   CALL FindBitExcitLevel(iLutHF,iLutnI,NIfD,WalkExcitLevel,2)
               ENDIF
!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info
!               WRITE(6,'(A,3G,2I)') "PPP",dProb,dProbNorm,dProb*dProbNorm,WalkExcitLevel,iCompositeSize
!               CALL SumEContrib(DetCurr,WalkExcitLevel,iSgn,iLutnI,HDiagCurr,(dProb*dProbNorm))



               IF(iDebug.gt.4) then
                  Write(6,'(A,I)',advance='no') "Killing at excitor: ",iPartDie
                  call WriteBitEx(6,iLutHF,CurrentDets(:,iPartDie),.false.)
                  Write(6,'(A,G)') "Prob: ",dProb*dProbNorm
               endif

   !Die with a probability as normal, but biased because we only selected this set
   !of excitors with dProb.
               iDie=AttemptDieProbPar(DetCurr,HDiagCurr,WalkExcitLevel,iSgn,1/(dProb*dProbNorm))
               IF(iDebug.gt.4.or.((iDebug.eq.3.or.iDebug.eq.4).and.iDie.ne.0)) then
                  Write(6,'(A,I)',advance='no') " Killing at excitor: ",iPartDie
                  Write(6,'(A)',advance='no') " chosen "
                  call WriteBitEx(6,iLutHF,CurrentDets(:,iPartDie),.false.)
                  WRITE(6,'(A,I)',advance='no') "Number died ",iDie
                  if(iCompositeSize.gt.1) then
                     WRITE(6,'(A)',advance='no') " from "
                     do i=1,iCompositeSize
                        call WriteBitEx(6,iLutHF,SelectedExcitors(:,i),.false.)
                     enddo
                  endif
                  WRITE(6,*)
                  if(iDebug.eq.4) Write(6,'(A,G)') "Prob: ",dProb*dProbNorm
               endif
               NoDied=NoDied+iDie          !Update death counter
               
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
        write(89,*)
!#if 0
        dT1Sq=(dT1Sq/HFCount)
        dT1SqCuml=dT1SqCuml+dT1Sq
        if (iMaxExcitorSelections.gt.1) then
!            ENumCyc=ENumCyc+dT1Sq
!            IF(Iter.gt.NEquilSteps) SumENum=SumENum+dT1Sq
        endif
        write(78,'(I,4G)') Iter, dT1Sq/HFCount,AJWTProjE/HFCount,AJWTProjE,dT1SqCuml/Iter
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
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,NEl)
            ELSE
                CALL FindBitExcitLevel(iLutHF,CurrentDets(:,j),NIfD,WalkExcitLevel,2)
            ENDIF
!HDiags are stored.
            HDiagCurr=CurrentH(j)
            CALL DecodeBitDet(DetCurr,CurrentDets(:,j),NEl,NIfD)

!Sum in any energy contribution from the determinant, including other parameters, such as excitlevel info

            CALL SumEContrib(DetCurr,WalkExcitLevel,CurrentSign(j),CurrentDets(:,j),HDiagCurr,1.d0)
            CopySign=CurrentSign(j)
            IF(CopySign.ne.0) THEN
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


#endif
