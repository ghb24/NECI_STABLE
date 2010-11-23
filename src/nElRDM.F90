MODULE nElRDMMod
! This file contains the routines used to find the full n electron reduced density matrix (nElRDM).
! This is done on the fly to avoid having to histogram the full wavefunction which is extremely time and memory inefficient.
! In this way, these routines differ slightly from those in NatOrbsMod (which take a histogrammed wavefunction usually 
! truncated around double excitations and form the one electron RDM) but the basic formula is still the same.
! For example, he elements of the one electron reduced density matrix are given by:
! 1RDM_pq   = < Psi | a_p+ a_q | Psi > 
! where Psi is the full wavefunction, a_p+ is the creation operator and a_q is the annihilation operator.
!           = < sum_i c_i D_i | a_p+ a_q | sum_j c_j D_j >
! The elements 1RDM_pq therefore come from the sum of the contributions c_i*c_j from all pairs of determinants D_i and D_j which 
! are related by a single excitation between p and q.
! This can be generalised for the nElRDM, where n = 1 or 2.

! The algorithm for calculating the nElRDM on the fly will be similar in nature to the direct annihilation routines.
! Each set of processors has a list of determinants.  
! These each select the first of these D_i, generate all the allowed single or double excitations D_j, and order them in terms of the 
! processor they would be on if they are occupied.
! The excitations are then sent to the relevant processor along with the original determinant, and it's sign c_i.
! Each processor then receives nProcessors sets of excitations (D_j's).
! For each of these they binary search their list of occupied determinants.  
! If an excitation D_j is found, c_i.c_j is added to the matrix element corresponding to the orbitals involved in the excitation.

! NOTE: There will be possible speed ups considering the fact that the 1RDM is symmetrical.
! Can initially find all elements and average the two values pq and qp (more accurate??).
! But should put a condition into the excitaiton generator so that only single excitations with q > p are generated.

! By finding the full 1RDM, we have the ability to derive the natural orbitals as well as electron densities etc.
        
        USE Global_Utilities
        USE Parallel
        USE bit_reps , only : NIfTot
        USE SystemData , only : NEl,nBasis,tStoreSpinOrbs
        USE NatOrbsMod , only : NatOrbMat,NatOrbMatTag,Evalues,EvaluesTag
        USE CalcData , only : MemoryFacPart
        USE constants , only : n_int
        USE Logging , only : RDMExcitLevel
        USE RotateOrbsData , only : CoeffT1, CoeffT1Tag, tTurnStoreSpinOff
        IMPLICIT NONE
        INTEGER , ALLOCATABLE :: Sing_InitExcSlots(:),Sing_ExcList(:)
        INTEGER , ALLOCATABLE :: Doub_InitExcSlots(:),Doub_ExcList(:)
        INTEGER(kind=n_int) , ALLOCATABLE :: Sing_ExcDjs(:,:),Sing_ExcDjs2(:,:)
        INTEGER(kind=n_int) , ALLOCATABLE :: Doub_ExcDjs(:,:),Doub_ExcDjs2(:,:)
        INTEGER :: Sing_ExcDjsTag,Sing_ExcDjs2Tag,TwoElRDMTag,AllTwoElRDMTag
        INTEGER :: Doub_ExcDjsTag,Doub_ExcDjs2Tag,OneElRDMTag
        REAL*8 , ALLOCATABLE :: OneElRDM(:,:)
        REAL*8 , ALLOCATABLE :: TwoElRDM(:,:)
        REAL*8 , ALLOCATABLE :: AllTwoElRDM(:,:)
        REAL*8 :: OneEl_Gap,TwoEl_Gap,AllTotPartstemp
        type(timer), save :: nElRDM_Time

    contains

    SUBROUTINE InitRDM()
! This routine initialises any of the arrays needed to calculate the reduced density matrix.    
        USE NatOrbsMod , only : SetupNatOrbLabels 
        USE SystemData , only : tSeparateOccVirt
        USE RotateOrbsMod , only : SymLabelCounts2,SymLabelCounts2Tag,NoOrbs,SpatOrbs,NoRotOrbs
        USE RotateOrbsMod , only : SymLabelList2,SymLabelListInv,SymLabelList3,SymLabelList2Tag
        USE RotateOrbsMod , only : SymLabelList3Tag,SymLabelListInvTag
        USE Logging , only : tCalc_RDMEnergy, tDiagRDM
        INTEGER :: ierr,i
        CHARACTER(len=*), PARAMETER :: this_routine='InitRDM'

        IF(tCalc_RDMEnergy) THEN
            WRITE(6,*) 'Calculating the energy from the reduced density matrix, this required both the 1 and 2 electron RDM.'
            RDMExcitLevel = 3
        ENDIF

        IF(tDiagRDM.and.(RDMExcitLevel.eq.2)) THEN
            WRITE(6,*) 'WARNING : Requesting to diagonalise 1-RDM, but only the 2-RDM is being calculated.'
        ENDIF

! The stuff below here is so that we can set up the symmetry arrays according to the NatOrbs routines.
! This then allows us to keep the one electron reduced density matrix in symmetry blocks, and therefore diagonalise it still in those 
! symmetry blocks later on.
! Don't bother when we're getting the two electron reduced density matrix.
        tTurnStoreSpinOff=.false.
        IF(.not.tStoreSpinOrbs) THEN
            tTurnStoreSpinOff=.true.
            tStoreSpinOrbs=.true.
        ENDIF
! We need tStoreSpinOrbs to be temporarily turned on because we always want OneRDM to be in spin orbitals.        
        tSeparateOccVirt=.false.
        NoOrbs=nBasis
        SpatOrbs=nBasis/2
        NoRotOrbs=NoOrbs

        IF((RDMExcitLevel.eq.2).or.(RDMExcitLevel.eq.3)) THEN

! This array actually contains the excitations in blocks of the processor they will be sent to.        
            ALLOCATE(Doub_ExcDjs(0:NIfTot,NINT(((NEl*nBasis)**2)*MemoryFacPart)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Doub_ExcDjs array.')
            CALL LogMemAlloc('Doub_ExcDjs',NINT(((NEl*nBasis)**2)*MemoryFacPart)*(NIfTot+1),size_n_int,this_routine,Doub_ExcDjsTag,ierr)

            ALLOCATE(Doub_ExcDjs2(0:NIfTot,NINT(((NEl*nBasis)**2)*MemoryFacPart)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Doub_ExcDjs2 array.')
            CALL LogMemAlloc('Doub_ExcDjs2',NINT(((NEl*nBasis)**2)*MemoryFacPart)*(NIfTot+1),size_n_int,this_routine,Doub_ExcDjs2Tag,ierr)

! We also need to allocate the actual nElRDM on each processor, and an allnElRDM on only the root.
            ALLOCATE(TwoElRDM(((nBasis*(nBasis-1))/2),((nBasis*(nBasis-1))/2)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating TwoElRDM array,')
            CALL LogMemAlloc('TwoElRDM',(((nBasis*(nBasis-1))/2)**2),8,this_routine,TwoElRDMTag,ierr)
            TwoElRDM(:,:)=0.D0

!            ALLOCATE(TwoElRDM(nBasis,nBasis,nBasis,nBasis),stat=ierr)
!            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating TwoElRDM array,')
!            CALL LogMemAlloc('TwoElRDM',(nBasis**4),8,this_routine,TwoElRDMTag,ierr)
!            TwoElRDM(:,:,:,:)=0.D0

            IF(iProcIndex.eq.0) THEN
                ALLOCATE(AllTwoElRDM(((nBasis*(nBasis-1))/2),((nBasis*(nBasis-1))/2)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating AllTwoElRDM array,')
                CALL LogMemAlloc('AllTwoEleRDM',(((nBasis*(nBasis-1))/2)**2),8,this_routine,AllTwoElRDMTag,ierr)
                AllTwoElRDM(:,:)=0.D0
            ENDIF


!            IF(iProcIndex.eq.0) THEN
!                ALLOCATE(AllTwoElRDM(nBasis,nBasis,nBasis,nBasis),stat=ierr)
!                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating AllTwoElRDM array,')
!                CALL LogMemAlloc('AllTwoEleRDM',(nBasis**4),8,this_routine,AllTwoElRDMTag,ierr)
!                AllTwoElRDM(:,:,:,:)=0.D0
!            ENDIF

! We need room to potentially generate (N*M)^2 double excitations but these will be spread across each processor.        
            TwoEl_Gap=(((REAL(NEl)*REAL(nBasis))**2)*MemoryFacPart)/REAL(nProcessors)

            Doub_ExcDjs(:,:)=0
            Doub_ExcDjs2(:,:)=0

! This array contains the initial positions of the excitations for each processor.
            ALLOCATE(Doub_InitExcSlots(0:(nProcessors-1)),stat=ierr)
            do i=0,nProcessors-1
                Doub_InitExcSlots(i)=NINT(TwoEl_Gap*i)+1
            enddo

! This array contains the current position of the excitations as they're added.
            ALLOCATE(Doub_ExcList(0:(nProcessors-1)),stat=ierr)
            Doub_ExcList(:)=Doub_InitExcSlots(:)

        ENDIF            

        IF((RDMExcitLevel.eq.1).or.(RDMExcitLevel.eq.3)) THEN

! This array actually contains the excitations in blocks of the processor they will be sent to.        
            ALLOCATE(Sing_ExcDjs(0:NIfTot,NINT((NEl*nBasis)*MemoryFacPart)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Sing_ExcDjs array.')
            CALL LogMemAlloc('Sing_ExcDjs',NINT(NEl*nBasis*MemoryFacPart)*(NIfTot+1),size_n_int,this_routine,Sing_ExcDjsTag,ierr)

            ALLOCATE(Sing_ExcDjs2(0:NIfTot,NINT((NEl*nBasis)*MemoryFacPart)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Sing_ExcDjs2 array.')
            CALL LogMemAlloc('Sing_ExcDjs2',NINT(NEl*nBasis*MemoryFacPart)*(NIfTot+1),size_n_int,this_routine,Sing_ExcDjs2Tag,ierr)


! We also need to allocate the actual nElRDM on each processor, and an allnElRDM on only the root.
            ALLOCATE(OneElRDM(nBasis,nBasis),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating OneElRDM array,')
            CALL LogMemAlloc('nElRDM',nBasis**2,8,this_routine,OneElRDMTag,ierr)
            OneElRDM(:,:)=0.D0

            IF(iProcIndex.eq.0) THEN
! This is the AllnElRDM, called NatOrbMat simply because we use the natural orbital routines to diagonalise etc - am gonna change this so it's passed around).        
                ALLOCATE(NatOrbMat(nBasis,nBasis),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating NatOrbMat array,')
                CALL LogMemAlloc('NatOrbMat',nBasis**2,8,this_routine,NatOrbMatTag,ierr)
                NatOrbMat(:,:)=0.D0
            ENDIF

            Sing_ExcDjs(:,:)=0
            Sing_ExcDjs2(:,:)=0

! We need room to potentially generate N*M single excitations but these will be spread across each processor.        

            OneEl_Gap=(REAL(NEl)*REAL(nBasis)*MemoryFacPart)/REAL(nProcessors)

! This array contains the initial positions of the excitations for each processor.
            ALLOCATE(Sing_InitExcSlots(0:(nProcessors-1)),stat=ierr)
            do i=0,nProcessors-1
                Sing_InitExcSlots(i)=NINT(OneEl_Gap*i)+1
            enddo

! This array contains the current position of the excitations as they're added.
            ALLOCATE(Sing_ExcList(0:(nProcessors-1)),stat=ierr)
            Sing_ExcList(:)=Sing_InitExcSlots(:)


            ALLOCATE(SymLabelCounts2(2,32),stat=ierr)
            CALL LogMemAlloc('SymLabelCounts2',2*32,4,this_routine,SymLabelCounts2Tag,ierr)
            SymLabelCounts2(:,:)=0

            ALLOCATE(SymLabelList2(NoOrbs),stat=ierr)
            CALL LogMemAlloc('SymLabelList2',NoOrbs,4,this_routine,SymLabelList2Tag,ierr)
            SymLabelList2(:)=0                     
            ALLOCATE(SymLabelList3(NoOrbs),stat=ierr)
            CALL LogMemAlloc('SymLabelList3',NoOrbs,4,this_routine,SymLabelList3Tag,ierr)
            SymLabelList3(:)=0                     
     
            ALLOCATE(SymLabelListInv(NoOrbs),stat=ierr)
            CALL LogMemAlloc('SymLabelListInv',NoOrbs,4,this_routine,SymLabelListInvTag,ierr)
            SymLabelListInv(:)=0   

            IF(iProcIndex.eq.0) THEN
                ALLOCATE(Evalues(NoOrbs),stat=ierr)
                CALL LogMemAlloc('Evalues',NoOrbs,8,this_routine,EvaluesTag,ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for Evalues failed.")
                Evalues(:)=0.D0
            ENDIF

            CALL SetupNatOrbLabels() 

        ENDIF            

        !Need to now turn tStoreSpinOrbs back off for the rest of the calculation, but all the 
        !symmetry arrays etc are set up for tStoreSpinOrbs to be on.
        !Note, if this has been meddled with, the UMAT will be in spatial orbitals and everything else in spin.
        IF(tTurnStoreSpinOff) tStoreSpinOrbs=.false.

        nElRDM_Time%timer_name='nElRDMTime'


    END SUBROUTINE InitRDM


    SUBROUTINE FillRDMthisIter(TotWalkers)
        USE FciMCData , only : CurrentDets,TotParts 
        INTEGER(int64) , INTENT(IN) :: TotWalkers
        INTEGER(kind=n_int) :: iLutnI(0:NIfTot)
        INTEGER(int64) :: MaxTotWalkers,TotWalkIn(2),TotWalkOut(2)
        INTEGER :: i,error
        REAL*8 :: TempTotParts
        LOGICAL :: blank_det


        CALL set_timer(nElRDM_Time,30)

! Run through the current determinants.
! Find the max number of determinants on a processor - all need to run through this number so that the communication can be done at all stages.

        TotWalkIn(1)=TotWalkers
        TotWalkIn(2)=iProcIndex
        TempTotParts=REAL(TotParts(1))

        CALL MPIAllReduceDatatype(TotWalkIn,1,MPI_MAXLOC,MPI_2INTEGER,TotWalkOut)
        CALL MPIAllReduce(TempTotParts,MPI_SUM,AllTotPartstemp)

        MaxTotWalkers=TotWalkOut(1)

        do i=1,MaxTotWalkers

! But if the actual number of determinants on this processor is less than the number we're running through, feed in 0 determinants and 0 sign.
            IF(i.gt.TotWalkers) THEN
                iLutnI(:)=0
                blank_det=.true.
            ELSE
                iLutnI(:)=CurrentDets(:,i)
                blank_det=.false.
            ENDIF

            CALL AddRDMContrib(iLutnI,blank_det)

        enddo

        CALL halt_timer(nElRDM_Time)


    END SUBROUTINE FillRDMthisIter



    SUBROUTINE AddRDMContrib(iLutnI,blank_det)
! This is the general routine for taking a particular determinant in the spawned list, D_i and adding it's contribution to 
! the reduced density matrix.
        INTEGER(kind=n_int), INTENT(IN) :: iLutnI(0:NIfTot)
        LOGICAL, INTENT(IN) :: blank_det
        INTEGER :: i
        
! Set up excitation arrays.
! These are blocked according to the processor the excitation would be on if occupied.
! In each block, the first entry is the sign of determinant D_i and the second the bit string of the determinant (these need 
! to be sent along with the excitations).
! Each processor will have a different Di.

        IF((RDMExcitLevel.eq.1).or.(RDMExcitLevel.eq.3)) THEN
            Sing_ExcDjs(:,:)=0
            Sing_ExcList(:)=0
            Sing_ExcList(:) = Sing_InitExcSlots(:)

            do i=0,nProcessors-1
                Sing_ExcDjs(:,Sing_ExcList(i)) = iLutnI(:)
                Sing_ExcList(i) = Sing_ExcList(i)+1
            enddo
        ENDIF
        IF((RDMExcitLevel.eq.2).or.(RDMExcitLevel.eq.3)) THEN
            Doub_ExcDjs(:,:) = 0
            Doub_ExcList(:) = 0
            Doub_ExcList(:) = Doub_InitExcSlots(:)

            do i=0,nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(i)) = iLutnI(:)
                Doub_ExcList(i) = Doub_ExcList(i)+1
            enddo
        ENDIF

        IF(.not.blank_det) CALL GenExcDjs(iLutnI)
! Out of here we will get a filled ExcDjs array with all the single or double excitations from Dj, this will be done for each proc. 

! We then need to send the excitations to the relevant processors.
        CALL SendProcExcDjs()
! This routine then calls SearchOccDets which takes each excitation and and binary searches the occupied determinants for this.
! If found, we re-find the orbitals and parity involved in the excitation, and add the c_i*c_j contributions to 
! the corresponding matrix element.

    END SUBROUTINE AddRDMContrib



    SUBROUTINE GenExcDjs(iLutnI)
! This uses GenExcitations3 in symexcit3.F90 to generate all the possible either single or double excitations from D_i, 
! finds the processor they would be on if occupied, and puts them in the SingExcDjs array according to that processor.
        USE DetBitOps , only : EncodeBitDet
        USE AnnihilationMod , only : DetermineDetNode
        USE SymExcit3 , only : GenExcitations3
        USE RotateOrbsData , only : SymLabelListInv
        USE bit_reps , only : extract_bit_rep
        INTEGER(kind=n_int) , INTENT(IN) :: iLutnI(0:NIfTot)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: SignDi
        INTEGER :: ExcitMat3(2,2),nI(NEl),nJ(NEl),Proc,i,j,FlagsDi,Ind,a,b
        LOGICAL :: tParity,tAllExcitFound


        call extract_bit_rep (iLutnI, nI, SignDi, FlagsDi)
! Unfortunately uses the decoded determinant - might want to look at this.        

! Need to add in the diagonal elements.
! Do this now, while it is decoded.
! The RDM are always in spin orbitals, so just adding the orbital as is, is fine.
        IF((RDMExcitLevel.eq.1).or.(RDMExcitLevel.eq.3)) THEN
            do i=1,NEl
                OneElRDM(SymLabelListInv(nI(i)),SymLabelListInv(nI(i))) = &
                            OneElRDM(SymLabelListInv(nI(i)),SymLabelListInv(nI(i))) &
                             +((REAL(SignDi(1))*REAL(SignDi(1)))/AllTotPartstemp)
            enddo
        ENDIF
        IF((RDMExcitLevel.eq.2).or.(RDMExcitLevel.eq.3)) THEN
            do i=1,NEl
                do j=i+1,NEl
                    Ind=( ( (nI(j)-2) * (nI(j)-1) ) / 2 ) + nI(i)
                    TwoElRDM( Ind , Ind ) = TwoElRDM( Ind , Ind ) &
                             +( (REAL(SignDi(1))*REAL(SignDi(1)) ) / AllTotPartstemp )
                enddo
            enddo
        ENDIF
!        IF((RDMExcitLevel.eq.2).or.(RDMExcitLevel.eq.3)) THEN
!            do i = 1, NEl
!                do j = 1, NEl
!                    do a = 1, NEl
!                        do b = 1, NEl
!                            IF((i.eq.a).and.(j.eq.b)) THEN
!                                TwoElRDM(i,j,a,b) = TwoElRDM(i,j,a,b) &
!                                     +( (REAL(SignDi(1))*REAL(SignDi(1)) ) / AllTotPartstemp )
!                            ENDIF
!                            IF((i.eq.b).and.(j.eq.a)) THEN
!                                TwoElRDM(i,j,a,b) = TwoElRDM(i,j,a,b) &
!                                     -( (REAL(SignDi(1))*REAL(SignDi(1)) ) / AllTotPartstemp )
!                            ENDIF
!                        enddo
!                    enddo
!                enddo
!            enddo
!        ENDIF


        IF((RDMExcitLevel.eq.1).or.(RDMExcitLevel.eq.3)) THEN

            ExcitMat3(:,:)=0
! Zeros in ExcitMat3 starts off at the first single excitation.        
            tAllExcitFound=.false.
! This becomes true when all the excitations have been found.        

            do while (.not.tAllExcitFound)
                CALL GenExcitations3(nI,iLutnI,nJ,1,ExcitMat3(:,:),tParity,tAllExcitFound)            
! Passed out of here is the singly excited determinant, nJ.
! Information such as the orbitals involved in the excitation and the parity is also found in this step,
! we are not currently storing this, and it is re-calculated later on (after the determinants are passed
! to the relevant processor) - but the speed of sending this information vs recalculating it will be tested.
! RDMExcitLevel is passed through, if this is 1, only singles are generated, if it is 2 only doubles are found.

                IF(tAllExcitFound) EXIT

                iLutnJ(:)=0
                CALL EncodeBitDet(nJ,iLutnJ)

                Proc = DetermineDetNode(nJ,0)   !This will return a value between 0 -> nProcessors-1
                Sing_ExcDjs(:,Sing_ExcList(Proc)) = iLutnJ(:)
                Sing_ExcList(Proc) = Sing_ExcList(Proc)+1

! Want a quick test to see if arrays are getting full.            
                IF(Sing_ExcList(Proc).gt.NINT(OneEl_Gap*(Proc+1))) THEN
                    WRITE(6,*) 'Proc',Proc
                    WRITE(6,*) 'Sing_ExcList',Sing_ExcList
                    WRITE(6,*) 'No. spaces for each proc',NINT(OneEl_Gap)
                    CALL Stop_All('GenExcDjs','Too many excitations for space available.')
                ENDIF
            enddo
        ENDIF

        IF((RDMExcitLevel.eq.2).or.(RDMExcitLevel.eq.3)) THEN

            ExcitMat3(:,:)=0
! Zeros in ExcitMat3 starts off at the first single excitation.        
            tAllExcitFound=.false.
! This becomes true when all the excitations have been found.        

            do while (.not.tAllExcitFound)
                CALL GenExcitations3(nI,iLutnI,nJ,2,ExcitMat3(:,:),tParity,tAllExcitFound)            
! Passed out of here is the singly excited determinant, nJ.
! Information such as the orbitals involved in the excitation and the parity is also found in this step,
! we are not currently storing this, and it is re-calculated later on (after the determinants are passed
! to the relevant processor) - but the speed of sending this information vs recalculating it will be tested.
! RDMExcitLevel is passed through, if this is 1, only singles are generated, if it is 2 only doubles are found.

                IF(tAllExcitFound) EXIT

                iLutnJ(:)=0
                CALL EncodeBitDet(nJ,iLutnJ)

                Proc = DetermineDetNode(nJ,0)   !This will return a value between 0 -> nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(Proc)) = iLutnJ(:)
                Doub_ExcList(Proc) = Doub_ExcList(Proc)+1

! Want a quick test to see if arrays are getting full.            
                IF(Doub_ExcList(Proc).gt.NINT(TwoEl_Gap*(Proc+1))) THEN
                    WRITE(6,*) 'Proc',Proc
                    WRITE(6,*) 'Doub_ExcList',Doub_ExcList
                    WRITE(6,*) 'No. spaces for each proc',NINT(TwoEl_Gap)
                    CALL Stop_All('GenExcDjs','Too many excitations for space available.')
                ENDIF
            enddo
        ENDIF


    END SUBROUTINE GenExcDjs

    

    SUBROUTINE SendProcExcDjs()
! In this routine the excitations are sent to the relevant processors.
! Sent with them will be the Di they were excited from and its sign.
! Each processor will receive nProcessor number of lists with different Di determinants.
! The original Di's will (I think) still be in the original InitSingExcSlots positions.
! This follows the directannihilation algorithm closely.
        INTEGER :: i,j,sendcounts(nProcessors),disps(nProcessors),sing_recvcounts(nProcessors)
        INTEGER :: sing_recvdisps(nProcessors),error,MaxSendIndex,MaxIndex
        INTEGER :: doub_recvcounts(nProcessors),doub_recvdisps(nProcessors)

        IF((RDMExcitLevel.eq.1).or.(RDMExcitLevel.eq.3)) THEN
            do i=0,nProcessors-1
                sendcounts(i+1)=Sing_ExcList(i)-(NINT(OneEl_Gap*i)+1)
! Sendcounts is the number of singly excited determinants we want to send for each processor (but goes from 1, not 0).            
                disps(i+1)=NINT(OneEl_Gap*i)
! and I think disps is the first slot for each processor - 1.            
            enddo

            MaxSendIndex=Sing_ExcList(nProcessors-1)-1

! We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
            sing_recvcounts(1:nProcessors)=0
            CALL MPIAlltoAll(sendcounts,1,sing_recvcounts,1,error)
! I think recvcounts(i) is the number of determinants sent from processor i.        

! We can now get recvdisps from recvcounts, since we want the data to be contiguous after the move.
            sing_recvdisps(1)=0
            do i=2,nProcessors
                sing_recvdisps(i)=sing_recvdisps(i-1)+sing_recvcounts(i-1)
            enddo

            MaxIndex=sing_recvdisps(nProcessors)+sing_recvcounts(nProcessors)
! But the actual number of integers we need to send is the calculated values * NIfTot+1.
            do i=1,nProcessors
                sendcounts(i)=sendcounts(i)*(NIfTot+1)
                disps(i)=disps(i)*(NIfTot+1)
                sing_recvcounts(i)=sing_recvcounts(i)*(NIfTot+1)
                sing_recvdisps(i)=sing_recvdisps(i)*(NIfTot+1)
            enddo
#ifdef PARALLEL
            CALL MPIAlltoAllv(Sing_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,Sing_ExcDjs2,sing_recvcounts,sing_recvdisps,error)
#else
            Sing_ExcDjs2(0:NIfTot,1:MaxIndex)=Sing_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif

            CALL Sing_SearchOccDets(sing_recvcounts,sing_recvdisps)
        ENDIF


        IF((RDMExcitLevel.eq.2).or.(RDMExcitLevel.eq.3)) THEN
            do i=0,nProcessors-1
                sendcounts(i+1)=Doub_ExcList(i)-(NINT(TwoEl_Gap*i)+1)
! Sendcounts is the number of singly excited determinants we want to send for each processor (but goes from 1, not 0).            
                disps(i+1)=NINT(TwoEl_Gap*i)
! and I think disps is the first slot for each processor - 1.            
            enddo

            MaxSendIndex = Doub_ExcList(nProcessors-1)-1

! We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
            doub_recvcounts(1:nProcessors)=0
            CALL MPIAlltoAll(sendcounts,1,doub_recvcounts,1,error)
! I think recvcounts(i) is the number of determinants sent from processor i.        

! We can now get recvdisps from recvcounts, since we want the data to be contiguous after the move.
            doub_recvdisps(1)=0
            do i=2,nProcessors
                doub_recvdisps(i)=doub_recvdisps(i-1)+doub_recvcounts(i-1)
            enddo

            MaxIndex=doub_recvdisps(nProcessors)+doub_recvcounts(nProcessors)
! But the actual number of integers we need to send is the calculated values * NIfTot+1.
            do i=1,nProcessors
                sendcounts(i)=sendcounts(i)*(NIfTot+1)
                disps(i)=disps(i)*(NIfTot+1)
                doub_recvcounts(i)=doub_recvcounts(i)*(NIfTot+1)
                doub_recvdisps(i)=doub_recvdisps(i)*(NIfTot+1)
            enddo

! This is the main send of all the single excitations to the corresponding processors.        
#ifdef PARALLEL
            CALL MPIAlltoAllv(Doub_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,Doub_ExcDjs2,doub_recvcounts,doub_recvdisps,error)
#else
            Doub_ExcDjs2(0:NIfTot,1:MaxIndex)=Doub_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif
            CALL Doub_SearchOccDets(doub_recvcounts,doub_recvdisps)

        ENDIF

        
    END SUBROUTINE SendProcExcDjs


    SUBROUTINE Sing_SearchOccDets(recvcounts,recvdisps)
! We now have arrays SingExcDjs2 which contain all the single excitations from each processor.
! These number sent from processor i is recvcounts(i), and the first 2 have information about the determinant Di from which
! the Dj's are single excitations (and it's sign).
        USE AnnihilationMod , only : BinSearchParts
        USE FciMCData , only : TotWalkers,CurrentDets
        USE RotateOrbsData , only : SymLabelListInv
        USE bit_reps , only : extract_bit_rep
        INTEGER, INTENT(IN) :: recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: SignDi,SignDj
        INTEGER :: i,j,NoDets,StartDets,PartInd,Indij,Indab 
        INTEGER :: nI(NEl),nJ(NEl),Ex(2,2),FlagsDi,FlagsDj
        LOGICAL :: tDetFound,tParity
        REAL*8 :: ParityFactor

! Take each Dj, and binary search the CurrentDets to see if it is occupied.

        do i=1,nProcessors
! Doing determinants from each processor separately because each has a different D_i it was excited from.

            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            IF(NoDets.gt.1) THEN
                call extract_bit_rep (Sing_ExcDjs2(:,StartDets), nI, SignDi, FlagsDi)

                do j=StartDets+1,(NoDets+StartDets-1)
! D_i is in the first spot - start from the second.                
                
                    iLutnJ(:)=Sing_ExcDjs2(:,j)
! This binary searches CurrentDets between 1 and TotWalkers for determinant iLutnJ.
! If found, tDetFound will be true, and PartInd the index in CurrentDets where the determinant is.
                    CALL BinSearchParts(iLutnJ,1,int(TotWalkers,int32),PartInd,tDetFound)
                    IF(tDetFound) THEN
! Determinant occupied; add c_i*c_j to the relevant element of nElRDM.                    
! Need to first find the orbitals involved in the excitation from D_i -> D_j and the parity.

                        ParityFactor=1.D0
                        Ex(:,:)=0
! Ex(1,1) goes in as the max number of excitations - we know this is an excitation of level RDMExcitLevel.                    
                        Ex(1,1)=1

                        call extract_bit_rep (CurrentDets(:,PartInd), nJ, SignDj, FlagsDj)

! Ex(1,:) comes out as the orbital(s) excited from, Ex(2,:) comes out as the orbital(s) excited to.                    
                        CALL GetExcitation(nI,nJ,NEl,Ex,tParity)

                        IF(Ex(1,1).le.0) CALL Stop_All('Sing_SearchOccDets','nJ is not the correct excitation of nI.')

                        IF(tParity) ParityFactor=-1.D0

                        Indij = SymLabelListInv(Ex(1,1)) 
                        Indab = SymLabelListInv(Ex(2,1)) 
                        
                        OneElRDM( Indij , Indab ) = OneElRDM( Indij , Indab ) + &
                            ( (ParityFactor * (REAL(SignDj(1)) * REAL(SignDi(1))) ) / AllTotPartstemp )

! No normalisation factor just yet - possibly need to revise.                    
                    ENDIF

                enddo
            ENDIF

        enddo
      
    END SUBROUTINE Sing_SearchOccDets


    SUBROUTINE Doub_SearchOccDets(recvcounts,recvdisps)
! We now have arrays SingExcDjs2 which contain all the single excitations from each processor.
! These number sent from processor i is recvcounts(i), and the first 2 have information about the determinant Di from which
! the Dj's are single excitations (and it's sign).
        USE AnnihilationMod , only : BinSearchParts
        USE FciMCData , only : TotWalkers,CurrentDets
        USE RotateOrbsData , only : SymLabelListInv
        USE bit_reps , only : extract_bit_rep
        INTEGER, INTENT(IN) :: recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: SignDi,SignDj
        INTEGER :: i,j,NoDets,StartDets,PartInd,Indij,Indab 
        INTEGER :: nI(NEl),nJ(NEl),Ex(2,2),FlagsDi,FlagsDj
        LOGICAL :: tDetFound,tParity
        REAL*8 :: ParityFactor

! Take each Dj, and binary search the CurrentDets to see if it is occupied.

        do i=1,nProcessors
! Doing determinants from each processor separately because each has a different D_i it was excited from.

            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            IF(NoDets.gt.1) THEN
                call extract_bit_rep (Doub_ExcDjs2(:,StartDets), nI, SignDi, FlagsDi)

                do j=StartDets+1,(NoDets+StartDets-1)
! D_i is in the first spot - start from the second.                
                
                    iLutnJ(:)=Doub_ExcDjs2(:,j)
! This binary searches CurrentDets between 1 and TotWalkers for determinant iLutnJ.
! If found, tDetFound will be true, and PartInd the index in CurrentDets where the determinant is.
                    CALL BinSearchParts(iLutnJ,1,int(TotWalkers,int32),PartInd,tDetFound)
                    IF(tDetFound) THEN
! Determinant occupied; add c_i*c_j to the relevant element of nElRDM.                    
! Need to first find the orbitals involved in the excitation from D_i -> D_j and the parity.

                        ParityFactor=1.D0
                        Ex(:,:)=0
! Ex(1,1) goes in as the max number of excitations - we know this is an excitation of level RDMExcitLevel.                    
                        Ex(1,1)=2

                        call extract_bit_rep (CurrentDets(:,PartInd), nJ, SignDj, FlagsDj)

! Ex(1,:) comes out as the orbital(s) excited from, Ex(2,:) comes out as the orbital(s) excited to.                    
                        CALL GetExcitation(nI,nJ,NEl,Ex,tParity)

                        IF(Ex(1,1).le.0) CALL Stop_All('SearchOccDets','nJ is not the correct excitation of nI.')

                        IF(tParity) ParityFactor=-1.D0

                        !Find unique index for the pairs of orbitals we are exciting from (ij) and to (ab).
!                        Indij = ( ( (Ex(1,2)-2) * (Ex(1,2)-1) ) / 2 ) + Ex(1,1)
!                        Indab = ( ( (Ex(2,2)-2) * (Ex(2,2)-1) ) / 2 ) + Ex(2,1)
                        IF((Ex(1,1).gt.Ex(1,2)).or.(Ex(2,1).gt.Ex(2,2))) THEN
                            WRITE(6,*) 'Ex(1,1)',Ex(1,1),'Ex(1,2)',Ex(1,2),'Ex(2,1)',Ex(2,1),'Ex(2,2)',Ex(2,2)
                            CALL Stop_All('SearchOccDets','The orbitals involved in excitation are not in the expected order.')
                        ENDIF

                        TwoElRDM( Indij , Indab ) = TwoElRDM( Indij , Indab ) + &
                            ( (ParityFactor * (REAL(SignDj(1)) * REAL(SignDi(1))) ) / AllTotPartstemp )

                        !2RDM(p,q,r,s)  = < DI | Epq Ers - delta_qr Eps | DJ>
                        !               = < DI | E(p -> q) E(r -> s) - delta_qr E(p -> s) | DJ >
!                        TwoElRDM( Ex(1,1) , Ex(2,1), Ex(1,2), Ex(2,2) ) = TwoElRDM( Ex(1,1), Ex(2,1), Ex(1,2), Ex(2,2) ) + &
!                            ( (ParityFactor * (REAL(SignDj(1)) * REAL(SignDi(1))) ) / AllTotPartstemp )

! No normalisation factor just yet - possibly need to revise.                    
                    ENDIF

                enddo
            ENDIF

        enddo
      
    END SUBROUTINE Doub_SearchOccDets


    SUBROUTINE FinaliseRDM()
! This routine finalises the one electron reduced density matrix stuff.
! This includes summing each of the individual matrices from each processor.
! Normalisation of some sort?
! Calling the diagonalisation routines if we want to get the occupation numbers.
        USE Logging , only : tDiagRDM,tCalc_RDMEnergy 
        USE NatOrbsMod , only : DiagNatOrbMat,OrderCoeffT1 
        USE NatOrbsMod , only : FillCoeffT1
        USE SystemData , only : tSeparateOccVirt,tRotateVirtOnly,tRotateOccOnly
        USE SystemData , only : tNoRODump, ARR, BRR, G1
        USE RotateOrbsMod , only : FourIndInts, FourIndIntsTag,Transform2ElIntsMemSave
        USE RotateOrbsMod , only : CalcFOCKMatrix,RefillUMATandTMAT2D
        USE RotateOrbsData , only : NoOrbs
        INTEGER :: error,i,j,ierr
        REAL*8 :: SumDiag
        CHARACTER(len=*), PARAMETER :: this_routine='FinaliseRDM'

!        WRITE(6,*) 'nElRDM'
!        do i=1,nBasis
!            do j=1,nBasis
!                WRITE(6,'(F15.1)',advance='no') nElRDM(j,i)
!            enddo
!            WRITE(6,*) ''
!        enddo

        WRITE(6,*) 'here 01'
        CALL FLUSH(6)

        IF((RDMExcitLevel.eq.1).or.(RDMExcitLevel.eq.3)) CALL MPIReduce(OneElRDM,MPI_SUM,NatOrbMat)

        WRITE(6,*) 'here 02'
        CALL FLUSH(6)


        IF((RDMExcitLevel.eq.2).or.(RDMExcitLevel.eq.3)) CALL MPIReduce(TwoElRDM,MPI_SUM,AllTwoElRDM)

        WRITE(6,*) 'here 03'
        CALL FLUSH(6)

!        IF(iProcIndex.eq.0) THEN
!            WRITE(6,*) 'NatOrbMat'
!            do i=1,nBasis
!                do j=1,nBasis
!                    WRITE(6,'(F15.1)',advance='no') NatOrbMat(j,i)
!                enddo
!                WRITE(6,*) ''
!            enddo
!        ENDIF

!        CALL Stop_All('','')

! Getting back into the rotate routines, which use symlabellist2 etc etc, so we need to tell these routines 
! they've been set up as spin orbitals.
! The only time this causes a problem is in UMAT, which will be spatial orbitals.
        tTurnStoreSpinOff=.false.
        IF(.not.tStoreSpinOrbs) THEN
            tTurnStoreSpinOff=.true.
            tStoreSpinOrbs=.true.
        ENDIF

        IF(tCalc_RDMEnergy) CALL Calc_Energy_from_RDM()

        call MPIBarrier(error)

!        CALL Stop_all('','')

        IF(tDiagRDM.and.((RDMExcitLevel.eq.1).or.(RDMExcitLevel.eq.3))) THEN
! Call the routines from NatOrbs that diagonalise the one electron reduced density matrix.

            IF(iProcIndex.eq.0) THEN
        
!                SumDiag=0.D0
!                do i=1,nBasis
!                    SumDiag=SumDiag+ABS(NatOrbMat(i,i))
!                enddo
!                do i=1,nBasis
!                    do j=1,nBasis
!                        NatOrbMat(j,i)=NatOrbMat(j,i)/SumDiag
!                    enddo
!                enddo

                tRotateVirtOnly=.true.
                tRotateOccOnly=.false.
                tSeparateOccVirt=.false.
                CALL DiagNatOrbMat()

                CALL OrderCoeffT1()

!                IF(tTurnStoreSpinOff) tStoreSpinOrbs=.false.

                SumDiag=0.D0
                do i=1,nBasis
                    SumDiag=SumDiag+Evalues(i)
                enddo

                do i=1,nBasis
                    Evalues(i)=Evalues(i)/(SumDiag/REAL(NEl))
                enddo
                WRITE(6,*) 'Normalised Evalues:'
                do i=1,nBasis
                    WRITE(6,*) Evalues(i)
                enddo

                IF(.not.tNoRODump) THEN

                    ALLOCATE(CoeffT1(NoOrbs,NoOrbs),stat=ierr)
                    CALL LogMemAlloc(this_routine,NoOrbs*NoOrbs,8,this_routine,CoeffT1Tag,ierr)
                    CoeffT1(:,:)=0.D0

                    CALL FillCoeffT1()

                    ALLOCATE(FourIndInts(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
                    CALL LogMemAlloc('FourIndInts',(NoOrbs**4),8,this_routine,FourIndIntsTag,ierr)

! Then, transform2ElInts
                    WRITE(6,*) 'Transforming the four index integrals'
                    CALL Transform2ElIntsMemSave()

                    WRITE(6,*) 'Re-calculating the fock matrix'
                    CALL CalcFOCKMatrix()

                    WRITE(6,*) 'Refilling the UMAT and TMAT2D'
! The ROFCIDUMP is also printed out in here.        
                    CALL RefillUMATandTMAT2D()        

                    CALL FLUSH(6)

                    CALL WRITEBASIS(6,G1,nBasis,ARR,BRR)

                ENDIF
            ENDIF

        ELSEIF(tDiagRDM) THEN
        
            WRITE(6,*) 'WARNING: Request to diagonalise two electron reduced density matrix.'

        ENDIF

! Print out some stuff about the one electron reduced density matrix.

! This is where we would likely call any further calculations of force etc.


    
    END SUBROUTINE FinaliseRDM



    SUBROUTINE DeallocateRDM()
        USE RotateOrbsMod , only : SymLabelList2,SymLabelListInv,SymLabelListInvTag,SymLabelList2Tag
        USE RotateOrbsMod , only : FourIndInts, FourIndIntsTag, SymLabelCounts2, SymLabelCounts2Tag
        USE RotateOrbsMod , only : SymLabelList3, SymLabelList3Tag
        USE SystemData , only : tNoRODump
        USE Logging , only : tDiagRDM
        CHARACTER(len=*), PARAMETER :: this_routine='DeallocateRDM'

        IF((RDMExcitLevel.eq.1).or.(RDMExcitLevel.eq.3)) THEN
! This array contains the initial positions of the single excitations for each processor.
            DEALLOCATE(Sing_InitExcSlots)
 
! This array contains the current position of the single excitations as they're added.
            DEALLOCATE(Sing_ExcList)

! This array actually contains the single excitations in blocks of the processor they will be sent to.        
            DEALLOCATE(Sing_ExcDjs)
            CALL LogMemDeAlloc(this_routine,Sing_ExcDjsTag)
 
            DEALLOCATE(Sing_ExcDjs2)
            CALL LogMemDeAlloc(this_routine,Sing_ExcDjs2Tag)

! We also need to allocate the actual 1RDM on each processor, and an all1RDM on only the root.        
            DEALLOCATE(OneElRDM)
            CALL LogMemDeAlloc(this_routine,OneElRDMTag)

            IF(iProcIndex.eq.0) THEN
                DEALLOCATE(Evalues)
                CALL LogMemDeAlloc(this_routine,EvaluesTag)

                DEALLOCATE(NatOrbMat)
                CALL LogMemDeAlloc(this_routine,NatOrbMatTag)
            ENDIF

            DEALLOCATE(SymLabelCounts2)
            CALL LogMemDeAlloc(this_routine,SymLabelCounts2Tag)
 
            DEALLOCATE(SymLabelList2)
            CALL LogMemDeAlloc(this_routine,SymLabelList2Tag)
 
            DEALLOCATE(SymLabelList3)
            CALL LogMemDeAlloc(this_routine,SymLabelList3Tag)
 
            DEALLOCATE(SymLabelListInv)
            CALL LogMemDeAlloc(this_routine,SymLabelListInvTag)
 
            IF(tDiagRDM.and.(.not.tNoRODump)) THEN
                DEALLOCATE(CoeffT1)
                CALL LogMemDeAlloc(this_routine,CoeffT1Tag)
     
                DEALLOCATE(FourIndInts)
                CALL LogMemDeAlloc(this_routine,FourIndIntsTag)
     
            ENDIF
        ENDIF

        IF((RDMExcitLevel.eq.2).or.(RDMExcitLevel.eq.3)) THEN
! This array contains the initial positions of the single excitations for each processor.
            DEALLOCATE(Doub_InitExcSlots)
 
! This array contains the current position of the single excitations as they're added.
            DEALLOCATE(Doub_ExcList)

! This array actually contains the single excitations in blocks of the processor they will be sent to.        
            DEALLOCATE(Doub_ExcDjs)
            CALL LogMemDeAlloc(this_routine,Doub_ExcDjsTag)
 
            DEALLOCATE(Doub_ExcDjs2)
            CALL LogMemDeAlloc(this_routine,Doub_ExcDjs2Tag)

! We also need to allocate the actual 1RDM on each processor, and an all1RDM on only the root.        
            WRITE(6,*) 'Having issues deallocating TwoElRDM (possibly because of array bounds issues)'
            CALL FLUSH(6)

            DEALLOCATE(TwoElRDM)
            CALL LogMemDeAlloc(this_routine,TwoElRDMTag)

            WRITE(6,*) 'Nope, deallocation of TwoElRDM is o.k'
            CALL FLUSH(6)

            IF(iProcIndex.eq.0) THEN
                DEALLOCATE(AllTwoElRDM)
                CALL LogMemDeAlloc(this_routine,AllTwoElRDMTag)
            ENDIF

        ENDIF

    END SUBROUTINE DeallocateRDM


    SUBROUTINE Calc_Energy_from_RDM()
!This routine takes the 1 electron and 2 electron reduced density matrices and calculated the energy they give.    
!The equation for the energy is as follows:
!
!   E = Tr(h1 1RDM) + 1/2 Tr(h2 2RDM)
!
!where h1 are the 2 index integrals, and h2 the 4 index integrals.  The traces, Tr, are given by:
!   Tr(h1 1RDM) = Sum_i,j [ h1(i,j) 1RDM(j,i) ]
!   Tr(h2 2RDM) = Sum_i,j;k,l [ h2(i,j;k,l) 2RDM(k,l;i,j) ]
        USE IntegralsData , only : UMAT
        USE UMatCache , only : UMatInd
        USE OneEInts , only : TMAT2D
        USE RotateOrbsMod , only : SymLabelList2
        INTEGER :: i,j,k,l,Ind2,Ind1,i2,j2,k2,l2
        REAL*8 :: RDMEnergy, Trace_1RDM

!First of all 'normalise' the reduced density matrices.
!These are not even close to being normalised at the moment, because of the way they are calculated on the fly.
!They should be calculated from a normalised wavefunction.
!But we know that the trace of the one electron reduced density matrix must be equal to the number of the electrons.
!We can use this to find the factor we must divide the 1RDM through by and apply this to the 2RDM as well.

        WRITE(6,*) 'inside this routine'
        CALL FLUSH(6)

        Trace_1RDM = 0.D0

        RDMEnergy = 0.D0

        IF(iProcIndex.eq.0) THEN
            do i = 1, nBasis
                Trace_1RDM = Trace_1RDM + NatOrbMat(i,i)
            enddo

            WRITE(6,*) 'trace o.k'
            WRITE(6,*) 'Trace_1RDM',Trace_1RDM

            WRITE(6,*) '((nBasis*(nBasis-1))/2)',((nBasis*(nBasis-1))/2)
            WRITE(6,*) 'tTurnStoreSpinOff',tTurnStoreSpinOff
            CALL FLUSH(6)

            !Need to multiply each element of the reduced density matrices by NEl / Trace_1RDM,
            !and then add it's contribution to the energy.
            do i = 1, nBasis 

                IF(tTurnStoreSpinOff) THEN
                    i2 = CEILING(REAL(i)/2.D0)
                ELSE
                    i2 = i
                ENDIF
                do j = i, nBasis

                    NatOrbMat(SymLabelList2(j),SymLabelList2(i)) = &
                            NatOrbMat(SymLabelList2(j),SymLabelList2(i)) * ( REAL(NEl,8) / Trace_1RDM )  
                    
    !                WRITE(6,*) 'NatOrbMat(i,j)',NatOrbMat(j,i)
    !                CALL FLUSH(6)

                    !TMAT2D always stored in spin orbitals so o.k. 
                    RDMEnergy = RDMEnergy + ( REAL(TMAT2D(i,j),8)  &
                                                * NatOrbMat(SymLabelList2(j),SymLabelList2(i)) )

    !                WRITE(6,*) 'RDMEnergy 01',RDMEnergy

                    IF(i.ne.j) THEN
                        NatOrbMat(SymLabelList2(i),SymLabelList2(j)) = &
                                NatOrbMat(SymLabelList2(i),SymLabelList2(j)) * ( REAL(NEl,8) / Trace_1RDM )  

                        RDMEnergy = RDMEnergy + ( REAL(TMAT2D(j,i),8) * &
                                                        NatOrbMat(SymLabelList2(i),SymLabelList2(j)) )
                    ENDIF

    !                WRITE(6,*) 'RDMEnergy 02',RDMEnergy

                    Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i

                    IF(tTurnStoreSpinOff) THEN
                        j2 = CEILING(REAL(j)/2.D0)
                    ELSE
                        j2 = j
                    ENDIF

                    do k = 1, nBasis

                        IF(tTurnStoreSpinOff) THEN
                            k2 = CEILING(REAL(k)/2.D0)
                        ELSE
                            k2 = k
                        ENDIF
                        do l = k, nBasis

                            IF(tTurnStoreSpinOff) THEN
                                l2 = CEILING(REAL(l)/2.D0)
                            ELSE
                                l2 = l
                            ENDIF

                            Ind2 = ( ( (l-2) * (l-1) ) / 2 ) + k

    !                        WRITE(6,*) 'Ind1',Ind1,'Ind2',Ind2

                            AllTwoElRDM(Ind2,Ind1) = AllTwoElRDM(Ind2,Ind1) * ( REAL(NEl,8) / Trace_1RDM )  

    !                        WRITE(6,*) 'RDMEnergy 03',RDMEnergy
    !                        WRITE(6,*) 'i',i,'j',j,'k',k,'l',l

                                ! AllTwoElRDM etc will be in spin orbitals, but UMAT is in spatial, need to correct for this.

                                !UmatInd uses chemical notation.
                            RDMEnergy = RDMEnergy + ( 0.50 * REAL(UMAT(UMatInd(i2,k2,j2,l2,0,0)),8) * AllTwoElRDM(Ind2,Ind1) )

                            IF(i2.ne.j2) THEN

                                RDMEnergy = RDMEnergy + ( 0.50 * REAL(UMAT(UMatInd(j2,k2,i2,l2,0,0)),8) * AllTwoElRDM(Ind2,Ind1) )

                            ENDIF

                            IF(k2.ne.l2) THEN

                                RDMEnergy = RDMEnergy + ( 0.50 * REAL(UMAT(UMatInd(i2,l2,j2,k2,0,0)),8) * AllTwoElRDM(Ind2,Ind1) )

                            ENDIF

                            IF((i2.ne.j2).and.(k2.ne.l2)) THEN

                                RDMEnergy = RDMEnergy + ( 0.50 * REAL(UMAT(UMatInd(j2,l2,i2,k2,0,0)),8) * AllTwoElRDM(Ind2,Ind1) )

                            ENDIF

    !                        WRITE(6,*) 'RDMEnergy 04',RDMEnergy

                        enddo
                    enddo
                enddo
            enddo

            WRITE(6,*) '******'
            WRITE(6,*) 'The ENERGY calculated using the REDUCED DENSITY MATRICES is:',RDMEnergy

        ENDIF

    END SUBROUTINE Calc_Energy_from_RDM


! =============================== OLD ROUTINES =============================================================================    

    SUBROUTINE NormandDiagRDM()
!This routine takes the OneRDM with non-normalised off-diagonal elements, adds the diagonal elements from the HF determinant
!and normalises it according to the number of walkers on the HF determinant.
!It then diagonalises the 1-RDM to find linear combinations of the HF orbitals that are closer to the natural orbitals,
!and the occupation numbers of these new orbitals (e-values).
        USE NatOrbsMod , only : FindNatOrbsOld
        IMPLICIT NONE
        INTEGER :: i,j,error,ierr
        REAL*8 :: TempSumNoatHF
        REAL*8 , ALLOCATABLE :: TempRDM(:,:)

        INTEGER :: SumNoatHF,AllSumNoatHF,RDM(100,100),HFDet(10)

        TempSumNoatHF=real(SumNoatHF)
!Sum TempSumNoatHF over all processors and then send to all processes
        CALL MPI_AllReduce(TempSumNoatHF,AllSumNoatHF,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
        ALLOCATE(TempRDM(nBasis,nBasis),stat=ierr)
        IF(ierr.ne.0) THEN
            CALL Stop_All("NormandDiagRDM","Could not allocate TempRDM")
        ENDIF

        CALL MPI_AllReduce(RDM(:,:),TempRDM(:,:),nBasis*nBasis,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
        RDM(:,:)=TempRDM(:,:)
        DEALLOCATE(TempRDM)
        
!Normalise the off-diag OneRDM elements using the number of walkers on the HF determinant
        do i=1,nBasis
            do j=i+1,nBasis
                RDM(i,j)=RDM(i,j)/AllSumNoatHF
                RDM(j,i)=RDM(i,j)
            enddo
        enddo
!Using the decoded version of the HF determinant, place values of 1.D0 in the diagonal elements
!of the 1-RDM corresponding to the occupied orbitals.
        do i=1,NEl
            RDM(HFDet(i),HFDet(i))=1.D0
        enddo
    
        CALL FindNatOrbsOld()           !Diagonalise the 1-RDM

    END SUBROUTINE NormandDiagRDM



!          IF(tConstructNOs) THEN

!!Fill the 1-RDM to find natural orbital on-the-fly.
!!Find the orbitals that are involved in the excitation (those that differ in occupation to the ref orbital).
!                CALL FindSingleOrbs(iLutHF,iLutCurr,NIfD,Orbs)
!!Add 1.D0 (or -1.D0) to the off-diagonal element connecting the relevent orbitals.
!                IF(Iter.gt.NEquilSteps) THEN
!                    OneRDM(Orbs(1),Orbs(2))=OneRDM(Orbs(1),Orbs(2))+REAL(WSign,dp)
!                    OneRDM(Orbs(2),Orbs(1))=OneRDM(Orbs(1),Orbs(2))
!                ENDIF
!!At the end of all iterations, this OneRDM will contain only the unnormalised off-diagonal elements.
!          ENDIF
 
END MODULE nElRDMMod


