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
        USE IntegralsData , only : UMAT
        USE UMatCache , only : UMatInd
        USE SystemData , only : NEl,nBasis,tStoreSpinOrbs, G1, BRR, lNoSymmetry, ARR
        USE SystemData , only : tUseMP2VarDenMat, Ecore, LMS 
        USE NatOrbsMod , only : NatOrbMat,NatOrbMatTag,Evalues,EvaluesTag
        USE CalcData , only : MemoryFacPart
        USE constants , only : n_int
        USE OneEInts , only : TMAT2D
        USE Logging , only : RDMExcitLevel, tROFciDUmp, NoDumpTruncs
        USE RotateOrbsData , only : CoeffT1, CoeffT1Tag, tTurnStoreSpinOff, NoFrozenVirt
        USE RotateOrbsData , only : SymLabelCounts2,SymLabelList2,SymLabelListInv,NoOrbs
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
        USE RotateOrbsMod , only : SymLabelCounts2Tag,SpatOrbs,NoRotOrbs
        USE RotateOrbsMod , only : SymLabelList2Tag,SymLabelListInvTag
        USE RotateOrbsData , only : SymLabelList3, SymLabelList3Tag
        USE Logging , only : tCalc_RDMEnergy, tDiagRDM
        INTEGER :: ierr,i
        CHARACTER(len=*), PARAMETER :: this_routine='InitRDM'

        IF(tCalc_RDMEnergy) THEN
            WRITE(6,*) 'Calculating the energy from the reduced density matrix, this requires both the 1 and 2 electron RDM.'
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
        NoDumpTruncs = 1
        NoFrozenVirt = 0

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

            CALL SetUpSymLabels_RDM() 

        ENDIF            

        !Need to now turn tStoreSpinOrbs back off for the rest of the calculation, but all the 
        !symmetry arrays etc are set up for tStoreSpinOrbs to be on.
        !Note, if this has been meddled with, the UMAT will be in spatial orbitals and everything else in spin.
        !i.e. everything is in spin orbitals, except for UMAT which will be in spatial if tStoreSpinOrbs is false.
        IF(tTurnStoreSpinOff) tStoreSpinOrbs=.false.

        nElRDM_Time%timer_name='nElRDMTime'


    END SUBROUTINE InitRDM


    SUBROUTINE SetUpSymLabels_RDM() 
!We always want the RDM's to be stored in spin orbitals (for now), so this routine just sets up the symmetry labels so that
!the orbitals are ordered according to symmetry (all beta then all alpha).
!We also always want to mix both the occupied and virtual orbitals, so there doesn't need to be any separation of these.
        USE sort_mod
        IMPLICIT NONE
        INTEGER , ALLOCATABLE :: LabOrbs(:), SymOrbs(:)
        INTEGER :: LabOrbsTag, SymOrbsTag, ierr, spin, i , j 
        INTEGER :: StartFill, Prev, lo, hi, Symi, SymCurr
        CHARACTER(len=*) , PARAMETER :: this_routine = 'SetUpSymLabels_RDM'

        ALLOCATE(LabOrbs(nBasis/2),stat=ierr)
        CALL LogMemAlloc('LabOrbs',nBasis/2,4,this_routine,LabOrbsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for LabOrbs failed.")
        ALLOCATE(SymOrbs(nBasis/2),stat=ierr)
        CALL LogMemAlloc('SymOrbs',nBasis/2,4,this_routine,SymOrbsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for SymOrbs failed.")

! Brr has the orbital numbers in order of energy... i.e Brr(2) = the orbital index with the second lowest energy.
!            do i=1,nBasis
!                WRITE(6,*) BRR(i)
!            enddo
!            CALL FLUSH(6)
!            CALL Stop_All('','')

! We separate the alphas from the betas (because there'll never be mixing between the two - so all diagonalisation etc 
! can be done separately).  e.g SymLabelList has all beta orbitals followed by all alpha orbitals.
        do spin = 1,2       !do beta then alpha

            LabOrbs(:)=0
            SymOrbs(:)=0

! *** STEP 1 *** Fill SymLabelList2.
! we first find all the beta orbitals - and order them in terms of symmetry, and then the alpha.
            do i=1,nBasis/2
                IF(spin.eq.1) THEN
                    LabOrbs(i)=BRR((2*i))-1
                    SymOrbs(i)=INT(G1(LabOrbs(i))%sym%S,4)
                ELSEIF(spin.eq.2) THEN
                    LabOrbs(i)=BRR((2*i))
                    SymOrbs(i)=INT(G1(LabOrbs(i))%sym%S,4)
                ENDIF
            enddo

!            WRITE(6,*) 'before sort'
!            do i=1,nBasis/2
!                WRITE(6,*) LabOrbs(i),SymOrbs(i)
!            enddo

            call sort (SymOrbs, LabOrbs)
            ! Sorts LabOrbs according to the order of SymOrbs (i.e. in terms of symmetry). 

!            WRITE(6,*) 'after sort'
!            do i=1,nBasis/2
!                WRITE(6,*) LabOrbs(i),SymOrbs(i)
!            enddo
!            CALL FLUSH(6)
!            stop

 
! SymLabelList2 is then filled with the symmetry ordered beta then alpha arrays.        
            IF(spin.eq.1) THEN
                j = 0
                do i = 1, nBasis/2
                    j = j + 1
                    SymLabelList2(i) = LabOrbs(j)
                enddo
            ELSEIF(spin.eq.2) THEN
                j = 0
                do i = (nBasis/2) + 1 , nBasis
                    j = j + 1
                    SymLabelList2(i) = LabOrbs(j)
                enddo
            ENDIF


!*** STEP 2 *** Fill SymLabelCounts2.
!This again has all beta indices followed by all alpha.
!SymLabelCounts(1,:) contains the position in SymLabelList2 where the symmetry index starts, and
!SymLabelCounts(2,:) contains the number of orbitals in that symmetry index.

            IF(lNoSymmetry) THEN
                ! if we are ignoring symmetry, all orbitals essentially have symmetry 0.
                IF(spin.eq.1) THEN
                    SymLabelCounts2(1,1)=1
                    SymLabelCounts2(2,1)=nBasis/2
                ELSEIF(spin.eq.2) THEN
                    SymLabelCounts2(1,9)=1
                    SymLabelCounts2(2,9)=nBasis/2
                ENDIF
                
            ELSE 
                ! otherwise we run through the occupied orbitals, counting the number with each symmetry
                ! and noting where in SymLabelList2 each symmetry block starts.
                IF(spin.eq.1) THEN
                    StartFill=1
                    Prev=0
                ELSEIF(spin.eq.2) THEN
                    StartFill=9
                    Prev=nBasis/2
                ENDIF
                SymCurr=0
                SymLabelCounts2(1,StartFill)=1+Prev
                do i=1,nBasis/2
                    Symi=INT(G1(SymLabelList2(i+Prev))%sym%S,4)
                    SymLabelCounts2(2,(Symi+StartFill))=SymLabelCounts2(2,(Symi+StartFill))+1
                    IF(Symi.ne.SymCurr) THEN
                        SymLabelCounts2(1,(Symi+StartFill))=i+Prev
                        SymCurr=Symi
                    ENDIF
                enddo
            ENDIF
     
            ! Go through each symmetry group, making sure the orbital pairs are ordered lowest to highest.
            IF(spin.eq.1) THEN
                do i=1,8
                    IF(SymLabelCounts2(2,i).ne.0) THEN
                        lo = SymLabelCounts2(1, i)
                        hi = lo + SymLabelCounts2(2, i) - 1
                        call sort (SymLabelList2 (lo:hi))
                    ENDIF
                enddo
            ELSEIF(spin.eq.2) THEN
                do i=9,16
                    IF(SymLabelCounts2(2,i).ne.0) THEN
                        lo = SymLabelCounts2(1, i)
                        hi = lo + SymLabelCounts2(2, i) - 1
                        call sort (SymLabelList2 (lo:hi))
                    ENDIF
                enddo
            ENDIF

        enddo

        do i=1,nBasis
            SymLabelListInv(SymLabelList2(i))=i
        enddo

! Deallocate the arrays just used in this routine.

        DEALLOCATE(SymOrbs)
        CALL LogMemDealloc(this_routine,SymOrbsTag)

        DEALLOCATE(LabOrbs)
        CALL LogMemDealloc(this_routine,LabOrbsTag)


!        WRITE(6,*) 'Sym Label Counts'
!        do i=1,16
!            WRITE(6,*) i,SymLabelCounts2(1,i),SymLabelCounts2(2,i)
!        enddo
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and their symmetries according to G1'
!        do i=1,nBasis
!            WRITE(6,*) i,SymLabelList2(i),INT(G1(SymLabelList2(i))%sym%S,4)
!        enddo
!        WRITE(6,*) 'i','ARR(SymLabelList2(i),1)','ARR(SymLabelList2(i),2)','Sym'
!        do i=1,NoOrbs
!            IF(tStoreSpinOrbs) THEN
!                WRITE(6,*) i,ARR(SymLabelList2(i),1),ARR(SymLabelList2(i),2),INT(G1(SymLabelList2(i))%sym%S,4)
!            ENDIF
!        enddo
!
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and its inverse'
!        do i=1,NoOrbs
!            WRITE(6,*) SymLabelList2(i),SymLabelListInv(i)
!        enddo
!        CALL FLUSH(6)
!        CALL Stop_All('SetUpSymLabels_RDM','Checking orbital labelling.')


    END SUBROUTINE SetUpSymLabels_RDM


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

! There is no need to use the SymLabelList arrays for the 2 el RDM because we are not diagonalising 
! or anything.
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

!            WRITE(6,*) 'Generating from nI',nI
!            WRITE(6,*) 'bit rep',iLutnI

            do while (.not.tAllExcitFound)
                CALL GenExcitations3(nI,iLutnI,nJ,2,ExcitMat3(:,:),tParity,tAllExcitFound)            
! Passed out of here is the doubly excited determinant, nJ.
! Information such as the orbitals involved in the excitation and the parity is also found in this step,
! we are not currently storing this, and it is re-calculated later on (after the determinants are passed
! to the relevant processor) - but the speed of sending this information vs recalculating it will be tested.
! RDMExcitLevel is passed through, if this is 1, only singles are generated, if it is 2 only doubles are found.

                IF(tAllExcitFound) EXIT

                iLutnJ(:)=0
                CALL EncodeBitDet(nJ,iLutnJ)

                Proc = DetermineDetNode(nJ,0)   !This will return a value between 0 -> nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(Proc)) = iLutnJ(:)
                ! All the double excitations from this particular nI are being stored in Doub_ExcDjs.

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
        INTEGER :: i,j,k,NoDets,StartDets,PartInd,Indij,Indab 
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

                            
                        IF(RDMExcitLevel.eq.3) THEN
                            do k = 1, NEl                            
                                IF(nI(k).ne.Ex(1,1)) THEN
                                    Indij = ( ( (MAX(nI(k),Ex(1,1))-2) * (MAX(nI(k),Ex(1,1))-1) ) / 2 ) + MIN(nI(k),Ex(1,1))
                                    Indab = ( ( (MAX(nI(k),Ex(2,1))-2) * (MAX(nI(k),Ex(2,1))-1) ) / 2 ) + MIN(nI(k),Ex(2,1))

                                    TwoElRDM( Indij , Indab ) = TwoElRDM( Indij , Indab ) + &
                                        ( (ParityFactor * (REAL(SignDj(1)) * REAL(SignDi(1))) ) / AllTotPartstemp )
                                ENDIF
                            enddo
                        ENDIF

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
    
!        WRITE(6,*) 'Searching for generated nJs in occupied list'

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

!                        WRITE(6,*) 'DetFound'
!                        WRITE(6,*) 'ilutnJ',iLutnJ
!                        WRITE(6,*) 'CurrentDets(:,PartInd)',CurrentDets(:,PartInd)
!                        WRITE(6,*) 'nJ',nJ
!                        WRITE(6,*) 'From',Ex(1,1),Ex(1,2)
!                        WRITE(6,*) 'To',Ex(2,1),Ex(2,2)
!                        WRITE(6,*) 'Adding in ',SignDj(1),'times',SignDi(1)

                        !Find unique index for the pairs of orbitals we are exciting from (ij) and to (ab).
                        Indij = ( ( (Ex(1,2)-2) * (Ex(1,2)-1) ) / 2 ) + Ex(1,1)
                        Indab = ( ( (Ex(2,2)-2) * (Ex(2,2)-1) ) / 2 ) + Ex(2,1)
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
        USE SystemData , only : tSeparateOccVirt,tRotateVirtOnly,tRotateOccOnly
        USE SystemData , only : tNoRODump, ARR, BRR, G1
        USE RotateOrbsMod , only : FourIndInts, FourIndIntsTag
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

        IF((RDMExcitLevel.eq.1).or.(RDMExcitLevel.eq.3)) CALL MPIReduce(OneElRDM,MPI_SUM,NatOrbMat)

        IF((RDMExcitLevel.eq.2).or.(RDMExcitLevel.eq.3)) CALL MPIReduce(TwoElRDM,MPI_SUM,AllTwoElRDM)

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

                CALL DiagRDM()

                tRotateVirtOnly=.true.
                tRotateOccOnly=.false.
                tSeparateOccVirt=.false.
                WRITE(6,*) 'NatOrbMat'
                do i = 1, nBasis
                    do j = 1, nBasis
                        WRITE(6,'(F20.10)',advance='no') NatOrbMat(j,i)
                    enddo
                    WRITE(6,*) ''
                enddo

                CALL OrderRDM()

                IF(tTurnStoreSpinOff) tStoreSpinOrbs=.false.

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

 !                   WRITE(6,*) 'NoOrbs',NoOrbs

                    ALLOCATE(CoeffT1(NoOrbs,NoOrbs),stat=ierr)
                    CALL LogMemAlloc(this_routine,NoOrbs*NoOrbs,8,this_routine,CoeffT1Tag,ierr)
                    CoeffT1(:,:)=0.D0

                    CALL FillCoeffT1_RDM()

!                    WRITE(6,*) 'CoeffT1'
!                    do i = 1, nBasis
!                        do j = 1, nBasis
!                            WRITE(6,'(F10.6)',advance='no') CoeffT1(j,i)
!                        enddo
!                        WRITE(6,*) ''
!                    enddo
!                    WRITE(6,*) 'SymLabelList2'
!                    do i = 1, nBasis
!                        WRITE(6,*) SymLabelList2(i),INT(G1(SymLabelList2(i))%sym%S,4) 
!                    enddo

                    ALLOCATE(FourIndInts(nBasis,nBasis,nBasis,nBasis),stat=ierr)
                    CALL LogMemAlloc('FourIndInts',(nBasis**4),8,this_routine,FourIndIntsTag,ierr)

! Then, transform2ElInts
                    WRITE(6,*) 'Transforming the four index integrals'
                    CALL Transform2ElIntsMemSave_RDM()

                    WRITE(6,*) 'Re-calculating the fock matrix'
                    CALL CalcFOCKMatrix_RDM()

                    WRITE(6,*) 'Refilling the UMAT and TMAT2D'
! The ROFCIDUMP is also printed out in here.        
                    CALL RefillUMATandTMAT2D_RDM()        

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


    SUBROUTINE DiagRDM()
! The diagonalisation routine reorders the orbitals in such a way that the corresponding orbital labels are lost.
! In order to keep the spin and spatial symmetries, each symmetry must be fed into the diagonalisation routine separately.
! The best way to do this is to order the orbitals so that all the alpha orbitals follow all the beta orbitals, with the 
! occupied orbitals first, in terms of symmetry, and the virtual second, also ordered by symmetry.
! This gives us flexibility w.r.t rotating only the occupied or only virtual and looking at high spin states.
        IMPLICIT NONE
        REAL*8 :: SumTrace,SumDiagTrace
        REAL*8 , ALLOCATABLE :: WORK2(:),EvaluesSym(:),NOMSym(:,:)
        INTEGER :: ierr,i,j,spin,Sym,LWORK2,WORK2Tag,SymStartInd,NoSymBlock,PrevSym
        INTEGER :: EvaluesSymTag,NOMSymTag
        CHARACTER(len=*), PARAMETER :: this_routine='DiagRDM'

 
!        DiagNatOrbMat_Time%timer_name='DiagNatOrb'
!        CALL set_timer(DiagNatOrbMat_Time,30)

! Test that we're not breaking symmetry.
        do i=1,nBasis
            do j=1,nBasis
!                WRITE(6,*) INT(G1(SymLabelList2(i))%sym%S,4),INT(G1(SymLabelList2(j))%sym%S,4),NatOrbMat(i,j)
                IF((INT(G1(SymLabelList2(i))%sym%S,4).ne.INT(G1(SymLabelList2(j))%sym%S,4))) THEN
                    IF(ABS(NatOrbMat(i,j)).ge.1.0E-15) THEN
                        WRITE(6,'(6A8,A20)') 'i','j','Label i','Label j','Sym i','Sym j','Matrix value'
                        WRITE(6,'(6I3,F40.20)') i,j,SymLabelList2(i),SymLabelList2(j),INT(G1(SymLabelList2(i))%sym%S,4),INT(G1(SymLabelList2(j))%sym%S,4),NatOrbMat(i,j)
                        IF(tUseMP2VarDenMat) THEN
                            WRITE(6,*) '**WARNING** - There is a non-zero NatOrbMat value between orbitals of different symmetry.'
                            WRITE(6,*) 'These elements will be ignored, and the symmetry maintained in the final transformation matrix.'
                        ELSE
                            CALL Stop_All(this_routine,'Non-zero NatOrbMat value between different symmetries.')
                        ENDIF
                    ENDIF
                    NatOrbMat(i,j)=0.D0
                ENDIF
            enddo
        enddo

        SumTrace=0.D0
        do i=1,nBasis
            SumTrace=SumTrace+NatOrbMat(i,i)
        enddo

        WRITE(6,*) 'Calculating eigenvectors and eigenvalues of NatOrbMat'
        CALL FLUSH(6)

        ! We need to feed in the alpha and beta spins separately.
        ! Otherwise these jumble up and the final ordering is uncorrect. 
        ! There should be no non-zero values between these, but can put a check in for this.

        do spin=1,2

! If we want to maintain the symmetry, we cannot have all the orbitals jumbled up when the diagonaliser reorders the eigenvectors.
! Must instead feed each symmetry block in separately.
! This means that although the transformed orbitals are jumbled within the symmetry blocks, the symmetry labels are all that are relevant and these are unaffected.
            IF(spin.eq.1) THEN
                PrevSym=1
            ELSEIF(spin.eq.2) THEN
                PrevSym=9
            ENDIF

            Sym=0
            LWORK2=-1
            do while (Sym.le.7)

                NoSymBlock=SymLabelCounts2(2,Sym+PrevSym)

                SymStartInd=SymLabelCounts2(1,Sym+PrevSym)-1
                ! This is one less than the index that the symmetry starts, so that when we run through i=1,..., we can
                ! start at SymStartInd+i.

                IF(NoSymBlock.gt.1) THEN

                    ALLOCATE(NOMSym(NoSymBlock,NoSymBlock),stat=ierr)
                    CALL LogMemAlloc('NOMSym',NoSymBlock**2,8,this_routine,NOMSymTag,ierr)
                    IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating NOMSym.")
                    ALLOCATE(EvaluesSym(NoSymBlock),stat=ierr)
                    CALL LogMemAlloc('EvaluesSym',NoSymBlock,8,this_routine,EvaluesSymTag,ierr)
                    IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating EvaluesSym.")

                    LWORK2=3*NoSymBlock+1
                    ALLOCATE(WORK2(LWORK2),stat=ierr)
                    CALL LogMemAlloc('WORK2',LWORK2,8,this_routine,WORK2Tag,ierr)
                    IF(ierr.ne.0) CALL Stop_All(this_routine,"Problem allocating WORK2.")

                    do j=1,NoSymBlock
                        do i=1,NoSymBlock
                            NOMSym(i,j)=NatOrbMat(SymStartInd+i,SymStartInd+j)
                        enddo
                    enddo


                    WRITE(6,*) '*****'
                    WRITE(6,*) 'Symmetry ',Sym, 'with spin ',spin,' has ',NoSymBlock,' orbitals.'
                    WRITE(6,*) 'The NatOrbMat for this symmetry block is '
                    do i=1,NoSymBlock
                        do j=1,NoSymBlock
                            WRITE(6,'(F20.10)',advance='no') NOMSym(j,i)
                        enddo
                        WRITE(6,*) ''
                    enddo

                    CALL DSYEV('V','L',NoSymBlock,NOMSym,NoSymBlock,EvaluesSym,WORK2,LWORK2,ierr)
                    ! NOMSym goes in as the original NOMSym, comes out as the eigenvectors (Coefficients).
                    ! EvaluesSym comes out as the eigenvalues in ascending order.

                    WRITE(6,*) 'After diagonalization, the e-vectors (diagonal elements) of this matrix are ,'
                    do i=1,NoSymBlock
                        WRITE(6,'(F20.10)',advance='no') EvaluesSym(i)
                    enddo
                    WRITE(6,*) ''
                    WRITE(6,*) 'These go from orbital ,',SymStartInd+1,' to ',SymStartInd+NoSymBlock

                    do i=1,NoSymBlock
                        Evalues(SymStartInd+i)=EvaluesSym(i)
                    enddo

                    ! CAREFUL if eigenvalues are put in ascending order, this may not be correct, with the labelling system.
                    ! may be better to just take coefficients and transform TMAT2DRot in transform2elints.
                    ! a check that comes out as diagonal is a check of this routine anyway.

                    WRITE(6,*) 'The eigenvectors (coefficients) for symmetry block ',Sym
                    do i=1,NoSymBlock
                        do j=1,NoSymBlock
                            WRITE(6,'(F20.10)',advance='no') NOMSym(j,i)
                        enddo
                        WRITE(6,*) ''
                    enddo
             
                    do j=1,NoSymBlock
                        do i=1,NoSymBlock
                            NatOrbMat(SymStartInd+i,SymStartInd+j)=NOMSym(i,j)
                        enddo
                    enddo
                    ! Directly fill the coefficient matrix with the eigenvectors from the diagonalization.

                    DEALLOCATE(WORK2)
                    CALL LogMemDealloc(this_routine,WORK2Tag)

                    DEALLOCATE(NOMSym)
                    CALL LogMemDealloc(this_routine,NOMSymTag)

                    DEALLOCATE(EvaluesSym)
                    CALL LogMemDealloc(this_routine,EvaluesSymTag)

                ELSEIF(NoSymBlock.eq.1) THEN
                    ! The eigenvalue is the lone value, while the eigenvector is 1.

                    Evalues(SymStartInd+1)=NatOrbMat(SymStartInd+1,SymStartInd+1)
                    NatOrbMat(SymStartInd+1,SymStartInd+1)=1.D0
                    WRITE(6,*) '*****'
                    WRITE(6,*) 'Symmetry ',Sym,' has only one orbital.'
                    WRITE(6,*) 'Copying diagonal element ,',SymStartInd+1,'to NatOrbMat'
                ENDIF

                Sym=Sym+1
            enddo
        enddo

        WRITE(6,*) 'Matrix diagonalised'
        CALL FLUSH(6)

        SumDiagTrace=0.D0
        do i=1,nBasis
            SumDiagTrace=SumDiagTrace+Evalues(i)
        enddo
        IF((ABS(SumDiagTrace-SumTrace)).gt.1E-5) THEN
            WRITE(6,*) 'Sum of diagonal NatOrbMat elements : ',SumTrace
            WRITE(6,*) 'Sum of eigenvalues : ',SumDiagTrace
            CALL Stop_All(this_routine,'The trace of the 1RDM matrix before diagonalisation is not equal to that after.')
        ENDIF

!        CALL halt_timer(DiagNatOrbMat_Time)

    END SUBROUTINE DiagRDM


    SUBROUTINE OrderRDM()
        USE sort_mod
        USE Logging , only : tTruncRODump
        USE RotateOrbsData , only : SymLabelList3
        IMPLICIT NONE
        INTEGER :: spin,i,j,ierr,StartSort,EndSort
        CHARACTER(len=*), PARAMETER :: this_routine='OrderRDM'
        INTEGER , ALLOCATABLE :: SymOrbsTemp(:)
        INTEGER :: SymOrbsTempTag
        

! Here, if symmetry is kept, we are going to have to reorder the eigenvectors according to the size of the eigenvalues, while taking
! the orbital labels (and therefore symmetries) with them. This will be put back into MP2VDM from MP2VDMTemp.

! Want to reorder the eigenvalues from largest to smallest, taking the eigenvectors with them and the symmetry as well.  
! If using spin orbitals, do this for the alpha spin and then the beta.
 
!        OrderCoeff_Time%timer_name='OrderCoeff'
!        CALL set_timer(OrderCoeff_Time,30)

!        WRITE(6,*) 'before sorting'
!        do i = 1, nBasis
!            write(6,*) i, SymLabelList2(i), INT(G1(SymLabelList2(i))%sym%S,4), evalues(i)
!            do j = 1, nBasis
!                write(6,'(F15.10)',advance='no') natorbmat(j,i)
!            enddo
!            WRITE(6,*) ''
!        enddo

        do i = 1, nBasis
            SymLabelList3(i) = SymLabelList2(i)
        enddo


        IF(tTruncRODump) THEN
            ! If we are truncating, the orbitals stay in this order, so we want to take their symmetries with them.
            ALLOCATE(SymOrbsTemp(nBasis),stat=ierr)
            CALL LogMemAlloc('SymOrbsTemp',nBasis,4,this_routine,SymOrbsTempTag,ierr)
            SymOrbsTemp(:)=0

            do i=1,nBasis
                SymOrbsTemp(i)=INT(G1(SymLabelList2(i))%sym%S,4)
            enddo

            do spin=1,2

                IF(spin.eq.1) THEN
                    StartSort=1
                    EndSort=(nBasis/2)
                ELSEIF(spin.eq.2) THEN
                    StartSort=(nBasis/2)+1
                    EndSort=nBasis
                ENDIF

                call sort (Evalues(startSort:endSort), &
                           natOrbMat(startSort:endSort, startSort:endSort), &
                           symOrbsTemp(startSort:endSort))

            enddo
               
        ELSE
            ! If we are not truncating, the orbitals get put back into their original order, so the symmetry information is still 
            ! correct, no need for the SymOrbs array.
            ! Instead, just take the labels of SymLabelList2 with them.

            do spin=1,2

                IF(spin.eq.1) THEN
                    StartSort=1
                    EndSort=nBasis/2
                ELSEIF(spin.eq.2) THEN
                    StartSort=(nBasis/2)+1
                    EndSort=nBasis
                ENDIF

                call sort (EValues(startSort:endSort), &
                           NatOrbMat(1:nBasis, startSort:endSort), &
                           symLabelList2(startSort:endSort))
            enddo 
            
        ENDIF

!        CALL halt_timer(OrderCoeff_Time)
!        WRITE(6,*) 'after sorting'
!        do i = 1, nBasis
!            write(6,*) i, SymLabelList2(i), INT(G1(SymLabelList2(i))%sym%S,4), evalues(i)
!            do j = 1, nBasis
!                write(6,'(F15.10)',advance='no') natorbmat(j,i)
!            enddo
!            WRITE(6,*) ''
!        enddo

        WRITE(6,*) 'Eigen-values: '
        do i=1,nBasis
            WRITE(6,*) Evalues(i)
        enddo

    END SUBROUTINE OrderRDM


    SUBROUTINE FillCoeffT1_RDM
        USE RotateOrbsData , only : CoeffT1,SymLabelList3,SymOrbs,SymOrbsTag,TruncEval,NoRotOrbs,EvaluesTrunc,EvaluesTruncTag
        USE Logging , only : tTruncRODump,tTruncDumpbyVal
        use util_mod, only: get_free_unit
        IMPLICIT NONE
        INTEGER :: l,k,i,j,NoRotAlphBet, io1, io2, SymLabelListInvNewTag
        CHARACTER(len=*), PARAMETER :: this_routine='FillCoeffT1_RDM'
        CHARACTER(len=5) :: Label
        CHARACTER(len=20) :: LabelFull
        REAL*8 :: OccEnergies(1:NoRotOrbs)
        INTEGER , ALLOCATABLE :: SymLabelListInvNew(:)
        INTEGER :: ierr, Orb, New_Pos
  
!        FillCoeff_Time%timer_name='FillCoeff'
!        CALL set_timer(FillCoeff_Time,30)

!        WRITE(6,*) 'NatOrbMat'
!        do i = 1, nBasis
!            do j = 1, nBasis
!                WRITE(6,'(F10.6)',advance='no') NatOrbMat(j,i)
!            enddo
!            WRITE(6,*) ''
!        enddo


! run through i, find out what orbital this corresponds to.  the position of this orbital is given by SymLabelList2Inv.
! the new position of this orbital is given by SymLabelList2Inv2.

        ALLOCATE(SymLabelListInvNew(nBasis),stat=ierr)
        CALL LogMemAlloc('SymLabelListInvNew',nBasis,4,this_routine,SymLabelListInvNewTag,ierr)
        SymLabelListInvNew(:)=0   

        do i=1,nBasis
            SymLabelListInvNew(SymLabelList2(i))=i
        enddo

        do i=1,nBasis
            do j = 1, nBasis
                Orb = SymLabelList3(j)      !This is the orbital the old position corresponds to.
                New_Pos = SymLabelListInvNew(Orb)   !This is the new position that orbital should go into.
                CoeffT1(New_Pos,i)=NatOrbMat(j,i)
            enddo
        enddo

        DEALLOCATE(SymLabelListInvNew)

!        do i=1,NoOrbs
!            WRITE(6,*) Evalues(i)
!        enddo
!        do i=1,NoOrbs
!            WRITE(6,*) NatOrbMat(:,i)
!        enddo
!        CALL FLUSH(6)
!        stop

        io2 = get_free_unit()
        OPEN(io2,FILE='EVALUES',status='unknown')
        WRITE(io2,*) nBasis
        k=0
        do i=1,nBasis,2
            k=k+1
            WRITE(io2,'(2I5,ES20.10,I5,A5,I5,ES20.10,I5)') (nBasis-i+1),i,Evalues(k),INT(G1(SymLabelList2(k))%Sym%S,4),'  *  ',&
                                                     &i+1,Evalues(k+(nBasis/2)),INT(G1(SymLabelList2(k+(nBasis/2)))%Sym%S,4)
        enddo
        CLOSE(io2)

!        CoeffT1(:,:) = 0.D0
!        do i = 1, nBasis
!            CoeffT1(i,i) = 1.D0
!        enddo

!        CALL HistNatOrbEvalues()

!        WRITE(6,*) 'NatOrbMat matrix'
!        do i=1,NoOrbs
!            WRITE(6,*) NatOrbMat(:,i)
!        enddo

!        OPEN(io1,FILE='TRANSFORMMAT',status='unknown')
!        do i=1,NoOrbs
!            do j=1,NoOrbs-NoFrozenVirt
!                WRITE(io1,*) i,j,CoeffT1(i,j)
!            enddo
!        enddo
!        CLOSE(io1)

!        CALL halt_timer(FillCoeff_Time)
 

    ENDSUBROUTINE FillCoeffT1_RDM

    
!This is an M^5 transform, which transforms all the two-electron integrals into the new basis described by the Coeff matrix.
!This is v memory inefficient and currently does not use any spatial symmetry information.
    SUBROUTINE Transform2ElIntsMemSave_RDM()
        USE RotateOrbsMod , only : FourIndInts
        INTEGER :: i,j,k,l,a,b,g,d,ierr,Temp4indintsTag,a2,b2,g2,d2
        REAL*8 , ALLOCATABLE :: Temp4indints(:,:)
#ifdef __CMPLX
        call stop_all('Transform2ElIntsMemSave_RDM', 'Rotating orbitals not implemented for complex orbitals.')
#endif
        
!        Transform2ElInts_Time%timer_name='Transform2ElIntsTime'
!        CALL set_timer(Transform2ElInts_time,30)

!Zero arrays from previous transform

 
        ALLOCATE(Temp4indints(nBasis,nBasis),stat=ierr)
        CALL LogMemAlloc('Temp4indints',nBasis**2,8,'Transform2ElIntsMemSave_RDM',Temp4indintsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All('Transform2ElIntsMemSave_RDM','Problem allocating memory to Temp4indints.')
 
        FourIndInts(:,:,:,:)=0.D0

!        WRITE(6,*) 'to here 01'
!        WRITE(6,*) tStoreSpinOrbs
!        CALL FLUSH(6)

! **************
! Calculating the two-transformed, four index integrals.

! The untransformed <alpha beta | gamma delta> integrals are found from UMAT(UMatInd(i,j,k,l,0,0)
! All our arrays are in spin orbitals - if tStoreSpinOrbs is false, UMAT will be in spatial orbitals - need to account for this.

! Running through 1,nBasis - the actual orbitals corresponding to that index are given by SymLabelList2

        do b=1,nBasis
            IF(.not.tStoreSpinOrbs) THEN
                b2=CEILING(REAL(SymLabelList2(b))/2.D0)
            ELSE
                b2=SymLabelList2(b)
            ENDIF
            do d=1,b
                IF(.not.tStoreSpinOrbs) THEN
                    d2=CEILING(REAL(SymLabelList2(d))/2.D0)
                ELSE
                    d2=SymLabelList2(d)
                ENDIF
                do a=1,nBasis
                    IF(.not.tStoreSpinOrbs) THEN
                        a2=CEILING(REAL(SymLabelList2(a))/2.D0)
                    ELSE
                        a2=SymLabelList2(a)
                    ENDIF
                    do g=1,a
                        IF(.not.tStoreSpinOrbs) THEN
                            g2=CEILING(REAL(SymLabelList2(g))/2.D0)
                        ELSE
                            g2=SymLabelList2(g)
                        ENDIF

!UMatInd in physical notation, but FourIndInts in chemical (just to make it more clear in these transformations).
!This means that here, a and g are interchangable, and so are b and d.
                        FourIndInts(a,g,b,d)=REAL(UMAT(UMatInd(a2,b2,g2,d2,0,0)),8)
                        FourIndInts(g,a,b,d)=REAL(UMAT(UMatInd(a2,b2,g2,d2,0,0)),8)
                        FourIndInts(a,g,d,b)=REAL(UMAT(UMatInd(a2,b2,g2,d2,0,0)),8)
                        FourIndInts(g,a,d,b)=REAL(UMAT(UMatInd(a2,b2,g2,d2,0,0)),8)
                    enddo
                enddo

                Temp4indints(:,:)=0.D0
                CALL DGEMM('T','N',nBasis,nBasis,nBasis,1.0,CoeffT1(:,:),nBasis,FourIndInts(1:nBasis,1:nBasis,b,d),nBasis,0.0,Temp4indints(1:nBasis,1:nBasis),nBasis)
                ! Temp4indints(i,g) comes out of here, so to transform g to k, we need the transpose of this.

                CALL DGEMM('T','T',nBasis,nBasis,nBasis,1.0,CoeffT1(:,:),nBasis,Temp4indints(1:nBasis,1:nBasis),nBasis,0.0,FourIndInts(1:nBasis,1:nBasis,b,d),nBasis)
                ! Get Temp4indits02(i,k)

                do i=1,nBasis
                    do k=1,i
                        FourIndInts(i,k,d,b)=FourIndInts(i,k,b,d)
                        FourIndInts(k,i,d,b)=FourIndInts(i,k,b,d)
                        FourIndInts(i,k,b,d)=FourIndInts(i,k,b,d)
                        FourIndInts(k,i,b,d)=FourIndInts(i,k,b,d)
                    enddo
                enddo
            enddo
        enddo
        
! Calculating the 3 transformed, 4 index integrals. 01=a untransformed,02=b,03=g,04=d
        do i=1,nBasis
            do k=1,i

                Temp4indints(:,:)=0.D0
                CALL DGEMM('T','N',nBasis,nBasis,nBasis,1.0,CoeffT1(:,:),nBasis,FourIndInts(i,k,1:nBasis,1:nBasis),nBasis,0.0,Temp4indints(1:nBasis,1:nBasis),nBasis)

                CALL DGEMM('T','T',nBasis,nBasis,nBasis,1.0,CoeffT1(:,:),nBasis,Temp4indints(1:nBasis,1:nBasis),nBasis,0.0,FourIndInts(i,k,1:nBasis,1:nBasis),nBasis)
                do l=1,nBasis
                    do j=1,l
                        FourIndInts(k,i,j,l)=FourIndInts(i,k,j,l)
                        FourIndInts(k,i,l,j)=FourIndInts(i,k,j,l)
                        FourIndInts(i,k,j,l)=FourIndInts(i,k,j,l)
                        FourIndInts(i,k,l,j)=FourIndInts(i,k,j,l)
                    enddo
                enddo
            enddo
        enddo

        DEALLOCATE(Temp4indints)
        CALL LogMemDeAlloc('Transform2ElIntsMemSave_RDM',Temp4indintsTag)
 
!        CALL halt_timer(Transform2ElInts_Time)

    END SUBROUTINE Transform2ElIntsMemSave_RDM


    SUBROUTINE CalcFOCKMatrix_RDM()
        USE SystemData , only : nBasis
        USE Logging , only : tRDMonfly
        INTEGER :: i,j,k,l,a,b,ierr,ArrDiagNewTag
        REAL*8 :: FOCKDiagSumHF,FOCKDiagSumNew
        CHARACTER(len=*) , PARAMETER :: this_routine='CalcFOCKMatrix_RDM'
        REAL*8 , ALLOCATABLE :: ArrDiagNew(:)
        !NEED TO FIX THIS!

! This subroutine calculates and writes out the fock matrix for the transformed orbitals.
! ARR is originally the fock matrix in the HF basis.
! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

! When transforming the orbitals into approximate natural orbitals, we want to save memory, so don't bother
! calculating the whole matrix, just the diagonal elements that we actually need.

    
        ALLOCATE(ArrDiagNew(nBasis),stat=ierr)
        CALL LogMemAlloc('ArrDiagNew',nBasis,8,this_routine,ArrDiagNewTag,ierr)
        ArrDiagNew(:)=0.D0                     

!        WRITE(6,*) 'The diagonal fock elements in the HF basis set'
!        do a=1,nBasis
!            WRITE(6,'(F20.10)',advance='no') Arr(a,2)
!            WRITE(6,*) Arr(:,2)
!        enddo


! First calculate the sum of the diagonal elements, ARR.
! Check if this is already being done.
        FOCKDiagSumHF=0.D0
        do a=1,nBasis        
            FOCKDiagSumHF=FOCKDiagSumHF+Arr(a,2)
        enddo

        WRITE(6,*) 'Sum of the fock matrix diagonal elements in the HF basis set = ',FOCKDiagSumHF

!        WRITE(6,*) 'Coeffs'
!        do i=1,NoOrbs
!            do j=1,NoOrbs
!                WRITE(6,'(F20.10)',advance='no') CoeffT1(j,i)
!            enddo
!            WRITE(6,*) ''
!        enddo

! Then calculate the fock matrix in the transformed basis, and the sum of the new diagonal elements.
! Our Arr in spin orbitals.
!        do j=1,NoOrbs
!            ArrNew(j,j)=Arr(2*j,2)
!        enddo

        FOCKDiagSumNew=0.D0
        do j=1,nBasis
            l=SymLabelList2(j)
            ArrDiagNew(l) = 0.D0
            do a=1,nBasis
                b=SymLabelList2(a)
                ArrDiagNew(l)=ArrDiagNew(l)+(CoeffT1(a,j)*ARR(b,2)*CoeffT1(a,j))
            enddo
            FOCKDiagSumNew=FOCKDiagSumNew+(ArrDiagNew(l))
        enddo
        ! If we are truncation the virtual space, only the unfrozen entries will be transformed.
        

        WRITE(6,*) 'Sum of the fock matrix diagonal elements in the transformed basis set = ',FOCKDiagSumNew

!        WRITE(6,*) 'The fock matrix for the transformed orbitals'
!        do j=1,NoOrbs
!            do i=1,NoOrbs
!                WRITE(6,'(F20.10)',advance='no') ArrNew(i,j)
!            enddo
!            WRITE(6,*) ''
!        enddo

!        WRITE(6,*) 'BRR then ARR before being changed',nBasis
!        do i=1,nBasis
!            WRITE(6,*) i,BRR(i),ARR(i,1),ARR(BRR(i),2),ArrDiagNew(i)
!        enddo
       
!        WRITE(6,*) 'to here',NoDumpTruncs,NoOrbs 
!        CALL FLUSH(6)

! Refill ARR(:,1) (ordered in terms of energies), and ARR(:,2) (ordered in terms of orbital number).
! ARR(:,2) needs to be ordered in terms of symmetry and then energy (like SymLabelList), so currently this ordering will not be 
! correct when reading in qchem INTDUMPS as the orbital number ordering is by energy.

!        IF(NoDumpTruncs.le.1) THEN
! If we are only writing out 1 ROFCIDUMP or we are not truncating at all - can refill ARR etc.            

        do j=1,nBasis
            ARR(j,2)=ArrDiagNew(j)
            ARR(j,1)=ArrDiagNew(BRR(j))
        enddo

!        WRITE(6,*) 'BRR then ARR after being changed'
!        do i=1,nBasis
!            WRITE(6,*) i,BRR(i),ARR(i,1),ARR(BRR(i),2)
!        enddo
!        CALL FLUSH(6)
!        stop       

        DEALLOCATE(ArrDiagNew)
        CALL LogMemDealloc(this_routine,ArrDiagNewTag)

!        WRITE(6,*) 'end of calcfockmatrix'
!        call flush(6)

    ENDSUBROUTINE CalcFOCKMatrix_RDM


    SUBROUTINE RefillUMATandTMAT2D_RDM()
        USE RotateOrbsMod , only : FourIndInts, PrintROFCIDUMP, PrintRepeatROFCIDUMP
        INTEGER :: l,k,j,i,a,b,g,d,c,nBasis2,TMAT2DPartTag,ierr
        REAL*8 :: NewTMAT
        REAL*8 , ALLOCATABLE :: TMAT2DPart(:,:)
#ifdef __CMPLX
        call stop_all('RefillUMATandTMAT2D_RDM', 'Rotating orbitals not implemented for complex orbitals.')
#endif

        IF(tStoreSpinOrbs) THEN
            ALLOCATE(TMAT2DPart(nBasis,nBasis),stat=ierr)
            CALL LogMemAlloc('TMAT2DPart',nBasis*nBasis,8,'RefillUMAT_RDM',TMAT2DPartTag,ierr)
        ELSE
            ALLOCATE(TMAT2DPart(nBasis,nBasis),stat=ierr)
            CALL LogMemAlloc('TMAT2DPart',nBasis*nBasis,8,'RefillUMAT_RDM',TMAT2DPartTag,ierr)
        ENDIF
        TMAT2DPart(:,:)=0.D0

!        RefillUMAT_Time%timer_name='RefillUMATandTMAT'
!        CALL set_timer(RefillUMAT_Time,30)

!        do i = 1,nBasis
!            WRITE(6,*) SymLabelList2(i)
!        enddo


! Make the UMAT elements the four index integrals.  These are calculated by transforming the HF orbitals using the
! coefficients that have been found
        IF(NoDumpTruncs.le.1) THEN
            do l=1,nBasis
                IF(.not.tStoreSpinOrbs) THEN
                    d=CEILING(REAL(SymLabelList2(l))/2.D0)
                ELSE
                    d=SymLabelList2(l)
                ENDIF
                do k=1,nBasis

                    IF(.not.tStoreSpinOrbs) THEN
                        b=CEILING(REAL(SymLabelList2(k))/2.D0)
                    ELSE
                        b=SymLabelList2(k)
                    ENDIF
     
                    do j=1,l

                        IF(.not.tStoreSpinOrbs) THEN
                            g=CEILING(REAL(SymLabelList2(j))/2.D0)
                        ELSE
                            g=SymLabelList2(j)
                        ENDIF
                        do i=1,k

                            IF(.not.tStoreSpinOrbs) THEN
                                a=CEILING(REAL(SymLabelList2(i))/2.D0)
                            ELSE
                                a=SymLabelList2(i)
                            ENDIF

!The FourIndInts are in chemical notation, the UMatInd in physical.                            
                            UMAT(UMatInd(a,b,g,d,0,0))=FourIndInts(i,j,k,l)
                        enddo
                    enddo
                enddo
            enddo
        ENDIF

! Also calculate the 2 index integrals, and make these the elements of the TMAT2D matrix.
! TMAT2D is in spin orbitals.

!        WRITE(6,*) 'TMAT2D before transformation' 
!        do l=1,nBasis
!            do k=1,nBasis
!                WRITE(6,'(F10.6)',advance='no') REAL(TMAT2D(k,l),8)
!                WRITE(6,*) REAL(TMAT2D(:,:),8)
!            enddo
!            WRITE(6,*) ''
!        enddo

!        WRITE(6,*) 'SpatOrbs',SpatOrbs
!        WRITE(6,*) 'NoRotOrbs',NoRotOrbs
!        CALL FLUSH(6)


        do a=1,nBasis
            do k=1,nBasis
                i=SymLabelList2(k)
                NewTMAT=0.D0
                do b=1,nBasis
                    d=SymLabelList2(b)
                    NewTMAT=NewTMAT+(CoeffT1(b,k)*REAL(TMAT2D(d,a),8))
                enddo
                TMAT2DPart(i,a)=NewTMAT
            enddo
        enddo

        do k=1,nBasis
            do l=1,nBasis
                j=SymLabelList2(l)
                NewTMAT=0.D0
                do a=1,nBasis
                    c=SymLabelList2(a)
                    NewTMAT=NewTMAT+(CoeffT1(a,l)*TMAT2DPart(k,c))
                enddo
                TMAT2D(k,j)=NewTMAT
            enddo
        enddo

    
!        WRITE(6,*) 'TMAT2D after transformation'
!        do l=1,nBasis
!            do k=1,nBasis
!                WRITE(6,'(F10.6)',advance='no') REAL(TMAT2D(k,l),8)
!            enddo
!            WRITE(6,*) ''
!        enddo
!        CALL FLUSH(6)
!        stop

        DEALLOCATE(TMAT2DPart)
        CALL LogMemDeAlloc('RefillUMAT_RDM',TMAT2DPartTag)


!        CALL set_timer(RefillUMAT_Time,30)

        WRITE(6,'(A,I5,A)') ' Printing the new ROFCIDUMP file for a truncation of ',NoFrozenVirt,' orbitals.'
        IF(tROFciDump.and.(NoDumpTruncs.gt.1)) THEN
            CALL PrintRepeatROFCIDUMP()
        ELSEIF(tROFciDUmp) THEN
            CALL PrintROFCIDUMP_RDM()
        ENDIF


    ENDSUBROUTINE RefillUMATandTMAT2D_RDM


    SUBROUTINE PrintROFCIDUMP_RDM()
!This prints out a new FCIDUMP file in the same format as the old one.
        use util_mod, only: get_free_unit
        INTEGER :: i,j,k,l,iunit
        CHARACTER(len=5) :: Label
        CHARACTER(len=20) :: LabelFull

        IF(tStoreSpinOrbs) THEN
            NoOrbs = nBasis
        ELSE
            NoOrbs = nBasis/2
        ENDIF

!        PrintROFCIDUMP_Time%timer_name='PrintROFCIDUMP'
!        CALL set_timer(PrintROFCIDUMP_Time,30)

        Label=''
        LabelFull=''
        WRITE(Label,'(I5)') NoFrozenVirt
        LabelFull='ROFCIDUMP-'//adjustl(Label)

        iunit = get_free_unit()
        OPEN(iunit,FILE=LabelFull,STATUS='unknown')
        
        WRITE(iunit,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',NoOrbs,',NELEC=',NEl,',MS2=',LMS,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,NoOrbs
            IF(tStoreSpinOrbs) THEN
                WRITE(iunit,'(I1,A1)',advance='no') (INT(G1(i)%sym%S,4)+1),','
            ELSE
                WRITE(iunit,'(I1,A1)',advance='no') (INT(G1(i*2)%sym%S,4)+1),','
            ENDIF
        enddo
        WRITE(iunit,*) ''
        IF(tStoreSpinOrbs) THEN
            WRITE(iunit,'(A7,I1,A11)') 'ISYM=',1,' UHF=.TRUE.'
        ELSE
            WRITE(iunit,'(A7,I1)') 'ISYM=',1
        ENDIF
        WRITE(iunit,'(A5)') '&END'
       
        do i=1,NoOrbs
            do j=1,i
                do l=1,NoOrbs
!                    Sym=IEOR(INT(G1(j*2)%sym%S,4),IEOR(INT(G1(k*2)%sym%S,4),INT(G1(i*2)%sym%S,4)))
                    ! Potential to put symmetry in here, have currently taken it out, because when we're only printing non-zero values,
                    ! it is kind of unnecessary - although it may be used to speed things up.
                    do k=1,NoOrbs
!                        Syml=INT(G1(l*2)%sym%S,4)
!                        IF((Syml.eq.Sym).and.((REAL(UMat(UMatInd(i,j,k,l,0,0)),8)).ne.0.D0)) &

!UMatInd is in physical notation <ij|kl>, but the indices printed in the FCIDUMP are in chemical notation (ik|jl).
                        IF((ABS(REAL(UMat(UMatInd(i,j,k,l,0,0)),8))).ne.0.D0) &
                                        &WRITE(iunit,'(F21.12,4I3)') REAL(UMat(UMatInd(i,j,k,l,0,0)),8),i,k,j,l 
 
                    enddo
                enddo
           enddo
        enddo

! TMAT2D stored as spin orbitals
        do i=1,NoOrbs
            ! Symmetry?
            do k=1,NoOrbs
                IF(tStoreSpinOrbs) THEN
                    IF((REAL(TMAT2D(i,k),8)).ne.0.D0) WRITE(iunit,'(F21.12,4I3)') REAL(TMAT2D(i,k),8),i,k,0,0
                ELSE
                    IF((REAL(TMAT2D(2*i,2*k),8)).ne.0.D0) WRITE(iunit,'(F21.12,4I3)') REAL(TMAT2D(2*i,2*k),8),i,k,0,0
                ENDIF
            enddo
        enddo

! ARR has the energies of the orbitals (eigenvalues).  ARR(:,2) has ordering we want.
! ARR is stored as spin orbitals.

        do k=1,NoOrbs
            IF(tStoreSpinOrbs) THEN
                WRITE(iunit,'(F21.12,4I3)') Arr(k,2),k,0,0,0
            ELSE
                WRITE(iunit,'(F21.12,4I3)') Arr(2*k,2),k,0,0,0
            ENDIF
        enddo

        WRITE(iunit,'(F21.12,4I3)') ECore,0,0,0,0
        
        CALL FLUSH(iunit)

        CLOSE(iunit)

!        CALL halt_timer(PrintROFCIDUMP_Time)


    ENDSUBROUTINE PrintROFCIDUMP_RDM


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

            write(6,*) 'here o.k 05'
            call flush(6)

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
        USE UMatCache, only: GTID
        INTEGER :: i,j,k,l,Ind2,Ind1,i2,j2,k2,l2
        REAL*8 :: RDMEnergy, Trace_1RDM, Trace_2RDM, Coul, Exch 
        REAL*8 :: Trace_2RDM_New, RDMEnergy1El, RDMEnergy2El


        IF(tTurnStoreSpinOff) tStoreSpinOrbs=.false.

!First of all 'normalise' the reduced density matrices.
!These are not even close to being normalised at the moment, because of the way they are calculated on the fly.
!They should be calculated from a normalised wavefunction.
!But we know that the trace of the one electron reduced density matrix must be equal to the number of the electrons.
!We can use this to find the factor we must divide the 1RDM through by.

!We also know that the trace of the two electron reduced density matrix must be equal to the number of electron pairs 
!in the system = 1/2 N ( N - 1), so we can do the same for the 2RDM.

        WRITE(6,*) 'inside this routine'
        CALL FLUSH(6)

        Trace_1RDM = 0.D0
        Trace_2RDM = 0.D0

        RDMEnergy = 0.D0
        RDMEnergy1El = 0.D0
        RDMEnergy2El = 0.D0

        IF(iProcIndex.eq.0) THEN
            do i = 1, nBasis
                Trace_1RDM = Trace_1RDM + NatOrbMat(i,i)
                Trace_2RDM = Trace_2RDM + AllTwoElRDM(i,i)
            enddo

            WRITE(6,*) 'trace o.k'
            WRITE(6,*) 'Trace_1RDM',Trace_1RDM
            WRITE(6,*) 'Trace_2RDM',Trace_2RDM

            WRITE(6,*) '((nBasis*(nBasis-1))/2)',((nBasis*(nBasis-1))/2)
            WRITE(6,*) 'tTurnStoreSpinOff',tTurnStoreSpinOff
            CALL FLUSH(6)

            do i = 1, nBasis
                NatOrbMat(SymLabelListInv(i),SymLabelListInv(i)) = &
                        NatOrbMat(SymLabelListInv(i),SymLabelListInv(i)) * ( REAL(NEl,8) / Trace_1RDM )  
                do j = 1, nBasis                            
                    IF(i.ne.j) NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) = &
                            NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) * ( REAL(NEl,8) / Trace_1RDM )  
 
                enddo
            enddo

!            WRITE(6,*) 'before'
!            do i = 1, ((nBasis*(nBasis-1))/2)
!                do j = 1, ((nBasis*(nBasis-1))/2)
!                    IF(AllTwoElRDM(i,j).ne.0) WRITE(6,*) i,j,AllTwoElRDM(i,j)
!                enddo
!            enddo

            Trace_2RDM_New = 0.D0
            do i = 1, ((nBasis*(nBasis-1))/2)
                AllTwoElRDM(i,i) = AllTwoElRDM(i,i) * ( (0.5 * REAL((NEl * (NEl - 1)) ,8)) / Trace_2RDM )  
                Trace_2RDM_New = Trace_2RDM_New + AllTwoElRDM(i,i)
                do j = 1, ((nBasis*(nBasis-1))/2)
                    IF(i.ne.j) AllTwoElRDM(i,j) = AllTwoElRDM(i,j) * ( (0.5 * REAL((NEl * (NEl - 1)) ,8)) / Trace_2RDM )  
                enddo
            enddo

            WRITE(6,*) 'Trace_2RDM_New',Trace_2RDM_New

!            WRITE(6,*) 'after'
!            do i = 1, ((nBasis*(nBasis-1))/2)
!                do j = 1, i
!                    WRITE(6,*) i,j,AllTwoElRDM(i,j),AllTwoElRDM(j,i)
!                enddo
!            enddo

            !Need to multiply each element of the 1 electron reduced density matrices by NEl / Trace_1RDM,
            !and then add it's contribution to the energy.
            do i = 1, nBasis 

                i2 = gtID(i)
 
                do j = 1, nBasis
                    !i and j correspond to *orbitals*.

                    j2 = gtID(j)

!SymLabelListInv(j) = x, gives the position of orbital j in NatOrbMat (orbital j is in position x).
!We want to find orbital j, because we're multiplying it by TMAT2D(i,j) where i and j are the *orbitals* not the position.
!                    IF(i.ne.j) NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) = &
!                            NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) * ( REAL(NEl,8) / Trace_1RDM )  
                    
!                    WRITE(6,*) 'NatOrbMat(i,j)',NatOrbMat(j,i)
!                    CALL FLUSH(6)

                    !TMAT2D always stored in spin orbitals so o.k. 
                    RDMEnergy1El = RDMEnergy1El + ( REAL(TMAT2D(i,j),8)  &
                                                * NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) )

!                    WRITE(6,*) '01 Adding:',( REAL(TMAT2D(i,j),8)  &
!                                                * NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) )

!                    IF(i.ne.j) THEN
!                        NatOrbMat(SymLabelList2(i),SymLabelList2(j)) = &
!                                NatOrbMat(SymLabelList2(i),SymLabelList2(j)) * ( REAL(NEl,8) / Trace_1RDM )  
!
!                        RDMEnergy = RDMEnergy + ( REAL(TMAT2D(j,i),8) * &
!                                                        NatOrbMat(SymLabelList2(i),SymLabelList2(j)) )
!                    ENDIF

!                    WRITE(6,*) 'RDMEnergy 01',RDMEnergy


                    do k = 1, i-1

                        IF(k.eq.i) CYCLE

                        k2 = gtID(k)

                        Ind1 = ( ( (MAX(i,k)-2) * (MAX(i,k)-1) ) / 2 ) + MIN(i,k)
!                            WRITE(6,*) 'Ind1',Ind1,'i',i,'k',k

                        do l = 1, j - 1

                            IF(l.eq.j) CYCLE

                            l2 = gtID(l)

                            Ind2 = ( ( (MAX(j,l)-2) * (MAX(j,l)-1) ) / 2 ) + MIN(j,l)
!                                    WRITE(6,*) 'Ind2',Ind2,'j',j,'l',l

!                                    WRITE(6,*) 'Ind1',Ind1,'Ind2',Ind2

!                                    IF(Ind2.ge.Ind1) AllTwoElRDM(Ind1,Ind2) = AllTwoElRDM(Ind1,Ind2) * ( REAL( (0.5 * NEl * (NEl - 1)) ,8) / Trace_2RDM )  
!                                    WRITE(6,*) 'RDMEnergy 02',RDMEnergy
!                                    WRITE(6,*) 'i',i,'j',j,'k',k,'l',l

                            ! AllTwoElRDM etc will be in spin orbitals, but UMAT is in spatial, need to correct for this.

                            !UmatInd uses *physical* notation.
!                                    RDMEnergy2El = RDMEnergy2El + ( 0.50 * REAL(UMAT(UMatInd(i2,k2,j2,l2,0,0)),8) * AllTwoElRDM(Ind1,Ind2) )
!                                    RDMEnergy2El = RDMEnergy2El + ( 0.50 * (REAL(UMAT(UMatInd(i2,k2,j2,l2,0,0)),8) &
!                                                                            - REAL(UMAT(UMatInd(i2,k2,l2,j2,0,0)),8)) * AllTwoElRDM(Ind1,Ind2) )

                            IF( (G1(i)%Ms .eq. G1(j)%Ms) .and. &
                                 (G1(k)%Ms .eq. G1(l)%Ms) ) THEN
                                Coul = REAL(UMAT(UMatInd(i2,k2,j2,l2,0,0)),8)
                            ELSE
                                Coul = 0.D0
                            ENDIF

                            IF( (G1(i)%Ms .eq. G1(l)%Ms) .and. &
                                 (G1(k)%Ms .eq. G1(j)%Ms) ) THEN
                                Exch = REAL(UMAT(UMatInd(i2,k2,l2,j2,0,0)),8)
                            ELSE
                                Exch = 0.D0
                            ENDIF

!                            RDMEnergy2El = RDMEnergy2El + ( (REAL(UMAT(UMatInd(i2,k2,j2,l2,0,0)),8) &
!                                                                    - REAL(UMAT(UMatInd(i2,k2,l2,j2,0,0)),8)) * AllTwoElRDM(Ind1,Ind2) )
                            RDMEnergy2El = RDMEnergy2El + ( 0.5 * ( Coul - Exch ) * AllTwoElRDM(Ind1,Ind2) )
 
                            ! We are just running from k = 1, i-1 (for example) so we are adding this contribution in twice (simply omit the 
                            ! factor of 0.5).
                            ! This is so that we don't accidentally count the 'diagonal' (i.e. Ind1 = Ind2) terms twice.


!                                    WRITE(6,*) '02 Adding:',( 0.50 * REAL(UMAT(UMatInd(i2,k2,j2,l2,0,0)),8) * AllTwoElRDM(Ind1,Ind2) )
!                                    IF(( 0.50 * REAL(UMAT(UMatInd(i2,k2,j2,l2,0,0)),8) * AllTwoElRDM(Ind1,Ind2) ).gt.1E100) THEN
!                                        WRITE(6,*) 'Ind1',Ind1,'Ind2',Ind2
!                                        WRITE(6,*) 'UMAT',REAL(UMAT(UMatInd(i2,k2,j2,l2,0,0)),8)
!                                        WRITE(6,*) 'AllTwoElRDM',AllTwoElRDM(Ind1,Ind2)
!                                    ENDIF

!                                    IF(i2.ne.j2) THEN
!                                        RDMEnergy = RDMEnergy + ( 0.50 * REAL(UMAT(UMatInd(j2,k2,i2,l2,0,0)),8) * AllTwoElRDM(Ind2,Ind1) )
!                                    ENDIF

!                                    IF(k2.ne.l2) THEN
!                                        RDMEnergy = RDMEnergy + ( 0.50 * REAL(UMAT(UMatInd(i2,l2,j2,k2,0,0)),8) * AllTwoElRDM(Ind2,Ind1) )
!                                    ENDIF

!                                    IF((i2.ne.j2).and.(k2.ne.l2)) THEN
!                                        RDMEnergy = RDMEnergy + ( 0.50 * REAL(UMAT(UMatInd(j2,l2,i2,k2,0,0)),8) * AllTwoElRDM(Ind2,Ind1) )
!                                    ENDIF

!                                    WRITE(6,*) 'RDMEnergy 04',RDMEnergy

                        enddo
                    enddo
                enddo
            enddo
            RDMEnergy = RDMEnergy1El + RDMEnergy2El

            WRITE(6,*) '******'
            WRITE(6,*) 'The ENERGY calculated using the REDUCED DENSITY MATRICES is:',RDMEnergy
            WRITE(6,*) 'The contribution from the 1 el RDM:',RDMEnergy1El
            WRITE(6,*) 'The contribution from the 2 el RDM:',RDMEnergy2El

            CALL Test_Energy_Calc()

        ENDIF

        tTurnStoreSpinOff=.false.
        IF(.not.tStoreSpinOrbs) THEN
            tTurnStoreSpinOff=.true.
            tStoreSpinOrbs=.true.
        ENDIF


    END SUBROUTINE Calc_Energy_from_RDM


    SUBROUTINE Test_Energy_Calc()
! This routine calculates the energy based on the simple expression Energy = Sum_IJ C_I C_J H_IJ where I and J are determinants    
! (rather than occupied orbitals), and H_IJ is the Hamiltonian element between them.
        USE FciMCData , only : TotWalkers,CurrentDets
        USE Determinants, only : get_helement
        USE bit_reps , only : extract_bit_rep, extract_sign
        USE UMatCache, only: GTID
        INTEGER :: I, J, nI(NEl), nJ(NEl), FlagsI, FlagsJ, IC, Ex(2,2)
        INTEGER :: k,l,k2,l2,a2,b2
        LOGICAL :: tParity
        INTEGER, DIMENSION(lenof_sign) :: SignI, SignJ
        REAL*8 :: Test_Energy, Sum_Coeffs, SignIreal, SignJreal 
        REAL*8 :: Test_Energy_1El, Test_Energy_2El, ParityFactor
        REAL*8 :: Coul, Exch
        HElement_t :: H_IJ

        WRITE(6,*) '****************'
        WRITE(6,*) 'TESTING ENERGY CALCULATION **** '

        IF(tTurnStoreSpinOff) tStoreSpinOrbs=.false.

        Sum_Coeffs = 0.D0
        do i = 1, TotWalkers
            call extract_sign (CurrentDets(:,I), SignI)
            Sum_Coeffs = Sum_Coeffs + ( REAL(SignI(1)) * REAL(SignI(1)) )
        enddo
        Sum_Coeffs = SQRT( Sum_Coeffs )

! Just need to run over all occupied determinants, find the matrix element between them, and sum them all up.
        Test_Energy = 0.D0
        Test_Energy_1El = 0.D0
        Test_Energy_2El = 0.D0

        do I = 1, TotWalkers
            call extract_bit_rep (CurrentDets(:,I), nI, SignI, FlagsI)
            SignIreal = REAL(SignI(1)) * ( 1.D0 / Sum_Coeffs )

            do J = 1, TotWalkers
                call extract_bit_rep (CurrentDets(:,J), nJ, SignJ, FlagsJ)

                SignJreal = REAL(SignJ(1)) * ( 1.D0 / Sum_Coeffs )

!                WRITE(6,*) 'I',I,'J',J
!                WRITE(6,*) 'nI',nI
!                WRITE(6,*) 'nJ',nJ
!                WRITE(6,*) 'ilutI',CurrentDets(:,I)
!                WRITE(6,*) 'ilutj',CurrentDets(:,J)
!                CALL FLUSH(6)

                H_IJ = get_helement(nI, nJ, CurrentDets(:,I), CurrentDets(:,J), IC)

!                WRITE(6,*) 'H_IJ',H_IJ
!                CALL FLUSH(6)

                IF(IC.eq.0) THEN
                    !Add to both 1El and 2El contributions.

                    do k = 1, NEl
                        k2 = gtID(nI(k))

                        Test_Energy_1El = Test_Energy_1El + ( SignIreal * SignIreal * REAL( TMAT2D(nI(k),nI(k)) ) )

                        ! For i,j -> i,j type doubles, UMat(i,j,i,j) = UMat(i,j,j,i) (so <ij||ij> = 0.
                        do l = 1, NEl
                            IF(k.eq.l) CYCLE

                            l2 = gtID(nI(l))

                            IF( G1(nI(k))%Ms .eq. G1(nI(l))%Ms ) THEN
                                Exch = REAL(UMAT(UMatInd(k2,l2,l2,k2,0,0)),8)
                            ELSE
                                Exch = 0.D0
                            ENDIF

                            Test_Energy_2El = Test_Energy_2El + (0.5 * SignIreal * SignIreal &
                                                                    * (REAL(UMAT(UMatInd(k2,l2,k2,l2,0,0)),8) - Exch) )
!!                                                                          * (REAL(UMAT(UMatInd(k2,l2,k2,l2,0,0)),8)) )
                        enddo
                    enddo

                ELSEIF(IC.eq.1) THEN

                    Ex(:,:) = 0
                    Ex(1,1) = 1
                    tParity = .false.
                    ParityFactor = 1.D0
                    CALL GetExcitation(nI,nJ,NEl,Ex,tParity)
! Ex(1,:) comes out as the orbital(s) excited from, Ex(2,:) comes out as the orbital(s) excited to.                    
                    IF(tParity) ParityFactor=-1.D0

                    Test_Energy_1El = Test_Energy_1El + (ParityFactor * SignIreal * SignJreal * REAL( TMAT2D(Ex(1,1),Ex(2,1)) ) )

                    IF(G1(Ex(1,1))%Ms.eq.G1(Ex(2,1))%Ms) THEN

                        do k = 1, NEl
                            IF(nI(k).eq.Ex(1,1)) CYCLE

                            k2 = gtID(nI(k))

                            IF(G1(nI(k))%Ms.eq.G1(Ex(1,1))%Ms) THEN
                                Exch = REAL(UMAT(UMatInd(k2,Ex(1,1),Ex(2,1),k2,0,0)),8)
                            ELSE
                                Exch = 0.D0
                            ENDIF

                            Test_Energy_2El = Test_Energy_2El + ( ParityFactor * SignIreal * SignJreal &
                                                                    * (REAL(UMAT(UMatInd(k2,Ex(1,1),k2,Ex(2,1),0,0)),8) - Exch)) 
!                                                                    * (REAL(UMAT(UMatInd(k2,Ex(1,1),k2,Ex(2,1),0,0)),8)) )
!                            WRITE(6,*) '*'
!                            WRITE(6,*) k2,Ex(1,1),k2,Ex(2,1),REAL(UMAT(UMatInd(k2,Ex(1,1),k2,Ex(2,1),0,0)),8)
!                            WRITE(6,*) k2,Ex(1,1),Ex(2,1),k2,REAL(UMAT(UMatInd(k2,Ex(1,1),Ex(2,1),k2,0,0)),8)
                        enddo
                    ENDIF

                ELSEIF(IC.eq.2) THEN

                    Ex(:,:) = 0
                    Ex(1,1) = 2
                    tParity = .false.
                    ParityFactor = 1.D0
                    CALL GetExcitation(nI,nJ,NEl,Ex,tParity)
! Ex(1,:) comes out as the orbital(s) excited from, Ex(2,:) comes out as the orbital(s) excited to.                    
                    IF(tParity) ParityFactor=-1.D0

                    k2 = gtID(Ex(1,1))
                    l2 = gtID(Ex(1,2))
                    a2 = gtID(Ex(2,1))
                    b2 = gtID(Ex(2,2))

                    IF( (G1(Ex(1,1))%Ms .eq. G1(Ex(2,1))%Ms) .and. &
                         (G1(Ex(1,2))%Ms .eq. G1(Ex(2,2))%Ms) ) THEN
                        Coul = REAL(UMAT(UMatInd(k2,l2,a2,b2,0,0)),8)
                    ELSE
                        Coul = 0.D0
                    ENDIF

                    IF( (G1(Ex(1,1))%Ms .eq. G1(Ex(2,2))%Ms) .and. &
                         (G1(Ex(1,2))%Ms .eq. G1(Ex(2,1))%Ms) ) THEN
                        Exch = REAL(UMAT(UMatInd(k2,l2,b2,a2,0,0)),8)
                    ELSE
                        Exch = 0.D0
                    ENDIF

                    Test_Energy_2El = Test_Energy_2El + ( ParityFactor * SignIreal * SignJreal * (Coul - Exch) )
!                    WRITE(6,*) '*'
!                    WRITE(6,*) k2,l2,a2,b2,REAL(UMAT(UMatInd(k2,l2,a2,b2,0,0)),8)
!                    WRITE(6,*) k2,l2,a2,b2,REAL(UMAT(UMatInd(k2,l2,b2,a2,0,0)),8)

                ENDIF

                Test_Energy = Test_Energy + ( SignIreal * SignJreal * REAL(H_IJ) )

!                WRITE(6,*) 'Test_Energy',Test_Energy
!                CALL FLUSH(6)

            enddo

        enddo

        WRITE(6,*) 'The test Energy contribution to the 1 electron part = ',Test_Energy_1El
        WRITE(6,*) 'The test Energy contribution to the 2 electron part = ',Test_Energy_2El
        WRITE(6,*) 'The sum of these = ',Test_Energy_1El + Test_Energy_2El

        WRITE(6,*) 'The TEST ENERGY calculated from Sum_IJ C_I C_J H_IJ  = ',Test_Energy
        WRITE(6,*) ''
        WRITE(6,*) 'The 2 electron part should be : ',Test_Energy - Test_Energy_1El

        tTurnStoreSpinOff=.false.
        IF(.not.tStoreSpinOrbs) THEN
            tTurnStoreSpinOff=.true.
            tStoreSpinOrbs=.true.
        ENDIF

    END SUBROUTINE Test_Energy_Calc


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


