MODULE nElRDMMod
! This file contains the routines used to find the full n electron reduced density matrix (nElRDM).
! This is done on the fly to avoid having to histogram the full wavefunction which is extremely 
! time and memory inefficient.
! In this way, these routines differ slightly from those in NatOrbsMod (which take a histogrammed 
! wavefunction usually truncated around double excitations and form the one electron RDM) but the basic 
! formula is still the same.
! For example, he elements of the one electron reduced density matrix are given by:
! 1RDM_pq   = < Psi | a_p+ a_q | Psi > 
! where Psi is the full wavefunction, a_p+ is the creation operator and a_q is the annihilation operator.
!           = < sum_i c_i D_i | a_p+ a_q | sum_j c_j D_j >
! The elements 1RDM_pq therefore come from the sum of the contributions c_i*c_j from all pairs of 
! determinants D_i and D_j which are related by a single excitation between p and q.
! This can be generalised for the nElRDM, where n = 1 or 2.

! The algorithm for calculating the nElRDM on the fly will be similar in nature to the direct 
! annihilation routines.
! Each set of processors has a list of determinants.  
! These each select the first of these D_i, generate all the allowed single or double excitations D_j, 
! and order them in terms of the processor they would be on if they are occupied.
! The excitations are then sent to the relevant processor along with the original determinant, 
! and it's sign c_i.
! Each processor then receives nProcessors sets of excitations (D_j's).
! For each of these they binary search their list of occupied determinants.  
! If an excitation D_j is found, c_i.c_j is added to the matrix element corresponding to the orbitals 
! involved in the excitation.

! NOTE: There will be possible speed ups considering the fact that the 1RDM is symmetrical.
! Can initially find all elements and average the two values pq and qp (more accurate??).
! But should put a condition into the excitaiton generator so that only single excitations with q > p 
! are generated.

! By finding the full 1RDM, we have the ability to derive the natural orbitals as well as electron 
! densities etc.
        
        USE Global_Utilities
        USE Parallel
        USE bit_reps , only : NIfTot, NIfDBO, decode_bit_det
        USE IntegralsData , only : UMAT
        USE UMatCache , only : UMatInd
        USE SystemData , only : NEl,nBasis,tStoreSpinOrbs, G1, BRR, lNoSymmetry, ARR
        USE SystemData , only : tUseMP2VarDenMat, Ecore, LMS, tHPHF
        USE NatOrbsMod , only : NatOrbMat,NatOrbMatTag,Evalues,EvaluesTag
        USE CalcData , only : MemoryFacPart
        USE constants , only : n_int, dp, Root2, sizeof_int
        USE OneEInts , only : TMAT2D
        USE FciMCData , only : MaxWalkersPart, MaxSpawned, Spawned_Parents, PreviousCycles,&
                               Spawned_Parents_Index, Spawned_ParentsTag, AccumRDMNorm_Inst,&
                               Spawned_Parents_IndexTag, Iter, AccumRDMNorm, AvNoatHF,&
                               iLutRef, tSinglePartPhase, AllAccumRDMNorm
        USE Logging , only : RDMExcitLevel, tROFciDUmp, NoDumpTruncs, tExplicitAllRDM, &
                             tHF_S_D_Ref, tHF_Ref_Explicit, tHF_S_D, tPrint1RDM
        USE RotateOrbsData , only : CoeffT1Tag, tTurnStoreSpinOff, NoFrozenVirt
        USE RotateOrbsData , only : SymLabelCounts2,SymLabelList2,SymLabelListInv,NoOrbs, SpatOrbs
        USE util_mod , only : get_free_unit
        IMPLICIT NONE
        INTEGER , ALLOCATABLE :: Sing_InitExcSlots(:),Sing_ExcList(:)
        INTEGER , ALLOCATABLE :: Doub_InitExcSlots(:),Doub_ExcList(:)
        INTEGER(kind=n_int) , ALLOCATABLE :: Sing_ExcDjs(:,:),Sing_ExcDjs2(:,:)
        INTEGER(kind=n_int) , ALLOCATABLE :: Doub_ExcDjs(:,:),Doub_ExcDjs2(:,:)
        INTEGER :: Sing_ExcDjsTag,Sing_ExcDjs2Tag,aaaa_RDMTag,All_aaaa_RDMTag
        INTEGER :: Doub_ExcDjsTag,Doub_ExcDjs2Tag,UMATTempTag
        INTEGER :: Energies_unit, ActualStochSign_unit, abab_RDMTag, All_abab_RDMTag
        INTEGER :: abba_RDMTag, All_abba_RDMTag, NoSymLabelCounts
        REAL(dp) , ALLOCATABLE :: aaaa_RDM(:,:), abab_RDM(:,:), abba_RDM(:,:)
        REAL(dp) , ALLOCATABLE :: All_aaaa_RDM(:,:),All_abab_RDM(:,:), All_abba_RDM(:,:)
        REAL(dp) , ALLOCATABLE :: UMATTemp(:,:)
        REAL(dp) :: OneEl_Gap,TwoEl_Gap, Normalisation,Trace_2RDM_Inst, Trace_2RDM, Trace_1RDM, norm
        LOGICAL :: tFinalRDMEnergy, tCalc_RDMEnergy
        type(timer), save :: nElRDM_Time, FinaliseRDM_time, RDMEnergy_time

    contains

    SUBROUTINE InitRDM()
! This routine initialises any of the arrays needed to calculate the reduced density matrix.    
! It is used for both the explicit and stochastic RDMs.
        USE NatOrbsMod , only : SetupNatOrbLabels 
        USE RotateOrbsMod , only : SymLabelCounts2Tag,NoRotOrbs
        USE RotateOrbsMod , only : SymLabelList2Tag,SymLabelListInvTag
        USE Logging , only : tDo_Not_Calc_RDMEnergy, tDiagRDM, tReadRDMs, &
                            TPopsFile, tno_RDMs_to_read, twrite_RDMs_to_read,&
                            tWriteMultRDMs
        USE CalcData , only : tRegenDiagHEls
        use DetBitOps , only : TestClosedShellDet
        implicit none
        INTEGER :: ierr,i, MemoryAlloc, MemoryAlloc_Root
        CHARACTER(len=*), PARAMETER :: this_routine='InitRDM'

! First thing is to check we're not trying to fill the RDMs in a way that is 
! not compatible with the code (not every case has been accounted for yet).
#ifdef __CMPLX
        CAll Stop_All(this_routine,'Filling of reduced density matrices not working with &
                                    &complex walkers yet.')
#endif
        ! Only spatial orbitals for the 2-RDMs (and F12).
        if((.not.TestClosedShellDet(iLutRef)).and.(RDMExcitLevel.ne.1)) &
            call stop_all(this_routine,'2-RDM calculations not set up for open shell systems.')
                
        if(tStoreSpinOrbs.and.(RDMExcitLevel.ne.1)) &
            call stop_all(this_routine,'2-RDM calculations not set up for systems stored &
                                        &as spin orbitals.')
    
        ! The averaged coefficients used for calculating the RDMs are stored with the CurrentH 
        ! array (which stores the diagonal H elements).  Will need to set up a new array or something
        ! if we're not storing these Kii values.
        if(tRegenDiagHEls) &
            call stop_all(this_routine,'RDMs not currently set up for regenerating the &
                                    &diagonal H elements.  This should not be difficult though.')

        if(tExplicitAllRDM) then
            write(6,'(A)') " Explicitly calculating the reduced density matrices from the &
                                                        &FCIQMC wavefunction."
        else
            write(6,'(A)') " Stochastically calculating the reduced density matrices from the &
                            &FCIQMC wavefunction" 
            write(6,'(A)') " (incl. explicit connections to the reference determinant)."
        endif

        IF(RDMExcitLevel.eq.1) THEN
            tCalc_RDMEnergy = .false.

        ELSE
! If the RDMExcitLevel is 2 or 3 - and we're calculating the 2-RDM, 
! then we automatically calculate the energy unless we specifically say not to.
            IF(tDo_Not_Calc_RDMEnergy) THEN
                tCalc_RDMEnergy = .false.            
            ELSE
                tCalc_RDMEnergy = .true.
                WRITE(6,'(A)') ' Calculating the energy from the reduced &
                &density matrix, this requires the 2 electron RDM from which the 1-RDM can also be constructed.'
            ENDIF
        ENDIF

        ! Have not got HPHF working with the explicit or truncated methods yet.
        ! Neither of these would be too difficult to implement.
        if(tHPHF.and.tExplicitAllRDM) CALL Stop_All('InitRDM',&
                'HPHF not set up with the explicit calculation of the RDM.')

        if(tHPHF.and.(tHF_S_D_Ref.or.tHF_S_D)) CALL Stop_All('InitRDM',&
                'HPHF not set up when doing a HF, S, D calculation.')

        ! Can't diagonalise the non-hermitian matrix.                
        if(tDiagRDM.and.(tHF_S_D_Ref.or.tHF_Ref_Explicit)) then
            write(6,*) 'Ignoring request to diagonalise the 1-RDM calculated using the HF or HF, S, D &
                &as a reference - this is not an appropriate matrix for natural orbitals.'
            tDiagRDM = .false.
        endif

        SpatOrbs=nBasis/2
        if(tStoreSpinOrbs) then
            NoOrbs=nBasis
        else
            NoOrbs=SpatOrbs
        endif

! Here we're allocating arrays for the actual calculation of the RDM.

        MemoryAlloc = 0
        MemoryAlloc_Root = 0            ! Memory allocated in bytes.

! First for the storage of the actual 1- or 2-RMD.

        IF(RDMExcitLevel.eq.1) THEN
! This is the AllnElRDM, called NatOrbMat simply because we use the natural 
! orbital routines to diagonalise etc.        
! We don't have an instantaneous 1-RDM.
            ALLOCATE(NatOrbMat(NoOrbs,NoOrbs),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating NatOrbMat array,')
            CALL LogMemAlloc('NatOrbMat',NoOrbs**2,8,this_routine,NatOrbMatTag,ierr)
            NatOrbMat(:,:)=0.D0

            MemoryAlloc = MemoryAlloc + ( NoOrbs * NoOrbs * 8 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( NoOrbs * NoOrbs * 8 ) 
        ELSE
            ! If we're calculating the 2-RDM, the 1-RDM does not need to be calculated as well 
            ! because all its info is in the 2-RDM anyway.

            ! The 2-RDM of the type alpha alpha alpha alpha ( = beta beta beta beta).
            ! These *do not* include any 2-RDM(i,j,a,b) terms where i=j or a=b (if they're the same 
            ! spin this can't happen).
            ALLOCATE(aaaa_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating aaaa_RDM array,')
            CALL LogMemAlloc('aaaa_RDM',(((SpatOrbs*(SpatOrbs-1))/2)**2),8,this_routine,aaaa_RDMTag,ierr)
            aaaa_RDM(:,:)=0.D0

            ! The 2-RDM of the type alpha beta beta alpha ( = beta alpha alpha beta).
            ! These also *do not* also include 2-RDM(i,j,a,b) terms where i=j or a=b (these are the same as the abab elements).
            ALLOCATE(abba_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating abba_RDM array,')
            CALL LogMemAlloc('abba_RDM',(((SpatOrbs*(SpatOrbs-1))/2)**2),8,this_routine,abba_RDMTag,ierr)
            abba_RDM(:,:)=0.D0

            MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 

            ! The 2-RDM of the type alpha beta alpha beta ( = beta alpha beta alpha).
            ! These *do* include 2-RDM(i,j,a,b) terms where i=j or a=b, if they're different spin this 
            ! is possible - hence the slightly different size to the aaaa array.
            ALLOCATE(abab_RDM(((SpatOrbs*(SpatOrbs+1))/2),((SpatOrbs*(SpatOrbs+1))/2)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating abab_RDM array,')
            CALL LogMemAlloc('abab_RDM',(((SpatOrbs*(SpatOrbs+1))/2)**2),8,this_routine,abab_RDMTag,ierr)
            abab_RDM(:,:)=0.D0

            MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 ) * 8 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 ) * 8 ) 

            IF(iProcIndex.eq.0) THEN
                ! Each of these currently need to be stored on the root as well as each node.
                ! This allows us to separately calculate the instantaneous energy.
                ! TODO : Option to only calculate the accumulated RDMs - cut storage on the root in half.
                ALLOCATE(All_aaaa_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating All_aaaa_RDM array,')
                CALL LogMemAlloc('All_aaaa_RDM',(((SpatOrbs*(SpatOrbs-1))/2)**2),8,this_routine,All_aaaa_RDMTag,ierr)
                All_aaaa_RDM(:,:)=0.D0

                ALLOCATE(All_abab_RDM(((SpatOrbs*(SpatOrbs+1))/2),((SpatOrbs*(SpatOrbs+1))/2)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating All_abab_RDM array,')
                CALL LogMemAlloc('All_abab_RDM',(((SpatOrbs*(SpatOrbs+1))/2)**2),8,this_routine,All_abab_RDMTag,ierr)
                All_abab_RDM(:,:)=0.D0

                ALLOCATE(All_abba_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating All_abba_RDM array,')
                CALL LogMemAlloc('All_abba_RDM',(((SpatOrbs*(SpatOrbs-1))/2)**2),8,this_routine,All_abba_RDMTag,ierr)
                All_abba_RDM(:,:)=0.D0

                MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
                MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 ) * 8 ) 

                if(tDiagRDM.or.tPrint1RDM) then
                    ! Still need to allocate 1-RDM to get nat orb occupation numbers.
                    ALLOCATE(NatOrbMat(SpatOrbs,SpatOrbs),stat=ierr)
                    IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating NatOrbMat array,')
                    CALL LogMemAlloc('NatOrbMat',SpatOrbs**2,8,this_routine,NatOrbMatTag,ierr)
                    NatOrbMat(:,:)=0.D0

                    MemoryAlloc_Root = MemoryAlloc_Root + ( SpatOrbs * SpatOrbs * 8 ) 
                endif
            ENDIF
        ENDIF            

! We then need to allocate the arrays for excitations etc when doing the explicit all calculation.        
        IF(tExplicitAllRDM) THEN            

            ! We always calculate the single stuff - and if RDMExcitLevel is 1, this is all,
            ! otherwise calculate the double stuff too.

! This array actually contains the excitations in blocks of the processor they will be sent to.        
! Only needed if the 1-RDM is the only thing being calculated.
            ALLOCATE(Sing_ExcDjs(0:NIfTot,NINT((NEl*nBasis)*MemoryFacPart)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Sing_ExcDjs array.')
            CALL LogMemAlloc('Sing_ExcDjs',NINT(NEl*nBasis*MemoryFacPart)*(NIfTot+1),&
                                            size_n_int,this_routine,Sing_ExcDjsTag,ierr)

            ALLOCATE(Sing_ExcDjs2(0:NIfTot,NINT((NEl*nBasis)*MemoryFacPart)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Sing_ExcDjs2 array.')
            CALL LogMemAlloc('Sing_ExcDjs2',NINT(NEl*nBasis*MemoryFacPart)*(NIfTot+1),&
                                            size_n_int,this_routine,Sing_ExcDjs2Tag,ierr)

            Sing_ExcDjs(:,:)=0
            Sing_ExcDjs2(:,:)=0

            MemoryAlloc = MemoryAlloc + ( (NIfTot + 1) * NINT((NEl*nBasis)*MemoryFacPart) * size_n_int * 2 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( (NIfTot + 1) * NINT((NEl*nBasis)*MemoryFacPart) * size_n_int * 2 ) 

! We need room to potentially generate N*M single excitations but these will be 
! spread across each processor.        

            OneEl_Gap=(REAL(NEl)*REAL(nBasis)*MemoryFacPart)/REAL(nProcessors)

! This array contains the initial positions of the excitations for each processor.
            ALLOCATE(Sing_InitExcSlots(0:(nProcessors-1)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Sing_InitExcSlots array,')
            do i=0,nProcessors-1
                Sing_InitExcSlots(i)=NINT(OneEl_Gap*i)+1
            enddo

! This array contains the current position of the excitations as they're added.
            ALLOCATE(Sing_ExcList(0:(nProcessors-1)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Sing_ExcList array,')
            Sing_ExcList(:)=Sing_InitExcSlots(:)

            IF(RDMExcitLevel.ne.1) THEN
! This array actually contains the excitations in blocks of the processor 
! they will be sent to.        
                ALLOCATE(Doub_ExcDjs(0:NIfTot,NINT(((NEl*nBasis)**2)*MemoryFacPart)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Doub_ExcDjs array.')
                CALL LogMemAlloc('Doub_ExcDjs',NINT(((NEl*nBasis)**2)*MemoryFacPart)&
                                *(NIfTot+1),size_n_int,this_routine,Doub_ExcDjsTag,ierr)

                ALLOCATE(Doub_ExcDjs2(0:NIfTot,NINT(((NEl*nBasis)**2)*MemoryFacPart)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Doub_ExcDjs2 array.')
                CALL LogMemAlloc('Doub_ExcDjs2',NINT(((NEl*nBasis)**2)*MemoryFacPart)&
                                *(NIfTot+1),size_n_int,this_routine,Doub_ExcDjs2Tag,ierr)

                MemoryAlloc = MemoryAlloc + ( (NIfTot + 1) * NINT(((NEl*nBasis)**2)*MemoryFacPart) * size_n_int * 2 ) 
                MemoryAlloc_Root = MemoryAlloc_Root + ( (NIfTot + 1) * NINT(((NEl*nBasis)**2)*MemoryFacPart) * size_n_int * 2 ) 

! We need room to potentially generate (N*M)^2 double excitations but these 
! will be spread across each processor.        
                TwoEl_Gap=(((REAL(NEl)*REAL(nBasis))**2)*MemoryFacPart)/REAL(nProcessors)

                Doub_ExcDjs(:,:)=0
                Doub_ExcDjs2(:,:)=0

! This array contains the initial positions of the excitations for each processor.
                ALLOCATE(Doub_InitExcSlots(0:(nProcessors-1)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Doub_InitExcSlots array,')
                do i=0,nProcessors-1
                    Doub_InitExcSlots(i)=NINT(TwoEl_Gap*i)+1
                enddo

! This array contains the current position of the excitations as they're added.
                ALLOCATE(Doub_ExcList(0:(nProcessors-1)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Doub_ExcList array,')
                Doub_ExcList(:)=Doub_InitExcSlots(:)
            ENDIF

        ELSEIF(.not.tHF_Ref_Explicit) THEN

! Finally, we need to hold onto the parents of the spawned particles.            
! This is not necessary if we're doing completely explicit calculations.
            ALLOCATE(Spawned_Parents(0:(NIfDBO+1),MaxSpawned),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Spawned_Parents array,')
            CALL LogMemAlloc('Spawned_Parents',MaxSpawned*(NIfDBO+2),size_n_int,&
                                                this_routine,Spawned_ParentsTag,ierr)
            ALLOCATE(Spawned_Parents_Index(2,MaxSpawned),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Spawned_Parents_Index array,')
            CALL LogMemAlloc('Spawned_Parents_Index',MaxSpawned*2,4,this_routine,&
                                                        Spawned_Parents_IndexTag,ierr)

            MemoryAlloc = MemoryAlloc + ( (NIfTot + 2) * MaxSpawned * size_n_int ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( (NIfTot + 2) * MaxSpawned * size_n_int ) 

            MemoryAlloc = MemoryAlloc + ( 2 * MaxSpawned * 4 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( 2 * MaxSpawned * 4 ) 

        ENDIF

        if(iProcIndex.eq.0) then
            write(6,"(A,F14.6,A,F14.6,A)") " Main RDM memory arrays consists of : ", &
                    & REAL(MemoryAlloc_Root,dp)/1048576.D0," Mb/Processor on the root, and ", &
                    & REAL(MemoryAlloc,dp)/1048576.D0," Mb/Processor on other processors."
        endif

        ! These parameters are set for the set up of the symmetry arrays, which are later used 
        ! for the diagonalisation / rotation of the 1-RDMs.

        if(tStoreSpinOrbs) then
            NoSymLabelCounts = 16
        else
            NoSymLabelCounts = 8
        endif

        IF((RDMExcitLevel.eq.1).or.tDiagRDM.or.tPrint1RDM) THEN
            ! These arrays contain indexing systems to order the 1-RDM orbitals in terms of 
            ! symmetry.
            ! This allows the diagonalisation of the RDMs to be done in symmetry blocks (a lot 
            ! quicker/easier).
            ! The 2-RDM does not need to be reordered as it's never diagonalised. 

            ALLOCATE(SymLabelCounts2(2,NoSymLabelCounts),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating SymLabelCounts2 array,')
            CALL LogMemAlloc('SymLabelCounts2',2*NoSymLabelCounts,4,this_routine,SymLabelCounts2Tag,ierr)
            SymLabelCounts2(:,:)=0

            ALLOCATE(SymLabelList2(NoOrbs),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating SymLabelList2 array,')
            CALL LogMemAlloc('SymLabelList2',NoOrbs,4,this_routine,SymLabelList2Tag,ierr)
            SymLabelList2(:)=0                     
     
            ALLOCATE(SymLabelListInv(NoOrbs),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating SymLabelListInv array,')
            CALL LogMemAlloc('SymLabelListInv',NoOrbs,4,this_routine,SymLabelListInvTag,ierr)
            SymLabelListInv(:)=0   

            IF((iProcIndex.eq.0).and.tDiagRDM) THEN
                ALLOCATE(Evalues(NoOrbs),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Evalues array,')
                CALL LogMemAlloc('Evalues',NoOrbs,8,this_routine,EvaluesTag,ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for Evalues failed.")
                Evalues(:)=0.D0
            ENDIF

            ! This routine actually sets up the symmetry labels for the 1-RDM.
            ! TODO : Merge this routine (and rotations later) with the NatOrbs file.
            CALL SetUpSymLabels_RDM() 

        ENDIF            

        if(iProcIndex.eq.0) write(6,'(A)') " RDM memory allocation successful... "                    

        ! Open file to keep track of RDM Energies (if they're being calculated). 
        IF((iProcIndex.eq.0).and.tCalc_RDMEnergy) THEN
            Energies_unit = get_free_unit()
            OPEN(Energies_unit,file='RDMEnergies',status='unknown')

            WRITE(Energies_unit, "(A1,3A30)") '#','Iteration','RDM Energy - Inst','RDM Energy - Accum'
        ENDIF
        tFinalRDMEnergy = .false.

        Trace_2RDM = 0.0_dp
        Trace_2RDM_Inst = 0.0_dp
        AccumRDMNorm = 0.0_dp
        AccumRDMNorm_Inst = 0.0_dp
        ! AccumRDMNorm is the normalisation for the case where we're using a limited reference to calculate
        ! the RDM.

        ! Reads in the RDMs from a previous calculation, sets the accumulating normalisations, 
        ! writes out the starting energy.
        if(tReadRDMs) then
            if(tSinglePartPhase) then
                write(6,'(A)') 'WARNING - Asking to read in the RDMs, but not varying shift from &
                                & the beginning of the calculation.'
                write(6,'(A)') 'Ignoring the request to read in the RDMs and starting again.'
                tReadRDMs = .false.
            else
                call Read_In_RDMs()
            endif
        endif

        ! By default, if we're writing out a popsfile (and doing an RDM calculation), we also 
        ! write out the unnormalised RDMs that can be read in when restarting a calculation.
        ! If the NORDMSTOREAD option is on, these wont be printed.  
        if(TPopsFile.and.(.not.tno_RDMs_to_read)) twrite_RDMs_to_read = .true.

        nElRDM_Time%timer_name='nElRDMTime'
        FinaliseRDM_Time%timer_name='FinaliseRDMTime'
        RDMEnergy_Time%timer_name='RDMEnergyTime'

    END SUBROUTINE InitRDM

    subroutine zero_rdms()
        implicit none

        if(RDMExcitLevel.eq.1) then
            NatOrbMat(:,:) = 0.0_dp
        else
            aaaa_RDM(:,:) = 0.0_dp
            abab_RDM(:,:) = 0.0_dp
            abba_RDM(:,:) = 0.0_dp
            All_aaaa_RDM(:,:) = 0.0_dp
            All_abab_RDM(:,:) = 0.0_dp
            All_abba_RDM(:,:) = 0.0_dp
        endif

        Trace_2RDM = 0.0_dp
        Trace_2RDM_Inst = 0.0_dp
        AccumRDMNorm = 0.0_dp
        AccumRDMNorm_Inst = 0.0_dp

    end subroutine

    subroutine Read_In_RDMs()
! Reads in the arrays to restart the RDM calculation (and continue accumulating).
! These arrays are not normalised, so the trace is also calculated.
! The energy is then calculated (if required) from the RDMs read in only.
        use Logging , only : IterRDMonFly
        implicit none
        logical :: exists_aaaa,exists_abab,exists_abba,exists_one
        integer :: RDM_unit, FileEnd
        integer :: i,j,a,b,Ind1,Ind2
        real(dp) :: Temp_RDM_Element

        if(iProcIndex.eq.0) then 

            if(RDMExcitLevel.eq.1) then

                write(6,'(A)') ' Reading in the 1-RDM'

                ! The OneRDM will have been printed exactly as is.  Without having been 
                ! made hermitian, without being normalised, and in spatial orbitals if 
                ! tStoreSpinOrbs is false.

                INQUIRE(FILE='OneRDM_POPS',EXIST=exists_one)
                if(exists_one) then
                    RDM_unit = get_free_unit()
                    open(RDM_unit,FILE='OneRDM_POPS',status='old',form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i,j,Temp_RDM_Element
                        if(FileEnd.gt.0) call stop_all("Read_In_RDMs","Error reading OneRDM_POPS")
                        if(FileEnd.lt.0) exit

                        NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) = Temp_RDM_Element
                    enddo
                    close(RDM_unit)
                else
                    call stop_all('Read_In_RDMs','Attempting to read in the OneRDM, but &
                                    &the OneRDM_POPS file does not exist.')
                endif                                    

            else

                write(6,'(A)') ' Reading in the 2-RDMs'

                ! The TwoRDMs will have been printed exactly as they were.  Without having been 
                ! made hermitian, without being normalised, and in spatial orbitals. 

                ! Only read in the 2-RDMs (the 1-RDM becomes redundant).
                INQUIRE(FILE='TwoRDM_POPS_aaaa',EXIST=exists_aaaa)
                INQUIRE(FILE='TwoRDM_POPS_abab',EXIST=exists_abab)
                INQUIRE(FILE='TwoRDM_POPS_abba',EXIST=exists_abba)
                if(exists_aaaa.and.exists_abab.and.exists_abba) THEN
                    ! All TOREAD RDM files are present - read in.
                    RDM_unit = get_free_unit()
                    open(RDM_unit,FILE='TwoRDM_POPS_aaaa',status='old',form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i,j,a,b,Temp_RDM_Element 
                        if(FileEnd.gt.0) call stop_all("Read_In_RDMs","Error reading TwoRDM_POPS_aaaa")
                        if(FileEnd.lt.0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        All_aaaa_RDM(Ind1,Ind2) = Temp_RDM_Element
                    enddo
                    close(RDM_unit)

                    open(RDM_unit,FILE='TwoRDM_POPS_abab',status='old',form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i,j,a,b,Temp_RDM_Element 
                        if(FileEnd.gt.0) call stop_all("Read_In_RDMs","Error reading TwoRDM_POPS_abab")
                        if(FileEnd.lt.0) exit

                        Ind1 = ( ( (j-1) * j ) / 2 ) + i
                        Ind2 = ( ( (b-1) * b ) / 2 ) + a
                        All_abab_RDM(Ind1,Ind2) = Temp_RDM_Element
                    enddo
                    close(RDM_unit)

                    open(RDM_unit,FILE='TwoRDM_POPS_abba',status='old',form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i,j,a,b,Temp_RDM_Element 
                        if(FileEnd.gt.0) call stop_all("Read_In_RDMs","Error reading TwoRDM_POPS_abba")
                        if(FileEnd.lt.0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        All_abba_RDM(Ind1,Ind2) = Temp_RDM_Element
                    enddo
                    close(RDM_unit)

                else
                    write(6,*) 'exists_aaaa',exists_aaaa
                    write(6,*) 'exists_abab',exists_abab
                    write(6,*) 'exists_abba',exists_abba
                    call flush(6)
                    CALL Stop_All('Read_in_RDMs',"Attempting to read in the RDMs, &
                                    &but at least one of the TwoRDM_a***_TOREAD files are missing.")
                endif

            endif
        endif

        ! Calculate the energy for the matrices read in (if we're calculating more than the 1-RDM).
        if(tCalc_RDMEnergy) then
            call Calc_Energy_from_RDM()
        endif

! Continue calculating the RDMs from the first iteration when the POPSFILES (and RDMs) are read in.
! This overwrites the iteration number put in the input.
        IterRDMonFly = 0

    end subroutine Read_In_RDMs

    SUBROUTINE SetUpSymLabels_RDM() 
! This routine just sets up the symmetry labels so that
! the orbitals are ordered according to symmetry (all beta then all alpha if spin orbs).
        USE sort_mod
        USE UMatCache , only : GTID
        IMPLICIT NONE
        INTEGER , ALLOCATABLE :: SymOrbs(:)
        INTEGER :: LabOrbsTag, SymOrbsTag, ierr, i , j 
        INTEGER :: lo, hi, Symi, SymCurr, Symi2, SymCurr2
        CHARACTER(len=*) , PARAMETER :: this_routine = 'SetUpSymLabels_RDM'

        ! This is only allocated temporarily to be used to order the orbitals by.
        ALLOCATE(SymOrbs(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymOrbs',NoOrbs,4,this_routine,SymOrbsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for SymOrbs failed.")

! Brr has the orbital numbers in order of energy... i.e Brr(2) = the orbital index 
! with the second lowest energy.
! Brr is always in spin orbitals.
!        write(6,*) 'BRR'
!        do i=1,nBasis
!            WRITE(6,*) i, BRR(i), INT(G1(BRR(i))%sym%S,4)
!        enddo
!        CALL FLUSH(6)
!        CALL Stop_All('','')

! Now we want to put the spatial orbital index, followed by the symmetry.        
        SymLabelList2(:) = 0
        SymOrbs(:)=0

! *** STEP 1 *** Fill SymLabelList2.
! find the orbitals and order them in terms of symmetry.
        do i=1,SpatOrbs
            if(tStoreSpinOrbs) then
                ! for open shell systems, all alpha are followed by all beta.
                SymLabelList2(i)=BRR(2*i)
                SymOrbs(i)=INT(G1(BRR(2*i))%sym%S,4)
                SymLabelList2(i+SpatOrbs)=BRR((2*i)-1)
                SymOrbs(i+SpatOrbs)=INT(G1(BRR((2*i)-1))%sym%S,4)
            else
                SymLabelList2(i)=gtID(BRR(2*i))
                SymOrbs(i)=INT(G1(BRR(2*i))%sym%S,4)
                ! Orbital BRR(2*i) for i = 1 will be the beta orbital with the 
                ! second lowest energy - want the spatial orbital index to go with this.
                ! G1 also in spin orbitals - get symmetry of this beta orbital, will 
                ! be the same as the spatial orbital.
            endif
        enddo

        call sort (SymOrbs(1:SpatOrbs), SymLabelList2(1:SpatOrbs))
        ! Sorts SymLabelList2 according to the order of SymOrbs (i.e. in terms of symmetry). 
        if(tStoreSpinOrbs) &
            call sort (SymOrbs(SpatOrbs+1:nBasis), SymLabelList2(SpatOrbs+1:nBasis))
            ! Also do this for the beta set if spin orbitals.

!*** STEP 2 *** Fill SymLabelCounts2.
!SymLabelCounts(1,:) contains the position in SymLabelList2 where the symmetry index starts,
!SymLabelCounts(2,:) contains the number of orbitals in that symmetry index.
!Again if spin orbs, all alpha are followed by all beta - i.e. first 8 refer to alpha, second 8 to beta.

        IF(lNoSymmetry) THEN
            ! if we are ignoring symmetry, all orbitals essentially have symmetry 0.
            SymLabelCounts2(1,1) = 1
            SymLabelCounts2(2,1) = SpatOrbs
            if(tStoreSpinOrbs) then
                SymLabelCounts2(1,9) = SpatOrbs+1
                SymLabelCounts2(2,9) = SpatOrbs
            endif
        ELSE 
            ! otherwise we run through the occupied orbitals, counting the number with 
            ! each symmetry and noting where in SymLabelList2 each symmetry block starts.
            SymCurr = 0
            SymLabelCounts2(1,1) = 1
            if(tStoreSpinOrbs) then
                SymCurr2 = 0
                SymLabelCounts2(1,9) = SpatOrbs + 1
            endif
            do i=1,SpatOrbs
                if(tStoreSpinOrbs) then
                    Symi=SymOrbs(i)
                    Symi2=SymOrbs(i + SpatOrbs)
                else
                    Symi=SymOrbs(i)
                endif
                SymLabelCounts2(2,(Symi+1))= SymLabelCounts2(2,(Symi+1))+1
                IF(Symi.ne.SymCurr) THEN
                    do j = SymCurr + 1, Symi
                        SymLabelCounts2(1,(j+1))=i
                    enddo
                    SymCurr=Symi
                ENDIF
                if(tStoreSpinOrbs) then
                    SymLabelCounts2(2,(Symi2+9))= SymLabelCounts2(2,(Symi2+9))+1
                    IF(Symi2.ne.SymCurr2) THEN
                        do j = SymCurr2 + 1, Symi2
                            SymLabelCounts2(1,(j+9))=i + SpatOrbs
                        enddo
                        SymCurr2=Symi2
                    ENDIF
                endif
            enddo
        ENDIF

        ! Go through each symmetry group, making sure the orbitals are 
        ! ordered lowest to highest within each symmetry.
        do i=1,NoSymLabelCounts
            IF(SymLabelCounts2(2,i).ne.0) THEN
                lo = SymLabelCounts2(1, i)
                hi = lo + SymLabelCounts2(2, i) - 1
                call sort (SymLabelList2 (lo:hi))
            ENDIF
        enddo

        ! Construct the inverse matrix.  While SymLabelList2 takes a position and tells us 
        ! what orbital is in it, we also might need to take an orbital and find out what 
        ! position to put its contribution in.
        do i=1,NoOrbs
            SymLabelListInv(SymLabelList2(i))=i
        enddo

! Deallocate the arrays just used in this routine.
        DEALLOCATE(SymOrbs)
        CALL LogMemDealloc(this_routine,SymOrbsTag)

!        WRITE(6,*) 'Sym Label Counts'
!        do i=1,8
!            WRITE(6,*) i,SymLabelCounts2(1,i),SymLabelCounts2(2,i)
!        enddo
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), &
!                               & and their symmetries according to G1'
!        do i=1,SpatOrbs
!            WRITE(6,*) i,SymLabelList2(i),INT(G1(2*SymLabelList2(i))%sym%S,4)
!        enddo
!        WRITE(6,*) 'i','ARR(SymLabelList2(i),1)','ARR(SymLabelList2(i),2)','Sym'
!        do i=1,SpatOrbs
!             WRITE(6,*) i,ARR(2*SymLabelList2(i),1),ARR(2*SymLabelList2(i),2),&
!                                               INT(G1(2*SymLabelList2(i))%sym%S,4)
!        enddo
!
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and its inverse'
!        do i=1,SpatOrbs
!            WRITE(6,*) SymLabelList2(i),SymLabelListInv(i)
!        enddo
!        CALL FLUSH(6)
!        CALL Stop_All('SetUpSymLabels_RDM','Checking orbital labelling.')

    END SUBROUTINE SetUpSymLabels_RDM


    subroutine DeAlloc_Alloc_SpawnedParts()
! When calculating the RDMs, we need to store the parent from which a child is spawned along with the 
! children in the spawned array.
! This means a slightly larger array is communicated between processors, which there is no point in doing 
! for the first part of the calculation.
! When we start calculating the RDMs this routine is called and the SpawnedParts array is made larger to 
! accommodate the parents.
        USE FciMCData , only : SpawnVec, SpawnVec2, SpawnVecTag, SpawnVec2Tag, &
                               SpawnedParts, SpawnedParts2
        implicit none                               
        INTEGER :: ierr                               
        CHARACTER(len=*), PARAMETER :: this_routine='DeAlloc_Alloc_SpawnedParts'

        DEALLOCATE(SpawnVec)
        CALL LogMemDealloc(this_routine,SpawnVecTag)
        DEALLOCATE(SpawnVec2)
        CALL LogMemDealloc(this_routine,SpawnVec2Tag)
 
        ALLOCATE(SpawnVec(0:(NIftot+NIfDBO+2),MaxSpawned),stat=ierr)
        CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+NIfDBO+3),size_n_int,this_routine,SpawnVecTag,ierr)
        ALLOCATE(SpawnVec2(0:(NIfTot+NIfDBO+2),MaxSpawned),stat=ierr)
        CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+NIfDBO+3),size_n_int,this_routine,SpawnVec2Tag,ierr)

!        SpawnVec(:,:)=0
!        SpawnVec2(:,:)=0

!Point at correct spawning arrays
        SpawnedParts=>SpawnVec
        SpawnedParts2=>SpawnVec2

        WRITE(6,'(A54,F10.4,A4,F10.4,A13)') 'Memory requirement for spawned arrays increased from ',&
                                        REAL(((NIfTot+1)*MaxSpawned*2*size_n_int),dp)/1048576.D0,' to ',&
                                        REAL(((NIfTot+NIfDBO+3)*MaxSpawned*2*size_n_int),dp)/1048576.D0, ' Mb/Processor'

    end subroutine DeAlloc_Alloc_SpawnedParts

! The following extract_bit_rep_rdm routines both extract the bit representation 
! of the current determinant, and at the same time add in the current determinants contributions 
! to the diagonal elements of the RDMs.

    subroutine extract_bit_rep_rdm_diag_no_rdm(iLutnI, CurrH_I, nI, SignI, &
                                                FlagsI, IterRDMStartI, AvSignI, Store)
! This is just the standard extract_bit_rep routine for when we're not calculating the RDMs.    
        use constants , only : dp, n_int, lenof_sign
        use SystemData , only : NEl
        use bit_reps , only : NIfTot, extract_bit_rep, extract_bit_rep_rdm
        use FciMCData , only : excit_gen_store_type, NCurrH
        use DetBitOps , only : TestClosedShellDet
        implicit none
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp) , intent(in) :: CurrH_I(NCurrH)
        integer, intent(out) :: nI(nel), FlagsI
        integer, dimension(lenof_sign), intent(out) :: SignI
        real(dp) , intent(out) :: IterRDMStartI, AvSignI
        type(excit_gen_store_type), intent(inout), optional :: Store

        call extract_bit_rep (iLutnI, nI, SignI, FlagsI, Store)

        IterRDMStartI = 0.0_dp
        AvSignI = 0.0_dp

    end subroutine extract_bit_rep_rdm_diag_no_rdm

    subroutine extract_bit_rep_rdm_diag_norm(iLutnI, CurrH_I, nI, SignI, &
                                                FlagsI, IterRDMStartI, AvSignI, Store)
! While extracting the orbitals from the bit representation of the determinant, we 
! simulaneously add in the contribution of each orbital to the diagonal elements of the RDMs.
! This is used for a standard full space RDM calculation without HPHF.
        use constants , only : dp, n_int, lenof_sign
        use SystemData , only : NEl
        use bit_reps , only : NIfTot, extract_bit_rep, extract_bit_rep_rdm
        use FciMCData , only : excit_gen_store_type, NCurrH
        use DetBitOps , only : TestClosedShellDet
        implicit none
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp) , intent(in) :: CurrH_I(NCurrH)
        integer, intent(out) :: nI(nel), FlagsI
        integer, dimension(lenof_sign), intent(out) :: SignI
        real(dp) , intent(out) :: IterRDMStartI, AvSignI
        type(excit_gen_store_type), intent(inout), optional :: Store

        ! This is the iteration from which this determinant has been occupied.
        IterRDMStartI = CurrH_I(3)
        ! If there is nothing stored there yet, the first iteration is this one.
        IF(IterRDMStartI.eq.0.0_dp) IterRDMStartI = real(Iter, dp)
 
        ! This extracts everything, calculates the average sign, and adds in the 
        ! diagonal elements.
        call extract_bit_rep_rdm (iLutnI, IterRDMStartI, CurrH_I(2), 1.0_dp, nI, SignI, FlagsI, AvSignI)

    end subroutine extract_bit_rep_rdm_diag_norm

    subroutine extract_bit_rep_rdm_diag_hf_s_d(iLutnI, CurrH_I, nI, SignI, FlagsI, IterRDMStartI, AvSignI, Store)
! This routine is used when we're doing some truncated RDM calculation.    
! In this case, all the diagonal elements are calculated later on.
! Here we just need to calculate the average sign and extract the info as usual.
        use constants , only : dp, n_int, lenof_sign
        use SystemData , only : NEl
        use bit_reps , only : NIfTot, extract_bit_rep, extract_bit_rep_rdm
        use FciMCData , only : excit_gen_store_type, NCurrH
        use DetBitOps , only : TestClosedShellDet
        implicit none
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp) , intent(in) :: CurrH_I(NCurrH)
        integer, intent(out) :: nI(nel), FlagsI
        integer, dimension(lenof_sign), intent(out) :: SignI
        real(dp) , intent(out) :: IterRDMStartI, AvSignI
        type(excit_gen_store_type), intent(inout), optional :: Store

        call extract_bit_rep (iLutnI, nI, SignI, FlagsI, Store)

        if(tHF_Ref_Explicit) then
            IterRDMStartI = 0.0_dp
            AvSignI = real(SignI(1),dp)

        else
            IterRDMStartI = CurrH_I(3)
            IF(IterRDMStartI.eq.0.0_dp) IterRDMStartI = real(Iter, dp)
            AvSignI = ( ((real(Iter,dp) - IterRDMStartI) * CurrH_I(2)) &
                        + real(SignI(1),dp) ) / ( real(Iter,dp) - IterRDMStartI + 1.0_dp )
        endif

    end subroutine extract_bit_rep_rdm_diag_hf_s_d

    subroutine extract_bit_rep_rdm_diag_hphf(iLutnI, CurrH_I, nI, SignI, FlagsI, IterRDMStartI, AvSignI, Store)
! While extracting the orbitals from the bit representation of the determinant, we 
! simulaneously add in the contribution of each orbital to the diagonal elements of the RDMs.
! But we also need to add in the contribution from it's spin coupled partner if using HPHF, and possibly 
! connections between the two.
! This is not done in the best way at the moment, but is fine for now.
        use constants , only : dp, n_int, lenof_sign
        use SystemData , only : NEl
        use bit_reps , only : NIfTot, extract_bit_rep, extract_bit_rep_rdm
        use FciMCData , only : excit_gen_store_type, NCurrH
        use DetBitOps , only : TestClosedShellDet
        implicit none
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp) , intent(in) :: CurrH_I(NCurrH)
        integer, intent(out) :: nI(nel), FlagsI
        integer, dimension(lenof_sign), intent(out) :: SignI
        real(dp) , intent(out) :: IterRDMStartI, AvSignI
        type(excit_gen_store_type), intent(inout), optional :: Store

        IterRDMStartI = CurrH_I(3)
        IF(IterRDMStartI.eq.0.0_dp) IterRDMStartI = real(Iter, dp)

        if(.not.TestClosedShellDet(iLutnI)) then
            ! The 0.5 factor is because each coefficient needs to be divided by SQRT(2) (when it has 
            ! a spin coupled buddy that is not the same), but 
            ! we then square the coefficients - multiply the squared value by 0.5.
            call extract_bit_rep_rdm (iLutnI, IterRDMStartI, CurrH_I(2), 0.5_dp, nI, SignI, FlagsI, AvSignI)

            ! If HPHF is on, we need to consider I' as well as I.
            ! i.e. only one of the spin coupled determinants will be in the list, need to flip 
            ! the spins of all determinants to generate the other, and add in the C_I to this I' too.
            ! Only need to do this if I is open shell.
            ! TODO : Make this better
            call Add_StochRDM_Diag_HPHF(iLutnI, nI, AvSignI)
        else
            call extract_bit_rep_rdm (iLutnI, IterRDMStartI, CurrH_I(2), 1.0_dp, nI, SignI, FlagsI, AvSignI)
        endif

    end subroutine extract_bit_rep_rdm_diag_hphf


    subroutine Add_StochRDM_Diag_HPHF(iLutCurr,DetCurr,AvSignCurr)
! This is called when we run over all TotWalkers in CurrentDets.    
! It is only called if HF is being used and we've encountered an open shell det.
! The decoding routine would have added in the diagonal elements of the original HF pair,
! need to take care of the spin coupled one, and any connection between them.
        use FciMCData , only : HFDet
        use hphf_integrals , only : hphf_sign
        use HPHFRandExcitMod , only : FindExcitBitDetSym
        use DetBitOps , only : FindBitExcitLevel
        implicit none
        integer(kind=n_int), intent(in) :: iLutCurr(0:NIfTot)
        integer , intent(in) :: DetCurr(NEl)
        real(dp) , intent(in) :: AvSignCurr
        integer(kind=n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(NEl), SignFac, HPHFExcitLevel

! Add diagonal elements to reduced density matrices.

! C_X D_X = C_X / SQRT(2) [ D_I +/- D_I'] - for open shell dets, divide stored C_X by SQRT(2). 
! Add in I.
        call FindExcitBitDetSym(iLutCurr, SpinCoupDet)
        call decode_bit_det (nSpinCoup, SpinCoupDet)
        ! Find out if it's + or - in the above expression.                
        SignFac = hphf_sign(iLutCurr)

        call Fill_Diag_RDM(nSpinCoup, (real(SignFac,dp)*AvSignCurr)/SQRT(2.D0))

! For HPHF we're considering < D_I + D_I' | a_a+ a_b+ a_j a_i | D_I + D_I' >
! Not only do we have diagonal < D_I | a_a+ a_b+ a_j a_i | D_I > terms, but also cross terms
! < D_I | a_a+ a_b+ a_j a_i | D_I' > if D_I and D_I' can be connected by a single or double 
! excitation.
! Find excitation level between D_I and D_I' and add in the contribution if connected.
        HPHFExcitLevel = FindBitExcitLevel (iLutCurr, SpinCoupDet, 2)
        if(HPHFExcitLevel.le.2) & 
            call Add_RDM_From_IJ_Pair(DetCurr,nSpinCoup,&
                                        AvSignCurr/SQRT(2.D0), &
                                        (real(SignFac,dp)*AvSignCurr)/SQRT(2.D0),.true.)

    end subroutine Add_StochRDM_Diag_HPHF

! These Add_RDM_HFConnections routines take the current determinant and if it is 
! a single or double excitation of the HF, they explicitly add in the contribution to the RDM 
! from the current and hf determinant.

    subroutine Add_RDM_HFConnections_Norm(iLutJ,nJ,AvSignJ,walkExcitLevel)
! This is called when we run over all TotWalkers in CurrentDets.    
! It is called for each CurrentDet which is a single or double of the HF.
! It explicitly adds in the HF - S/D connection, as if the HF were D_i and 
! the single or double D_j.
! This is the standard full space RDM calc (No HPHF).
! In this case the diagonal elements wll already be taken care of.
        use constants , only : n_int, lenof_sign, dp
        use SystemData , only : NEl
        use bit_reps , only : NIfTot
        use FciMCData , only : HFDet, AvNoatHF
        use hphf_integrals , only : hphf_sign
        use HPHFRandExcitMod , only : FindExcitBitDetSym
        use DetBitOps , only : FindBitExcitLevel, TestClosedShellDet
        implicit none
        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        integer , intent(in) :: nJ(NEl)
        real(dp) , intent(in) :: AvSignJ
        integer , intent(in) :: walkExcitLevel
        integer(kind=n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(NEl), HPHFExcitLevel

! Quick check that the HF population is being calculated correctly.
        if(walkExcitLevel.eq.0) then
            if(AvSignJ.ne.AvNoatHF) then
                write(6,*) 'AvSignJ',AvSignJ
                write(6,*) 'AvNoatHF',AvNoatHF
                CALL Stop_All('PerformFCIMCycPar','Incorrect instantaneous HF population.')
            endif
        endif

! If we have a single or double, add in the connection to the HF, symmetrically.        
        if((walkExcitLevel.eq.1).or.(walkExcitLevel.eq.2)) &
            call Add_RDM_From_IJ_Pair(HFDet,nJ,AvNoatHF,AvSignJ,.true.)

    end subroutine Add_RDM_HFConnections_Norm


    subroutine Add_RDM_HFConnections_HPHF(iLutJ,nJ,AvSignJ,walkExcitLevel)
! This is called when we run over all TotWalkers in CurrentDets.    
! It is called for each CurrentDet which is a single or double of the HF.
! It adds in the HF - S/D connection.
! The diagonal elements will already have been taken care of by the extract routine.
        use constants , only : n_int, lenof_sign, dp
        use SystemData , only : NEl
        use bit_reps , only : NIfTot
        use FciMCData , only : HFDet, AvNoatHF
        use hphf_integrals , only : hphf_sign
        use HPHFRandExcitMod , only : FindExcitBitDetSym
        use DetBitOps , only : FindBitExcitLevel, TestClosedShellDet
        implicit none
        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        integer , intent(in) :: nJ(NEl)
        real(dp) , intent(in) :: AvSignJ
        integer , intent(in) :: walkExcitLevel
        integer(kind=n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(NEl), HPHFExcitLevel

        if(walkExcitLevel.eq.0) then
            if(AvSignJ.ne.AvNoatHF) then
                write(6,*) 'AvSignJ',AvSignJ
                write(6,*) 'AvNoatHF',AvNoatHF
                CALL Stop_All('PerformFCIMCycPar','Incorrect instantaneous HF population.')
            endif
        endif

! Now if the determinant is connected to the HF (i.e. single or double), add in the diagonal elements
! of this connection as well - symmetrically because no probabilities are involved.
        if((walkExcitLevel.eq.1).or.(walkExcitLevel.eq.2)) &
            call Fill_Spin_Coupled_RDM_v2(iLutRef,iLutJ,HFDet,nJ,&
                                            AvNoatHF,AvSignJ,.true.)

    end subroutine Add_RDM_HFConnections_HPHF

    subroutine Add_RDM_HFConnections_HF_S_D(iLutJ,nJ,AvSignJ,walkExcitLevel)
! This is called when we run over all TotWalkers in CurrentDets.    
! This finds all the connections to the HF when doing some sort of truncated RDM 
! calculation.
! Here, the diagonal elements will not have been added in by the extract routines.
! In the case of HF_Ref_Explicit, this routine does all the work.
        use constants , only : n_int, lenof_sign, dp
        use SystemData , only : NEl
        use bit_reps , only : NIfTot
        use FciMCData , only : HFDet, AvNoatHF
        use hphf_integrals , only : hphf_sign
        use HPHFRandExcitMod , only : FindExcitBitDetSym
        use DetBitOps , only : FindBitExcitLevel, TestClosedShellDet
        implicit none
        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        integer , intent(in) :: nJ(NEl)
        real(dp) , intent(in) :: AvSignJ
        integer , intent(in) :: walkExcitLevel
        integer(kind=n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(NEl), HPHFExcitLevel

! Add diagonal elements to reduced density matrices.

! If HF_S_D_Ref, we are only considering determinants connected to the HF, 
! doubles and singles (so theoretically up to quadruples).
! But for the diagonal elements - only consider doubles and singles (and HF).
            
        ! In all of these cases the HF is a diagonal element.
        if(walkExcitLevel.eq.0) then

            call Fill_Diag_RDM(nJ, AvSignJ)
            AccumRDMNorm = AccumRDMNorm + (AvSignJ * AvSignJ)
            AccumRDMNorm_Inst = AccumRDMNorm_Inst + (AvSignJ * AvSignJ)

            if(AvSignJ.ne.AvNoatHF) then
                write(6,*) 'AvSignJ',AvSignJ
                write(6,*) 'HF Sign',AvNoatHF
                call stop_all('Add_RDM_HFConnections_HF_S_D','HF population is incorrect.')
            endif

            ! The HF is always closed shell (at the moment), 
            ! so don't need to account for HPHF here.

        elseif(walkExcitLevel.le.2) then

            if(tHF_Ref_Explicit) then
                
                if(tHPHF) then

                    ! Now if the determinant is connected to the HF (i.e. single or double), 
                    ! add in the elements of this connection as well - symmetrically 
                    ! because no probabilities are involved.
                    call Fill_Spin_Coupled_RDM_v2(iLutRef,iLutJ,HFDet,nJ,&
                                AvNoatHF,AvSignJ,.false.)

                else

                    ! The singles and doubles are connected and explicitly calculated 
                    ! - but not symmetrically.
                    call Add_RDM_From_IJ_Pair(HFDet, nJ, AvNoatHF, &
                                                AvSignJ,.false.)

                endif

            else
                ! For the HF,S,D symmetric case, and the HF,S,D reference, the S and D
                ! are diagonal terms too.
                ! These options are not set up for HPHF.
                call Fill_Diag_RDM(nJ, AvSignJ)
                AccumRDMNorm = AccumRDMNorm + (AvSignJ * AvSignJ)
                AccumRDMNorm_Inst = AccumRDMNorm_Inst + (AvSignJ * AvSignJ)

                call Add_RDM_From_IJ_Pair(HFDet,nJ,AvNoatHF,AvSignJ,.true.)

            endif

        endif

    end subroutine Add_RDM_HFConnections_HF_S_D


    subroutine Add_RDM_HFConnections_null(iLutJ,nJ,AvSignJ,walkExcitLevel)
! This is called when we are not filling the density matrices.
        use constants , only : n_int, lenof_sign, dp
        use SystemData , only : NEl
        use bit_reps , only : NIfTot
        use FciMCData , only : HFDet, AvNoatHF
        use hphf_integrals , only : hphf_sign
        use HPHFRandExcitMod , only : FindExcitBitDetSym
        use DetBitOps , only : FindBitExcitLevel, TestClosedShellDet
        implicit none
        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        integer , intent(in) :: nJ(NEl)
        real(dp) , intent(in) :: AvSignJ
        integer , intent(in) :: walkExcitLevel
        integer(kind=n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(NEl), HPHFExcitLevel

    end subroutine Add_RDM_HFConnections_null

! This routine does the same as Fill_Spin_Coupled_RDM, but hopefully more efficiently!
! It takes to HPHF functions, and calculate what needs to be summed into the RDMs
    subroutine Fill_Spin_Coupled_RDM_v2(iLutnI,iLutnJ,nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)
        use systemData, only: tOddS_hphf
        use HPHFRandExcitMod, only: FindExcitBitDetSym,FindDetSpinSym
        use HPHF_Integrals , only : hphf_sign
        USE DetBitOps , only : FindBitExcitLevel, TestClosedShellDet
        implicit none
        integer(n_int), intent(in) :: iLutnI(0:NIfTot),iLutnJ(0:NIfTot)
        real(dp) , intent(in) :: realSignI, realSignJ
        integer, intent(in) :: nI(NEl),nJ(NEl)
        logical, intent(in) :: tFill_CiCj_Symm
        integer(n_int) :: iLutnI2(0:NIfTot)
        integer :: nI2(NEl),nJ2(NEl)
        real(dp) :: NewSignJ,NewSignI,PermSignJ,PermSignI
        integer :: I_J_ExcLevel,ICoup_J_ExcLevel
        character(*), parameter :: t_r='Fill_Spin_Coupled_RDM_v2'

        if(TestClosedShellDet(iLutnI)) then
            if(tOddS_HPHF) then
                call stop_all(t_r,"Should not be any closed shell determinants in high S states")
            endif

            if(TestClosedShellDet(iLutnJ)) then
                !Closed shell -> Closed shell - just as in determinant case
!                write(6,*) "CS -> CS "
                call Add_RDM_From_IJ_Pair(nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)
!                write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),realSignI,realSignJ
            else
                !Closed shell -> open shell.
!                write(6,*) "CS -> OS "
                call FindDetSpinSym(nJ,nJ2,NEl)
                NewSignJ = realSignJ/Root2
                call Add_RDM_From_IJ_Pair(nI,nJ,realSignI,NewSignJ,tFill_CiCj_Symm)
!                write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),realSignI,NewSignJ
                !What is the permutation between Di and Dj'
                NewSignJ = NewSignJ * hphf_sign(iLutnJ)
                call Add_RDM_From_IJ_Pair(nI,nJ2,realSignI,NewSignJ,tFill_CiCj_Symm)
!                write(6,"(A,4I4,2F12.6)") "I -> J' : ",nI(:),nJ2(:),realSignI,NewSignJ
            endif
        elseif(TestClosedShellDet(iLutnJ)) then
            !Open shell -> closed shell
!            write(6,*) "OS -> CS "
            call FindDetSpinSym(nI,nI2,NEl)
            NewSignI = realSignI/Root2
            call Add_RDM_From_IJ_Pair(nI,nJ,NewSignI,realSignJ,tFill_CiCj_Symm)
!            write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),NewSignI,realSignJ
            !What is the permutation between Di' and Dj 
            NewSignI = NewSignI * hphf_sign(iLutnI)
            call Add_RDM_From_IJ_Pair(nI2,nJ,NewSignI,realSignJ,tFill_CiCj_Symm)
!            write(6,"(A,4I4,2F12.6)") "I' -> J : ",nI2(:),nJ(:),NewSignI,realSignJ
        else
!            write(6,*) "OS -> OS "
            !Open shell -> open shell
            NewSignI = realSignI/Root2
            NewSignJ = realSignJ/Root2
            PermSignJ = NewSignJ * real(hphf_sign(iLutnJ),dp)
            PermSignI = NewSignI * real(hphf_sign(iLutnI),dp)
            call FindExcitBitDetSym(iLutnI, iLutnI2)
            call FindDetSpinSym(nI,nI2,NEl)
            call FindDetSpinSym(nJ,nJ2,NEl)
            I_J_ExcLevel = FindBitExcitLevel(iLutnI, iLutnJ,2)
            ICoup_J_ExcLevel = FindBitExcitLevel(iLutnI2,iLutnJ,2)
            if(I_J_ExcLevel.le.2) then
                call Add_RDM_From_IJ_Pair(nI,nJ,NewSignI,NewSignJ,tFill_CiCj_Symm)      !Di -> Dj
!                write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),NewSignI,NewSignJ
                call Add_RDM_From_IJ_Pair(nI2,nJ2,PermSignI,PermSignJ,tFill_CiCj_Symm)   !Di' -> Dj'  (both permuted sign)
!                write(6,"(A,4I4,2F12.6)") "I' -> J' : ",nI2(:),nJ2(:),PermSignI,PermSignJ
            endif
            if(ICoup_J_ExcLevel.le.2) then
                call Add_RDM_From_IJ_Pair(nI2,nJ,PermSignI,NewSignJ,tFill_CiCj_Symm)    !Di' -> Dj  (i permuted sign)
!                write(6,"(A,4I4,2F12.6)") "I' -> J : ",nI2(:),nJ(:),PermSignI,NewSignJ
                call Add_RDM_From_IJ_Pair(nI,nJ2,NewSignI,PermSignJ,tFill_CiCj_Symm)     !Di  -> Dj'  (j permuted sign)
!                write(6,"(A,4I4,2F12.6)") "I -> J' : ",nI(:),nJ2(:),NewSignI,PermSignJ
            endif
        endif

    end subroutine Fill_Spin_Coupled_RDM_v2

    subroutine Fill_Spin_Coupled_RDM(iLutnI,iLutnJ,nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)
!If the two HPHF determinants we're considering consist of I + I' and J + J', 
!where X' is the spin coupled (all spins flipped) version of X,
!then we have already considered the I -> J excitation.
!And if I and J are connected by a double excitation, tDoubleConnection is true and we have 
!also considered I' -> J'.
!But we need to also account for I -> J' and I' -> J.
        use HPHFRandExcitMod, only: FindExcitBitDetSym
        use HPHF_Integrals , only : hphf_sign
        USE DetBitOps , only : FindBitExcitLevel, TestClosedShellDet
        implicit none
        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot),iLutnJ(0:NIfTot)
        integer , intent(in) :: nI(NEl), nJ(NEl)
        real(dp) , intent(in) :: realSignI, realSignJ
        logical , intent(in) :: tFill_CiCj_Symm
        integer(kind=n_int) :: iLutnI2(0:NIfTot),iLutnJ2(0:NIfTot)
        integer :: Ex(2,2), SpinCoupI_J_ExcLevel, nI2(NEl), nJ2(NEl)
        integer :: SignFacI, SignFacJ, I_J_ExcLevel
        logical :: tParity
        real(dp) :: realSignFacI, realSignFacJ

!First we flip the spin of both determinants, and store I' and J'.
!Actually if I and J are related by a double excitation, we don't need J'.        

!First we flip the spin of I', and find out the excitation level between I' and J.
!If this is a double excitation, we don't actually need to find J' - we can just invert 
!the excitation matrix of the I' -> J transition.
!If this is anything above a double, we likewise don't need to find J', because I -> J' 
!will also have a 0 matrix element.

!        write(6,*) '***'
!        write(6,'(A5)',advance='no') 'nI'
!        do i = 1,4
!            write(6,'(I5)',advance='no') nI(i)
!        enddo
!        write(6,*) ''
!        write(6,'(A5)',advance='no') 'nJ'
!        do i = 1,4
!            write(6,'(I5)',advance='no') nJ(i)
!        enddo
!        write(6,*) ''
!        write(6,*) 'realSignI',realSignI
!        write(6,*) 'realSignJ',realSignJ

        I_J_ExcLevel = FindBitExcitLevel (iLutnI, iLutnJ, 2)

        if (.not. TestClosedShellDet(iLutnI)) then

            ! I is open shell, and so a spin coupled determinant I' exists.
!            write(6,*) 'I open shell'

            !Find I'.
            call FindExcitBitDetSym(iLutnI, iLutnI2)
            call decode_bit_det (nI2, iLutnI2)
            SignFacI = hphf_sign(iLutnI)
            realSignFacI = real(SignFacI,dp) / SQRT(2.0)

!            write(6,*) 'spin coupled nI'
!            do i = 1,4
!                write(6,'(I5)',advance='no') nI2(i)
!            enddo
!            write(6,*) ''

            !Find excitation level between I' and J - not necessarily the same as 
            !that between I and J.
            SpinCoupI_J_ExcLevel = FindBitExcitLevel (iLutnI2, iLutnJ, 2)

            IF( (.not.(I_J_ExcLevel.le.2)) .and. (.not.(SpinCoupI_J_ExcLevel.le.2)) ) &
                call Stop_All('Fill_Spin_Coupled_RDM','No spin combination are connected.')
                
            if ( .not. TestClosedShellDet(iLutnJ) ) then
                
                ! Both I and J are open shell, need all 4 combinations.

                !Find J'.
                call FindExcitBitDetSym(iLutnJ, iLutnJ2)
                call decode_bit_det (nJ2, iLutnJ2)
                SignFacJ = hphf_sign(iLutnJ)
                realSignFacJ = real(SignFacJ,dp) / SQRT(2.D0)

!                write(6,*) "OS -> OS "

                if (I_J_ExcLevel.le.2) then

                    ! I -> J.
                    call Add_RDM_From_IJ_Pair(nI,nJ,(realSignI/SQRT(2.D0)),&
                                                (realSignJ/SQRT(2.D0)),tFill_CiCj_Symm)
!                    write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),realSignI/SQRT(2.D0),realSignJ/SQRT(2.D0)
 
                    ! I' -> J'.
                    call Add_RDM_From_IJ_Pair(nI2,nJ2,(realSignFacI*realSignI),&
                                              (realSignFacJ*realSignJ),tFill_CiCj_Symm)
!                    write(6,"(A,4I4,2F12.6)") "I' -> J' : ",nI2(:),nJ2(:),(realSignFacI*realSignI),(realSignFacJ*realSignJ)

                endif

                if (SpinCoupI_J_ExcLevel.le.2) then

                    ! I' -> J.
                    call Add_RDM_From_IJ_Pair(nI2,nJ,(realSignFacI*realSignI),&
                                                (realSignJ/SQRT(2.D0)),tFill_CiCj_Symm)

!                    write(6,"(A,4I4,2F12.6)") "I' -> J : ",nI2(:),nJ(:),realSignFacI*realSignI,realSignJ/SQRT(2.D0)
                    ! I -> J'.
                    call Add_RDM_From_IJ_Pair(nI, nJ2,(realSignI/SQRT(2.D0)),&
                                                 (realSignFacJ*realSignJ),tFill_CiCj_Symm)

!                    write(6,"(A,4I4,2F12.6)") "I -> J' : ",nI(:),nJ2(:),realSignI/SQRT(2.D0),realSignFacJ*realSignJ

                endif

            else
                
                ! I is open shell, but J is not.
                ! Need I -> J and I' -> J.

!                write(6,*) "OS -> CS "
                ! I -> J.
                call Add_RDM_From_IJ_Pair(nI,nJ,(realSignI/SQRT(2.D0)),&
                                                        realSignJ,tFill_CiCj_Symm)
!                write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),realSignI/SQRT(2.D0),realSignJ
                ! I' -> J.
                call Add_RDM_From_IJ_Pair(nI2,nJ,(realSignFacI*realSignI),&
                                                        realSignJ,tFill_CiCj_Symm)
!                write(6,"(A,4I4,2F12.6)") "I' -> J : ",nI2(:),nJ(:),realSignFacI*realSignI,realSignJ

            endif

        elseif( .not. TestClosedShellDet(iLutnJ) ) then
            ! This is the case where I is closed shell, but J is not.
            ! Need I -> J and I -> J'. 

!            write(6,*) "CS -> OS "
            ! I -> J.
            if (I_J_ExcLevel.le.2) call Add_RDM_From_IJ_Pair(nI,nJ,realSignI,&
                                               (realSignJ/SQRT(2.D0)),tFill_CiCj_Symm)
!            write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),realSignI,realSignJ/SQRT(2.D0)

            ! Find J'.
            call FindExcitBitDetSym(iLutnJ, iLutnJ2)
            SignFacJ = hphf_sign(iLutnJ)
            realSignFacJ = real(SignFacJ,dp) / SQRT(2.D0)

            !Find excitation level between I and J'.
            SpinCoupI_J_ExcLevel = FindBitExcitLevel (iLutnI, iLutnJ2, 2)

            if (SpinCoupI_J_ExcLevel.le.2) then
                call decode_bit_det (nJ2, iLutnJ2)
                
                ! I -> J'.
                call Add_RDM_From_IJ_Pair(nI,nJ2,realSignI,&
                                            (realSignFacJ*realSignJ),tFill_CiCj_Symm)
!                write(6,"(A,4I4,2F12.6)") "I -> J' : ",nI(:),nJ2(:),realSignI,realSignFacJ*realSignJ

           endif

       elseif(I_J_ExcLevel.le.2) then

            ! I and J are both closed shell.

            ! Just I -> J.
            call Add_RDM_From_IJ_Pair(nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)
!            write(6,*) "CS -> CS "
!            write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),realSignI,realSignJ

       endif

    end subroutine Fill_Spin_Coupled_RDM

    subroutine Add_RDM_From_IJ_Pair(nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)
! This routine takes a pair of different determinants Di and Dj, and figures out which type 
! of elements need to be added in to the RDM.
        implicit none
        integer , intent(in) :: nI(NEl), nJ(NEl)
        real(dp) , intent(in) :: realSignI, realSignJ
        logical , intent(in) :: tFill_CiCj_Symm
        integer :: Ex(2,2),j
        logical :: tParity

        Ex(:,:) = 0
        Ex(1,1) = 2         ! Maximum excitation level - we know they are connected by 
                            ! a double or single.
        tParity = .false.

        call GetExcitation(nI,nJ,NEl,Ex,tParity)
! Ex(1,:) comes out as the orbital(s) excited from, i.e. i,j 
! Ex(2,:) comes out as the orbital(s) excited to, i.e. a,b.         

        IF(Ex(1,1).le.0) THEN
            ! Error.
            write(6,*) '*'
            WRITE(6,*) 'nI',nI
            WRITE(6,*) 'nJ',nJ
            WRITE(6,*) 'Ex(:,:)',Ex(1,1),Ex(1,2),Ex(2,1),Ex(2,2)
            WRITE(6,*) 'tParity',tParity
            WRITE(6,*) 'realSignI',realSignI
            WRITE(6,*) 'realSignJ',realSignJ
            write(6,*) '*'
            call flush(6)
            CALL Stop_All('Add_RDM_From_IJ_Pair',&
                    'Excitation level between pair not 1 or 2 as it should be.')
        ENDIF

        if((Ex(1,2).eq.0).and.(Ex(2,2).eq.0)) then

            ! Di and Dj are separated by a single excitation.
            ! Add in the contribution from this pair into the 1- and 2-RDM.
            call Fill_Sings_RDM(nI,Ex,tParity,realSignI,realSignJ,tFill_CiCj_Symm)
    
        elseif(RDMExcitLevel.ne.1) then

            ! Otherwise Di and Dj are connected by a double excitation.
            ! Add in this contribution to the 2-RDM (as long as we're calculating this obv).
            call Fill_Doubs_RDM(Ex,tParity,realSignI,realSignJ,tFill_CiCj_Symm)

        endif

    end subroutine Add_RDM_From_IJ_Pair

! EXPLICIT ROUTINES    
    SUBROUTINE Fill_ExplicitRDM_this_Iter(TotWalkers)
        USE FciMCData , only : CurrentDets,TotParts 
        USE bit_reps , only : encode_sign, extract_sign
        implicit none
        INTEGER(int64) , INTENT(IN) :: TotWalkers
        INTEGER(kind=n_int) :: iLutnI(0:NIfTot)
        INTEGER(int64) :: MaxTotWalkers,TotWalkIn(2),TotWalkOut(2)
        INTEGER :: i,error
        REAL(dp) :: TempTotParts, NormalisationTemp, Sum_Coeffs
        LOGICAL :: blank_det
        INTEGER, DIMENSION(lenof_sign) :: SignI, SignI2

! Run through the current determinants.
! Find the max number of determinants on a processor - all need to run through this 
! number so that the communication can be done at all stages.

        TotWalkIn(1)=TotWalkers
        TotWalkIn(2)=iProcIndex

        CALL MPIAllReduceDatatype(TotWalkIn,1,MPI_MAXLOC,MPI_2INTEGER,TotWalkOut)

        MaxTotWalkers=TotWalkOut(1)

! This little commented routine calculates the normalisation factor for the coefficients 
! at this iteration.
! It appears this is not necessary, and in reality it is fine to just add in 
! C_i * AllTotPartsCurr / AllTotParts , and then scale the density matrices at the end so 
! that their traces are correct.
! This also means the contributions are weighted by the number of walkers in the system 
! at the time.
!        Normalisation = 0.D0
!        NormalisationTemp = 0.D0
!        do i = 1, TotWalkers
!            call extract_sign (CurrentDets(:,I), SignI)
!            NormalisationTemp = NormalisationTemp + ( REAL(SignI(1)) * REAL(SignI(1)) )
!        enddo
!        WRITE(6,*) 'NormalisationTemp',NormalisationTemp

!        CALL MPIAllReduce(NormalisationTemp,MPI_SUM,Normalisation)

!        Normalisation = SQRT( Normalisation )
!        WRITE(6,*) 'Normalisation',Normalisation
!        Sum_Coeffs = SQRT( Normalisation )

!        do I = 1, TotWalkers
!            SignI2 = 0
!            call extract_sign (CurrentDets(:,I), SignI)
!            SignI2(1) = NINT(REAL(SignI(1)) * ( 1.D0 / Sum_Coeffs ))
!            call encode_sign(CurrentDets(:,I), SignI2)
!        enddo

        CALL set_timer(nElRDM_Time,30)

        do i=1,MaxTotWalkers

! But if the actual number of determinants on this processor is less than the number 
! we're running through, feed in 0 determinants and 0 sign.
            IF(i.gt.TotWalkers) THEN
                iLutnI(:)=0
                blank_det=.true.
            ELSE
                iLutnI(:)=CurrentDets(:,i)
                blank_det=.false.
            ENDIF

            CALL Add_ExplicitRDM_Contrib(iLutnI,blank_det)

        enddo
        CALL halt_timer(nElRDM_Time)

    END SUBROUTINE Fill_ExplicitRDM_this_Iter


    SUBROUTINE Fill_Hist_ExplicitRDM_this_Iter(TotWalkers)
        USE FciMCData , only : CurrentDets,TotParts 
        USE bit_reps , only : encode_sign, extract_sign
        USE DetCalcData , only : Det
        use hist_data, only: AllHistogram, Histogram
        use DetCalcData , only : FCIDets
        implicit none
        INTEGER(int64) , INTENT(IN) :: TotWalkers
        INTEGER(kind=n_int) :: iLutnI(0:NIfTot)
        INTEGER :: i,error
        REAL(dp) :: TempTotParts, NormalisationTemp, Sum_Coeffs
        LOGICAL :: blank_det
        INTEGER, DIMENSION(lenof_sign) :: TempSign

! Run through the current determinants.
! Find the max number of determinants on a processor - all need to run through this 
! number so that the communication can be done at all stages.

!        TotWalkIn(1)=TotWalkers
!        TotWalkIn(2)=iProcIndex
!
!        CALL MPIAllReduceDatatype(TotWalkIn,1,MPI_MAXLOC,MPI_2INTEGER,TotWalkOut)
!
!        MaxTotWalkers=TotWalkOut(1)

! This little commented routine calculates the normalisation factor for the coefficients 
! at this iteration.
! It appears this is not necessary, and in reality it is fine to just add in 
! C_i * AllTotPartsCurr / AllTotParts , and then scale the density matrices at the end so 
! that their traces are correct.
! This also means the contributions are weighted by the number of walkers in the system 
! at the time.
!        Normalisation = 0.D0
!        NormalisationTemp = 0.D0
!        do i = 1, TotWalkers
!            call extract_sign (CurrentDets(:,I), SignI)
!            NormalisationTemp = NormalisationTemp + ( REAL(SignI(1)) * REAL(SignI(1)) )
!        enddo
!        WRITE(6,*) 'NormalisationTemp',NormalisationTemp

!        CALL MPIAllReduce(NormalisationTemp,MPI_SUM,Normalisation)

!        Normalisation = SQRT( Normalisation )
!        WRITE(6,*) 'Normalisation',Normalisation
!        Sum_Coeffs = SQRT( Normalisation )

!        do I = 1, TotWalkers
!            SignI2 = 0
!            call extract_sign (CurrentDets(:,I), SignI)
!            SignI2(1) = NINT(REAL(SignI(1)) * ( 1.D0 / Sum_Coeffs ))
!            call encode_sign(CurrentDets(:,I), SignI2)
!        enddo

        CALL set_timer(nElRDM_Time,30)

        CALL MPISumAll(Histogram,AllHistogram)

        norm=0.D0
        if(iProcIndex.eq.0) then
            do i=1,Det
                norm=norm+AllHistogram(1,i)**2
            enddo
            norm=SQRT(norm)
        endif

        CALL MPISumAll_inplace(norm)
 
        do i=1,Det

! But if the actual number of determinants on this processor is less than the number 
! we're running through, feed in 0 determinants and 0 sign.
!            IF(AllHistogram(1,i).eq.0.0_dp) THEN
            IF(Histogram(1,i).eq.0.0_dp) THEN
                iLutnI(:)=0
                blank_det=.true.
            ELSE
                iLutnI(:)=FCIDets(:,i)
                blank_det=.false.
            ENDIF

            TempSign(1) = i
            call encode_sign(iLutnI,TempSign)

            CALL Add_Hist_ExplicitRDM_Contrib(iLutnI,blank_det)

        enddo
        CALL halt_timer(nElRDM_Time)

    END SUBROUTINE Fill_Hist_ExplicitRDM_this_Iter


    SUBROUTINE Add_ExplicitRDM_Contrib(iLutnI,blank_det)
! This is the general routine for taking a particular determinant in the spawned list, 
! D_i and adding it's contribution to the reduced density matrix.
        implicit none
        INTEGER(kind=n_int), INTENT(IN) :: iLutnI(0:NIfTot)
        LOGICAL, INTENT(IN) :: blank_det
        INTEGER :: i
        
! Set up excitation arrays.
! These are blocked according to the processor the excitation would be on if occupied.
! In each block, the first entry is the sign of determinant D_i and the second the bit 
! string of the determinant (these need to be sent along with the excitations).
! Each processor will have a different Di.

        Sing_ExcDjs(:,:)=0
        Sing_ExcList(:)=0
        Sing_ExcList(:) = Sing_InitExcSlots(:)

        do i=0,nProcessors-1
            Sing_ExcDjs(:,Sing_ExcList(i)) = iLutnI(:)
            Sing_ExcList(i) = Sing_ExcList(i)+1
        enddo
        IF(RDMExcitLevel.ne.1) THEN
            Doub_ExcDjs(:,:) = 0
            Doub_ExcList(:) = 0
            Doub_ExcList(:) = Doub_InitExcSlots(:)

            do i=0,nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(i)) = iLutnI(:)
                Doub_ExcList(i) = Doub_ExcList(i)+1
            enddo
        ENDIF

        IF(.not.blank_det) CALL GenExcDjs(iLutnI)
! Out of here we will get a filled ExcDjs array with all the single or double excitations 
! from Dj, this will be done for each proc. 

! We then need to send the excitations to the relevant processors.
        CALL SendProcExcDjs()
! This routine then calls SearchOccDets which takes each excitation and and binary 
! searches the occupied determinants for this.
! If found, we re-find the orbitals and parity involved in the excitation, and add the 
! c_i*c_j contributions to the corresponding matrix element.

    END SUBROUTINE Add_ExplicitRDM_Contrib

    SUBROUTINE Add_Hist_ExplicitRDM_Contrib(iLutnI,blank_det)
! This is the general routine for taking a particular determinant in the spawned list, 
! D_i and adding it's contribution to the reduced density matrix.
        implicit none
        INTEGER(kind=n_int), INTENT(IN) :: iLutnI(0:NIfTot)
        LOGICAL, INTENT(IN) :: blank_det
        INTEGER :: i
        
! Set up excitation arrays.
! These are blocked according to the processor the excitation would be on if occupied.
! In each block, the first entry is the sign of determinant D_i and the second the bit 
! string of the determinant (these need to be sent along with the excitations).
! Each processor will have a different Di.

        Sing_ExcDjs(:,:)=0
        Sing_ExcList(:)=0
        Sing_ExcList(:) = Sing_InitExcSlots(:)

        do i=0,nProcessors-1
            Sing_ExcDjs(:,Sing_ExcList(i)) = iLutnI(:)
            Sing_ExcList(i) = Sing_ExcList(i)+1
        enddo
        IF(RDMExcitLevel.ne.1) THEN
            Doub_ExcDjs(:,:) = 0
            Doub_ExcList(:) = 0
            Doub_ExcList(:) = Doub_InitExcSlots(:)

            do i=0,nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(i)) = iLutnI(:)
                Doub_ExcList(i) = Doub_ExcList(i)+1
            enddo
        ENDIF

        IF(.not.blank_det) CALL Gen_Hist_ExcDjs(iLutnI)
! Out of here we will get a filled ExcDjs array with all the single or double excitations 
! from Dj, this will be done for each proc. 

! We then need to send the excitations to the relevant processors.
        CALL Send_Hist_ProcExcDjs()
! This routine then calls SearchOccDets which takes each excitation and and binary 
! searches the occupied determinants for this.
! If found, we re-find the orbitals and parity involved in the excitation, and add the 
! c_i*c_j contributions to the corresponding matrix element.

    END SUBROUTINE Add_Hist_ExplicitRDM_Contrib


    SUBROUTINE GenExcDjs(iLutnI)
! This uses GenExcitations3 in symexcit3.F90 to generate all the possible either 
! single or double excitations from D_i, finds the processor they would be on if occupied, 
! and puts them in the SingExcDjs array according to that processor.
        USE DetBitOps , only : EncodeBitDet
        USE AnnihilationMod , only : DetermineDetNode
        USE SymExcit3 , only : GenExcitations3
        USE RotateOrbsData , only : SymLabelListInv
        USE bit_reps , only : extract_bit_rep
        implicit none
        INTEGER(kind=n_int) , INTENT(IN) :: iLutnI(0:NIfTot)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: SignDi, SignDi2
        INTEGER :: ExcitMat3(2,2),nI(NEl),nJ(NEl),Proc,FlagsDi,a,b,CountTemp
        LOGICAL :: tParity,tAllExcitFound

        call extract_bit_rep (iLutnI, nI, SignDi, FlagsDi)
! Unfortunately uses the decoded determinant - might want to look at this.        

        call Fill_Diag_RDM(nI,real(SignDi(1),dp))

!        CountTemp = 0

        ExcitMat3(:,:)=0
! Zeros in ExcitMat3 starts off at the first single excitation.        
        tAllExcitFound=.false.
! This becomes true when all the excitations have been found.        

        do while (.not.tAllExcitFound)
!                write(6,*) 'generating singles'
!                call flush(6)
            CALL GenExcitations3(nI,iLutnI,nJ,1,ExcitMat3(:,:),tParity,&
                                                        tAllExcitFound,.true.)            
! Passed out of here is the singly excited determinant, nJ.
! Information such as the orbitals involved in the excitation and the parity is also found 
! in this step, we are not currently storing this, and it is re-calculated later on 
! (after the determinants are passed to the relevant processor) - but the speed of sending 
! this information vs recalculating it will be tested.
! RDMExcitLevel is passed through, if this is 1, only singles are generated, 
! if it is 2 only doubles are found.

            IF(tAllExcitFound) EXIT

            iLutnJ(:)=0
            CALL EncodeBitDet(nJ,iLutnJ)

            Proc = DetermineDetNode(nJ,0)   
            !This will return a value between 0 -> nProcessors-1
            Sing_ExcDjs(:,Sing_ExcList(Proc)) = iLutnJ(:)
            Sing_ExcList(Proc) = Sing_ExcList(Proc)+1
!                CountTemp = CountTemp + 1

!                IF((nI(1).eq.5).and.(nI(2).eq.6)) WRITE(6,*) CountTemp,': nJ',nJ

! Want a quick test to see if arrays are getting full.            
            IF(Sing_ExcList(Proc).gt.NINT(OneEl_Gap*(Proc+1))) THEN
                WRITE(6,*) 'Proc',Proc
                WRITE(6,*) 'Sing_ExcList',Sing_ExcList
                WRITE(6,*) 'No. spaces for each proc',NINT(OneEl_Gap)
                CALL Stop_All('GenExcDjs',&
                            'Too many excitations for space available.')
            ENDIF
        enddo

        IF(RDMExcitLevel.ne.1) THEN            

            ExcitMat3(:,:)=0
! Zeros in ExcitMat3 starts off at the first single excitation.        
            tAllExcitFound=.false.
! This becomes true when all the excitations have been found.        

!            WRITE(6,*) 'Generating from nI',nI
!            WRITE(6,*) 'bit rep',iLutnI

            do while (.not.tAllExcitFound)
!                write(6,*) 'generating doubles'
!                call flush(6)
                CALL GenExcitations3(nI,iLutnI,nJ,2,ExcitMat3(:,:),tParity,&
                                                            tAllExcitFound,.true.)            
! Passed out of here is the doubly excited determinant, nJ.
! Information such as the orbitals involved in the excitation and the parity is 
! also found in this step, we are not currently storing this, and it is re-calculated 
! later on (after the determinants are passed to the relevant processor) - but the speed 
! of sending this information vs recalculating it will be tested.
! RDMExcitLevel is passed through, if this is 1, only singles are generated, 
! if it is 2 only doubles are found.

                IF(tAllExcitFound) EXIT

                iLutnJ(:)=0
                CALL EncodeBitDet(nJ,iLutnJ)

                Proc = DetermineDetNode(nJ,0)   
                !This will return a value between 0 -> nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(Proc)) = iLutnJ(:)
                ! All the double excitations from this particular nI are being 
                ! stored in Doub_ExcDjs.

                Doub_ExcList(Proc) = Doub_ExcList(Proc)+1

!                IF((nI(1).eq.5).and.(nI(2).eq.6)) WRITE(6,*) CountTemp,': nJ',nJ
!                IF((nI(1).eq.5).and.(nI(2).eq.6)) WRITE(6,*) 'Ind', &
!                                                (( (nJ(2)-2) * (nJ(2)-1) ) / 2 ) + nJ(1)

! Want a quick test to see if arrays are getting full.            
                IF(Doub_ExcList(Proc).gt.NINT(TwoEl_Gap*(Proc+1))) THEN
                    WRITE(6,*) 'Proc',Proc
                    WRITE(6,*) 'Doub_ExcList',Doub_ExcList
                    WRITE(6,*) 'No. spaces for each proc',NINT(TwoEl_Gap)
                    CALL Stop_All('GenExcDjs','Too many excitations for space available.')
                ENDIF
            enddo
        ENDIF

!        IF((nI(1).eq.5).and.(nI(2).eq.6)) CALL Stop_All('','')

    END SUBROUTINE GenExcDjs


    SUBROUTINE Gen_Hist_ExcDjs(iLutnI)
! This uses GenExcitations3 in symexcit3.F90 to generate all the possible either 
! single or double excitations from D_i, finds the processor they would be on if occupied, 
! and puts them in the SingExcDjs array according to that processor.
        USE DetBitOps , only : EncodeBitDet
        USE AnnihilationMod , only : DetermineDetNode
        USE SymExcit3 , only : GenExcitations3
        USE RotateOrbsData , only : SymLabelListInv
        USE bit_reps , only : extract_bit_rep
        use hist_data, only: AllHistogram
        implicit none
        INTEGER(kind=n_int) , INTENT(IN) :: iLutnI(0:NIfTot)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: HistPos
        INTEGER :: ExcitMat3(2,2),nI(NEl),nJ(NEl),Proc,FlagsDi,a,b,CountTemp
        LOGICAL :: tParity,tAllExcitFound
        real(dp) :: realSignDi

        call extract_bit_rep (iLutnI, nI, HistPos, FlagsDi)
! Unfortunately uses the decoded determinant - might want to look at this.        

        realSignDi = AllHistogram(1,HistPos(1))/norm
        
        call Fill_Diag_RDM(nI,realSignDi)

!        CountTemp = 0

        ExcitMat3(:,:)=0
! Zeros in ExcitMat3 starts off at the first single excitation.        
        tAllExcitFound=.false.
! This becomes true when all the excitations have been found.        

        do while (.not.tAllExcitFound)
!                write(6,*) 'generating singles'
!                call flush(6)
            CALL GenExcitations3(nI,iLutnI,nJ,1,ExcitMat3(:,:),tParity,&
                                                        tAllExcitFound,.true.)            
! Passed out of here is the singly excited determinant, nJ.
! Information such as the orbitals involved in the excitation and the parity is also found 
! in this step, we are not currently storing this, and it is re-calculated later on 
! (after the determinants are passed to the relevant processor) - but the speed of sending 
! this information vs recalculating it will be tested.
! RDMExcitLevel is passed through, if this is 1, only singles are generated, 
! if it is 2 only doubles are found.

            IF(tAllExcitFound) EXIT

            iLutnJ(:)=0
            CALL EncodeBitDet(nJ,iLutnJ)

            Proc = DetermineDetNode(nJ,0)   
            !This will return a value between 0 -> nProcessors-1
            Sing_ExcDjs(:,Sing_ExcList(Proc)) = iLutnJ(:)
            Sing_ExcList(Proc) = Sing_ExcList(Proc)+1
!                CountTemp = CountTemp + 1

!                IF((nI(1).eq.5).and.(nI(2).eq.6)) WRITE(6,*) CountTemp,': nJ',nJ

! Want a quick test to see if arrays are getting full.            
            IF(Sing_ExcList(Proc).gt.NINT(OneEl_Gap*(Proc+1))) THEN
                WRITE(6,*) 'Proc',Proc
                WRITE(6,*) 'Sing_ExcList',Sing_ExcList
                WRITE(6,*) 'No. spaces for each proc',NINT(OneEl_Gap)
                CALL Stop_All('GenExcDjs',&
                            'Too many excitations for space available.')
            ENDIF
        enddo

        IF(RDMExcitLevel.ne.1) THEN            

            ExcitMat3(:,:)=0
! Zeros in ExcitMat3 starts off at the first single excitation.        
            tAllExcitFound=.false.
! This becomes true when all the excitations have been found.        

!            WRITE(6,*) 'Generating from nI',nI
!            WRITE(6,*) 'bit rep',iLutnI

            do while (.not.tAllExcitFound)
!                write(6,*) 'generating doubles'
!                call flush(6)
                CALL GenExcitations3(nI,iLutnI,nJ,2,ExcitMat3(:,:),tParity,&
                                                            tAllExcitFound,.true.)            
! Passed out of here is the doubly excited determinant, nJ.
! Information such as the orbitals involved in the excitation and the parity is 
! also found in this step, we are not currently storing this, and it is re-calculated 
! later on (after the determinants are passed to the relevant processor) - but the speed 
! of sending this information vs recalculating it will be tested.
! RDMExcitLevel is passed through, if this is 1, only singles are generated, 
! if it is 2 only doubles are found.

                IF(tAllExcitFound) EXIT

                iLutnJ(:)=0
                CALL EncodeBitDet(nJ,iLutnJ)

                Proc = DetermineDetNode(nJ,0)   
                !This will return a value between 0 -> nProcessors-1
                Doub_ExcDjs(:,Doub_ExcList(Proc)) = iLutnJ(:)
                ! All the double excitations from this particular nI are being 
                ! stored in Doub_ExcDjs.

                Doub_ExcList(Proc) = Doub_ExcList(Proc)+1

!                IF((nI(1).eq.5).and.(nI(2).eq.6)) WRITE(6,*) CountTemp,': nJ',nJ
!                IF((nI(1).eq.5).and.(nI(2).eq.6)) WRITE(6,*) 'Ind', &
!                                                (( (nJ(2)-2) * (nJ(2)-1) ) / 2 ) + nJ(1)

! Want a quick test to see if arrays are getting full.            
                IF(Doub_ExcList(Proc).gt.NINT(TwoEl_Gap*(Proc+1))) THEN
                    WRITE(6,*) 'Proc',Proc
                    WRITE(6,*) 'Doub_ExcList',Doub_ExcList
                    WRITE(6,*) 'No. spaces for each proc',NINT(TwoEl_Gap)
                    CALL Stop_All('GenExcDjs','Too many excitations for space available.')
                ENDIF
            enddo
        ENDIF

!        IF((nI(1).eq.5).and.(nI(2).eq.6)) CALL Stop_All('','')

    END SUBROUTINE Gen_Hist_ExcDjs


    SUBROUTINE SendProcExcDjs()
! In this routine the excitations are sent to the relevant processors.
! Sent with them will be the Di they were excited from and its sign.
! Each processor will receive nProcessor number of lists with different Di determinants.
! The original Di's will (I think) still be in the original InitSingExcSlots positions.
! This follows the directannihilation algorithm closely.
        implicit none
        INTEGER :: i,j,sendcounts(nProcessors),disps(nProcessors)
        INTEGER :: sing_recvcounts(nProcessors)
        INTEGER :: sing_recvdisps(nProcessors),error,MaxSendIndex,MaxIndex
        INTEGER :: doub_recvcounts(nProcessors),doub_recvdisps(nProcessors)

        do i=0,nProcessors-1
            sendcounts(i+1)=Sing_ExcList(i)-(NINT(OneEl_Gap*i)+1)
! Sendcounts is the number of singly excited determinants we want to send for 
! each processor (but goes from 1, not 0).            
            disps(i+1)=NINT(OneEl_Gap*i)
! and I think disps is the first slot for each processor - 1.            
        enddo

        MaxSendIndex=Sing_ExcList(nProcessors-1)-1

! We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
        sing_recvcounts(1:nProcessors)=0
        CALL MPIAlltoAll(sendcounts,1,sing_recvcounts,1,error)
! I think recvcounts(i) is the number of determinants sent from processor i.        

! We can now get recvdisps from recvcounts, since we want the data to be 
! contiguous after the move.
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
        CALL MPIAlltoAllv(Sing_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,&
                            Sing_ExcDjs2,sing_recvcounts,sing_recvdisps,error)
#else
        Sing_ExcDjs2(0:NIfTot,1:MaxIndex)=Sing_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif

        CALL Sing_SearchOccDets(sing_recvcounts,sing_recvdisps)


        IF(RDMExcitLevel.ne.1) THEN
            do i=0,nProcessors-1
                sendcounts(i+1)=Doub_ExcList(i)-(NINT(TwoEl_Gap*i)+1)
! Sendcounts is the number of singly excited determinants we want to send for 
! each processor (but goes from 1, not 0).            
                disps(i+1)=NINT(TwoEl_Gap*i)
! and I think disps is the first slot for each processor - 1.            
            enddo

            MaxSendIndex = Doub_ExcList(nProcessors-1)-1

! We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
            doub_recvcounts(1:nProcessors)=0
            CALL MPIAlltoAll(sendcounts,1,doub_recvcounts,1,error)
! I think recvcounts(i) is the number of determinants sent from processor i.        

! We can now get recvdisps from recvcounts, since we want the data to be contiguous 
! after the move.
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
            CALL MPIAlltoAllv(Doub_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,&
                                    Doub_ExcDjs2,doub_recvcounts,doub_recvdisps,error)
#else
            Doub_ExcDjs2(0:NIfTot,1:MaxIndex)=Doub_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif
            CALL Doub_SearchOccDets(doub_recvcounts,doub_recvdisps)

        ENDIF

        
    END SUBROUTINE SendProcExcDjs

    SUBROUTINE Send_Hist_ProcExcDjs()
! In this routine the excitations are sent to the relevant processors.
! Sent with them will be the Di they were excited from and its sign.
! Each processor will receive nProcessor number of lists with different Di determinants.
! The original Di's will (I think) still be in the original InitSingExcSlots positions.
! This follows the directannihilation algorithm closely.
        implicit none
        INTEGER :: i,j,sendcounts(nProcessors),disps(nProcessors)
        INTEGER :: sing_recvcounts(nProcessors)
        INTEGER :: sing_recvdisps(nProcessors),error,MaxSendIndex,MaxIndex
        INTEGER :: doub_recvcounts(nProcessors),doub_recvdisps(nProcessors)

        do i=0,nProcessors-1
            sendcounts(i+1)=Sing_ExcList(i)-(NINT(OneEl_Gap*i)+1)
! Sendcounts is the number of singly excited determinants we want to send for 
! each processor (but goes from 1, not 0).            
            disps(i+1)=NINT(OneEl_Gap*i)
! and I think disps is the first slot for each processor - 1.            
        enddo

        MaxSendIndex=Sing_ExcList(nProcessors-1)-1

! We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
        sing_recvcounts(1:nProcessors)=0
        CALL MPIAlltoAll(sendcounts,1,sing_recvcounts,1,error)
! I think recvcounts(i) is the number of determinants sent from processor i.        

! We can now get recvdisps from recvcounts, since we want the data to be 
! contiguous after the move.
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
        CALL MPIAlltoAllv(Sing_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,&
                            Sing_ExcDjs2,sing_recvcounts,sing_recvdisps,error)
#else
        Sing_ExcDjs2(0:NIfTot,1:MaxIndex)=Sing_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif

        CALL Sing_Hist_SearchOccDets(sing_recvcounts,sing_recvdisps)


        IF(RDMExcitLevel.ne.1) THEN            
            do i=0,nProcessors-1
                sendcounts(i+1)=Doub_ExcList(i)-(NINT(TwoEl_Gap*i)+1)
! Sendcounts is the number of singly excited determinants we want to send for 
! each processor (but goes from 1, not 0).            
                disps(i+1)=NINT(TwoEl_Gap*i)
! and I think disps is the first slot for each processor - 1.            
            enddo

            MaxSendIndex = Doub_ExcList(nProcessors-1)-1

! We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
            doub_recvcounts(1:nProcessors)=0
            CALL MPIAlltoAll(sendcounts,1,doub_recvcounts,1,error)
! I think recvcounts(i) is the number of determinants sent from processor i.        

! We can now get recvdisps from recvcounts, since we want the data to be contiguous 
! after the move.
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
            CALL MPIAlltoAllv(Doub_ExcDjs(:,1:MaxSendIndex),sendcounts,disps,&
                                    Doub_ExcDjs2,doub_recvcounts,doub_recvdisps,error)
#else
            Doub_ExcDjs2(0:NIfTot,1:MaxIndex)=Doub_ExcDjs(0:NIfTot,1:MaxSendIndex)
#endif
            CALL Doub_Hist_SearchOccDets(doub_recvcounts,doub_recvdisps)

        ENDIF

        
    END SUBROUTINE Send_Hist_ProcExcDjs


    SUBROUTINE Sing_SearchOccDets(recvcounts,recvdisps)
! We now have arrays SingExcDjs2 which contain all the single excitations from 
! each processor.
! These number sent from processor i is recvcounts(i), and the first 2 have information 
! about the determinant Di from which the Dj's are single excitations (and it's sign).
        USE AnnihilationMod , only : BinSearchParts
        USE FciMCData , only : TotWalkers,CurrentDets
        USE RotateOrbsData , only : SymLabelListInv
        USE bit_reps , only : extract_bit_rep
        implicit none
        INTEGER, INTENT(IN) :: recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: SignDi,SignDj, SignDi2,SignDj2
        INTEGER :: PartInd
        INTEGER :: i,j,NoDets,StartDets
        INTEGER :: nI(NEl),nJ(NEl),Ex(2,2),FlagsDi,FlagsDj
        LOGICAL :: tDetFound,tParity
        REAL(dp) :: realSignDi, realSignDj

! Take each Dj, and binary search the CurrentDets to see if it is occupied.

        do i=1,nProcessors
! Doing determinants from each processor separately because each has a different 
! D_i it was excited from.

            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            IF(NoDets.gt.1) THEN
                call extract_bit_rep (Sing_ExcDjs2(:,StartDets), nI, SignDi, FlagsDi)

                realSignDi = real(SignDi(1))

                do j=StartDets+1,(NoDets+StartDets-1)
! D_i is in the first spot - start from the second.                
                
                    iLutnJ(:)=Sing_ExcDjs2(:,j)
! This binary searches CurrentDets between 1 and TotWalkers for determinant iLutnJ.
! If found, tDetFound will be true, and PartInd the index in CurrentDets where the 
! determinant is.
                    CALL BinSearchParts(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)
                    IF(tDetFound) THEN
! Determinant occupied; add c_i*c_j to the relevant element of nElRDM.                    
! Need to first find the orbitals involved in the excitation from D_i -> D_j and the parity.

                        Ex(:,:)=0
! Ex(1,1) goes in as the max number of excitations - we know this is an excitation of 
! level RDMExcitLevel. 
                        Ex(1,1)=1
                        tParity = .false.

                        call extract_bit_rep (CurrentDets(:,PartInd), nJ, SignDj, FlagsDj)

                        realSignDj = real(SignDj(1))

! Ex(1,:) comes out as the orbital(s) excited from, Ex(2,:) comes out as the orbital(s) 
! excited to.    
                        CALL GetExcitation(nI,nJ,NEl,Ex,tParity)

                        IF(Ex(1,1).le.0) CALL Stop_All('Sing_SearchOccDets',&
                                            'nJ is not the correct excitation of nI.')

                        call Fill_Sings_RDM(nI,Ex,tParity,realSignDi,realSignDj,.true.)

! No normalisation factor just yet - possibly need to revise.                    
                    ENDIF

                enddo
            ENDIF
        enddo
      
    END SUBROUTINE Sing_SearchOccDets


    SUBROUTINE Doub_SearchOccDets(recvcounts,recvdisps)
! We now have arrays SingExcDjs2 which contain all the single excitations 
! from each processor.
! These number sent from processor i is recvcounts(i), and the first 2 have information 
! about the determinant Di from which the Dj's are single excitations (and it's sign).
        USE AnnihilationMod , only : BinSearchParts
        USE FciMCData , only : TotWalkers,CurrentDets
        USE RotateOrbsData , only : SymLabelListInv
        USE bit_reps , only : extract_bit_rep
        implicit none
        INTEGER, INTENT(IN) :: recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: SignDi,SignDj, SignDi2, SignDj2
        INTEGER :: PartInd
        INTEGER :: i,j,NoDets,StartDets
        INTEGER :: nI(NEl),nJ(NEl),Ex(2,2),FlagsDi,FlagsDj
        LOGICAL :: tDetFound,tParity
        REAL(dp) :: realSignDi,realSignDj

! Take each Dj, and binary search the CurrentDets to see if it is occupied.
    
!        WRITE(6,*) 'Searching for generated nJs in occupied list'

        do i=1,nProcessors
! Doing determinants from each processor separately because each has a different D_i 
! it was excited from.

            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            IF(NoDets.gt.1) THEN
                call extract_bit_rep (Doub_ExcDjs2(:,StartDets), nI, SignDi, FlagsDi)

                realSignDi = real(SignDi(1))

                do j=StartDets+1,(NoDets+StartDets-1)
! D_i is in the first spot - start from the second.                
                
                    iLutnJ(:)=Doub_ExcDjs2(:,j)

! This binary searches CurrentDets between 1 and TotWalkers for determinant iLutnJ.
! If found, tDetFound will be true, and PartInd the index in CurrentDets where the 
! determinant is.
                    CALL BinSearchParts(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

                    IF(tDetFound) THEN
! Determinant occupied; add c_i*c_j to the relevant element of nElRDM.                    
! Need to first find the orbitals involved in the excitation from D_i -> D_j and 
! the parity.

                        Ex(:,:)=0
! Ex(1,1) goes in as the max number of excitations - we know this is an excitation of 
! level RDMExcitLevel. 
                        Ex(1,1)=2
                        tParity = .false.

                        call extract_bit_rep (CurrentDets(:,PartInd), nJ, SignDj, FlagsDj)

                        realSignDj = real(SignDj(1))

! Ex(1,:) comes out as the orbital(s) excited from, Ex(2,:) comes out as the orbital(s) 
! excited to. 
                        CALL GetExcitation(nI,nJ,NEl,Ex,tParity)

                        IF(Ex(1,1).le.0) CALL Stop_All('SearchOccDets',&
                                            'nJ is not the correct excitation of nI.')

                        call Fill_Doubs_RDM(Ex,tParity,realSignDi,realSignDj,.true.)
                        
                        
                    ENDIF
                enddo
            ENDIF
        enddo
      
    END SUBROUTINE Doub_SearchOccDets


    SUBROUTINE Sing_Hist_SearchOccDets(recvcounts,recvdisps)
! We now have arrays SingExcDjs2 which contain all the single excitations from 
! each processor.
! These number sent from processor i is recvcounts(i), and the first 2 have information 
! about the determinant Di from which the Dj's are single excitations (and it's sign).
        USE AnnihilationMod , only : BinSearchParts
        USE FciMCData , only : TotWalkers,CurrentDets
        USE RotateOrbsData , only : SymLabelListInv
        USE bit_reps , only : extract_bit_rep
        USE FciMCData , only : iluthf
        use DetBitOps , only : FindBitExcitLevel
        use hist_data, only: AllHistogram
        use hist , only : find_hist_coeff_explicit
        implicit none
        INTEGER, INTENT(IN) :: recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: HistPos
        INTEGER :: PartInd, ExcitLevel
        INTEGER :: i,j,NoDets,StartDets
        INTEGER :: nI(NEl),nJ(NEl),Ex(2,2),FlagsDi,FlagsDj
        LOGICAL :: tDetFound,tParity
        REAL(dp) :: realSignDi, realSignDj

! Take each Dj, and binary search the CurrentDets to see if it is occupied.

        do i=1,nProcessors
! Doing determinants from each processor separately because each has a different 
! D_i it was excited from.

            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            IF(NoDets.gt.1) THEN
                call extract_bit_rep (Sing_ExcDjs2(:,StartDets), nI, HistPos, FlagsDi)

                realSignDi = AllHistogram(1,HistPos(1))/norm

                do j=StartDets+1,(NoDets+StartDets-1)
! D_i is in the first spot - start from the second.                
                
                    iLutnJ(:)=Sing_ExcDjs2(:,j)
! This binary searches CurrentDets between 1 and TotWalkers for determinant iLutnJ.
! If found, tDetFound will be true, and PartInd the index in CurrentDets where the 
! determinant is.
!                    CALL BinSearchParts(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

                    ExcitLevel = FindBitExcitLevel (iLutHF, iLutnJ, NEl)
                    call find_hist_coeff_explicit (iLutnJ, ExcitLevel,PartInd,tDetFound)

                    IF(tDetFound) THEN
! Determinant occupied; add c_i*c_j to the relevant element of nElRDM.                    
! Need to first find the orbitals involved in the excitation from D_i -> D_j and the parity.

                        Ex(:,:)=0
! Ex(1,1) goes in as the max number of excitations - we know this is an excitation of 
! level RDMExcitLevel. 
                        Ex(1,1)=1
                        tParity = .false.

!                        call extract_bit_rep (CurrentDets(:,PartInd), nJ, SignDj, FlagsDj)
!                        realSignDj = real(SignDj(1))

                        call decode_bit_det(nJ,iLutnJ)

                        realSignDj = AllHistogram(1,PartInd)/norm

! Ex(1,:) comes out as the orbital(s) excited from, Ex(2,:) comes out as the orbital(s) 
! excited to.    
                        CALL GetExcitation(nI,nJ,NEl,Ex,tParity)

                        IF(Ex(1,1).le.0) CALL Stop_All('Sing_SearchOccDets',&
                                            'nJ is not the correct excitation of nI.')

                        call Fill_Sings_RDM(nI,Ex,tParity,realSignDi,realSignDj,.true.)

! No normalisation factor just yet - possibly need to revise.                    
                    ENDIF

                enddo
            ENDIF
        enddo
      
    END SUBROUTINE Sing_Hist_SearchOccDets


    SUBROUTINE Doub_Hist_SearchOccDets(recvcounts,recvdisps)
! We now have arrays SingExcDjs2 which contain all the single excitations 
! from each processor.
! These number sent from processor i is recvcounts(i), and the first 2 have information 
! about the determinant Di from which the Dj's are single excitations (and it's sign).
        USE AnnihilationMod , only : BinSearchParts
        USE FciMCData , only : TotWalkers,CurrentDets
        USE RotateOrbsData , only : SymLabelListInv
        USE bit_reps , only : extract_bit_rep
        USE FciMCData , only : iluthf
        use DetBitOps , only : FindBitExcitLevel
        use hist_data, only: AllHistogram
        use hist , only : find_hist_coeff_explicit
        implicit none
        INTEGER, INTENT(IN) :: recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: HistPos
        INTEGER :: PartInd, ExcitLevel
        INTEGER :: i,j,NoDets,StartDets
        INTEGER :: nI(NEl),nJ(NEl),Ex(2,2),FlagsDi,FlagsDj
        LOGICAL :: tDetFound,tParity
        REAL(dp) :: realSignDi,realSignDj

! Take each Dj, and binary search the CurrentDets to see if it is occupied.
    
!        WRITE(6,*) 'Searching for generated nJs in occupied list'

        do i=1,nProcessors
! Doing determinants from each processor separately because each has a different D_i 
! it was excited from.

            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            IF(NoDets.gt.1) THEN
                call extract_bit_rep (Doub_ExcDjs2(:,StartDets), nI, HistPos, FlagsDi)

                realSignDi = AllHistogram(1,HistPos(1))/norm

                do j=StartDets+1,(NoDets+StartDets-1)
! D_i is in the first spot - start from the second.                
                
                    iLutnJ(:)=Doub_ExcDjs2(:,j)

! This binary searches CurrentDets between 1 and TotWalkers for determinant iLutnJ.
! If found, tDetFound will be true, and PartInd the index in CurrentDets where the 
! determinant is.
!                    CALL BinSearchParts(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

                    ExcitLevel = FindBitExcitLevel (iLutHF, iLutnJ, NEl)
                    call find_hist_coeff_explicit (iLutnJ, ExcitLevel, PartInd, tDetFound)

                    IF(tDetFound) THEN
! Determinant occupied; add c_i*c_j to the relevant element of nElRDM.                    
! Need to first find the orbitals involved in the excitation from D_i -> D_j and 
! the parity.

                        Ex(:,:)=0
! Ex(1,1) goes in as the max number of excitations - we know this is an excitation of 
! level RDMExcitLevel. 
                        Ex(1,1)=2
                        tParity = .false.

!                        call extract_bit_rep (CurrentDets(:,PartInd), nJ, SignDj, FlagsDj)
!                        realSignDj = real(SignDj(1))

                        call decode_bit_det(nJ,iLutnJ)
                        realSignDj = AllHistogram(1,PartInd)/norm

! Ex(1,:) comes out as the orbital(s) excited from, Ex(2,:) comes out as the orbital(s) 
! excited to. 
                        CALL GetExcitation(nI,nJ,NEl,Ex,tParity)

                        IF(Ex(1,1).le.0) CALL Stop_All('SearchOccDets',&
                                            'nJ is not the correct excitation of nI.')

                        call Fill_Doubs_RDM(Ex,tParity,realSignDi,realSignDj,.true.)
                        
                        
                    ENDIF
                enddo
            ENDIF
        enddo
      
    END SUBROUTINE Doub_Hist_SearchOccDets


! THESE NEXT ROUTINES ARE GENERAL TO BOTH STOCHASTIC AND EXPLICIT    

    subroutine Fill_Diag_RDM(nI,realSignDi)
! Fill diagonal elements of 1- and 2-RDM.
! These are < Di | a_i+ a_i | Di > and < Di | a_i+ a_j+ a_j a_i | Di >.
        USE UMatCache , only : GTID
        implicit none
        integer , intent(in) :: nI(NEl)
        real(dp) , intent(in) :: realSignDi
        integer :: i, j, iSpat, jSpat, Ind, iInd

! Need to add in the diagonal elements.
        
!        WRITE(6,*) realSignDi

        if(RDMExcitLevel.eq.1) then
            do i=1,NEl
                if(tStoreSpinOrbs) then
                    iInd = SymLabelListInv(nI(i))
                else
                    ! SymLabelListInv will be in spat orbitals too.
                    iInd = SymLabelListInv(gtID(nI(i)))
                endif
                NatOrbMat(iInd,iInd) = NatOrbMat(iInd,iInd) &
                                          + ( realSignDi * realSignDi ) 
            enddo
        else
            ! Only calculating 2-RDM.
            do i=1,NEl - 1
                iSpat = gtID(nI(i))

                ! Orbitals in nI ordered lowest to highest so nI(j) > nI(i), 
                ! and jSpat >= iSpat (can only be equal if different spin).
                do j=i+1,NEl
                    jSpat = gtID(nI(j))

                    ! either alpha alpha or beta beta -> aaaa array.
                    if( ((mod(nI(i),2).ne.0).and.(mod(nI(j),2).ne.0)) .or. &
                        ((mod(nI(i),2).eq.0).and.(mod(nI(j),2).eq.0)) ) then

                        ! Ind doesn't include diagonal terms (when iSpat = jSpat).
                        Ind=( ( (jSpat-2) * (jSpat-1) ) / 2 ) + iSpat
                        aaaa_RDM( Ind , Ind ) = aaaa_RDM( Ind , Ind ) &
                                          + ( realSignDi * realSignDi ) 

                    ! either alpha beta or beta alpha -> abab array.                                              
                    else

                        ! Ind does include diagonal terms (when iSpat = jSpat)
                        Ind=( ( (jSpat-1) * jSpat ) / 2 ) + iSpat
                        abab_RDM( Ind , Ind ) = abab_RDM( Ind , Ind ) &
                                          + ( realSignDi * realSignDi ) 
                    endif

                enddo
            enddo
        endif

    end subroutine Fill_Diag_RDM

    subroutine Fill_Sings_RDM(nI,Ex,tParity,realSignDi,realSignDj,tFill_CiCj_Symm)
! This routine adds in the contribution to the 1- and 2-RDM from determinants connected
! by a single excitation.
        USE UMatCache , only : GTID
        implicit none
        integer , intent(in) :: nI(NEl), Ex(2,2)
        logical , intent(in) :: tParity
        real(dp) , intent(in) :: realSignDi, realSignDj
        logical , intent(in) :: tFill_CiCj_Symm
        integer :: k, Indik, Indak, iSpat, aSpat, kSpat, iInd, aInd
        real(dp) :: ParityFactor, ParityFactor2

!        WRITE(6,*) '* In singles'
!        call flush(6)
!        WRITE(6,*) 'Ex(1,:)',Ex(1,:)
!        WRITE(6,*) 'Ex(2,:)',Ex(2,:)
!        WRITE(6,*) 'tParity',tParity
!        WRITE(6,*) 'nI',nI

        ParityFactor=1.D0
        IF(tParity) ParityFactor=-1.D0

        if(RDMExcitLevel.eq.1) then
            ! SymLabelList2(i), gives the orbital in position i
            ! SymLabelListInv(i), gives the position orbital i should go in.

            if(tStoreSpinOrbs) then
                iInd = Ex(1,1)
                aInd = Ex(2,1)
            else
                iInd = gtID(Ex(1,1))
                aInd = gtID(Ex(2,1))   ! These two must have the same spin.
            endif
            Indik = SymLabelListInv(iInd)    ! Position of i 
            Indak = SymLabelListInv(aInd)    ! Position of a.
            
            ! Adding to 1-RDM(i,a), ci.cj effectively.
            NatOrbMat( Indik , Indak ) = NatOrbMat( Indik , Indak ) + (ParityFactor * &
                                                             realSignDi * realSignDj )

            if(tFill_CiCj_Symm) then                                
                NatOrbMat( Indak , Indik ) = NatOrbMat( Indak , Indik ) + (ParityFactor * &
                                                             realSignDi * realSignDj )

            endif
        else
            ! Looking at elements of the type Gamma(i,k,a,k)

            ! The two determinants Di and Dj will have the same occupations except for the i and a.
            ! Any of the N-1 other electrons can be annihilated and created in the same orbital.
            ! So we run over all k = all N-1 other occupied orbitals.
            
            iSpat = gtID(Ex(1,1))
            aSpat = gtID(Ex(2,1))   ! These two must have the same spin.

            do k = 1, NEl                            

                kSpat = gtID(nI(k))

                IF(nI(k).ne.Ex(1,1)) THEN

                    if((iSpat.eq.kSpat).or.(aSpat.eq.kSpat)) then
                        ! the only array with i = j and a = b is abab.

                        Indik=( ( (max(iSpat,kSpat)-1) * max(iSpat,kSpat) ) / 2 ) + min(iSpat,kSpat)
                        Indak=( ( (max(aSpat,kSpat)-1) * max(aSpat,kSpat) ) / 2 ) + min(aSpat,kSpat)

                        ! This could be either abab or abba, but in both cases, add into the abab.
                        ! Kind of pretent the abba is of the form abab.
                        abab_RDM( Indik , Indak ) = abab_RDM( Indik , Indak ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                        if(tFill_CiCj_Symm) then
                            abab_RDM( Indak , Indik ) = abab_RDM( Indak , Indik ) + ( ParityFactor * &
                                                                         realSignDi * realSignDj )
                        endif


                    else
                        ! Checking spins of i and k.
                        ! If same, i.e alpha alpha or beta beta -> aaaa array.
                        if( ((mod(Ex(1,1),2).ne.0).and.(mod(nI(k),2).ne.0)) .or. &
                            ((mod(Ex(1,1),2).eq.0).and.(mod(nI(k),2).eq.0)) ) then

                            ! 2-RDM(i,j,a,b) is defined to have i < j and a < b, as that is how the unique 
                            ! indices are defined for i,j and a,b.
                            ! But the parity is defined so that the i -> a excitation is aligned.

                            ! I.e. we're adding these as nI(k),Ex(1,1) -> nI(k), Ex(2,1)
                            ! So if Ex(1,1) < nI(k), or Ex(2,1) < nI(k) then we need 
                            ! to switch the parity.
                            ParityFactor2 = ParityFactor
                            IF((Ex(1,1).lt.nI(k)).and.(Ex(2,1).gt.nI(k))) &
                                                    ParityFactor2 = ParityFactor * (-1.D0)
                            IF((Ex(1,1).gt.nI(k)).and.(Ex(2,1).lt.nI(k))) &
                                                    ParityFactor2 = ParityFactor * (-1.D0)

                            ! Ind doesn't include diagonal terms (when iSpat = jSpat).
                            Indik=( ( (max(iSpat,kSpat)-2) * (max(iSpat,kSpat)-1) ) / 2 ) + min(iSpat,kSpat)
                            Indak=( ( (max(aSpat,kSpat)-2) * (max(aSpat,kSpat)-1) ) / 2 ) + min(aSpat,kSpat)

                            aaaa_RDM( Indik , Indak ) = aaaa_RDM( Indik , Indak ) + ( ParityFactor2 * &
                                                                                 realSignDi * realSignDj )

                            if(tFill_CiCj_Symm) then
                                aaaa_RDM( Indak , Indik ) = aaaa_RDM( Indak , Indik ) + ( ParityFactor2 * &
                                                                                 realSignDi * realSignDj )
                            endif

                        ! either abab or abba array. 
                        ! we distinguish between these because i<j and a<b.
                        else

                            ! ordered with k's aligned, i k j k -> abab array.
                            if( ((Ex(1,1).lt.nI(k)).and.(Ex(2,1).lt.nI(k))).or. &
                                ((Ex(1,1).gt.nI(k)).and.(Ex(2,1).gt.nI(k))) ) then

                                ! It is possible for i = k or j = k if they are spat orbitals 
                                ! and have different spins.
                                Indik=( ( (max(iSpat,kSpat)-1) * max(iSpat,kSpat) ) / 2 ) + min(iSpat,kSpat)
                                Indak=( ( (max(aSpat,kSpat)-1) * max(aSpat,kSpat) ) / 2 ) + min(aSpat,kSpat)

                                abab_RDM( Indik , Indak ) = abab_RDM( Indik , Indak ) + ( ParityFactor * &
                                                                                 realSignDi * realSignDj )

                                if(tFill_CiCj_Symm) then
                                    abab_RDM( Indak , Indik ) = abab_RDM( Indak , Indik ) + ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                endif

                            ! ordered with k's not aligned, i k k j -> abba array
                            else

                                ! Ind doesn't include diagonal terms (when iSpat = jSpat).
                                Indik=( ( (max(iSpat,kSpat)-2) * (max(iSpat,kSpat)-1) ) / 2 ) + min(iSpat,kSpat)
                                Indak=( ( (max(aSpat,kSpat)-2) * (max(aSpat,kSpat)-1) ) / 2 ) + min(aSpat,kSpat)

                                abba_RDM( Indik , Indak ) = abba_RDM( Indik , Indak ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )

                                if(tFill_CiCj_Symm) then
                                    abba_RDM( Indak , Indik ) = abba_RDM( Indak , Indik ) - ( ParityFactor * &
                                                                                 realSignDi * realSignDj )
                                endif

                            endif
                        endif
                    endif
                ENDIF
            enddo

        endif

    end subroutine Fill_Sings_RDM

    subroutine Fill_Doubs_RDM(Ex,tParity,realSignDi,realSignDj,tFill_CiCj_Symm)
! This routine adds in the contribution to the 2-RDM from determinants connected
! by a double excitation.
        USE UMatCache , only : GTID
        implicit none
        integer , intent(in) :: Ex(2,2)
        logical , intent(in) :: tParity
        real(dp) , intent(in) :: realSignDi, realSignDj
        logical , intent(in) :: tFill_CiCj_Symm
        integer :: Indij, Indab, iSpat, jSpat, aSpat, bSpat
        real(dp) :: ParityFactor

        ! Adding to elements Gamma(i,j,a,b)

        ParityFactor=1.D0
        IF(tParity) ParityFactor=-1.D0

!        WRITE(6,*) '* In doubles'
!        call flush(6)
!        WRITE(6,*) 'Ex(1,:)',Ex(1,:)
!        WRITE(6,*) 'Ex(2,:)',Ex(2,:)
!        WRITE(6,*) 'tParity',tParity
!        WRITE(6,*) 'Adding realSignDi, realSignDj',realSignDi,realSignDj

        iSpat = gtID(Ex(1,1))
        jSpat = gtID(Ex(1,2))       ! Ex(1,1) < Ex(1,2)
        aSpat = gtID(Ex(2,1)) 
        bSpat = gtID(Ex(2,2))       ! Ex(2,1) < Ex(2,2)

        if((iSpat.eq.jSpat).or.(aSpat.eq.bSpat)) then

            ! i and a are different spin -> abba (but adding as abab - mult by -1).
            if( ((mod(Ex(1,1),2).eq.0).and.(mod(Ex(2,1),2).ne.0)) .or. &
                ((mod(Ex(1,1),2).ne.0).and.(mod(Ex(2,1),2).eq.0)) ) &
                    ParityFactor = ParityFactor * (-1.0_dp)

            Indij=( ( (jSpat-1) * jSpat ) / 2 ) + iSpat
            Indab=( ( (bSpat-1) * bSpat ) / 2 ) + aSpat


            abab_RDM( Indij , Indab ) = abab_RDM( Indij , Indab ) + ( ParityFactor * &
                                                         realSignDi * realSignDj )

            if(tFill_CiCj_Symm) then
                abab_RDM( Indab , Indij ) = abab_RDM( Indab , Indij ) + ( ParityFactor * &
                                                         realSignDi * realSignDj )
            endif

        else
            ! Checking spins of i and j (these must be same combination as a and b).
            ! If alpha alpha or beta beta -> aaaa array.
            if( ((mod(Ex(1,1),2).ne.0).and.(mod(Ex(1,2),2).ne.0)) .or. &
                ((mod(Ex(1,1),2).eq.0).and.(mod(Ex(1,2),2).eq.0)) ) then

                ! Don't need to worry about diagonal terms, i can't equal j.
                ! jSpat > iSpat and bSpat > aSpat
                Indij=( ( (jSpat-2) * (jSpat-1) ) / 2 ) + iSpat
                Indab=( ( (bSpat-2) * (bSpat-1) ) / 2 ) + aSpat

                aaaa_RDM( Indij , Indab ) = aaaa_RDM( Indij , Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )

                if(tFill_CiCj_Symm) then
                    aaaa_RDM( Indab , Indij ) = aaaa_RDM( Indab , Indij ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                endif

            ! Either alpha beta or beta alpha -> abab array.
            else
                
                ! if when ordering i < j and a < b, is it abab or abba.

                ! i and a are the same spin -> abab
                if( ((mod(Ex(1,1),2).eq.0).and.(mod(Ex(2,1),2).eq.0)) .or. &
                    ((mod(Ex(1,1),2).ne.0).and.(mod(Ex(2,1),2).ne.0)) ) then

                    Indij=( ( (jSpat-1) * jSpat ) / 2 ) + iSpat
                    Indab=( ( (bSpat-1) * bSpat ) / 2 ) + aSpat


                    abab_RDM( Indij , Indab ) = abab_RDM( Indij , Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )

                    if(tFill_CiCj_Symm) then
                        abab_RDM( Indab , Indij ) = abab_RDM( Indab , Indij ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                    endif

                ! i and a are different spin -> abba
                ! the only double excitation case with Indij = Indab will go in here.
                else

                    ! Don't need to worry about diagonal terms, i can't equal j.
                    ! jSpat > iSpat and bSpat > aSpat
                    Indij=( ( (jSpat-2) * (jSpat-1) ) / 2 ) + iSpat
                    Indab=( ( (bSpat-2) * (bSpat-1) ) / 2 ) + aSpat

                    abba_RDM( Indij , Indab ) = abba_RDM( Indij , Indab ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )

                    if(tFill_CiCj_Symm) then
                        abba_RDM( Indab , Indij ) = abba_RDM( Indab , Indij ) + ( ParityFactor * &
                                                                 realSignDi * realSignDj )
                    endif
                endif
            endif
        endif

    end subroutine Fill_Doubs_RDM

    subroutine FinaliseRDM()
! This routine finalises the one electron reduced density matrix stuff at the point of a softexit.
! This includes summing each of the individual matrices from each processor,
! and calling the diagonalisation routines if we want to get the occupation numbers.
        USE Logging , only : tDiagRDM
        implicit none
        INTEGER :: error
        real(dp) :: Norm_2RDM, Norm_2RDM_Inst
        CHARACTER(len=*), PARAMETER :: this_routine='FinaliseRDM'

        CALL set_timer(FinaliseRDM_Time)

        write(6,*) ''
        if(tExplicitAllRDM) then
            write(6,*) '**** RDMs CALCULATED EXPLICITLY **** '
        elseif(tHF_Ref_Explicit) then
            write(6,'(A)') ' **** RDMs CALCULATED EXPLICITLY USING THE HF AS A REFERENCE**** '
        else
            write(6,*) '**** RDMs CALCULATED STOCHASTICALLY **** '
        endif
        write(6,*) ''

        ! Combine the 1- or 2-RDM from all processors etc.

        if(RDMExcitLevel.eq.1) then

            call Finalise_1e_RDM()  

        else
            ! We always want to calculate one final RDM energy, whether or not we're 
            ! calculating the energy throughout the calculation.
            ! Unless of course, only the 1-RDM is being calculated.

            ! Calculate the energy one last time - and write out everything we need.
            tFinalRDMEnergy = .true.
            CALL Calc_Energy_from_RDM()

!            CALL Test_Energy_Calc()

            if(tPrint1RDM) call Finalise_1e_RDM()

        endif
        call MPIBarrier(error)

! Call the routines from NatOrbs that diagonalise the one electron reduced 
! density matrix.
        IF(tDiagRDM) call find_nat_orb_occ_numbers()

! This is where we would likely call any further calculations of force etc.

        CALL halt_timer(FinaliseRDM_Time)

    
    end subroutine FinaliseRDM

    subroutine Finalise_1e_RDM() 
! This routine takes the 1-RDM (NatOrbMat), normalises it, makes it 
! hermitian if required, and prints out the versions we're interested in.    
! This is only ever called at the very end of a calculation.
        use Logging , only : twrite_RDMs_to_read, twrite_normalised_RDMs
        implicit none
        integer :: i
        real(dp) :: Norm_1RDM, Trace_1RDM

        AllAccumRDMNorm = 0.D0
        IF(tHF_S_D_Ref.or.tHF_Ref_Explicit.or.tHF_S_D) &
            CALL MPIReduce(AccumRDMNorm,MPI_SUM,AllAccumRDMNorm)

        if(RDMExcitLevel.eq.1) CALL MPISum_inplace(NatOrbMat)
        
        if(iProcIndex.eq.0) then 

            ! Find the normalisation.
            call calc_1e_norms(Trace_1RDM, Norm_1RDM)

            ! Write out the unnormalised, non-hermitian OneRDM_POPS.
            if(twrite_RDMs_to_read) call Write_out_1RDM(Norm_1RDM,.false.)

            ! Enforce the hermiticity condition.  If the RDMExcitLevel is not 1, the 
            ! 1-RDM has been constructed from the hermitian 2-RDM, so this will not 
            ! be necessary.
            ! The HF_Ref and HF_S_D_Ref cases are not hermitian by definition.
            if((RDMExcitLevel.eq.1).and.(.not.(tHF_Ref_Explicit.or.tHF_S_D_Ref))) &
                call make_1e_rdm_hermitian(Norm_1RDM)

            ! Write out the final, normalised, hermitian OneRDM.                
            if(twrite_normalised_RDMs) call Write_out_1RDM(Norm_1RDM,.true.)

        endif

    end subroutine Finalise_1e_RDM

    subroutine calc_1e_norms(Trace_1RDM, Norm_1RDM)
! We want to 'normalise' the reduced density matrices.
! These are not even close to being normalised at the moment, because of the way they are 
! calculated on the fly.
! They should be calculated from a normalised wavefunction.
! But we know that the trace of the one electron reduced density matrix must be equal to 
! the number of the electrons.
! We can use this to find the factor we must divide the 1RDM through by.
        implicit none                            
        real(dp) , intent(out) :: Trace_1RDM, Norm_1RDM
        integer :: i

        Trace_1RDM = 0.D0
        Norm_1RDM = 0.D0

        do i = 1, NoOrbs
            Trace_1RDM = Trace_1RDM + NatOrbMat(i,i)
        enddo

        IF(tHF_S_D_Ref.or.tHF_Ref_Explicit.or.tHF_S_D) THEN
            Norm_1RDM = 1.D0 / AllAccumRDMNorm
        ELSE
            ! Sum of diagonal elements of 1 electron RDM must equal NEl, 
            ! number of electrons.
            Norm_1RDM = ( REAL(NEl,8) / Trace_1RDM )
        ENDIF

!        if(tFinalRDMEnergy) then
!            WRITE(6,*) 'AllAccumRDMNorm',AllAccumRDMNorm
!            WRITE(6,*) 'Norm_1RDM',Norm_1RDM
!            WRITE(6,*) 'Trace_1RDM',Trace_1RDM
!        endif

        !Need to multiply each element of the 1 electron reduced density matrices 
        !by NEl / Trace_1RDM,
        !and then add it's contribution to the energy.

    end subroutine calc_1e_norms

    subroutine make_1e_rdm_hermitian(Norm_1RDM)
! Simply average the 1-RDM(i,j) and 1-RDM(j,i) elements which should be equal in a perfect world.    
        implicit none 
        real(dp) , intent(in) :: Norm_1RDM
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity 
        integer :: i, j
        real(dp) :: Temp

        Max_Error_Hermiticity = 0.D0
        Sum_Error_Hermiticity = 0.D0
        do i = 1, NoOrbs
            do j = i, NoOrbs
                IF((abs((NatOrbMat(SymLabelListInv(i),SymLabelListInv(j))*Norm_1RDM) - &
                        (NatOrbMat(SymLabelListInv(j),SymLabelListInv(i))*Norm_1RDM))).gt.Max_Error_Hermiticity) &
                    Max_Error_Hermiticity = abs((NatOrbMat(SymLabelListInv(i),SymLabelListInv(j))*Norm_1RDM) - &
                                                (NatOrbMat(SymLabelListInv(j),SymLabelListInv(i))*Norm_1RDM))

                Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                        abs((NatOrbMat(SymLabelListInv(i),SymLabelListInv(j))*Norm_1RDM) - &
                                            (NatOrbMat(SymLabelListInv(j),SymLabelListInv(i))*Norm_1RDM))

                Temp = (NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) + &
                        NatOrbMat(SymLabelListInv(j),SymLabelListInv(i)))/2.D0

                NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) = Temp
                NatOrbMat(SymLabelListInv(j),SymLabelListInv(i)) = Temp
            enddo
        enddo

        ! Output the hermiticity errors.
        write(6,'(A29,F30.20)') ' MAX ABS ERROR IN HERMITICITY', Max_Error_Hermiticity
        write(6,'(A29,F30.20)') ' SUM ABS ERROR IN HERMITICITY', Sum_Error_Hermiticity

    end subroutine make_1e_rdm_hermitian

    subroutine Write_out_1RDM(Norm_1RDM,tNormalise)
! This routine writes out the OneRDM.
! If tNormalise is true, we are printing the normalised, hermitian matrix.
! Otherwise, Norm_1RDM is ignored and we print both 1-RDM(i,j) and 1-RDM(j,i) (in binary) 
! for the OneRDM_POPS file to be read in in a restart calculation.
        USE UMatCache, only: GTID
        implicit none
        real(dp) , intent(in) :: Norm_1RDM
        logical , intent(in) :: tNormalise
        integer :: i, j, iSpat, jSpat
        integer :: OneRDM_unit
 
        if(tNormalise) then
            ! Haven't got the capabilities to produce multiple 1-RDMs yet.
            write(6,*) 'Writing out the *normalised* 1 electron density matrix to file'
            call flush(6)
            OneRDM_unit = get_free_unit()
            OPEN(OneRDM_unit,file='OneRDM',status='unknown')
        else
            ! Only every write out 1 of these at the moment.
            write(6,*) 'Writing out the *unnormalised* 1 electron density matrix to file for reading in'
            call flush(6)
            OneRDM_unit = get_free_unit()
            OPEN(OneRDM_unit,file='OneRDM_POPS',status='unknown',form='unformatted')
        endif
 
        ! Currently always printing 1-RDM in spin orbitals.
        do i = 1, nBasis
            do j = 1, nBasis
                if(tStoreSpinOrbs) then
                    if(NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)).ne.0.D0) then 
                        if(tNormalise.and.((i.le.j).or.tHF_Ref_Explicit.or.tHF_S_D_Ref)) then
                            write(OneRDM_unit,"(2I6,G25.17)") i,j, & 
                                NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) * Norm_1RDM
                        elseif(.not.tNormalise) then
                            ! For the pops, we haven't made the 1-RDM hermitian yet, 
                            ! so print both the 1-RDM(i,j) and 1-RDM(j,i) elements.
                            ! This is written in binary.
                            write(OneRDM_unit) i,j,NatOrbMat(SymLabelListInv(i),SymLabelListInv(j))
                        endif
                    endif
                else
                    iSpat = gtID(i)
                    jSpat = gtID(j)
                    if(NatOrbMat(SymLabelListInv(iSpat),SymLabelListInv(jSpat)).ne.0.D0) then 
                        if(tNormalise.and.((i.le.j).or.tHF_Ref_Explicit.or.tHF_S_D_Ref)) then
                            if(((mod(i,2).eq.0).and.(mod(j,2).eq.0)).or.&
                                ((mod(i,2).ne.0).and.(mod(j,2).ne.0))) then
                                write(OneRDM_unit,"(2I6,G25.17)") i,j, & 
                                    ( NatOrbMat(SymLabelListInv(iSpat),SymLabelListInv(jSpat)) &
                                                                    * Norm_1RDM ) / 2.0_dp
                            endif
                        elseif(.not.tNormalise) then
                            ! The popsfile can be printed in spatial orbitals.
                            if((mod(i,2).eq.0).and.(mod(j,2).eq.0)) then
                                write(OneRDM_unit) iSpat,jSpat, & 
                                    NatOrbMat(SymLabelListInv(iSpat),SymLabelListInv(jSpat)) 
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
                
        close(OneRDM_unit)

    end subroutine Write_out_1RDM


    subroutine Finalise_2e_RDM(Norm_2RDM_Inst, Norm_2RDM) 
! This routine sums, normalises, hermitian-ises, and prints the 2-RDMs.    
! This may be called multiple times if we want to print multiple 2-RDMs.
        use Logging , only : twrite_RDMs_to_read, twrite_normalised_RDMs, &
                                RDMEnergyIter, IterWriteRDMs, tWriteMultRDMs
        use FciMCData , only : IterRDMStart
        implicit none
        real(dp) , intent(out) :: Norm_2RDM_Inst, Norm_2RDM
        real(dp) :: AllAccumRDMNorm_Inst
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity
        integer :: RDM_Cycles

        ! If Iter = 0, this means we have just read in the TwoRDM_POPS_a*** matrices into All_a***_RDM, and 
        ! just want to calculate the old energy.
        ! Don't need to do all this stuff here, because a***_RDM will be empty.
        if(Iter.ne.0) then

            ! All the arrays are summed into the one on processor 0.
            CALL MPISum_inplace(aaaa_RDM(:,:))
            CALL MPISum_inplace(abab_RDM(:,:))
            CALL MPISum_inplace(abba_RDM(:,:))

            ! The TwoElRDM on the root is now the sum of all 'instantaneous' RDMs (summed over 
            ! the energy update cycle).
            ! Whereas AllTwoElRDM is accumulated over the entire run.

            ! The AllTwoElRDM's are actually averaged over the iterations.
            ! This is to keep the trace's etc not too large, not sure if it's the best or not.
            RDM_Cycles =  ( Iter - IterRDMStart ) / RDMEnergyIter
            
            if(iProcIndex.eq.0) then
                All_aaaa_RDM(:,:) = ( All_aaaa_RDM(:,:) * &
                                    ( real(RDM_Cycles,dp) / ( real(RDM_Cycles,dp) + 1.0_dp ) ) ) &
                                    + ( aaaa_RDM(:,:) / ( real(RDM_Cycles,dp) + 1.0_dp ) )
                All_abab_RDM(:,:) = ( All_abab_RDM(:,:) * &
                                    ( real(RDM_Cycles,dp) / ( real(RDM_Cycles,dp) + 1.0_dp ) ) ) &
                                    + ( abab_RDM(:,:) / ( real(RDM_Cycles,dp) + 1.0_dp ) )
                All_abba_RDM(:,:) = ( All_abba_RDM(:,:) * &
                                    ( real(RDM_Cycles,dp) / ( real(RDM_Cycles,dp) + 1.0_dp ) ) ) &
                                    + ( abba_RDM(:,:) / ( real(RDM_Cycles,dp) + 1.0_dp ) )
            endif

            AllAccumRDMNorm_Inst = 0.D0
            if(tHF_S_D_Ref.or.tHF_Ref_Explicit.or.tHF_S_D) then
                CALL MPIReduce(AccumRDMNorm_Inst,MPI_SUM,AllAccumRDMNorm_Inst)
                AllAccumRDMNorm = ( AllAccumRDMNorm * &
                                    ( real(RDM_Cycles,dp) / ( real(RDM_Cycles,dp) + 1.0_dp ) ) ) &
                                    + ( AllAccumRDMNorm_Inst / ( real(RDM_Cycles,dp) + 1.0_dp ) )
            endif
        endif

        if(iProcIndex.eq.0) then

            ! Calculate the normalisations.
            call calc_2e_norms(AllAccumRDMNorm_Inst, Norm_2RDM_Inst, Norm_2RDM)

            ! Print out the relevant 2-RDMs.
            if( tFinalRDMEnergy .or. &
                ( tWriteMultRDMs .and. (mod((Iter - IterRDMStart)+1,IterWriteRDMs).eq.0) ) ) then

                if(tFinalRDMEnergy) then
                    ! Only ever want to print the POPS 2-RDMs (for reading in) at the end.
                    if(twrite_RDMs_to_read) call Write_out_2RDM(Norm_2RDM,.false.)

                    ! We also don't want to make the 2-RDMs hermitian until the end, so that we can 
                    ! get the hermiticity error from the final matrix.
                    if(.not.(tHF_Ref_Explicit.or.tHF_S_D_Ref)) then
                        call make_2e_rdm_hermitian(Norm_2RDM, Max_Error_Hermiticity, Sum_Error_Hermiticity)

                        write(6,'(A29,F30.20)') ' MAX ABS ERROR IN HERMITICITY', Max_Error_Hermiticity
                        write(6,'(A29,F30.20)') ' SUM ABS ERROR IN HERMITICITY', Sum_Error_Hermiticity
                    endif
                endif

                ! This writes out the normalised, hermitian 2-RDMs.
                if(twrite_normalised_RDMs) call Write_out_2RDM(Norm_2RDM,.true.)

            endif
        endif

    end subroutine 

    subroutine calc_2e_norms(AllAccumRDMNorm_Inst, Norm_2RDM_Inst, Norm_2RDM)
! We want to 'normalise' the reduced density matrices.
! These are not even close to being normalised at the moment, because of the way they are 
! calculated on the fly.
! They should be calculated from a normalised wavefunction.
!
! We also know that the trace of the two electron reduced density matrix must be equal to the 
! number of electron pairs in the system = 1/2 N ( N - 1), so we can do the same for the 2RDM.
        implicit none                            
        real(dp) , intent(in) :: AllAccumRDMNorm_Inst
        real(dp) , intent(out) :: Norm_2RDM_Inst, Norm_2RDM
        integer :: i

        ! Find the current, unnormalised trace of each matrix.
        ! TODO: This can be merged into the spin averaging when everything is working.

        Trace_2RDM_Inst = 0.D0
        Trace_2RDM = 0.D0

        do i = 1, ((SpatOrbs*(SpatOrbs+1))/2)
            if(i.le.((SpatOrbs*(SpatOrbs-1))/2)) then
                Trace_2RDM_Inst = Trace_2RDM_Inst + aaaa_RDM(i,i)
                Trace_2RDM = Trace_2RDM + All_aaaa_RDM(i,i)
            endif
            Trace_2RDM_Inst = Trace_2RDM_Inst + abab_RDM(i,i)
            Trace_2RDM = Trace_2RDM + All_abab_RDM(i,i)
        enddo

        Norm_2RDM_Inst = 0.D0
        Norm_2RDM = 0.D0

        IF(tHF_S_D_Ref.or.tHF_Ref_Explicit.or.tHF_S_D) THEN
            Norm_2RDM_Inst = 1.D0 / AllAccumRDMNorm_Inst
            Norm_2RDM = 1.D0 / AllAccumRDMNorm
        ELSE
            ! Sum of diagonal elements of 2 electron RDM must equal number of 
            ! pairs of electrons, = NEl ( NEl - 1 ) / 2
            Norm_2RDM_Inst = ( (0.50 * (REAL(NEl) * (REAL(NEl) - 1.D0))) / Trace_2RDM_Inst )
            Norm_2RDM = ( (0.50 * (REAL(NEl) * (REAL(NEl) - 1.D0))) / Trace_2RDM )
        ENDIF

!        if(tFinalRDMEnergy) then
!            WRITE(6,*) 'AllAccumRDMNorm',AllAccumRDMNorm
!            WRITE(6,*) 'Norm_2RDM',Norm_2RDM
!            WRITE(6,*) 'Trace_2RDM',Trace_2RDM
!        endif

        !Need to multiply each element of the 1 electron reduced density matrices 
        !by NEl / Trace_1RDM,
        !and then add it's contribution to the energy.

    end subroutine calc_2e_norms

    subroutine make_2e_rdm_hermitian(Norm_2RDM, Max_Error_Hermiticity, Sum_Error_Hermiticity)
! This averages 2-RDM(i,j;a,b) and 2-RDM(a,b;i,j) or equivalently 2-RDM(Ind1,Ind2) and 2-RDM(Ind2,Ind1).
        implicit none 
        real(dp) , intent(in) :: Norm_2RDM
        real(dp) , intent(out) :: Max_Error_Hermiticity, Sum_Error_Hermiticity 
        integer :: i, j
        real(dp) :: Temp

        Max_Error_Hermiticity = 0.D0
        Sum_Error_Hermiticity = 0.D0

        do i = 1, ((SpatOrbs*(SpatOrbs+1))/2)
            do j = i+1, ((SpatOrbs*(SpatOrbs+1))/2)

                if((i.le.((SpatOrbs*(SpatOrbs-1))/2)).and.(j.le.((SpatOrbs*(SpatOrbs-1))/2))) then

                    IF((abs((All_aaaa_RDM(i,j)*Norm_2RDM)-(All_aaaa_RDM(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                        Max_Error_Hermiticity = abs((All_aaaa_RDM(i,j)*Norm_2RDM)-(All_aaaa_RDM(j,i)*Norm_2RDM))

                    Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                            abs((All_aaaa_RDM(i,j)*Norm_2RDM)-(All_aaaa_RDM(j,i)*Norm_2RDM))

                    Temp = (All_aaaa_RDM(i,j) + All_aaaa_RDM(j,i)) / 2.D0

                    All_aaaa_RDM(i,j) = Temp
                    All_aaaa_RDM(j,i) = Temp

                    IF((abs((All_abba_RDM(i,j)*Norm_2RDM)-(All_abba_RDM(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                        Max_Error_Hermiticity = abs((All_abba_RDM(i,j)*Norm_2RDM)-(All_abba_RDM(j,i)*Norm_2RDM))

                    Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                            abs((All_abba_RDM(i,j)*Norm_2RDM)-(All_abba_RDM(j,i)*Norm_2RDM))

                    Temp = (All_abba_RDM(i,j) + All_abba_RDM(j,i)) / 2.D0

                    All_abba_RDM(i,j) = Temp
                    All_abba_RDM(j,i) = Temp
                endif

                IF((abs((All_abab_RDM(i,j)*Norm_2RDM)-(All_abab_RDM(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                    Max_Error_Hermiticity = abs((All_abab_RDM(i,j)*Norm_2RDM)-(All_abab_RDM(j,i)*Norm_2RDM))

                Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                        abs((All_abab_RDM(i,j)*Norm_2RDM)-(All_abab_RDM(j,i)*Norm_2RDM))

                Temp = (All_abab_RDM(i,j) + All_abab_RDM(j,i)) / 2.D0

                All_abab_RDM(i,j) = Temp
                All_abab_RDM(j,i) = Temp


            enddo
        enddo

    end subroutine make_2e_rdm_hermitian


    subroutine Write_out_2RDM(Norm_2RDM,tNormalise)
! Writes out the 2-RDMs.  If tNormalise is true, we print the normalised (hermitian) matrix.    
! Otherwise we print the unnormalised 2-RDMs, and we print (in binary) both 2-RDM(Ind1,Ind2) 
! and 2-RDM(Ind2,Ind1) because this matrix wont be hermitian.

! While, for instance, the TwoRDM_aaaa so far has actually been a sum of the aaaa elements and 
! the bbbb elements.  We only want to print the aaaa elements.
        use util_mod , only : get_unique_filename
        use Logging , only : tWriteMultRDMs
        implicit none
        real(dp) , intent(in) :: Norm_2RDM
        logical , intent(in) :: tNormalise
        real(dp) :: Tot_Spin_Projection, SpinPlus, SpinMinus
        real(dp) :: ParityFactor,Divide_Factor 
        integer :: i, j, a, b, Ind1_aa, Ind1_ab, Ind2_aa, Ind2_ab
        integer :: aaaa_RDM_unit, abab_RDM_unit, abba_RDM_unit
        character(255) :: TwoRDM_aaaa_name, TwoRDM_abab_name, TwoRDM_abba_name

        if(tNormalise) then
            write(6,*) 'Writing out the *normalised* 2 electron density matrix to file'
            call flush(6)
            ! This takes the TwoRDM_aaaa, and if tWriteMultPops is true (and given that we've put 
            ! .true. in the 3rd position, it'll find the next unused TwoRDM_aaaa.X file name.
            call get_unique_filename('TwoRDM_aaaa',tWriteMultRDMs,.true.,1,TwoRDM_aaaa_name)
            aaaa_RDM_unit = get_free_unit()
            OPEN(aaaa_RDM_unit,file=TwoRDM_aaaa_name,status='unknown')

            call get_unique_filename('TwoRDM_abab',tWriteMultRDMs,.true.,1,TwoRDM_abab_name)
            abab_RDM_unit = get_free_unit()
            OPEN(abab_RDM_unit,file=TwoRDM_abab_name,status='unknown')

            call get_unique_filename('TwoRDM_abba',tWriteMultRDMs,.true.,1,TwoRDM_abba_name)
            abba_RDM_unit = get_free_unit()
            OPEN(abba_RDM_unit,file=TwoRDM_abba_name,status='unknown')
        else
            write(6,*) 'Writing out the *unnormalised* 2 electron density matrix to file for reading in'
            call flush(6)
            aaaa_RDM_unit = get_free_unit()
            OPEN(aaaa_RDM_unit,file='TwoRDM_POPS_aaaa',status='unknown',form='unformatted')
            abab_RDM_unit = get_free_unit()
            OPEN(abab_RDM_unit,file='TwoRDM_POPS_abab',status='unknown',form='unformatted')
            abba_RDM_unit = get_free_unit()
            OPEN(abba_RDM_unit,file='TwoRDM_POPS_abba',status='unknown',form='unformatted')
        endif
        
        Tot_Spin_Projection = 0.D0
        do i = 1, SpatOrbs

            do j = i, SpatOrbs

                Ind1_aa = ( ( (j-2) * (j-1) ) / 2 ) + i
                Ind1_ab = ( ( (j-1) * j ) / 2 ) + i

                ! TODO : Fix this.
!                    call sum_in_spin_proj(i,j,Ind1,Norm_2RDM,Tot_Spin_Projection)

                do a = 1, SpatOrbs

                    do b = a, SpatOrbs

                        Ind2_aa = ( ( (b-2) * (b-1) ) / 2 ) + a
                        Ind2_ab = ( ( (b-1) * b ) / 2 ) + a

                        ! usually each element will have two contributions (from aaaa and bbbb).
                        ! we then need to divide each by 2.
                        ! but in cases where i and j, and a and b, are in the same spatial 
                        ! orbital, there will be only one contribution.
                        if((i.eq.j).and.(a.eq.b)) then
                            Divide_Factor = 1.0_dp
                        else
                            Divide_Factor = 2.0_dp
                        endif
                        
                        if((i.ne.j).and.(a.ne.b)) then

                            if( All_aaaa_RDM(Ind1_aa,Ind2_aa).ne.0.0_dp) then
                                ! If we're normalising (and have made the matrix hermitian) we only 
                                ! need to write out Ind1 < Ind2.
                                ! Otherwise we print out Ind1, Ind2 and Ind2, Ind1 so we can 
                                ! find the hermiticity error in the final matrix (after all runs).
                                if(tNormalise.and.((Ind1_aa.le.Ind2_aa).or.tHF_Ref_Explicit.or.tHF_S_D_Ref)) then
                                    if(tFinalRDMEnergy) then
                                        ! For the final calculation, the 2-RDMs will have been made hermitian.
                                        write(aaaa_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                                ( All_aaaa_RDM(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                    else
                                        ! If we're printing the 2-RDMs early (using WRITERDMSEVERY), the actual 
                                        ! matrix will not be hermitian, but we want to print a hermitian version.
                                        ! Average the values here.
                                        write(aaaa_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                            ( ((All_aaaa_RDM(Ind1_aa,Ind2_aa) + All_aaaa_RDM(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                * Norm_2RDM ) / Divide_Factor
                                    endif
                                elseif(.not.tNormalise) then
                                    ! for the POPS files, print everything to binary.
                                    ! no divide factor, we just read them in as is.
                                    write(aaaa_RDM_unit) i,j,a,b, &
                                            All_aaaa_RDM(Ind1_aa,Ind2_aa) 
                                endif
                            endif

                            if( All_abba_RDM(Ind1_aa,Ind2_aa).ne.0.0_dp) then
                                if(tNormalise.and.((Ind1_aa.le.Ind2_aa).or.tHF_Ref_Explicit.or.tHF_S_D_Ref)) then
                                    if(tFinalRDMEnergy) then
                                        write(abba_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                            ( All_abba_RDM(Ind1_aa,Ind2_aa) * Norm_2RDM ) / Divide_Factor
                                    else
                                        write(abba_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                            ( ((All_abba_RDM(Ind1_aa,Ind2_aa) + All_abba_RDM(Ind2_aa,Ind1_aa))/2.0_dp) &
                                                                    * Norm_2RDM ) / Divide_Factor
                                    endif
                                elseif(.not.tNormalise) then
                                    write(abba_RDM_unit) i,j,a,b, &
                                        All_abba_RDM(Ind1_aa,Ind2_aa) 
                                endif
                            endif

                        endif

                        if( All_abab_RDM(Ind1_ab,Ind2_ab).ne.0.0_dp) then
                            if(tNormalise.and.((Ind1_ab.le.Ind2_ab).or.tHF_Ref_Explicit.or.tHF_S_D_Ref)) then
                                if(tFinalRDMEnergy) then
                                    write(abab_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                        ( All_abab_RDM(Ind1_ab,Ind2_ab) * Norm_2RDM ) / Divide_Factor
                                else
                                    write(abab_RDM_unit,"(4I6,G25.17)") i,j,a,b, &
                                        ( ((All_abab_RDM(Ind1_ab,Ind2_ab) + All_abab_RDM(Ind2_ab,Ind1_ab))/2.0_dp) &
                                                                * Norm_2RDM ) / Divide_Factor
                                endif
                            elseif(.not.tNormalise) then
                                write(abab_RDM_unit) i,j,a,b, &
                                    All_abab_RDM(Ind1_ab,Ind2_ab) 
                            endif
                        endif

                    enddo
                enddo

            enddo

        enddo
        close(aaaa_RDM_unit)
        close(abab_RDM_unit)
        close(abba_RDM_unit)

!        Tot_Spin_Projection = Tot_Spin_Projection + (3.D0 * real(NEl,dp))
!! Tot_Spin_Projection is now equal to 4S(S+1) - find S. 
!        Tot_Spin_Projection = Tot_Spin_Projection/4.D0
!        if((1.D0 + 4.D0*Tot_Spin_Projection).le.0) then
!            call Warning_neci('Write_out_1and_2RDM',"Complex spin calculated from density matrices!")
!        else
!            SpinPlus = (-1.D0 + sqrt(1.D0 + 4.D0*Tot_Spin_Projection))/2.D0
!            SpinMinus = (-1.D0 + sqrt(1.D0 + 4.D0*Tot_Spin_Projection))/2.D0
!!            write(6,*) 'SpinPlus',SpinPlus
!!            write(6,*) 'SpinMinus',SpinMinus
!        endif

!        if(RDMExcitLevel.eq.3) then
!            write(6,*) ''
!            write(6,'(A22,F30.20)') ' TOTAL SPIN PROJECTION', Max(SpinPlus,SpinMinus) 
!            write(6,*) ''
!        endif

    end subroutine Write_out_2RDM


    SUBROUTINE Calc_Energy_from_RDM()
! This routine takes the 1 electron and 2 electron reduced density matrices 
! and calculated the energy they give.    
! The equation for the energy is as follows:
!
!   E = Tr(h1 1RDM) + 1/2 Tr(h2 2RDM)
!
! where h1 are the 2 index integrals, and h2 the 4 index integrals.  The traces, Tr, 
! are given by:
!   Tr(h1 1RDM) = Sum_i,j [ h1(i,j) 1RDM(j,i) ]
!   Tr(h2 2RDM) = Sum_i,j;k,l [ h2(i,j;k,l) 2RDM(k,l;i,j) ]
        USE IntegralsData , only : UMAT
        USE UMatCache , only : UMatInd
        USE RotateOrbsMod , only : SymLabelList2
        USE UMatCache , only : GTID
        implicit none
        real(dp) :: Norm_2RDM, Norm_2RDM_Inst
        INTEGER :: i,j,a,b,Ind1_aa,Ind1_ab,Ind2_aa,Ind2_ab,ierr
        INTEGER :: iSpin, jSpin, error
        REAL(dp) :: RDMEnergy_Inst, RDMEnergy, Coul, Exch, Parity_Factor 
        REAL(dp) :: Trace_2RDM_New, RDMEnergy1, RDMEnergy2

        CALL set_timer(RDMEnergy_Time,30)

        Trace_2RDM_New = 0.D0

        RDMEnergy_Inst = 0.D0
        RDMEnergy1 = 0.D0
        RDMEnergy2 = 0.D0
        RDMEnergy = 0.D0
    
        ! Normalise, make hermitian, print etc.
        call Finalise_2e_RDM(Norm_2RDM_Inst, Norm_2RDM)

        if(tFinalRDMEnergy) then
            write(6,*) ''
            write(6,*) 'Calculating the final RDM energy'
        endif

        if(iProcIndex.eq.0) then

            do i = 1, SpatOrbs
                iSpin = 2 * i

                do j = i, SpatOrbs
                    jSpin = 2 * j

                    Ind1_aa = ( ( (j-2) * (j-1) ) / 2 ) + i
                    Ind1_ab = ( ( (j-1) * j ) / 2 ) + i

                    do a = 1, SpatOrbs

                        ! Adding in contributions effectively from the 1-RDM (although these are calculated 
                        ! from the 2-RDM.
                        call calc_1RDM_energy(i,j,a,iSpin,jSpin, Norm_2RDM, Norm_2RDM_Inst, &
                                                    RDMEnergy_Inst, RDMEnergy1)

                        do b = a, SpatOrbs

                            Ind2_aa = ( ( (b-2) * (b-1) ) / 2 ) + a
                            Ind2_ab = ( ( (b-1) * b ) / 2 ) + a

                            ! UMAT in chemical notation.
                            ! In spin or spatial orbitals.
                            Coul = REAL(UMAT(UMatInd(i,j,a,b,0,0)),8)
                            Exch = REAL(UMAT(UMatInd(i,j,b,a,0,0)),8)

                            if((i.ne.j).and.(a.ne.b)) then
                                ! Cannot get i=j or a=b contributions in aaaa.
                                RDMEnergy_Inst = RDMEnergy_Inst + ( aaaa_RDM(Ind1_aa,Ind2_aa) &
                                                                    * Norm_2RDM_Inst * ( Coul - Exch ) )
                                RDMEnergy2 = RDMEnergy2 + ( All_aaaa_RDM(Ind1_aa,Ind2_aa) &
                                                            * Norm_2RDM * ( Coul - Exch ) ) 
                                if(Ind1_aa.eq.Ind2_aa) &
                                    Trace_2RDM_New = Trace_2RDM_New + &
                                                        All_aaaa_RDM(Ind1_aa,Ind2_aa) * Norm_2RDM

                                ! For abab cases, coul element will be non-zero, exchange zero.
                                RDMEnergy_Inst = RDMEnergy_Inst + ( abab_RDM(Ind1_ab,Ind2_ab) &
                                                                        * Norm_2RDM_Inst * Coul )
                                RDMEnergy2 = RDMEnergy2 + ( All_abab_RDM(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * Coul ) 

                                if(Ind1_ab.eq.Ind2_ab) &
                                    Trace_2RDM_New = Trace_2RDM_New + &
                                                        All_abab_RDM(Ind1_ab,Ind2_ab) * Norm_2RDM

                                ! For abba cases, coul element will be zero, exchange non-zero.
                                RDMEnergy_Inst = RDMEnergy_Inst - ( abba_RDM(Ind1_aa,Ind2_aa) &
                                                                        * Norm_2RDM_Inst * Exch )
                                RDMEnergy2 = RDMEnergy2 - ( All_abba_RDM(Ind1_aa,Ind2_aa) &
                                                            * Norm_2RDM * Exch ) 


                            else
                                ! i = j or a = b
                                ! abab has both abab and abba elements in them effectively.
                                ! half will have non-zero coul, and half non-zero exchange.

                                RDMEnergy_Inst = RDMEnergy_Inst + ( 0.5_dp * abab_RDM(Ind1_ab,Ind2_ab) &
                                                                        * Norm_2RDM_Inst * Coul ) &
                                                                + ( 0.5_dp * abab_RDM(Ind1_ab,Ind2_ab) &
                                                                        * Norm_2RDM_Inst * Exch )
                                RDMEnergy2 = RDMEnergy2 + ( 0.5_dp * All_abab_RDM(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * Coul ) &
                                                        + ( 0.5_dp * All_abab_RDM(Ind1_ab,Ind2_ab) &
                                                            * Norm_2RDM * Exch )
                                if(Ind1_ab.eq.Ind2_ab) &
                                    Trace_2RDM_New = Trace_2RDM_New + &
                                                        All_abab_RDM(Ind1_ab,Ind2_ab) * Norm_2RDM

                            endif

                       enddo

                    enddo
                enddo
            enddo

            RDMEnergy_Inst = RDMEnergy_Inst + Ecore
            RDMEnergy = RDMEnergy1 + RDMEnergy2 + Ecore 

            ! Obviously this 'instantaneous' energy is actually accumulated between energy 
            ! print outs.
            WRITE(Energies_unit, "(I31,2F30.15)") Iter+PreviousCycles, RDMEnergy_Inst, RDMEnergy

            if(tFinalRDMEnergy) then
                write(6,*) 'Trace of 2-el-RDM before normalisation : ',Trace_2RDM
                write(6,*) 'Trace of 2-el-RDM after normalisation : ',Trace_2RDM_New
                write(6,*) 'Energy contribution from the 1-RDM: ',RDMEnergy1
                write(6,*) 'Energy contribution from the 2-RDM: ',RDMEnergy2
                write(6,'(A64,F30.20)') ' *TOTAL ENERGY* CALCULATED USING THE *REDUCED &
                                            &DENSITY MATRICES*:',RDMEnergy
                write(6,*) ''
                CLOSE(Energies_unit) 
            endif

        endif

        ! Zero all the 'instantaneous' stuff.
        aaaa_RDM(:,:) = 0.0_dp
        abab_RDM(:,:) = 0.0_dp
        abba_RDM(:,:) = 0.0_dp
        AccumRDMNorm_Inst = 0.0_dp
        Trace_2RDM_Inst = 0.0_dp

        CALL halt_timer(RDMEnergy_Time)

    END SUBROUTINE Calc_Energy_from_RDM


    subroutine calc_1RDM_energy(i,j,a,iSpin,jSpin,Norm_2RDM,Norm_2RDM_Inst,&
                                                    RDMEnergy_Inst,RDMEnergy1)
! This routine calculates the 1-RDM part of the RDM energy, and constructs the 
! 1-RDM if required for diagonalisation or something.
        ! gamma(i,j) = [1/(NEl - 1)] * SUM_a Gamma(i,a,j,a) 
        ! want to calculate:    gamma(i,j) * h_ij
        ! h_ij => TMAT2D(iSpin,jSpin)
        USE OneEInts , only : TMAT2D
        USE Logging , only : tDiagRDM
        implicit none
        integer , intent(in) :: i,j,a,iSpin,jSpin
        real(dp) , intent(in) :: Norm_2RDM, Norm_2RDM_Inst
        real(dp) , intent(inout) :: RDMEnergy_Inst, RDMEnergy1
        real(dp) :: Parity_Factor
        integer :: Ind1_1e_ab, Ind2_1e_ab
        integer :: Ind1_1e_aa, Ind2_1e_aa


        Ind1_1e_ab = ( ( (max(i,a)-1) * max(i,a) ) / 2 ) + min(i,a)
        Ind2_1e_ab = ( ( (max(j,a)-1) * max(j,a) ) / 2 ) + min(j,a)

        ! for i a -> j a excitation, when lined up as min max -> min max, 
        ! if a's are aligned, only a b a b arrays contain single excitations, 
        ! if a's not aligned, a b b a.
        ! all a a a a will contain single excitations.
        if(((i.le.a).and.(j.le.a)).or.((i.ge.a).and.(j.ge.a))) then
            RDMEnergy_Inst = RDMEnergy_Inst + ( (abab_RDM(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM_Inst) &
                                                * REAL(TMAT2D(iSpin,jSpin),8) &
                                                * (1.0_dp / real(NEl - 1,dp)) )
            RDMEnergy1 = RDMEnergy1 + ( (All_abab_RDM(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM) &
                                                * REAL(TMAT2D(iSpin,jSpin),8) &
                                                * (1.0_dp / real(NEl - 1,dp)) )

            if((tDiagRDM.or.tPrint1RDM).and.tFinalRDMEnergy) then                                                
                NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) = &
                            NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) &
                                        + ( All_abab_RDM(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) ) 
            endif

            ! For Gamma elements corresponding to 1-RDMs ( Gamma(i,a,j,a) ), we're only considering 
            ! i =< j and therefore we need to sum in the opposite contribution too.
            if(Ind1_1e_ab.ne.Ind2_1e_ab) then                                                                
                RDMEnergy_Inst = RDMEnergy_Inst + ( (abab_RDM(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM_Inst) &
                                                    * REAL(TMAT2D(jSpin,iSpin),8) &
                                                    * (1.0_dp / real(NEl - 1,dp)) )
                RDMEnergy1 = RDMEnergy1 + ( (All_abab_RDM(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM) &
                                                    * REAL(TMAT2D(jSpin,iSpin),8) &
                                                    * (1.0_dp / real(NEl - 1,dp)) )

                if((tDiagRDM.or.tPrint1RDM).and.tFinalRDMEnergy) then                                                
                    NatOrbMat(SymLabelListInv(j),SymLabelListInv(i)) = &
                            NatOrbMat(SymLabelListInv(j),SymLabelListInv(i)) &
                                        + ( All_abab_RDM(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) ) 
                endif

            endif

            ! But since we're running over all a, i a and a i will both be counted, but i i only once 
            ! (whereas it should be counted twice).
            if((i.eq.j).and.(i.eq.a)) then                                                                
                RDMEnergy_Inst = RDMEnergy_Inst + ( (abab_RDM(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM_Inst) &
                                            * REAL(TMAT2D(iSpin,jSpin),8) &
                                            * (1.0_dp / real(NEl - 1,dp)) )
                RDMEnergy1 = RDMEnergy1 + ( (All_abab_RDM(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM) &
                                            * REAL(TMAT2D(iSpin,jSpin),8) &
                                            * (1.0_dp / real(NEl - 1,dp)) )

                if((tDiagRDM.or.tPrint1RDM).and.tFinalRDMEnergy) then                                                
                    NatOrbMat(SymLabelListInv(j),SymLabelListInv(i)) = &
                            NatOrbMat(SymLabelListInv(j),SymLabelListInv(i)) &
                                        + ( All_abab_RDM(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) ) 
                endif

            endif

        endif

        if((i.ne.a).and.(j.ne.a)) then
            Ind1_1e_aa = ( ( (max(i,a)-2) * (max(i,a)-1) ) / 2 ) + min(i,a)
            Ind2_1e_aa = ( ( (max(j,a)-2) * (max(j,a)-1) ) / 2 ) + min(j,a)

            if((i.ne.j).and.((i.lt.a).and.(j.gt.a)).or.((i.gt.a).and.(j.lt.a))) then
                RDMEnergy_Inst = RDMEnergy_Inst - ( (abba_RDM(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM_Inst) &
                                                    * REAL(TMAT2D(iSpin,jSpin),8) &
                                                    * (1.0_dp / real(NEl - 1,dp)) )
                RDMEnergy1 = RDMEnergy1 - ( (All_abba_RDM(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                    * REAL(TMAT2D(iSpin,jSpin),8) &
                                                    * (1.0_dp / real(NEl - 1,dp)) )

                if((tDiagRDM.or.tPrint1RDM).and.tFinalRDMEnergy) then                                                
                    NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) = &
                                NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) &
                                            - ( All_abba_RDM(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                    * (1.0_dp / real(NEl - 1,dp)) ) 
                endif

                if(Ind1_1e_aa.ne.Ind2_1e_aa) then
                    RDMEnergy_Inst = RDMEnergy_Inst - ( (abba_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM_Inst) &
                                                        * REAL(TMAT2D(iSpin,jSpin),8) &
                                                        * (1.0_dp / real(NEl - 1,dp)) )
                    RDMEnergy1 = RDMEnergy1 - ( (All_abba_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                        * REAL(TMAT2D(iSpin,jSpin),8) &
                                                        * (1.0_dp / real(NEl - 1,dp)) )

                    if((tDiagRDM.or.tPrint1RDM).and.tFinalRDMEnergy) then                                                
                        NatOrbMat(SymLabelListInv(j),SymLabelListInv(i)) = &
                                NatOrbMat(SymLabelListInv(j),SymLabelListInv(i)) &
                                            - ( ( All_abba_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                    * (1.0_dp / real(NEl - 1,dp)) ) / 2.0_dp )
                    endif
                endif
            endif


            if(((i.lt.a).and.(j.lt.a)).or.((i.gt.a).and.(j.gt.a))) then
                Parity_Factor = 1.0_dp
            else
                Parity_Factor = -1.0_dp
            endif
!
            RDMEnergy_Inst = RDMEnergy_Inst + ( (aaaa_RDM(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM_Inst) &
                                                * REAL(TMAT2D(iSpin,jSpin),8) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
            RDMEnergy1 = RDMEnergy1 + ( (All_aaaa_RDM(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                * REAL(TMAT2D(iSpin,jSpin),8) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )

            if((tDiagRDM.or.tPrint1RDM).and.tFinalRDMEnergy) then                                                
                NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) = &
                            NatOrbMat(SymLabelListInv(i),SymLabelListInv(j)) &
                                        + ( All_aaaa_RDM(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 
            endif

            if(Ind1_1e_aa.ne.Ind2_1e_aa) then
                RDMEnergy_Inst = RDMEnergy_Inst + ( (aaaa_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM_Inst) &
                                                * REAL(TMAT2D(jSpin,iSpin),8) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
                RDMEnergy1 = RDMEnergy1 + ( (All_aaaa_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                * REAL(TMAT2D(jSpin,iSpin),8) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )

                if((tDiagRDM.or.tPrint1RDM).and.tFinalRDMEnergy) then                                                
                    NatOrbMat(SymLabelListInv(j),SymLabelListInv(i)) = &
                            NatOrbMat(SymLabelListInv(j),SymLabelListInv(i)) &
                                        + ( All_aaaa_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
                endif
            endif
        endif

    end subroutine calc_1RDM_energy


    subroutine find_nat_orb_occ_numbers()
! Diagonalises the 1-RDM (NatOrbMat), so that after this routine NatOrbMat is the 
! eigenfunctions of the 1-RDM (the matrix transforming the MO's into the NOs).
! This also gets the NO occupation numbers (evaluse) and correlation entropy.
        USE Logging , only : tPrintRODump, tNoNOTransform
        USE SystemData , only : ARR, BRR, G1
        USE RotateOrbsMod , only : FourIndInts, FourIndIntsTag
        USE RotateOrbsData , only : NoOrbs
        USE UMatCache , only : GTID
        implicit none
        integer :: i, j, ierr, Evalues_unit, NatOrbs_unit, jSpat, jInd
        REAL(dp) :: SumDiag, Corr_Entropy, Norm_Evalues
        logical :: tNegEvalue
        CHARACTER(len=*), PARAMETER :: this_routine='find_nat_orb_occ_numbers'

        IF(iProcIndex.eq.0) THEN
            
            ! Diagonalises the 1-RDM.  NatOrbMat goes in as the 1-RDM, comes out as the 
            ! eigenvector of the 1-RDM (the matrix transforming the MO's into the NOs).
            CALL DiagRDM(SumDiag)
            Norm_Evalues = SumDiag/REAL(NEl)

            ! Write out normalised evalues to file and calculate the correlation entropy.
            Corr_Entropy = 0.D0
            Evalues_unit = get_free_unit()
            OPEN(Evalues_unit,file='NO_OCC_NUMBERS',status='unknown')
            WRITE(Evalues_unit,'(A)') '# NORMALISED 1RDM EVALUES (NATURAL ORBITAL OCCUPATION NUMBERS):'
            tNegEvalue = .false.
            do i=1,NoOrbs
                if(tStoreSpinOrbs) then
                    WRITE(Evalues_unit,'(I6,G25.17)') i,Evalues(i)/Norm_Evalues
                    if(Evalues(i).gt.0.D0) then
                        Corr_Entropy = Corr_Entropy - ( abs(Evalues(i)/ Norm_Evalues) &
                                                        * LOG(abs(Evalues(i)/ Norm_Evalues)) )
                    else
                        tNegEvalue = .true.
                    endif
                else
                    WRITE(Evalues_unit,'(I6,G25.17)') i,Evalues(i)/(2.0_dp*Norm_Evalues)
                    if(Evalues(i).gt.0.D0) then
                        Corr_Entropy = Corr_Entropy - (2.0_dp * ( abs(Evalues(i)/(2.0_dp*Norm_Evalues)) &
                                                        * LOG(abs(Evalues(i)/(2.0_dp*Norm_Evalues))) ) )
                    else
                        tNegEvalue = .true.
                    endif
                endif
            enddo
            if(.not.tStoreSpinOrbs) then
                do i=1,SpatOrbs
                    WRITE(Evalues_unit,'(I6,G25.17)') i+SpatOrbs,Evalues(i)/(2.0_dp*Norm_Evalues)
                enddo
            endif
            close(Evalues_unit)
            WRITE(6,'(A20,F30.20)') ' CORRELATION ENTROPY', Corr_Entropy
            WRITE(6,'(A33,F30.20)') ' CORRELATION ENTROPY PER ELECTRON', Corr_Entropy / real(NEl,dp) 
            if(tNegEvalue) write(6,'(A)') ' WARNING: Negative NO occupation numbers found.'

            ! Write out the evectors to file.
            ! This is the matrix that transforms the molecular orbitals into the natural orbitals.
            ! Evalue(i) corresponds to Evector NatOrbsMat(1:nBasis,i)
            ! We just want the Evalues in the same order as above, but the 1:nBasis part (corresponding 
            ! to the molecular orbitals), needs to refer to the actual orbital labels.
            ! Want these orbitals to preferably be in order, run through the orbital, need the position 
            ! to find the corresponding NatOrbs element, use SymLabelListInv
            if(.not.tNoNOTransform) then
                NatOrbs_unit = get_free_unit()
                OPEN(NatOrbs_unit,file='NO_TRANSFORM',status='unknown')
                write(NatOrbs_unit,'(2A6,A20)') '#   MO','NO','Transform Coeff'
                ! write out in terms of spin orbitals, all alpha then all beta.
                do i = 1, NoOrbs
                    do j = 1, nBasis
                        ! Here i corresponds to the natural orbital, and j to the molecular orbital.
                        ! i is actually the spin orbital in this case.
                        if(tStoreSpinOrbs) then
                            jInd = j
                        else
                            if(mod(j,2).ne.0) then
                                jInd = gtID(j)
                            else
                                CYCLE
                            endif
                        endif
                        if(NatOrbMat(SymLabelListInv(jInd),i).ne.0.0_dp) &
                            write(NatOrbs_unit,'(2I6,G25.17)') j,i,NatOrbMat(SymLabelListInv(jInd),i)
                    enddo
                enddo
                if(.not.tStoreSpinOrbs) then
                    do i = 1, SpatOrbs
                        do j = 2, nBasis, 2
                            ! Here i corresponds to the natural orbital, and j to the molecular orbital.
                            ! i is actually the spin orbital in this case.
                            jSpat = gtID(j)
                            if(NatOrbMat(SymLabelListInv(jSpat),i).ne.0.0_dp) &
                                write(NatOrbs_unit,'(2I6,G25.17)') j,i+SpatOrbs,&
                                                                NatOrbMat(SymLabelListInv(jSpat),i)
                        enddo
                    enddo
                endif
                close(NatOrbs_unit)
            endif

            if(tPrintRODump) then                

                ALLOCATE(FourIndInts(NoOrbs,NoOrbs,NoOrbs,NoOrbs),stat=ierr)
                CALL LogMemAlloc('FourIndInts',(NoOrbs**4),8,this_routine,&
                                                        FourIndIntsTag,ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating FourIndInts array,')

! Then, transform2ElInts
                WRITE(6,*) ''
                WRITE(6,*) 'Transforming the four index integrals'
                CALL Transform2ElIntsMemSave_RDM()

                WRITE(6,*) 'Re-calculating the fock matrix'
                CALL CalcFOCKMatrix_RDM()

                WRITE(6,*) 'Refilling the UMAT and TMAT2D'
! The ROFCIDUMP is also printed out in here.        
                CALL RefillUMATandTMAT2D_RDM()        

                CALL FLUSH(6)

                CALL WRITEBASIS(6,G1,NoOrbs,ARR,BRR)

            ENDIF
        ENDIF

    end subroutine find_nat_orb_occ_numbers


    SUBROUTINE DiagRDM(SumTrace)
! The diagonalisation routine reorders the orbitals in such a way that the 
! corresponding orbital labels are lost.
! In order to keep the spin and spatial symmetries, each symmetry must be fed into 
! the diagonalisation routine separately.
! The best way to do this is to order the orbitals so that all the alpha orbitals 
! follow all the beta orbitals, with the occupied orbitals first, in terms of symmetry, 
! and the virtual second, also ordered by symmetry.
! This gives us flexibility w.r.t rotating only the occupied or only virtual and 
! looking at high spin states.
        IMPLICIT NONE
        REAL(dp) , intent(out) :: SumTrace
        REAL(dp) :: SumDiagTrace
        REAL(dp) , ALLOCATABLE :: WORK2(:),EvaluesSym(:),NOMSym(:,:)
        INTEGER :: ierr,i,j,spin,Sym,LWORK2,WORK2Tag,SymStartInd,NoSymBlock
        INTEGER :: EvaluesSymTag,NOMSymTag,k,MaxSym
        LOGICAL :: tDiffSym
        CHARACTER(len=*), PARAMETER :: this_routine='DiagRDM'

! Test that we're not breaking symmetry.
! And calculate the trace at the same time.
        SumTrace=0.D0
        do i=1,NoOrbs
            do j=1,NoOrbs
                tDiffSym = .false.
                if(tStoreSpinOrbs) then
                    IF((INT(G1(SymLabelList2(i))%sym%S,4).ne.&
                        INT(G1(SymLabelList2(j))%sym%S,4))) tDiffSym = .true.
                else
                    IF((INT(G1(2*SymLabelList2(i))%sym%S,4).ne.&
                        INT(G1(2*SymLabelList2(j))%sym%S,4))) tDiffSym = .true.
                endif
                if(tDiffSym) then
                    IF(ABS(NatOrbMat(i,j)).ge.1.0E-15) THEN
                        WRITE(6,'(6A8,A20)') 'i','j','Label i','Label j','Sym i',&
                                                                'Sym j','Matrix value'
                        if(tStoreSpinOrbs) then                                                              
                            WRITE(6,'(6I3,F40.20)') i,j,SymLabelList2(i),SymLabelList2(j),&
                                    INT(G1(SymLabelList2(i))%sym%S,4),&
                                    INT(G1(SymLabelList2(j))%sym%S,4),NatOrbMat(i,j)
                        else
                            WRITE(6,'(6I3,F40.20)') i,j,SymLabelList2(i),SymLabelList2(j),&
                                    INT(G1(2*SymLabelList2(i))%sym%S,4),&
                                    INT(G1(2*SymLabelList2(j))%sym%S,4),NatOrbMat(i,j)
                        endif
                        IF(tUseMP2VarDenMat) THEN
                            WRITE(6,*) '**WARNING** - There is a non-zero NatOrbMat &
                            &value between orbitals of different symmetry.'
                            WRITE(6,*) 'These elements will be ignored, and the symmetry &
                            &maintained in the final transformation matrix.'
                        ELSE
                            write(6,*) 'k,SymLabelList2(k),SymLabelListInv(k)'
                            do k = 1,NoOrbs
                                write(6,*) k,SymLabelList2(k),SymLabelListInv(k)
                            enddo
                            call flush(6)
                            CALL Stop_All(this_routine,'Non-zero NatOrbMat value between &
                            &different symmetries.')
                        ENDIF
                    ENDIF
                    NatOrbMat(i,j)=0.D0
                ENDIF
            enddo
            SumTrace=SumTrace+NatOrbMat(i,i)
        enddo

        write(6,*) ''
        WRITE(6,*) 'Calculating eigenvectors and eigenvalues of NatOrbMat'
        CALL FLUSH(6)

! If we want to maintain the symmetry, we cannot have all the orbitals jumbled up when the 
! diagonaliser reorders the eigenvectors.
! Must instead feed each symmetry block in separately.
! This means that although the transformed orbitals are jumbled within the symmetry blocks, 
! the symmetry labels are all that are relevant and these are unaffected.
        Sym=0
        LWORK2=-1
        if(tStoreSpinOrbs) then
            MaxSym = 15
        else
            MaxSym = 7
        endif
        do while (Sym.le.MaxSym)

            NoSymBlock=SymLabelCounts2(2,Sym+1)

            SymStartInd=SymLabelCounts2(1,Sym+1)-1
            ! This is one less than the index that the symmetry starts, so that when we 
            ! run through i=1,..., we can start at SymStartInd+i.

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

!                WRITE(6,*) '*****'
!                WRITE(6,*) 'Symmetry ',Sym, 'with spin ',spin,' has ',&
!                                                        NoSymBlock,' orbitals.'
!                WRITE(6,*) 'The NatOrbMat for this symmetry block is '
!                do i=1,NoSymBlock
!                    do j=1,NoSymBlock
!                        WRITE(6,'(F20.10)',advance='no') NOMSym(j,i)
!                    enddo
!                    WRITE(6,*) ''
!                enddo

                CALL DSYEV('V','L',NoSymBlock,NOMSym,NoSymBlock,EvaluesSym,&
                                                                WORK2,LWORK2,ierr)
                ! NOMSym goes in as the original NOMSym, comes out as the 
                ! eigenvectors (Coefficients).
                ! EvaluesSym comes out as the eigenvalues in ascending order.

!                WRITE(6,*) 'After diagonalization, the e-vectors (diagonal elements) of this matrix are ,'
!                do i=1,NoSymBlock
!                    WRITE(6,'(F20.10)',advance='no') EvaluesSym(NoSymBlock-i+1)
!                enddo
!                WRITE(6,*) ''
!                WRITE(6,*) 'These go from orbital ,',SymStartInd+1,' to '&
!                                                            ,SymStartInd+NoSymBlock

                do i=1,NoSymBlock
                    Evalues(SymStartInd+i)=EvaluesSym(NoSymBlock-i+1)
                enddo

                ! CAREFUL if eigenvalues are put in ascending order, this may not be 
                ! correct, with the labelling system.
                ! may be better to just take coefficients and transform TMAT2DRot 
                ! in transform2elints.
                ! a check that comes out as diagonal is a check of this routine anyway.

!                WRITE(6,*) 'The eigenvectors (coefficients) for symmetry block ',Sym
!                WRITE(6,'(I20)',advance='no') 0
!                do i = 1, NoSymBlock
!                    WRITE(6,'(I20)',advance='no') SymLabelList2(SymStartInd+i)
!                enddo
!                write(6,*) ''
!                do i=1,NoSymBlock
!                    write(6,'(I20)',advance='no') SymLabelList2(SymStartInd+i)
!                    do j=1,NoSymBlock
!                        WRITE(6,'(F20.10)',advance='no') NOMSym(j,NoSymBlock-i+1)
!                    enddo
!                    WRITE(6,*) ''
!                enddo
         
                do j=1,NoSymBlock
                    do i=1,NoSymBlock
                        NatOrbMat(SymStartInd+i,SymStartInd+j)=NOMSym(i,NoSymBlock-j+1)
                    enddo
                enddo
                ! Directly fill the coefficient matrix with the eigenvectors from 
                ! the diagonalization.

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
!                WRITE(6,*) '*****'
!                WRITE(6,*) 'Symmetry ',Sym,' has only one orbital.'
!                WRITE(6,*) 'Copying diagonal element ,',SymStartInd+1,'to NatOrbMat'
            ENDIF

            Sym=Sym+1
        enddo

        WRITE(6,*) 'Matrix diagonalised'
        CALL FLUSH(6)

        SumDiagTrace=0.D0
        do i=1,NoOrbs
            SumDiagTrace=SumDiagTrace+Evalues(i)
        enddo
        IF((ABS(SumDiagTrace-SumTrace)).gt.1.D0) THEN
            WRITE(6,*) 'Sum of diagonal NatOrbMat elements : ',SumTrace
            WRITE(6,*) 'Sum of eigenvalues : ',SumDiagTrace
            WRITE(6,*) 'WARNING : &
            &The trace of the 1RDM matrix before diagonalisation is '
            write(6,*) 'not equal to that after.'
        ENDIF

! The MO's still correspond to SymLabelList2, but the NO's are slightly jumbled now.

!        call OrderNatOrbMat()

    END SUBROUTINE DiagRDM


    SUBROUTINE OrderNatOrbMat()
        USE sort_mod
        USE Logging , only : tTruncRODump
        IMPLICIT NONE
        INTEGER :: spin,i,j,ierr,StartSort,EndSort
        CHARACTER(len=*), PARAMETER :: this_routine='OrderRDM'
        INTEGER , ALLOCATABLE :: NatOrbMatTemp(:,:), SymLabelList3(:), SymLabelListTemp(:), EvaluesTemp(:)
        INTEGER :: NatOrbMatTempTag, SymLabelList3Tag, Orb, New_Pos
        
! Here, if symmetry is kept, we are going to have to reorder the eigenvectors 
! according to the size of the eigenvalues, while taking the orbital labels 
! (and therefore symmetries) with them. This will be put back into MP2VDM from MP2VDMTemp.

! Want to reorder the eigenvalues from largest to smallest, taking the eigenvectors 
! with them and the symmetry as well.  
! If using spin orbitals, do this for the alpha spin and then the beta.

        ALLOCATE(NatOrbMatTemp(NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('NatOrbMatTemp',NoOrbs**2,8,&
                            'OrderNatOrbMat',NatOrbMatTempTag,ierr)
        ALLOCATE(SymLabelList3(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList3',NoOrbs,8,&
                            'OrderNatOrbMat',SymLabelList3Tag,ierr)

        ALLOCATE(SymLabelListTemp(NoOrbs),stat=ierr)
        ALLOCATE(EvaluesTemp(NoOrbs),stat=ierr)

        do i = 1, NoOrbs
            SymLabelList3(i) = SymLabelList2(i)
        enddo

        StartSort=1
        EndSort=SpatOrbs

        call sort (EValues(startSort:endSort), &
                   NatOrbMat(1:NoOrbs, startSort:endSort), &
                   symLabelList2(startSort:endSort))

        if(tStoreSpinOrbs) then                  
            StartSort=SpatOrbs + 1
            EndSort=nBasis

            call sort (EValues(startSort:endSort), &
                       NatOrbMat(1:NoOrbs, startSort:endSort), &
                       symLabelList2(startSort:endSort))

        endif                       

        write(6,*) 'after sort'
        do i = 1, NoOrbs
            write(6,*) i,SymLabelList2(i),SymLabelList3(i)
        enddo

        SymLabelListTemp(:) = SymLabelList2(:)
        EvaluesTemp(:) = Evalues(:)
        do i = 1, NoOrbs
            SymLabelList2(NoOrbs - i + 1) = SymLabelListTemp(i)
            Evalues(NoOrbs - i + 1) = EvaluesTemp(i)
        enddo

        SymLabelListInv(:)=0   
        do i=1,NoOrbs
            SymLabelListInv(SymLabelList2(i))=i
        enddo

        NatOrbMatTemp(:,:) = NatOrbMat(:,:)
        NatOrbMat(:,:) = 0.0_dp

        do i=1,NoOrbs
            do j = 1, NoOrbs
                Orb = SymLabelList3(j)      
                !This is the orbital the old position corresponds to.
                New_Pos = SymLabelListInv(Orb)   
                !This is the new position that orbital should go into.
                NatOrbMat(NoOrbs-New_Pos+1,NoOrbs-i+1)=NatOrbMatTemp(j,i)
            enddo
        enddo


        DEALLOCATE(NatOrbMatTemp)
        DEALLOCATE(SymLabelList3)


    END SUBROUTINE OrderNatOrbMat

   
! This is an M^5 transform, which transforms all the two-electron integrals 
! into the new basis described by the Coeff matrix.
! This is v memory inefficient and currently does not use any spatial 
! symmetry information.
    SUBROUTINE Transform2ElIntsMemSave_RDM()
        USE RotateOrbsMod , only : FourIndInts
        implicit none
        INTEGER :: i,j,k,l,a,b,g,d,ierr,Temp4indintsTag,a2,d2,b2,g2
        REAL(dp) , ALLOCATABLE :: Temp4indints(:,:)
        CHARACTER(len=*), PARAMETER :: this_routine='Transform2ElIntsMemSave_RDM'
#ifdef __CMPLX
        call stop_all('Transform2ElIntsMemSave_RDM', &
                    'Rotating orbitals not implemented for complex orbitals.')
#endif
        
!Zero arrays from previous transform

        ALLOCATE(Temp4indints(NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('Temp4indints',NoOrbs**2,8,&
                            'Transform2ElIntsMemSave_RDM',Temp4indintsTag,ierr)
        IF(ierr.ne.0) CALL Stop_All('Transform2ElIntsMemSave_RDM',&
                                    'Problem allocating memory to Temp4indints.')
 
        FourIndInts(:,:,:,:)=0.D0

! **************
! Calculating the two-transformed, four index integrals.

! The untransformed <alpha beta | gamma delta> integrals are found from 
! UMAT(UMatInd(i,j,k,l,0,0)
! All our arrays are in spin orbitals - if tStoreSpinOrbs is false, 
! UMAT will be in spatial orbitals - need to account for this.

! Running through 1,NoOrbs - the actual orbitals corresponding to that index 
! are given by SymLabelList2

        do b=1,NoOrbs
            b2 = SymLabelList2(b)
            do d=1,NoOrbs
                d2 = SymLabelList2(d)
                do a=1,NoOrbs
                    a2 = SymLabelList2(a)
                    do g=1,NoOrbs
                        g2 = SymLabelList2(g)

! UMatInd in physical notation, but FourIndInts in chemical 
! (just to make it more clear in these transformations).
! This means that here, a and g are interchangable, and so are b and d.
                        FourIndInts(a,g,b,d)=REAL(UMAT(UMatInd(a2,b2,g2,d2,0,0)),8)
                    enddo
                enddo

                Temp4indints(:,:)=0.D0
                CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,NatOrbMat(:,:),NoOrbs,&
                            FourIndInts(1:NoOrbs,1:NoOrbs,b,d),NoOrbs,0.0,&
                            Temp4indints(1:NoOrbs,1:NoOrbs),NoOrbs)
                ! Temp4indints(i,g) comes out of here, so to transform g to k, 
                ! we need the transpose of this.

                CALL DGEMM('T','T',NoOrbs,NoOrbs,NoOrbs,1.0,NatOrbMat(:,:),NoOrbs,&
                            Temp4indints(1:NoOrbs,1:NoOrbs),NoOrbs,0.0,&
                            FourIndInts(1:NoOrbs,1:NoOrbs,b,d),NoOrbs)
                ! Get Temp4indits02(i,k)
            enddo
        enddo
        
! Calculating the 3 transformed, 4 index integrals. 01=a untransformed,02=b,03=g,04=d
        do i=1,NoOrbs
            do k=1,NoOrbs

                Temp4indints(:,:)=0.D0
                CALL DGEMM('T','N',NoOrbs,NoOrbs,NoOrbs,1.0,NatOrbMat(:,:),NoOrbs,&
                            FourIndInts(i,k,1:NoOrbs,1:NoOrbs),NoOrbs,0.0,&
                            Temp4indints(1:NoOrbs,1:NoOrbs),NoOrbs)

                CALL DGEMM('T','T',NoOrbs,NoOrbs,NoOrbs,1.0,NatOrbMat(:,:),&
                            NoOrbs,Temp4indints(1:NoOrbs,1:NoOrbs),NoOrbs,0.0,&
                            FourIndInts(i,k,1:NoOrbs,1:NoOrbs),NoOrbs)
            enddo
        enddo

        DEALLOCATE(Temp4indints)
        CALL LogMemDeAlloc('Transform2ElIntsMemSave_RDM',Temp4indintsTag)
 
    END SUBROUTINE Transform2ElIntsMemSave_RDM


    SUBROUTINE CalcFOCKMatrix_RDM()
! Calculate the fock matrix in the nat orb basis.    
        USE SystemData , only : nBasis
        USE Logging , only : tRDMonfly
        implicit none
        INTEGER :: i,j,k,l,a,b,ierr,ArrDiagNewTag
        REAL(dp) :: FOCKDiagSumHF,FOCKDiagSumNew
        CHARACTER(len=*) , PARAMETER :: this_routine='CalcFOCKMatrix_RDM'
        REAL(dp) , ALLOCATABLE :: ArrDiagNew(:)

! This subroutine calculates and writes out the fock matrix for the transformed orbitals.
! ARR is originally the fock matrix in the HF basis.
! ARR(:,1) - ordered by energy, ARR(:,2) - ordered by spin-orbital index.

! When transforming the orbitals into approximate natural orbitals, we want to save memory, 
! so don't bother calculating the whole matrix, just the diagonal elements 
! that we actually need.

        ALLOCATE(ArrDiagNew(nBasis),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating ArrDiagNew array,')
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

!        WRITE(6,*) 'Sum of the fock matrix diagonal elements in the HF basis set = ',&
!                                                                        FOCKDiagSumHF
! Then calculate the fock matrix in the transformed basis, 
! and the sum of the new diagonal elements.
! Our Arr in spin orbitals.
!        do j=1,NoOrbs
!            ArrNew(j,j)=Arr(2*j,2)
!        enddo

        FOCKDiagSumNew=0.D0
        do j=1,NoOrbs
            l=SymLabelList2(j)
            if(tStoreSpinOrbs) then
                ArrDiagNew(l) = 0.D0
            else
                ArrDiagNew(2*l) = 0.D0
                ArrDiagNew((2*l)-1) = 0.D0
            endif
            do a=1,NoOrbs
                b=SymLabelList2(a)
                if(tStoreSpinOrbs) then
                    ArrDiagNew(l)=ArrDiagNew(l)+(NatOrbMat(a,j)*ARR(b,2)*NatOrbMat(a,j))
                else
                    ArrDiagNew(2*l)=ArrDiagNew(2*l)+(NatOrbMat(a,j)*ARR(2*b,2)*NatOrbMat(a,j))
                    ArrDiagNew((2*l)-1)=ArrDiagNew((2*l)-1)+(NatOrbMat(a,j)*ARR((2*b)-1,2)*NatOrbMat(a,j))
                endif
            enddo
            if(tStoreSpinOrbs) then
                FOCKDiagSumNew=FOCKDiagSumNew+(ArrDiagNew(l))
            else
                FOCKDiagSumNew=FOCKDiagSumNew+(ArrDiagNew(2*l))
                FOCKDiagSumNew=FOCKDiagSumNew+(ArrDiagNew((2*l)-1))
            endif
        enddo
        ! If we are truncation the virtual space, only the unfrozen entries will 
        ! be transformed.

!        WRITE(6,*) 'Sum of the fock matrix diagonal elements in the transformed &
!                        &basis set = ',FOCKDiagSumNew

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
       
! Refill ARR(:,1) (ordered in terms of energies), and ARR(:,2) 
! (ordered in terms of orbital number).
! ARR(:,2) needs to be ordered in terms of symmetry and then energy (like SymLabelList), 
! so currently this ordering will not be correct when reading in qchem INTDUMPS as the 
! orbital number ordering is by energy.

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

    ENDSUBROUTINE CalcFOCKMatrix_RDM

    SUBROUTINE RefillUMATandTMAT2D_RDM()
! UMat is in spin or spatial orbitals, TMAT2D only spin.
! This routine refills these to more easily write out the ROFCIDUMP, and originally 
! to be able to continue a calculation (although I doubt this works at the moment).
        USE RotateOrbsMod , only : FourIndInts, PrintROFCIDUMP, PrintRepeatROFCIDUMP
        implicit none
        INTEGER :: l,k,j,i,a,b,g,d,c,nBasis2,TMAT2DPartTag,ierr
        REAL(dp) :: NewTMAT
        REAL(dp) , ALLOCATABLE :: TMAT2DPart(:,:)
        CHARACTER(len=*), PARAMETER :: this_routine='RefillUMATandTMAT2D_RDM'
#ifdef __CMPLX
        call stop_all('RefillUMATandTMAT2D_RDM', &
                    'Rotating orbitals not implemented for complex orbitals.')
#endif

        ! TMAT2D is always in spin orbitals.
        ALLOCATE(TMAT2DPart(nBasis,nBasis),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating TMAT2DPart array,')
        CALL LogMemAlloc('TMAT2DPart',nBasis*nBasis,8,&
                                    'RefillUMAT_RDM',TMAT2DPartTag,ierr)
        TMAT2DPart(:,:)=0.D0

! Make the UMAT elements the four index integrals.  
! These are calculated by transforming the HF orbitals using the coefficients 
! that have been found
        do l=1,NoOrbs
            d=SymLabelList2(l)
            do k=1,NoOrbs
                b=SymLabelList2(k)
                do j=1,NoOrbs
                    g=SymLabelList2(j)
                    do i=1,NoOrbs
                        a=SymLabelList2(i)
!The FourIndInts are in chemical notation, the UMatInd in physical.                            
                        UMAT(UMatInd(a,b,g,d,0,0))=FourIndInts(i,j,k,l)
                    enddo
                enddo
            enddo
        enddo

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

        do a=1,nBasis
            do k=1,NoOrbs
                i=SymLabelList2(k)
                NewTMAT=0.D0
                do b=1,NoOrbs
                    d=SymLabelList2(b)
                    if(tStoreSpinOrbs) then
                        NewTMAT=NewTMAT+(NatOrbMat(b,k)*REAL(TMAT2D(d,a),8))
                    else
                        NewTMAT=NewTMAT+(NatOrbMat(b,k)*REAL(TMAT2D(2*d,a),8))
                    endif
                enddo
                if(tStoreSpinOrbs) then
                    TMAT2DPart(i,a)=NewTMAT
                else
                    if(mod(a,2).eq.0) then
                        TMAT2DPart(2*i,a)=NewTMAT
                    else
                        TMAT2DPart((2*i)-1,a)=NewTMAT
                    endif
                endif
            enddo
        enddo

        do k=1,nBasis
            do l=1,NoOrbs
                j=SymLabelList2(l)
                NewTMAT=0.D0
                do a=1,NoOrbs
                    c=SymLabelList2(a)
                    if(tStoreSpinOrbs) then
                        NewTMAT=NewTMAT+(NatOrbMat(a,l)*TMAT2DPart(k,c))
                    else
                        NewTMAT=NewTMAT+(NatOrbMat(a,l)*TMAT2DPart(k,2*c))
                    endif
                enddo
                if(tStoreSpinOrbs) then
                    TMAT2D(k,j)=NewTMAT
                else
                    if(mod(k,2).eq.0) then
                        TMAT2D(k,2*j)=NewTMAT
                    else
                        TMAT2D(k,(2*j)-1)=NewTMAT
                    endif
                endif
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

        CALL PrintROFCIDUMP_RDM()

    ENDSUBROUTINE RefillUMATandTMAT2D_RDM


    SUBROUTINE PrintROFCIDUMP_RDM()
!This prints out a new FCIDUMP file in the same format as the old one.
        implicit none
        INTEGER :: i,j,k,l,iunit

!        PrintROFCIDUMP_Time%timer_name='PrintROFCIDUMP'
!        CALL set_timer(PrintROFCIDUMP_Time,30)

        iunit = get_free_unit()
        OPEN(iunit,FILE='ROFCIDUMP',STATUS='unknown')
        
        WRITE(iunit,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',NoOrbs,&
                                                ',NELEC=',NEl,',MS2=',LMS,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,NoOrbs
            IF(tStoreSpinOrbs) THEN
                WRITE(iunit,'(I1,A1)',advance='no') INT(G1(i)%sym%S,4)+1,','
            ELSE
                WRITE(iunit,'(I1,A1)',advance='no') INT(G1(i*2)%sym%S,4)+1,','
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
            do j=1,NoOrbs
                do l=1,j
                    ! Potential to put symmetry in here, have currently taken it out, 
                    ! because when we're only printing non-zero values, it is kind 
                    ! of unnecessary - although it may be used to speed things up.
                    do k=1,i
! UMatInd is in physical notation <ij|kl>, but the indices printed in the FCIDUMP 
! are in chemical notation (ik|jl).
                        IF((ABS(REAL(UMat(UMatInd(i,j,k,l,0,0)),8))).ne.0.D0) &
                                WRITE(iunit,'(F21.12,4I3)') &
                                REAL(UMat(UMatInd(i,j,k,l,0,0)),8),i,k,j,l 
 
                    enddo
                enddo
           enddo
        enddo

! TMAT2D stored as spin orbitals
        do i=1,NoOrbs
            ! Symmetry?
            do k=1,NoOrbs
                IF(tStoreSpinOrbs) THEN
                    IF((REAL(TMAT2D(i,k),8)).ne.0.D0) WRITE(iunit,'(F21.12,4I3)') &
                                                        REAL(TMAT2D(i,k),8),i,k,0,0
                ELSE
                    IF((REAL(TMAT2D(2*i,2*k),8)).ne.0.D0) WRITE(iunit,'(F21.12,4I3)') &
                                                        REAL(TMAT2D(2*i,2*k),8),i,k,0,0
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
! This routine just deallocates the arrays allocated in InitRDM.
! If the NECI calculation softexits before the RDMs start to fill, this is all that 
! is called at the end.
        USE RotateOrbsMod , only : SymLabelList2,SymLabelListInv,&
                                   SymLabelListInvTag,SymLabelList2Tag,&
                                   FourIndInts, FourIndIntsTag, &
                                   SymLabelCounts2, SymLabelCounts2Tag
        USE Logging , only : tDiagRDM, tPrintRODump
        implicit none
        CHARACTER(len=*), PARAMETER :: this_routine='DeallocateRDM'

        IF(tExplicitAllRDM) THEN

! This array contains the initial positions of the single excitations for 
! each processor.
            DEALLOCATE(Sing_InitExcSlots)
 
! This array contains the current position of the single excitations as 
! they're added.
            DEALLOCATE(Sing_ExcList)

! This array actually contains the single excitations in blocks of the 
! processor they will be sent to.        
            DEALLOCATE(Sing_ExcDjs)
            CALL LogMemDeAlloc(this_routine,Sing_ExcDjsTag)

            DEALLOCATE(Sing_ExcDjs2)
            CALL LogMemDeAlloc(this_routine,Sing_ExcDjs2Tag)


            IF(RDMExcitLevel.ne.1) THEN
! This array contains the initial positions of the single excitations for 
! each processor.
                DEALLOCATE(Doub_InitExcSlots)
 
! This array contains the current position of the single excitations as 
! they're added.
                DEALLOCATE(Doub_ExcList)

! This array actually contains the single excitations in blocks of the 
! processor they will be sent to.        
                DEALLOCATE(Doub_ExcDjs)
                CALL LogMemDeAlloc(this_routine,Doub_ExcDjsTag)
     
                DEALLOCATE(Doub_ExcDjs2)
                CALL LogMemDeAlloc(this_routine,Doub_ExcDjs2Tag)
            ENDIF

        ELSEIF(.not.tHF_Ref_Explicit) THEN

            DEALLOCATE(Spawned_Parents)
            CALL LogMemDeAlloc(this_routine,Spawned_ParentsTag)

            DEALLOCATE(Spawned_Parents_Index)
            CALL LogMemDeAlloc(this_routine,Spawned_Parents_IndexTag)

        ENDIF

        IF(RDMExcitLevel.eq.1) THEN

            DEALLOCATE(NatOrbMat)
            CALL LogMemDeAlloc(this_routine,NatOrbMatTag)

            IF((iProcIndex.eq.0).and.tDiagRDM) THEN
                DEALLOCATE(Evalues)
                CALL LogMemDeAlloc(this_routine,EvaluesTag)

                IF(tPrintRODump) THEN
                    DEALLOCATE(FourIndInts)
                    CALL LogMemDeAlloc(this_routine,FourIndIntsTag)
                ENDIF

            ENDIF

            DEALLOCATE(SymLabelCounts2)
            CALL LogMemDeAlloc(this_routine,SymLabelCounts2Tag)

            DEALLOCATE(SymLabelList2)
            CALL LogMemDeAlloc(this_routine,SymLabelList2Tag)

            DEALLOCATE(SymLabelListInv)
            CALL LogMemDeAlloc(this_routine,SymLabelListInvTag)

        ELSE

            DEALLOCATE(aaaa_RDM)
            CALL LogMemDeAlloc(this_routine,aaaa_RDMTag)
            DEALLOCATE(abab_RDM)
            CALL LogMemDeAlloc(this_routine,abab_RDMTag)
            DEALLOCATE(abba_RDM)
            CALL LogMemDeAlloc(this_routine,abba_RDMTag)

            IF(iProcIndex.eq.0) THEN
                DEALLOCATE(All_aaaa_RDM)
                CALL LogMemDeAlloc(this_routine,All_aaaa_RDMTag)
                DEALLOCATE(All_abab_RDM)
                CALL LogMemDeAlloc(this_routine,All_abab_RDMTag)
                DEALLOCATE(All_abba_RDM)
                CALL LogMemDeAlloc(this_routine,All_abba_RDMTag)
 
                if(tDiagRDM.or.tPrint1RDM) then
                    DEALLOCATE(NatOrbMat)
                    CALL LogMemDeAlloc(this_routine,NatOrbMatTag)
                endif
            ENDIF

        ENDIF

    END SUBROUTINE DeallocateRDM


!    subroutine sum_in_spin_proj(i,j,Ind1,Norm_2RDM,Tot_Spin_Projection)
!        implicit none
!        integer , intent(in) :: i, j, Ind1
!        real(dp) , intent(in) :: Norm_2RDM
!        real(dp) , intent(inout) :: Tot_Spin_Projection
!        integer :: Ind2
!        real(dp) :: ParityFactor, Mult_Fac1, Mult_Fac2
!
!        ! Total Spin Projection: 
!        ! Sum_i_j [ 2RDM(ia,ja;ia,ja) + 2RDM(ib,jb;ib,jb) 
!        !           - 2 * 2RDM(ia,jb;ia,jb) - 4 * 2RDM(ia,jb,ja,ib) ] + (3 * NEl) = 4S(S+1)
!
!        ! i and j both alpha (even i and j).
!        if((mod(i,2).eq.0).and.(mod(j,2).eq.0)) then
!            Tot_Spin_Projection = Tot_Spin_Projection & 
!                                + (2.D0 * AllTwoElRDM(Ind1, Ind1) * Norm_2RDM )  
!
!        ! i and j both beta (odd i and j).
!        elseif((mod(i,2).ne.0).and.(mod(j,2).ne.0)) then
!            Tot_Spin_Projection = Tot_Spin_Projection & 
!                                + (2.D0 * AllTwoElRDM(Ind1, Ind1) * Norm_2RDM )  
!
!        ! i alpha (even), j beta (odd).                                                              
!        elseif((mod(i,2).eq.0).and.(mod(j,2).ne.0)) then
!
!            if(j.eq.(i-1)) then
!                ! this is the case where spat i = spat j
!                Mult_Fac1 = 2.D0
!                Mult_Fac2 = 4.D0
!            else
!                Mult_Fac1 = 4.D0
!                Mult_Fac2 = 8.D0
!            endif
!
!            ! adding in ialpha jbeta ialpha jbeta
!            Tot_Spin_Projection = Tot_Spin_Projection &
!                                    - (Mult_Fac1 *  AllTwoElRDM(Ind1, Ind1) * Norm_2RDM)  
!
!            ParityFactor = 1.D0
!            ! we find the index based on i < j, in reality this will not always be the case.
!            ! if i > j, times the parity by -1.
!            if(max(i,j).eq.i) ParityFactor = ParityFactor * (-1.D0)
!
!            ! adding in ialpha jbeta jalpha ibeta                                            
!            ! Ind2 is for jalpha ibeta
!
!            ! i is even (alpha) - want same orbital but beta = i-1                                    
!            ! j is odd (beta) - want same orbital but alpha = j+1                                    
! 
!            ! i currently alpha (2), j currently beta (1)
!            ! want i beta (1), j alpha (2)
!            ! but the parity will be wrong
!
!            Ind2 = ( ( (max(i-1,j+1)-2) * (max(i-1,j+1)-1) ) / 2 ) + min(i-1,j+1)
!
!            ! Expect now that j+1 < i-1, check if this is actually true.
!            if(max(i-1,j+1).eq.(j+1)) ParityFactor = ParityFactor * (-1.D0)
!
!            Tot_Spin_Projection = Tot_Spin_Projection &
!                                    - (Mult_Fac2 * AllTwoElRDM(Ind1, Ind2) * Norm_2RDM * ParityFactor)
!        endif
!
!    end subroutine

!
!    SUBROUTINE Average_Spins_and_Sum_1e_Norms(Trace_1RDM)
!! This will only be called by the root.
!! It takes the NatOrbMat and averages each of the elements that should be equal by spin.
!! If HPHF is on, these things should already be true and this is not called.
!        implicit none
!        real(dp) , intent(out) :: Trace_1RDM
!        real(dp) :: Entry_bb, Entry_aa
!        integer :: i, a
!
!        Trace_1RDM = 0.D0
!
!        do i = 1, nBasis - 1, 2
!
!            ! ialpha -> ialpha = ibeta -> ibeta
!            Entry_bb = NatOrbMat(SymLabelListInv(i),SymLabelListInv(i))
!            Entry_aa = NatOrbMat(SymLabelListInv(i+1),SymLabelListInv(i+1))
!
!            NatOrbMat(SymLabelListInv(i),SymLabelListInv(i)) = &
!                                        ( Entry_bb + Entry_aa ) / 2.D0 
!
!            NatOrbMat(SymLabelListInv(i+1),SymLabelListInv(i+1)) = &
!                                        ( Entry_bb + Entry_aa ) / 2.D0 
!
!            Trace_1RDM = Trace_1RDM + NatOrbMat(SymLabelListInv(i),SymLabelListInv(i)) &
!                                    + NatOrbMat(SymLabelListInv(i+1),SymLabelListInv(i+1))
!
!            do a = 1, nBasis - 1, 2
!
!                ! ialpha -> jalpha = ibeta -> jbeta
!                Entry_bb = NatOrbMat(SymLabelListInv(i),SymLabelListInv(a))
!                Entry_aa = NatOrbMat(SymLabelListInv(i+1),SymLabelListInv(a+1))
!
!                NatOrbMat(SymLabelListInv(i),SymLabelListInv(a)) = &
!                                        ( Entry_bb + Entry_aa ) / 2.D0 
!
!                NatOrbMat(SymLabelListInv(i+1),SymLabelListInv(a+1)) = &
!                                        ( Entry_bb + Entry_aa ) / 2.D0 
!            enddo
!        enddo
!
!!    
!    END SUBROUTINE Average_Spins_and_Sum_1e_Norms
!
!
!
!    SUBROUTINE Average_Spins_and_Sum_2e_Norms(Trace_2RDM_Inst, Trace_2RDM)
!! This will only be called by the root.
!! It takes the AllTwoElRDM, and averages each of the elements that should be equal by spin.
!! If HPHF is on, these things should already be true and this is not called.
!
!! Elements that should be equal;
!! < a a | a a > = < b b | b b > 
!! < a b | a b > = < b a | b a >
!! < a b | b a > = < b a | a b >
!        implicit none
!        real(dp) , intent(out) :: Trace_2RDM_Inst, Trace_2RDM
!        real(dp) :: Entry_abab, Entry_baba, Sign_abab, Sign_baba
!        real(dp) :: Entry_abba, Entry_baab, Sign_abba, Sign_baab
!        real(dp) :: Entry_aaaa, Entry_bbbb
!        integer :: i, j, a, b
!        integer :: Indij_ab, Indij_aa, Indij_bb, Indij_ba
!        integer :: Indab_ab, Indab_aa, Indab_bb, Indab_ba
!        real(dp) :: AllTwoElRDM(1000,1000)
!        real(dp) :: TwoElRDM(1000,1000)
!
!        Trace_2RDM_Inst = 0.D0
!        Trace_2RDM = 0.D0
!
!        do i = 1, nBasis - 1, 2
!
!            do a = 1, nBasis - 1, 2
!
!                do j = i, nBasis - 1, 2
!
!                    do b = a, nBasis - 1, 2
!
!                        ! ialpha ialpha -> jalpha jalpha
!                        ! = ibeta ibeta -> jbeta jbeta
!                        if((i.ne.j).and.(a.ne.b)) then
!                            Indij_bb = ( ( (j-2) * (j-1) ) / 2 ) + i
!                            Indab_bb = ( ( (b-2) * (b-1) ) / 2 ) + a
!                            Entry_bbbb = AllTwoElRDM(Indij_bb,Indab_bb)
!
!                            Indij_aa = ( ( ((j+1)-2) * ((j+1)-1) ) / 2 ) + (i+1)
!                            Indab_aa = ( ( ((b+1)-2) * ((b+1)-1) ) / 2 ) + (a+1)
!                            Entry_aaaa = AllTwoElRDM(Indij_aa,Indab_aa)
!
!                            AllTwoElRDM(Indij_bb,Indab_bb) = &
!                                        ( Entry_bbbb + Entry_aaaa ) / 2.D0 
!
!!                            AllTwoElRDM(Indij_aa,Indab_aa) = &
!                                        ( Entry_bbbb + Entry_aaaa ) / 2.D0 
!
!                            if(Indij_aa.eq.Indab_aa) then
!                                Trace_2RDM_Inst = Trace_2RDM_Inst + TwoElRDM(Indij_aa,Indab_aa)  
!                                Trace_2RDM = Trace_2RDM + AllTwoElRDM(Indij_aa,Indab_aa)  
!                            endif
!                            if(Indij_bb.eq.Indab_bb) then
!                                Trace_2RDM_Inst = Trace_2RDM_Inst + TwoElRDM(Indij_bb,Indab_bb)  
!                                Trace_2RDM = Trace_2RDM + AllTwoElRDM(Indij_bb,Indab_bb)  
!                            endif
!
!                        endif
!
!                        ! alpha beta -> alpha beta
!                        ! = beta alpha -> beta alpha
!                        Indij_ba = ( ( (max(i,j+1)-2) * (max(i,j+1)-1) ) / 2 ) + min(i,j+1)
!                        Indab_ba = ( ( (max(a,b+1)-2) * (max(a,b+1)-1) ) / 2 ) + min(a,b+1)
!                        Entry_baba = AllTwoElRDM(Indij_ba, Indab_ba)
!
!                        ! assume i < j+1 and a < b+1 - actually this must be true.
!                        ! keep this in for now anyway.
!                        Sign_baba = 1.D0
!                        if( ((max(i,j+1).eq.i).and.(max(a,b+1).ne.a)) .or. & 
!                            ((max(i,j+1).eq.(j+1)).and.(max(a,b+1).ne.(b+1))) ) Sign_baba = -1.D0
!
!                        Indij_ab = ( ( (max(i+1,j)-2) * (max(i+1,j)-1) ) / 2 ) + min(i+1,j)
!                        Indab_ab = ( ( (max(a+1,b)-2) * (max(a+1,b)-1) ) / 2 ) + min(a+1,b)
!                        Entry_abab = AllTwoElRDM(Indij_ab,Indab_ab)
!
!                        ! assume i+1 < j and a+1 < b - not always true if i and j are from 
!                        ! the same spatial orbital.
!                        Sign_abab = 1.D0
!                        if( ((max(i+1,j).ne.j).and.(max(a+1,b).eq.b)) .or. &
!                            ((max(i+1,j).ne.(i+1)).and.(max(a+1,b).eq.(a+1))) ) Sign_abab = -1.D0
!
!                        
!                        AllTwoElRDM(Indij_ba, Indab_ba) = &
!                            Sign_baba * ( ( ( Entry_baba * Sign_baba ) + &
!                                            ( Entry_abab * Sign_abab) ) / 2.D0 )
!
!
!                        AllTwoElRDM(Indij_ab,Indab_ab) = &
!                            Sign_abab * ( ( ( Entry_baba * Sign_baba ) + &
!                                                ( Entry_abab * Sign_abab ) )  / 2.D0 )
!
!                        if(Indij_ba.eq.Indab_ba) then
!                            Trace_2RDM_Inst = Trace_2RDM_Inst + TwoElRDM(Indij_ba,Indab_ba)  
!                            Trace_2RDM = Trace_2RDM + AllTwoElRDM(Indij_ba,Indab_ba)  
!                        endif
!!                        if(Indij_ab.eq.Indab_ab) then
!                            if(i.ne.j) then
!                                Trace_2RDM_Inst = Trace_2RDM_Inst + TwoElRDM(Indij_ab,Indab_ab)  
!                                Trace_2RDM = Trace_2RDM + AllTwoElRDM(Indij_ab,Indab_ab)  
!                            endif
!                        endif
!
!                        ! beta alpha -> alpha beta
!                        ! = alpha beta -> beta alpha
!                        Indij_ba = ( ( (max(i,j+1)-2) * (max(i,j+1)-1) ) / 2 ) + min(i,j+1)
!                        Indab_ab = ( ( (max(a+1,b)-2) * (max(a+1,b)-1) ) / 2 ) + min(a+1,b)
!                        Entry_baab = AllTwoElRDM(Indij_ba,Indab_ab)
!
!                        Sign_baab = 1.D0
!                        if( ((max(i,j+1).eq.i).and.(max(a+1,b).eq.b)) .or. &
!                            ((max(i,j+1).eq.(j+1)).and.(max(a+1,b).eq.(a+1))) ) Sign_baab = -1.D0
!
!                        Indij_ab = ( ( (max(i+1,j)-2) * (max(i+1,j)-1) ) / 2 ) + min(i+1,j)
!                        Indab_ba = ( ( (max(a,b+1)-2) * (max(a,b+1)-1) ) / 2 ) + min(a,b+1)
!                        Entry_abba = AllTwoElRDM(Indij_ab,Indab_ba)
!
!                        Sign_abba = 1.D0
!                        if( ((max(i+1,j).eq.(i+1)).and.(max(a,b+1).eq.(b+1))) .or. &
!                            ((max(i+1,j).eq.j).and.(max(a,b+1).eq.a)) ) Sign_abba = -1.D0
!
!                        
!                        AllTwoElRDM(Indij_ba,Indab_ab) = &
!                            Sign_baab * ( ( (Entry_baab * Sign_baab) + &
!                                                (Entry_abba * Sign_abba) ) / 2.D0 )
!
!                        AllTwoElRDM(Indij_ab,Indab_ba) = & 
!                            Sign_abba * ( ( (Entry_baab * Sign_baab) + &
!                                                (Entry_abba * Sign_abba) ) / 2.D0 )
!
!                    enddo
!                enddo
!            enddo
!        enddo
!
!    
!    END SUBROUTINE Average_Spins_and_Sum_2e_Norms
!
    SUBROUTINE Test_Energy_Calc()
! This routine calculates the energy based on the simple expression Energy = Sum_IJ C_I C_J H_IJ 
! where I and J are determinants    
! (rather than occupied orbitals), and H_IJ is the Hamiltonian element between them.
        USE FciMCData , only : TotWalkers,CurrentDets,iluthf
        USE Determinants, only : get_helement
        USE bit_reps , only : extract_bit_rep, extract_sign,nifdbo
        USE DetBitOps, only : detbiteq
        USE UMatCache, only: GTID
        implicit none
        INTEGER :: I, J, nI(NEl), nJ(NEl), FlagsI, FlagsJ, IC, Ex(2,2)
        INTEGER :: k,l,k2,l2,a2,b2,i2,j2, AllCurrentDetsTag
        INTEGER(int64) :: AllTotWalkers_local
        LOGICAL :: tParity
        INTEGER, DIMENSION(lenof_sign) :: SignI, SignJ
        REAL(dp) :: Test_Energy, Sum_Coeffs, SignIreal, SignJreal 
        REAL(dp) :: Test_Energy_1El, Test_Energy_2El, ParityFactor
        REAL(dp) :: Coul, Exch 
        REAL(dp) , ALLOCATABLE :: TestRDM(:,:)
        INTEGER(n_int) , ALLOCATABLE :: AllCurrentDets(:,:)
        HElement_t :: H_IJ
        INTEGER :: Ind1,Ind2,TestRDMTag,ierr,comm
        INTEGER :: lengthsout(0:nProcessors-1), disp(0:nProcessors-1)
        CHARACTER(len=*), PARAMETER :: this_routine='Test_Energy_Calc'

        WRITE(6,*) '****************'
        WRITE(6,*) '**** TESTING ENERGY CALCULATION **** '

        AllTotWalkers_local = 0
        CALL MPIReduce(TotWalkers,MPI_SUM,AllTotWalkers_local)

        IF(iProcIndex.eq.0) THEN
!            ALLOCATE(TestRDM(((nBasis*(nBasis-1))/2),((nBasis*(nBasis-1))/2)),stat=ierr)
!            IF(ierr.ne.0) CALL Stop_All('test_energy_calc','Problem allocating TestRDM array,')
!            TestRDM(:,:)=0.D0

            ALLOCATE(AllCurrentDets(0:NIfTot,AllTotWalkers_local),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating AllCurrentDets array,')
            CALL LogMemAlloc('AllCurrentDets',int(AllTotWalkers_local*(NIfTot+1),sizeof_int),size_n_int,&
                                'Test_Energy_Calc',AllCurrentDetsTag,ierr)
            AllCurrentDets(0:NIfTot,1:AllTotWalkers_local)=0
        ENDIF

        lengthsout(0:nProcessors-1) = 0
        CALL MPIAllGather(int(TotWalkers,sizeof_int)*(NIfTot+1),lengthsout,ierr)

        disp(:) = 0
        do i = 1, nProcessors-1
            disp(i) = disp(i-1) + lengthsout(i-1)
        enddo
        CALL MPIGatherv(CurrentDets(0:NIfTot,1:TotWalkers), AllCurrentDets, &
                                                lengthsout, Disp, ierr)

        IF(iProcIndex.eq.0) THEN

            Sum_Coeffs = 0.D0
            do i = 1, AllTotWalkers_local
                call extract_sign (AllCurrentDets(:,I), SignI)
                Sum_Coeffs = Sum_Coeffs + ( REAL(SignI(1)) * REAL(SignI(1)) )
            enddo
!            WRITE(6,*) 'Sum_Coeffs',Sum_Coeffs
            Sum_Coeffs = SQRT( Sum_Coeffs )
            
! Just need to run over all occupied determinants, find the matrix element between them, 
! and sum them all up.
            Test_Energy = 0.D0
            Test_Energy_1El = 0.D0
            Test_Energy_2El = 0.D0

!            WRITE(6,*) 'SignIreal'

            do I = 1, AllTotWalkers_local
                call extract_bit_rep (AllCurrentDets(:,I), nI, SignI, FlagsI)
                SignIreal = REAL(SignI(1)) * ( 1.D0 / Sum_Coeffs )

!                WRITE(6,*) SignIreal

                do J = 1, AllTotWalkers_local
                    call extract_bit_rep (AllCurrentDets(:,J), nJ, SignJ, FlagsJ)

                    SignJreal = REAL(SignJ(1)) * ( 1.D0 / Sum_Coeffs )

                    H_IJ = get_helement(nI, nJ, AllCurrentDets(:,I), AllCurrentDets(:,J), IC)

!                    WRITE(6,*) 'H_IJ',H_IJ
!                    CALL FLUSH(6)

                    IF(IC.eq.0) THEN
                        !Add to both 1El and 2El contributions.

!                        WRITE(6,*) 'nI',nI

                        Test_Energy_1El = Test_Energy_1El + &
                                ( SignIreal * SignIreal * REAL( TMAT2D(nI(NEl),nI(NEl)) ) )

!                        TestRDM(SymLabelListInv(nI(NEl)),SymLabelListInv(nI(NEl))) &
!                                = TestRDM(SymLabelListInv(nI(NEl)),SymLabelListInv(nI(NEl))) &
!                                                                + (SignIreal * SignIreal)

                        do k = 1, NEl - 1
                            k2 = gtID(nI(k))

                            Test_Energy_1El = Test_Energy_1El + &
                                        ( SignIreal * SignIreal * REAL( TMAT2D(nI(k),nI(k)) ) )

!                            TestRDM(SymLabelListInv(nI(k)),SymLabelListInv(nI(k))) &
!                                    = TestRDM(SymLabelListInv(nI(k)),SymLabelListInv(nI(k))) &
!                                                                       + (SignIreal * SignIreal)

                            ! For i,j -> i,j type doubles, UMat(i,j,i,j) &
!                                                = UMat(i,j,j,i) (so <ij||ij> = 0.
                            do l = k + 1, NEl 

                                l2 = gtID(nI(l))

                                IF( G1(nI(k))%Ms .eq. G1(nI(l))%Ms ) THEN
                                    Exch = REAL(UMAT(UMatInd(k2,l2,l2,k2,0,0)),8)
                                ELSE
                                    Exch = 0.D0
                                ENDIF

                                Test_Energy_2El = Test_Energy_2El + ( SignIreal * SignIreal &
                                            * (REAL(UMAT(UMatInd(k2,l2,k2,l2,0,0)),8) - Exch) )

!                                Ind1 = ( ( (nI(l)-2) * (nI(l)-1) ) / 2 ) + nI(k)
!                                TestRDM(Ind1,Ind1) = (REAL(UMAT(UMatInd(k2,l2,k2,l2,0,0)),8) - Exch)  
!                                TestRDM(Ind1,Ind1) = TestRDM(Ind1,Ind1) + (SignIreal * SignIreal)

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

                        Test_Energy_1El = Test_Energy_1El + &
                            ( ParityFactor * SignIreal * SignJreal * REAL( TMAT2D(Ex(1,1),Ex(2,1)) ) )

!                        TestRDM(SymLabelListInv(Ex(1,1)),SymLabelListInv(Ex(2,1))) &
!                            = TestRDM(SymLabelListInv(Ex(1,1)),SymLabelListInv(Ex(2,1))) + &
!                                                    (ParityFactor * SignIreal * SignJreal)
!                        TestRDM(Ex(1,1),Ex(2,1)) = TestRDM(Ex(1,1),Ex(2,1)) + &
!                                                        (ParityFactor * SignIreal * SignJreal)
!                        TestRDM(Ind1,Ind2) = TestRDM(Ind1,Ind2) + ( ParityFactor * SignIreal * SignJreal) 

                        IF(G1(Ex(1,1))%Ms.eq.G1(Ex(2,1))%Ms) THEN

                            do k = 1, NEl
                                IF(nI(k).eq.Ex(1,1)) CYCLE

                                k2 = gtID(nI(k))
                                l2 = gtID(Ex(1,1))
                                a2 = gtID(Ex(2,1))

                                Coul =  REAL(UMAT(UMatInd(k2,l2,k2,a2,0,0)),8)

                                IF(G1(nI(k))%Ms.eq.G1(Ex(1,1))%Ms) THEN
                                    Exch = REAL(UMAT(UMatInd(k2,l2,a2,k2,0,0)),8)
                                ELSE
                                    Exch = 0.D0
                                ENDIF

                                Test_Energy_2El = Test_Energy_2El + ( ParityFactor * SignIreal * &
                                                                    SignJreal * (Coul - Exch) )

!                                Ind1 = ( ( (MAX(nI(k),Ex(1,1))-2) * (MAX(nI(k),Ex(1,1))-1) ) / 2 ) &
!                                                                            + MIN(nI(k),Ex(1,1))
!                                Ind2 = ( ( (MAX(nI(k),Ex(2,1))-2) * (MAX(nI(k),Ex(2,1))-1) ) / 2 ) &
!                                                                            + MIN(nI(k),Ex(2,1))
!                                TestRDM(Ind1,Ind2) = ( REAL(UMAT(UMatInd(k2,l2,k2,a2,0,0)),8) - Exch)
!                                TestRDM(Ind1,Ind2) = TestRDM(Ind1,Ind2) + &
!                                                        ( ParityFactor * SignIreal * SignJreal) 

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

                        Test_Energy_2El = Test_Energy_2El + &
                                    ( ParityFactor * SignIreal * SignJreal * (Coul - Exch) )

!                        Ind1 = ( ( (Ex(1,2)-2) * (Ex(1,2)-1) ) / 2 ) + Ex(1,1)
!                        Ind2 = ( ( (Ex(2,2)-2) * (Ex(2,2)-1) ) / 2 ) + Ex(2,1)
!                        TestRDM(Ind1,Ind2) = ( (Coul - Exch))
!                        TestRDM(Ind1,Ind2) = TestRDM(Ind1,Ind2) + &
!                                                (ParityFactor * SignIreal * SignJreal)

                    ENDIF

                    Test_Energy = Test_Energy + ( SignIreal * SignJreal * REAL(H_IJ) )
                enddo

            enddo

            WRITE(6,*) 'Test energy contribution from the 1-el-RDM :',Test_Energy_1El
            WRITE(6,*) 'Test energy contribution from the 2-el-RDM :',Test_Energy_2El
            WRITE(6,*) 'The sum of these plus the core energy = ',&
                                                    Test_Energy_1El + Test_Energy_2El + Ecore
            WRITE(6,*) '       ********        '
            WRITE(6,*) 'The TEST ENERGY calculated from Sum_IJ C_I C_J H_IJ  = ',Test_Energy
            WRITE(6,*) '       ********        '
!            WRITE(6,*) 'The 2 electron part should be : ',Test_Energy - Test_Energy_1El

!            WRITE(6,*) 'SymLabelListInv',SymLabelListInv

!            do i = 1, nBasis
!                i2 = gtID(i)
!
!                WRITE(6,'(A8,4I5,2F30.15)') '** diag',i,i,i2,i2,TestRDM(SymLabelListInv(i),&
!                                            SymLabelListInv(i)),NatOrbMat(SymLabelListInv(i),&
!                                            SymLabelListInv(i)) * ( REAL(NEl,8) / Trace_1RDM ) 
!
!                do k = i + 1, nBasis
!                    k2 = gtID(k)
!                    Ind1 = ( ( (k-2) * (k-1) ) / 2 ) + i
!
!                    IF( G1(i)%Ms .eq. G1(k)%Ms ) THEN
!                        Exch = REAL(UMAT(UMatInd(i2,k2,k2,i2,0,0)),8)
!                    ELSE
!                        Exch = 0.D0
!                    ENDIF
!
!                    WRITE(6,'(6I5,2F30.15)') i,k,i2,k2,SymLabelListInv(i),SymLabelListInv(k), &
!                                            TestRDM(SymLabelListInv(i),SymLabelListInv(k)),&
!                                            NatOrbMat(SymLabelListInv(i),SymLabelListInv(k)) * &
!                                            ( REAL(NEl,8) / Trace_1RDM ) 
!    
!                    WRITE(6,*) '** diag',i,k,TestRDM(Ind1,Ind1),UMATTemp(Ind1,Ind1)
!!!                    WRITE(6,'(A8,6I5,2F30.15)') '** diag',i,k,i,k,Ind1,Ind1,TestRDM(Ind1,Ind1),&
!!!                    (AllTwoElRDM(Ind1,Ind1)*( (0.50 * (REAL(NEl) * (REAL(NEl) - 1.D0))) &
!!!                     / Trace_2RDM )) !,&
!                       REAL(UMAT(UMatInd(i2,k2,i2,k2,0,0)),8),Exch
!!                    do j = 1, nBasis
!!                        j2 = gtID(j)
!!                        do l = j + 1, nBasis
!!                            l2 = gtID(l)
!!                            Ind2 = ( ( (l-2) * (l-1) ) / 2 ) + j
!    
!                             IF( (G1(i)%Ms .eq. G1(j)%Ms ).and.(G1(k)%Ms .eq. G1(l)%Ms )) THEN
!                                 Coul = REAL(UMAT(UMatInd(i2,k2,j2,l2,0,0)),8)
!                             ELSE
!                                 Coul = 0.D0
!                             ENDIF
!    
!                             IF( (G1(i)%Ms .eq. G1(l)%Ms ).and.(G1(k)%Ms .eq. G1(j)%Ms )) THEN
!                                 Exch = REAL(UMAT(UMatInd(i2,k2,l2,j2,0,0)),8)
!                             ELSE
!                                 Exch = 0.D0
!                             ENDIF
!    
!!                            IF(Ind1.ne.Ind2) THEN
!!                                IF(AllTwoElRDM(Ind1,Ind2).ne.0) THEN
!                                    IF((TestRDM(Ind1,Ind2).ne.0).or.(UMATTemp(Ind1,Ind2).ne.0)) &
!                                            WRITE(6,*) i,k,j,l,Ind1,Ind2,&
!                                            TestRDM(Ind1,Ind2),UMATTemp(Ind1,Ind2)
!                                    IF(ABS(TestRDM(Ind1,Ind2)-UMATTemp(Ind1,Ind2)).gt.1E-15) &
!                                            WRITE(6,*) 'diff',TestRDM(Ind1,Ind2)-UMATTemp(Ind1,Ind2)
!!                                    IF((TestRDM(Ind1,Ind2).ne.0).or.(AllTwoElRDM(Ind1,Ind2).ne.0)) & 
!!                                           WRITE(6,'(10I5,2F30.15)') i,k,j,l,i2,k2,j2,l2,Ind1,Ind2,&
!!                                                                       TestRDM(Ind1,Ind2),&
!!                                                                    (AllTwoElRDM(Ind1,Ind2) * &
!!                                      ( (0.50 * (REAL(NEl) * (REAL(NEl) - 1.D0))) / Trace_2RDM )) !, &
!                                                                                Coul,Exch
!!                                    IF(ABS(TestRDM(Ind1,Ind2)-(AllTwoElRDM(Ind1,Ind2)* &
!!                                            ( (0.50 * (REAL(NEl) * (REAL(NEl) - 1.D0))) &
!!                                              / Trace_2RDM )) ).gt.1E-15) &
!!                                                   WRITE(6,*) 'diff',TestRDM(Ind1,Ind2)- &
!!                                                   (AllTwoElRDM(Ind1,Ind2) *( (0.50 * (REAL(NEl) * &
!!                                                   (REAL(NEl) - 1.D0))) / Trace_2RDM )) 
!     
!!                                ENDIF
!!                            ENDIF
!!                        enddo
!!                    enddo
!                enddo
!            enddo

!            DEALLOCATE(TestRDM)

            DEALLOCATE(AllCurrentDets)
            CALL LogMemDeAlloc('test_energy_calc',AllCurrentDetsTag)


        ENDIF


    END SUBROUTINE Test_Energy_Calc


END MODULE nElRDMMod


    SUBROUTINE DiDj_Found_FillRDM(Spawned_No,iLutJ,realSignJ)
! This routine is called when we have found a Di (or multiple Di's) spawning onto a Dj 
! with sign /= 0 (i.e. occupied).
! We then want to run through all the Di, Dj pairs and add their coefficients 
! (with appropriate de-biasing factors) into the 1 and 2 electron RDM.
        USE FciMCData , only : Spawned_Parents, Spawned_Parents_Index, iLutHF
        USE bit_reps , only : NIfTot, NIfDBO, decode_bit_det
        USE nElRDMMod , only : Add_RDM_From_IJ_Pair, Fill_Spin_Coupled_RDM, Fill_Spin_Coupled_RDM_v2
        USE Logging , only : tHF_S_D_Ref, tHF_S_D
        USE SystemData , only : NEl,tHPHF
        USE constants , only : n_int, dp, lenof_sign
        USE DetBitOps , only : DetBitEQ, FindBitExcitLevel
        IMPLICIT NONE
        integer , intent(in) :: Spawned_No
        integer(kind=n_int) , intent(in) :: iLutJ(0:NIfTot)
        real(dp) , intent(in) :: realSignJ
        integer :: i, j, nI(NEl), nJ(NEl), walkExcitLevel
        real(dp) :: realSignI
        logical :: tParity, tDetAdded

! Spawning from multiple parents, to iLutJ, which has SignJ.        

! We are at position Spawned_No in the SpawnedParts array.
! Spawned_Parents_Index(1,Spawned_No) is therefore the start position of the list of parents (Di's) 
! which spawned on the Dj in SpawnedParts(Spawned_No)
! There are Spawned_Parents_Index(2,Spawned_No) of these parent Di's.
! Spawned_Parents(0:NIfDBO,x) is the determinant Di, Spawned_Parents(NIfDBO+1,x) is the un-biased ci.

        ! Run through all Di's.
        do i = Spawned_Parents_Index(1,Spawned_No), &
                Spawned_Parents_Index(1,Spawned_No) + Spawned_Parents_Index(2,Spawned_No) - 1 

            IF(tHF_S_D_Ref.or.tHF_S_D) THEN
                ! In the case of the HF_S_D_Ref option, we'll only be in 
                ! this loop if the Dj is le 4.
                ! And for HF_S_D if Dj has excitation level le 2.
                ! Calc excitation level of Di - this needs to be 1 or 2 in 
                ! both cases (connections to the HF have already been included).
                walkExcitLevel = FindBitExcitLevel (iLutHF, Spawned_Parents(0:NIfDBO,i), 2)
                IF(walkExcitLevel.gt.2) CYCLE
                IF(walkExcitLevel.eq.0) CYCLE
            ELSEIF(DetBitEQ(iLutHF,Spawned_Parents(0:NIfDBO,i),NIfDBO)) then
                ! We've already added HF - S, and HF - D symmetrically.
                ! Any connection with the HF has therefore already been added.
                CYCLE
            ENDIF
            
            call decode_bit_det (nI, Spawned_Parents(0:NIfDBO,i))
            call decode_bit_det (nJ, iLutJ)

!            write(6,*) 'nI',nI
!            write(6,*) 'nJ',nJ
!            write(6,*) 'nI walkexcitlevel',walkexcitlevel

            ! Ci and Cj.
            realSignI = transfer( Spawned_Parents(NIfDBO+1,i), realSignI )

            !SignJ passed in as real (realSignJ)

            ! Given the Di,Dj and Ci,Cj - find the orbitals involved in the excitation, 
            ! and therefore the RDM elements we want to add the Ci.Cj to.
            IF(tHPHF) THEN
                call Fill_Spin_Coupled_RDM_v2(Spawned_Parents(0:NIfDBO,i),iLutJ,nI,nJ,&
                                                    realSignI,realSignJ,.false.)
            ELSE
                call Add_RDM_From_IJ_Pair(nI,nJ,realSignI,realSignJ,.false.)
            ENDIF

        enddo

    END SUBROUTINE DiDj_Found_FillRDM

    subroutine Fill_Diag_RDM_FromOrbs(j, nI_Prev, Prev, AvSignDi, SignFac)
! Given an occupied orbital, and the sign of the determinant it is occupied in,
! add in the contribution to the 1 and 2-electron RDMs.
        USE SystemData , only : NEl, tStoreSpinOrbs
        USE constants , only : dp, lenof_sign
        USE nElRDMMod , only : aaaa_RDM, abab_RDM
        USE NatOrbsMod , only : NatOrbMat
        USE Logging , only : RDMExcitLevel
        USE RotateOrbsData , only : SymLabelListInv
        USE UMatCache , only : GTID
        implicit none
        integer , intent(in) :: j, nI_Prev(NEl), Prev
        real(dp) , intent(in) :: AvSignDi, SignFac
        integer :: i, Ind, iSpat, jSpat, jInd

! Need to add in the diagonal elements.
! The RDM are always in spin orbitals, so just adding the orbital as is, is fine.

!        WRITE(6,*) realSignDi

        IF(RDMExcitLevel.eq.1) THEN
            if(tStoreSpinOrbs) then
                jInd = SymLabelListInv(j)
            else
                jInd = SymLabelListInv(gtID(j))
            endif
            NatOrbMat(jInd,jInd) = NatOrbMat(jInd,jInd) &
                                      + ( AvSignDi * AvSignDi * SignFac )
        ELSE

! There is no need to use the SymLabelList arrays for the 2 el RDM because we are 
! not diagonalising or anything.

            jSpat = gtID(j)
            ! effectively running from i = 1 -> j-1.
            ! so jSpat >= iSpat (can only be equal if different spin).
            do i=1,Prev
                iSpat = gtID(nI_Prev(i))

                ! either alpha alpha or beta beta -> aaaa array.
                if( ((mod(j,2).ne.0).and.(mod(nI_Prev(i),2).ne.0)) .or. &
                    ((mod(j,2).eq.0).and.(mod(nI_Prev(i),2).eq.0)) ) then

                    ! Ind doesn't include diagonal terms (when iSpat = jSpat).
                    Ind=( ( (jSpat-2) * (jSpat-1) ) / 2 ) + iSpat
                    aaaa_RDM( Ind , Ind ) = aaaa_RDM( Ind , Ind ) &
                                            + ( AvSignDi * AvSignDi * SignFac )

                ! either alpha beta or beta alpha -> abab array.                                              
                else

                    ! Ind does include diagonal terms (when iSpat = jSpat)
                    Ind=( ( (jSpat-1) * jSpat ) / 2 ) + iSpat
                    abab_RDM( Ind , Ind ) = abab_RDM( Ind , Ind ) &
                                            + ( AvSignDi * AvSignDi * SignFac )

                endif
            enddo

        ENDIF

!        write(6,*) 'adding to diagonal',AvSignDi * AvSignDi * SignFac

    end subroutine Fill_Diag_RDM_FromOrbs
 

