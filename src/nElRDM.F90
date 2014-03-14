! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
MODULe nElRDMMod
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
! Can initially find all elements and average the two values pq and qp (more accurate?).
! But should put a condition into the excitaiton generator so that only single excitations with q > p 
! are generated.

! By finding the full 1RDM, we have the ability to derive the natural orbitals as well as electron 
! densities etc.
        
    use Global_Utilities
    use Parallel_neci
    use bit_reps, only: NIfTot, NIfDBO, decode_bit_det, extract_bit_rep, &
                        encode_sign, extract_sign
    use bit_rep_data, only: flag_deterministic, test_flag
    use IntegralsData, only: UMAT
    use UMatCache, only: UMatInd, GTID
    use SystemData, only: NEl, nBasis, tStoreSpinOrbs, G1, BRR, lNoSymmetry, &
                          ARR, tUseMP2VarDenMat, Ecore, LMS, tHPHF, tFixLz, &
                          iMaxLz, tRef_Not_HF, tOddS_hphf
    use NatOrbsMod, only: NatOrbMat, NatOrbMatTag, Evalues, EvaluesTag, &
                          SetupNatOrbLabels
    use CalcData, only: MemoryFacPart, tRegenDiagHEls, NMCyc, InitiatorWalkNo
    use OneEInts, only: TMAT2D
    use FciMCData, only: MaxWalkersPart, MaxSpawned, Spawned_Parents, &
                         PreviousCycles, Spawned_Parents_Index, &
                         Spawned_ParentsTag, AccumRDMNorm_Inst, &
                         Spawned_Parents_IndexTag, Iter, AccumRDMNorm, &
                         AvNoatHF, tSinglePartPhase, AllAccumRDMNorm, iLutRef,&
                         HFDet_True, NCurrH, ilutHF_True, SpawnVec, &
                         SpawnVec2, SpawnVecTag, SpawnVec2Tag, SpawnedParts, &
                         SpawnedParts2, excit_gen_store_type, CurrentDets, &
                         CurrentH, IterRDMStart, ValidSpawnedList, &
                         TempSpawnedPartsInd, TempSpawnedParts, TotParts, &
                         TotWalkers, iLutHF, core_space, IterLastRDMFill, &
                         determ_proc_sizes,determ_proc_indices, partial_determ_vector, &
                         full_determ_vector, full_determ_vector_av, tHashWalkerList
    use LoggingData, only: RDMExcitLevel, tROFciDump, NoDumpTruncs, tHF_S_D, &
                       tExplicitAllRDM, tHF_S_D_Ref, tHF_Ref_Explicit, &
                       tHF_S_D, tPrint1RDM, tInitiatorRDM, RDMEnergyIter, &
                       tDo_Not_Calc_RDMEnergy, tDiagRDM, tReadRDMs, &
                       tPopsfile, tNo_RDMs_to_read, twrite_RDMs_to_read, &
                       tWriteMultRDMs, tDumpForcesInfo, IterRDMonFly, &
                       tWrite_normalised_RDMs, IterWriteRDMs, tPrintRODump, &
                       tNoNOTransform, tTruncRODump, tRDMonfly, tInitiatorRDMDiag, &
                       tTaperDiagRDM, tTaperSQDiagRDM, tCorrectRDMErf, erf_factor1, &
                       erf_factor2, ThreshOccRDM, tThreshOccRDMDiag,tDipoles, &
                       tBrokenSymNOs,occ_numb_diff
    use RotateOrbsData, only: CoeffT1Tag, tTurnStoreSpinOff, NoFrozenVirt, &
                              SymLabelCounts2_rot,SymLabelList2_rot, &
                              SymLabelListInv_rot,NoOrbs, SpatOrbs, &
                              SymLabelCounts2_rotTag, SymLabelList2_rotTag, &
                              NoRotOrbs, SymLabelListInv_rotTag, NoOrbs
    use DeterminantData, only: write_det
    use hphf_integrals, only: hphf_sign
    use HPHFRandExcitMod, only: FindExcitBitDetSym, FindDetSpinSym
    use DetBitOps, only: TestClosedShellDet, FindBitExcitLevel, DetBitEQ, &
                         EncodeBitDet, DetBitLT, get_bit_excitmat
    use DetCalcData, only: Det, FCIDets
    use hist_data, only: AllHistogram, Histogram
    use RotateOrbsMod, only: FourIndInts, FourIndIntsTag, PrintROFCIDUMP, &
                             PrintRepeatROFCIDUMP
    use hist, only: find_hist_coeff_explicit
    use sparse_arrays, only: sparse_core_ham, SparseCoreHamilTags, deallocate_sparse_ham, &
                            core_connections
    use SymExcit3, only: GenExcitations3
    use OneEInts, only: TMAT2D
    use hash, only: DetermineDetNode
    use constants
    use util_mod
    use sort_mod

        IMPLICIT NONE
        INTEGER , ALLOCATABLE :: Sing_InitExcSlots(:),Sing_ExcList(:)
        INTEGER , ALLOCATABLE :: Doub_InitExcSlots(:),Doub_ExcList(:)
        INTEGER(kind=n_int) , ALLOCATABLE :: Sing_ExcDjs(:,:),Sing_ExcDjs2(:,:)
        INTEGER(kind=n_int) , ALLOCATABLE :: Doub_ExcDjs(:,:),Doub_ExcDjs2(:,:)
        INTEGER :: Sing_ExcDjsTag,Sing_ExcDjs2Tag,aaaa_RDMTag,All_aaaa_RDMTag
        INTEGER :: Doub_ExcDjsTag,Doub_ExcDjs2Tag,UMATTempTag
        INTEGER :: Energies_unit, ActualStochSign_unit, abab_RDMTag, All_abab_RDMTag
        INTEGER :: abba_RDMTag, All_abba_RDMTag, NoSymLabelCounts, Rho_iiTag
        REAL(dp) , ALLOCATABLE :: aaaa_RDM(:,:), abab_RDM(:,:), abba_RDM(:,:)
        REAL(dp) , ALLOCATABLE :: All_aaaa_RDM(:,:),All_abab_RDM(:,:), All_abba_RDM(:,:)
        REAL(dp) , ALLOCATABLE :: UMATTemp(:,:), Rho_ii(:)
        REAL(dp) , ALLOCATABLE :: Lagrangian(:,:)
        REAL(dp) :: OneEl_Gap,TwoEl_Gap, Normalisation,Trace_2RDM_Inst, Trace_2RDM, Trace_1RDM, norm
        LOGICAL :: tFinalRDMEnergy, tCalc_RDMEnergy
        type(timer), save :: nElRDM_Time, FinaliseRDM_time, RDMEnergy_time
        logical :: trotatedNOs=.false.

    contains

    SUBROUTINE InitRDM()
! This routine initialises any of the arrays needed to calculate the reduced density matrix.    
! It is used for both the explicit and stochastic RDMs.
        INTEGER :: ierr,i, MemoryAlloc, MemoryAlloc_Root
        CHARACTER(len=*), PARAMETER :: this_routine='InitRDM'

! First thing is to check we're not trying to fill the RDMs in a way that is 
! not compatible with the code (not every case has been accounted for yet).
#ifdef __CMPLX
        CAll Stop_All(this_routine,'Filling of reduced density matrices not working with &
                                    &complex walkers yet.')
#endif
        if(tRDMonFly.and.tHashWalkerList) &
            call stop_all("FciMCPar", "Linear scaling + RDMs doesn't give the correct &
            & solution yet. If realcoeffs are turned on, it runs but gives the wrong RDM energy, &
            & but in a very non-obvious way.  No particular elements or groups of elements &
            & seem especially wrong, but the RDM energy is clearly incorrect.  If using integer &
            & coefficients it doesn't even run -- this is a problem somewhere in annihilation &
            & that comes about from the sign and flag being stored together")
        
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
            write(6,'(A)',advance='no') " incl. explicit connections to the following HF determinant:"
            call write_det (6, HFDet_True, .true.)
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
            NatOrbMat(:,:)=0.0_dp

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
            aaaa_RDM(:,:)=0.0_dp

            ! The 2-RDM of the type alpha beta beta alpha ( = beta alpha alpha beta).
            ! These also *do not* also include 2-RDM(i,j,a,b) terms where i=j or a=b (these are the same as the abab elements).
            ALLOCATE(abba_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating abba_RDM array,')
            CALL LogMemAlloc('abba_RDM',(((SpatOrbs*(SpatOrbs-1))/2)**2),8,this_routine,abba_RDMTag,ierr)
            abba_RDM(:,:)=0.0_dp

            MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 

            ! The 2-RDM of the type alpha beta alpha beta ( = beta alpha beta alpha).
            ! These *do* include 2-RDM(i,j,a,b) terms where i=j or a=b, if they're different spin this 
            ! is possible - hence the slightly different size to the aaaa array.
            ALLOCATE(abab_RDM(((SpatOrbs*(SpatOrbs+1))/2),((SpatOrbs*(SpatOrbs+1))/2)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating abab_RDM array,')
            CALL LogMemAlloc('abab_RDM',(((SpatOrbs*(SpatOrbs+1))/2)**2),8,this_routine,abab_RDMTag,ierr)
            abab_RDM(:,:)=0.0_dp

            MemoryAlloc = MemoryAlloc + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 ) * 8 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 ) * 8 ) 

            IF(iProcIndex.eq.0) THEN
                ! Each of these currently need to be stored on the root as well as each node.
                ! This allows us to separately calculate the instantaneous energy.
                ! TODO : Option to only calculate the accumulated RDMs - cut storage on the root in half.
                ALLOCATE(All_aaaa_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating All_aaaa_RDM array,')
                CALL LogMemAlloc('All_aaaa_RDM',(((SpatOrbs*(SpatOrbs-1))/2)**2),8,this_routine,All_aaaa_RDMTag,ierr)
                All_aaaa_RDM(:,:)=0.0_dp

                ALLOCATE(All_abab_RDM(((SpatOrbs*(SpatOrbs+1))/2),((SpatOrbs*(SpatOrbs+1))/2)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating All_abab_RDM array,')
                CALL LogMemAlloc('All_abab_RDM',(((SpatOrbs*(SpatOrbs+1))/2)**2),8,this_routine,All_abab_RDMTag,ierr)
                All_abab_RDM(:,:)=0.0_dp

                ALLOCATE(All_abba_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating All_abba_RDM array,')
                CALL LogMemAlloc('All_abba_RDM',(((SpatOrbs*(SpatOrbs-1))/2)**2),8,this_routine,All_abba_RDMTag,ierr)
                All_abba_RDM(:,:)=0.0_dp

                MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs-1))/2 ) ** 2 ) * 2 * 8 ) 
                MemoryAlloc_Root = MemoryAlloc_Root + ( ( ( (SpatOrbs*(SpatOrbs+1))/2 ) ** 2 ) * 8 ) 

                if(tDiagRDM.or.tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) then
                    ! Still need to allocate 1-RDM to get nat orb occupation numbers.
                    ALLOCATE(NatOrbMat(NoOrbs,NoOrbs),stat=ierr)
                    IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating NatOrbMat array,')
                    CALL LogMemAlloc('NatOrbMat',NoOrbs**2,8,this_routine,NatOrbMatTag,ierr)
                    NatOrbMat(:,:)=0.0_dp
                    MemoryAlloc_Root = MemoryAlloc_Root + ( NoOrbs * NoOrbs * 8 ) 
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

            OneEl_Gap=(REAL(NEl,dp)*REAL(nBasis,dp)*MemoryFacPart)/REAL(nProcessors,dp)

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
                TwoEl_Gap=(((REAL(NEl,dp)*REAL(nBasis,dp))**2)*MemoryFacPart)/REAL(nProcessors,dp)

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
                    & REAL(MemoryAlloc_Root,dp)/1048576.0_dp," Mb/Processor on the root, and ", &
                    & REAL(MemoryAlloc,dp)/1048576.0_dp," Mb/Processor on other processors."
        endif

        ! These parameters are set for the set up of the symmetry arrays, which are later used 
        ! for the diagonalisation / rotation of the 1-RDMs.

        if(tStoreSpinOrbs) then
            if(tFixLz) then
                NoSymLabelCounts = 16 * ( (2 * iMaxLz) + 1 )
            else
                NoSymLabelCounts = 16 
            endif
        else
            if(tFixLz) then
                NoSymLabelCounts = 8 * ( (2 * iMaxLz) + 1 )
            else
                NoSymLabelCounts = 8
            endif
        endif

        IF((RDMExcitLevel.eq.1).or.tDiagRDM.or.tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) THEN
            ! These arrays contain indexing systems to order the 1-RDM orbitals in terms of 
            ! symmetry.
            ! This allows the diagonalisation of the RDMs to be done in symmetry blocks (a lot 
            ! quicker/easier).
            ! The 2-RDM does not need to be reordered as it's never diagonalised. 

            ALLOCATE(SymLabelCounts2_rot(2,NoSymLabelCounts),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating SymLabelCounts2_rot array,')
            CALL LogMemAlloc('SymLabelCounts2_rot',2*NoSymLabelCounts,4,this_routine,SymLabelCounts2_rotTag,ierr)
            SymLabelCounts2_rot(:,:)=0

            ALLOCATE(SymLabelList2_rot(NoOrbs),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating SymLabelList2_rot array,')
            CALL LogMemAlloc('SymLabelList2_rot',NoOrbs,4,this_routine,SymLabelList2_rotTag,ierr)
            SymLabelList2_rot(:)=0                     
     
            ALLOCATE(SymLabelListInv_rot(NoOrbs),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating SymLabelListInv_rot array,')
            CALL LogMemAlloc('SymLabelListInv_rot',NoOrbs,4,this_routine,SymLabelListInv_rotTag,ierr)
            SymLabelListInv_rot(:)=0   

            IF((iProcIndex.eq.0).and.tDiagRDM) THEN
                ALLOCATE(Evalues(NoOrbs),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Evalues array,')
                CALL LogMemAlloc('Evalues',NoOrbs,8,this_routine,EvaluesTag,ierr)
                Evalues(:)=0.0_dp

                ALLOCATE(Rho_ii(NoOrbs),stat=ierr)
                IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating Rho_ii array,')
                CALL LogMemAlloc('Rho_ii',NoOrbs,8,this_routine,Rho_iiTag,ierr)
                Rho_ii(:)=0.0_dp

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
        logical :: exists_aaaa,exists_abab,exists_abba,exists_one
        integer :: RDM_unit, FileEnd
        integer :: i,j,a,b,Ind1,Ind2
        real(dp) :: Temp_RDM_Element, Norm_2RDM

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

                        NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = Temp_RDM_Element
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
                    call neci_flush(6)
                    CALL Stop_All('Read_in_RDMs',"Attempting to read in the RDMs, &
                                    &but at least one of the TwoRDM_a***_TOREAD files are missing.")
                endif

            endif
        endif

        ! Calculate the energy for the matrices read in (if we're calculating more than the 1-RDM).
        if(tCalc_RDMEnergy) then
            call Calc_Energy_from_RDM(Norm_2RDM)
        endif

! Continue calculating the RDMs from the first iteration when the POPSFILES (and RDMs) are read in.
! This overwrites the iteration number put in the input.
        IterRDMonFly = 0

    end subroutine Read_In_RDMs

    SUBROUTINE SetUpSymLabels_RDM() 
! This routine just sets up the symmetry labels so that
! the orbitals are ordered according to symmetry (all beta then all alpha if spin orbs).
        INTEGER , ALLOCATABLE :: SymOrbs_rot(:)
        INTEGER :: LabOrbsTag, SymOrbs_rotTag, ierr, i , j, SpatSym, LzSym 
        INTEGER :: lo, hi, Symi, SymCurr, Symi2, SymCurr2
        CHARACTER(len=*) , PARAMETER :: this_routine = 'SetUpSymLabels_RDM'

        ! This is only allocated temporarily to be used to order the orbitals by.
        ALLOCATE(SymOrbs_rot(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymOrbs_rot',NoOrbs,4,this_routine,SymOrbs_rotTag,ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,"Mem allocation for SymOrbs_rot failed.")

! Brr has the orbital numbers in order of energy... i.e Brr(2) = the orbital index 
! with the second lowest energy.
! Brr is always in spin orbitals.
!        write(6,*) 'BRR'
!        do i=1,nBasis
!            WRITE(6,*) i, BRR(i), INT(G1(BRR(i))%sym%S)
!        enddo
!        CALL neci_flush(6)
!        CALL Stop_All('','')

! Now we want to put the spatial orbital index, followed by the symmetry.        
        SymLabelList2_rot(:) = 0
        SymOrbs_rot(:)=0

! *** STEP 1 *** Fill SymLabelList2_rot.
! find the orbitals and order them in terms of symmetry.
        do i=1,SpatOrbs
            if(tStoreSpinOrbs) then
                ! for open shell systems, all alpha are followed by all beta.
                SymLabelList2_rot(i) = BRR(2*i)
                SymLabelList2_rot(i+SpatOrbs) = BRR((2*i)-1)

                if(tFixLz) then
                    SpatSym = INT(G1(BRR(2*i))%sym%S)
                    LzSym = INT(G1(BRR(2*i))%Ml)
                    SymOrbs_rot(i) = ( SpatSym * ((2 * iMaxLz) + 1) ) + ( LzSym + iMaxLz )

                    SpatSym = INT(G1(BRR((2*i)-1))%sym%S)
                    LzSym = INT(G1(BRR((2*i)-1))%Ml)
                    SymOrbs_rot(i+SpatOrbs) = ( SpatSym * ((2 * iMaxLz) + 1) ) + ( LzSym + iMaxLz )
                else
                    SymOrbs_rot(i) = INT(G1(BRR(2*i))%sym%S) 
                    SymOrbs_rot(i+SpatOrbs) = INT(G1(BRR((2*i)-1))%sym%S) 
                endif
            else
                SymLabelList2_rot(i) = gtID(BRR(2*i))
                if(tFixLz) then
                    SpatSym = INT(G1(BRR(2*i))%sym%S)
                    LzSym = INT(G1(BRR(2*i))%Ml)
                    SymOrbs_rot(i) = ( SpatSym * ((2 * iMaxLz) + 1) ) + ( LzSym + iMaxLz )
                else
                    SymOrbs_rot(i) = INT(G1(BRR(2*i))%sym%S) 
                endif
                ! Orbital BRR(2*i) for i = 1 will be the beta orbital with the 
                ! second lowest energy - want the spatial orbital index to go with this.
                ! G1 also in spin orbitals - get symmetry of this beta orbital, will 
                ! be the same as the spatial orbital.
            endif
        enddo

        call sort (SymOrbs_rot(1:SpatOrbs), SymLabelList2_rot(1:SpatOrbs))
        ! Sorts SymLabelList2_rot according to the order of SymOrbs_rot (i.e. in terms of symmetry). 
        if(tStoreSpinOrbs) &
            call sort (SymOrbs_rot(SpatOrbs+1:nBasis), SymLabelList2_rot(SpatOrbs+1:nBasis))
            ! Also do this for the beta set if spin orbitals.

!*** STEP 2 *** Fill SymLabelCounts2_rot_rot. This is like SymLabelCounts2_rot, but has a different ordering - BEWARE
!SymLabelCounts(1,:) contains the position in SymLabelList2_rot where the symmetry index starts,
!SymLabelCounts(2,:) contains the number of orbitals in that symmetry index.
!Again if spin orbs, all alpha are followed by all beta - i.e. first 8 refer to alpha, second 8 to beta.

        IF(lNoSymmetry) THEN
            ! if we are ignoring symmetry, all orbitals essentially have symmetry 0.
            SymLabelCounts2_rot(1,1) = 1
            SymLabelCounts2_rot(2,1) = SpatOrbs
            if(tStoreSpinOrbs) then
                SymLabelCounts2_rot(1,9) = SpatOrbs+1
                SymLabelCounts2_rot(2,9) = SpatOrbs
            endif
        ELSE 
            ! otherwise we run through the occupied orbitals, counting the number with 
            ! each symmetry (spatial and Lz) and noting where in SymLabelList2_rot each symmetry block starts.
            SymCurr = 0
            SymLabelCounts2_rot(1,1) = 1
            if(tStoreSpinOrbs) then
                SymCurr2 = 0
                SymLabelCounts2_rot(1,9) = SpatOrbs + 1
            endif
            do i=1,SpatOrbs
                if(tStoreSpinOrbs) then
                    Symi=SymOrbs_rot(i)
                    Symi2=SymOrbs_rot(i + SpatOrbs)
                else
                    Symi=SymOrbs_rot(i)
                endif
                SymLabelCounts2_rot(2,(Symi+1))= SymLabelCounts2_rot(2,(Symi+1))+1
                IF(Symi.ne.SymCurr) THEN
                    do j = SymCurr + 1, Symi
                        SymLabelCounts2_rot(1,(j+1))=i
                    enddo
                    SymCurr=Symi
                ENDIF
                if(tStoreSpinOrbs) then
                    SymLabelCounts2_rot(2,(Symi2+9))= SymLabelCounts2_rot(2,(Symi2+9))+1
                    IF(Symi2.ne.SymCurr2) THEN
                        do j = SymCurr2 + 1, Symi2
                            SymLabelCounts2_rot(1,(j+9))=i + SpatOrbs
                        enddo
                        SymCurr2=Symi2
                    ENDIF
                endif
            enddo
        ENDIF

        ! Go through each symmetry group, making sure the orbitals are 
        ! ordered lowest to highest within each symmetry.
        do i=1,NoSymLabelCounts
            IF(SymLabelCounts2_rot(2,i).ne.0) THEN
                lo = SymLabelCounts2_rot(1, i)
                hi = lo + SymLabelCounts2_rot(2, i) - 1
                call sort (SymLabelList2_rot (lo:hi))
            ENDIF
        enddo

        ! Construct the inverse matrix.  While SymLabelList2_rot takes a position and tells us 
        ! what orbital is in it, we also might need to take an orbital and find out what 
        ! position to put its contribution in.
        do i=1,NoOrbs
            SymLabelListInv_rot(SymLabelList2_rot(i))=i
        enddo

! Deallocate the arrays just used in this routine.
        DEALLOCATE(SymOrbs_rot)
        CALL LogMemDealloc(this_routine,SymOrbs_rotTag)

!        WRITE(6,*) 'Sym Label Counts'
!        do i=1,8
!            WRITE(6,*) i,SymLabelCounts2_rot(1,i),SymLabelCounts2_rot(2,i)
!        enddo
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), &
!                               & and their symmetries according to G1'
!        do i=1,SpatOrbs
!            WRITE(6,*) i,SymLabelList2_rot(i),INT(G1(2*SymLabelList2_rot(i))%sym%S)
!        enddo
!        WRITE(6,*) 'i','ARR(SymLabelList2_rot(i),1)','ARR(SymLabelList2_rot(i),2)','Sym'
!        do i=1,SpatOrbs
!             WRITE(6,*) i,ARR(2*SymLabelList2_rot(i),1),ARR(2*SymLabelList2_rot(i),2),&
!                                               INT(G1(2*SymLabelList2_rot(i))%sym%S)
!        enddo
!
!        WRITE(6,*) 'Sym label list (i.e the orbitals in symm order), and its inverse'
!        do i=1,SpatOrbs
!            WRITE(6,*) SymLabelList2_rot(i),SymLabelListInv_rot(i)
!        enddo
!        CALL neci_flush(6)
!        CALL Stop_All('SetUpSymLabels_RDM','Checking orbital labelling.')

    END SUBROUTINE SetUpSymLabels_RDM

    subroutine DeAlloc_Alloc_SpawnedParts()
! When calculating the RDMs, we need to store the parent from which a child is spawned along with the 
! children in the spawned array.
! This means a slightly larger array is communicated between processors, which there is no point in doing 
! for the first part of the calculation.
! When we start calculating the RDMs this routine is called and the SpawnedParts array is made larger to 
! accommodate the parents.
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
                                        REAL(((NIfTot+1)*MaxSpawned*2*size_n_int),dp)/1048576.0_dp,' to ',&
                                        REAL(((NIfTot+NIfDBO+3)*MaxSpawned*2*size_n_int),dp)/1048576.0_dp, ' Mb/Processor'

    end subroutine DeAlloc_Alloc_SpawnedParts

    subroutine extract_bit_rep_avsign_no_rdm(iLutnI, CurrH_I, nI, SignI, &
                                                FlagsI, IterRDMStartI, AvSignI, Store)
! This is just the standard extract_bit_rep routine for when we're not calculating the RDMs.    
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp) , intent(in) :: CurrH_I(NCurrH)
        integer, intent(out) :: nI(nel), FlagsI
        real(dp), dimension(lenof_sign), intent(out) :: SignI
        real(dp) , intent(out) :: IterRDMStartI, AvSignI
        type(excit_gen_store_type), intent(inout), optional :: Store

        call extract_bit_rep (iLutnI, nI, SignI, FlagsI, Store)

        IterRDMStartI = 0.0_dp
        AvSignI = 0.0_dp

    end subroutine extract_bit_rep_avsign_no_rdm

    subroutine extract_bit_rep_avsign_norm(iLutnI, CurrH_I, nI, SignI, &
                                                FlagsI, IterRDMStartI, AvSignI, Store)
! The following extract_bit_rep_avsign routine extracts the bit representation 
! of the current determinant, and calculate the average sign since this determinant became 
! occupied. 
! This is called for each determinant in the occupied list at the beginning of its FCIQMC cycle. 
! It is used if we're calculating the RDMs with or without HPHF. 
! Input:    iLutnI (bit rep of current determinant).
!           CurrH_I(1) - diagonal H_ii element for current det i.
!           CurrH_I(2) - previous AvSign, CurrH_I(3) - previous IterRDMStartI
! Output:   nI, SignI, FlagsI after extract.                                              
!           IterRDMStartI - new iteration the determinant became occupied (as a real).
!           AvSignI - the new average walker population during this time (also real).
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp) , intent(in) :: CurrH_I(NCurrH)
        integer, intent(out) :: nI(nel), FlagsI
        real(dp), dimension(lenof_sign), intent(out) :: SignI
        real(dp) , intent(out) :: IterRDMStartI, AvSignI
        type(excit_gen_store_type), intent(inout), optional :: Store

        ! This is the iteration from which this determinant has been occupied.
        IterRDMStartI = CurrH_I(3)
        ! If there is nothing stored there yet, the first iteration the determinant 
        ! became occupied is this one.
        IF(IterRDMStartI.eq.0.0_dp) IterRDMStartI = real(Iter, dp)

        ! This extracts everything.
        call extract_bit_rep (iLutnI, nI, SignI, FlagsI)

        ! Update the average population.
        ! This just comes out as the current population (SignI) if this is the first 
        ! time the determinant has become occupied.
        AvSignI = ( ((real(Iter,dp) - IterRDMStartI) * CurrH_I(2)) &
                        + SignI(1) ) / ( real(Iter,dp) - IterRDMStartI + 1.0_dp )

    end subroutine extract_bit_rep_avsign_norm

    subroutine fill_rdm_diag_currdet_norm(iLutnI, nI, CurrH_I, ExcitLevelI, tCoreSpaceDet, IterLastRDMFill) 
! This routine calculates the diagonal RDM contribution, and explicit connections to the HF, from the 
! current determinant. 
! Each determinant (iLutnI/nI), has Hii element CurrH_I(1), average 
! population CurrH_I(2), and has been occupied since iteration CurrH_I(3).
! IterLastRDMFill is the number of iterations since the last time the RDM contributions were added in 
! (often the frequency of the RDM energy calculation). 
! We need to multiply the RDM contributions by either this, or the number of iterations 
! the determinant has been occupied, which ever is fewer.
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp) , intent(in) :: CurrH_I(NCurrH)
        logical, intent(in), optional :: tCoreSpaceDet
        integer, intent(in) :: nI(nel), ExcitLevelI, IterLastRDMFill
        real(dp) :: IterDetOcc, IterRDM
        integer(n_int) :: SpinCoupDet(0:nIfTot)
        integer :: nSpinCoup(nel), SignFac, HPHFExcitLevel

        ! This is the number of iterations this determinant has been occupied.
        IterDetOcc = real(Iter,dp) - CurrH_I(3) + 1.0_dp

        ! IterRDM is then the number of iterations we want to multiply the contributions by.
        IterRDM = min(IterDetOcc,real(IterLastRDMFill,dp))
        
        if(tHPHF) then
            if(.not.TestClosedShellDet(iLutnI)) then
                call Fill_Diag_RDM(nI, CurrH_I(2)/SQRT(2.0_dp), IterRDM,tCoreSpaceDet)
                
! C_X D_X = C_X / SQRT(2) [ D_I +/- D_I'] - for open shell dets, divide stored C_X by SQRT(2). 
! Add in I.
                call FindExcitBitDetSym(iLutnI, SpinCoupDet)
                call decode_bit_det (nSpinCoup, SpinCoupDet)
                ! Find out if it's + or - in the above expression.                
                SignFac = hphf_sign(iLutnI)

                call Fill_Diag_RDM(nSpinCoup, (real(SignFac,dp)*CurrH_I(2))/SQRT(2.0_dp), IterRDM,tCoreSpaceDet)

! For HPHF we're considering < D_I + D_I' | a_a+ a_b+ a_j a_i | D_I + D_I' >
! Not only do we have diagonal < D_I | a_a+ a_b+ a_j a_i | D_I > terms, but also cross terms
! < D_I | a_a+ a_b+ a_j a_i | D_I' > if D_I and D_I' can be connected by a single or double 
! excitation.
! Find excitation level between D_I and D_I' and add in the contribution if connected.
                HPHFExcitLevel = FindBitExcitLevel (iLutnI, SpinCoupDet, 2)
                if(HPHFExcitLevel.le.2) & 
                    call Add_RDM_From_IJ_Pair(nI, nSpinCoup, &
                                                IterRDM * CurrH_I(2)/SQRT(2.0_dp), &
                                                (real(SignFac,dp)*CurrH_I(2))/SQRT(2.0_dp), .true.)

            else

                ! HPHF on, but determinant closed shell.
                call Fill_Diag_RDM(nI, CurrH_I(2), IterRDM,tCoreSpaceDet)

            endif

            call Add_RDM_HFConnections_HPHF(iLutnI, nI, CurrH_I(2), ExcitLevelI, IterRDM)   

        else
            ! No HPHF

            call Fill_Diag_RDM(nI, CurrH_I(2), IterRDM, tCoreSpaceDet)
            
            call Add_RDM_HFConnections_Norm(iLutnI, nI, CurrH_I(2), ExcitLevelI, IterRDM)   

        endif

    end subroutine fill_rdm_diag_currdet_norm

    subroutine fill_rdm_diag_currdet_hfsd(iLutnI, nI, CurrH_I, ExcitLevelI, tCoreSpaceDet, IterLastRDMFill) 
! This routine calculates the diagonal RDM contribution, and explicit connections to the HF, from the 
! current determinant. 
! Each determinant (iLutnI/nI), has Hii element CurrH_I(1), average 
! population CurrH_I(2), and has been occupied since iteration CurrH_I(3).
! IterLastRDMFill is the number of iterations since the last time the RDM contributions were added in 
! (often the frequency of the RDM energy calculation). 
! We need to multiply the RDM contributions by either this, or the number of iterations 
! the determinant has been occupied, which ever is fewer.
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp) , intent(in) :: CurrH_I(NCurrH)
        logical, intent(in), optional :: tCoreSpaceDet
        integer, intent(in) :: nI(nel), ExcitLevelI, IterLastRDMFill
        real(dp) :: IterDetOcc, IterRDM
        integer(n_int) :: SpinCoupDet(0:nIfTot)
        integer :: nSpinCoup(nel), SignFac, HPHFExcitLevel

        ! This is the number of iterations this determinant has been occupied.
        IterDetOcc = real(Iter,dp) - CurrH_I(3) + 1.0_dp

        ! IterRDM is then the number of iterations we want to multiply the contributions by.
        IterRDM = min(IterDetOcc,real(IterLastRDMFill,dp))

        call Add_RDM_HFConnections_HF_S_D(iLutnI, nI, CurrH_I(2), ExcitLevelI, IterRDM)

    end subroutine fill_rdm_diag_currdet_hfsd

    subroutine fill_rdm_softexit(nDets)
! !This routine is called if a softexit command is used.  
! !In this case, the diagonal RDM elements, and explicit connections to the HF will not have been 
! !calculated in the final iteration.  
! !Need to run over the occupied determinants and do so now. 
! !This is clearly not all that efficient, and so could be incorporated into the popsfile write out 
! !if it becomes an issue.
        integer(int64) , intent(in) :: nDets
  !      integer :: nI(nel), ExcitLevel, i

  !      ! If it happens to be an iteration where the diagonal RDM elements are already being 
  !      ! calculated, then no need to do this.
  !      if(.not.((Iter.eq.NMCyc).or.(mod((Iter - IterRDMStart + 1),RDMEnergyIter).eq.0))) then

  !          ! IterLastRDMFill is the number of iterations from the last time the RDM elements 
  !          ! were included.
  !          IterLastRDMFill = mod((Iter - IterRDMStart + 1),RDMEnergyIter)

  !          ! Run over all determinants in the occupied (CurrentDets) list.
  !          do i = 1, int(nDets,sizeof_int)

  !              ! If we are only using initiators to calculate the RDMs, only add in the diagonal and 
  !              ! explicit contributions if the average population is greater than n_add = InitiatorWalkNo.
  !              if((abs(CurrentH(2,i)).gt.real(InitiatorWalkNo,dp)).or.(.not.tInitiatorRDM)) then

  !                  if(tRef_Not_HF) then
  !                      ExcitLevel = FindBitExcitLevel (iLutHF_true, CurrentDets(:,i), 2) 
  !                  else
  !                      ExcitLevel = FindBitExcitLevel (iLutRef, CurrentDets(:,i), 2)
  !                  endif

  !                  call decode_bit_det (nI, CurrentDets(:,i))


  !                  if(tHF_Ref_Explicit.or.tHF_S_D.or.tHF_S_D_Ref) then
  !                      call fill_rdm_diag_currdet_hfsd(CurrentDets(:,i), nI, CurrentH(:,i), &
  !                                                                      ExcitLevel)
  !                  else
  !                      call fill_rdm_diag_currdet_norm(CurrentDets(:,i), nI, CurrentH(:,i), &
  !                                                                      ExcitLevel)
  !                  endif
  !                  
  !              endif

  !          enddo

  !          !Now would need to add in off-diagonal contributions from the deterministic space

  !      endif

    end subroutine fill_rdm_softexit

    subroutine det_removed_fill_diag_rdm( iLutnI, CurrH_I )
! This routine is called if a determinant is removed from the list of currently occupied.  
! At this point we need to add in its diagonal contribution for the number of iterations it has 
! been occupied (or since the contribution was last included). 
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        real(dp) , intent(in) :: CurrH_I(NCurrH)
        integer :: nI(nel), ExcitLevel

        ! If the determinant is removed on an iteration that the diagonal RDM elements are 
        ! already being calculated, it will already have been counted.
        if(.not.((Iter.eq.NMCyc).or.(mod((Iter - IterRDMStart + 1),RDMEnergyIter).eq.0))) then

            if((abs(CurrH_I(2)).gt.InitiatorWalkNo).or.(.not.tInitiatorRDM)) then

                ! IterLastRDMFill is the number of iterations from the last time the RDM elements 
                ! were included.
                IterLastRDMFill = mod((Iter - IterRDMStart + 1),RDMEnergyIter)

                call decode_bit_det (nI, iLutnI)
                if(tRef_Not_HF) then
                    ExcitLevel = FindBitExcitLevel (iLutHF_true, iLutnI, 2)
                else
                    ExcitLevel = FindBitExcitLevel (iLutRef, iLutnI, 2)
                endif

                if(tHF_Ref_Explicit.or.tHF_S_D.or.tHF_S_D_Ref) then
                    call fill_rdm_diag_currdet_hfsd(iLutnI, nI, CurrH_I, ExcitLevel, .false., IterLastRDMFill)
                else
                    call fill_rdm_diag_currdet_norm(iLutnI, nI, CurrH_I, ExcitLevel, .false., IterLastRDMFill)
                endif

            endif

        endif

    end subroutine det_removed_fill_diag_rdm

! These Add_RDM_HFConnections routines take the current determinant and if it is 
! a single or double excitation of the HF, they explicitly add in the contribution to the RDM 
! from the current and hf determinant.

    subroutine Add_RDM_HFConnections_Norm(iLutJ,nJ,AvSignJ,walkExcitLevel,IterRDM)
! This is called when we run over all TotWalkers in CurrentDets.    
! It is called for each CurrentDet which is a single or double of the HF.
! It explicitly adds in the HF - S/D connection, as if the HF were D_i and 
! the single or double D_j.
! This is the standard full space RDM calc (No HPHF).
! In this case the diagonal elements wll already be taken care of.
        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        integer , intent(in) :: nJ(NEl)
        real(dp) , intent(in) :: AvSignJ, IterRDM
        integer , intent(in) :: walkExcitLevel
        integer(kind=n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(NEl), HPHFExcitLevel

! The two measures are calculated and updated at different places in the code, so it's not necessarily
! an issue that they don't agree.  E.g. whether or not you rezero the stats if the sign goes to 
! zero instantaneously during annihilation or not.
! Quick check that the HF population is being calculated correctly.
        if(walkExcitLevel.eq.0) then
            if(AvSignJ.ne.AvNoatHF) then
                write(6,*) 'HFDet_True',HFDet_True
                write(6,*) 'nJ',nJ
                write(6,*) 'iLutJ',iLutJ
                write(6,*) 'AvSignJ',AvSignJ
                write(6,*) 'AvNoatHF',AvNoatHF
                CALL Stop_All('Add_RDM_HFConnections_Norm','Incorrect instantaneous HF population.')
            endif
        endif

! If we have a single or double, add in the connection to the HF, symmetrically.       
        if((walkExcitLevel.eq.1).or.(walkExcitLevel.eq.2)) &
            call Add_RDM_From_IJ_Pair(HFDet_True,nJ,AvNoatHF,IterRDM * AvSignJ,.true.)

    end subroutine Add_RDM_HFConnections_Norm


    subroutine Add_RDM_HFConnections_HPHF(iLutJ,nJ,AvSignJ,walkExcitLevel,IterRDM)
! This is called when we run over all TotWalkers in CurrentDets.    
! It is called for each CurrentDet which is a single or double of the HF.
! It adds in the HF - S/D connection.
! The diagonal elements will already have been taken care of by the extract routine.
        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        integer , intent(in) :: nJ(NEl)
        real(dp) , intent(in) :: AvSignJ, IterRDM
        integer , intent(in) :: walkExcitLevel
        integer(kind=n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(NEl), HPHFExcitLevel

        if(walkExcitLevel.eq.0) then
            if(AvSignJ.ne.AvNoatHF) then
                write(6,*) 'AvSignJ',AvSignJ
                write(6,*) 'AvNoatHF',AvNoatHF
                CALL Stop_All('Add_RDM_HFConnections_HPHF','Incorrect instantaneous HF population.')
            endif
        endif

! Now if the determinant is connected to the HF (i.e. single or double), add in the diagonal elements
! of this connection as well - symmetrically because no probabilities are involved.
        if((walkExcitLevel.eq.1).or.(walkExcitLevel.eq.2)) &
            call Fill_Spin_Coupled_RDM_v2(iLutHF_True,iLutJ,HFDet_True,nJ,&
                                            AvNoatHF,IterRDM * AvSignJ,.true.)

    end subroutine Add_RDM_HFConnections_HPHF

    subroutine Add_RDM_HFConnections_HF_S_D(iLutJ,nJ,AvSignJ,walkExcitLevel, IterRDM)
! This is called when we run over all TotWalkers in CurrentDets.    
! This finds all the connections to the HF when doing some sort of truncated RDM 
! calculation.
! Here, the diagonal elements will not have been added in by the extract routines.
! In the case of HF_Ref_Explicit, this routine does all the work.
        integer(kind=n_int), intent(in) :: iLutJ(0:NIfTot)
        integer , intent(in) :: nJ(NEl)
        real(dp) , intent(in) :: AvSignJ, IterRDM
        integer , intent(in) :: walkExcitLevel
        integer(kind=n_int) :: SpinCoupDet(0:niftot)
        integer :: nSpinCoup(NEl), HPHFExcitLevel

! Add diagonal elements to reduced density matrices.

! If HF_S_D_Ref, we are only considering determinants connected to the HF, 
! doubles and singles (so theoretically up to quadruples).
! But for the diagonal elements - only consider doubles and singles (and HF).

        ! In all of these cases the HF is a diagonal element.
        if(walkExcitLevel.eq.0) then

            call Fill_Diag_RDM(nJ, AvSignJ, IterRDM)
            AccumRDMNorm = AccumRDMNorm + (AvSignJ * AvSignJ * IterRDM)
            AccumRDMNorm_Inst = AccumRDMNorm_Inst + (AvSignJ * AvSignJ * IterRDM)
    
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
                    call Fill_Spin_Coupled_RDM_v2(iLutHF_True,iLutJ,HFDet_True,nJ,&
                                AvNoatHF,AvSignJ * IterRDM,.false.)

                else

                    ! The singles and doubles are connected and explicitly calculated 
                    ! - but not symmetrically.
                    call Add_RDM_From_IJ_Pair(HFDet_True, nJ, AvNoatHF, &
                                                AvSignJ * IterRDM,.false.)

                endif

            else
                ! For the HF,S,D symmetric case, and the HF,S,D reference, the S and D
                ! are diagonal terms too.
                ! These options are not set up for HPHF.
                call Fill_Diag_RDM(nJ, AvSignJ, IterRDM)
                AccumRDMNorm = AccumRDMNorm + (AvSignJ * AvSignJ * IterRDM)
                AccumRDMNorm_Inst = AccumRDMNorm_Inst + (AvSignJ * AvSignJ * IterRDM)

                call Add_RDM_From_IJ_Pair(HFDet_True,nJ,AvNoatHF,AvSignJ * IterRDM,.true.)

            endif

        endif

    end subroutine Add_RDM_HFConnections_HF_S_D


    subroutine calc_rdmbiasfac(p_spawn_rdmfac,p_gen,AvSignCurr,SignCurr,RDMBiasFacCurr)
        real(dp), intent(in) :: p_gen, AvSignCurr
        real(dp), dimension(lenof_sign), intent(in) :: SignCurr
        real(dp) , intent(out) :: RDMBiasFacCurr
        real(dp), intent(in) :: p_spawn_rdmfac
        real(dp) :: p_notlist_rdmfac, p_spawn, p_not_spawn, p_max_walktospawn

        ! We eventually turn this real bias factor into an integer to be passed around 
        ! with the spawned children and their parents - this only works with 64 bit at the mo.
        if(n_int.eq.4) CALL Stop_All('attempt_create_normal', &
                        'the bias factor currently does not work with 32 bit integers.')

        ! Otherwise calculate the 'sign' of Di we are eventually going to add in as Di.Dj.                                
        ! Because we only add in Di.Dj when we successfully spawn from Di.Dj, we need to unbias (scale up) 
        ! Di by the probability of this happening.
        ! We need the probability that the determinant i, with population n_i, will spawn on j.
        ! We only consider one instance of a pair Di,Dj, so just want the probability of any of the n_i 
        ! walkers spawning at least once on Dj.
        ! P_successful_spawn(j | i)[n_i] =  1 - P_not_spawn(j | i)[n_i]
        ! P_not_spawn(j | i )[n_i] is the probability of none of the n_i walkers spawning on j from i.
        ! This requires either not generating j, or generating j and not succesfully spawning, n_i times.
        ! P_not_spawn(j | i )[n_i] = [(1 - P_gen(j | i)) + ( P_gen( j | i ) * (1 - P_spawn(j | i))]^n_i

        p_notlist_rdmfac = ( 1.0_dp - p_gen ) + ( p_gen * (1.0_dp - p_spawn_rdmfac) )

        ! The bias fac is now n_i / P_successful_spawn(j | i)[n_i]
       

        if(real(int(SignCurr(1)),dp).ne.SignCurr(1)) then
            !There's a non-integer population on this determinant
            !We need to consider both possibilities - whether we attempted to spawn 
            !int(SignCurr) times or int(SignCurr)+1 times
            p_max_walktospawn=abs(SignCurr(1)-real(int(SignCurr(1)),dp))
            p_not_spawn = (1.0_dp - p_max_walktospawn)*(p_notlist_rdmfac**abs(int(SignCurr(1)))) + &
                        p_max_walktospawn*(p_notlist_rdmfac**(abs(int(SignCurr(1)))+1))

        else
            p_not_spawn=p_notlist_rdmfac**(abs(SignCurr(1)))
        endif

        p_spawn=abs(1.0_dp - p_not_spawn)
        
        if(tInitiatorRDM) then
            if(abs(AvSignCurr).gt.InitiatorWalkNo) then
                ! The Di is an initiator (on average) - keep passing around its sign until we 
                ! know if the Dj is an initiator.
                RDMBiasFacCurr = AvSignCurr / p_spawn   
            else
                ! The determinant we are spawning from is not an initiator (on average) 
                ! - do not want to add this Di.Dj contribution into the RDM.
                RDMBiasFacCurr = 0.0_dp
            endif
        else
            RDMBiasFacCurr = AvSignCurr / p_spawn   
        endif
        
    end subroutine calc_rdmbiasfac

    subroutine store_parent_with_spawned(RDMBiasFacCurr, WalkerNumber, iLutI, DetSpawningAttempts, iLutJ, procJ)
    !We are spawning from iLutI to SpawnedParts(:,ValidSpawnedList(proc)).
    !This routine stores the parent (D_i) with the spawned child (D_j) so that we can
    !add in Ci.Cj to the RDM later on.
    !The parent is NIfDBO integers long, and stored in the second part of the SpawnedParts array 
    !from NIfTot+1 -> NIfTot+1 + NIfDBO.
        real(dp), intent(in) :: RDMBiasFacCurr
        integer, intent(in) :: WalkerNumber, procJ
        integer, intent(in) :: DetSpawningAttempts
        integer(kind=n_int), intent(in) :: iLutI(0:niftot),iLutJ(0:niftot)
        logical :: tRDMStoreParent
        integer :: j

        if(RDMBiasFacCurr.eq.0.0_dp) then
            ! If RDMBiasFacCurr is exactly zero, any contribution from Ci.Cj will be zero 
            ! so it is not worth carrying on. 
            ! This happens when we are only using initiators for the RDM, and Di is a non-initiator.

            SpawnedParts(niftot+1:niftot+nifdbo+2, ValidSpawnedList(procJ)) = 0
        else

            ! First we want to check if this Di.Dj pair has already been accounted for.
            ! This means searching the Dj's that have already been spawned from this Di, to make sure 
            ! the new Di being spawned on here is not the same.
            ! The Dj children spawned by the current Di are being stored in the array TempSpawnedParts, 
            ! so that the reaccurance of a Di.Dj pair may be monitored.

            ! Store the Di parent with the spawned child, unless we find this Dj has already been spawned on.
            tRDMStoreParent = .true.

            ! Run through the Dj walkers that have already been spawned from this particular Di.
            ! If this is the first to be spawned from Di, TempSpawnedPartsInd will be zero, so we 
            ! just wont run over anything.
            do j = 1,TempSpawnedPartsInd
                if(DetBitEQ(iLutJ(0:NIfDBO),TempSpawnedParts(0:NIfDBO,j),NIfDBO)) then
                    ! If this Dj is found, we do not want to store the parent with this spawned walker.
                    tRDMStoreParent = .false.
                    EXIT
                endif
            enddo

            if(tRDMStoreParent) then
                ! This is a new Dj that has been spawned from this Di.
                ! We want to store it in the temporary list of spawned parts which have come from this Di.
                if(WalkerNumber.ne.DetSpawningAttempts) then
                    ! Don't bother storing these if we're on the last walker, or if we only have one 
                    ! walker on Di.
                    TempSpawnedPartsInd = TempSpawnedPartsInd + 1
                    TempSpawnedParts(0:NIfDBO,TempSpawnedPartsInd) = iLutJ(0:NIfDBO)
                endif

                ! We also want to make sure the parent Di is stored with this Dj.
                SpawnedParts(niftot+1:niftot+nifdbo+1, ValidSpawnedList(procJ)) = iLutI(0:nifdbo) 

                ! We need to carry with the child (and the parent), the sign of the parent.
                ! In actual fact this is the sign of the parent divided by the probability of generating 
                ! that pair Di and Dj, to account for the 
                ! fact that Di and Dj are not always added to the RDM, but only when Di spawns on Dj.
                ! This RDMBiasFacCurr factor is turned into an integer to pass around to the relevant processors.
                SpawnedParts(niftot+nifdbo+2, ValidSpawnedList(procJ)) = &
                    transfer(RDMBiasFacCurr,SpawnedParts(niftot+nifdbo+2, ValidSpawnedList(procJ)))

            else
                ! This Di has already spawned on this Dj - don't store the Di parent with this child, 
                ! so that the pair is not double counted.  
                ! We are using the probability that Di spawns onto Dj *at least once*, so we don't want to 
                ! double count this pair.
                SpawnedParts(niftot+1:niftot+nifdbo+2, ValidSpawnedList(procJ)) = 0
            endif
        endif

    end subroutine store_parent_with_spawned


    subroutine check_fillRDM_DiDj(Spawned_No,iLutJ,realSignJ)
! The spawned parts contain the Dj's spawned by the Di's in CurrentDets.
! If the SpawnedPart is found in the CurrentDets list, it means that the Dj has a non-zero 
! cj - and therefore the Di.Dj pair will have a non-zero ci.cj to contribute to the RDM.
! The index i tells us where to look in the parent array, for the Di's to go with this Dj.
        integer , intent(in) :: Spawned_No
        integer(kind=n_int) , intent(in) :: iLutJ(0:NIfTot)
        real(dp) , intent(in) :: realSignJ
        integer :: ExcitLevel
    
        ! In all cases, we've already symmetrically added in 
        ! connections to the HF, so we don't want to re add any pair 
        ! containing the HF.
        if(tHF_S_D) then 
            ! In the case of the HF S D matrix (symmetric), Di and Dj can both 
            ! be the HF, singles or doubles.
            ! This is the excitation level of Dj.
            ExcitLevel = FindBitExcitLevel (iLutHF_True, iLutJ, 2)
            if((ExcitLevel.eq.2).or.(ExcitLevel.eq.1)) &
                CALL DiDj_Found_FillRDM(Spawned_No,iLutJ,realSignJ)
        elseif(tHF_S_D_Ref) then
            ! In the case of the HF and singles and doubles Ref, 
            ! Di is only ever the HF, and Dj is 
            ! anything connected - i.e. up to quadruples.
            ExcitLevel = FindBitExcitLevel (iLutHf_True, iLutJ, 4)
            if((ExcitLevel.le.4).and.(ExcitLevel.ne.0)) &
                CALL DiDj_Found_FillRDM(Spawned_No,iLutJ,realSignJ)
        elseif(.not.DetBitEQ(iLutHF_True,iLutJ,NIfDBO)) then
            if(tInitiatorRDM) then
                if(abs(realSignJ).gt.InitiatorWalkNo) &
                    CALL DiDj_Found_FillRDM(Spawned_No,iLutJ,realSignJ)
            else
                CALL DiDj_Found_FillRDM(Spawned_No,iLutJ,realSignJ)
            endif
        endif

    end subroutine check_fillRDM_DiDj
 

    SUBROUTINE DiDj_Found_FillRDM(Spawned_No,iLutJ,realSignJ)
! This routine is called when we have found a Di (or multiple Di's) spawning onto a Dj 
! with sign /= 0 (i.e. occupied).
! We then want to run through all the Di, Dj pairs and add their coefficients 
! (with appropriate de-biasing factors) into the 1 and 2 electron RDM.
        use bit_rep_data, only : flag_deterministic, test_flag
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
                walkExcitLevel = FindBitExcitLevel (iLutHF_True, Spawned_Parents(0:NIfDBO,i), 2)
                IF(walkExcitLevel.gt.2) CYCLE
                IF(walkExcitLevel.eq.0) CYCLE
            ELSEIF(DetBitEQ(iLutHF_True,Spawned_Parents(0:NIfDBO,i),NIfDBO)) then
                ! We've already added HF - S, and HF - D symmetrically.
                ! Any connection with the HF has therefore already been added.
                CYCLE
            ENDIF
            
            call decode_bit_det (nI, Spawned_Parents(0:NIfDBO,i))
            call decode_bit_det (nJ, iLutJ)

            realSignI = transfer( Spawned_Parents(NIfDBO+1,i), realSignI )
            
            !SignJ passed in as real (realSignJ)

            ! Given the Di,Dj and Ci,Cj - find the orbitals involved in the excitation, 
            ! and therefore the RDM elements we want to add the Ci.Cj to.
            IF(tHPHF) THEN
                call Fill_Spin_Coupled_RDM_v2(Spawned_Parents(0:NIfDBO,i), iLutJ, nI, nJ, &
                                                    realSignI, realSignJ, .false.)
            ELSE
                call Add_RDM_From_IJ_Pair(nI, nJ, realSignI, realSignJ, .false.)
            ENDIF

        enddo

    END SUBROUTINE DiDj_Found_FillRDM

! This routine does the same as Fill_Spin_Coupled_RDM, but hopefully more efficiently!
! It takes to HPHF functions, and calculate what needs to be summed into the RDMs
    subroutine Fill_Spin_Coupled_RDM_v2(iLutnI,iLutnJ,nI,nJ,realSignI,realSignJ,tFill_CiCj_Symm)
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
!Above Fill_Spin_Coupled_RDM_v2 is more efficient version of this routine.
!If the two HPHF determinants we're considering consist of I + I' and J + J', 
!where X' is the spin coupled (all spins flipped) version of X,
!then we have already considered the I -> J excitation.
!And if I and J are connected by a double excitation, tDoubleConnection is true and we have 
!also considered I' -> J'.
!But we need to also account for I -> J' and I' -> J.
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
                realSignFacJ = real(SignFacJ,dp) / SQRT(2.0_dp)

!                write(6,*) "OS -> OS "

                if (I_J_ExcLevel.le.2) then

                    ! I -> J.
                    call Add_RDM_From_IJ_Pair(nI,nJ,(realSignI/SQRT(2.0_dp)),&
                                                (realSignJ/SQRT(2.0_dp)),tFill_CiCj_Symm)
!                    write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),realSignI/SQRT(2.0_dp),realSignJ/SQRT(2.0_dp)
 
                    ! I' -> J'.
                    call Add_RDM_From_IJ_Pair(nI2,nJ2,(realSignFacI*realSignI),&
                                              (realSignFacJ*realSignJ),tFill_CiCj_Symm)
!                    write(6,"(A,4I4,2F12.6)") "I' -> J' : ",nI2(:),nJ2(:),(realSignFacI*realSignI),(realSignFacJ*realSignJ)

                endif

                if (SpinCoupI_J_ExcLevel.le.2) then

                    ! I' -> J.
                    call Add_RDM_From_IJ_Pair(nI2,nJ,(realSignFacI*realSignI),&
                                                (realSignJ/SQRT(2.0_dp)),tFill_CiCj_Symm)

!                    write(6,"(A,4I4,2F12.6)") "I' -> J : ",nI2(:),nJ(:),realSignFacI*realSignI,realSignJ/SQRT(2.0_dp)
                    ! I -> J'.
                    call Add_RDM_From_IJ_Pair(nI, nJ2,(realSignI/SQRT(2.0_dp)),&
                                                 (realSignFacJ*realSignJ),tFill_CiCj_Symm)

!                    write(6,"(A,4I4,2F12.6)") "I -> J' : ",nI(:),nJ2(:),realSignI/SQRT(2.0_dp),realSignFacJ*realSignJ

                endif

            else
                
                ! I is open shell, but J is not.
                ! Need I -> J and I' -> J.

!                write(6,*) "OS -> CS "
                ! I -> J.

                call Add_RDM_From_IJ_Pair(nI,nJ,(realSignI/SQRT(2.0_dp)),&
                                                        realSignJ,tFill_CiCj_Symm)
!                write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),realSignI/SQRT(2.0_dp),realSignJ
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
                                               (realSignJ/SQRT(2.0_dp)),tFill_CiCj_Symm)
!            write(6,"(A,4I4,2F12.6)") "I -> J : ",nI(:),nJ(:),realSignI,realSignJ/SQRT(2.0_dp)

            ! Find J'.
            call FindExcitBitDetSym(iLutnJ, iLutnJ2)
            SignFacJ = hphf_sign(iLutnJ)
            realSignFacJ = real(SignFacJ,dp) / SQRT(2.0_dp)

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
            call neci_flush(6)
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

! =======================================================================================    
! EXPLICIT ROUTINES    
! =======================================================================================    
    SUBROUTINE Fill_ExplicitRDM_this_Iter(TotWalkers)
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
!        Normalisation = 0.0_dp
!        NormalisationTemp = 0.0_dp
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
!            SignI2(1) = NINT(REAL(SignI(1)) * ( 1.0_dp / Sum_Coeffs ))
!            call encode_sign(CurrentDets(:,I), SignI2)
!        enddo

        CALL set_timer(nElRDM_Time,30)

        do i=1,int(MaxTotWalkers,sizeof_int)

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
        INTEGER(int64) , INTENT(IN) :: TotWalkers
        INTEGER(kind=n_int) :: iLutnI(0:NIfTot)
        INTEGER :: i,error
        REAL(dp) :: TempTotParts, NormalisationTemp, Sum_Coeffs,AllNode_norm
        LOGICAL :: blank_det
        REAL(dp), DIMENSION(lenof_sign) :: TempSign

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
!        Normalisation = 0.0_dp
!        NormalisationTemp = 0.0_dp
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
!            SignI2(1) = NINT(REAL(SignI(1)) * ( 1.0_dp / Sum_Coeffs ))
!            call encode_sign(CurrentDets(:,I), SignI2)
!        enddo

        CALL set_timer(nElRDM_Time,30)

        CALL MPISumAll(Histogram,AllHistogram)

        norm=0.0_dp
        if(iProcIndex.eq.0) then
            do i=1,Det
                norm=norm+AllHistogram(1,i)**2
            enddo
            norm=SQRT(norm)
        endif

        CALL MPISumAll(norm,allNode_norm)
        norm=allNode_norm
 
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

            TempSign(1) = real(i,dp)
            call encode_sign(iLutnI,TempSign)

            CALL Add_Hist_ExplicitRDM_Contrib(iLutnI,blank_det)

        enddo
        CALL halt_timer(nElRDM_Time)

    END SUBROUTINE Fill_Hist_ExplicitRDM_this_Iter


    SUBROUTINE Add_ExplicitRDM_Contrib(iLutnI,blank_det)
! This is the general routine for taking a particular determinant in the spawned list, 
! D_i and adding it's contribution to the reduced density matrix.
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
        INTEGER(kind=n_int) , INTENT(IN) :: iLutnI(0:NIfTot)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        REAL(dp), dimension(lenof_sign) :: SignDi, SignDi2
        integer :: ExcitMat3(2,2), nI(NEl), nJ(NEl), Proc, FlagsDi
        integer :: a, b, CountTemp, exflag
        logical :: tAllExcitFound, tParity

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
!                call neci_flush(6)
            exflag = 1
            CALL GenExcitations3(nI,iLutnI,nJ,exflag,ExcitMat3(:,:),tParity,&
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
!                call neci_flush(6)
                exflag = 2
                CALL GenExcitations3(nI,iLutnI,nJ,exflag,ExcitMat3(:,:),tParity,&
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
        INTEGER(kind=n_int) , INTENT(IN) :: iLutnI(0:NIfTot)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: HistPos
        REAL(dp), dimension(lenof_sign) :: RealHistPos
        integer :: ExcitMat3(2,2), nI(NEl), nJ(NEl), Proc, FlagsDi
        integer :: a, b, CountTemp, exflag
        logical :: tAllExcitFound, tParity
        real(dp) :: realSignDi

        call extract_bit_rep (iLutnI, nI, RealHistPos, FlagsDi)
! Unfortunately uses the decoded determinant - might want to look at this.        
        HistPos=int(RealHistPos)

        realSignDi = AllHistogram(1,HistPos(1))/norm
        
        call Fill_Diag_RDM(nI,realSignDi)

!        CountTemp = 0

        ExcitMat3(:,:)=0
! Zeros in ExcitMat3 starts off at the first single excitation.        
        tAllExcitFound=.false.
! This becomes true when all the excitations have been found.        

        do while (.not.tAllExcitFound)
!                write(6,*) 'generating singles'
!                call neci_flush(6)
            exflag = 1
            CALL GenExcitations3(nI,iLutnI,nJ,exflag,ExcitMat3(:,:),tParity,&
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
!                call neci_flush(6)
                exflag = 2
                CALL GenExcitations3(nI,iLutnI,nJ,exflag,ExcitMat3(:,:),tParity,&
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
        INTEGER :: i,j
        INTEGER(MPIArg) :: sendcounts(nProcessors),disps(nProcessors)
        INTEGER(MPIArg) :: sing_recvcounts(nProcessors)
        INTEGER :: error,MaxSendIndex,MaxIndex
        INTEGER(MPIArg) :: sing_recvdisps(nProcessors)
        INTEGER(MPIArg) :: doub_recvcounts(nProcessors),doub_recvdisps(nProcessors)

        do i=0,nProcessors-1
            sendcounts(i+1)=int(Sing_ExcList(i)-(NINT(OneEl_Gap*i)+1),MPIArg)
! Sendcounts is the number of singly excited determinants we want to send for 
! each processor (but goes from 1, not 0).            
            disps(i+1)=NINT(OneEl_Gap*i,MPIArg)
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
            sendcounts(i)=sendcounts(i)*(int(NIfTot+1,MPIArg))
            disps(i)=disps(i)*(int(NIfTot+1,MPIArg))
            sing_recvcounts(i)=sing_recvcounts(i)*(int(NIfTot+1,MPIArg))
            sing_recvdisps(i)=sing_recvdisps(i)*(int(NIfTot+1,MPIArg))
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
                sendcounts(i+1)=int(Doub_ExcList(i)-(NINT(TwoEl_Gap*i)+1),MPIArg)
! Sendcounts is the number of singly excited determinants we want to send for 
! each processor (but goes from 1, not 0).            
                disps(i+1)=NINT(TwoEl_Gap*i,MPIArg)
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
                sendcounts(i)=sendcounts(i)*(int(NIfTot+1,MPIArg))
                disps(i)=disps(i)*(int(NIfTot+1,MPIArg))
                doub_recvcounts(i)=doub_recvcounts(i)*(int(NIfTot+1,MPIArg))
                doub_recvdisps(i)=doub_recvdisps(i)*(int(NIfTot+1,MPIArg))
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
        INTEGER :: i,j
        INTEGER(MPIArg) :: sendcounts(nProcessors),disps(nProcessors)
        INTEGER(MPIArg) :: sing_recvcounts(nProcessors)
        INTEGER :: error,MaxSendIndex,MaxIndex
        INTEGER(MPIArg) :: sing_recvdisps(nProcessors)
        INTEGER(MPIArg) :: doub_recvcounts(nProcessors),doub_recvdisps(nProcessors)

        do i=0,nProcessors-1
            sendcounts(i+1)=int(Sing_ExcList(i)-(NINT(OneEl_Gap*i)+1),MPIArg)
! Sendcounts is the number of singly excited determinants we want to send for 
! each processor (but goes from 1, not 0).            
            disps(i+1)=NINT(OneEl_Gap*i,MPIArg)
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
            sendcounts(i)=sendcounts(i)*(int(NIfTot+1,MPIArg))
            disps(i)=disps(i)*(int(NIfTot+1,MPIArg))
            sing_recvcounts(i)=sing_recvcounts(i)*(int(NIfTot+1,MPIArg))
            sing_recvdisps(i)=sing_recvdisps(i)*(int(NIfTot+1,MPIArg))
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
                sendcounts(i+1)=int(Doub_ExcList(i)-(NINT(TwoEl_Gap*i)+1),MPIArg)
! Sendcounts is the number of singly excited determinants we want to send for 
! each processor (but goes from 1, not 0).            
                disps(i+1)=NINT(TwoEl_Gap*i,MPIArg)
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
                sendcounts(i)=sendcounts(i)*(int(NIfTot+1,MPIArg))
                disps(i)=disps(i)*(int(NIfTot+1,MPIArg))
                doub_recvcounts(i)=doub_recvcounts(i)*(int(NIfTot+1,MPIArg))
                doub_recvdisps(i)=doub_recvdisps(i)*(int(NIfTot+1,MPIArg))
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
        INTEGER(MPIArg), INTENT(IN) :: recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        REAL(dp), dimension(lenof_sign) :: SignDi,SignDj, SignDi2,SignDj2
        integer :: i, j, NoDets, StartDets, PartInd
        integer :: nI(NEl), nJ(NEl), Ex(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity
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
                    CALL BinSearchParts_rdm(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)
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
        INTEGER(MPIArg), INTENT(IN) :: recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        REAL(dp), dimension(lenof_sign) :: SignDi,SignDj, SignDi2, SignDj2
        integer :: i, j, NoDets, StartDets, PartInd
        integer :: nI(NEl), nJ(NEl), Ex(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity
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
                    CALL BinSearchParts_rdm(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

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
!        USE AnnihilationMod , only : BinSearchParts
        INTEGER(MPIArg), INTENT(IN) :: recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: HistPos
        REAL(dp), dimension(lenof_sign) :: RealHistPos

        integer :: i, j, NoDets, StartDets, PartInd, ExcitLevel
        integer :: nI(NEl), nJ(NEl), Ex(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity
        REAL(dp) :: realSignDi, realSignDj

! Take each Dj, and binary search the CurrentDets to see if it is occupied.

        do i=1,nProcessors
! Doing determinants from each processor separately because each has a different 
! D_i it was excited from.

            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            IF(NoDets.gt.1) THEN
                call extract_bit_rep (Sing_ExcDjs2(:,StartDets), nI, RealHistPos, FlagsDi)

                HistPos=int(RealHistPos)

                realSignDi = AllHistogram(1,HistPos(1))/norm

                do j=StartDets+1,(NoDets+StartDets-1)
! D_i is in the first spot - start from the second.                
                
                    iLutnJ(:)=Sing_ExcDjs2(:,j)
! This binary searches CurrentDets between 1 and TotWalkers for determinant iLutnJ.
! If found, tDetFound will be true, and PartInd the index in CurrentDets where the 
! determinant is.
                    CALL BinSearchParts_rdm(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

                    ExcitLevel = FindBitExcitLevel (iLutHF_true, iLutnJ, NEl)
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
        INTEGER(MPIArg), INTENT(IN) :: recvcounts(nProcessors),recvdisps(nProcessors)
        INTEGER(kind=n_int) :: iLutnJ(0:NIfTot)
        INTEGER, dimension(lenof_sign) :: HistPos
        REAL(dp), dimension(lenof_sign) :: RealHistPos
        integer :: i, j, NoDets, StartDets,PartInd, ExcitLevel
        integer :: nI(NEl), nJ(NEl), Ex(2,2), FlagsDi, FlagsDj
        logical :: tDetFound, tParity
        REAL(dp) :: realSignDi,realSignDj

! Take each Dj, and binary search the CurrentDets to see if it is occupied.
    
!        WRITE(6,*) 'Searching for generated nJs in occupied list'

        do i=1,nProcessors
! Doing determinants from each processor separately because each has a different D_i 
! it was excited from.

            NoDets=recvcounts(i)/(NIfTot+1)
            StartDets=(recvdisps(i)/(NIfTot+1))+1

            IF(NoDets.gt.1) THEN
                call extract_bit_rep (Doub_ExcDjs2(:,StartDets), nI, RealHistPos, FlagsDi)

                HistPos=int(RealHistPos)

                realSignDi = AllHistogram(1,HistPos(1))/norm

                do j=StartDets+1,(NoDets+StartDets-1)
! D_i is in the first spot - start from the second.                
                
                    iLutnJ(:)=Doub_ExcDjs2(:,j)

! This binary searches CurrentDets between 1 and TotWalkers for determinant iLutnJ.
! If found, tDetFound will be true, and PartInd the index in CurrentDets where the 
! determinant is.
                    CALL BinSearchParts_rdm(iLutnJ,1,int(TotWalkers,sizeof_int),PartInd,tDetFound)

                    ExcitLevel = FindBitExcitLevel (iLutHF_True, iLutnJ, NEl)
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


! =======================================================================================    
! THESE NEXT ROUTINES ARE GENERAL TO BOTH STOCHASTIC AND EXPLICIT    
! =======================================================================================    

    subroutine Fill_Diag_RDM(nI,realSignDi,RDMItersIn,tCoreSpaceDet)
! Fill diagonal elements of 1- and 2-RDM.
! These are < Di | a_i+ a_i | Di > and < Di | a_i+ a_j+ a_j a_i | Di >.
        integer , intent(in) :: nI(NEl)
        real(dp) , intent(in) :: realSignDi
        real(dp) , intent(in) , optional :: RDMItersIn
        logical , intent(in) , optional :: tCoreSpaceDet
        integer :: i, j, iSpat, jSpat, Ind, iInd
        real(dp) :: RDMIters, ScaleContribFac

! Need to add in the diagonal elements.
        
        if(.not.present(RDMItersIn)) then
            RDMIters=1.0_dp
        else
            RDMIters=RDMItersIn
        endif

        ScaleContribFac=1.0

        if ((.not. tCoreSpaceDet) .or. .not.present(tCoreSpaceDet)) then
            !Dets in the core space are never removed from main list, so strictly do not require corrections
            if (tThreshOccRDMDiag .and. (abs(RealSignDi) .le. ThreshOccRDM)) ScaleContribFac=0.0_dp
           
        endif
        
        if(RDMExcitLevel.eq.1) then
            do i=1,NEl
                if(tStoreSpinOrbs) then
                    iInd = SymLabelListInv_rot(nI(i))
                else
                    ! SymLabelListInv_rot will be in spat orbitals too.
                    iInd = SymLabelListInv_rot(gtID(nI(i)))
                endif
                NatOrbMat(iInd,iInd) = NatOrbMat(iInd,iInd) &
                                          + ( realSignDi * realSignDi * RDMIters)*ScaleContribFac 
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
                                          + ( realSignDi * realSignDi * RDMIters)*ScaleContribFac
                    ! either alpha beta or beta alpha -> abab array.                                              
                    else

                        ! Ind does include diagonal terms (when iSpat = jSpat)
                        Ind=( ( (jSpat-1) * jSpat ) / 2 ) + iSpat
                        abab_RDM( Ind , Ind ) = abab_RDM( Ind , Ind ) &
                                          + ( realSignDi * realSignDi * RDMIters)*ScaleContribFac
                    endif

                enddo
            enddo
        endif

    end subroutine Fill_Diag_RDM

    subroutine Fill_Sings_RDM(nI,Ex,tParity,realSignDi,realSignDj,tFill_CiCj_Symm)
! This routine adds in the contribution to the 1- and 2-RDM from determinants connected
! by a single excitation.
        integer , intent(in) :: nI(NEl), Ex(2,2)
        logical , intent(in) :: tParity
        real(dp) , intent(in) :: realSignDi, realSignDj
        logical , intent(in) :: tFill_CiCj_Symm
        integer :: k, Indik, Indak, iSpat, aSpat, kSpat, iInd, aInd
        real(dp) :: ParityFactor, ParityFactor2

!        WRITE(6,*) '* In singles'
!        call neci_flush(6)
!        WRITE(6,*) 'Ex(1,:)',Ex(1,:)
!        WRITE(6,*) 'Ex(2,:)',Ex(2,:)
!        WRITE(6,*) 'tParity',tParity
!        WRITE(6,*) 'nI',nI

        ParityFactor=1.0_dp
        IF(tParity) ParityFactor=-1.0_dp

        if(RDMExcitLevel.eq.1) then
            ! SymLabelList2_rot(i), gives the orbital in position i
            ! SymLabelListInv_rot(i), gives the position orbital i should go in.

            if(tStoreSpinOrbs) then
                iInd = Ex(1,1)
                aInd = Ex(2,1)
            else
                iInd = gtID(Ex(1,1))
                aInd = gtID(Ex(2,1))   ! These two must have the same spin.
            endif
            Indik = SymLabelListInv_rot(iInd)    ! Position of i 
            Indak = SymLabelListInv_rot(aInd)    ! Position of a.
            
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
                                                    ParityFactor2 = ParityFactor * (-1.0_dp)
                            IF((Ex(1,1).gt.nI(k)).and.(Ex(2,1).lt.nI(k))) &
                                                    ParityFactor2 = ParityFactor * (-1.0_dp)

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
        integer , intent(in) :: Ex(2,2)
        logical , intent(in) :: tParity
        real(dp) , intent(in) :: realSignDi, realSignDj
        logical , intent(in) :: tFill_CiCj_Symm
        integer :: Indij, Indab, iSpat, jSpat, aSpat, bSpat
        real(dp) :: ParityFactor

        ! Adding to elements Gamma(i,j,a,b)

        ParityFactor=1.0_dp
        IF(tParity) ParityFactor=-1.0_dp

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
        INTEGER :: error
        real(dp) :: Norm_2RDM, Norm_2RDM_Inst
        real(dp) :: Norm_1RDM, Trace_1RDM, SumN_Rho_ii
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

            call Finalise_1e_RDM(Norm_1RDM)  

        else
            ! We always want to calculate one final RDM energy, whether or not we're 
            ! calculating the energy throughout the calculation.
            ! Unless of course, only the 1-RDM is being calculated.

            ! Calculate the energy one last time - and write out everything we need.
            tFinalRDMEnergy = .true.

            !1RDM is contructed here (in calc_1RDM_energy)
            CALL Calc_Energy_from_RDM(Norm_2RDM)

            if(tPrint1RDM) then
                call Finalise_1e_RDM(Norm_1RDM)
            elseif(tDiagRDM.and.(iProcIndex.eq.0)) then
                call calc_1e_norms(Trace_1RDM, Norm_1RDM, SumN_Rho_ii)
                write(6,*) ''
                write(6,'(A55,F30.20)') ' SUM OF 1-RDM(i,i) FOR THE N LOWEST ENERGY &
                                            &HF ORBITALS: ',SumN_Rho_ii
            endif
            if (tDumpForcesInfo) then
                if (.not. tPrint1RDM) call Finalise_1e_RDM(Norm_1RDM)
                CALL Calc_Lagrangian_from_RDM(Norm_1RDM, Norm_2RDM)
                call convert_mats_Molpforces(Norm_1RDM, Norm_2RDM)
            endif

        endif
        call MPIBarrier(error)

! Call the routines from NatOrbs that diagonalise the one electron reduced 
! density matrix.
        tRotatedNOs = .false. ! Needed for BrokenSymNo routine
        IF(tDiagRDM) call find_nat_orb_occ_numbers()

! This is where we would likely call any further calculations of force etc.
        if(tDipoles) then
            if (.not. tPrint1RDM) call Finalise_1e_RDM(Norm_1RDM)
            call CalcDipoles(Norm_1RDM)
        endif

! After all the NO calculations are finished we'd like to do another rotation
! to obtain symmetry-broken natural orbitals
        if (tBrokenSymNOs) then
            call BrokenSymNO(occ_numb_diff)
        endif

        CALL halt_timer(FinaliseRDM_Time)

    
    end subroutine FinaliseRDM

    subroutine Finalise_1e_RDM(Norm_1RDM) 
! This routine takes the 1-RDM (NatOrbMat), normalises it, makes it 
! hermitian if required, and prints out the versions we're interested in.    
! This is only ever called at the very end of a calculation.
        use Logging , only : twrite_RDMs_to_read, twrite_normalised_RDMs
                             
        !implicit none
        integer :: i, ierr
        real(dp), intent(out) :: Norm_1RDM
        real(dp) :: Trace_1RDM, SumN_Rho_ii
        real(dp), allocatable :: AllNode_NatOrbMat(:,:)

        Norm_1RDM = 0.0_dp
        AllAccumRDMNorm = 0.0_dp
        IF(tHF_S_D_Ref.or.tHF_Ref_Explicit.or.tHF_S_D) &
            CALL MPIReduce(AccumRDMNorm,MPI_SUM,AllAccumRDMNorm)

        if(RDMExcitLevel.eq.1) then

            ALLOCATE(AllNode_NatOrbMat(NoOrbs,NoOrbs),stat=ierr)
            
            CALL MPISumAll(NatOrbMat,AllNode_NatOrbMat)
            NatOrbMat=AllNode_NatOrbMat
            
            DEALLOCATE(AllNode_NatOrbMat)

        endif

        if(iProcIndex.eq.0) then 

            ! Find the normalisation.
            call calc_1e_norms(Trace_1RDM, Norm_1RDM, SumN_Rho_ii)

            ! Write out the unnormalised, non-hermitian OneRDM_POPS.
            if(twrite_RDMs_to_read) call Write_out_1RDM(Norm_1RDM,.false.)

            ! Enforce the hermiticity condition.  If the RDMExcitLevel is not 1, the 
            ! 1-RDM has been constructed from the hermitian 2-RDM, so this will not 
            ! be necessary.
            ! The HF_Ref and HF_S_D_Ref cases are not hermitian by definition.
            if((RDMExcitLevel.eq.1).and.(.not.(tHF_Ref_Explicit.or.tHF_S_D_Ref))) then
                call make_1e_rdm_hermitian(Norm_1RDM)
            endif

            ! Write out the final, normalised, hermitian OneRDM.                
            if(twrite_normalised_RDMs) call Write_out_1RDM(Norm_1RDM,.true.)

            write(6,'(A55,F30.20)') ' SUM OF 1-RDM(i,i) FOR THE N LOWEST ENERGY HF ORBITALS: ',SumN_Rho_ii

        endif

    end subroutine Finalise_1e_RDM

    subroutine calc_1e_norms(Trace_1RDM, Norm_1RDM, SumN_Rho_ii)
! We want to 'normalise' the reduced density matrices.
! These are not even close to being normalised at the moment, because of the way they are 
! calculated on the fly.
! They should be calculated from a normalised wavefunction.
! But we know that the trace of the one electron reduced density matrix must be equal to 
! the number of the electrons.
! We can use this to find the factor we must divide the 1RDM through by.
        real(dp) , intent(out) :: Trace_1RDM, Norm_1RDM, SumN_Rho_ii
        integer :: i, HFDet_ID, BRR_ID

        Trace_1RDM = 0.0_dp
        Norm_1RDM = 0.0_dp

        do i = 1, NoOrbs
            Trace_1RDM = Trace_1RDM + NatOrbMat(i,i)
        enddo

        IF(tHF_S_D_Ref.or.tHF_Ref_Explicit.or.tHF_S_D) THEN
            Norm_1RDM = 1.0_dp / AllAccumRDMNorm
        ELSE
            ! Sum of diagonal elements of 1 electron RDM must equal NEl, 
            ! number of electrons.
            Norm_1RDM = ( REAL(NEl,dp) / Trace_1RDM )
        ENDIF

!        if(tFinalRDMEnergy) then
!            WRITE(6,*) 'AllAccumRDMNorm',AllAccumRDMNorm
!            WRITE(6,*) 'Norm_1RDM',Norm_1RDM
!            WRITE(6,*) 'Trace_1RDM',Trace_1RDM
!        endif

        !Need to multiply each element of the 1 electron reduced density matrices 
        !by NEl / Trace_1RDM,
        !and then add it's contribution to the energy.
        
! Want to sum the diagonal elements of the 1-RDM for the HF orbitals.
! Given the HF orbitals, SymLabelListInv_rot tells us their position in the 1-RDM.
        SumN_Rho_ii = 0.0_dp
        do i = 1, NoOrbs
            ! Rho_ii is the diagonal elements of the 1-RDM.
            ! Want this ordered according to the energy of the orbitals.
            ! Brr has the orbital numbers in order of energy... i.e Brr(2) = the orbital index 
            ! with the second lowest energy.
            ! Brr is always in spin orbitals.
            ! i gives the energy level, BRR gives the orbital, SymLabelListInv_rot gives the position of 
            ! this orbital in NatOrbMat.
            if(tDiagRDM) then
                if(tStoreSpinOrbs) then
                    Rho_ii(i) = NatOrbMat(SymLabelListInv_rot(BRR(i)),SymLabelListInv_rot(BRR(i))) * Norm_1RDM
                else
                    BRR_ID = gtID(BRR(2*i))
                    Rho_ii(i) = NatOrbMat(SymLabelListInv_rot(BRR_ID),SymLabelListInv_rot(BRR_ID)) * Norm_1RDM
                endif
            endif
    
            if(i.le.NEl) then
                if(tStoreSpinOrbs) then
                    SumN_Rho_ii = SumN_Rho_ii + &
                            ( NatOrbMat(SymLabelListInv_rot(HFDet_True(i)),SymLabelListInv_rot(HFDet_True(i))) &
                                * Norm_1RDM )
                else
                    HFDet_ID = gtID(HFDet_True(i))
                    SumN_Rho_ii = SumN_Rho_ii + &
                            ( NatOrbMat(SymLabelListInv_rot(HFDet_ID),SymLabelListInv_rot(HFDet_ID)) &
                                * Norm_1RDM ) / 2.0_dp
                endif
            endif
        enddo

    end subroutine calc_1e_norms

    subroutine make_1e_rdm_hermitian(Norm_1RDM)
! Simply average the 1-RDM(i,j) and 1-RDM(j,i) elements which should be equal in a perfect world.    
        real(dp) , intent(in) :: Norm_1RDM
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity 
        integer :: i, j
        real(dp) :: Temp

        Max_Error_Hermiticity = 0.0_dp
        Sum_Error_Hermiticity = 0.0_dp
        do i = 1, NoOrbs
            do j = i, NoOrbs
                IF((abs((NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM) - &
                        (NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i))*Norm_1RDM))).gt.Max_Error_Hermiticity) &
                    Max_Error_Hermiticity = abs((NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM) - &
                                                (NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i))*Norm_1RDM))

                Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                        abs((NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM) - &
                                            (NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i))*Norm_1RDM))

                Temp = (NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) + &
                        NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)))/2.0_dp

                NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = Temp
                NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = Temp
            enddo
        enddo

        ! Output the hermiticity errors.
        write(6,'(A29,F30.20)') ' MAX ABS ERROR IN 1RDM HERMITICITY', Max_Error_Hermiticity
        write(6,'(A29,F30.20)') ' SUM ABS ERROR IN 1RDM HERMITICITY', Sum_Error_Hermiticity

    end subroutine make_1e_rdm_hermitian

    subroutine Write_out_1RDM(Norm_1RDM,tNormalise)
! This routine writes out the OneRDM.
! If tNormalise is true, we are printing the normalised, hermitian matrix.
! Otherwise, Norm_1RDM is ignored and we print both 1-RDM(i,j) and 1-RDM(j,i) (in binary) 
! for the OneRDM_POPS file to be read in in a restart calculation.
        real(dp) , intent(in) :: Norm_1RDM
        logical , intent(in) :: tNormalise
        integer :: i, j, iSpat, jSpat
        integer :: OneRDM_unit
 
        if(tNormalise) then
            ! Haven't got the capabilities to produce multiple 1-RDMs yet.
            write(6,*) 'Writing out the *normalised* 1 electron density matrix to file'
            call neci_flush(6)
            OneRDM_unit = get_free_unit()
            OPEN(OneRDM_unit,file='OneRDM',status='unknown')
        else
            ! Only every write out 1 of these at the moment.
            write(6,*) 'Writing out the *unnormalised* 1 electron density matrix to file for reading in'
            call neci_flush(6)
            OneRDM_unit = get_free_unit()
            OPEN(OneRDM_unit,file='OneRDM_POPS',status='unknown',form='unformatted')
        endif
 
        ! Currently always printing 1-RDM in spin orbitals.
        do i = 1, nBasis
            do j = 1, nBasis
                if(tStoreSpinOrbs) then
                    if(NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)).ne.0.0_dp) then 
                        if(tNormalise.and.((i.le.j).or.tHF_Ref_Explicit.or.tHF_S_D_Ref)) then
                            WRITE(6,*) "Written * 1.75"
                            write(OneRDM_unit,"(2I6,G25.17)") i,j, & 
                                NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) * Norm_1RDM
                        elseif(.not.tNormalise) then
                            ! For the pops, we haven't made the 1-RDM hermitian yet, 
                            ! so print both the 1-RDM(i,j) and 1-RDM(j,i) elements.
                            ! This is written in binary.
                            write(OneRDM_unit) i,j,NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))
                        endif
                    endif
                else
                    iSpat = gtID(i)
                    jSpat = gtID(j)
                    if(NatOrbMat(SymLabelListInv_rot(iSpat),SymLabelListInv_rot(jSpat)).ne.0.0_dp) then 
                        if(tNormalise.and.((i.le.j).or.tHF_Ref_Explicit.or.tHF_S_D_Ref)) then
                            if(((mod(i,2).eq.0).and.(mod(j,2).eq.0)).or.&
                                ((mod(i,2).ne.0).and.(mod(j,2).ne.0))) then
                                write(OneRDM_unit,"(2I6,G25.17)") i,j, & 
                                    ( NatOrbMat(SymLabelListInv_rot(iSpat),SymLabelListInv_rot(jSpat)) &
                                                                    * Norm_1RDM ) / 2.0_dp
                            endif
                        elseif(.not.tNormalise) then
                            ! The popsfile can be printed in spatial orbitals.
                            if((mod(i,2).eq.0).and.(mod(j,2).eq.0)) then
                                write(OneRDM_unit) iSpat,jSpat, & 
                                    NatOrbMat(SymLabelListInv_rot(iSpat),SymLabelListInv_rot(jSpat)) 
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
        real(dp) , intent(out) :: Norm_2RDM_Inst, Norm_2RDM
        real(dp) :: AllAccumRDMNorm_Inst
        logical :: tmake_herm
        real(dp), allocatable :: AllNodes_aaaa_RDM(:,:)
        real(dp), allocatable :: AllNodes_abab_RDM(:,:)
        real(dp), allocatable :: AllNodes_abba_RDM(:,:)
        integer :: ierr

!        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity

        ! If Iter = 0, this means we have just read in the TwoRDM_POPS_a*** matrices into All_a***_RDM, and 
        ! just want to calculate the old energy.
        ! Don't need to do all this stuff here, because a***_RDM will be empty.
        if(Iter.ne.0) then

            ALLOCATE(AllNodes_aaaa_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
            ALLOCATE(AllNodes_abba_RDM(((SpatOrbs*(SpatOrbs-1))/2),((SpatOrbs*(SpatOrbs-1))/2)),stat=ierr)
            ALLOCATE(AllNodes_abab_RDM(((SpatOrbs*(SpatOrbs+1))/2),((SpatOrbs*(SpatOrbs+1))/2)),stat=ierr)

            CALL MPISumAll(aaaa_RDM(:,:),AllNodes_aaaa_RDM(:,:))
            aaaa_RDM(:,:) = AllNodes_aaaa_RDM(:,:)

            CALL MPISumAll(abab_RDM(:,:),AllNodes_abab_RDM(:,:))
            abab_RDM(:,:) = AllNodes_abab_RDM(:,:)

            CALL MPISumAll(abba_RDM(:,:),AllNodes_abba_RDM(:,:))
            abba_RDM(:,:) = AllNodes_abba_RDM(:,:)

            DEALLOCATE(AllNodes_aaaa_RDM)
            DEALLOCATE(AllNodes_abab_RDM)
            DEALLOCATE(AllNodes_abba_RDM)

            !CALL MPISum_inplace(aaaa_RDM(:,:))
            !CALL MPISum_inplace(abab_RDM(:,:))
            !CALL MPISum_inplace(abba_RDM(:,:))

            ! The TwoElRDM on the root is now the sum of all 'instantaneous' RDMs (summed over 
            ! the energy update cycle).
            ! Whereas AllTwoElRDM is accumulated over the entire run.

            if(iProcIndex.eq.0) then
                All_aaaa_RDM(:,:) = All_aaaa_RDM(:,:) + aaaa_RDM(:,:) 
                All_abab_RDM(:,:) = All_abab_RDM(:,:) + abab_RDM(:,:) 
                All_abba_RDM(:,:) = All_abba_RDM(:,:) + abba_RDM(:,:)
            endif

            AllAccumRDMNorm_Inst = 0.0_dp
            if(tHF_S_D_Ref.or.tHF_Ref_Explicit.or.tHF_S_D) then
                CALL MPIReduce(AccumRDMNorm_Inst,MPI_SUM,AllAccumRDMNorm_Inst)
                AllAccumRDMNorm = AllAccumRDMNorm + AllAccumRDMNorm_Inst
            endif
        endif

        if(iProcIndex.eq.0) then

            ! Calculate the normalisations.
            call calc_2e_norms(AllAccumRDMNorm_Inst, Norm_2RDM_Inst, Norm_2RDM)

            ! Print out the relevant 2-RDMs.
            if( tFinalRDMEnergy .or. &
                ( tWriteMultRDMs .and. (mod((Iter - IterRDMStart)+1,IterWriteRDMs).eq.0) ) ) then

                tmake_herm = .false.

                if(tFinalRDMEnergy) then
                    ! Only ever want to print the POPS 2-RDMs (for reading in) at the end.
                    if(twrite_RDMs_to_read) call Write_out_2RDM(Norm_2RDM,.false.,.false.)

!                    ! We also don't want to make the 2-RDMs hermitian until the end, so that we can 
!                    ! get the hermiticity error from the final matrix.
!                    if(.not.(tHF_Ref_Explicit.or.tHF_S_D_Ref)) then
!                        call make_2e_rdm_hermitian(Norm_2RDM, Max_Error_Hermiticity, Sum_Error_Hermiticity)
!
!                        write(6,'(A29,F30.20)') ' MAX ABS ERROR IN HERMITICITY', Max_Error_Hermiticity
!                        write(6,'(A29,F30.20)') ' SUM ABS ERROR IN HERMITICITY', Sum_Error_Hermiticity
!                    endif

                    if(.not.(tHF_Ref_Explicit.or.tHF_S_D_Ref)) tmake_herm = .true.

                endif

                ! This writes out the normalised, hermitian 2-RDMs.
                if(twrite_normalised_RDMs) call Write_out_2RDM(Norm_2RDM,.true.,tmake_herm)

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
        real(dp) , intent(in) :: AllAccumRDMNorm_Inst
        real(dp) , intent(out) :: Norm_2RDM_Inst, Norm_2RDM
        integer :: i

        ! Find the current, unnormalised trace of each matrix.
        ! TODO: This can be merged into the spin averaging when everything is working.

        Trace_2RDM_Inst = 0.0_dp
        Trace_2RDM = 0.0_dp

        do i = 1, ((SpatOrbs*(SpatOrbs+1))/2)
            if(i.le.((SpatOrbs*(SpatOrbs-1))/2)) then
                Trace_2RDM_Inst = Trace_2RDM_Inst + aaaa_RDM(i,i)
                Trace_2RDM = Trace_2RDM + All_aaaa_RDM(i,i)
            endif
            Trace_2RDM_Inst = Trace_2RDM_Inst + abab_RDM(i,i)
            Trace_2RDM = Trace_2RDM + All_abab_RDM(i,i)
        enddo

        Norm_2RDM_Inst = 0.0_dp
        Norm_2RDM = 0.0_dp

        IF(tHF_S_D_Ref.or.tHF_Ref_Explicit.or.tHF_S_D) THEN
            Norm_2RDM_Inst = 1.0_dp / AllAccumRDMNorm_Inst
            Norm_2RDM = 1.0_dp / AllAccumRDMNorm
        ELSE
            ! Sum of diagonal elements of 2 electron RDM must equal number of 
            ! pairs of electrons, = NEl ( NEl - 1 ) / 2
            Norm_2RDM_Inst = ( (0.50_dp * (REAL(NEl,dp) * (REAL(NEl,dp) - 1.0_dp))) / Trace_2RDM_Inst )
            Norm_2RDM = ( (0.50_dp * (REAL(NEl,dp) * (REAL(NEl,dp) - 1.0_dp))) / Trace_2RDM )
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
        real(dp) , intent(in) :: Norm_2RDM
        real(dp) , intent(out) :: Max_Error_Hermiticity, Sum_Error_Hermiticity 
        integer :: i, j
        real(dp) :: Temp

        Max_Error_Hermiticity = 0.0_dp
        Sum_Error_Hermiticity = 0.0_dp

        do i = 1, ((SpatOrbs*(SpatOrbs+1))/2)
            do j = i+1, ((SpatOrbs*(SpatOrbs+1))/2)

                if((i.le.((SpatOrbs*(SpatOrbs-1))/2)).and.(j.le.((SpatOrbs*(SpatOrbs-1))/2))) then

                    IF((abs((All_aaaa_RDM(i,j)*Norm_2RDM)-(All_aaaa_RDM(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                        Max_Error_Hermiticity = abs((All_aaaa_RDM(i,j)*Norm_2RDM)-(All_aaaa_RDM(j,i)*Norm_2RDM))

                    Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                            abs((All_aaaa_RDM(i,j)*Norm_2RDM)-(All_aaaa_RDM(j,i)*Norm_2RDM))

                    Temp = (All_aaaa_RDM(i,j) + All_aaaa_RDM(j,i)) / 2.0_dp

                    All_aaaa_RDM(i,j) = Temp
                    All_aaaa_RDM(j,i) = Temp

                    IF((abs((All_abba_RDM(i,j)*Norm_2RDM)-(All_abba_RDM(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                        Max_Error_Hermiticity = abs((All_abba_RDM(i,j)*Norm_2RDM)-(All_abba_RDM(j,i)*Norm_2RDM))

                    Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                            abs((All_abba_RDM(i,j)*Norm_2RDM)-(All_abba_RDM(j,i)*Norm_2RDM))

                    Temp = (All_abba_RDM(i,j) + All_abba_RDM(j,i)) / 2.0_dp

                    All_abba_RDM(i,j) = Temp
                    All_abba_RDM(j,i) = Temp
                endif

                IF((abs((All_abab_RDM(i,j)*Norm_2RDM)-(All_abab_RDM(j,i)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                    Max_Error_Hermiticity = abs((All_abab_RDM(i,j)*Norm_2RDM)-(All_abab_RDM(j,i)*Norm_2RDM))

                Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                        abs((All_abab_RDM(i,j)*Norm_2RDM)-(All_abab_RDM(j,i)*Norm_2RDM))

                Temp = (All_abab_RDM(i,j) + All_abab_RDM(j,i)) / 2.0_dp

                All_abab_RDM(i,j) = Temp
                All_abab_RDM(j,i) = Temp


            enddo
        enddo

    end subroutine make_2e_rdm_hermitian


    subroutine Write_out_2RDM(Norm_2RDM,tNormalise,tmake_herm)
! Writes out the 2-RDMs.  If tNormalise is true, we print the normalised (hermitian) matrix.    
! Otherwise we print the unnormalised 2-RDMs, and we print (in binary) both 2-RDM(Ind1,Ind2) 
! and 2-RDM(Ind2,Ind1) because this matrix wont be hermitian.

! While, for instance, the TwoRDM_aaaa so far has actually been a sum of the aaaa elements and 
! the bbbb elements.  We only want to print the aaaa elements.
        real(dp) , intent(in) :: Norm_2RDM
        logical , intent(in) :: tNormalise, tmake_herm
        real(dp) :: Tot_Spin_Projection, SpinPlus, SpinMinus
        real(dp) :: ParityFactor,Divide_Factor 
        integer :: i, j, a, b, Ind1_aa, Ind1_ab, Ind2_aa, Ind2_ab
        integer :: aaaa_RDM_unit, abab_RDM_unit, abba_RDM_unit, No_Herm_Elements
        character(255) :: TwoRDM_aaaa_name, TwoRDM_abab_name, TwoRDM_abba_name
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity, Sum_Herm_Percent 
        real(dp) :: Temp


        if(tNormalise) then
            write(6,*) 'Writing out the *normalised* 2 electron density matrix to file'
            call neci_flush(6)
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
            call neci_flush(6)
            aaaa_RDM_unit = get_free_unit()
            OPEN(aaaa_RDM_unit,file='TwoRDM_POPS_aaaa',status='unknown',form='unformatted')
            abab_RDM_unit = get_free_unit()
            OPEN(abab_RDM_unit,file='TwoRDM_POPS_abab',status='unknown',form='unformatted')
            abba_RDM_unit = get_free_unit()
            OPEN(abba_RDM_unit,file='TwoRDM_POPS_abba',status='unknown',form='unformatted')
        endif
        
        Max_Error_Hermiticity = 0.0_dp
        Sum_Error_Hermiticity = 0.0_dp
        Sum_Herm_Percent = 0.0_dp
        No_Herm_Elements = 0
        Tot_Spin_Projection = 0.0_dp
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

                            if( (All_aaaa_RDM(Ind1_aa,Ind2_aa).ne.0.0_dp).or.&
                                (All_aaaa_RDM(Ind2_aa,Ind1_aa).ne.0.0_dp) )then
                                ! If we're normalising (and have made the matrix hermitian) we only 
                                ! need to write out Ind1 < Ind2.
                                ! Otherwise we print out Ind1, Ind2 and Ind2, Ind1 so we can 
                                ! find the hermiticity error in the final matrix (after all runs).
                                if(tNormalise.and.((Ind1_aa.le.Ind2_aa).or.tHF_Ref_Explicit.or.tHF_S_D_Ref)) then
                                    
                                    IF((abs((All_aaaa_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                - (All_aaaa_RDM(Ind2_aa,Ind1_aa)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                        Max_Error_Hermiticity = abs((All_aaaa_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (All_aaaa_RDM(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                    Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                            abs((All_aaaa_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (All_aaaa_RDM(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                    Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                        (abs((All_aaaa_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                - (All_aaaa_RDM(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                        (abs((All_aaaa_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                + (All_aaaa_RDM(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                    No_Herm_Elements = No_Herm_Elements + 1                                                        

                                    if(tmake_herm) then                                                            
                                        Temp = (All_aaaa_RDM(Ind1_aa,Ind2_aa) + All_aaaa_RDM(Ind2_aa,Ind1_aa)) / 2.0_dp

                                        All_aaaa_RDM(Ind1_aa,Ind2_aa) = Temp
                                        All_aaaa_RDM(Ind2_aa,Ind1_aa) = Temp
                                    endif

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

                            if( (All_abba_RDM(Ind1_aa,Ind2_aa).ne.0.0_dp).or.&
                                (All_abba_RDM(Ind2_aa,Ind1_aa).ne.0.0_dp) ) then
                                if(tNormalise.and.((Ind1_aa.le.Ind2_aa).or.tHF_Ref_Explicit.or.tHF_S_D_Ref)) then

                                    IF((abs((All_abba_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                - (All_abba_RDM(Ind2_aa,Ind1_aa)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                        Max_Error_Hermiticity = abs((All_abba_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                - (All_abba_RDM(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                    Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                            abs((All_abba_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                                    - (All_abba_RDM(Ind2_aa,Ind1_aa)*Norm_2RDM))

                                    Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                        (abs((All_abba_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                        - (All_abba_RDM(Ind2_aa,Ind1_aa)*Norm_2RDM)) / &
                                                        (abs((All_abba_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM) &
                                                        + (All_abba_RDM(Ind2_aa,Ind1_aa)*Norm_2RDM)) / 2.0_dp) )
                                    No_Herm_Elements = No_Herm_Elements + 1                                                        

                                    if(tmake_herm) then                                                            
                                        Temp = (All_abba_RDM(Ind1_aa,Ind2_aa) + All_abba_RDM(Ind2_aa,Ind1_aa)) / 2.0_dp

                                        All_abba_RDM(Ind1_aa,Ind2_aa) = Temp
                                        All_abba_RDM(Ind2_aa,Ind1_aa) = Temp
                                    endif

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

                        if( (All_abab_RDM(Ind1_ab,Ind2_ab).ne.0.0_dp).or.&
                            (All_abab_RDM(Ind2_ab,Ind1_ab).ne.0.0_dp) ) then

                            IF((abs((All_abab_RDM(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                        - (All_abab_RDM(Ind2_ab,Ind1_ab)*Norm_2RDM))).gt.Max_Error_Hermiticity) &
                                Max_Error_Hermiticity = abs((All_abab_RDM(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                        - (All_abab_RDM(Ind2_ab,Ind1_ab)*Norm_2RDM))

                            Sum_Error_Hermiticity = Sum_Error_Hermiticity +     &
                                                    abs((All_abab_RDM(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                        - (All_abab_RDM(Ind2_ab,Ind1_ab)*Norm_2RDM))

                            Sum_Herm_Percent = Sum_Herm_Percent +   &
                                                (abs((All_abab_RDM(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                    - (All_abab_RDM(Ind2_ab,Ind1_ab)*Norm_2RDM)) / &
                                                (abs((All_abab_RDM(Ind1_ab,Ind2_ab)*Norm_2RDM) &
                                                    + (All_abab_RDM(Ind2_ab,Ind1_ab)*Norm_2RDM)) / 2.0_dp) )
                            No_Herm_Elements = No_Herm_Elements + 1                                                        

                            if(tmake_herm) then                                                            
                                Temp = (All_abab_RDM(Ind1_ab,Ind2_ab) + All_abab_RDM(Ind2_ab,Ind1_ab)) / 2.0_dp

                                All_abab_RDM(Ind1_ab,Ind2_ab) = Temp
                                All_abab_RDM(Ind2_ab,Ind1_ab) = Temp
                            endif

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

        if(tNormalise.and.(.not.(tHF_Ref_Explicit.or.tHF_S_D_Ref))) then
            write(6,'(I15,F30.20,A20,A39)') Iter+PreviousCycles, Max_Error_Hermiticity, &
                                            '( Iteration,',' MAX ABS ERROR IN HERMITICITY )'
            write(6,'(I15,F30.20,A20,A39)') Iter+PreviousCycles, Sum_Error_Hermiticity, &
                                            '( Iteration,',' SUM ABS ERROR IN HERMITICITY )'
            write(6,'(I15,F30.20,A20,A51)') Iter+PreviousCycles, Sum_Herm_Percent/real(No_Herm_Elements,dp), &
                                            '( Iteration,',' AVERAGE ABS PERCENTAGE HERMITICITY ERROR )'
        endif

!        Tot_Spin_Projection = Tot_Spin_Projection + (3.0_dp * real(NEl,dp))
!! Tot_Spin_Projection is now equal to 4S(S+1) - find S. 
!        Tot_Spin_Projection = Tot_Spin_Projection/4.0_dp
!        if((1.0_dp + 4.0_dp*Tot_Spin_Projection).le.0) then
!            call Warning_neci('Write_out_1and_2RDM',"Complex spin calculated from density matrices!")
!        else
!            SpinPlus = (-1.0_dp + sqrt(1.0_dp + 4.0_dp*Tot_Spin_Projection))/2.0_dp
!            SpinMinus = (-1.0_dp + sqrt(1.0_dp + 4.0_dp*Tot_Spin_Projection))/2.0_dp
!!            write(6,*) 'SpinPlus',SpinPlus
!!            write(6,*) 'SpinMinus',SpinMinus
!        endif

!        if(RDMExcitLevel.eq.3) then
!            write(6,*) ''
!            write(6,'(A22,F30.20)') ' TOTAL SPIN PROJECTION', Max(SpinPlus,SpinMinus) 
!            write(6,*) ''
!        endif

    end subroutine Write_out_2RDM


    SUBROUTINE Calc_Energy_from_RDM(Norm_2RDM)
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
        real(dp), intent(out) :: Norm_2RDM
        real(dp) :: Norm_2RDM_Inst
        INTEGER :: i,j,a,b,Ind1_aa,Ind1_ab,Ind2_aa,Ind2_ab,ierr
        INTEGER :: iSpin, jSpin, error
        REAL(dp) :: RDMEnergy_Inst, RDMEnergy, Coul, Exch, Parity_Factor 
        REAL(dp) :: Trace_2RDM_New, RDMEnergy1, RDMEnergy2

        CALL set_timer(RDMEnergy_Time,30)

        Trace_2RDM_New = 0.0_dp

        RDMEnergy_Inst = 0.0_dp
        RDMEnergy1 = 0.0_dp
        RDMEnergy2 = 0.0_dp
        RDMEnergy = 0.0_dp
    
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
                            Coul = REAL(UMAT(UMatInd(i,j,a,b,0,0)),dp)
                            Exch = REAL(UMAT(UMatInd(i,j,b,a,0,0)),dp)

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
            call neci_flush(Energies_unit)

            if(tFinalRDMEnergy) then
                write(6,*) 'Trace of 2-el-RDM before normalisation : ',Trace_2RDM
                write(6,*) 'Trace of 2-el-RDM after normalisation : ',Trace_2RDM_New
                write(6,*) 'Energy contribution from the 1-RDM: ',RDMEnergy1
                write(6,*) 'Energy contribution from the 2-RDM: ',RDMEnergy2
                write(6,'(A64,F30.20)') ' *TOTAL ENERGY* CALCULATED USING THE *REDUCED &
                                            &DENSITY MATRICES*:',RDMEnergy
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
    
    SUBROUTINE Calc_Lagrangian_from_RDM(Norm_1RDM,Norm_2RDM)
! CMO
! This routine takes the 1 electron and 2 electron reduced density matrices 
! and calculated the Lagrangian term, X, required for the calculation of forces.    
! The equation for X is as follows:
!
!   X_pq = Sum_r[h_pr 1RDM_qr] + 0.5*Sum_rst[(pr|st)[2RDM_qrst + 2RDM_rqst]]
!        where 2RDM is defined in chemical notation sense: 2RDM_ijkl = <Psi| a_i+ a_k+ a_l a_j|Psi>
!
        USE IntegralsData , only : UMAT
        USE UMatCache , only : UMatInd
        USE RotateOrbsMod , only : SymLabelList2_rot
        USE UMatCache , only : GTID
        USE Logging, only : tDumpForcesInfo, tPrintLagrangian
        implicit none
        real(dp), intent(in) :: Norm_2RDM
        real(dp), intent(in) :: Norm_1RDM
        real(dp) :: Norm_2RDM_Inst
        INTEGER :: p,q,r,s,t,ierr,stat
        INTEGER :: pSpin, qSpin, rSpin, error
        REAL(dp) :: RDMEnergy_Inst, RDMEnergy, Coul, Exch, Parity_Factor 
        REAL(dp) :: Trace_2RDM_New, RDMEnergy1, RDMEnergy2
        REAL(dp) :: qrst, rqst
        REAL(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity, Temp

        ALLOCATE(Lagrangian(SpatOrbs,SpatOrbs),stat=ierr)
        Lagrangian(:,:)=0.D0
        
        ! We will begin by calculating the Lagrangian in chemical notation - we will explicitely calculate
        ! both halves (X_pq and X_qp) in order to check if it is symmetric or not.
        !  - a symmetric Lagrangian is important in allowing us to use the derivative overlap matrix rather than the 
        !    coupled-perturbed coefficients when calculating the Forces later on 
        !    (see Sherrill, Analytic Gradients of CI Energies eq38)

        write(6,*) ''
        write(6,*) 'Calculating the Lagrangian X from the final density matrices'

        !Calculating the Lagrangian X in terms of spatial orbitals
        if(iProcIndex.eq.0) then

            do p = 1, SpatOrbs  !Run over spatial orbitals
                pSpin = 2*p    ! Picks out beta component
                do q = 1, SpatOrbs  !Want to calculate X(p,q) separately from X(q,p) for now to see if we're symmetric
                    do r = 1, SpatOrbs
                        rSpin=2*r
                        ! Adding in contributions effectively from the 1-RDM
                        ! We made sure earlier that the 1RDM is contructed, so we can call directly from this
                        if(tStoreSpinOrbs) then
                            ! Include both aa and bb contributions
                            Lagrangian(p,q)=Lagrangian(p,q)+2.0_dp*(NatOrbMat(SymLabelListInv_rot(2*q),SymLabelListInv_rot(2*r)))* &
                                                                REAL(TMAT2D(pSpin,rSpin),8)*Norm_1RDM
                        else
                            !We will be here most often (?)
                            Lagrangian(p,q)=Lagrangian(p,q)+NatOrbMat(SymLabelListInv_rot(q),SymLabelListInv_rot(r))* &
                                                                REAL(TMAT2D(pSpin,rSpin),8)*Norm_1RDM
                        endif

                        do s = 1, SpatOrbs
                
                            ! Here we're looking to start adding on the contributions from the 2-RDM and the 2-el integrals
                            ! For X(p,q), these have the form 0.5*Sum_rst[(pr|st)[2RDM_qrst + 2RDM_rqst]]
                            
                            ! NOTE: In some notations, the factor of a half goes *inside* the 2RDM - 
                            ! consistent with Yamaguchi etc (see eq 11 in Sherrill, Analytic Gradients
                            ! of CI energies).  We will keep it outside for now, consistent with the storage
                            ! within neci

                            do t=1, SpatOrbs
                                
                                !Integral (pr|st) = <ps|rt>
                                !Give indices in PHYSICAL NOTATION
                                !NB, FCIDUMP is labelled in chemical notation
                                Coul = REAL(UMAT(UMatInd(p,s,r,t,0,0)),8)

                                qrst=Find_Spatial_2RDM_Chem(q,r,s,t, Norm_2RDM)
                                rqst=Find_Spatial_2RDM_Chem(r,q,s,t, Norm_2RDM)
                                
                                Lagrangian(p,q)=Lagrangian(p,q) + 0.5_dp*Coul*(qrst+rqst)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
       
            !! Now symmetrise (make hermitian, such that X_pq = X_qp) the Lagrangian X

            Max_Error_Hermiticity = 0.D0
            Sum_Error_Hermiticity = 0.D0
            do p = 1, SpatOrbs
                do q = p, SpatOrbs
                    IF(abs(Lagrangian(p,q) - Lagrangian(q,p)).gt.Max_Error_Hermiticity) &
                        Max_Error_Hermiticity = abs(Lagrangian(p,q)-Lagrangian(q,p))

                    Sum_Error_Hermiticity = Sum_Error_Hermiticity+abs(Lagrangian(p,q) - Lagrangian(q,p))

                    Temp = (Lagrangian(p,q) + Lagrangian(q,p))/2.D0
                    Lagrangian(p,q) = Temp
                    Lagrangian(q,p) = Temp
                enddo
            enddo

            ! Output the hermiticity errors.
            write(6,'(A40,F30.20)') ' MAX ABS ERROR IN Lagrangian HERMITICITY', Max_Error_Hermiticity
            write(6,'(A40,F30.20)') ' SUM ABS ERROR IN Lagrangian HERMITICITY', Sum_Error_Hermiticity

        endif

    end subroutine

    subroutine calc_1RDM_energy(i,j,a,iSpin,jSpin,Norm_2RDM,Norm_2RDM_Inst,&
                                                        RDMEnergy_Inst,RDMEnergy1)
    ! This routine calculates the 1-RDM part of the RDM energy, and constructs the 
    ! 1-RDM if required for diagonalisation or something.
    ! gamma(i,j) = [1/(NEl - 1)] * SUM_a Gamma(i,a,j,a) 
    ! want to calculate:    gamma(i,j) * h_ij
    ! h_ij => TMAT2D(iSpin,jSpin)
    
        USE OneEInts , only : TMAT2D
        USE Logging , only : tDiagRDM, tDumpForcesInfo, tDipoles
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
                                                * REAL(TMAT2D(iSpin,jSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) )
            RDMEnergy1 = RDMEnergy1 + ( (All_abab_RDM(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM) &
                                                * REAL(TMAT2D(iSpin,jSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) )

            if((tDiagRDM.or.tPrint1RDM.or.tDumpForcesInfo.or.tDipoles)  &
                .and.tFinalRDMEnergy) then                                                
                NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = &
                            NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) &
                                        + ( All_abab_RDM(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) ) 
            endif

            ! For Gamma elements corresponding to 1-RDMs ( Gamma(i,a,j,a) ), we're only considering 
            ! i =< j and therefore we need to sum in the opposite contribution too.
            if(Ind1_1e_ab.ne.Ind2_1e_ab) then                                                                
                RDMEnergy_Inst = RDMEnergy_Inst + ( (abab_RDM(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM_Inst) &
                                                    * REAL(TMAT2D(jSpin,iSpin),dp) &
                                                    * (1.0_dp / real(NEl - 1,dp)) )
                RDMEnergy1 = RDMEnergy1 + ( (All_abab_RDM(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM) &
                                                    * REAL(TMAT2D(jSpin,iSpin),dp) &
                                                    * (1.0_dp / real(NEl - 1,dp)) )

                if((tDiagRDM.or.tPrint1RDM.or.tDumpForcesInfo.or.tDipoles).and.tFinalRDMEnergy) then  
                    NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = &
                            NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) &
                                        + ( All_abab_RDM(Ind2_1e_ab,Ind1_1e_ab) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) ) 
                endif

            endif

            ! But since we're running over all a, i a and a i will both be counted, but i i only once 
            ! (whereas it should be counted twice).
            if((i.eq.j).and.(i.eq.a)) then                                                                
                RDMEnergy_Inst = RDMEnergy_Inst + ( (abab_RDM(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM_Inst) &
                                            * REAL(TMAT2D(iSpin,jSpin),dp) &
                                            * (1.0_dp / real(NEl - 1,dp)) )
                RDMEnergy1 = RDMEnergy1 + ( (All_abab_RDM(Ind1_1e_ab,Ind2_1e_ab) * Norm_2RDM) &
                                            * REAL(TMAT2D(iSpin,jSpin),dp) &
                                            * (1.0_dp / real(NEl - 1,dp)) )

                if((tDiagRDM.or.tPrint1RDM.or.tDumpForcesInfo.or.tDipoles).and.tFinalRDMEnergy) then 
                    NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = &
                            NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) &
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
                                                    * REAL(TMAT2D(iSpin,jSpin),dp) &
                                                    * (1.0_dp / real(NEl - 1,dp)) )
                RDMEnergy1 = RDMEnergy1 - ( (All_abba_RDM(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                    * REAL(TMAT2D(iSpin,jSpin),dp) &
                                                    * (1.0_dp / real(NEl - 1,dp)) )

                if((tDiagRDM.or.tPrint1RDM.or.tDumpForcesInfo.or.tDipoles).and.tFinalRDMEnergy) then 
                    NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = &
                                NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) &
                                            - ( All_abba_RDM(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                    * (1.0_dp / real(NEl - 1,dp)) ) 
                endif

                if(Ind1_1e_aa.ne.Ind2_1e_aa) then
                    RDMEnergy_Inst = RDMEnergy_Inst - ( (abba_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM_Inst) &
                                                        * REAL(TMAT2D(iSpin,jSpin),dp) &
                                                        * (1.0_dp / real(NEl - 1,dp)) )
                    RDMEnergy1 = RDMEnergy1 - ( (All_abba_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                        * REAL(TMAT2D(iSpin,jSpin),dp) &
                                                        * (1.0_dp / real(NEl - 1,dp)) )

                    if((tDiagRDM.or.tPrint1RDM.or.tDumpForcesInfo.or.tDipoles).and.tFinalRDMEnergy) then 
                        NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = &
                                NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) &
                                            - ( ( All_abba_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM &
                                                    * (1.0_dp / real(NEl - 1,dp)) ))
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
                                                * REAL(TMAT2D(iSpin,jSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
            RDMEnergy1 = RDMEnergy1 + ( (All_aaaa_RDM(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM) &
                                                * REAL(TMAT2D(iSpin,jSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )

            if((tDiagRDM.or.tPrint1RDM.or.tDumpForcesInfo.or.tDipoles).and.tFinalRDMEnergy) then 
                NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = &
                            NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) &
                                        + ( All_aaaa_RDM(Ind1_1e_aa,Ind2_1e_aa) * Norm_2RDM &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor ) 
            endif

            if(Ind1_1e_aa.ne.Ind2_1e_aa) then
                RDMEnergy_Inst = RDMEnergy_Inst + ( (aaaa_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM_Inst) &
                                                * REAL(TMAT2D(jSpin,iSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )
                RDMEnergy1 = RDMEnergy1 + ( (All_aaaa_RDM(Ind2_1e_aa,Ind1_1e_aa) * Norm_2RDM) &
                                                * REAL(TMAT2D(jSpin,iSpin),dp) &
                                                * (1.0_dp / real(NEl - 1,dp)) * Parity_Factor )

                if((tDiagRDM.or.tPrint1RDM.or.tDumpForcesInfo.or.tDipoles).and.tFinalRDMEnergy) then 
                    NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) = &
                            NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(i)) &
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
        integer :: ierr
        REAL(dp) :: SumDiag
        CHARACTER(len=*), PARAMETER :: this_routine='find_nat_orb_occ_numbers'

        IF(iProcIndex.eq.0) THEN
            
            ! Diagonalises the 1-RDM.  NatOrbMat goes in as the 1-RDM, comes out as the 
            ! eigenvector of the 1-RDM (the matrix transforming the MO's into the NOs).
            CALL DiagRDM(SumDiag)

            ! Writes out the NO occupation numbers and evectors to files.
            call write_evales_and_transform_mat(SumDiag)

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

                CALL neci_flush(6)

                CALL WRITEBASIS(6,G1,NoOrbs,ARR,BRR)

            ENDIF
        ENDIF

    end subroutine find_nat_orb_occ_numbers

    subroutine write_evales_and_transform_mat(SumDiag)
        real(dp) , intent(in) :: SumDiag
        integer :: i, j, Evalues_unit, NatOrbs_unit, jSpat, jInd, NO_Number
        REAL(dp) :: Corr_Entropy, Norm_Evalues, SumN_NO_Occ
        logical :: tNegEvalue, tWrittenEvalue

        if(tStoreSpinOrbs) then
            Norm_Evalues = SumDiag/REAL(NEl)
        else
            Norm_Evalues = 2.0_dp*(SumDiag/REAL(NEl))
        endif

        ! Write out normalised evalues to file and calculate the correlation entropy.
        Corr_Entropy = 0.0_dp
        Evalues_unit = get_free_unit()
        OPEN(Evalues_unit,file='NO_OCC_NUMBERS',status='unknown')

        WRITE(Evalues_unit,'(A)') '# NOs (natural orbitals) ordered by occupation number' 
        WRITE(Evalues_unit,'(A)') '# MOs (HF orbitals) ordered by energy' 
        WRITE(Evalues_unit,'(A1,A5,A30,A20,A30)') '#','NO','NO OCCUPATION NUMBER','MO','MO OCCUPATION NUMBER'
        tNegEvalue = .false.
        SumN_NO_Occ = 0.0_dp
        NO_Number = 1
        do i=1,NoOrbs
            if(tStoreSpinOrbs) then
                WRITE(Evalues_unit,'(I6,G35.17,I15,G35.17)') i,Evalues(i)/Norm_Evalues, &
                                                                BRR(i), Rho_ii(i)
                if(Evalues(i).gt.0.0_dp) then
                    Corr_Entropy = Corr_Entropy - ( abs(Evalues(i)/ Norm_Evalues) &
                                                    * LOG(abs(Evalues(i)/ Norm_Evalues)) )
                else
                    tNegEvalue = .true.
                endif
                if(i.le.NEl) SumN_NO_Occ = SumN_NO_Occ + (Evalues(i)/Norm_Evalues)
            else
                WRITE(Evalues_unit,'(I6,G35.17,I15,G35.17)') (2*i)-1,Evalues(i)/Norm_Evalues, &
                                                            BRR((2*i)-1), Rho_ii(i)/2.0_dp
                if(Evalues(i).gt.0.0_dp) then
                    Corr_Entropy = Corr_Entropy - (2.0_dp * ( abs(Evalues(i)/Norm_Evalues) &
                                                    * LOG(abs(Evalues(i)/Norm_Evalues)) ) )
                else
                    tNegEvalue = .true.
                endif
                WRITE(Evalues_unit,'(I6,G35.17,I15,G35.17)') 2*i,Evalues(i)/Norm_Evalues, &
                                                            BRR(2*i), Rho_ii(i)/2.0_dp
                if(i.le.(NEl/2)) SumN_NO_Occ = SumN_NO_Occ + (2.0_dp * (Evalues(i)/Norm_Evalues))
            endif
        enddo
        close(Evalues_unit)

        write(6,'(A45,F30.20)') ' SUM OF THE N LARGEST NO OCCUPATION NUMBERS: ',SumN_NO_Occ
   
        WRITE(6,'(A20,F30.20)') ' CORRELATION ENTROPY', Corr_Entropy
        WRITE(6,'(A33,F30.20)') ' CORRELATION ENTROPY PER ELECTRON', Corr_Entropy / real(NEl,dp) 
        if(tNegEvalue) write(6,'(A)') ' WARNING: Negative NO occupation numbers found.'

        ! Write out the evectors to file.
        ! This is the matrix that transforms the molecular orbitals into the natural orbitals.
        ! Evalue(i) corresponds to Evector NatOrbsMat(1:nBasis,i)
        ! We just want the Evalues in the same order as above, but the 1:nBasis part (corresponding 
        ! to the molecular orbitals), needs to refer to the actual orbital labels.
        ! Want these orbitals to preferably be in order, run through the orbital, need the position 
        ! to find the corresponding NatOrbs element, use SymLabelListInv_rot
        if(.not.tNoNOTransform) then
            NatOrbs_unit = get_free_unit()
            OPEN(NatOrbs_unit,file='NO_TRANSFORM',status='unknown')
            write(NatOrbs_unit,'(2A6,2A30)') '#   MO','NO','Transform Coeff','NO OCC NUMBER'
            ! write out in terms of spin orbitals, all alpha then all beta.
            NO_Number = 1
            do i = 1, NoOrbs
                tWrittenEvalue = .false.
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
                    if(tWrittenEvalue) then
                        if(NatOrbMat(SymLabelListInv_rot(jInd),i).ne.0.0_dp) &
                            write(NatOrbs_unit,'(2I6,G35.17)') j,NO_Number,NatOrbMat(SymLabelListInv_rot(jInd),i)
                    else
                        if(NatOrbMat(SymLabelListInv_rot(jInd),i).ne.0.0_dp) then 
                            write(NatOrbs_unit,'(2I6,2G35.17)') j,NO_Number,NatOrbMat(SymLabelListInv_rot(jInd),i),&
                                                                                Evalues(i)/Norm_Evalues
                            tWrittenEvalue = .true.
                        endif
                    endif
                enddo
                NO_Number = NO_Number + 1
                if(.not.tStoreSpinOrbs) then
                    tWrittenEvalue = .false.
                    do j = 2, nBasis, 2
                        ! Here i corresponds to the natural orbital, and j to the molecular orbital.
                        ! i is actually the spin orbital in this case.
                        jSpat = gtID(j)
                        if(tWrittenEvalue) then
                            if(NatOrbMat(SymLabelListInv_rot(jSpat),i).ne.0.0_dp) &
                                write(NatOrbs_unit,'(2I6,G35.17)') j,NO_Number,&
                                                                NatOrbMat(SymLabelListInv_rot(jSpat),i)
                        else
                            if(NatOrbMat(SymLabelListInv_rot(jSpat),i).ne.0.0_dp) then
                                write(NatOrbs_unit,'(2I6,2G35.17)') j,NO_Number,&
                                                        NatOrbMat(SymLabelListInv_rot(jSpat),i), &
                                                        Evalues(i)/Norm_Evalues
                                tWrittenEvalue = .true.
                            endif
                        endif
                    enddo
                    NO_Number = NO_Number + 1
                endif
            enddo
            close(NatOrbs_unit)
        endif

    end subroutine write_evales_and_transform_mat

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
        REAL(dp) , intent(out) :: SumTrace
        REAL(dp) :: SumDiagTrace
        REAL(dp) , ALLOCATABLE :: WORK2(:),EvaluesSym(:),NOMSym(:,:)
        INTEGER :: ierr,i,j,spin,Sym,LWORK2,WORK2Tag,SymStartInd,NoSymBlock
        INTEGER :: EvaluesSymTag,NOMSymTag,k,MaxSym
        LOGICAL :: tDiffSym, tDiffLzSym
        CHARACTER(len=*), PARAMETER :: this_routine='DiagRDM'

! Test that we're not breaking symmetry.
! And calculate the trace at the same time.
        SumTrace=0.0_dp
        do i=1,NoOrbs
            do j=1,NoOrbs
                tDiffSym = .false.
                tDiffLzSym = .false.
                if(tStoreSpinOrbs) then
                    IF((INT(G1(SymLabelList2_rot(i))%sym%S).ne.&
                        INT(G1(SymLabelList2_rot(j))%sym%S))) tDiffSym = .true.
                    IF((INT(G1(SymLabelList2_rot(i))%Ml).ne.&
                        INT(G1(SymLabelList2_rot(j))%Ml))) tDiffLzSym = .true.
                else
                    IF((INT(G1(2*SymLabelList2_rot(i))%sym%S).ne.&
                        INT(G1(2*SymLabelList2_rot(j))%sym%S))) tDiffSym = .true.
                    IF((INT(G1(2*SymLabelList2_rot(i))%Ml).ne.&
                        INT(G1(2*SymLabelList2_rot(j))%Ml))) tDiffLzSym = .true.
                endif
                if(tDiffSym) then
                    IF(ABS(NatOrbMat(i,j)).ge.1.0E-15) THEN
                        WRITE(6,'(6A8,A20)') 'i','j','Label i','Label j','Sym i',&
                                                                'Sym j','Matrix value'
                        if(tStoreSpinOrbs) then                                                              
                            WRITE(6,'(6I3,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    INT(G1(SymLabelList2_rot(i))%sym%S),&
                                    INT(G1(SymLabelList2_rot(j))%sym%S),NatOrbMat(i,j)
                        else
                            WRITE(6,'(6I3,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    INT(G1(2*SymLabelList2_rot(i))%sym%S),&
                                    INT(G1(2*SymLabelList2_rot(j))%sym%S),NatOrbMat(i,j)
                        endif
                        IF(tUseMP2VarDenMat) THEN
                            WRITE(6,*) '**WARNING** - There is a non-zero NatOrbMat &
                            &value between orbitals of different symmetry.'
                            WRITE(6,*) 'These elements will be ignored, and the symmetry &
                            &maintained in the final transformation matrix.'
                        ELSE
                            write(6,*) 'k,SymLabelList2_rot(k),SymLabelListInv_rot(k)'
                            do k = 1,NoOrbs
                                write(6,*) k,SymLabelList2_rot(k),SymLabelListInv_rot(k)
                            enddo
                            call neci_flush(6)
                            CALL Stop_All(this_routine,'Non-zero NatOrbMat value between &
                            &different symmetries.')
                        ENDIF
                    ENDIF
                    NatOrbMat(i,j)=0.0_dp
                ENDIF
                if(tDiffLzSym) then
                    IF(ABS(NatOrbMat(i,j)).ge.1.0E-15) THEN
                        WRITE(6,'(6A8,A40)') 'i','j','Label i','Label j','Lz i',&
                                                                'Lz j','Matrix value'
                        if(tStoreSpinOrbs) then                                                              
                            WRITE(6,'(6I8,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    INT(G1(SymLabelList2_rot(i))%Ml),&
                                    INT(G1(SymLabelList2_rot(j))%Ml),NatOrbMat(i,j)
                        else
                            WRITE(6,'(6I8,F40.20)') i,j,SymLabelList2_rot(i),SymLabelList2_rot(j),&
                                    INT(G1(2*SymLabelList2_rot(i))%Ml),&
                                    INT(G1(2*SymLabelList2_rot(j))%Ml),NatOrbMat(i,j)
                        endif
                        write(6,'(A)') ' **WARNING** - There is a non-zero NatOrbMat element &
                        &between orbitals of different Lz symmetry.'
                    endif
                endif
            enddo
            SumTrace=SumTrace+NatOrbMat(i,i)
        enddo

        write(6,*) ''
        WRITE(6,*) 'Calculating eigenvectors and eigenvalues of NatOrbMat'
        CALL neci_flush(6)

! If we want to maintain the symmetry, we cannot have all the orbitals jumbled up when the 
! diagonaliser reorders the eigenvectors.
! Must instead feed each symmetry block in separately.
! This means that although the transformed orbitals are jumbled within the symmetry blocks, 
! the symmetry labels are all that are relevant and these are unaffected.
        Sym=0
        LWORK2=-1
        if(tStoreSpinOrbs) then
            if(tFixLz) then
                MaxSym = (16 * ( ( 2 * iMaxLz ) + 1 ) ) - 1
            else
                MaxSym = 15
            endif
        else
            if(tFixLz) then
                MaxSym = (8 * ( ( 2 * iMaxLz ) + 1 ) ) - 1
            else
                MaxSym = 7
            endif
        endif
        do while (Sym.le.MaxSym)

            NoSymBlock=SymLabelCounts2_rot(2,Sym+1)

            SymStartInd=SymLabelCounts2_rot(1,Sym+1)-1
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
!                    WRITE(6,'(I20)',advance='no') SymLabelList2_rot(SymStartInd+i)
!                enddo
!                write(6,*) ''
!                do i=1,NoSymBlock
!                    write(6,'(I20)',advance='no') SymLabelList2_rot(SymStartInd+i)
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
                NatOrbMat(SymStartInd+1,SymStartInd+1)=1.0_dp
!                WRITE(6,*) '*****'
!                WRITE(6,*) 'Symmetry ',Sym,' has only one orbital.'
!                WRITE(6,*) 'Copying diagonal element ,',SymStartInd+1,'to NatOrbMat'
            ENDIF

            Sym=Sym+1
        enddo

        WRITE(6,*) 'Matrix diagonalised'
        CALL neci_flush(6)

        SumDiagTrace=0.0_dp
        do i=1,NoOrbs
            SumDiagTrace=SumDiagTrace+Evalues(i)
        enddo
        IF((ABS(SumDiagTrace-SumTrace)).gt.1.0_dp) THEN
            WRITE(6,*) 'Sum of diagonal NatOrbMat elements : ',SumTrace
            WRITE(6,*) 'Sum of eigenvalues : ',SumDiagTrace
            WRITE(6,*) 'WARNING : &
            &The trace of the 1RDM matrix before diagonalisation is '
            write(6,*) 'not equal to that after.'
        ENDIF

! The MO's still correspond to SymLabelList2_rot.
! Although the NO's are slightly jumbled, they are only jumbled within their symmetry blocks.  
! They still correspond to the symmetries of SymLabelList2_rot, which is the important part.

! But in order to look at the output, it is easier to consider them in terms of highest 
! occupied to lowest occupied - i.e. in terms of the NO eigenvalues (occupation numbers).
        call OrderNatOrbMat()

    END SUBROUTINE DiagRDM


    SUBROUTINE OrderNatOrbMat()
        INTEGER :: spin,i,j,ierr,StartSort,EndSort
        CHARACTER(len=*), PARAMETER :: this_routine='OrderRDM'
        INTEGER , ALLOCATABLE :: SymLabelList3_rot(:)
        REAL(dp) , ALLOCATABLE :: NatOrbMatTemp(:,:), EvaluesTemp(:)
        INTEGER :: NatOrbMatTempTag, SymLabelList3_rotTag, EvaluesTempTag, Orb, New_Pos
        
! Here, if symmetry is kept, we are going to have to reorder the eigenvectors 
! according to the size of the eigenvalues, while taking the orbital labels 
! (and therefore symmetries) with them. This will be put back into MP2VDM from MP2VDMTemp.

! Want to reorder the eigenvalues from largest to smallest, taking the eigenvectors 
! with them and the symmetry as well.  
! If using spin orbitals, do this for the alpha spin and then the beta.

        ALLOCATE(NatOrbMatTemp(NoOrbs,NoOrbs),stat=ierr)
        CALL LogMemAlloc('NatOrbMatTemp',NoOrbs**2,8,&
                            'OrderNatOrbMat',NatOrbMatTempTag,ierr)
        ALLOCATE(SymLabelList3_rot(NoOrbs),stat=ierr)
        CALL LogMemAlloc('SymLabelList3_rot',NoOrbs,4,&
                            'OrderNatOrbMat',SymLabelList3_rotTag,ierr)
        ALLOCATE(EvaluesTemp(NoOrbs),stat=ierr)
        CALL LogMemAlloc('EvaluesTemp',NoOrbs,4,&
                            'OrderNatOrbMat',EvaluesTempTag,ierr)

! Want to remember the original orbital ordering, as after the sort, the MO's 
! will still have this ordering.
        SymLabelList3_rot(:) = SymLabelList2_rot(:)

        StartSort=1
        EndSort=SpatOrbs

! Unfortunately this sort routine orders the orbitals in ascending order... which is not quite 
! what we want.  Just remember this when printing out the Evalues.
        call sort (EValues(startSort:endSort), &
                   NatOrbMat(1:NoOrbs, startSort:endSort), &
                   SymLabelList2_rot(startSort:endSort))

        if(tStoreSpinOrbs) then                  
            StartSort=SpatOrbs + 1
            EndSort=nBasis

            call sort (EValues(startSort:endSort), &
                       NatOrbMat(1:NoOrbs, startSort:endSort), &
                       SymLabelList2_rot(startSort:endSort))

        endif                       

! We now have the NO's ordered according to the size of their Evalues (occupation 
! numbers).  This will have jumbled up their symmetries.  Want to reorder the 
! MO's to match this ordering (so that we only have one SymLabelList array).

! Need a new SymLabelListInv_rot too.        
        SymLabelListInv_rot(:)=0   
        do i=1,NoOrbs
            SymLabelListInv_rot(SymLabelList2_rot(NoOrbs-i+1))=i
        enddo

        NatOrbMatTemp(:,:) = NatOrbMat(:,:)
        NatOrbMat(:,:) = 0.0_dp

        do i=1,NoOrbs
            do j = 1, NoOrbs

                ! In position j, the MO orbital Orb is currently there.
                Orb = SymLabelList3_rot(j)      

                ! Want to move it to the position the NO's are in.
                New_Pos = SymLabelListInv_rot(Orb)   

                ! But we also want to reverse the order of everything... 
                NatOrbMat(New_Pos,NoOrbs - i + 1)=NatOrbMatTemp(j,i)
            enddo
        enddo

        SymLabelList3_rot(:) = SymLabelList2_rot(:)
        EvaluesTemp(:) = Evalues(:)
        do i = 1, NoOrbs
            SymLabelList2_rot(i) = SymLabelList3_rot(NoOrbs - i + 1)
            Evalues(i) = EvaluesTemp(NoOrbs - i + 1)
        enddo

        DEALLOCATE(NatOrbMatTemp)
        DEALLOCATE(SymLabelList3_rot)
        DEALLOCATE(EvaluesTemp)

    END SUBROUTINE OrderNatOrbMat

   
! This is an M^5 transform, which transforms all the two-electron integrals 
! into the new basis described by the Coeff matrix.
! This is v memory inefficient and currently does not use any spatial 
! symmetry information.
    SUBROUTINE Transform2ElIntsMemSave_RDM()
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
 
        FourIndInts(:,:,:,:)=0.0_dp

! **************
! Calculating the two-transformed, four index integrals.

! The untransformed <alpha beta | gamma delta> integrals are found from 
! UMAT(UMatInd(i,j,k,l,0,0)
! All our arrays are in spin orbitals - if tStoreSpinOrbs is false, 
! UMAT will be in spatial orbitals - need to account for this.

! Running through 1,NoOrbs - the actual orbitals corresponding to that index 
! are given by SymLabelList2_rot

        do b=1,NoOrbs
            b2 = SymLabelList2_rot(b)
            do d=1,NoOrbs
                d2 = SymLabelList2_rot(d)
                do a=1,NoOrbs
                    a2 = SymLabelList2_rot(a)
                    do g=1,NoOrbs
                        g2 = SymLabelList2_rot(g)

! UMatInd in physical notation, but FourIndInts in chemical 
! (just to make it more clear in these transformations).
! This means that here, a and g are interchangable, and so are b and d.
                        FourIndInts(a,g,b,d)=REAL(UMAT(UMatInd(a2,b2,g2,d2,0,0)),dp)
                    enddo
                enddo

                Temp4indints(:,:)=0.0_dp
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

                Temp4indints(:,:)=0.0_dp
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
        ArrDiagNew(:)=0.0_dp                     

!        WRITE(6,*) 'The diagonal fock elements in the HF basis set'
!        do a=1,nBasis
!            WRITE(6,'(F20.10)',advance='no') Arr(a,2)
!            WRITE(6,*) Arr(:,2)
!        enddo

! First calculate the sum of the diagonal elements, ARR.
! Check if this is already being done.
        FOCKDiagSumHF=0.0_dp
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

        FOCKDiagSumNew=0.0_dp
        do j=1,NoOrbs
            l=SymLabelList2_rot(j)
            if(tStoreSpinOrbs) then
                ArrDiagNew(l) = 0.0_dp
            else
                ArrDiagNew(2*l) = 0.0_dp
                ArrDiagNew((2*l)-1) = 0.0_dp
            endif
            do a=1,NoOrbs
                b=SymLabelList2_rot(a)
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
!        CALL neci_flush(6)
!        stop       

        DEALLOCATE(ArrDiagNew)
        CALL LogMemDealloc(this_routine,ArrDiagNewTag)

    ENDSUBROUTINE CalcFOCKMatrix_RDM

    SUBROUTINE RefillUMATandTMAT2D_RDM()
! UMat is in spin or spatial orbitals, TMAT2D only spin.
! This routine refills these to more easily write out the ROFCIDUMP, and originally 
! to be able to continue a calculation (although I doubt this works at the moment).
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
        TMAT2DPart(:,:)=0.0_dp

! Make the UMAT elements the four index integrals.  
! These are calculated by transforming the HF orbitals using the coefficients 
! that have been found
        do l=1,NoOrbs
            d=SymLabelList2_rot(l)
            do k=1,NoOrbs
                b=SymLabelList2_rot(k)
                do j=1,NoOrbs
                    g=SymLabelList2_rot(j)
                    do i=1,NoOrbs
                        a=SymLabelList2_rot(i)
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
                i=SymLabelList2_rot(k)
                NewTMAT=0.0_dp
                do b=1,NoOrbs
                    d=SymLabelList2_rot(b)
                    if(tStoreSpinOrbs) then
                        NewTMAT=NewTMAT+(NatOrbMat(b,k)*REAL(TMAT2D(d,a),dp))
                    else
                        NewTMAT=NewTMAT+(NatOrbMat(b,k)*REAL(TMAT2D(2*d,a),dp))
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
                j=SymLabelList2_rot(l)
                NewTMAT=0.0_dp
                do a=1,NoOrbs
                    c=SymLabelList2_rot(a)
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
!        CALL neci_flush(6)
!        stop

        DEALLOCATE(TMAT2DPart)
        CALL LogMemDeAlloc('RefillUMAT_RDM',TMAT2DPartTag)

        if (trotatedNOs) then
            call PrintROFCIDUMP_RDM("BSFCIDUMP")
        else
            CALL PrintROFCIDUMP_RDM("ROFCIDUMP")
        endif

    ENDSUBROUTINE RefillUMATandTMAT2D_RDM


    SUBROUTINE PrintROFCIDUMP_RDM(filename)
!This prints out a new FCIDUMP file in the same format as the old one.
        INTEGER :: i,j,k,l,iunit
        character(len=9) :: filename

!        PrintROFCIDUMP_Time%timer_name='PrintROFCIDUMP'
!        CALL set_timer(PrintROFCIDUMP_Time,30)

        iunit = get_free_unit()
        OPEN(iunit,FILE=filename,STATUS='unknown')!'ROFCIDUMP',STATUS='unknown')
        
        WRITE(iunit,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',NoOrbs,&
                                                ',NELEC=',NEl,',MS2=',LMS,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,NoOrbs
            IF(tStoreSpinOrbs) THEN
                WRITE(iunit,'(I1,A1)',advance='no') INT(G1(i)%sym%S)+1,','
            ELSE
                if (tRotatedNOs.and.tBrokenSymNOs) then
                    write(iunit,'(I1,A1)',advance='no') 1,','
                else
                    WRITE(iunit,'(I1,A1)',advance='no') INT(G1(i*2)%sym%S)+1,','
                endif
            ENDIF
        enddo
        WRITE(iunit,*) ''
        IF(tStoreSpinOrbs) THEN
            WRITE(iunit,'(A7,I1,A11)') 'ISYM=',1,' UHF=.TRUE.'
        ELSE
            WRITE(iunit,'(A7,I1,A12)') 'ISYM=',1,' UHF=.FALSE.'
        ENDIF
        if(tFixLz) then
            WRITE(iunit,'(A7)',advance='no') 'SYML='
            do i=1,NoOrbs
                if(i.eq.NoOrbs) then
                    WRITE(iunit,'(I3,A1)') -20,','
                else
                    WRITE(iunit,'(I3,A1)',advance='no') -20,','
                endif
            enddo
            WRITE(iunit,'(A8)',advance='no') 'SYMLZ='
            do i=1,NoOrbs
                if(i.eq.NoOrbs) then
                    IF(tStoreSpinOrbs) THEN
                        WRITE(iunit,'(I2,A1)') INT(G1(i)%Ml),','
                    ELSE
                        WRITE(iunit,'(I2,A1)') INT(G1(i*2)%Ml),','
                    ENDIF
                else
                    IF(tStoreSpinOrbs) THEN
                        WRITE(iunit,'(I2,A1)',advance='no') INT(G1(i)%Ml),','
                    ELSE
                        WRITE(iunit,'(I2,A1)',advance='no') INT(G1(i*2)%Ml),','
                    ENDIF
                endif
            enddo
        endif
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
                        IF((ABS(REAL(UMat(UMatInd(i,j,k,l,0,0)),dp))).ne.0.0_dp) &
                                WRITE(iunit,'(F21.12,4I3)') &
                                REAL(UMat(UMatInd(i,j,k,l,0,0)),dp),i,k,j,l 
 
                    enddo
                enddo
           enddo
        enddo

! TMAT2D stored as spin orbitals
        do i=1,NoOrbs
            ! Symmetry?
            do k=1,NoOrbs
                IF(tStoreSpinOrbs) THEN
                    IF((REAL(TMAT2D(i,k),dp)).ne.0.0_dp) WRITE(iunit,'(F21.12,4I3)') &
                                                        REAL(TMAT2D(i,k),dp),i,k,0,0
                ELSE
                    IF((REAL(TMAT2D(2*i,2*k),dp)).ne.0.0_dp) WRITE(iunit,'(F21.12,4I3)') &
                                                        REAL(TMAT2D(2*i,2*k),dp),i,k,0,0
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
        
        CALL neci_flush(iunit)

        CLOSE(iunit)

!        CALL halt_timer(PrintROFCIDUMP_Time)

    ENDSUBROUTINE PrintROFCIDUMP_RDM

    subroutine BrokenSymNO(occ_numb_diff)

        ! We'd like to keep both the original NOs (Natural Orbitals) and the broken maximally
        ! localised NOs for checking
        real(dp), intent(in) :: occ_numb_diff
        real(dp) :: diffnorm,SumDiag,sum_old,sum_new
        real(dp), allocatable :: no_store(:,:)
        real(dp), allocatable :: trans_2orbs_coeffs(:,:)
        real(dp), allocatable :: selfint(:)
        integer, allocatable :: rotate_list(:,:)!,orbs(:,:)
        integer :: l1,l2,l3,l4,m,n
        logical :: partnerfound

        allocate(trans_2orbs_coeffs(2,2))
        !allocate(orbs(6,2))
        allocate(rotate_list((NoOrbs*(NoOrbs-1)),4))
        allocate(selfint(NoOrbs))
        allocate(no_store(NoOrbs,NoOrbs))

        ! need to store NO coefficients since these are overwritten during 
        ! the orbital rotation
        no_store = NatOrbMat

        if (iProcIndex.eq.0) then

            ! Normalisation
            SumDiag = sum(Evalues)

            if (tStoreSpinOrbs) then
                diffnorm = SumDiag/dble(NEl)
            else
                diffnorm = 2.0_dp*(SumDiag/dble(NEl))
            endif
            
            trotatedNOs=.false.
     
            if (tStoreSpinorbs) then
                call Stop_all("BrokenSymNO","Broken symmetry NOs currently not implemented for UHF")
            endif

            write(6,*) '------------------------------------------------------------------------------'
            write(6,*) 'Localising NOs whose occupation numbers differ by less than&
                & threshold'
            write(6,*) '------------------------------------------------------------------------------'
            write(6,*) 'Threshold for orbitals to rotate:',occ_numb_diff

            ! self-interactions
            selfint(:) = 0.0_dp
            do l1=1,NoOrbs
                selfint(l1) = Umat(UmatInd(l1,l1,l1,l1,0,0))
            enddo

            write(6,*) 'Self-interactions for NOs:'
            do l1=1,NoOrbs
                write(6,'(I3,3X,G25.12)') l1, selfint(l1)
            enddo
            write(6,*) 'Sum of NO selfinteractions:',sum(selfint)

            NatOrbMat = 0.0_dp
            ! Generate the list of orbitals which are rotated amongst each
            ! other
            rotate_list(:,:) = 0
            ! Need to account for spatial and spin orbital representations
            ! since orbitals of different spin cannot be mixed
            ! List contains the NOs which are rotated
            ! It can deal with a maximum of four NOs which are mixed
            m = 0
            n = 1
            do l1=1,NoOrbs
                if ((m.ne.0).and.(l1.le.rotate_list(m,n))) cycle
                    partnerfound = .false.
                    n = 1
                    do l2=(l1+1),NoOrbs
                    !write(6,*) (dabs((Evalues(l1)/diffnorm))-(Evalues(l2)/diffnorm))/dabs((Evalues(l2)/diffnorm))
                        if ((abs((Evalues(l1)/diffnorm)-(Evalues(l2)/diffnorm))/abs((Evalues(l2)/diffnorm)))&
                            &.lt.occ_numb_diff) then
                            if (.not.partnerfound) then
                                m = m + 1
                                n = n + 1
                                rotate_list(m,1) = l1
                                rotate_list(m,2) = l2
                                partnerfound = .true.
                            elseif (partnerfound) then
                                n = n + 1
                                if (n.gt.2) then
                                    n = 2
                                    write(6,*) '***Warning***'
                                    write(6,*) 'Threshold generated more than 2-fold degeneracy'
                                    write(6,*) 'NOs around:',l2
                                    cycle  ! don't want to rotate more than 2 orbitals
                                endif
                                ! this is for up to four-fold degeneracy
                             !       if (n.gt.4) then
                             !           n = 4
                             !           write(6,*) '***Warning***'
                             !           write(6,*) 'Threshold generated more than 4-fold degeneracy'
                             !           write(6,*) 'NOs around:',l2
                             !           cycle  ! don't want to rotate more than 4 orbitals
                             !       endif
                                rotate_list(m,n) = l2
                            endif
                        endif
                    enddo
                if (.not.partnerfound) then
                    ! don't rotate orbital
                    NatOrbMat(l1,l1) = 1.0_dp
                endif
            enddo

            write(6,*) 'The following pairs of orbitals will be rotated:'
            do l1=1,m
                write(6,'(I3,3X,4(I3))') l1,rotate_list(l1,:)
            enddo
            ! rotate two-fold degeneracate pairs first
            do l1=1,m
                ! if only two orbitals have same occupation numbers
                if (rotate_list(l1,3).eq.0) then
                    write(6,'(A20,4(I3))') 'Rotating NOs:',rotate_list(l1,:)
                    call Rotate2Orbs(rotate_list(l1,1),rotate_list(l1,2),trans_2orbs_coeffs,selfint(rotate_list(l1,1)),&
                        &selfint(rotate_list(l1,2)))
                    ! The new NOs are 
                    ! phi_{i'} = cos a p_{i} + sin a p_{j}
                    ! phi_{j'} = -sin a p_{i} + cos a p_{j}
                    NatorbMat(rotate_list(l1,1),rotate_list(l1,1)) = trans_2orbs_coeffs(1,1)
                    NatorbMat(rotate_list(l1,2),rotate_list(l1,1)) = trans_2orbs_coeffs(2,1)
                    NatorbMat(rotate_list(l1,1),rotate_list(l1,2)) = trans_2orbs_coeffs(1,2)
                    NatorbMat(rotate_list(l1,2),rotate_list(l1,2)) = trans_2orbs_coeffs(2,2)
                    write(6,*) 'Sum of rotated NO self-interactions:',sum(selfint)
                endif
           enddo
           ! this is for dealing with up to four-fold degeneracy
           ! it needs some further work to work pro[erly
     !      do l1=1,m
     !            if ((rotate_list(l1,3).ne.0).and.(rotate_list(l1,4).eq.0)) then
     !                   sum_new = sum(selfint)
     !                   orbs(1,1) = 1
     !                   orbs(1,2) = 2
     !                   orbs(2,1) = 1
     !                   orbs(2,2) = 3
     !                   orbs(3,1) = 2
     !                   orbs(3,2) = 3
     !       !            do
     !                       sum_old = sum_new
     !                       write(6,'(A20,4(I3))') 'Rotating NOs:',rotate_list(l1,:)
     !                       do l3=1,1!3
     !                           NatOrbMat(:,:) = 0.0_dp
     !                           do l4=1,NoOrbs
     !                               NatOrbMat(l1,l1) = 1.0_dp
     !                           enddo
     !                           NatOrbMat(rotate_list(l1,orbs(l3,1)),rotate_list(l1,orbs(l3,1))) = 0.0_dp
     !                           NatOrbMat(rotate_list(l1,orbs(l3,2)),rotate_list(l1,orbs(l3,2))) = 0.0_dp
     !                           call Rotate2Orbs(rotate_list(l1,orbs(l3,1)),rotate_list(l1,orbs(l3,2)),&
     !                               &trans_2orbs_coeffs,selfint(rotate_list(l1,orbs(l3,1))),&
     !                               &selfint(rotate_list(l1,orbs(l3,2))))
     !                           ! The new NOs are 
     !                           ! phi_{i'} = cos a p_{i} + sin a p_{j}
     !                           ! phi_{j'} = -sin a p_{i} + cos a p_{j}
     !                           NatorbMat(rotate_list(l1,orbs(l3,1)),rotate_list(l1,orbs(l3,1)))&
     !                               & = trans_2orbs_coeffs(1,1)
     !                           NatorbMat(rotate_list(l1,orbs(l3,2)),rotate_list(l1,orbs(l3,1)))&
     !                               & = trans_2orbs_coeffs(2,1)
     !                           NatorbMat(rotate_list(l1,orbs(l3,1)),rotate_list(l1,orbs(l3,2)))&
     !                               & = trans_2orbs_coeffs(1,2)
     !                           NatorbMat(rotate_list(l1,orbs(l3,2)),rotate_list(l1,orbs(l3,2)))&
     !                               & = trans_2orbs_coeffs(2,2)
     !                       enddo
     !        !               sum_new = sum(selfint)
     !                       !write(6,*) sum_new,sum_old
     !        !               if (abs(sum_new).le.abs(sum_old)) then
     !        !                   exit
     !        !               endif
     !        !           enddo
     !               !endif
     !               write(6,*) 'Sum of rotated NO self-interactions:',sum(selfint)
     !           endif
     !      enddo
     !      do l1=1,m
     !           ! if four orbitals have the same occupation numbers
     !           if (rotate_list(l1,4).ne.0) then
     !               sum_new = sum(selfint)
     !               orbs(1,1) = 1
     !               orbs(1,2) = 2
     !               orbs(2,1) = 3
     !               orbs(2,2) = 4
     !               orbs(3,1) = 1
     !               orbs(3,2) = 3
     !               orbs(4,1) = 2
     !               orbs(4,2) = 4
     !               orbs(5,1) = 1
     !               orbs(5,2) = 4
     !               orbs(6,1) = 2
     !               orbs(6,2) = 3
     !         !      do
     !                   sum_old = sum_new
     !                   write(6,'(A20,4(I3))') 'Rotating NOs:',rotate_list(l1,:)
     !                   do l3=1,1!3
     !                       NatOrbMat(:,:) = 0.0_dp
     !                       do l4=1,NoOrbs
     !                           NatOrbMat(l1,l1) = 1.0_dp
     !                       enddo
     !                       NatOrbMat(rotate_list(l1,orbs(((2*l3)-1),1)),&
     !                           &rotate_list(l1,orbs(((2*l3)-1),1))) = 0.0_dp
     !                       NatOrbMat(rotate_list(l1,orbs(((2*l3)-1),2)),&
     !                           &rotate_list(l1,orbs(((2*l3)-1),2))) = 0.0_dp
     !                       NatOrbMat(rotate_list(l1,orbs(((2*l3)),1)),&
     !                           &rotate_list(l1,orbs(((2*l3)),1))) = 0.0_dp
     !                       NatOrbMat(rotate_list(l1,orbs(((2*l3)),2)),&
     !                           &rotate_list(l1,orbs(((2*l3)),2))) = 0.0_dp
     !                       call Rotate2Orbs(rotate_list(l1,orbs(((2*l3)-1),1)),&
     !                           &rotate_list(l1,orbs(((2*l3)-1),2)),&
     !                           &trans_2orbs_coeffs,selfint(rotate_list(l1,orbs(((2*l3)-1),1))),&
     !                           &selfint(rotate_list(l1,orbs(((2*l3)-1),2))))
     !                           ! The new NOs are 
     !                           ! phi_{i'} = cos a p_{i} + sin a p_{j}
     !                           ! phi_{j'} = -sin a p_{i} + cos a p_{j}
     !                           NatorbMat(rotate_list(l1,orbs(((2*l3)-1),1)),&
     !                               &rotate_list(l1,orbs(((2*l3)-1),1))) = trans_2orbs_coeffs(1,1)
     !                           NatorbMat(rotate_list(l1,orbs(((2*l3)-1),2)),&
     !                               &rotate_list(l1,orbs(((2*l3)-1),1))) = trans_2orbs_coeffs(2,1)
     !                           NatorbMat(rotate_list(l1,orbs(((2*l3)-1),1)),&
     !                               &rotate_list(l1,orbs(((2*l3)-1),2))) = trans_2orbs_coeffs(1,2)
     !                           NatorbMat(rotate_list(l1,orbs(((2*l3)-1),2)),&
     !                               &rotate_list(l1,orbs(((2*l3)-1),2))) = trans_2orbs_coeffs(2,2)
     !                           call Rotate2Orbs(rotate_list(l1,orbs(((2*l3)),1)),&
     !                               &rotate_list(l1,orbs(((2*l3)),2)),&
     !                               &trans_2orbs_coeffs,selfint(rotate_list(l1,orbs(((2*l3)),1))),&
     !                               &selfint(rotate_list(l1,orbs(((2*l3)),2))))
     !                           NatorbMat(rotate_list(l1,orbs(((2*l3)),1)),&
     !                               &rotate_list(l1,orbs(((2*l3)),1))) = trans_2orbs_coeffs(1,1)
     !                           NatorbMat(rotate_list(l1,orbs(((2*l3)),2)),&
     !                               &rotate_list(l1,orbs(((2*l3)),1))) = trans_2orbs_coeffs(2,1)
     !                           NatorbMat(rotate_list(l1,orbs(((2*l3)),1)),&
     !                               &rotate_list(l1,orbs(((2*l3)),2))) = trans_2orbs_coeffs(1,2)
     !                           NatorbMat(rotate_list(l1,orbs(((2*l3)),2)),&
     !                               &rotate_list(l1,orbs(((2*l3)),2))) = trans_2orbs_coeffs(2,2)
     !                       enddo
     !             !          sum_new = sum(selfint)
     !             !          !write(6,*) sum_new,sum_old
     !             !          if (abs(sum_new).le.abs(sum_old)) then
     !             !              exit
     !             !          endif
     !             !      enddo
     !              ! endif
     !                   write(6,*) 'Sum of rotated NO self-interactions:',sum(selfint)
     !               endif
     !           enddo

            write(6,*) 'Final self-interactions for rotated NOs:'
            do l1=1,NoOrbs
                write(6,'(I3,3X,G25.12)') l1, selfint(l1)
            enddo
            write(6,*) 'Sum of rotated NO self-interactions:',sum(selfint)

            !if (tstillrotating) then
            write(6,*) '------------------------------------------------------'
            write(6,*) 'Writing out BSFCIDUMP...'
            write(6,*) '------------------------------------------------------'
            tRotatedNOs=.true.
            call Transform2ElIntsMemSave_RDM()
            call RefillUMATandTMAT2D_RDM()
            !call PrintROFCIDUMP_RDM("BSFCIDUMP")

 
        endif
        
        ! restore original NOs
        NatOrbMat(:,:) = 0.0_dp
        NatOrbMat = no_store

        deallocate(trans_2orbs_coeffs)
        deallocate(rotate_list)
        !deallocate(orbs)
        deallocate(selfint)
        deallocate(no_store)


    endsubroutine BrokenSymNO

    subroutine Rotate2Orbs(i,j,trans_2orbs_coeffs,selfintorb1,selfintorb2)

        ! This routine takes two orbitals i,j, and rotates them in order to maximally localise these
        ! It employs an Edminston-Ruedenberg type localisation which maximises
        ! \sum_{i=1}^{2} \sum_{r,s,u,v} (c_{ir})*(c_{is})*c_{iu}c_{iv} <p_{i}p_{i}|u|p_{i}p_{i}>
        ! where p_{i} are the original NOs
        ! The the coefficients c are given by the following matrix:
        ! c =  cos a   sin a
        !      -sin a  cos a
        ! Then angle a is found by differentiating and setting the sum equal to 0 which gives
        ! the following analytical expression (notation <ii|u|ii> =  <p_{i}p_{i}|u|p_{i}p_{i}>)
        ! tan a = -x/y
        ! where x and y are sums of the original NO four index inegrals
        real(dp), allocatable, intent(inout) :: trans_2orbs_coeffs(:,:)
        real(dp), intent(inout) :: selfintorb1,selfintorb2
        real(dp) :: alpha(2)
        !real(dp) :: secondderiv(2)
        real(dp) :: selfinteractions(2)
        real(dp) :: coeffcos,coeffsin
        integer :: indicesij(2)
        integer, intent(in) :: i,j
        integer :: l1,l2,l3,l4,l5

        ! Umat(UMatInd(i,j,k,l,0,0)) contains the four-index integrals
        ! <ij|kl> (physical notation) in the NO basis

        indicesij(1) = i
        indicesij(2) = j
        trans_2orbs_coeffs(:,:) = 0.0_dp

        coeffcos = Umat(UmatInd(i,i,i,j,0,0)) + Umat(UmatInd(i,i,j,i,0,0)) + Umat(UmatInd(i,j,i,i,0,0)) &
            & - Umat(UmatInd(i,j,j,j,0,0)) + Umat(UmatInd(j,i,i,i,0,0)) - Umat(UmatInd(j,i,j,j,0,0)) &
            & - Umat(UmatInd(j,j,i,j,0,0)) - Umat(UmatInd(j,j,j,i,0,0))

        coeffsin = -Umat(UmatInd(i,i,i,i,0,0)) + Umat(UmatInd(i,i,j,j,0,0)) + Umat(UmatInd(i,j,i,j,0,0)) &
            & + Umat(UmatInd(i,j,j,i,0,0)) + Umat(UmatInd(j,i,i,j,0,0)) + Umat(UmatInd(j,i,j,i,0,0)) &
            & + Umat(UmatInd(j,j,i,i,0,0)) - Umat(UmatInd(j,j,j,j,0,0))

        ! atan return a value in [-pi/2,pi/2]
        ! there are two values of alpha in the range 0,2pi
        ! which satisfy the equation
        alpha(1) = atan((-coeffsin/coeffcos))
        ! if alpha < 0 shift it in the positive range
        if (alpha(1).lt.0.0_dp) then
            alpha(1) = alpha(1) + pi
        endif
        !q the second value is alpha+pi
        alpha(2) = alpha(1) + pi
        
        alpha = alpha/4.0_dp

        !! second derivatives to find maximum (necessary since the minimum, i.e. fully delocalised
        !! orbitals satisfy the same conditions
        !secondderiv(1) = (4.0_dp*coeffsin*cos(4.0_dp*alpha(1))) - (4.0_dp*coeffcos*sin(4.0_dp*alpha(1)))
        !secondderiv(2) = (4.0_dp*coeffsin*cos(4.0_dp*alpha(2))) - (4.0_dp*coeffcos*sin(4.0_dp*alpha(2)))

        ! compute selfinteractions to check which one is largest
        ! this is a better measure than the second derivatives
        selfinteractions(:) = 0.0_dp 

        do l1=1,2
            trans_2orbs_coeffs(1,1) = cos(alpha(l1))
            trans_2orbs_coeffs(2,1) = sin(alpha(l1))
            trans_2orbs_coeffs(1,2) = -sin(alpha(l1))
            trans_2orbs_coeffs(2,2) = cos(alpha(l1))

            do l2=1,2
                do l3=1,2
                    do l4=1,2
                        do l5=1,2
                            selfinteractions(l1) = selfinteractions(l1) + trans_2orbs_coeffs(l2,1)&
                                &*trans_2orbs_coeffs(l3,1)*&
                                &trans_2orbs_coeffs(l4,1)*trans_2orbs_coeffs(l5,1)*&
                                &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5),0,0))
                        enddo
                    enddo
                enddo
            enddo
            do l2=1,2
                do l3=1,2
                    do l4=1,2
                        do l5=1,2
                            selfinteractions(l1) = selfinteractions(l1) + trans_2orbs_coeffs(l2,2)&
                                &*trans_2orbs_coeffs(l3,2)*&
                                &trans_2orbs_coeffs(l4,2)*trans_2orbs_coeffs(l5,2)*&
                                &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5),0,0))
                        enddo
                    enddo
                enddo
            enddo
        enddo

 
        if (selfinteractions(1).ge.selfinteractions(2)) then
            trans_2orbs_coeffs(1,1) = cos(alpha(1))
            trans_2orbs_coeffs(2,1) = sin(alpha(1))
            trans_2orbs_coeffs(1,2) = -sin(alpha(1))
            trans_2orbs_coeffs(2,2) = cos(alpha(1))
        elseif (selfinteractions(2).gt.selfinteractions(1)) then
            trans_2orbs_coeffs(1,1) = cos(alpha(2))
            trans_2orbs_coeffs(2,1) = sin(alpha(2))
            trans_2orbs_coeffs(1,2) = -sin(alpha(2))
            trans_2orbs_coeffs(2,2) = cos(alpha(2))
        endif

        ! new sefl-interactions
        selfintorb1 = 0.0_dp
        selfintorb2 = 0.0_dp
        do l2=1,2
            do l3=1,2
                do l4=1,2
                    do l5=1,2
                        selfintorb1 = selfintorb1 + trans_2orbs_coeffs(l2,1)&
                            &*trans_2orbs_coeffs(l3,1)*&
                            &trans_2orbs_coeffs(l4,1)*trans_2orbs_coeffs(l5,1)*&
                            &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5),0,0))
                        selfintorb2 = selfintorb2 + trans_2orbs_coeffs(l2,2)&
                            &*trans_2orbs_coeffs(l3,2)*&
                            &trans_2orbs_coeffs(l4,2)*trans_2orbs_coeffs(l5,2)*&
                            &Umat(UmatInd(indicesij(l2),indicesij(l3),indicesij(l4),indicesij(l5),0,0))
                    enddo
                enddo
            enddo
        enddo


    endsubroutine Rotate2Orbs

    SUBROUTINE DeallocateRDM()
! This routine just deallocates the arrays allocated in InitRDM.
! If the NECI calculation softexits before the RDMs start to fill, this is all that 
! is called at the end.
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

            if(allocated(NatOrbMat)) then
                DEALLOCATE(NatOrbMat)
                CALL LogMemDeAlloc(this_routine,NatOrbMatTag)
            endif

            IF((iProcIndex.eq.0).and.tDiagRDM) THEN
                if(allocated(Evalues)) then
                    DEALLOCATE(Evalues)
                    CALL LogMemDeAlloc(this_routine,EvaluesTag)
                endif

                if(allocated(Rho_ii)) then
                    DEALLOCATE(Rho_ii)
                    CALL LogMemDeAlloc(this_routine,Rho_iiTag)
                endif

                IF(tPrintRODump) THEN
                    if(allocated(FourIndInts)) then
                        DEALLOCATE(FourIndInts)
                        CALL LogMemDeAlloc(this_routine,FourIndIntsTag)
                    endif
                ENDIF

            ENDIF

            DEALLOCATE(SymLabelCounts2_rot)
            CALL LogMemDeAlloc(this_routine,SymLabelCounts2_rotTag)

            DEALLOCATE(SymLabelList2_rot)
            CALL LogMemDeAlloc(this_routine,SymLabelList2_rotTag)

            DEALLOCATE(SymLabelListInv_rot)
            CALL LogMemDeAlloc(this_routine,SymLabelListInv_rotTag)

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
                    if(allocated(NatOrbMat)) then
                        DEALLOCATE(NatOrbMat)
                        CALL LogMemDeAlloc(this_routine,NatOrbMatTag)
                    endif
                endif
            ENDIF

        ENDIF

    END SUBROUTINE DeallocateRDM

    subroutine convert_mats_Molpforces(Norm_1RDM, Norm_2RDM)
    use SystemData, only : nEl,nbasis,LMS,ECore
    use IntegralsData, only : nFrozen
    use NatOrbsMod, only : NatOrbMat
    use SymData, only : Sym_Psi, SymLabelCounts,nSymLabels
    use SymExcitDataMod, only: SpinOrbSymLabel,SymLabelCounts2
    use GenRandSymExcitNUMod , only : ClassCountInd, RandExcitSymLabelProd
    use sym_mod
    use UMatCache, only: GTID
    
    integer :: iblkq, iseccr, istat1, isyref, ms2
    integer :: posn1, posn2
    integer :: i, j, k, l
    integer :: myname, ifil, intrel, iout
    integer :: orb1, orb2, Sym_i, Sym_j, Sym_ij
    integer :: Sym_k, Sym_l, Sym_kl
    integer, dimension(8) :: iact, ldact !iact(:) # of active orbs per sym, 
                                         !ldact(:) - # Pairs of orbs that multiply to give given sym
    integer, dimension(8) :: icore, iclos 
    integer, dimension(nSymLabels):: blockstart1, blockstart2
    integer, dimension(nSymLabels) :: elements_assigned1, elements_assigned2
    integer :: FC_Lag_Len  !Length of the Frozen Core Lagrangian
    integer :: Len_1RDM, Len_2RDM, FCLag_Len
    real(dp) :: ijkl, jikl
    real(dp), intent(in) :: Norm_1RDM, Norm_2RDM
    real(dp), allocatable :: SymmetryPacked2RDM(:), SymmetryPacked1RDM(:) !SymPacked2RDM in CHEMICAL Notation
    real(dp), allocatable :: SymmetryPackedLagrangian(:)
    real(dp), allocatable :: FC_Lagrangian(:)
    integer :: iWfRecord, iWfSym
    real(dp) :: WfRecord
    character(len=*), parameter :: t_r='convert_mats_Molpforces'

#ifdef __INT64
    intrel=1
#else
    intrel=2
#endif

    if (iProcIndex .eq. 0) then
        !Calculating header information needed for the molpro dump routine
        !Considering all electron for now
        icore(:)=0  !Number of frozen orbitals per symmetry
        iclos(:)=0  !icore(i) + Number of 'closed' orbitals per symmetry (iclos is the same as icore for us)
!        iblkq=21402
        call molpro_get_reference_info(iWfRecord, iWfSym)
        iblkq = iWfRecord ! record number for orbitals
        iseccr=0          ! record number for core orbitals
        istat1=1        !ground state
        isyref=Sym_Psi+1  !spatial symmetry of the wavefunction
#ifdef MOLPRO
        if (isyref.ne.iWfSym) call stop_all(t_r,"NECI and common/cref do not agree on irrep of wave function")
#endif
        ms2=LMS  !2 * M_s
        myname=5001 !Arbitrary file names
        ifil=1
        ! ^- cgk: might want to use nexfre(igrsav) for these two.
        iout=molpro_get_iout()
        ! ^- thise one should come from common/tapes. 
        ldact(:)=0
        iact(:)=0
        Len_1RDM=0
        Len_2RDM=0
        blockstart1(:)=0
        blockstart2(:)=0
        elements_assigned1(:)=0
        elements_assigned2(:)=0
        
        
        ! Find out the number of orbital pairs that multiply to a given sym (ldact)
        do i=1,SpatOrbs  
            !run over spatial orbitals
            do j=1,  i    ! i .ge. j
                Sym_i=SpinOrbSymLabel(2*i)  !Consider only alpha orbitals
                Sym_j=SpinOrbSymLabel(2*j)
                Sym_ij=RandExcitSymLabelProd(Sym_i, Sym_j)
                ldact(Sym_ij+1)=ldact(Sym_ij+1)+1
            enddo
        enddo

        !Calculate lengths of arrays, and where each sym block starts
        do i=0,nSymLabels-1

            !CMO: Check if Sym_i goes 0-->7 or 1-->8
            !Find position of each symmetry block in sym-packed forms of RDMS 1 & 2
            blockstart1(i+1)=Len_1RDM+1 !N.B. Len_1RDM still being updated in this loop
            blockstart2(i+1)=Len_2RDM+1 !N.B. Len_2RDM still being updated in this loop
            
            ! Count the number of active orbitals of the given symmetry
            iact(i+1)=SymLabelCounts2(2,ClassCountInd(1,i,0)) !Count the number of active orbitals of the given symmetry
            Len_1RDM=Len_1RDM+(iact(i+1)*(iact(i+1)+1)/2) ! add on # entries in sym-packed 1RDM for sym i
            
            Len_2RDM=Len_2RDM+(ldact(i+1))**2 !Assumes no frozen orbitals
        enddo

        FCLag_Len=SpatOrbs**2  !! Arbitrarily set this for now - we will not be printing it whilst nfrozen=0
        
        !Allocate arrays accordingly
        allocate(SymmetryPacked1RDM(Len_1RDM))
        allocate(SymmetryPackedLagrangian(Len_1RDM))
        allocate(SymmetryPacked2RDM(Len_2RDM))
        allocate(FC_Lagrangian(FCLag_Len))

        FC_Lagrangian(:)=0  !Frozen-core Lagrandian -- whilst we do all electron calcs
        
        ! Constructing the Symmetry Packed arrays
        ! We convert our 1RDM, Lagrangian and  2RDM into the required Molpro symmetry-packed format
        ! 2RDM is stored as a full square matrix, separated into symmetry blocks
        ! We store D_ijkl (chemical notation, spatial orbs) where (i .ge. j) and (k .ge. l)
        ! For each symmetry X, there will be a block where (i,j) and (k,l) both have symmetry X
        ! making (ij,kl) totally symmetric, and D_ijkl (potentially) non-zero.
        ! 1RDM and Lagrangian are stored as upper triangles, separated by symmetry block

        SymmetryPacked2RDM(:)=0.0_dp
        SymmetryPacked1RDM(:)=0.0_dp
        SymmetryPackedLagrangian(:) = 0.0_dp

        do i=1, SpatOrbs  !run over spatial orbitals, ALL ELECTRON ONLY
            do j=1,  i    ! i .ge. j
                Sym_i=SpinOrbSymLabel(2*i)  !Consider only alpha orbitals
                Sym_j=SpinOrbSymLabel(2*j)
                Sym_ij=RandExcitSymLabelProd(Sym_i, Sym_j)
                if (Sym_ij .eq. 0) then
                    posn1=blockstart1(Sym_i+1)+elements_assigned1(Sym_i+1)
                    ! Add pre-symmetrised contribution to the symmetry-packed 1RDM
                    if(tStoreSpinOrbs) then
                        !Include both aa and bb contributions
                        SymmetryPacked1RDM(posn1)=2.0_dp*&
                                (NatOrbMat(SymLabelListInv_rot(2*i),SymLabelListInv_rot(2*j)))*Norm_1RDM
                    else
                        SymmetryPacked1RDM(posn1)=NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM
                    endif

                    !Add in the symmetrised Lagrangian contribution to the sym-packed Lagrangian
                    SymmetryPackedLagrangian(posn1)=Lagrangian(i,j)
                    elements_assigned1(Sym_i+1)=elements_assigned1(Sym_i+1)+1
                endif
                do k=1, SpatOrbs
                    do l=1, k
                        Sym_k=SpinOrbSymLabel(2*k)  !Consider only alpha orbitals
                        Sym_l=SpinOrbSymLabel(2*l)
                        Sym_kl=RandExcitSymLabelProd(Sym_k, Sym_l)

                        if(Sym_kl .eq. Sym_ij) then
                            posn2=blockstart2(Sym_ij+1)+elements_assigned2(Sym_ij+1)
                        
                            !Extracts spat orb chemical notation term from spin separated physical notation RDMs 
                            ijkl=Find_Spatial_2RDM_Chem(i,j,k,l, Norm_2RDM)
                            jikl=Find_Spatial_2RDM_Chem(j,i,k,l, Norm_2RDM)

                            SymmetryPacked2RDM(posn2)=0.5*(ijkl+jikl)
                            elements_assigned2(Sym_ij+1)=elements_assigned2(Sym_ij+1)+1
                        endif
                    enddo
                enddo
            enddo
        enddo
        
        call molpro_dump_mcscf_dens_for_grad(myname,ifil, &
            icore, iclos, iact, nEL, isyref, ms2, iblkq,&
            iseccr, istat1, SymmetryPacked1RDM, Len_1RDM, &
            SymmetryPacked2RDM, Len_2RDM, SymmetryPackedLagrangian, &
            Len_1RDM, FC_Lagrangian, FC_Lag_Len, &
            iout, intrel)
        
    endif
        
    end subroutine

   function Find_Spatial_2RDM_Chem(p,q,r,s, Norm_2RDM) result(pqrs)

   !This routine calculates the spatial orbital, chemical notation 2RDM component
   !                 D_pqrs = <Psi | p+ r+ s q | Psi>
   !This is achieved by looking up the D_pr,qs component in the various spin-separated
   !versions of the 2RDM that are currently stored with spatial orb numbering, in 
   !physical notation (e.g. All_aaaa_RDM etc).

   !To convert from spin to spatial orbitals, we need to apply the following:
   !D_pr,qs = D_pr,qs(aaaa) + D_pr,qs(bbbb) + D_pr,qs(abab) + D_pr,qs(baba) (Eq. ***)

   !We note now the following quirks of the All_aaaa_RDM-type arrays for the manner in
   !which they store these components
   !    1. In most cases the current RDMs store the *sum* of the spin-inverted terms
   !         - ie, All_aaaa_RDM(pr,qs) contains the sum of the aaaa and bbbb contributions
   !    2. When p=r and q=s, there is only one contribution generated in NECI
   !         - ie, All_abab_RDM(pp,qq) contains only one of the two identical abab and baba contributions
   !         - Terms of this kind but be explicitly multiplied by two to satisfy Eq. *** above
   !         - This is stored in the "Mult_Factor"
   !    3. The existing 2RDMs only store terms with r>=p and s>=q
   !         - If we wish to look up a term with a different order to this, we must swap the
   !           order of the indices, considering the swapped spin and introducing appropriate signs
   !         - ie if p>r and s>q, D_pr,qs(abab) is found by looking up -D_rp,qs(abba)
    
   integer, intent(in) :: p,q,r,s
   real(dp), intent(in) :: Norm_2RDM
   real(dp) :: pqrs
   real(dp) :: Mult_Factor
    integer :: Ind1_aa, Ind1_ab, Ind2_aa, Ind2_ab 
   pqrs=0.0_dp

    if((p.eq.r).and.(q.eq.s)) then
        Mult_Factor = 2.0_dp
    else
        Mult_Factor = 1.0_dp
    endif
    
    if ((r.ge.p) .and. (s.ge.q)) then !D_pr,qs correctly ordered
        Ind1_aa = ( ( (r-2) * (r-1) ) / 2 ) + p
        Ind1_ab = ( ( (r-1) * r ) / 2 ) + p
        Ind2_aa = ( ( (s-2) * (s-1) ) / 2 ) + q
        Ind2_ab = ( ( (s-1) * s ) / 2 ) + q
        pqrs=pqrs+Mult_Factor*All_abab_RDM(Ind1_ab,Ind2_ab)*Norm_2RDM
        if ((p.ne.r) .and. (q.ne.s)) pqrs=pqrs+Mult_Factor*All_aaaa_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM
    elseif ((p.gt.r) .and. (s.gt.q)) then !Need to reorder D_pr,qs to -D_rp,qs
        Ind1_aa = ( ( (p-2) * (p-1) ) / 2 ) + r
        Ind2_aa = ( ( (s-2) * (s-1) ) / 2 ) + q
        pqrs=pqrs-Mult_Factor*All_abba_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM
        pqrs=pqrs-Mult_Factor*All_aaaa_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM
    elseif ((r.gt.p) .and. (q.gt.s)) then !Need to reorder D_pr,qs to -D_pr,sq
        Ind1_aa = ( ( (r-2) * (r-1) ) / 2 ) + p
        Ind2_aa = ( ( (q-2) * (q-1) ) / 2 ) + s
        pqrs=pqrs-Mult_Factor*All_abba_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM
        pqrs=pqrs-Mult_Factor*All_aaaa_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM
    else !Must need to reorder D_pr,qs to D_rp,sq
        Ind1_aa = ( ( (p-2) * (p-1) ) / 2 ) + r
        Ind1_ab = ( ( (p-1) * p ) / 2 ) + r
        Ind2_aa = ( ( (q-2) * (q-1) ) / 2 ) + s
        Ind2_ab = ( ( (q-1) * q ) / 2 ) + s
        pqrs=pqrs+Mult_Factor*All_abab_RDM(Ind1_ab,Ind2_ab)*Norm_2RDM
        if ((p.ne.r) .and. (q.ne.s)) pqrs=pqrs+Mult_Factor*All_aaaa_RDM(Ind1_aa,Ind2_aa)*Norm_2RDM
    endif
        
    end function

    subroutine molpro_set_igrdsav(igrsav_)
        implicit double precision(a-h,o-z)
        implicit integer(i-n)
        character(len=*), parameter :: t_r='molpro_set_igrdsav'
#ifdef MOLPRO
        include "common/cwsave"
        ! tell gradient programs where gradient information is stored
        igrsav = igrsav_
#else
        call stop_all(t_r,'Should be being called from within MOLPRO')
#endif
    end subroutine
    
    subroutine molpro_get_reference_info(iWfRecord, iWfSym)
        integer :: iWfSym,iWfRecord
#ifdef MOLPRO
        include "common/code"
        include "common/cref"
        ! wf(1) should contain the record number, as integer, of the last used orbital set.
        iWfRecord = wf(1)
        iWfSym = isyref
#else
        iWfRecord = 21402
#endif
    end subroutine
    
    function molpro_get_iout() result(iiout)
        integer :: iiout
#ifdef MOLPRO
        include "common/tapes"
        ! output i/o unit for main output (usually 6, but might be different
        ! if in logfile mode, e.g., during geometry optimization)
        iiout = iout
#else
        iiout = 6
#endif
    end function
    

    subroutine molpro_dump_mcscf_dens_for_grad(name,ifil, &
          icore,iclos,iact,nelec,isyref,ms2,iblkq,iseccr,istat1, &
          den1,lden1, &
          den2,lden2, &
          eps,leps, &
          epsc,lepsc,iout,intrel)
      implicit double precision(a-h,o-z)
    
    !  Use the following subroutine to dump the information needed for molpro
    !  forces calculations, such that it is in the exact format that molpro
    !  would have dumped from a MCSCF calculation.  When interfaced with Molpro,
    !  I understand that the molpro equivalent of this routine should be called. 
    
    !  implicit double precision (a-h,o-z)
      integer,                            intent(inout) :: name
      integer,                            intent(inout) :: ifil
      integer, dimension(8),              intent(in)    :: icore
      integer, dimension(8),              intent(in)    :: iclos
      integer, dimension(8),              intent(in)    :: iact
      integer,                            intent(in)    :: nelec
      integer,                            intent(in)    :: isyref
      integer,                            intent(in)    :: ms2
      integer,                            intent(in)    :: iblkq
      integer,                            intent(in)    :: iseccr
      integer,                            intent(in)    :: istat1
      integer,                            intent(in)    :: lden1
      real(dp), dimension(lden1), intent(in)    :: den1
      integer,                            intent(in)    :: lden2
      real(dp), dimension(lden2), intent(in)    :: den2
      integer,                            intent(in)    :: leps
      real(dp), dimension(leps),  intent(in)    :: eps
      integer,                            intent(in)    :: lepsc
      real(dp), dimension(lepsc), intent(in)    :: epsc
      integer,                            intent(in)    :: iout
      integer,                            intent(in)    :: intrel
      integer                                           :: igrsav

    !> Write mcscf density in format needed for gradient program.
    !> \param[in,out] name record number to be written
    !> \param[in,out] ifil file number to be written
    !> \param[in] icore numbers of frozen core orbitals in each symmetry
    !> \param[in] iclos plus numbers of closed-shell orbitals in each symmetry
    !> \param[in] iact numbers of active orbitals in each symmetry
    !> \param[in] nelec number of active electrons -- DONE
    !> \param[in] isyref spatial symmetry of wavefunction -- DONE
    !> \param[in] ms2 spin quauntum number times 2 -- DONE
    !> \param[in] iblkq record number * 10 + file number for orbitals (typically 21402) -- DONE
    !> \param[in] iseccr record number * 10 + file number for frozen orbitals (typically 21002) -- DONE
    !> \param[in] istat1 state number (1 for ground state) -- DONE
    !> \param[in] den1 1-particle density matrix
    !> \param[in] lden1 size of den1
    !> \param[in] den2 2-particle density matrix
    !> \param[in] lden2 size of den2
    !> \param[in] eps Lagrangian
    !> \param[in] leps size of leps
    !> \param[in] epsc Lagrangian for frozen core
    !> \param[in] lepsc size of epsc
    !> \param[in] iout unit for output, eg. 6
    !> \param[in] intrel byte size ratio of integer to double precision, normally 1 or 2
    !> \param[out] igrsav where the record was written (needed in common/cwsave)

      integer, dimension(30)      :: header
      integer :: i, ncore
      character(len=6), parameter :: label='MCGRAD'
      integer :: lhead

#ifndef MOLPRO
          call stop_all(t_r,'Should not be here if not running through molpro')
#endif

      ncore = 0
      do i=1,8
       ncore = ncore + icore(i)
       header(1-1+i)=icore(i)
       header(1+7+i)=iclos(i)
       header(1+15+i)=iact(i)
      enddo
      header(1+24)=nelec
      header(1+25)=isyref
      header(1+26)=ms2
      header(1+27)=iabs(iblkq)
      header(1+28)=iabs(iseccr)
      header(1+29)=iabs(istat1)
      lhead=30/intrel
      igrsav=10*(name-1)+ifil
      
      name=igrsav/10
      ifil=igrsav-10*name
#ifdef MOLPRO
      call reserv(lhead+lden1+lden2+leps+lepsc,ifil,name,-1)
      call writem(header,lhead,ifil,name,0,label)
      call writem(den1,lden1,ifil,name,lhead,label)
      call writem(den2,lden2,ifil,name,lhead+lden1,label)
      call writem(eps,leps,ifil,name,lhead+lden1+lden2,label)
      if(ncore.ne.0) call writem(epsc,lepsc,ifil,name,&
                    lhead+lden1+lden2+leps,label)
#endif
      write(iout,20) istat1,isyref,name,ifil
!       call setf(2,recmc,1)
      ! ^- cgk: not sure what this does.
20    format(/' Gradient information for state',i2,'.',i1,&
                  ' saved on record  ',i8,'.',i1)
      call molpro_set_igrdsav(igrsav)

      !igrsav=nexfre(igrsav)
      !name=igrsav/10
      !ifil=igrsav-10*name
      !call reserv(lhead+lden1+lden2+leps+lepsc,ifil,name,-1)
!      open (unit = ifil, file = "fciqmc_forces_info", form='UNFORMATTED', access='sequential')
!      write(ifil) header, den1, den2, eps
      !call writem(den2,lden2,ifil,name,lhead+lden1,label)
      !call writem(eps,leps,ifil,name,lhead+lden1+lden2,label)
      !if(ncore.ne.0) call writem(epsc,lepsc,ifil,name,
      !>                          lhead+lden1+lden2+leps,label)
!      write(iout,*) "istat1, isyref, name, ifil", istat1,isyref,name,ifil
!      write(iout,*) "header, den1, den2, eps", header, den1, den2, eps
!20    format(/' Gradient information for state',i2,'.',i1,
     !>        ' saved on record  ',i8,'.',i1)
     ! return
      end subroutine

!Do a binary search in CurrentDets, between the indices of MinInd and MaxInd. If successful, tSuccess will be true and 
!PartInd will be a coincident determinant. If there are multiple values, the chosen one may be any of them...
!If failure, then the index will be one less than the index that the particle would be in if it was present in the list.
!(or close enough!)
    SUBROUTINE BinSearchParts_rdm(iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER(KIND=n_int) :: iLut(0:NIfTot)
        INTEGER :: MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(CurrentDets(:,MaxInd),iLut(:),NIfDBO)
            IF(Comp.eq.0) THEN
                tSuccess=.true.
                PartInd=MaxInd
                RETURN
            ELSE
                tSuccess=.false.
                PartInd=MinInd
            ENDIF
        ENDIF

        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(CurrentDets(:,N),iLut(:),NIfDBO)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower 
!bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final 
!entry of the list and leave
!We need to check the last index.
                    Comp=DetBitLT(CurrentDets(:,i+1),iLut(:),NIfDBO)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN
                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchParts_rdm
    
    subroutine fill_RDM_offdiag_deterministic() 
        integer :: i, j
        integer, dimension(2) :: SingEx
        integer, dimension(2,2) :: Ex
        real(dp) :: InstSignI, InstSignJ
        real(dp) :: AvSignI, AvSignJ
        logical :: tParity
        integer(kind=n_int) :: iLutI(0:niftot), iLutJ(0:niftot)
        integer :: nI(nel), nJ(nel), IC
        integer :: IterRDM, connect_elem
   
        if(mod((Iter - IterRDMStart + 1),RDMEnergyIter).eq.0) then
            IterRDM=RDMEnergyIter
        else
            !This must be the final iteration, as we've got tFill_RDM=.true.
            !for an iteration where we wouldn't normally need the energy
            IterRDM=mod((Iter - IterRDMStart + 1),RDMEnergyIter)
        endif

        Ex(:,:)=0

        do i = 1, determ_proc_sizes(iProcIndex) !Core dets on this proc
            iLutI=core_space(:,determ_proc_indices(iProcIndex)+i)
                         
            if(DetBitEq(iLutI,iLutHF_True,NifDBO)) cycle  !Connections to HF done elsewhere.
           
!            InstSignI=partial_determ_vector(i)
            AvSignI=full_determ_vector_av(determ_proc_indices(iProcIndex)+i)

            CALL decode_bit_det(nI,iLutI)
                        
            do j = 1, sparse_core_ham(i)%num_elements-1   
                 !Running over all non-zero off-diag matrix elements
                 !Connections to whole space (1 row), excluding diagonal elements

                 !Note: determ_proc_indices holds sum(determ_proc_sizes(0:proc-1))
                 !Core space holds all the core determinants on every processor, so we need
                 !to shuffle up to the range of indices corresponding to this proc
                 ! (using determ_proc_indices) and then select the correct one, i
                 
                 iLutJ=core_space(:,core_connections(i)%positions(j))
                 if(DetBitEq(iLutJ,iLutHF_True,NifDBO)) cycle
                 !Connections to HF done elsewhere.
                 
!                 InstSignJ=full_determ_vector(core_connections(i)%positions(j))
                 AvSignJ=full_determ_vector_av(core_connections(i)%positions(j))
                 
                 connect_elem=core_connections(i)%elements(j)

                 IC=abs(connect_elem)

                 if(sign(1, connect_elem).gt.0) then
                     tParity=.false.
                 else
                     tParity=.true.
                 endif

                 if(tHPHF) then

                     call decode_bit_det(nJ, iLutJ)

                     call Fill_Spin_Coupled_RDM_v2(iLutI, iLutJ, nI, nJ, AvSignI*IterRDM, AvSignJ, .false.)
                     !if(IC.eq.1) then

                 else
                     if(IC.eq.1) then
                         !Single excitation - contributes to 1- and 2-RDM (if calculated)
                          
                         !Note: get_bit_excitmat may be buggy (DetBitOps), but will do for now as we need the Ex...
                         call get_bit_excitmat(iLutI,iLutJ,SingEx,IC)
                         Ex(:,1)=SingEx(:)
                        
                         !No need to explicitly fill symmetrically as we'll generate pairs of 
                         ! determinants both ways around using the connectivity matrix.
                         call Fill_Sings_RDM(nI,Ex,tParity,AvSignI*IterRDM,AvSignJ,.false.)

                     elseif((IC.eq.2).and.(RDMExcitLevel.ne.1)) then
                         
                         !Note: get_bit_excitmat may be buggy (DetBitOps), but will do for now as we need the Ex...
                         call get_bit_excitmat(iLutI,iLutJ,Ex,IC)
                                     
                         call Fill_Doubs_RDM(Ex,tParity,AvSignI*IterRDM,AvSignJ,.false.)
                     endif
                 endif
             end do
        end do
    end subroutine fill_RDM_offdiag_deterministic 

    !The only thing needed is the 1RDM (normalized)
    subroutine CalcDipoles(Norm_1RDM)
#ifdef MOLPRO
        use outputResult
        use SymData, only : Sym_Psi,nSymLabels
        use GenRandSymExcitNUMod , only : RandExcitSymLabelProd, ClassCountInd
        use SymExcitDataMod, only: SpinOrbSymLabel,SymLabelCounts2
        implicit none
        integer, dimension(nSymLabels) :: elements_assigned1, blockstart1
        real(dp) :: dipmom(3),znuc,zcor
        real(dp), allocatable :: SymmetryPacked1RDM(:),ints(:)
        integer :: i,j,ipr,Sym_i,Sym_j,Sym_ij,posn1,isize,isyref,mxv,iout
!        include "common/maxatm"
!        include "common/tapes"
!        include "common/dumpinfow"
!        include "common/cstate"
!        include "common/maxbfn"
!        include "common/corbdim"
!        include "common/casscf"
!        include "common/syminf"
!        include "common/jobopt"
!        include "common/big"
!        include "common/cbas"
!        include "common/clseg"
!        include "common/cref"
!        include "common/ctran2"
!        include "common/code"
!        include "common/cmpp"
!        include "common/d2gen_cvb"
#else
        implicit none
#endif
        real(dp), intent(in) :: Norm_1RDM
        character(len=*), parameter :: t_r='CalcDipoles'

#ifdef MOLPRO

        if(iProcIndex.eq.0) then
            iout=molpro_get_iout()
            !We need to work out a) how molpro symmetry-packs UHF integrals (ROHF would be fine though)
            !b) Ensure that the 1RDM is correctly calculated for UHF (It is always allocated as spatorbs)
            !c) Modify this routine for contracting over spin-orbitals
            if(tStoreSpinOrbs) call stop_all(t_r,'Not working for ROHF/UHF')

            isyref=Sym_Psi+1  !spatial symmetry of the wavefunction

            !Size of symmetry packed arrays (spatial)
            isize = 0
            blockstart1(:) = 0
            do i = 0,nSymLabels-1
                !Find position of each symmetry block in sym-packed forms of RDMS 1 & 2
                blockstart1(i+1)=isize+1 !N.B. Len_1RDM still being updated in this loop

                isize = isize + (SymLabelCounts2(2,ClassCountInd(1,i,0))*   &
                    (SymLabelCounts2(2,ClassCountInd(1,i,0))+1))/2 !Counting alpha orbitals
            enddo
!            write(6,*) "Size of symmetry packed 1-electron array",isize
            allocate(ints(isize))

            elements_assigned1(:) = 0
            allocate(SymmetryPacked1RDM(isize))
            SymmetryPacked1RDM(:) = 0.0_dp
            do i = 1,SpatOrbs !run over spatial orbitals, ALL ELECTRON ONLY
                do j = 1,i ! i .ge. j
                    Sym_i=SpinOrbSymLabel(2*i)  !Consider only alpha orbitals
                    Sym_j=SpinOrbSymLabel(2*j)
                    Sym_ij=RandExcitSymLabelProd(Sym_i, Sym_j)
                    if (Sym_ij .eq. 0) then
                        posn1=blockstart1(Sym_i+1)+elements_assigned1(Sym_i+1)

                        SymmetryPacked1RDM(posn1)=NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM
                        if(i.ne.j) then
                            !Double the off-diagonal elements of the 1RDM, so that when we contract over the 
                            !symmetry packed representation of the 1RDM, it is as if we are also including the other half of the matrix
                            SymmetryPacked1RDM(posn1) = 2.0_dp*SymmetryPacked1RDM(posn1)
                        endif
                    endif
                    elements_assigned1(Sym_i+1)=elements_assigned1(Sym_i+1)+1

                enddo
            enddo
                
!            !Create and dump integrals? 
!Not currently working. Currently requires them to be dumped from casscf calc already
!            n1elec = 5
!            i1elec(1) = 2
!            i1elec(2) = 4
!            i1elec(3) = 5
!            i1elec(4) = 6
!            i1elec(5) = 3
!            write(6,*) "Calling mukint..."
!            call flush(6)
!            call mukint(0,0)
!            write(6,*) "Exiting mukint..."
!            call flush(6)
!            itrsfm=1
!            write(6,*) "Symmetry packed 1RDM: ",SymmetryPacked1RDM(:)
            dipmom(:) = 0.0_dp

            call clearvar('DMX')
            call clearvar('DMY')
            call clearvar('DMZ')

            do ipr=4,6
            
                ints(isize) = 0.0_dp

                call pget(ints,ipr,znuc,zcor)

                !write(6,*) "Integrals for property ",ipr-3
                !write(6,*) ints(1:isize)

                !Now, contract
                do i = 1,isize
                    dipmom(ipr-3) = dipmom(ipr-3) - ints(i)*SymmetryPacked1RDM(i)
                enddo
                dipmom(ipr-3) = dipmom(ipr-3) + znuc - zcor
            enddo
            write(iout,"(A)") ""
            write(iout,"(A,3f15.8)") "DIPOLE MOMENT: ",dipmom(1:3)
            write(iout,"(A)") ""
            call output_result('FCIQMC','Dipole moment',dipmom(1:3),1,isyref,numberformat='3f15.8',debye=.TRUE.)
            mxv=1
            call setvar('DMX',dipmom(1),'AU',1,1,mxv,-1)
            call setvar('DMY',dipmom(2),'AU',1,1,mxv,-1)
            call setvar('DMZ',dipmom(3),'AU',1,1,mxv,-1)
            deallocate(ints,SymmetryPacked1RDM)
        endif

#else
        call warning_neci(t_r,'Cannot compute dipole moments if not running within molpro. Exiting...')
#endif

    end subroutine CalcDipoles


END MODULE nElRDMMod


