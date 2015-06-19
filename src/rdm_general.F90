module rdm_general

    ! This module contains general routines related to RDM calculation in
    ! FCIQMC. This includes the initialisation and end routines, and also a
    ! few routines (such as calc_rdmbiasfac and store_parent_with_spawned) used
    ! during the main simulation.

    use bit_rep_data, only: NIfTot, NIfDBO
    use SystemData, only: NEl, nBasis
    use constants

    implicit none

contains

    ! Initialisation routines.

    subroutine InitRDMs(nrdms)

        ! This routine initialises any of the arrays needed to calculate the
        ! reduced density matrix. It is used for both the explicit and
        ! stochastic RDMs.

        use DeterminantData, only: write_det
        use CalcData, only: MemoryFacPart
        use FciMCData, only: MaxSpawned, Spawned_Parents, Spawned_Parents_Index
        use FciMCData, only: Spawned_ParentsTag, Spawned_Parents_IndexTag
        use FciMCData, only: tSinglePartPhase, HFDet_True, tFinalRDMEnergy
        use LoggingData, only: tDo_Not_Calc_RDMEnergy, tNo_RDMs_To_Read, tWrite_RDMs_to_read
        use LoggingData, only: tRDMInstEnergy, RDMExcitLevel, tExplicitAllRDM, tPrint1RDM
        use LoggingData, only: tDiagRDM, tReadRDMs, tPopsfile, tDumpForcesInfo, tDipoles
        use NatOrbsMod, only: NatOrbMat, NatOrbMatTag, Evalues, EvaluesTag
        use Parallel_neci, only: iProcIndex, nProcessors
        use rdm_data, only: rdms, tOpenShell, tCalc_RDMEnergy, Sing_ExcDjs, Doub_ExcDjs
        use rdm_data, only: Sing_ExcDjs2, Doub_ExcDjs2, Sing_ExcDjsTag, Doub_ExcDjsTag
        use rdm_data, only: Sing_ExcDjs2Tag, Doub_ExcDjs2Tag, OneEl_Gap, TwoEl_Gap
        use rdm_data, only: Sing_InitExcSlots, Doub_InitExcSlots, Sing_ExcList, Doub_ExcList
        use rdm_data, only: rdm_estimates_unit, nElRDM_Time, FinaliseRDMs_time, RDMEnergy_time
        use RotateOrbsData, only: SymLabelCounts2_rot,SymLabelList2_rot, SymLabelListInv_rot
        use RotateOrbsData, only: SymLabelCounts2_rotTag, SymLabelList2_rotTag, NoOrbs
        use RotateOrbsData, only: SymLabelListInv_rotTag, SpatOrbs, NoSymLabelCounts
        use SystemData, only: tStoreSpinOrbs, tHPHF, tFixLz, iMaxLz, tROHF
        use util_mod, only: get_free_unit, LogMemAlloc

        integer, intent(in) :: nrdms

        integer :: ierr, i, iproc, rdm_size_1, rdm_size_2
        integer :: MemoryAlloc, MemoryAlloc_Root
        character(len=*), parameter :: t_r = 'InitRDMs'

        ! First thing is to check we're not trying to fill the RDMs in a way
        ! that is not compatible with the code (not every case has been
        ! accounted for yet).

#ifdef __CMPLX
        call stop_all(t_r, 'Filling of reduced density matrices not working with complex walkers yet.')
#endif

        if (nrdms > 1 .and. (tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles .or. RDMExcitLevel == 1)) then
            call stop_all(t_r, 'The RDM feature you have requested is not currently implemented &
                               &with multiple RDMs. In particular, forces, dipole moments, printing &
                               &1-RDMs, diagonalising 1-RDMs, or anything to do with natural orbitals &
                               & has not yet been implemented for more than one RDM.')
        end if

        ! For now, just allocate one rdm.
        allocate(rdms(nrdms))

        ! Only spatial orbitals for the 2-RDMs (and F12).
        if (tStoreSpinOrbs .and. (RDMExcitLevel .ne. 1)) &
            call stop_all(t_r, '2-RDM calculations not set up for systems stored as spin orbitals.')

        if (tROHF .or. tStoreSpinOrbs) then
            tOpenShell = .true.
        else
            tOpenShell = .false.
        end if

        if (tExplicitAllRDM) then
            write(6,'(1X,"Explicitly calculating the reduced density matrices from the FCIQMC wavefunction.")')
        else
            write(6,'(1X,"Stochastically calculating the reduced density matrices from the FCIQMC wavefunction")')
            write(6,'(1X,"incl. explicit connections to the following HF determinant:")', advance='no')
            call write_det (6, HFDet_True, .true.)
        end if

        if (RDMExcitLevel .eq. 1) then
            tCalc_RDMEnergy = .false.
        else
            ! If the RDMExcitLevel is 2 or 3 - and we're calculating the 2-RDM, 
            ! then we automatically calculate the energy unless we specifically
            ! say not to.
            if (tDo_Not_Calc_RDMEnergy) then
                tCalc_RDMEnergy = .false.            
            else
                tCalc_RDMEnergy = .true.
                write(6,'(1X,"Calculating the energy from the reduced density matrix. &
                              &This requires the 2 electron RDM from which the 1-RDM can also be constructed.")')
            end if
        end if

        ! Have not got HPHF working with the explicit or truncated methods yet.
        ! Neither of these would be too difficult to implement.
        if (tHPHF .and. tExplicitAllRDM) call stop_all(t_r, 'HPHF not set up with the explicit calculation of the RDM.')

        SpatOrbs = nBasis/2
        if (tOpenShell) then
            NoOrbs = nBasis
        else
            NoOrbs = SpatOrbs
        end if

        ! There are two different size arrays to allocate, as given by these
        ! array sizes (the RDMs are there sizes squared).
        rdm_size_1 = (SpatOrbs*(SpatOrbs-1))/2
        rdm_size_2 = (SpatOrbs*(SpatOrbs+1))/2

        ! Here we're allocating arrays for the actual calculation of the RDM.
        MemoryAlloc = 0
        MemoryAlloc_Root = 0 ! Memory allocated in bytes.

        ! First for the storage of the actual 1- or 2-RMD.
        if (RDMExcitLevel .eq. 1) then

            ! This is the AllnElRDM, called NatOrbMat simply because we use the
            ! natural orbital routines to diagonalise etc. We don't have an
            ! instantaneous 1-RDM.
            allocate(NatOrbMat(NoOrbs, NoOrbs), stat=ierr)
            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating NatOrbMat array,')
            call LogMemAlloc('NatOrbMat', NoOrbs**2, 8, t_r, NatOrbMatTag, ierr)
            NatOrbMat(:,:) = 0.0_dp

            MemoryAlloc = MemoryAlloc + ( NoOrbs * NoOrbs * 8 )
            MemoryAlloc_Root = MemoryAlloc_Root + ( NoOrbs * NoOrbs * 8 )
        else
            ! If we're calculating the 2-RDM, the 1-RDM does not need to be
            ! calculated as well because all its info is in the 2-RDM anyway.

            ! The 2-RDM of the type alpha alpha alpha alpha (= beta beta beta beta).
            ! These *do not* include any 2-RDM(i,j,a,b) terms where i=j or a=b (if
            ! they're the same spin this can't happen).

            do i = 1, nrdms
                if (tRDMInstEnergy) then
                    ! We will be filling up rdms(i)%aaaa_inst as we go along, which need
                    ! to be allocated per core.
                    ! When calculating the energy, these will be summed over cores
                    ! using an _inplace type command.
                    ! To calculate the full energy of the RDM (i.e. over full accum.
                    ! period), we need to allocate rdms(i)%aaaa_full on the head nodes.

                    allocate(rdms(i)%aaaa_inst(rdm_size_1, rdm_size_1), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r,'Problem allocating aaaa_inst RDM array,')
                    call LogMemAlloc('rdms(i)%aaaa_inst', (rdm_size_1**2), 8, t_r, rdms(i)%aaaa_instTag, ierr)
                    rdms(i)%aaaa_inst(:,:) = 0.0_dp

                    ! The 2-RDM of the type alpha beta beta alpha (= beta alpha alpha beta).
                    ! These also *do not* also include 2-RDM(i,j,a,b) terms where i=j or a=b
                    ! (these are the same as the abab elements).
                    allocate(rdms(i)%abba_inst(rdm_size_1, rdm_size_1), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating abba_inst RDM array,')
                    call LogMemAlloc('rdms(i)%abba_inst', (rdm_size_1**2), 8, t_r, rdms(i)%abba_instTag, ierr)
                    rdms(i)%abba_inst(:,:) = 0.0_dp

                    MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                    MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8

                    ! The 2-RDM of the type alpha beta alpha beta ( = beta alpha beta alpha).
                    ! These *do* include 2-RDM(i,j,a,b) terms where i=j or a=b, if they're
                    ! different spin this is possible - hence the slightly different size to
                    ! the aaaa array.
                    allocate(rdms(i)%abab_inst(rdm_size_2, rdm_size_2), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating abab_inst RDM array,')
                    call LogMemAlloc('rdms(i)%abab_inst', (rdm_size_2**2), 8, t_r, rdms(i)%abab_instTag, ierr)
                    rdms(i)%abab_inst(:,:) = 0.0_dp

                    MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8
                    MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8 

                    ! Extra arrays for open shell systems.
                    if (tOpenShell) then
                        allocate(rdms(i)%bbbb_inst(rdm_size_1, rdm_size_1), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating bbbb_inst RDM array,')
                        call LogMemAlloc('rdms(i)%bbbb_inst', (rdm_size_1**2), 8, t_r, rdms(i)%bbbb_instTag, ierr)
                        rdms(i)%bbbb_inst(:,:) = 0.0_dp

                        allocate(rdms(i)%baab_inst(rdm_size_1, rdm_size_1), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating baab_inst RDM array,')
                        call LogMemAlloc('rdms(i)%baab_inst', (rdm_size_1**2), 8, t_r, rdms(i)%baab_instTag, ierr)
                        rdms(i)%baab_inst(:,:) = 0.0_dp
                        MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8

                        allocate(rdms(i)%baba_inst(rdm_size_2, rdm_size_2), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating baba_inst RDM array,')
                        call LogMemAlloc('rdms(i)%baba_inst', (rdm_size_2**2), 8, t_r, rdms(i)%baba_instTag, ierr)
                        rdms(i)%baba_inst(:,:) = 0.0_dp
                        MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8
                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8
                    end if

                    if (iProcIndex .eq. 0) then
                        allocate(rdms(i)%aaaa_full(rdm_size_1, rdm_size_1), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating aaaa_full RDM array,')
                        call LogMemAlloc('rdms(i)%aaaa_full', (rdm_size_1**2), 8, t_r, rdms(i)%aaaa_fullTag, ierr)
                        rdms(i)%aaaa_full(:,:) = 0.0_dp

                        allocate(rdms(i)%abab_full(rdm_size_2, rdm_size_2), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating abab_full RDM array,')
                        call LogMemAlloc('rdms(i)%abab_full', (rdm_size_2**2), 8, t_r, rdms(i)%abab_fullTag, ierr)
                        rdms(i)%abab_full(:,:) = 0.0_dp

                        allocate(rdms(i)%abba_full(rdm_size_1, rdm_size_1), stat = ierr)
                        if (ierr .ne. 0) call stop_all(t_r,'Problem allocating abba_full RDM array,')
                        call LogMemAlloc('rdms(i)%abba_full', (rdm_size_1**2), 8, t_r, rdms(i)%abba_fullTag, ierr)
                        rdms(i)%abba_full(:,:) = 0.0_dp

                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8
                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8
                        MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                        MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8

                        if (tOpenShell) then
                            allocate(rdms(i)%bbbb_full(rdm_size_1, rdm_size_1), stat=ierr)
                            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating bbbb_full RDM array,')
                            call LogMemAlloc('rdms(i)%bbbb_full', (rdm_size_1**2), 8, t_r,rdms(i)%bbbb_fullTag,ierr)
                            rdms(i)%bbbb_full(:,:) = 0.0_dp

                            allocate(rdms(i)%baba_full(rdm_size_2,rdm_size_2),stat=ierr)
                            if (ierr .ne. 0) call stop_all(t_r,'Problem allocating baba_full RDM array,')
                            call LogMemAlloc('rdms(i)%baba_full', (rdm_size_2**2), 8,t_r, rdms(i)%baba_fullTag, ierr)
                            rdms(i)%baba_full(:,:) = 0.0_dp

                            allocate(rdms(i)%baab_full(rdm_size_1, rdm_size_1), stat=ierr)
                            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating baab_full array,')
                            call LogMemAlloc('rdms(i)%baab_full', (rdm_size_1**2), 8,t_r, rdms(i)%baab_fullTag, ierr)
                            rdms(i)%baab_full(:,:) = 0.0_dp

                            MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8
                            MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8
                            MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                            MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8

                        end if

                    end if
                    
                    rdms(i)%aaaa => rdms(i)%aaaa_inst
                    rdms(i)%abba => rdms(i)%abba_inst
                    rdms(i)%abab => rdms(i)%abab_inst

                    if (tOpenShell) then
                        rdms(i)%bbbb => rdms(i)%bbbb_inst
                        rdms(i)%baab => rdms(i)%baab_inst
                        rdms(i)%baba => rdms(i)%baba_inst
                    end if

                else
                    ! We're not calculating an instantaneous RDM energy.
                    ! Put RDM contributions directly into 'full' arrays, which are
                    ! now allocated every core
                    allocate(rdms(i)%aaaa_full(rdm_size_1,rdm_size_1), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating aaaa_full RDM array,')
                    call LogMemAlloc('rdms(i)%aaaa_full', (rdm_size_1**2), 8, t_r,rdms(i)%aaaa_fullTag,ierr)
                    rdms(i)%aaaa_full(:,:) = 0.0_dp

                    ! The 2-RDM of the type alpha beta beta alpha ( = beta alpha alpha beta).
                    ! These also *do not* also include 2-RDM(i,j,a,b) terms where i=j or a=b
                    ! (these are the same as the abab elements).
                    allocate(rdms(i)%abba_full(rdm_size_1,rdm_size_1),stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r,'Problem allocating abba_full RDM array,')
                    call LogMemAlloc('rdms(i)%abba_full', (rdm_size_1**2), 8, t_r, rdms(i)%abba_fullTag, ierr)
                    rdms(i)%abba_full(:,:) = 0.0_dp

                    MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                    MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8

                    ! The 2-RDM of the type alpha beta alpha beta ( = beta alpha beta alpha).
                    ! These *do* include 2-RDM(i,j,a,b) terms where i=j or a=b, if they're
                    ! different spin this is possible - hence the slightly different size
                    ! to the aaaa array.
                    allocate(rdms(i)%abab_full(rdm_size_2, rdm_size_2), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating abab_full RDM array,')
                    call LogMemAlloc('rdms(i)%abab_full', (rdm_size_2**2), 8, t_r, rdms(i)%abab_fullTag, ierr)
                    rdms(i)%abab_full(:,:) = 0.0_dp

                    MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8
                    MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8
                    
                    rdms(i)%aaaa => rdms(i)%aaaa_full
                    rdms(i)%abba => rdms(i)%abba_full
                    rdms(i)%abab => rdms(i)%abab_full

                    if (tOpenShell) then
                        allocate(rdms(i)%bbbb_full(rdm_size_1, rdm_size_1), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating bbbb_full RDM array,')
                        call LogMemAlloc('rdms(i)%bbbb_full', (rdm_size_1**2), 8, t_r, rdms(i)%bbbb_fullTag, ierr)
                        rdms(i)%bbbb_full(:,:) = 0.0_dp

                        allocate(rdms(i)%baab_full(rdm_size_1,rdm_size_1),stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating baab_full array,')
                        call LogMemAlloc('rdms(i)%baab_full', (rdm_size_1**2), 8, t_r, rdms(i)%baab_fullTag, ierr)
                        rdms(i)%baab_full(:,:) = 0.0_dp
                        MemoryAlloc = MemoryAlloc + (rdm_size_1**2)*2*8
                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_1**2)*2*8

                        allocate(rdms(i)%baba_full(rdm_size_2, rdm_size_2), stat=ierr)
                        if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating baba_full RDM array,')
                        call LogMemAlloc('rdms(i)%baba_full', (rdm_size_2**2), 8, t_r, rdms(i)%baba_fullTag, ierr)
                        rdms(i)%baba_full(:,:) = 0.0_dp
                        MemoryAlloc = MemoryAlloc + (rdm_size_2**2)*8
                        MemoryAlloc_Root = MemoryAlloc_Root + (rdm_size_2**2)*8
          
                        rdms(i)%bbbb => rdms(i)%bbbb_full
                        rdms(i)%baab => rdms(i)%baab_full
                        rdms(i)%baba => rdms(i)%baba_full
                    end if
                end if ! Not instantaneous
            end do ! Looping over all RDMs.

            if (iProcindex .eq. 0) then
                if (tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) then
                    ! Still need to allocate 1-RDM to get nat orb occupation numbers.
                    allocate(NatOrbMat(NoOrbs, NoOrbs), stat=ierr)
                    if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating NatOrbMat array,')
                    call LogMemAlloc('NatOrbMat', NoOrbs**2,8, t_r, NatOrbMatTag, ierr)
                    NatOrbMat(:,:) = 0.0_dp
                    MemoryAlloc_Root = MemoryAlloc_Root + ( NoOrbs * NoOrbs * 8 )
                end if
            end if

        end if            

        ! We then need to allocate the arrays for excitations etc when doing
        ! the explicit all calculation.
        if (tExplicitAllRDM) then

            ! We always calculate the single stuff - and if RDMExcitLevel is 1,
            ! this is all, otherwise calculate the double stuff too.

            ! This array actually contains the excitations in blocks of the
            ! processor they will be sent to. Only needed if the 1-RDM is the
            ! only thing being calculated.
            allocate(Sing_ExcDjs(0:NIfTot, nint((NEl*nBasis)*MemoryFacPart)), stat=ierr)
            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Sing_ExcDjs array.')
            call LogMemAlloc('Sing_ExcDjs', nint(NEl*nBasis*MemoryFacPart)*(NIfTot+1),&
                                            size_n_int, t_r, Sing_ExcDjsTag, ierr)

            allocate(Sing_ExcDjs2(0:NIfTot, nint((NEl*nBasis)*MemoryFacPart)), stat=ierr)
            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Sing_ExcDjs2 array.')
            call LogMemAlloc('Sing_ExcDjs2', nint(NEl*nBasis*MemoryFacPart)*(NIfTot+1),&
                                            size_n_int, t_r, Sing_ExcDjs2Tag, ierr)

            Sing_ExcDjs(:,:) = 0
            Sing_ExcDjs2(:,:) = 0

            MemoryAlloc = MemoryAlloc + ( (NIfTot + 1) * nint((NEl*nBasis)*MemoryFacPart) * size_n_int * 2 )
            MemoryAlloc_Root = MemoryAlloc_Root + ( (NIfTot + 1) * nint((NEl*nBasis)*MemoryFacPart) * size_n_int * 2 )

            ! We need room to potentially generate N*M single excitations but
            ! these will be spread across each processor.

            OneEl_Gap = (real(NEl,dp)*real(nBasis,dp)*MemoryFacPart)/real(nProcessors,dp)

            ! This array contains the initial positions of the excitations
            ! for each processor.
            allocate(Sing_InitExcSlots(0:(nProcessors-1)), stat=ierr)
            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Sing_InitExcSlots array,')
            do iproc = 0, nProcessors - 1
                Sing_InitExcSlots(iproc) = nint(OneEl_Gap*iproc) + 1
            end do

            ! This array contains the current position of the excitations as
            ! they're added.
            allocate(Sing_ExcList(0:(nProcessors-1)), stat=ierr)
            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Sing_ExcList array,')
            Sing_ExcList(:) = Sing_InitExcSlots(:)

            if (RDMExcitLevel .ne. 1) then
                ! This array actually contains the excitations in blocks of
                ! the processor they will be sent to.
                allocate(Doub_ExcDjs(0:NIfTot,nint(((NEl*nBasis)**2)*MemoryFacPart)), stat=ierr)
                if (ierr .ne. 0) call stop_all(t_r,'Problem allocating Doub_ExcDjs array.')
                call LogMemAlloc('Doub_ExcDjs', nint(((NEl*nBasis)**2)*MemoryFacPart)&
                                *(NIfTot+1), size_n_int, t_r, Doub_ExcDjsTag, ierr)

                allocate(Doub_ExcDjs2(0:NIfTot, nint(((NEl*nBasis)**2)*MemoryFacPart)), stat=ierr)
                if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Doub_ExcDjs2 array.')
                call LogMemAlloc('Doub_ExcDjs2',nint(((NEl*nBasis)**2)*MemoryFacPart)&
                                *(NIfTot+1), size_n_int, t_r, Doub_ExcDjs2Tag, ierr)

                MemoryAlloc = MemoryAlloc + ( (NIfTot + 1) * nint(((NEl*nBasis)**2)*MemoryFacPart) * size_n_int * 2 )
                MemoryAlloc_Root = MemoryAlloc_Root + ( (NIfTot + 1) * nint(((NEl*nBasis)**2)*MemoryFacPart) * size_n_int * 2 )

                ! We need room to potentially generate (N*M)^2 double excitations
                ! but these will be spread across each processor.        
                TwoEl_Gap = (((real(NEl,dp)*real(nBasis,dp))**2)*MemoryFacPart)/real(nProcessors,dp)

                Doub_ExcDjs(:,:) = 0
                Doub_ExcDjs2(:,:) = 0

                ! This array contains the initial positions of the excitations
                ! for each processor.
                allocate(Doub_InitExcSlots(0:(nProcessors-1)), stat=ierr)
                if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Doub_InitExcSlots array,')
                do iproc = 0, nProcessors-1
                    Doub_InitExcSlots(iproc) = nint(TwoEl_Gap*iproc) + 1
                end do

                ! This array contains the current position of the excitations
                ! as they're added.
                allocate(Doub_ExcList(0:(nProcessors-1)), stat=ierr)
                if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Doub_ExcList array,')
                Doub_ExcList(:) = Doub_InitExcSlots(:)
            end if

        else

            ! Finally, we need to hold onto the parents of the spawned particles.
            ! This is not necessary if we're doing completely explicit calculations.
            allocate(Spawned_Parents(0:(NIfDBO+2), MaxSpawned), stat=ierr)
            if (ierr .ne. 0) call stop_all(t_r,'Problem allocating Spawned_Parents array,')
            call LogMemAlloc('Spawned_Parents', MaxSpawned*(NIfDBO+3), size_n_int,&
                                                t_r,Spawned_ParentsTag, ierr)
            allocate(Spawned_Parents_Index(2,MaxSpawned),stat=ierr)
            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Spawned_Parents_Index array,')
            call LogMemAlloc('Spawned_Parents_Index', MaxSpawned*2,4, t_r,&
                                                        Spawned_Parents_IndexTag, ierr)

            MemoryAlloc = MemoryAlloc + ( (NIfTot + 2) * MaxSpawned * size_n_int ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( (NIfTot + 2) * MaxSpawned * size_n_int ) 

            MemoryAlloc = MemoryAlloc + ( 2 * MaxSpawned * 4 ) 
            MemoryAlloc_Root = MemoryAlloc_Root + ( 2 * MaxSpawned * 4 ) 

        end if

        if (iProcIndex .eq. 0) then
            write(6,"(A,F14.6,A,F14.6,A)") " Main RDM memory arrays consists of : ", &
                    & real(MemoryAlloc_Root,dp)/1048576.0_dp," Mb/Processor on the root, and ", &
                    & real(MemoryAlloc,dp)/1048576.0_dp," Mb/Processor on other processors."
        end if

        ! These parameters are set for the set up of the symmetry arrays, which
        ! are later used for the diagonalisation / rotation of the 1-RDMs.

        if (tOpenShell) then
            if (tFixLz) then
                NoSymLabelCounts = 16 * ( (2 * iMaxLz) + 1 )
            else
                NoSymLabelCounts = 16 
            end if
        else
            if (tFixLz) then
                NoSymLabelCounts = 8 * ( (2 * iMaxLz) + 1 )
            else
                NoSymLabelCounts = 8
            end if
        end if

        if ((RDMExcitLevel .eq. 1) .or. tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) then
            ! These arrays contain indexing systems to order the 1-RDM orbitals
            ! in terms of symmetry. This allows the diagonalisation of the RDMs
            ! to be done in symmetry blocks (a lot quicker/easier).
            ! The 2-RDM does not need to be reordered as it's never diagonalised. 

            allocate(SymLabelCounts2_rot(2, NoSymLabelCounts), stat=ierr)
            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating SymLabelCounts2_rot array,')
            call LogMemAlloc('SymLabelCounts2_rot', 2*NoSymLabelCounts, 4, t_r, SymLabelCounts2_rotTag, ierr)
            SymLabelCounts2_rot(:,:) = 0

            allocate(SymLabelList2_rot(NoOrbs), stat=ierr)
            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating SymLabelList2_rot array,')
            call LogMemAlloc('SymLabelList2_rot', NoOrbs, 4, t_r, SymLabelList2_rotTag, ierr)
            SymLabelList2_rot(:) = 0
     
            allocate(SymLabelListInv_rot(NoOrbs), stat=ierr)
            if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating SymLabelListInv_rot array,')
            call LogMemAlloc('SymLabelListInv_rot', NoOrbs, 4, t_r, SymLabelListInv_rotTag, ierr)
            SymLabelListInv_rot(:) = 0   

            if ((iProcIndex .eq. 0) .and. tDiagRDM) then
                allocate(Evalues(NoOrbs), stat=ierr)
                if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Evalues array,')
                call LogMemAlloc('Evalues', NoOrbs, 8, t_r, EvaluesTag, ierr)
                Evalues(:) = 0.0_dp

                allocate(rdms(1)%Rho_ii(NoOrbs), stat=ierr)
                if (ierr .ne. 0) call stop_all(t_r, 'Problem allocating Rho_ii array,')
                call LogMemAlloc('Rho_ii', NoOrbs, 8, t_r, rdms(1)%Rho_iiTag, ierr)
                rdms(1)%Rho_ii(:) = 0.0_dp
            end if

            ! This routine actually sets up the symmetry labels for the 1-RDM.
            ! TODO: Merge this routine (and rotations later) with the NatOrbs file.
            call SetUpSymLabels_RDM() 

        end if            

        if (iProcIndex .eq. 0) write(6,'(1X,"RDM memory allocation successful...")')

        ! Open file to keep track of RDM Energies (if they're being calculated). 
        if ((iProcIndex .eq. 0) .and. tCalc_RDMEnergy) then
            rdm_estimates_unit = get_free_unit()
            open(rdm_estimates_unit, file='RDMEstimates', status='unknown', position='append')

            write(rdm_estimates_unit, &
                "('#', 4X, 'Iteration', 6X, 'Energy numerator', 6X, 'Spin^2 numerator', 9X, 'Normalisation')")
        end if
        tFinalRDMEnergy = .false.

        ! Reads in the RDMs from a previous calculation, sets the accumulating
        ! normalisations, writes out the starting energy.
        if (tReadRDMs) then
            if (nrdms > 1) call stop_all(t_r, "Reading in multiple RDMs is not yet supported.")
            if (tSinglePartPhase(1) .or. tSinglePartPhase(inum_runs)) then
                write(6,'(A)') 'WARNING - Asking to read in the RDMs, but not varying shift from &
                                & the beginning of the calculation.'
                write(6,'(A)') 'Ignoring the request to read in the RDMs and starting again.'
                tReadRDMs = .false.
            else
                call Read_In_RDMs(rdms(1))
            end if
        end if

        ! By default, if we're writing out a popsfile (and doing an RDM
        ! calculation), we also write out the unnormalised RDMs that can be
        ! read in when restarting a calculation. If the NORDMSTOREAD option
        ! is on, these wont be printed.  
        if (tPopsfile .and. (.not. tno_RDMs_to_read)) twrite_RDMs_to_read = .true.

        nElRDM_Time%timer_name = 'nElRDMTime'
        FinaliseRDMs_Time%timer_name = 'FinaliseRDMsTime'
        RDMEnergy_Time%timer_name = 'RDMEnergyTime'

    end subroutine InitRDMs

    subroutine Read_In_RDMs(rdm)

        ! Reads in the arrays to restart the RDM calculation (and continue
        ! accumulating). These arrays are not normalised, so the trace is
        ! also calculated. The energy is then calculated (if required) from
        ! the RDMs read in only.

        use LoggingData, only: IterRDMonFly
        use LoggingData, only: RDMExcitLevel
        use NatOrbsMod, only: NatOrbMat
        use Parallel_neci, only: iProcIndex
        use rdm_data, only: rdm_t, tOpenShell, tCalc_RDMEnergy
        use rdm_estimators, only: rdm_output_wrapper
        use RotateOrbsData, only: SymLabelListInv_rot
        use SystemData, only: tStoreSpinOrbs
        use util_mod, only: get_free_unit

        type(rdm_t), intent(inout) :: rdm

        logical :: exists_one
        logical :: exists_aaaa, exists_abab, exists_abba
        logical :: exists_bbbb, exists_baba, exists_baab
        integer :: RDM_unit, FileEnd
        integer :: i, j, a, b, Ind1, Ind2
        real(dp) :: Temp_RDM_Element, Norm_2RDM
        character(len=*), parameter :: t_r = 'Read_In_RDMs'

        if (iProcIndex .eq. 0) then 

            if (RDMExcitLevel .eq. 1) then

                write(6,'(1X,"Reading in the 1-RDM")')

                ! The OneRDM will have been printed exactly as is. Without
                ! having been made hermitian, without being normalised, and in
                ! spatial orbitals if tStoreSpinOrbs is false.

                inquire(file='OneRDM_POPS', exist=exists_one)
                if (exists_one) then
                    RDM_unit = get_free_unit()
                    open(RDM_unit, file='OneRDM_POPS', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i, j, Temp_RDM_Element
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading OneRDM_POPS")
                        if (FileEnd .lt. 0) exit

                        NatOrbMat(SymLabelListInv_rot(i), SymLabelListInv_rot(j)) = Temp_RDM_Element
                    end do
                    close(RDM_unit)
                else
                    call stop_all(t_r, "Attempting to read in the OneRDM, but the OneRDM_POPS file does not exist.")
                end if                                    

            else

                write(6,'(1X,"Reading in the 2-RDMs")')

                ! The TwoRDMs will have been printed exactly as they were.
                ! Without having been made hermitian, without being
                ! normalised, and in spatial orbitals. 

                ! Only read in the 2-RDMs (the 1-RDM becomes redundant).
                inquire(file='TwoRDM_POPS_aaaa', exist=exists_aaaa)
                inquire(file='TwoRDM_POPS_abab', exist=exists_abab)
                inquire(file='TwoRDM_POPS_abba', exist=exists_abba)

                if (tOpenShell)then
                    inquire(file='TwoRDM_POPS_bbbb', exist=exists_bbbb)
                    inquire(file='TwoRDM_POPS_baba', exist=exists_baba)
                    inquire(file='TwoRDM_POPS_baab', exist=exists_baab)
                end if

                if (exists_aaaa .and. exists_abab .and. exists_abba) then
                    ! All TOREAD RDM files are present - read in.
                    RDM_unit = get_free_unit()
                    open(RDM_unit, file='TwoRDM_POPS_aaaa', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading TwoRDM_POPS_aaaa")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        rdm%aaaa_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit, file='TwoRDM_POPS_abab', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading TwoRDM_POPS_abab")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-1) * j ) / 2 ) + i
                        Ind2 = ( ( (b-1) * b ) / 2 ) + a
                        rdm%abab_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit, file='TwoRDM_POPS_abba', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit,iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading TwoRDM_POPS_abba")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        rdm%abba_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                else
                    write(6,*) 'exists_aaaa', exists_aaaa
                    write(6,*) 'exists_abab', exists_abab
                    write(6,*) 'exists_abba', exists_abba
                    call neci_flush(6)
                    call stop_all(t_r,"Attempting to read in the RDMs, &
                                    &but at least one of the TwoRDM_a***_TOREAD files are missing.")
                end if

                if (tOpenShell .and. exists_bbbb .and. exists_baba .and. exists_baab) then
                    ! All TOREAD RDM files for open shell RDMs are present - read in.
                    RDM_unit = get_free_unit()
                    open(RDM_unit, file='TwoRDM_POPS_bbbb', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat=FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading TwoRDM_POPS_bbbb")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        rdm%bbbb_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit, file='TwoRDM_POPS_baba', status='old', form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat = FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd.gt.0) call stop_all(t_r, "Error reading TwoRDM_POPS_baba")
                        if (FileEnd.lt.0) exit

                        Ind1 = ( ( (j-1) * j ) / 2 ) + i
                        Ind2 = ( ( (b-1) * b ) / 2 ) + a
                        rdm%baba_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                    open(RDM_unit,file='TwoRDM_POPS_baab',status='old',form='unformatted')
                    do while (.true.)
                        read(RDM_unit, iostat = FileEnd) i, j, a, b, Temp_RDM_Element 
                        if (FileEnd .gt. 0) call stop_all(t_r, "Error reading TwoRDM_POPS_baab")
                        if (FileEnd .lt. 0) exit

                        Ind1 = ( ( (j-2) * (j-1) ) / 2 ) + i
                        Ind2 = ( ( (b-2) * (b-1) ) / 2 ) + a
                        rdm%baab_full(Ind1,Ind2) = Temp_RDM_Element
                    end do
                    close(RDM_unit)

                else if (tOpenShell) then
                    write(6,*) 'exists_bbbb', exists_bbbb
                    write(6,*) 'exists_baba', exists_baba
                    write(6,*) 'exists_baab', exists_baab
                    call neci_flush(6)
                    call stop_all(t_r, "Attempting to read in the RDMs, &
                                  &but at least one of the TwoRDM_b***_TOREAD files are missing.")
                end if

            end if
        end if

        ! Calculate the energy for the matrices read in (if we're calculating more
        ! than the 1-RDM).
        if (tCalc_RDMEnergy) call rdm_output_wrapper(rdm, Norm_2RDM)

        ! Continue calculating the RDMs from the first iteration when the popsfiles
        ! (and RDMs) are read in. This overwrites the iteration number put in the input.
        IterRDMonFly = 0

    end subroutine Read_In_RDMs

    subroutine SetUpSymLabels_RDM() 

        ! This routine just sets up the symmetry labels so that the orbitals
        ! are ordered according to symmetry (all beta then all alpha if spin orbs).

        use rdm_data, only: tOpenShell
        use RotateOrbsData, only: SymLabelList2_rot, SymLabelCounts2_rot, SymLabelListInv_rot
        use RotateOrbsData, only: NoOrbs, SpatOrbs, NoSymLabelCounts
        use sort_mod, only: sort
        use SystemData, only: G1, BRR, lNoSymmetry, tFixLz, iMaxLz
        use UMatCache, only: gtID
        use util_mod, only: LogMemAlloc, LogMemDealloc

        integer, allocatable :: SymOrbs_rot(:)
        integer :: LabOrbsTag, SymOrbs_rotTag, ierr, i, j, SpatSym, LzSym 
        integer :: lo, hi, Symi, SymCurr, Symi2, SymCurr2
        character(len=*), parameter :: t_r = 'SetUpSymLabels_RDM'

        ! This is only allocated temporarily to be used to order the orbitals by.
        allocate(SymOrbs_rot(NoOrbs), stat=ierr)
        call LogMemAlloc('SymOrbs_rot', NoOrbs, 4, t_r, SymOrbs_rotTag, ierr)
        if (ierr .ne. 0) call stop_all(t_r, "Mem allocation for SymOrbs_rot failed.")

        ! Now we want to put the spatial orbital index, followed by the symmetry.
        SymLabelList2_rot(:) = 0
        SymOrbs_rot(:) = 0

        ! *** STEP 1 *** Fill SymLabelList2_rot.
        ! Find the orbitals and order them in terms of symmetry.
        do i=1, SpatOrbs
            if (tOpenShell) then
                ! For open shell systems, all alpha are followed by all beta.
                SymLabelList2_rot(i) = BRR(2*i)
                SymLabelList2_rot(i+SpatOrbs) = BRR((2*i)-1)

                if (tFixLz) then
                    SpatSym = int(G1(BRR(2*i))%sym%S)
                    LzSym = int(G1(BRR(2*i))%Ml)
                    SymOrbs_rot(i) = ( SpatSym * ((2 * iMaxLz) + 1) ) + ( LzSym + iMaxLz )

                    SpatSym = int(G1(BRR((2*i)-1))%sym%S)
                    LzSym = int(G1(BRR((2*i)-1))%Ml)
                    SymOrbs_rot(i+SpatOrbs) = ( SpatSym * ((2 * iMaxLz) + 1) ) + ( LzSym + iMaxLz )
                else
                    SymOrbs_rot(i) = int(G1(BRR(2*i))%sym%S) 
                    SymOrbs_rot(i+SpatOrbs) = int(G1(BRR((2*i)-1))%sym%S) 
                end if
            else
                SymLabelList2_rot(i) = gtID(BRR(2*i))
                if (tFixLz) then
                    SpatSym = int(G1(BRR(2*i))%sym%S)
                    LzSym = int(G1(BRR(2*i))%Ml)
                    SymOrbs_rot(i) = ( SpatSym * ((2 * iMaxLz) + 1) ) + ( LzSym + iMaxLz )
                else
                    SymOrbs_rot(i) = int(G1(BRR(2*i))%sym%S)
                end if
                ! Orbital BRR(2*i) for i = 1 will be the beta orbital with the 
                ! second lowest energy - want the spatial orbital index to go with this.
                ! G1 also in spin orbitals - get symmetry of this beta orbital, will 
                ! be the same as the spatial orbital.
            end if
        end do

        call sort (SymOrbs_rot(1:SpatOrbs), SymLabelList2_rot(1:SpatOrbs))
        ! Sorts SymLabelList2_rot according to the order of SymOrbs_rot
        ! (i.e. in terms of symmetry).
        if (tOpenShell) &
            call sort (SymOrbs_rot(SpatOrbs+1:nBasis), SymLabelList2_rot(SpatOrbs+1:nBasis))
            ! Also do this for the beta set if spin orbitals.

        ! *** STEP 2 *** Fill SymLabelCounts2_rot_rot. This is like
        ! SymLabelCounts2_rot, but has a different ordering - BEWARE.
        ! SymLabelCounts(1,:) contains the position in SymLabelList2_rot
        ! where the symmetry index starts, SymLabelCounts(2,:) contains the
        ! number of orbitals in that symmetry index. Again if spin orbs, all
        ! alpha are followed by all beta - i.e. first 8 refer to alpha, second
        ! 8 to beta.

        if (lNoSymmetry) then
            ! If we are ignoring symmetry, all orbitals essentially have
            ! symmetry 0.
            SymLabelCounts2_rot(1,1) = 1
            SymLabelCounts2_rot(2,1) = SpatOrbs
            if (tOpenShell) then
                SymLabelCounts2_rot(1,9) = SpatOrbs+1
                SymLabelCounts2_rot(2,9) = SpatOrbs
            end if
        else 
            ! Otherwise we run through the occupied orbitals, counting the
            ! number with each symmetry (spatial and Lz) and noting where in
            ! SymLabelList2_rot each symmetry block starts.
            SymCurr = 0
            SymLabelCounts2_rot(1,1) = 1
            if (tOpenShell) then
                SymCurr2 = 0
                SymLabelCounts2_rot(1,9) = SpatOrbs + 1
            end if
            do i = 1,SpatOrbs
                if (tOpenShell) then
                    Symi = SymOrbs_rot(i)
                    Symi2 = SymOrbs_rot(i + SpatOrbs)
                else
                    Symi = SymOrbs_rot(i)
                end if
                SymLabelCounts2_rot(2,(Symi+1)) = SymLabelCounts2_rot(2,(Symi+1)) + 1
                if (Symi .ne. SymCurr) then
                    do j = SymCurr + 1, Symi
                        SymLabelCounts2_rot(1,(j+1)) = i
                    end do
                    SymCurr = Symi
                end if
                if (tOpenShell) then
                    SymLabelCounts2_rot(2,(Symi2+9)) = SymLabelCounts2_rot(2,(Symi2+9))+1
                    if (Symi2 .ne. SymCurr2) then
                        do j = SymCurr2 + 1, Symi2
                            SymLabelCounts2_rot(1,(j+9)) = i + SpatOrbs
                        end do
                        SymCurr2 = Symi2
                    end if
                end if
            end do
        end if

        ! Go through each symmetry group, making sure the orbitals are 
        ! ordered lowest to highest within each symmetry.
        do i = 1, NoSymLabelCounts
            if (SymLabelCounts2_rot(2,i) .ne. 0) then
                lo = SymLabelCounts2_rot(1, i)
                hi = lo + SymLabelCounts2_rot(2, i) - 1
                call sort (SymLabelList2_rot (lo:hi))
            end if
        end do

        ! Construct the inverse matrix. While SymLabelList2_rot takes a
        ! position and tells us what orbital is in it, we also might need to
        ! take an orbital and find out what position to put its contribution in.
        do i = 1, NoOrbs
            SymLabelListInv_rot(SymLabelList2_rot(i)) = i
        end do

        ! Deallocate the arrays just used in this routine.
        deallocate(SymOrbs_rot)
        call LogMemDealloc(t_r,SymOrbs_rotTag)

    end subroutine SetUpSymLabels_RDM

    subroutine DeAlloc_Alloc_SpawnedParts()

        ! Routine called when RDM accumulation is turned on, usually midway
        ! through an FCIQMC simulation.

        ! When calculating the RDMs, we need to store the parent from which a
        ! child is spawned along with the children in the spawned array. This
        ! means a slightly larger array is communicated between processors,
        ! which there is no point in doing for the first part of the calculation.
        ! When we start calculating the RDMs this routine is called and the
        ! SpawnedParts array is made larger to accommodate the parents.

        use FciMCData, only: MaxSpawned, SpawnVec, SpawnVec2, SpawnVecTag, SpawnVec2Tag
        use FciMCData, only: SpawnedParts, SpawnedParts2
        use util_mod, only: LogMemAlloc, LogMemDealloc

        integer :: ierr                               
        character(len=*), parameter :: t_r = 'DeAlloc_Alloc_SpawnedParts'
        
        deallocate(SpawnVec)
        call LogMemDealloc(t_r,SpawnVecTag)
        deallocate(SpawnVec2)
        call LogMemDealloc(t_r,SpawnVec2Tag)
 
        allocate(SpawnVec(0:(NIftot+NIfDBO+2), MaxSpawned), stat=ierr)
        call LogMemAlloc('SpawnVec', MaxSpawned*(NIfTot+NIfDBO+3), size_n_int, t_r, SpawnVecTag, ierr)
        allocate(SpawnVec2(0:(NIfTot+NIfDBO+2), MaxSpawned),stat=ierr)
        call LogMemAlloc('SpawnVec2', MaxSpawned*(NIfTot+NIfDBO+3), size_n_int, t_r, SpawnVec2Tag, ierr)

        ! Point at correct spawning arrays
        SpawnedParts => SpawnVec
        SpawnedParts2 => SpawnVec2

        write(6,'(A54,F10.4,A4,F10.4,A13)') 'Memory requirement for spawned arrays increased from ',&
                                        real(((NIfTot+1)*MaxSpawned*2*size_n_int),dp)/1048576.0_dp,' to ',&
                                        real(((NIfTot+NIfDBO+3)*MaxSpawned*2*size_n_int),dp)/1048576.0_dp, ' Mb/Processor'

    end subroutine DeAlloc_Alloc_SpawnedParts

    ! Routines called at the end of a simulation.

    subroutine FinaliseRDMs(rdms)

        ! This routine performs some finalisation, including summing each of
        ! the individual matrices from each processor, and calling the
        ! diagonalisation routines if we want to get the occupation numbers.

        use FciMCData, only: tFinalRDMEnergy
        use LoggingData, only: tBrokenSymNOs, occ_numb_diff, RDMExcitLevel, tExplicitAllRDM
        use LoggingData, only: tPrint1RDM, tDiagRDM, tDumpForcesInfo, tDipoles
        use Parallel_neci, only: iProcIndex, MPIBarrier
        use rdm_data, only: rdm_t, tRotatedNos, FinaliseRDMs_Time
        use rdm_estimators, only: Calc_Lagrangian_from_RDM, convert_mats_Molpforces
        use rdm_estimators, only: rdm_output_wrapper, CalcDipoles
        use rdm_nat_orbs, only: find_nat_orb_occ_numbers, BrokenSymNo
        use util_mod, only: set_timer, halt_timer

        type(rdm_t), intent(inout) :: rdms(:)

        integer :: i, error
        real(dp) :: Norm_2RDM, Norm_2RDM_Inst
        real(dp) :: Norm_1RDM, Trace_1RDM, SumN_Rho_ii
        character(len=*), parameter :: t_r = 'FinaliseRDMs'

        call set_timer(FinaliseRDMs_Time)

        if (tExplicitAllRDM) then
            write(6,'(/,"**** RDMs CALCULATED EXPLICITLY ****",1X,/)')
        else
            write(6,'(/,"**** RDMs CALCULATED STOCHASTICALLY ****",1X,/)')
        end if

        ! Combine the 1- or 2-RDM from all processors, etc.

        do i = 1, size(rdms)

            if (RDMExcitLevel .eq. 1) then
                call Finalise_1e_RDM(rdms(i), Norm_1RDM)
            else
                ! We always want to calculate one final RDM energy, whether or not we're 
                ! calculating the energy throughout the calculation.
                ! Unless of course, only the 1-RDM is being calculated.

                ! Calculate the energy one last time - and write out everything we need.
                tFinalRDMEnergy = .true.

                ! 1-RDM is constructed here (in calc_1RDM_energy).
                call rdm_output_wrapper(rdms(i), Norm_2RDM)

                if (tPrint1RDM) then
                    call Finalise_1e_RDM(rdms(i), Norm_1RDM)
                else if (tDiagRDM .and. (iProcIndex .eq. 0)) then
                    call calc_1e_norms(rdms(i), Trace_1RDM, Norm_1RDM, SumN_Rho_ii)
                    write(6,'(/,1X,"SUM OF 1-RDM(i,i) FOR THE N LOWEST ENERGY HF ORBITALS:",1X,F20.13)') SumN_Rho_ii
                end if

                if (tDumpForcesInfo) then
                    if (.not. tPrint1RDM) call Finalise_1e_RDM(rdms(i), Norm_1RDM)
                    call Calc_Lagrangian_from_RDM(rdms(i), Norm_1RDM, Norm_2RDM)
                    call convert_mats_Molpforces(rdms(i), Norm_1RDM, Norm_2RDM)
                end if

            end if

            call MPIBarrier(error)

            ! Call the routines from NatOrbs that diagonalise the one electron
            ! reduced density matrix.
            tRotatedNOs = .false. ! Needed for BrokenSymNo routine
            if (tDiagRDM) call find_nat_orb_occ_numbers(rdms(i))

            ! This is where we would likely call any further calculations of
            ! forces, etc.
            if (tDipoles) then
                if (.not. tPrint1RDM) call Finalise_1e_RDM(rdms(i), Norm_1RDM)
                call CalcDipoles(Norm_1RDM)
            end if

            ! After all the NO calculations are finished we'd like to do another
            ! rotation to obtain symmetry-broken natural orbitals
            if (tBrokenSymNOs) then
                call BrokenSymNO(occ_numb_diff)
            end if

        end do

        call halt_timer(FinaliseRDMs_Time)
    
    end subroutine FinaliseRDMs

    subroutine Finalise_1e_RDM(rdm, Norm_1RDM) 

        ! This routine takes the 1-RDM (NatOrbMat), normalises it, makes it 
        ! hermitian if required, and prints out the versions we're interested
        ! in. This is only ever called at the very end of a calculation.

        use LoggingData, only: twrite_RDMs_to_read, twrite_normalised_RDMs, tForceCauchySchwarz
        use LoggingData, only: RDMExcitLevel
        use NatOrbsMod, only: NatOrbMat
        use Parallel_neci, only: iProcIndex, MPISumAll
        use rdm_data, only: rdm_t
        use RotateOrbsData, only: NoOrbs

        type(rdm_t), intent(inout) :: rdm
        real(dp), intent(out) :: Norm_1RDM
                             
        integer :: i, ierr
        real(dp) :: Trace_1RDM, SumN_Rho_ii
        real(dp), allocatable :: AllNode_NatOrbMat(:,:)

        Norm_1RDM = 0.0_dp

        if (RDMExcitLevel .eq. 1) then

            allocate(AllNode_NatOrbMat(NoOrbs, NoOrbs), stat=ierr)
            
            call MPISumAll(NatOrbMat, AllNode_NatOrbMat)
            NatOrbMat = AllNode_NatOrbMat
            
            deallocate(AllNode_NatOrbMat)

        end if

        if (iProcIndex .eq. 0) then 

            ! Find the normalisation.
            call calc_1e_norms(rdm, Trace_1RDM, Norm_1RDM, SumN_Rho_ii)

            ! Write out the unnormalised, non-hermitian OneRDM_POPS.
            if (twrite_RDMs_to_read) call Write_out_1RDM(Norm_1RDM, .false.)

            ! Enforce the hermiticity condition.  If the RDMExcitLevel is not 1, the 
            ! 1-RDM has been constructed from the hermitian 2-RDM, so this will not 
            ! be necessary.
            ! The HF_Ref and HF_S_D_Ref cases are not hermitian by definition.
            if (RDMExcitLevel .eq. 1) then
                call make_1e_rdm_hermitian(Norm_1RDM)
                
                if (tForceCauchySchwarz)then
                    call Force_Cauchy_Schwarz(Norm_1RDM)
                end if

            end if
            
            ! Write out the final, normalised, hermitian OneRDM.                
            if (twrite_normalised_RDMs) call Write_out_1RDM(Norm_1RDM, .true.)

            write(6,'(1X,"SUM OF 1-RDM(i,i) FOR THE N LOWEST ENERGY HF ORBITALS:",1X,F20.13)') SumN_Rho_ii

        end if

    end subroutine Finalise_1e_RDM

    subroutine calc_1e_norms(rdm, Trace_1RDM, Norm_1RDM, SumN_Rho_ii)

        ! We want to 'normalise' the reduced density matrices. These are not
        ! even close to being normalised at the moment, because of the way
        ! they are calculated on the fly. They should be calculated from a
        ! normalised wavefunction. But we know that the trace of the one
        ! electron reduced density matrix must be equal to the number of the
        ! electrons. We can use this to find the factor we must divide the
        ! 1-RDM through by.

        use FciMCData, only: HFDet_True
        use LoggingData, only: tDiagRDM
        use NatOrbsMod, only: NatOrbMat
        use rdm_data, only: rdm_t, tOpenShell
        use RotateOrbsData, only: SymLabelListInv_rot, NoOrbs
        use SystemData, only: BRR
        use UMatCache, only: gtID

        type(rdm_t), intent(inout) :: rdm
        real(dp), intent(out) :: Trace_1RDM, Norm_1RDM, SumN_Rho_ii

        integer :: i, HFDet_ID, BRR_ID

        Trace_1RDM = 0.0_dp
        Norm_1RDM = 0.0_dp

        do i = 1, NoOrbs
            Trace_1RDM = Trace_1RDM + NatOrbMat(i,i)
        end do

        Norm_1RDM = ( real(NEl, dp) / Trace_1RDM )
        
        ! Need to multiply each element of the 1 electron reduced density matrices 
        ! by NEl / Trace_1RDM,
        ! and then add it's contribution to the energy.
        
        ! Want to sum the diagonal elements of the 1-RDM for the HF orbitals.
        ! Given the HF orbitals, SymLabelListInv_rot tells us their position
        ! in the 1-RDM.
        SumN_Rho_ii = 0.0_dp

        do i = 1, NoOrbs

            ! Rho_ii is the diagonal elements of the 1-RDM. We want this
            ! ordered according to the energy of the orbitals. Brr has the
            ! orbital numbers in order of energy... i.e Brr(2) = the orbital
            ! index with the second lowest energy. Brr is always in spin
            ! orbitals. i gives the energy level, BRR gives the orbital,
            ! SymLabelListInv_rot gives the position of  this orbital in
            ! NatOrbMat.

            if (tDiagRDM) then
                if (tOpenShell) then
                    rdm%Rho_ii(i) = NatOrbMat(SymLabelListInv_rot(BRR(i)),SymLabelListInv_rot(BRR(i))) * Norm_1RDM
                else
                    BRR_ID = gtID(BRR(2*i))
                    rdm%Rho_ii(i) = NatOrbMat(SymLabelListInv_rot(BRR_ID),SymLabelListInv_rot(BRR_ID)) * Norm_1RDM
                end if
            end if
    
            if (i.le.NEl) then
                if (tOpenShell) then
                    SumN_Rho_ii = SumN_Rho_ii + &
                            ( NatOrbMat(SymLabelListInv_rot(HFDet_True(i)),SymLabelListInv_rot(HFDet_True(i))) &
                                * Norm_1RDM )
                else
                    HFDet_ID = gtID(HFDet_True(i))
                    SumN_Rho_ii = SumN_Rho_ii + &
                            ( NatOrbMat(SymLabelListInv_rot(HFDet_ID),SymLabelListInv_rot(HFDet_ID)) &
                                * Norm_1RDM ) / 2.0_dp
                end if
            end if
        end do

    end subroutine calc_1e_norms

    subroutine make_1e_rdm_hermitian(Norm_1RDM)

        ! Simply average the 1-RDM(i,j) and 1-RDM(j,i) elements which should
        ! be equal in a perfect world.

        use NatOrbsMod, only: NatOrbMat
        use RotateOrbsData, only: SymLabelListInv_rot, NoOrbs

        real(dp), intent(in) :: Norm_1RDM
        real(dp) :: Max_Error_Hermiticity, Sum_Error_Hermiticity 
        integer :: i, j
        real(dp) :: Temp

        Max_Error_Hermiticity = 0.0_dp
        Sum_Error_Hermiticity = 0.0_dp
        do i = 1, NoOrbs
            do j = i, NoOrbs
                if ((abs((NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))*Norm_1RDM) - &
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
            end do
        end do

        ! Output the hermiticity errors.
        write(6,'(1X,"MAX ABS ERROR IN 1RDM HERMITICITY",F20.13)') Max_Error_Hermiticity
        write(6,'(1X,"MAX ABS ERROR IN 1RDM HERMITICITY",F20.13)') Sum_Error_Hermiticity

    end subroutine make_1e_rdm_hermitian

    subroutine Force_Cauchy_Schwarz(Norm_1RDM)

        use NatOrbsMod, only: NatOrbMat
        use RotateOrbsData, only: SymLabelListInv_rot

        real(dp), intent(in) :: Norm_1RDM
        integer :: i, j
        real(dp) :: UpperBound

        write(6,'("Ensuring that Cauchy--Schwarz inequality holds.")')

        do i = 1, nBasis
            do j = 1, nBasis

                UpperBound = sqrt(NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(i))&
                    *NatOrbMat(SymLabelListInv_rot(j),SymLabelListInv_rot(j)))

                if (abs(NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j))) .gt. UpperBound)then

                    if (NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) .lt. 0.0_dp)then
                        NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = -UpperBound
                    else if (NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) .gt. 0.0_dp)then
                        NatOrbMat(SymLabelListInv_rot(i),SymLabelListInv_rot(j)) = UpperBound
                    end if

                    write(6,'("Changing element:")') i, j
                else
                    cycle
                end if
            end do
        end do

    end subroutine Force_Cauchy_Schwarz

    subroutine Write_out_1RDM(Norm_1RDM, tNormalise)

        ! This routine writes out the OneRDM. If tNormalise is true, we are
        ! printing the normalised, hermitian matrix. Otherwise, Norm_1RDM is
        ! ignored and we print both 1-RDM(i,j) and 1-RDM(j,i) (in binary) 
        ! for the OneRDM_POPS file to be read in in a restart calculation.

        use NatOrbsMod, only: NatOrbMat
        use rdm_data, only: tOpenShell
        use RotateOrbsData, only: SymLabelListInv_rot
        use UMatCache, only: gtID
        use util_mod, only: get_free_unit

        real(dp), intent(in) :: Norm_1RDM
        logical, intent(in) :: tNormalise
        integer :: i, j, iSpat, jSpat
        integer :: OneRDM_unit

        if (tNormalise) then
            ! Haven't got the capabilities to produce multiple 1-RDMs yet.
            write(6,'(1X,"Writing out the *normalised* 1 electron density matrix to file")')
            call neci_flush(6)
            OneRDM_unit = get_free_unit()
            open(OneRDM_unit,file='OneRDM',status='unknown')
        else
            ! Only every write out 1 of these at the moment.
            write(6,'(1X,"Writing out the *unnormalised* 1 electron density matrix to file for reading in")')
            call neci_flush(6)
            OneRDM_unit = get_free_unit()
            open(OneRDM_unit, file='OneRDM_POPS', status='unknown', form='unformatted')
        end if

        ! Currently always printing 1-RDM in spin orbitals.
        do i = 1, nBasis
            do j = 1, nBasis
                if (tOpenShell) then
                    if (NatOrbMat(SymLabelListInv_rot(i), SymLabelListInv_rot(j)) .ne. 0.0_dp) then 
                        if (tNormalise .and. (i .le. j)) then
                            write(OneRDM_unit,"(2I6,G25.17)") i, j, &
                                NatOrbMat(SymLabelListInv_rot(i), SymLabelListInv_rot(j)) * Norm_1RDM
                        else if (.not. tNormalise) then
                            ! For the pops, we haven't made the 1-RDM hermitian yet, 
                            ! so print both the 1-RDM(i,j) and 1-RDM(j,i) elements.
                            ! This is written in binary.
                            write(OneRDM_unit) i, j, NatOrbMat(SymLabelListInv_rot(i), SymLabelListInv_rot(j))
                        end if
                    end if
                else
                    iSpat = gtID(i)
                    jSpat = gtID(j)
                    if (NatOrbMat(SymLabelListInv_rot(iSpat), SymLabelListInv_rot(jSpat)) .ne. 0.0_dp) then
                        if (tNormalise .and. (i .le. j)) then
                            if (((mod(i,2).eq.0) .and. (mod(j,2) .eq. 0)) .or. &
                                ((mod(i,2).ne.0) .and. (mod(j,2) .ne. 0))) then
                                write(OneRDM_unit,"(2I6,G25.17)") i,j, & 
                                    ( NatOrbMat(SymLabelListInv_rot(iSpat),SymLabelListInv_rot(jSpat)) &
                                                                    * Norm_1RDM ) / 2.0_dp
                            end if
                        else if (.not. tNormalise) then
                            ! The popsfile can be printed in spatial orbitals.
                            if ((mod(i,2) .eq. 0) .and. (mod(j,2) .eq. 0)) then
                                write(OneRDM_unit) iSpat, jSpat, & 
                                    NatOrbMat(SymLabelListInv_rot(iSpat), SymLabelListInv_rot(jSpat)) 
                            end if
                        end if
                    end if
                end if
            end do
        end do

        close(OneRDM_unit)

    end subroutine Write_out_1RDM

    subroutine DeallocateRDM()

        ! This routine just deallocates the arrays allocated in InitRDMs.
        ! If the NECI calculation softexits before the RDMs start to fill,
        ! this is all that is called at the end.

        use FciMCData, only: Spawned_Parents, Spawned_Parents_Index
        use FciMCData, only: Spawned_ParentsTag, Spawned_Parents_IndexTag
        use LoggingData, only: RDMExcitLevel, tExplicitAllRDM
        use NatOrbsMod, only: NatOrbMat, NatOrbMatTag, Evalues, EvaluesTag
        use rdm_data, only: rdms, Sing_ExcDjs, Doub_ExcDjs
        use rdm_data, only: Sing_ExcDjs2, Doub_ExcDjs2, Sing_ExcDjsTag, Doub_ExcDjsTag
        use rdm_data, only: Sing_ExcDjs2Tag, Doub_ExcDjs2Tag
        use rdm_data, only: Sing_InitExcSlots, Doub_InitExcSlots, Sing_ExcList, Doub_ExcList
        use RotateOrbsData, only: SymLabelCounts2_rot,SymLabelList2_rot, SymLabelListInv_rot
        use RotateOrbsData, only: SymLabelCounts2_rotTag, SymLabelList2_rotTag
        use RotateOrbsData, only: SymLabelListInv_rotTag
        use RotateOrbsMod, only: FourIndInts, FourIndIntsTag
        use util_mod, only: LogMemDealloc

        integer :: i
        character(len=*), parameter :: t_r = 'DeallocateRDM'

        if (tExplicitAllRDM) then

            ! This array contains the initial positions of the single
            ! excitations for each processor.
            deallocate(Sing_InitExcSlots)
 
            ! This array contains the current position of the single
            ! excitations as they're added.
            deallocate(Sing_ExcList)

            ! This array actually contains the single excitations in blocks of
            ! the processor they will be sent to.        
            deallocate(Sing_ExcDjs)
            call LogMemDeAlloc(t_r,Sing_ExcDjsTag)

            deallocate(Sing_ExcDjs2)
            call LogMemDeAlloc(t_r,Sing_ExcDjs2Tag)


            if (RDMExcitLevel .ne. 1) then
                ! This array contains the initial positions of the
                ! single excitations for each processor.
                deallocate(Doub_InitExcSlots)
 
                ! This array contains the current position of the single
                ! excitations as they're added.
                deallocate(Doub_ExcList)

                ! This array actually contains the single excitations in
                ! blocks of the  processor they will be sent to.        
                deallocate(Doub_ExcDjs)
                call LogMemDeAlloc(t_r,Doub_ExcDjsTag)
     
                deallocate(Doub_ExcDjs2)
                call LogMemDeAlloc(t_r,Doub_ExcDjs2Tag)
            end if

        else

            if (allocated(Spawned_Parents)) then
                deallocate(Spawned_Parents)
                call LogMemDeAlloc(t_r,Spawned_ParentsTag)
            end if

            if (allocated(Spawned_Parents_Index)) then
                deallocate(Spawned_Parents_Index)
                call LogMemDeAlloc(t_r,Spawned_Parents_IndexTag)
            end if

        end if

        if (allocated(NatOrbMat)) then
            deallocate(NatOrbMat)
            call LogMemDeAlloc(t_r,NatOrbMatTag)
        end if

        if (allocated(Evalues)) then
            deallocate(Evalues)
            call LogMemDeAlloc(t_r,EvaluesTag)
        end if

        if (allocated(FourIndInts)) then
            deallocate(FourIndInts)
            call LogMemDeAlloc(t_r,FourIndIntsTag)
        end if

        if (allocated(SymLabelCounts2_rot)) then
            deallocate(SymLabelCounts2_rot)
            call LogMemDeAlloc(t_r,SymLabelCounts2_rotTag)
        end if

        if (allocated(SymLabelList2_rot)) then
            deallocate(SymLabelList2_rot)
            call LogMemDeAlloc(t_r,SymLabelList2_rotTag)
        end if

        if (allocated(SymLabelListInv_rot)) then
            deallocate(SymLabelListInv_rot)
            call LogMemDeAlloc(t_r,SymLabelListInv_rotTag)
        end if

        do i = 1, size(rdms)
            if (allocated(rdms(i)%Rho_ii)) then
                deallocate(rdms(i)%Rho_ii)
                call LogMemDeAlloc(t_r,rdms(i)%Rho_iiTag)
            end if

            if (associated(rdms(i)%aaaa_inst)) then
                deallocate(rdms(i)%aaaa_inst)
                call LogMemDeAlloc(t_r,rdms(i)%aaaa_instTag)
            end if

            if (associated(rdms(i)%abab_inst)) then
                deallocate(rdms(i)%abab_inst)
                call LogMemDeAlloc(t_r,rdms(i)%abab_instTag)
            end if

            if (associated(rdms(i)%abba_inst)) then
                deallocate(rdms(i)%abba_inst)
                call LogMemDeAlloc(t_r,rdms(i)%abba_instTag)
            end if

            if (associated(rdms(i)%bbbb_inst)) then
                deallocate(rdms(i)%bbbb_inst)
                call LogMemDeAlloc(t_r,rdms(i)%bbbb_instTag)
            end if

            if (associated(rdms(i)%baba_inst)) then
                deallocate(rdms(i)%baba_inst)
                call LogMemDeAlloc(t_r,rdms(i)%baba_instTag)
            end if

            if (associated(rdms(i)%baab_inst)) then
                deallocate(rdms(i)%baab_inst)
                call LogMemDeAlloc(t_r,rdms(i)%baab_instTag)
            end if

            if (associated(rdms(i)%aaaa_full)) then
                deallocate(rdms(i)%aaaa_full)
                call LogMemDeAlloc(t_r,rdms(i)%aaaa_fullTag)
            end if

            if (associated(rdms(i)%abab_full)) then
                deallocate(rdms(i)%abab_full)
                call LogMemDeAlloc(t_r,rdms(i)%abab_fullTag)
            end if

            if (associated(rdms(i)%abba_full)) then
                deallocate(rdms(i)%abba_full)
                call LogMemDeAlloc(t_r,rdms(i)%abba_fullTag)
            end if

            if (associated(rdms(i)%bbbb_full)) then
                deallocate(rdms(i)%bbbb_full)
                call LogMemDeAlloc(t_r,rdms(i)%bbbb_fullTag)
            end if

            if (associated(rdms(i)%baba_full)) then
                deallocate(rdms(i)%baba_full)
                call LogMemDeAlloc(t_r,rdms(i)%baba_fullTag)
            end if

            if (associated(rdms(i)%baab_full)) then
                deallocate(rdms(i)%baab_full)
                call LogMemDeAlloc(t_r,rdms(i)%baab_fullTag)
            end if
        end do

    end subroutine DeallocateRDM


    ! Some general routines used during the main simulation.


    subroutine extract_bit_rep_avsign_no_rdm(iLutnI, j, nI, SignI, FlagsI, IterRDMStartI, AvSignI, Store)

        ! This is just the standard extract_bit_rep routine for when we're not
        ! calculating the RDMs.    

        use bit_reps, only: extract_bit_rep
        use FciMCData, only: excit_gen_store_type

        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(in) :: j
        integer, intent(out) :: nI(nel), FlagsI
        real(dp), dimension(lenof_sign), intent(out) :: SignI
        real(dp), dimension(lenof_sign), intent(out) :: IterRDMStartI, AvSignI
        type(excit_gen_store_type), intent(inout), optional :: Store
        
        ! This extracts everything.
        call extract_bit_rep (iLutnI, nI, SignI, FlagsI, store)

        IterRDMStartI(:) = 0.0_dp
        AvSignI(:) = 0.0_dp

    end subroutine extract_bit_rep_avsign_no_rdm

    subroutine extract_bit_rep_avsign_norm(iLutnI, j, nI, SignI, FlagsI, IterRDMStartI, AvSignI, Store)

        ! The following extract_bit_rep_avsign routine extracts the bit
        ! representation of the current determinant, and calculates the average
        ! sign since this determinant became occupied. 

        ! In double run, we have to be particularly careful -- we need to start
        ! a new average when the determinant becomes newly occupied or
        ! unoccupied in either population (see CMO thesis). Additionally, we're
        ! also setting it up so that averages get restarted whenever we
        ! calculate the energy which saves a lot of faffing about, and storage
        ! of an extra set of RDMs, and is still unbiased. This is called for
        ! each determinant in the occupied list at the beginning of its FCIQMC
        ! cycle. It is used if we're calculating the RDMs with or without HPHF. 

        ! Input:    iLutnI (bit rep of current determinant).
        !           j - Which element in the CurrentDets array are we considering?
        ! Output:   nI, SignI, FlagsI after extract.                                              
        !           IterRDMStartI - new iteration the determinant became occupied (as a real).
        !           AvSignI - the new average walker population during this time (also real).

        use bit_reps, only: extract_bit_rep
        use FciMCData, only: PreviousCycles, Iter, IterRDMStart, excit_gen_store_type
        use global_det_data, only: get_iter_occ, get_av_sgn
        use LoggingData, only: RDMEnergyIter

        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(out) :: nI(nel), FlagsI
        integer, intent(in) :: j
        real(dp), dimension(lenof_sign), intent(out) :: SignI
        real(dp), dimension(lenof_sign), intent(out) :: IterRDMStartI, AvSignI
        type(excit_gen_store_type), intent(inout), optional :: Store
        integer :: part_ind

        ! This is the iteration from which this determinant has been occupied.
        IterRDMStartI(1:lenof_sign) = get_iter_occ(j)
        
        ! This extracts everything.
        call extract_bit_rep (iLutnI, nI, SignI, FlagsI)
            
        if (((Iter+PreviousCycles-IterRDMStart) .gt. 0) .and. &
            & (mod(((Iter-1)+PreviousCycles - IterRDMStart + 1), RDMEnergyIter) .eq. 0)) then 

            ! The previous iteration was one where we added in diagonal elements
            ! To keep things unbiased, we need to set up a new averaging block now.
            ! NB: if doing single run cutoff, note that doing things this way is now
            ! NOT the same as the technique described in CMO (and DMC's) thesis.
            ! Would expect diagonal elements to be slightly worse quality, improving
            ! as one calculates the RDM energy less frequently.  As this method is
            ! biased anyway, I'm not going to lose sleep over it.
            do part_ind = 1, lenof_sign
                AvSignI(part_ind) = SignI(part_ind)
                IterRDMStartI(part_ind) = real(Iter + PreviousCycles,dp)
            end do
        else
            ! Now let's consider other instances in which we need to start a new block:
            if (inum_runs .eq. 2) then
#if defined(__DOUBLERUN) || defined(__PROG_NUMRUNS) || defined(__CMPLX)
                if ((SignI(1) .eq. 0) .and. (IterRDMStartI(1) .ne. 0)) then
                    ! The population has just gone to zero on population 1.
                    ! Therefore, we need to start a new averaging block.
                    AvSignI(1) = 0
                    IterRDMStartI(1) = 0
                    AvSignI(2) = SignI(2)
                    IterRDMStartI(2) = real(Iter + PreviousCycles,dp)

                else if ((SignI(2) .eq. 0) .and. (IterRDMStartI(2) .ne. 0)) then
                    ! The population has just gone to zero on population 2.
                    ! Therefore, we need to start a new averaging block.
                    AvSignI(2) = 0
                    IterRDMStartI(2) = 0
                    AvSignI(1) = SignI(1)
                    IterRDMStartI(1) = real(Iter + PreviousCycles,dp)

                else if ((SignI(1) .ne. 0) .and. (IterRDMStartI(1) .eq. 0)) then
                    ! Population 1 has just become occupied.
                    IterRDMStartI(1) = real(Iter + PreviousCycles,dp)
                    IterRDMStartI(2) = real(Iter + PreviousCycles,dp)
                    AvSignI(1) = SignI(1)
                    AvSignI(2) = SignI(2)
                    if (SignI(2) .eq. 0) IterRDMStartI(2) = 0

                else if ((SignI(2) .ne. 0) .and. (IterRDMStartI(2) .eq. 0)) then
                    ! Population 2 has just become occupied.
                    IterRDMStartI(1) = real(Iter + PreviousCycles,dp)
                    IterRDMStartI(2) = real(Iter + PreviousCycles,dp)
                    AvSignI(1) = SignI(1)
                    AvSignI(2) = SignI(2)
                    if (SignI(1) .eq. 0) IterRDMStartI(1) = 0

                else
                    ! Nothing unusual has happened so update both populations
                    ! as normal.
                    do part_ind = 1, lenof_sign
                        ! Update the average population.
                        AvSignI(part_ind) = ( ((real(Iter+PreviousCycles,dp) - IterRDMStartI(part_ind)) * get_av_sgn(j, part_ind))&
                            + SignI(part_ind) ) / ( real(Iter+PreviousCycles,dp) - IterRDMStartI(part_ind) + 1.0_dp )
                    end do
                end if
#endif
            else
                do part_ind = 1, lenof_sign
                    ! If there is nothing stored there yet, the first iteration
                    ! the determinant became occupied is this one.
                    if (IterRDMStartI(part_ind) .eq. 0.0_dp) IterRDMStartI(part_ind) = real(Iter+PreviousCycles, dp)

                    ! Update the average population. This just comes out as the
                    ! current population (SignI) if this is the first  time the
                    ! determinant has become occupied.
                    AvSignI(part_ind) = ( ((real(Iter+PreviousCycles,dp) - IterRDMStartI(part_ind)) * get_av_sgn(j,part_ind)) &
                                    + SignI(part_ind) ) / ( real(Iter+PreviousCycles,dp) - IterRDMStartI(part_ind) + 1.0_dp )
                end do
            end if
        end if

    end subroutine extract_bit_rep_avsign_norm

    subroutine calc_rdmbiasfac(p_spawn_rdmfac, p_gen, SignCurr, RDMBiasFacCurr)

        real(dp), intent(in) :: p_gen
        real(dp), intent(in) :: SignCurr
        real(dp), intent(out) :: RDMBiasFacCurr
        real(dp), intent(in) :: p_spawn_rdmfac
        real(dp) :: p_notlist_rdmfac, p_spawn, p_not_spawn, p_max_walktospawn
        character(len=*), parameter :: t_r = 'attempt_create_normal'

        ! We eventually turn this real bias factor into an integer to be passed
        ! around with the spawned children and their parents - this only works
        ! with 64 bit at the moment.
        if (n_int .eq. 4) call stop_all(t_r, 'The bias factor currently does not work with 32 bit integers.')

        ! Otherwise calculate the 'sign' of Di we are eventually going to add
        ! in as Di.Dj. Because we only add in Di.Dj when we successfully spawn
        ! from Di.Dj, we need to unbias (scale up) Di by the probability of this
        ! happening. We need the probability that the determinant i, with
        ! population n_i, will spawn on j. We only consider one instance of a
        ! pair Di,Dj, so just want the probability of any of the n_i walkers
        ! spawning at least once on Dj.

        ! P_successful_spawn(j | i)[n_i] =  1 - P_not_spawn(j | i)[n_i]
        ! P_not_spawn(j | i )[n_i] is the probability of none of the n_i walkers spawning on j from i.
        ! This requires either not generating j, or generating j and not succesfully spawning, n_i times.
        ! P_not_spawn(j | i )[n_i] = [(1 - P_gen(j | i)) + ( P_gen( j | i ) * (1 - P_spawn(j | i))]^n_i

        p_notlist_rdmfac = ( 1.0_dp - p_gen ) + ( p_gen * (1.0_dp - p_spawn_rdmfac) )

        ! The bias fac is now n_i / P_successful_spawn(j | i)[n_i].

        if (real(int(SignCurr),dp) .ne. SignCurr) then
            ! There's a non-integer population on this determinant. We need to
            ! consider both possibilities - whether we attempted to spawn 
            ! int(SignCurr) times or int(SignCurr)+1 times.
            p_max_walktospawn = abs(SignCurr-real(int(SignCurr),dp))
            p_not_spawn = (1.0_dp - p_max_walktospawn)*(p_notlist_rdmfac**abs(int(SignCurr))) + &
                        p_max_walktospawn*(p_notlist_rdmfac**(abs(int(SignCurr))+1))

        else
            p_not_spawn = p_notlist_rdmfac**(abs(SignCurr))
        end if

        p_spawn = abs(1.0_dp - p_not_spawn)
        
        ! Always use instantaneous signs for stochastically sampled off-diag
        ! elements (see CMO thesis).
        RDMBiasFacCurr = SignCurr / p_spawn

    end subroutine calc_rdmbiasfac

    subroutine store_parent_with_spawned(RDMBiasFacCurr, WalkerNumber, iLutI, DetSpawningAttempts, iLutJ, procJ, part_type)

        ! We are spawning from iLutI to SpawnedParts(:,ValidSpawnedList(proc)).
        ! This routine stores the parent (D_i) with the spawned child (D_j) so
        ! that we can add in Ci.Cj to the RDM later on. The parent is NIfDBO
        ! integers long, and stored in the second part of the SpawnedParts array 
        ! from NIfTot+1 -> NIfTot+1 + NIfDBO.

        use DetBitOps, only: DetBitEQ
        use FciMCData, only: SpawnedParts, ValidSpawnedList, TempSpawnedParts, TempSpawnedPartsInd

        real(dp), intent(in) :: RDMBiasFacCurr
        integer, intent(in) :: WalkerNumber, procJ
        integer, intent(in) :: DetSpawningAttempts
        integer(kind=n_int), intent(in) :: iLutI(0:niftot),iLutJ(0:niftot)
        integer, intent(in) :: part_type
        logical :: tRDMStoreParent
        integer :: j

        if (RDMBiasFacCurr .eq. 0.0_dp) then
            ! If RDMBiasFacCurr is exactly zero, any contribution from Ci.Cj will be zero 
            ! so it is not worth carrying on. 
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

            do j = 1, TempSpawnedPartsInd
                if (DetBitEQ(iLutJ(0:NIfDBO), TempSpawnedParts(0:NIfDBO,j), NIfDBO)) then
                    ! If this Dj is found, we do not want to store the parent with this spawned walker.
                    tRDMStoreParent = .false.
                    exit
                end if
            end do

            if (tRDMStoreParent) then
                ! This is a new Dj that has been spawned from this Di.
                ! We want to store it in the temporary list of spawned parts which have come from this Di.
                if (WalkerNumber .ne. DetSpawningAttempts) then
                    ! Don't bother storing these if we're on the last walker, or if we only have one 
                    ! walker on Di.
                    TempSpawnedPartsInd = TempSpawnedPartsInd + 1
                    TempSpawnedParts(0:NIfDBO,TempSpawnedPartsInd) = iLutJ(0:NIfDBO)
                end if

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
            end if
        end if

    end subroutine store_parent_with_spawned

end module rdm_general
