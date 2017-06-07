#include "macros.h"

module rdm_general

    use bit_rep_data, only: NIfTot, NIfDBO
    use constants
    use SystemData, only: nel, nbasis

    implicit none

contains

    subroutine init_rdms(nrdms_standard, nrdms_transition)

        use DeterminantData, only: write_det
        use CalcData, only: MemoryFacPart
        use FciMCData, only: MaxSpawned, Spawned_Parents, Spawned_Parents_Index
        use FciMCData, only: Spawned_ParentsTag, Spawned_Parents_IndexTag
        use FciMCData, only: HFDet_True, tSinglePartPhase, AvNoatHF, IterRDM_HF
        use global_det_data, only: len_av_sgn_tot, len_iter_occ_tot
        use LoggingData, only: tDo_Not_Calc_2RDM_est, RDMExcitLevel, tExplicitAllRDM
        use LoggingData, only: tDiagRDM, tDumpForcesInfo, tDipoles, tPrint1RDM
        use LoggingData, only: tRDMInstEnergy, tReadRDMs, tPopsfile, tno_RDMs_to_read
        use LoggingData, only: twrite_RDMs_to_read, tPrint1RDMsFrom2RDMPops
        use LoggingData, only: tPrint1RDMsFromSpinfree
        use Parallel_neci, only: iProcIndex, nProcessors
        use rdm_data, only: rdm_estimates, one_rdms, two_rdm_spawn, two_rdm_main, two_rdm_recv
        use rdm_data, only: two_rdm_recv_2, tOpenShell, print_2rdm_est, Sing_ExcDjs, Doub_ExcDjs
        use rdm_data, only: Sing_ExcDjs2, Doub_ExcDjs2, Sing_ExcDjsTag, Doub_ExcDjsTag
        use rdm_data, only: Sing_ExcDjs2Tag, Doub_ExcDjs2Tag, OneEl_Gap, TwoEl_Gap
        use rdm_data, only: Sing_InitExcSlots, Doub_InitExcSlots, Sing_ExcList, Doub_ExcList
        use rdm_data, only: nElRDM_Time, FinaliseRDMs_time, RDMEnergy_time, states_for_transition_rdm
        use rdm_data, only: rdm_main_size_fac, rdm_spawn_size_fac, rdm_recv_size_fac
        use rdm_data, only: rdm_definitions
        use rdm_data_utils, only: init_rdm_spawn_t, init_rdm_list_t, init_one_rdm_t
        use rdm_data_utils, only: init_rdm_definitions_t, clear_one_rdms, clear_rdm_list_t
        use rdm_estimators, only: init_rdm_estimates_t, calc_2rdm_estimates_wrapper
        use rdm_reading
        use RotateOrbsData, only: SymLabelCounts2_rot,SymLabelList2_rot, SymLabelListInv_rot
        use RotateOrbsData, only: SymLabelCounts2_rotTag, SymLabelList2_rotTag, NoOrbs
        use RotateOrbsData, only: SymLabelListInv_rotTag, SpatOrbs, NoSymLabelCounts
        use SystemData, only: tStoreSpinOrbs, tHPHF, tFixLz, iMaxLz, tROHF, LMS
        use MemoryManager, only: LogMemAlloc

        integer, intent(in) :: nrdms_standard, nrdms_transition

        integer :: nrdms, rdm_nrows, nhashes_rdm_main, nhashes_rdm_spawn
        integer :: standard_spawn_size, min_spawn_size
        integer :: max_nelems_main, max_nelems_spawn, max_nelems_recv, max_nelems_recv_2
        integer :: memory_alloc, main_mem, spawn_mem, recv_mem
        integer :: irdm, iproc, ierr
        character(len=*), parameter :: t_r = 'init_rdms'

#ifdef __CMPLX
        call stop_all(t_r, 'Filling of reduced density matrices not working with complex walkers yet.')
#endif

        nrdms = nrdms_standard + nrdms_transition

        ! Only spatial orbitals for the 2-RDMs (and F12).
        if (tStoreSpinOrbs .and. RDMExcitLevel /= 1) then
            call stop_all(t_r, '2-RDM calculations not set up for systems stored as spin orbitals.')
        end if

        if (tROHF .or. tStoreSpinOrbs.or.LMS.ne.0) then
            tOpenShell = .true.
        else
            tOpenShell = .false.
        end if

        if (tExplicitAllRDM) then
            write(6,'(1X,"Explicitly calculating the reduced density matrices from the FCIQMC wavefunction.")')
        else
            write(6,'(1X,"Stochastically calculating the reduced density matrices from the FCIQMC wavefunction")')
            write(6,'(1X,"incl. explicit connections to the following HF determinant:")', advance='no')
            call write_det(6, HFDet_True, .true.)
        end if

        if (RDMExcitLevel == 1) then
            print_2rdm_est = .false.
        else
            ! If the RDMExcitLevel is 2 or 3 - and we're calculating the 2-RDM, 
            ! then we automatically calculate the energy (and other estimates!)
            ! unless we're specifically told not to.
            if (tDo_Not_Calc_2RDM_est) then
                print_2rdm_est = .false.
            else
                print_2rdm_est = .true.
                write(6,'(1X,"Calculating the energy from the reduced density matrix. &
                              &This requires the 2 electron RDM from which the 1-RDM can also be constructed.")')
            end if
        end if

        call init_rdm_definitions_t(rdm_definitions, nrdms_standard, nrdms_transition, states_for_transition_rdm)

        ! Allocate arrays for holding averaged signs and block lengths for the
        ! HF determinant.
        allocate(AvNoatHF(len_av_sgn_tot))
        AvNoatHF = 0.0_dp
        allocate(IterRDM_HF(len_iter_occ_tot))
        IterRDM_HF = 0

        ! Have not got HPHF working with the explicit or truncated methods yet.
        ! Neither of these would be too difficult to implement.
        if (tHPHF .and. tExplicitAllRDM) call stop_all(t_r, 'HPHF not set up with the explicit calculation of the RDM.')

        if (tDipoles) then
            write(6,'("WARNING - The calculation of dipole moments is currently not supported for the new RDM code. &
                      &Use the OLDRDMS option to use feature.")')
        end if

        SpatOrbs = nbasis/2
        if (tOpenShell) then
            NoOrbs = nbasis
        else
            NoOrbs = SpatOrbs
        end if

        ! The memory of (large) alloctaed arrays, per MPI process.
        memory_alloc = 0

        ! For now, create RDM array big enough so that *all* RDM elements on
        ! a particular processor can be stored, using the usual approximations
        ! to take symmetry into account. Include a factor of 1.5 to account for
        ! factors such as imperfect load balancing (which affects the spawned
        ! array).
        rdm_nrows = nbasis*(nbasis-1)/2
        max_nelems_main = 1.5*(rdm_nrows**2)/(8*nProcessors)*rdm_main_size_fac
        nhashes_rdm_main = 0.75*max_nelems_main*rdm_main_size_fac

        main_mem = max_nelems_main*(nrdms+1)*size_int_rdm
        write(6,'(/,1X,"About to allocate main RDM array, size per MPI process (MB):", f14.6)') real(main_mem,dp)/1048576.0_dp
        call init_rdm_list_t(two_rdm_main, nrdms, max_nelems_main, nhashes_rdm_main)
        write(6,'(1X,"Allocation of main RDM array complete.")')

        ! Factor of 10 over perfectly distributed size, for some safety.
        standard_spawn_size = 10.0*(rdm_nrows**2)/(8*nProcessors)
        ! For cases where we have a small number of orbitals but large number
        ! of processors (i.e., large CASSCF calculations), we may find the
        ! above standard_spawn_size is less than nProcessors. Thus, there
        ! would not be at least one spawning slot per processor. In such cases
        ! make sure that we have at least 50 per processor, for some safety.
        min_spawn_size = 50*nProcessors
        max_nelems_spawn = max(standard_spawn_size, min_spawn_size)*rdm_spawn_size_fac
        nhashes_rdm_spawn = 0.75*max_nelems_spawn*rdm_spawn_size_fac

        spawn_mem = max_nelems_spawn*(nrdms+1)*size_int_rdm
        write(6,'(1X,"About to allocate RDM spawning array, size per MPI process (MB):", f14.6)') real(spawn_mem,dp)/1048576.0_dp
        call init_rdm_spawn_t(two_rdm_spawn, rdm_nrows, nrdms, max_nelems_spawn, nhashes_rdm_spawn)
        write(6,'(1X,"Allocation of RDM spawning array complete.")')

        max_nelems_recv = 4.0*(rdm_nrows**2)/(8*nProcessors)*rdm_recv_size_fac
        max_nelems_recv_2 = 2.0*(rdm_nrows**2)/(8*nProcessors)*rdm_recv_size_fac

        recv_mem = (max_nelems_recv + max_nelems_recv_2)*(nrdms+1)*size_int_rdm
        write(6,'(1X,"About to allocate RDM receiving arrays, size per MPI process (MB):", f14.6)') real(recv_mem,dp)/1048576.0_dp
        ! Don't need the hash table for the received list, so pass 0 for nhashes.
        call init_rdm_list_t(two_rdm_recv, nrdms, max_nelems_recv, 0)
        call init_rdm_list_t(two_rdm_recv_2, nrdms, max_nelems_recv_2, 0)
        write(6,'(1X,"Allocation of RDM receiving arrays complete.",/)')

        ! Count the memory the various RDM lists (but this does *not* count
        ! the memory of the hash tables - this will increase dynamically
        ! throughout the simulation).
        memory_alloc = memory_alloc + main_mem + spawn_mem + recv_mem

        call init_rdm_estimates_t(rdm_estimates, nrdms_standard, nrdms_transition, print_2rdm_est)

        ! Initialise 1-RDM objects.
        if (RDMExcitLevel == 1 .or. RDMExcitLevel == 3 .or. &
              tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) then

            allocate(one_rdms(nrdms), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem allocating one_rdms array.')

            do irdm = 1, nrdms
                call init_one_rdm_t(one_rdms(irdm), NoOrbs)

                memory_alloc = memory_alloc + ( NoOrbs * NoOrbs * 8 )
            end do
        end if

        ! We then need to allocate the arrays for excitations etc when doing
        ! the explicit all calculation.
        if (tExplicitAllRDM) then
            ! We always calculate the single stuff - and if RDMExcitLevel is 1,
            ! this is all, otherwise calculate the double stuff too.

            ! This array actually contains the excitations in blocks of the
            ! processor they will be sent to. Only needed if the 1-RDM is the
            ! only thing being calculated.
            allocate(Sing_ExcDjs(0:NIfTot, nint((nel*nbasis)*MemoryFacPart)), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem allocating Sing_ExcDjs array.')
            call LogMemAlloc('Sing_ExcDjs', nint(nel*nbasis*MemoryFacPart)*(NIfTot+1),&
                                            size_n_int, t_r, Sing_ExcDjsTag, ierr)

            allocate(Sing_ExcDjs2(0:NIfTot, nint((nel*nbasis)*MemoryFacPart)), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem allocating Sing_ExcDjs2 array.')
            call LogMemAlloc('Sing_ExcDjs2', nint(nel*nbasis*MemoryFacPart)*(NIfTot+1),&
                                            size_n_int, t_r, Sing_ExcDjs2Tag, ierr)

            Sing_ExcDjs(:,:) = 0
            Sing_ExcDjs2(:,:) = 0

            memory_alloc = memory_alloc + ( (NIfTot + 1) * nint((nel*nbasis)*MemoryFacPart) * size_n_int * 2 )

            ! We need room to potentially generate N*M single excitations but
            ! these will be spread across each processor.
            OneEl_Gap = (real(nel,dp)*real(nbasis,dp)*MemoryFacPart)/real(nProcessors,dp)

            ! This array contains the initial positions of the excitations
            ! for each processor.
            allocate(Sing_InitExcSlots(0:(nProcessors-1)), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem allocating Sing_InitExcSlots array,')
            do iproc = 0, nProcessors - 1
                Sing_InitExcSlots(iproc) = nint(OneEl_Gap*iproc) + 1
            end do

            ! This array contains the current position of the excitations as
            ! they're added.
            allocate(Sing_ExcList(0:(nProcessors-1)), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem allocating Sing_ExcList array,')
            Sing_ExcList(:) = Sing_InitExcSlots(:)

            if (RDMExcitLevel /= 1) then
                ! This array actually contains the excitations in blocks of
                ! the processor they will be sent to.
                allocate(Doub_ExcDjs(0:NIfTot,nint(((nel*nbasis)**2)*MemoryFacPart)), stat=ierr)
                if (ierr /= 0) call stop_all(t_r,'Problem allocating Doub_ExcDjs array.')
                call LogMemAlloc('Doub_ExcDjs', nint(((nel*nbasis)**2)*MemoryFacPart)&
                                *(NIfTot+1), size_n_int, t_r, Doub_ExcDjsTag, ierr)

                allocate(Doub_ExcDjs2(0:NIfTot, nint(((nel*nbasis)**2)*MemoryFacPart)), stat=ierr)
                if (ierr /= 0) call stop_all(t_r, 'Problem allocating Doub_ExcDjs2 array.')
                call LogMemAlloc('Doub_ExcDjs2',nint(((nel*nbasis)**2)*MemoryFacPart)&
                                *(NIfTot+1), size_n_int, t_r, Doub_ExcDjs2Tag, ierr)

                memory_alloc = memory_alloc + ( (NIfTot + 1) * nint(((nel*nbasis)**2)*MemoryFacPart) * size_n_int * 2 )

                ! We need room to potentially generate (N*M)^2 double excitations
                ! but these will be spread across each processor.        
                TwoEl_Gap = (((real(nel,dp)*real(nbasis,dp))**2)*MemoryFacPart)/real(nProcessors,dp)

                Doub_ExcDjs(:,:) = 0
                Doub_ExcDjs2(:,:) = 0

                ! This array contains the initial positions of the excitations
                ! for each processor.
                allocate(Doub_InitExcSlots(0:(nProcessors-1)), stat=ierr)
                if (ierr /= 0) call stop_all(t_r, 'Problem allocating Doub_InitExcSlots array,')
                do iproc = 0, nProcessors-1
                    Doub_InitExcSlots(iproc) = nint(TwoEl_Gap*iproc) + 1
                end do

                ! This array contains the current position of the excitations
                ! as they're added.
                allocate(Doub_ExcList(0:(nProcessors-1)), stat=ierr)
                if (ierr /= 0) call stop_all(t_r, 'Problem allocating Doub_ExcList array,')
                Doub_ExcList(:) = Doub_InitExcSlots(:)
            end if

        else

            ! Finally, we need to hold onto the parents of the spawned particles.
            ! This is not necessary if we're doing completely explicit calculations.
            allocate(Spawned_Parents(0:(NIfDBO+2), MaxSpawned), stat=ierr)
            if (ierr /= 0) call stop_all(t_r,'Problem allocating Spawned_Parents array,')
            call LogMemAlloc('Spawned_Parents', MaxSpawned*(NIfDBO+3), size_n_int,&
                                                t_r,Spawned_ParentsTag, ierr)
            allocate(Spawned_Parents_Index(2,MaxSpawned),stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem allocating Spawned_Parents_Index array,')
            call LogMemAlloc('Spawned_Parents_Index', MaxSpawned*2,4, t_r,&
                                                        Spawned_Parents_IndexTag, ierr)

            memory_alloc = memory_alloc + ( (NIfTot + 2) * MaxSpawned * size_n_int ) 

            memory_alloc = memory_alloc + ( 2 * MaxSpawned * 4 ) 

        end if

        if (iProcIndex == 0) then
            write(6, "(A,F14.6,A,F14.6,A)") " Main RDM memory arrays consists of: ", &
                      real(memory_alloc,dp)/1048576.0_dp," MB per MPI process."
        end if

        ! These parameters are set for the set up of the symmetry arrays, which
        ! are later used for the diagonalisation / rotation of the 1-RDMs.

        if (tOpenShell) then
            if (tFixLz) then
                NoSymLabelCounts = 16*((2 * iMaxLz) + 1)
            else
                NoSymLabelCounts = 16 
            end if
        else
            if (tFixLz) then
                NoSymLabelCounts = 8*((2 * iMaxLz) + 1)
            else
                NoSymLabelCounts = 8
            end if
        end if

        if (RDMExcitLevel == 1 .or. RDMExcitLevel == 3 .or. tDiagRDM .or. tPrint1RDM .or. tDumpForcesInfo .or. tDipoles) then
            ! These arrays contain indexing systems to order the 1-RDM orbitals
            ! in terms of symmetry. This allows the diagonalisation of the RDMs
            ! to be done in symmetry blocks (a lot quicker/easier).

            if (.not. allocated(SymLabelCounts2_rot)) allocate(SymLabelCounts2_rot(2, NoSymLabelCounts), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem allocating SymLabelCounts2_rot array,')
            call LogMemAlloc('SymLabelCounts2_rot', 2*NoSymLabelCounts, 4, t_r, SymLabelCounts2_rotTag, ierr)
            SymLabelCounts2_rot = 0

            allocate(SymLabelList2_rot(NoOrbs), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem allocating SymLabelList2_rot array,')
            call LogMemAlloc('SymLabelList2_rot', NoOrbs, 4, t_r, SymLabelList2_rotTag, ierr)
            SymLabelList2_rot = 0
     
            allocate(SymLabelListInv_rot(NoOrbs), stat=ierr)
            if (ierr /= 0) call stop_all(t_r, 'Problem allocating SymLabelListInv_rot array,')
            call LogMemAlloc('SymLabelListInv_rot', NoOrbs, 4, t_r, SymLabelListInv_rotTag, ierr)
            SymLabelListInv_rot = 0

            ! This routine actually sets up the symmetry labels for the 1-RDM.
            call SetUpSymLabels_RDM()
        end if

        if (tPrint1RDMsFromSpinfree) then
            call read_spinfree_2rdm_files(rdm_definitions, two_rdm_main, two_rdm_spawn)
            call print_1rdms_from_sf2rdms_wrapper(rdm_definitions, one_rdms, two_rdm_main)
            ! now clear these objects before the main simulation.
            call clear_one_rdms(one_rdms)
            call clear_rdm_list_t(two_rdm_main)
        end if

        if (tReadRDMs) then
            if (RDMExcitLevel == 1) then
                do irdm = 1, size(one_rdms)
                    call read_1rdm(rdm_definitions, one_rdms(irdm), irdm)
                end do
            else
                call read_2rdm_popsfile(two_rdm_main, two_rdm_spawn)
                if (print_2rdm_est) call calc_2rdm_estimates_wrapper(rdm_definitions, rdm_estimates, two_rdm_main)

                if (tPrint1RDMsFrom2RDMPops) then
                    call print_1rdms_from_2rdms_wrapper(rdm_definitions, one_rdms, two_rdm_main, tOpenShell)
                end if
            end if

            if (any(tSinglePartPhase)) then
                write(6,'("WARNING - Asking to read in the RDMs, but not varying shift from the beginning of &
                          &the calculation. All RDMs just read in will be zeroed, to prevent invalid averaging.")')
                ! Clear these objects, before the main simulation, since we
                ! haven't started averaging RDMs yet.
                call clear_one_rdms(one_rdms)
                call clear_rdm_list_t(two_rdm_main)
                ! Turn off tReadRDMs, since the read in RDMs aren't being
                ! used. Leaving it on affects some other stuff later.
                tReadRDMs = .false.
            end if
        end if

        ! By default, if we're writing out a popsfile (and doing an RDM
        ! calculation), we also write out the unnormalised RDMs that can be
        ! read in when restarting a calculation. If the NORDMSTOREAD option
        ! is on, these wont be printed.
        if (tPopsfile .and. (.not. tno_RDMs_to_read)) twrite_RDMs_to_read = .true.

        if (iProcIndex == 0) write(6,'(1X,"RDM memory allocation complete.",/)')

        nElRDM_Time%timer_name = 'nElRDMTime'
        FinaliseRDMs_Time%timer_name = 'FinaliseRDMsTime'
        RDMEnergy_Time%timer_name = 'RDMEnergyTime'

    end subroutine init_rdms

    subroutine SetUpSymLabels_RDM()

        ! This routine just sets up the symmetry labels so that the orbitals
        ! are ordered according to symmetry (all beta then all alpha if spin orbs).

        use rdm_data, only: tOpenShell
        use RotateOrbsData, only: SymLabelList2_rot, SymLabelCounts2_rot, SymLabelListInv_rot
        use RotateOrbsData, only: NoOrbs, SpatOrbs, NoSymLabelCounts
        use sort_mod, only: sort
        use SystemData, only: G1, BRR, lNoSymmetry, tFixLz, iMaxLz
        use UMatCache, only: gtID
        use MemoryManager, only: LogMemAlloc, LogMemDealloc

        integer, allocatable :: SymOrbs_rot(:)
        integer :: SymOrbs_rotTag, ierr, i, j, SpatSym, LzSym 
        integer :: lo, hi, Symi, SymCurr, Symi2, SymCurr2
        character(len=*), parameter :: t_r = 'SetUpSymLabels_RDM'

        ! This is only allocated temporarily to be used to order the orbitals by.
        allocate(SymOrbs_rot(NoOrbs), stat=ierr)
        call LogMemAlloc('SymOrbs_rot', NoOrbs, 4, t_r, SymOrbs_rotTag, ierr)
        if (ierr /= 0) call stop_all(t_r, "Mem allocation for SymOrbs_rot failed.")

        ! Now we want to put the spatial orbital index, followed by the symmetry.
        SymLabelList2_rot(:) = 0
        SymOrbs_rot(:) = 0

        ! *** STEP 1 *** Fill SymLabelList2_rot.
        ! Find the orbitals and order them in terms of symmetry.
        do i = 1, SpatOrbs
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
                if (Symi /= SymCurr) then
                    do j = SymCurr + 1, Symi
                        SymLabelCounts2_rot(1,(j+1)) = i
                    end do
                    SymCurr = Symi
                end if
                if (tOpenShell) then
                    SymLabelCounts2_rot(2,(Symi2+9)) = SymLabelCounts2_rot(2,(Symi2+9))+1
                    if (Symi2 /= SymCurr2) then
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
            if (SymLabelCounts2_rot(2,i) /= 0) then
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

    subroutine realloc_SpawnedParts()

        ! Routine called when RDM accumulation is turned on, usually midway
        ! through an FCIQMC simulation.

        ! When calculating the RDMs, we need to store the parent from which a
        ! child is spawned along with the children in the spawned array. This
        ! means a slightly larger array is communicated between processors,
        ! which there is no point in doing for the first part of the calculation.
        ! When we start calculating the RDMs this routine is called and the
        ! SpawnedParts array is made larger to accommodate the parents.

        use bit_rep_data, only: nifbcast, NOffParent, bit_rdm_init
        use FciMCData, only: MaxSpawned, SpawnVec, SpawnVec2, SpawnVecTag, SpawnVec2Tag
        use FciMCData, only: SpawnedParts, SpawnedParts2
        use MemoryManager, only: LogMemAlloc, LogMemDealloc
        use util_mod_byte_size

        integer :: ierr, nifbcast_old
        character(len=*), parameter :: t_r = 'realloc_SpawnedParts'

        if (bit_rdm_init) &
            call stop_all(t_r, 'RDM broadcast representation already initialised')
        
        deallocate(SpawnVec)
        call LogMemDealloc(t_r,SpawnVecTag)
        deallocate(SpawnVec2)
        call LogMemDealloc(t_r,SpawnVec2Tag)

        ! Resize the RDM arrays
        NIfBCast_old = NIfBCast
        NOffParent = NIfBCast + 1

        NIfBCast = NIfBCast + NIfDBO + 2

        allocate(SpawnVec(0:NIfBCast, MaxSpawned), SpawnVec2(0:NIfBCast, MaxSpawned), stat=ierr)
        log_alloc(SpawnVec, SpawnVecTag, ierr)
        log_alloc(SpawnVec2, SpawnVec2Tag, ierr)

        ! Point at correct spawning arrays
        SpawnedParts => SpawnVec
        SpawnedParts2 => SpawnVec2

        write(6,'(A54,F10.4,A4,F10.4,A13)') &
            'Memory requirement for spawned arrays increased from ',&
            real(((NIfBCast_old+1)*MaxSpawned*2*size_n_int),dp)/1048576.0_dp, &
            ' to ', &
            real(((NIfBCast+1)*MaxSpawned*2*size_n_int),dp)/1048576.0_dp, &
            ' Mb/Processor'

        ! And we are done
        bit_rdm_init = .true.

    end subroutine realloc_SpawnedParts

    subroutine dealloc_global_rdm_data()

        ! This routine just deallocates the arrays allocated in InitRDMs.
        ! If the NECI calculation softexits before the RDMs start to fill,
        ! this is all that is called at the end.

        use FciMCData, only: Spawned_Parents, Spawned_Parents_Index
        use FciMCData, only: Spawned_ParentsTag, Spawned_Parents_IndexTag
        use FciMCData, only: AvNoatHF, IterRDM_HF
        use LoggingData, only: RDMExcitLevel, tExplicitAllRDM
        use rdm_data, only: two_rdm_main, two_rdm_recv, two_rdm_recv_2, two_rdm_spawn
        use rdm_data, only: rdm_estimates, one_rdms, Sing_ExcDjs, Doub_ExcDjs
        use rdm_data, only: Sing_ExcDjs2, Doub_ExcDjs2, Sing_ExcDjsTag, Doub_ExcDjsTag
        use rdm_data, only: Sing_ExcDjs2Tag, Doub_ExcDjs2Tag
        use rdm_data, only: Sing_InitExcSlots, Doub_InitExcSlots, Sing_ExcList, Doub_ExcList
        use rdm_data_old, only: rdms
        use rdm_data_utils, only: dealloc_rdm_list_t, dealloc_rdm_spawn_t, dealloc_one_rdm_t
        use rdm_estimators, only: dealloc_rdm_estimates_t
        use RotateOrbsData, only: SymLabelCounts2_rot, SymLabelList2_rot, SymLabelListInv_rot
        use RotateOrbsData, only: SymLabelCounts2_rotTag, SymLabelList2_rotTag
        use RotateOrbsData, only: SymLabelListInv_rotTag
        use RotateOrbsMod, only: FourIndInts, FourIndIntsTag
        use util_mod, only: LogMemDealloc

        character(len=*), parameter :: t_r = 'dealloc_global_rdm_data'

        integer :: irdm

        ! Deallocate global 2-RDM arrays, including the spawning object.
        call dealloc_rdm_list_t(two_rdm_main)
        call dealloc_rdm_list_t(two_rdm_recv)
        call dealloc_rdm_list_t(two_rdm_recv_2)
        call dealloc_rdm_spawn_t(two_rdm_spawn)

        ! Deallocate the RDM estimates object.
        call dealloc_rdm_estimates_t(rdm_estimates)

        ! Deallocate all 1-RDM objects.
        if (allocated(one_rdms)) then
            do irdm = 1, size(one_rdms)
                call dealloc_one_rdm_t(one_rdms(irdm))
            end do
            deallocate(one_rdms)
        end if

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

            if (RDMExcitLevel /= 1) then
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

        if (allocated(AvNoatHF)) deallocate(AvNoatHF)
        if (allocated(IterRDM_HF)) deallocate(IterRDM_HF)

    end subroutine dealloc_global_rdm_data

    ! Some general routines used during the main simulation.

    subroutine extract_bit_rep_avsign_no_rdm(rdm_defs, iLutnI, j, nI, SignI, FlagsI, IterRDMStartI, AvSignI, store)

        ! This is just the standard extract_bit_rep routine for when we're not
        ! calculating the RDMs.    

        use bit_reps, only: extract_bit_rep
        use FciMCData, only: excit_gen_store_type
        use global_det_data, only: len_av_sgn_tot, len_iter_occ_tot
        use rdm_data, only: rdm_definitions_t

        type(rdm_definitions_t), intent(in) :: rdm_defs
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(in) :: j
        integer, intent(out) :: nI(nel), FlagsI
        real(dp), intent(out) :: SignI(lenof_sign)
        real(dp), intent(out) :: IterRDMStartI(len_iter_occ_tot), AvSignI(len_av_sgn_tot)
        type(excit_gen_store_type), intent(inout), optional :: store

        integer :: iunused
        
        ! This extracts everything.
        call extract_bit_rep (iLutnI, nI, SignI, FlagsI, store)

        IterRDMStartI = 0.0_dp
        AvSignI = 0.0_dp

        ! Eliminate warnings
        iunused = j

    end subroutine extract_bit_rep_avsign_no_rdm

    subroutine extract_bit_rep_avsign_norm(rdm_defs, iLutnI, j, nI, SignI, FlagsI, IterRDMStartI, AvSignI, store)

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
        use CalcData, only: tPairedReplicas
        use FciMCData, only: PreviousCycles, Iter, IterRDMStart, excit_gen_store_type
        use global_det_data, only: get_iter_occ_tot, get_av_sgn_tot
        use global_det_data, only: len_av_sgn_tot, len_iter_occ_tot
        use rdm_data, only: rdm_definitions_t
        use LoggingData, only: RDMEnergyIter

        type(rdm_definitions_t), intent(in) :: rdm_defs
        integer(n_int), intent(in) :: iLutnI(0:nIfTot)
        integer, intent(out) :: nI(nel), FlagsI
        integer, intent(in) :: j
        real(dp), intent(out) :: SignI(lenof_sign)
        real(dp), intent(out) :: IterRDMStartI(len_iter_occ_tot), AvSignI(len_av_sgn_tot)
        type(excit_gen_store_type), intent(inout), optional :: store

        integer :: part_ind, irdm, iunused
        integer :: av_ind_1, av_ind_2

        ! This is the iteration from which this determinant has been occupied.
        IterRDMStartI = get_iter_occ_tot(j)

        ! This extracts everything.
        call extract_bit_rep (iLutnI, nI, SignI, FlagsI)

        associate(ind => rdm_defs%sim_labels)

        if (((Iter+PreviousCycles-IterRDMStart) > 0) .and. &
            & (mod(((Iter-1)+PreviousCycles - IterRDMStart + 1), RDMEnergyIter) == 0)) then 

            ! The previous iteration was one where we added in diagonal elements
            ! To keep things unbiased, we need to set up a new averaging block now.
            ! NB: if doing single run cutoff, note that doing things this way is now
            ! NOT the same as the technique described in CMO (and DMC's) thesis.
            ! Would expect diagonal elements to be slightly worse quality, improving
            ! as one calculates the RDM energy less frequently.  As this method is
            ! biased anyway, I'm not going to lose sleep over it.
            do irdm = 1, rdm_defs%nrdms
                av_ind_1 = irdm*2-1
                av_ind_2 = irdm*2

                AvSignI(av_ind_1) = SignI(ind(1,irdm))
                IterRDMStartI(av_ind_1) = real(Iter + PreviousCycles,dp)
                AvSignI(av_ind_2) = SignI(ind(2,irdm))
                IterRDMStartI(av_ind_2) = real(Iter + PreviousCycles,dp)
            end do
        else
            ! Now let's consider other instances in which we need to start a new block:
            do irdm = 1, rdm_defs%nrdms
                ! The indicies of the first and second replicas in this
                ! particular pair, in the *average* sign arrays.
                av_ind_1 = irdm*2-1
                av_ind_2 = irdm*2

                if ((SignI(ind(1,irdm)) == 0) .and. (IterRDMStartI(av_ind_1) /= 0)) then
                    ! The population has just gone to zero on population 1.
                    ! Therefore, we need to start a new averaging block.
                    AvSignI(av_ind_1) = 0
                    IterRDMStartI(av_ind_1) = 0
                    AvSignI(av_ind_2) = SignI(ind(2,irdm))
                    IterRDMStartI(av_ind_2) = real(Iter + PreviousCycles,dp)

                else if ((SignI(ind(2,irdm)) == 0) .and. (IterRDMStartI(av_ind_2) /= 0)) then
                    ! The population has just gone to zero on population 2.
                    ! Therefore, we need to start a new averaging block.
                    AvSignI(av_ind_2) = 0
                    IterRDMStartI(av_ind_2) = 0
                    AvSignI(av_ind_1) = SignI(ind(1,irdm))
                    IterRDMStartI(av_ind_1) = real(Iter + PreviousCycles,dp)

                else if ((SignI(ind(1,irdm)) /= 0) .and. (IterRDMStartI(av_ind_1) == 0)) then
                    ! Population 1 has just become occupied.
                    IterRDMStartI(av_ind_1) = real(Iter + PreviousCycles,dp)
                    IterRDMStartI(av_ind_2) = real(Iter + PreviousCycles,dp)
                    AvSignI(av_ind_1) = SignI(ind(1,irdm))
                    AvSignI(av_ind_2) = SignI(ind(2,irdm))
                    if (SignI(ind(2,irdm)) == 0) IterRDMStartI(av_ind_2) = 0

                else if ((SignI(ind(2,irdm)) /= 0) .and. (IterRDMStartI(av_ind_2) == 0)) then
                    ! Population 2 has just become occupied.
                    IterRDMStartI(av_ind_1) = real(Iter + PreviousCycles,dp)
                    IterRDMStartI(av_ind_2) = real(Iter + PreviousCycles,dp)
                    AvSignI(av_ind_1) = SignI(ind(1,irdm))
                    AvSignI(av_ind_2) = SignI(ind(2,irdm))
                    if (SignI(ind(1,irdm)) == 0) IterRDMStartI(av_ind_1) = 0

                else
                    ! Nothing unusual has happened so update both average
                    ! populations as normal.
                    AvSignI(av_ind_1) = &
                        ( ((real(Iter+PreviousCycles,dp) - IterRDMStartI(av_ind_1)) * get_av_sgn_tot(j, av_ind_1)) &
                          + SignI(ind(1,irdm)) ) / ( real(Iter+PreviousCycles,dp) - IterRDMStartI(av_ind_1) + 1.0_dp )

                    AvSignI(av_ind_2) = &
                        ( ((real(Iter+PreviousCycles,dp) - IterRDMStartI(av_ind_2)) * get_av_sgn_tot(j, av_ind_2)) &
                          + SignI(ind(2,irdm)) ) / ( real(Iter+PreviousCycles,dp) - IterRDMStartI(av_ind_2) + 1.0_dp )
                end if
            end do
        end if

        end associate

        ! Eliminate warnings
        iunused = store%nopen

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
        if (n_int == 4) call stop_all(t_r, 'The bias factor currently does not work with 32 bit integers.')

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

        if (abs(real(int(SignCurr), dp) - SignCurr) > 1.0e-12_dp) then
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

    subroutine store_parent_with_spawned(RDMBiasFacCurr, WalkerNumber, iLutI, DetSpawningAttempts, iLutJ, procJ)

        ! We are spawning from iLutI to SpawnedParts(:,ValidSpawnedList(proc)).
        ! This routine stores the parent (D_i) with the spawned child (D_j) so
        ! that we can add in Ci.Cj to the RDM later on. The parent is NIfDBO
        ! integers long, and stored in the second part of the SpawnedParts array 
        ! from NIfTot+1 -> NIfTot+1 + NIfDBO.

        use DetBitOps, only: DetBitEQ
        use FciMCData, only: SpawnedParts, ValidSpawnedList, TempSpawnedParts, TempSpawnedPartsInd
        use bit_reps, only: zero_parent, encode_parent

        real(dp), intent(in) :: RDMBiasFacCurr
        integer, intent(in) :: WalkerNumber, procJ
        integer, intent(in) :: DetSpawningAttempts
        integer(n_int), intent(in) :: iLutI(0:niftot), iLutJ(0:niftot)
        logical :: tRDMStoreParent
        integer :: j

        if (abs(RDMBiasFacCurr) < 1.0e-12_dp) then
            ! If RDMBiasFacCurr is exactly zero, any contribution from Ci.Cj will be zero 
            ! so it is not worth carrying on. 
            call zero_parent(SpawnedParts(:, ValidSpawnedList(procJ)))
        else

            ! First we want to check if this Di.Dj pair has already been accounted for.
            ! This means searching the Dj's that have already been spawned from this Di, to make sure 
            ! the new Dj being spawned on here is not the same.
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
                if (WalkerNumber /= DetSpawningAttempts) then
                    ! Don't bother storing these if we're on the last walker, or if we only have one 
                    ! walker on Di.
                    TempSpawnedPartsInd = TempSpawnedPartsInd + 1
                    TempSpawnedParts(0:NIfDBO,TempSpawnedPartsInd) = iLutJ(0:NIfDBO)
                end if

                ! We also want to make sure the parent Di is stored with this Dj.

                ! We need to carry with the child (and the parent), the sign of the parent.
                ! In actual fact this is the sign of the parent divided by the probability of generating
                ! that pair Di and Dj, to account for the 
                ! fact that Di and Dj are not always added to the RDM, but only when Di spawns on Dj.
                ! This RDMBiasFacCurr factor is turned into an integer to pass around to the relevant processors.
                call encode_parent(SpawnedParts(:, ValidSpawnedList(procJ)), &
                                   ilutI, RDMBiasFacCurr)

            else
                ! This Di has already spawned on this Dj - don't store the Di parent with this child, 
                ! so that the pair is not double counted.  
                ! We are using the probability that Di spawns onto Dj *at least once*, so we don't want to 
                ! double count this pair.
                call zero_parent(SpawnedParts(:, ValidSpawnedList(procJ)))
            end if
        end if

    end subroutine store_parent_with_spawned

end module rdm_general
