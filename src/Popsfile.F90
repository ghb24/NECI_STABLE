! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
#include "macros.h"

MODULE PopsfileMod

    use SystemData, only: nel, tHPHF, tFixLz, tCSF, nBasis, tNoBrillouin, &
                          AB_elec_pairs, par_elec_pairs
    use CalcData, only: tTruncInitiator, DiagSft, tWalkContGrow, nEquilSteps, &
                        ScaleWalkers, tReadPopsRestart, &
                        InitWalkers, tReadPopsChangeRef, nShiftEquilSteps, &
                        iWeightPopRead, iPopsFileNoRead, Tau, &
                        InitiatorWalkNo, MemoryFacPart, &
                        MemoryFacSpawn, tSemiStochastic, tTrialWavefunction, &
                        tCCMC, pops_norm, tWritePopsNorm
    use DetBitOps, only: DetBitLT, FindBitExcitLevel, DetBitEQ, EncodeBitDet, &
                         ilut_lt, ilut_gt
    use hash, only: DetermineDetNode, FindWalkerHash, clear_hash_table, &
                    fill_in_hash_table
    use Determinants, only: get_helement, write_det
    use hphf_integrals, only: hphf_diag_helement
    USE dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData
    use bit_rep_data, only: extract_sign, NOffSgn
    use bit_reps
    use constants
    use Parallel_neci
    use LoggingData, only: iWritePopsEvery, tPopsFile, iPopsPartEvery, tBinPops, &
                       tPrintPopsDefault, tIncrementPops, tPrintInitiators, &
                       tSplitPops, tZeroProjE, tRDMonFly, tExplicitAllRDM, &
                       binarypops_min_weight
    use sort_mod
    use util_mod, only: get_free_unit,get_unique_filename
    use tau_search, only: gamma_sing, gamma_doub, gamma_opp, gamma_par, &
        max_death_cpt
    use global_det_data, only: global_determinant_data, set_iter_occ, &
                               init_global_det_data, set_det_diagH

    implicit none

    logical :: tRealPOPSfile

    interface
        subroutine ChangeRefDet(det)
            use SystemData, only: nel
            implicit none
            integer :: det(nel)
        end subroutine
    end interface

contains

    !   V.3/4 POPSFILE ROUTINES   !
    !This routine reads in particle configurations from a POPSFILE v.3-4.
    !EndPopsList is the number of entries in the POPSFILE to read, and ReadBatch is the number of determinants
    !which can be read in in a single batch.
    subroutine ReadFromPopsfile(EndPopsList, ReadBatch, CurrWalkers, &
                                CurrParts, CurrHF, Dets, DetsLen, &
                                pops_nnodes, pops_walkers, PopNifSgn, &
                                PopNel, tCalcExtraInfo)

        use MemoryManager, only: TagIntType

        integer(int64) , intent(in) :: EndPopsList  !Number of entries in the POPSFILE.
        integer , intent(in) :: ReadBatch       !Size of the batch of determinants to read in in one go.
        integer(int64) , intent(out) :: CurrWalkers    !Number of determinants which end up on a given processor.
        real(dp), intent(out) :: CurrHF(lenof_sign)
        integer, intent(inout) :: PopNel
        ! If true, then calculate the diagonal Hamiltonian elements of the
        ! read-in popsfile, and, if using the linear scaling algorithm,
        ! fill in the walker hash table.
        logical, intent(in) :: tCalcExtraInfo
        integer, intent(in) :: pops_nnodes
        integer(int64), intent(in) :: pops_walkers(0:nProcessors-1)

        real(dp) :: CurrParts(lenof_sign)
        integer :: Slot, nJ(nel)
        integer :: iunit,j,ierr,PopsInitialSlots(0:nNodes-1), iunit_3
        integer(int64) :: i
        real(dp) :: BatchSize
        integer :: PopsSendList(0:nNodes-1),proc
        integer :: MaxSendIndex,err,DetHash
        integer(n_int) :: WalkerTemp(0:NIfTot)
        integer(int64) :: Det, AllCurrWalkers
        logical :: FormPops,BinPops,tReadAllPops,tStoreDet
        real(dp) , dimension(lenof_sign) :: SignTemp
        integer :: PopsVersion
        character(len=*) , parameter :: this_routine='ReadFromPopsfile'
        HElement_t :: HElemTemp
        character(255) :: popsfile
        !variables from header file
        logical :: tPop64Bit, tPopHPHF, tPopLz, tEOF
        integer :: iPopLenof_sign, iPopIter, PopNIfD, PopNIfY, PopNIfSgn, PopNIfFlag
        integer :: PopNIfTot, PopBlockingIter, read_nnodes, Popinum_runs
        integer :: PopRandomHash(1024)
        integer(int64) :: iPopAllTotWalkers
        real(dp) :: PopDiagSft, PopDiagSft2, read_tau, read_psingles
        real(dp) :: read_pparallel
        real(dp) , dimension(lenof_sign/inum_runs) :: PopSumNoatHF
        integer(int64) :: read_walkers_on_nodes(0:nProcessors-1)
        integer, intent(in) :: DetsLen
        INTEGER(kind=n_int), intent(out) :: Dets(0:nIfTot,DetsLen)
        character(12) :: tmp_num
        character(255) :: tmp_char
        HElement_t :: PopAllSumENum
        integer :: sgn(lenof_sign), flg, part_on_node = 0
        type(ll_node), pointer :: Temp
        integer, allocatable :: TempnI(:)
        real(dp) :: pops_norm_temp
        type(timer), save :: read_timer, process_timer

        ! Tme the overall popsfile read in
        read_timer%timer_name = 'POPS-read'
        process_timer%timer_name = 'POPS-process'
        call set_timer(read_timer)

        call open_pops_head(iunit,formpops,binpops)
        ! Determine version number.
        PopsVersion=FindPopsfileVersion(iunit)
        IF(FormPops) THEN
            if(PopsVersion.eq.3) then
                call ReadPopsHeadv3(iunit,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,PopNel, &
                    iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter, &
                    PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot)
            else
                call ReadPopsHeadv4(iunit,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,PopNel, &
                    iPopAllTotWalkers,PopDiagSft,PopDiagSft2,PopSumNoatHF,PopAllSumENum,iPopIter, &
                    PopNIfD,PopNIfY,PopNIfSgn,Popinum_runs,PopNIfFlag,PopNIfTot,read_tau, &
                    PopBlockingIter, PopRandomHash, read_psingles, read_pparallel, &
                    read_nnodes, read_walkers_on_nodes)
            endif

            if(EndPopsList.ne.iPopAllTotWalkers) then
                call stop_all(this_routine,"Error in assessing number of entries in POPSFILE")
            endif

        else if (BinPops) then
            ! If we are reading a binary popsfile, then we want to close the
            ! header file and open the file containing the determinants.
            ! We don't need to bother reading in the header, it has already
            ! been done (as it has been for non-binary popsfiles, too).
            if (iProcIndex == root .or. (bNodeRoot .and. tSplitPops)) then
                ! Close the header file.
                close(iunit)

                ! Use the correct pops file name.
                tmp_char = 'POPSFILEBIN'
                if (tSplitPops) then
                    write(tmp_num, '(i12)') iProcIndex
                    tmp_char = trim(tmp_char) // '-' // trim(adjustl(tmp_num))
                end if
                call get_unique_filename(trim(tmp_char), tIncrementPops, &
                                         .false., iPopsFileNoRead, popsfile)
                open(iunit, file=popsfile, status='old', form='unformatted')
            endif

            ! We need to consider the same parameters for particle
            ! distribution
            read_nnodes = pops_nnodes
            read_walkers_on_nodes = pops_walkers
        end if

        allocate(TempnI(PopNel))

        call mpibarrier(err)

        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,es17.10,A)") "Although ",EndPopsList, &
                " configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            else
                write(6,"(A,I15,A)") "Reading in a total of ",EndPopsList, " configurations from POPSFILE."
            ENDIF
            if(ScaleWalkers.ne.1) then
                call warning_neci(this_routine,"ScaleWalkers parameter found, but not implemented in POPSFILE v3 - ignoring.")
            endif
            call neci_flush(6)
        ENDIF

        ! Which of the POPSFILE reading routines are we going to make use of?
        if (tSplitPops) then
            CurrWalkers = read_pops_splitpops (iunit, PopNel, TempnI, BinPops, &
                                               Dets, DetsLen, iunit_3, PopNIfSgn)
        else if (pops_nnodes == nProcessors .and. (.not. tCCMC) .and. PopsVersion == 4) then
            CurrWalkers = read_pops_nnodes (iunit, PopNel, TempnI, BinPops, Dets, &
                                            DetsLen, read_walkers_on_nodes, &
                                            iunit_3, PopNIfSgn)
        else
            CurrWalkers = read_pops_general (iunit, PopNel, TempnI, BinPops, Dets, &
                                             DetsLen, ReadBatch, EndPopsList, &
                                             iunit_3, PopNIfSgn)
        end if

        ! Add the contributions to the norm of the popsfile wave function from
        ! all processes.
        pops_norm_temp = pops_norm
        call MPISumAll(pops_norm_temp, pops_norm)

        ! Close the popsfiles.
        if(iProcIndex == Root .or. (tSplitPops .and. bNodeRoot)) then
            close(iunit)
        endif

        ! Test we have still got all determinants
        write(6,*) "CurrWalkers: ", CurrWalkers
        call MPISum(CurrWalkers, 1, AllCurrWalkers)
        if (iProcIndex == Root) then
            if (iWeightPopRead == 0 .and. AllCurrWalkers /= EndPopsList) then
                write(6,*) "AllCurrWalkers: ", AllCurrWalkers
                write(6,*) "EndPopsList: ", EndPopsList

                if (tSplitPops) then
                    write(6,*)
                    write(6,*) '*******************************************'
                    write(6,*) 'Currently using: ', nProcessors, ' processors'
                    write(6,*)
                    write(6,*) 'Using pre-split popsfiles (SPLIT-POPS, &
                               &POPSFILEBIN-*).'
                    write(6,*) 'Ensure that the number of processors and &
                               &number of POPSFILEs match.'
                    write(6,*) '*******************************************'
                end if

                call stop_All(this_routine, "Not all walkers accounted for &
                                            &when reading in")
            endif
        endif

        ! Clear all deterministic and trial flags so that they can be changed later.
        if (tUseFlags) then
            do i = 1, CurrWalkers
                call clr_flag(Dets(:,i), flag_deterministic)
                call clr_flag(Dets(:,i), flag_determ_parent)
                call clr_flag(Dets(:,i), flag_trial)
                call clr_flag(Dets(:,i), flag_connected)
            end do
        end if

        call halt_timer(read_timer)
        call set_timer(process_timer)

        if (.not. tHashWalkerList) call sort(dets(:,1:CurrWalkers), ilut_lt, ilut_gt)

        if (tCalcExtraInfo) then

            if (tHashWalkerList) then
                call clear_hash_table(HashIndex)
                call fill_in_hash_table(HashIndex, nWalkerHashes, Dets, int(CurrWalkers,sizeof_int), .true.)
            endif

            ! Run through all determinants on each node, and calculate the total number of walkers, and noathf
            CurrHF = 0.0_dp
            CurrParts = 0.0_dp
            do i = 1, CurrWalkers
                call extract_sign(Dets(:,i),SignTemp)
                
                CurrParts=CurrParts+abs(SignTemp)
                if(DetBitEQ(Dets(:,i),iLutRef,NIfDBO)) then
                    if(CurrHF(1).ne.0) then
                        call stop_all(this_routine,"HF already found, but shouldn't have")
                    endif
                    CurrHF=CurrHF+SignTemp 
                    if (.not. tSemiStochastic) call set_det_diagH(int(i, sizeof_int), 0.0_dp)
                else
                    if (.not. tSemiStochastic) then
                    ! Calculate diagonal matrix element
                        call decode_bit_det(TempnI, Dets(:,i))
                        if (tHPHF) then
                            HElemTemp = hphf_diag_helement(TempnI,Dets(:,i))
                        else
                            HElemTemp = get_helement (TempnI, TempnI, 0)
                        endif
                        call set_det_diagH(int(i, sizeof_int), real(HElemTemp, dp) - Hii)
                    endif
                endif
            enddo

        end if

        call mpibarrier(err)

        if (allocated(TempnI)) deallocate(TempnI)

        call halt_timer(process_timer)

    end subroutine ReadFromPopsfile

    function read_pops_nnodes (iunit, PopNel, unused, binary_pops, det_list, &
                               max_dets, read_walkers_on_nodes, iunit_3, &
                               PopNIfSgn) &
                              result(CurrWalkers)

        ! A routine to read in popsfiles, making use of the stored information
        ! in the POPSFILE to tell us where each of the particles should end up
        ! (so that we don't need to calculate it).
        !
        ! --> Should be more efficient on communication overhead as well.

        integer, intent(in) :: iunit, PopNel, max_dets, iunit_3, PopNIfSgn
        integer, intent(out) :: unused(PopNel)
        integer(int64), intent(in) :: read_walkers_on_nodes(0:nProcessors-1)
        integer(n_int), intent(out) :: det_list(0:NIfTot, max_dets)
        logical, intent(in) :: binary_pops
        integer(int64) :: CurrWalkers
        character(*), parameter :: this_routine = 'read_pops_nnodes'

        ! The buffer is able to store the maximum number of particles on any
        ! determinant.
        integer(n_int), allocatable :: buffer(:,:)
        integer :: ndets, det, ierr, nelem, proc
        logical :: tEOF

        integer :: i

        ! A tag is used to identify this send/recv pair over any others
        integer, parameter :: mpi_tag = 123456  !z'beef'
        integer, parameter :: mpi_tag_h = 123457

        ! Initialise counters
        CurrWalkers = 0
        pops_norm = 0.0_dp
        det = 1

        ! A quick test for the initialisation of the walker arrays
        if (any(read_walkers_on_nodes > max_dets)) &
            call stop_all(this_routine, "Insufficient particle storage &
                         &allocated to store particles in POPSFILE")

        ! If we are on the root processor, then we do the reading in. Otherwise
        ! we just need to wait to have particles sent in!
        if (iProcIndex == root) then

            ! We need to allocate this, not put it on the stack, otherwise
            ! when using ifort we will get segfaults. gfortran sensibly
            ! allocates things bigger than the stacksize on the heap...
            allocate(buffer(0:NIfTot, max_dets), stat=ierr)
            if (ierr /= 0) &
                call stop_all(this_routine, 'Allocation of read buffer failed')

            ! Similarly, allocate the buffer for currentH
            if (ierr /= 0) &
                call stop_all(this_routine, 'Allocation of H-buffer failed')

            do proc = 0, nProcessors - 1

                ndets = 0
                do while (ndets < read_walkers_on_nodes(proc))

                    ! Read and store a particle for transmission
                    ndets = ndets + 1
                    tEOF = read_popsfile_det (iunit, PopNel, binary_pops, &
                                              buffer(:, ndets), unused, &
                                              PopNIfSgn, iunit_3, .false.)

                    ! Add the contribution from this determinant to the
                    ! norm of the popsfile wave function.
                    call add_pops_norm_contrib(buffer(:, ndets))

                    ! Catch a premature End-Of-File
                    if (tEOF) call stop_all(this_routine, &
                                            "Too few determinants found.")

                end do

                if (proc == root) then

                    ! If we have just read in the particles for the root
                    ! processor, just store them.
                    det_list(:, 1:ndets) = buffer(:, 1:ndets)
                    CurrWalkers = ndets

                else

                    ! Send the particles to the relevant processor. The counts
                    ! have already been transmitted.
                    nelem = ndets * (1 + NIfTot)
                    call MPISend(buffer(:,1:ndets), nelem, proc, mpi_tag, ierr)

                end if
                
            end do

            ! And deallocate the buffer
            deallocate(buffer)

        else if (bNodeRoot) then

            ! And receive the dets!
            ndets = int(read_walkers_on_nodes(iProcIndex))
            nelem = ndets * (1 + NIfTot)
            call MPIRecv(det_list, nelem, root, mpi_tag, ierr)


            ! Now we know how many particles are on this node
            CurrWalkers = ndets

        end if

    end function
        
    function read_pops_splitpops (iunit, PopNel, det_tmp, binary_pops, &
                                  det_list, max_dets, iunit_3, PopNIfSgn) &
                                  result(CurrWalkers)

        ! A routine to read in split popsfiles to each of the nodes.
        !
        ! In: iunit       - The popsfile being read from
        !     binary_pops - Is this a binary popsfile?
        
        integer, intent(in) :: iunit, PopNel, max_dets, iunit_3, PopNIfSgn
        integer, intent(out) :: det_tmp(PopNel)
        integer(n_int), intent(out) :: det_list(0:NifTot, max_dets)
        logical, intent(in) :: binary_pops
        integer(int64) :: CurrWalkers
        character(*), parameter :: this_routine = 'read_pops_splitpops'

        integer(n_int), allocatable :: BatchRead(:,:)
        integer(n_int) :: ilut_tmp(0:NIfTot)
        logical :: tEOF
        integer :: det
        integer :: proc

        allocate(BatchRead(0:NifTot, 1:MaxWalkersPart))

        write(6,*) 'Reading a maximum of ', MaxWalkersPart, ' particles to &
                   &each node from split POPSFILES'

        ! Initialise the relevant counters
        CurrWalkers = 0
        pops_norm = 0.0_dp

        ! If we are using pre-split popsfiles, then we need to do the
        ! reading on all of the nodes.
        if (bNodeRoot) then

            ! Get ready for reading in the next batch of walkers
            do while (.true.)

                ! Read the next entry, and store the walker in ilut_tmp.
                ! The decoded form is placed in det_tmp
                CurrWalkers = CurrWalkers + 1
                tEOF = read_popsfile_det (iunit, PopNel, binary_pops, &
                                          det_list(:, CurrWalkers), &
                                          det_tmp, PopNIfSgn, iunit_3, &
                                          .true.)

                ! Add the contribution from this determinant to the
                ! norm of the popsfile wave function.
                call add_pops_norm_contrib(det_list(:, CurrWalkers))

                ! And store the current H-values

                ! When we have got to the end of the file, we are done.
                if (tEOF) exit

                ! And a test that this split popsfile is somewhat valid...
                proc = DetermineDetNode(PopNel, det_tmp, 0)
                if (proc /= iProcIndex) &
                    call stop_all (this_routine, "Determinant in the &
                                   &wrong Split POPSFILE")

            enddo

        endif

        deallocate(BatchRead)

    end function read_pops_splitpops


    function read_pops_general (iunit, PopNel, det_tmp, binary_pops, det_list, &
                                max_dets, ReadBatch, EndPopsList, iunit_3, &
                                PopNIfSgn) &
                               result(CurrWalkers)

        ! General all-purpose pops reading routine.
        !
        ! In: iunit     - The popsfile being read from
        !     ReadBatch - The size of the buffer array to use. Normally
        !                 MaxSpawned, unless specified manually.

        integer, intent(in) :: iunit, PopNel, ReadBatch, max_dets, iunit_3
        integer, intent(out) :: det_tmp(PopNel)
        integer, intent(in) :: PopNIfSgn
        logical, intent(in) :: binary_pops
        integer(n_int), intent(out) :: det_list(0:NIfTot, max_dets)
        integer(int64), intent(in) :: EndPopsList
        integer(int64) :: CurrWalkers
        character(*), parameter :: this_routine = 'read_pops_general'

        logical :: tEOF, tReadAllPops
        integer(MPIArg) :: sendcounts(nNodes), disps(nNodes), recvcount
        integer(MPIArg) :: sendcounts2(nNodes), disps2(nNodes), recvcount2
        integer :: PopsInitialSlots(0:nNodes-1), PopsSendList(0:nNodes-1)
        integer :: batch_size, MaxSendIndex, i, j, det, nBatches, err, proc
        integer(n_int) :: ilut_tmp(0:NIfTot)

        integer(n_int), allocatable :: BatchRead(:,:)

        ! If we are on the root processor, we need a buffer array which stores
        ! the particles as they are read in. 
        !
        ! If we are on another processor, allocate a one-element array to avoid
        ! segfaults in MPIScatter.
        allocate(BatchRead(0:NIfTot, merge(ReadBatch, 1, iProcIndex == root)))

        if (iProcIndex == root) then

            ! This is the size of the batches in the above array (only has
            ! meaning on the root processor)
            batch_size = int(real(ReadBatch, dp) / real(nNodes, dp))

            ! What is the first element in the buffer array that the particles
            ! for ! each processor are placed in?
            forall(i = 0 : nNodes - 1) PopsInitialSlots(i) = batch_size * i + 1

            write(6, '(a,i12,a)') "Reading in a maximum of ", ReadBatch, &
                                  " determinants at a time from POPSFILE'"
            call neci_flush(6)

        end if

        ! Keep reading until all of the particles have been read in!
        det = 1
        tReadAllPops = .false.
        CurrWalkers = 0
        pops_norm = 0.0_dp
        nBatches = 0
        sendcounts = 0
        disps = 0
        MaxSendIndex = 1
r_loop: do while (.not. tReadAllPops)

            ! We read in the particles on the root node.
            ! --> Only do particle receiving on teh other nodes
            if (iProcIndex == root) then

                ! Get ready for reading in the next batch of walkers
                nBatches = nBatches + 1
                BatchRead(:,:) = 0
                PopsSendList(:) = PopsInitialSlots(:)

                do while (Det <= EndPopsList .or. tSplitPops)

                    ! Read the next entry, and store the walker in ilut_tmp.
                    ! The decoded form is placed in det_tmp
                    det = det + 1
                    tEOF = read_popsfile_det (iunit, PopNel, binary_pops, &
                                              ilut_tmp, det_tmp, PopNIfSgn, &
                                              iunit_3, .true.)
                    ! Add the contribution from this determinant to the
                    ! norm of the popsfile wave function.
                    call add_pops_norm_contrib(ilut_tmp)

                    ! When we have got to the end of the file, we are done.
                    if (tEOF) call stop_all (this_routine, &
                                             "Too few determinants found.")

                    ! Where should this particle be going?
                    proc = DetermineDetNode(PopNel, det_tmp, 0)

                    ! Store the found determinant in the temporary list, 
                    ! and if we have filled up the slot in the list then
                    ! distribute it when it is full.
                    BatchRead(:,PopsSendList(proc)) = ilut_tmp(:)
                    PopsSendList(proc) = PopsSendList(proc) + 1

                    ! If we have filled up the lists, exit the loop so that the
                    ! particles get distributed
                    if(proc /= nNodes - 1) then
                        if (PopsInitialSlots(proc+1) - PopsSendList(proc) < 2)&
                            exit
                    else
                        if (ReadBatch - PopsSendList(proc) < 2) exit
                    endif
                enddo

                ! Have we read in all of the particles?
                if (det > EndPopsList) tReadAllPops = .true.

                ! Prep the counts for transmitting the particles to all of
                ! the nodes.
                do j = 0, nNodes - 1
                    sendcounts(j+1) = &
                        int((PopsSendList(j) - PopsInitialSlots(j)) * &
                            (NIfTot + 1), MPIArg)
                    disps(j+1) = &
                        int((PopsInitialSlots(j) - 1) * (NIfTot + 1), MPIArg)
                enddo
                MaxSendIndex = (disps(nNodes) + sendcounts(nNodes)) &
                             / (nIfTot + 1)

            endif

            ! Now scatter the particles read in to their correct processors.
            if (bNodeRoot) then

                ! How much data goes to each processor?
                call MPIScatter (sendcounts, recvcount, err, roots)
                if (err /= 0) &
                    call stop_all (this_routine, "MPI scatter error")

                ! Transfer the walkers.
                call MPIScatterV (BatchRead(:,1:MaxSendIndex), sendcounts, &
                                  disps, &
                                  det_list(:,CurrWalkers+1:CurrWalkers+1+(recvcount/(NIfTot+1))), &
                                  recvcount, err, Roots)
                if (err /= 0) &
                    call stop_all (this_routine, "MPI scatterV error")


                CurrWalkers = CurrWalkers + recvcount / (NIfTot + 1)

            end if

            ! Are we done?
            call MPIBCast(tReadAllPops)
        end do r_loop

        write(6,"(a,i8)") "Number of batches required to distribute all &
                          &determinants in POPSFILE: ", nBatches
        write(6,*) "Number of configurations read in to this process: ", &
                   CurrWalkers 

        deallocate(BatchRead)

    end function read_pops_general

    ! This routine reads the next determinant entry from a popsfile and stores
    ! it in WalkerTemp, ready to be distributed.
    !
    ! This should only be called from root node.
    ! iunit = unit
    ! BinPops = Binary popsfile or formatted
    ! WalkerTemp = Determinant entry returned
    function read_popsfile_det (iunit, nel_loc, BinPops, WalkerTemp, nI, &
                                PopNifSgn, iunit_3, decode_det) result(tEOF)

        integer, intent(in) :: iunit
        integer, intent(in) :: nel_loc
        integer(n_int), intent(out) :: WalkerTemp(0:NIfTot)
        integer, intent(out) :: nI(nel_loc)
        integer, intent(in) :: PopNifSgn
        integer, intent(in), optional :: iunit_3
        logical, intent(in) :: BinPops, decode_det
        integer(n_int) :: WalkerTemp2(0:NIfTot)
        integer(n_int) :: sgn_int(PopNifSgn)
        integer :: elec, flg, i, j, stat, k
        real(dp) :: sgn(PopNifSgn)
        real(dp) :: new_sgn(lenof_sign)
        integer(n_int) :: flg_read
        logical :: tStoreDet, tEOF

        WalkerTemp = 0_n_int
        tStoreDet=.false.
        tEOF = .false.
r_loop: do while(.not.tStoreDet)

            ! All basis parameters match --> Read in directly.
            if (tRealPOPSfile) then
                if (BinPops) then
                    if (tUseFlags) then
                        read(iunit, iostat=stat) WalkerTemp(0:NIfDBO), sgn,&
                                                 flg_read
                    else
                        read(iunit, iostat=stat) WalkerTemp(0:NIfDBO), sgn
                    end if
                else
                    if (tUseFlags) then
                        read(iunit,*, iostat=stat) WalkerTemp(0:NIfDBO), &
                                                   sgn, flg_read
                    else
                        read(iunit,*, iostat=stat) WalkerTemp(0:NIfDBO), sgn
                    end if
                end if
            else
                if (BinPops) then
                    if (tUseFlags) then
                        read(iunit, iostat=stat) WalkerTemp(0:NIfDBO), &
                                                 sgn_int, flg_read
                    else
                        read(iunit, iostat=stat) WalkerTemp(0:NIfDBO), &
                                                 sgn_int
                    end if
                else
                    if (tUseFlags) then
                        read(iunit,*, iostat=stat) WalkerTemp(0:NIfDBO), &
                                                   sgn_int, flg_read
                    else
                        read(iunit,*, iostat=stat) WalkerTemp(0:NIfDBO), &
                                                   sgn_int
                    end if
                end if

                sgn = sgn_int
            end if
            if (stat < 0) then
                tEOF = .true. ! End of file reached.
                exit r_loop
            end if
            
            if((inum_runs.eq.2).and.(PopNifSgn.eq.1)) then
                !Read in pops from a single run. Distribute an identical set of walkers to each walker set
                !and then propagate the two independently
                new_sgn(1)=sgn(1)
                new_sgn(inum_runs)=sgn(1)
            else
                do k=1,lenof_sign
                    new_sgn(k)=sgn(k)
                enddo
            endif

            ! Store the sign and flag information in the determinant.
            flg = int(flg_read,sizeof_int)

            call encode_sign (WalkerTemp, new_sgn)
            if (tUseFlags) call encode_flags (WalkerTemp, flg)

            if((inum_runs.eq.2).and.(PopNifSgn.eq.1)) then
                if (test_flag(WalkerTemp, flag_is_initiator(1))) then
                    call set_flag(WalkerTemp, flag_is_initiator(2))
                else
                    call clr_flag(WalkerTemp, flag_is_initiator(2))
                endif

                if (test_flag(WalkerTemp, flag_make_initiator(1))) then
                    call set_flag(WalkerTemp, flag_make_initiator(2)) 
                else
                    call clr_flag(WalkerTemp, flag_make_initiator(2)) 
                endif
            endif
            
            ! Test if we actually want to store this walker...
            if (any(abs(sgn) >= iWeightPopRead)) then
                tStoreDet = .true.
                exit
            end if
        enddo r_loop

        ! Decode the determinant as required
        if (.not. tEOF .and. decode_det) then
            call decode_bit_det (nI, WalkerTemp)
        endif

    end function read_popsfile_det

    subroutine read_popsfile_wrapper(perturbs)

        type(perturbation), intent(in), allocatable, optional :: perturbs(:)

        integer :: iunithead, PopsVersion
        ! Variables from popsfile header...
        integer :: iPopLenof_sign, PopNel, iPopIter, PopNIfD, PopNIfY, WalkerListSize
        integer :: PopNIfSgn, PopNIfFlag, PopNIfTot, PopBlockingIter, read_nnodes
        integer :: Popinum_runs
        integer :: PopRandomHash(1024)
        logical :: tPop64Bit, tPopHPHF, tPopLz, formpops, binpops
        integer(int64) :: iPopAllTotWalkers
        integer(int64) :: read_walkers_on_nodes(0:nProcessors-1)
        real(dp) :: PopDiagSft, PopDiagSft2, read_tau
        real(dp) :: read_psingles, read_pparallel
        real(dp), dimension(lenof_sign/inum_runs) :: PopSumNoatHF
        HElement_t :: PopAllSumENum
        integer :: perturb_ncreate, perturb_nannihilate

        character(len=*), parameter :: t_r = "read_popsfile_wrapper"

        ! Read the header.
        call open_pops_head(iunithead,formpops,binpops)

        PopsVersion = FindPopsfileVersion(iunithead)

        if(PopsVersion == 4) then
            call ReadPopsHeadv4(iunithead,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,PopNel, &
                    iPopAllTotWalkers,PopDiagSft,PopDiagSft2,PopSumNoatHF,PopAllSumENum,iPopIter, &
                    PopNIfD,PopNIfY,PopNIfSgn,Popinum_runs,PopNIfFlag,PopNIfTot, &
                    read_tau,PopBlockingIter, PopRandomHash, read_psingles, &
                    read_pparallel, read_nnodes, read_walkers_on_nodes)
        else
            call stop_all(t_r, "Only version 4 popsfile are supported with kp-fciqmc.")
        endif

        ! Check the number of electrons created and annihilated by the
        ! perturbation operators.
        if (present(perturbs)) then
            perturb_ncreate = perturbs(1)%ncreate
            perturb_nannihilate = perturbs(1)%nannihilate
        else
            perturb_ncreate = 0
            perturb_nannihilate = 0
        end if

        call CheckPopsParams(tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,PopNel, &
                iPopAllTotWalkers,PopDiagSft,PopDiagSft2,PopSumNoatHF,PopAllSumENum,iPopIter, &
                PopNIfD,PopNIfY,PopNIfSgn,Popinum_runs,PopNIfFlag,PopNIfTot, &
                WalkerListSize,read_tau,PopBlockingIter, read_psingles, read_pparallel, &
                perturb_ncreate, perturb_nannihilate)

        if (iProcIndex == root) close(iunithead)

        call InitFCIMC_pops(iPopAllTotWalkers, PopNIfSgn, PopNel, read_nnodes, &
                            read_walkers_on_nodes, perturbs)

        ! If requested output the norm of the *unperturbed* walkers in the POPSFILE.
        if (tWritePopsNorm) call write_pops_norm()

    end subroutine read_popsfile_wrapper

    subroutine InitFCIMC_pops(iPopAllTotWalkers, PopNIfSgn, PopNel, pops_nnodes, &
                              pops_walkers, perturbs)

        use CalcData, only : iReadWalkersRoot
        use hash, only: clear_hash_table, fill_in_hash_table
        use perturbations, only: apply_perturbation_array
        use semi_stoch_procs, only: fill_in_diag_helements

        integer(int64), intent(in) :: iPopAllTotWalkers
        integer, intent(in) :: PopNIfSgn
        ! Note that PopNel might be recalculated in ReadFromPopsfile (and might
        ! not have even been set on input).
        integer, intent(inout) :: PopNel
        integer, intent(in) :: pops_nnodes
        integer(int64), intent(in) :: pops_walkers(0:nProcessors-1)
        ! Perturbation operators to apply to the determinants after they have
        ! been read in.
        type(perturbation), intent(in), allocatable, optional :: perturbs(:)
        integer :: run, ReadBatch
        integer :: nI(nel)
        logical :: apply_pert
        integer :: TotWalkersIn

        if (iReadWalkersRoot == 0) then
            ! ReadBatch is the number of walkers to read in from the 
            ! popsfile at one time. The larger it is, the fewer
            ! communictions will be needed to scatter the particles.
            !
            ! By default, the new array (which is only created on the 
            ! root processors) is the same length as the spawning 
            ! arrays.
            ReadBatch = MaxSpawned
        else
            ReadBatch = iReadWalkersRoot
        end if

        apply_pert = .false.
        if (present(perturbs)) then
            if (allocated(perturbs)) apply_pert = .true.
         end if

        ! If applying perturbations, read the popsfile into the array
        ! popsfile_dets and then apply the perturbations to the determinants
        ! in this array.
        ! If not, then read the popsfile straight into CurrentDets.
        if (apply_pert) then
            call ReadFromPopsfile(iPopAllTotWalkers, ReadBatch, TotWalkers, TotParts, NoatHF, &
                                  popsfile_dets, MaxWalkersPart, pops_nnodes, pops_walkers, PopNIfSgn, &
                                  PopNel, tCalcExtraInfo=.false.)

            TotWalkersIn = int(TotWalkers, sizeof_int)
            call apply_perturbation_array(perturbs, TotWalkersIn, popsfile_dets, CurrentDets)
            TotWalkers = int(TotWalkersIn, int64)
        else
            call ReadFromPopsfile(iPopAllTotWalkers, ReadBatch, TotWalkers, TotParts, NoatHF, &
                                  CurrentDets, MaxWalkersPart, pops_nnodes, pops_walkers, PopNIfSgn, &
                                  PopNel, tCalcExtraInfo=.false.)
        end if

        call fill_in_diag_helements()

        if (tHashWalkerList) then
            call clear_hash_table(HashIndex)
            call fill_in_hash_table(HashIndex, nWalkerHashes, CurrentDets, int(TotWalkers, sizeof_int), .true.)
        end if

        call set_initial_global_data(TotWalkers, CurrentDets)

    end subroutine InitFCIMC_pops
    
    subroutine CheckPopsParams(tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                    iPopAllTotWalkers,PopDiagSft,PopDiagSft2,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                    PopNIfD,PopNIfY,PopNIfSgn,Popinum_runs,PopNIfFlag,PopNIfTot, &
                    WalkerListSize,read_tau,PopBlockingIter, read_psingles, &
                    read_pparallel, perturb_ncreate, perturb_nann)
        use LoggingData , only : tZeroProjE
        logical , intent(in) :: tPop64Bit,tPopHPHF,tPopLz
        integer , intent(in) :: iPopLenof_sign,iPopNel,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot,Popinum_runs
        integer , intent(in) :: PopBlockingIter
        integer(int64) , intent(in) :: iPopAllTotWalkers
        real(dp) , intent(in) :: PopDiagSft,PopDiagSft2,read_tau
        real(dp), intent(in) :: read_psingles, read_pparallel
        real(dp) , dimension(lenof_sign/inum_runs) , intent(in) :: PopSumNoatHF
        HElement_t , intent(in) :: PopAllSumENum
        integer, intent(in) :: perturb_ncreate, perturb_nann
        integer , intent(out) :: WalkerListSize
        character(len=*) , parameter :: this_routine='CheckPopsParams'
        integer :: k

        !Ensure all NIF and symmetry options the same as when popsfile was written out.
#ifdef __INT64
        if(.not.tPop64Bit) call stop_all(this_routine,"Popsfile created with 32 bit walkers, but now using 64 bit.")
#else
        if(tPop64Bit) call stop_all(this_routine,"Popsfile created with 64 bit walkers, but now using 32 bit.")
#endif
        if(tPopHPHF.neqv.tHPHF) call stop_all(this_routine,"Popsfile HPHF and input HPHF not same")
        if(tPopLz.neqv.tFixLz) call stop_all(this_routine,"Popsfile Lz and input Lz not same")
        if(iPopNEl+perturb_ncreate-perturb_nann.ne.NEl) call stop_all(this_routine,"The number of electrons &
            &in the POPSFILE is not consistent with the number you have asked to run with.")
        if(PopNIfD.ne.NIfD) call stop_all(this_routine,"Popsfile NIfD and calculated NIfD not same")
        if(PopNIfY.ne.NIfY) call stop_all(this_routine,"Popsfile NIfY and calculated NIfY not same")
        if(inum_runs.eq.1) then
            !We allow these two values to be different if we're reading in a popsfile fine inum_runs=1 and we want to
            !continue the calculation with inum_runs=2
            if(iPopLenof_sign.ne.lenof_sign) call stop_all(this_routine,"Popsfile lenof_sign and input lenof_sign not same")
            if(PopNIfSgn.ne.NIfSgn) call stop_all(this_routine,"Popsfile NIfSgn and calculated NIfSgn not same")
        endif

        ! This test is no longer needed. NIfFlag depends on how we represent
        ! the flags in memory, PopsNIfFlag depends on tUseFlags. The are 
        ! allowed to differ.
!        if(PopNIfFlag.ne.NIfFlag) call stop_all(this_routine,"Popsfile NIfFlag and calculated NIfFlag not same")
        if (inum_runs.eq.1) then
            if (tUseFlags .and. NIfFlag == 0) then
                if (PopNIFTot /= NIfTot + 1) &
                    call stop_all(this_routine, "Popsfile NIfTot and &
                                 &calculated NIfTot don't match.")
            else
                if (PopNIfTot /= NIfTot) &
                    call stop_all(this_routine,"Popsfile NIfTot and calculated&
                                               & NIfTot not same")
            end if
        endif


        IF(.not.tWalkContGrow) THEN
!If we want the walker number to be stable, take the shift from the POPSFILE, otherwise, keep the input value.
            if (inum_runs.eq.2) then
                if(Popinum_runs.eq.2) then
                    DiagSft(1)=PopDiagSft
                    DiagSft(inum_runs)=PopDiagSft2
                else
                    !Previously we only had a single run, now we are restarting with double run
                    DiagSft(1)=PopDiagSft
                    DiagSft(inum_runs)=PopDiagSft
                endif
            else
                DiagSft=PopDiagSft
            endif
        ENDIF

        if(PopDiagSft.eq.0.0_dp) then
            !If the popsfile has a shift of zero, continue letting the population grow
            tWalkContGrow=.true.
            if (inum_runs.eq.2) then
                DiagSft(1)=PopDiagSft
                DiagSft(inum_runs)=PopDiagSft
            else
                DiagSft=PopDiagSft
            endif
        endif

        if(tWalkContGrow) then
            tSinglePartPhase=.true.
            !If continuing to grow, ensure we can allocate enough memory for what we hope to get the walker population to,
            !rather than the average number of determinants in the popsfile.
            WalkerListSize=int(max(int(initwalkers,int64),NINT(real(iPopAllTotWalkers,dp)/real(nNodes,dp),int64)),sizeof_int)
        else
            tSinglePartPhase=.false.
            WalkerListSize=NINT(real(iPopAllTotWalkers,dp)/real(nNodes,dp))
        endif

        if(inum_runs.eq.2) then
            AllSumNoatHF(1)=PopSumNoatHF(1)
            AllSumNoatHF(inum_runs)=PopSumNoatHF(1)
            AllSumENum(1)=PopAllSumENum
            AllSumENum(inum_runs)=PopAllSumENum
        elseif(lenof_sign.eq.2) then
            AllSumNoatHF(1)=PopSumNoatHF(1)
            AllSumNoatHF(lenof_sign)=PopSumNoatHF(lenof_sign)
            AllSumENum=PopAllSumENum
        else
            AllSumNoatHF(1)=PopSumNoatHF(1)
            AllSumENum=PopAllSumENum
        endif

        
        PreviousCycles=iPopIter
        WRITE(6,*) "Number of iterations in previous simulation: ",PreviousCycles
        IF(NEquilSteps.gt.0) THEN
            WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
            NEquilSteps=0
        ENDIF
        IF(TZeroProjE) THEN
            !Reset energy estimator
            WRITE(6,*) "Resetting projected energy counters to zero..."
            AllSumENum=0.0_dp
            AllSumNoatHF = 0
        ENDIF
        if(read_tau.eq.0.0_dp) then
            !Using popsfile v.3, where tau is not written out.
            !Exit if trying to dynamically search for timestep
            if(tSearchTau) then
                call stop_all(this_routine,"Cannot dynamically search for timestep if reading &
                    &in POPSFILE v.3. Manually specify timestep.")
            endif
            write(6,*) "Old popsfile detected."
            write(6,*) "Therefore automatic blocking will only start from current run"
            iBlockingIter = PreviousCycles
        else
            !Using popsfile v.4, where tau is written out and read in
            if(tSearchTau) then
                if((.not.tSinglePartPhase(1)).or.(.not.tSinglePartPhase(inum_runs))) then
                    tSearchTau=.false.
                endif
                Tau=read_tau
                write(6,"(A)") "Using timestep specified in POPSFILE, although continuing to dynamically adjust to optimise this"

                ! If we have been searching for tau, we may have been searching
                ! for psingles (it is done at the same time).
                if (read_psingles /= 0) then
                    if (tCSF) then ! .or. tSpinProjDets) then
                        call stop_all(this_routine, "pSingles storage not yet &
                                      &implemented for CSFs")
                    end if
                    pSingles = read_psingles
                    pDoubles = 1.0_dp - pSingles
                end if

                if (read_pparallel /= 0) then
                    pParallel = read_pparallel
                end if
            else
                !Tau specified. if it is different, write this here.
                if(abs(read_tau-Tau).gt.1.0e-5_dp) then
                    call warning_neci(this_routine,"Timestep specified in input file is different to that in the popsfile.")
                    write(6,"(A,F12.8)") "Old timestep: ",read_tau
                    write(6,"(A,F12.8)") "New timestep: ",tau
                endif
            endif
            iBlockingIter = PopBlockingIter
        endif
    
    end subroutine CheckPopsParams


!Routine for reading in from iunit the header information from a popsile v3 file.
    subroutine ReadPopsHeadv3(iunithead,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                iPopAllTotWalkers,PopDiagSft,PopSumNoatHF,PopAllSumENum,iPopIter,   &
                PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot)
        integer , intent(in) :: iunithead
        logical , intent(out) :: tPop64Bit,tPopHPHF,tPopLz
        integer , intent(out) :: iPopLenof_sign,iPopNel,iPopIter,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot
        integer(int64) , intent(out) :: iPopAllTotWalkers
        real(dp) , intent(out) :: PopDiagSft
        real(dp) , dimension(lenof_sign/inum_runs) , intent(out) :: PopSumNoatHF
        HElement_t, intent(out) :: PopAllSumENum
        character(len=24) :: junk,junk2,junk3,junk4,junk5
        integer :: PopsVersion

        PopsVersion=FindPopsfileVersion(iunithead)
        if(PopsVersion.ne.3) call stop_all("ReadPopsfileHeadv3","Wrong popsfile version for this routine.")
            
        if(iProcIndex.eq.root) then
            read(iunithead,'(A12,L5,A8,L5,A8,L5,A13,I5,A7,I6)') junk,tPop64Bit,junk2,tPopHPHF,junk3, &
                tPopLz,junk4,iPopLenof_sign,junk5,iPopNEl
            read(iunithead,*) iPopAllTotWalkers
            read(iunithead,*) PopDiagSft 
            read(iunithead,*) PopSumNoatHF 
            read(iunithead,*) PopAllSumENum 
            read(iunithead,*) iPopIter 
            read(iunithead,*) PopNIfD 
            read(iunithead,*) PopNIfY 
            read(iunithead,*) PopNIfSgn 
            read(iunithead,*) PopNIfFlag 
            read(iunithead,*) PopNIfTot
        endif
        !Broadcast the read in values from the header to all nodes.
        call MPIBCast(tPop64Bit)
        call MPIBCast(tPopHPHF)
        call MPIBCast(tPopLz)
        call MPIBCast(iPopLenof_sign)
        call MPIBCast(iPopNEl)
        call MPIBCast(iPopAllTotWalkers)
        call MPIBCast(PopDiagSft)
        call MPIBCast(PopSumNoatHF)
        call MPIBCast(PopAllSumENum)
        call MPIBCast(iPopIter)
        call MPIBCast(PopNIfD)
        call MPIBCast(PopNIfY)
        call MPIBCast(PopNIfSgn)
        call MPIBCast(PopNIfFlag)
        call MPIBCast(PopNIfTot)

    end subroutine ReadPopsHeadv3

    !Routine for reading in from iunit the header information from a popsile v4 file.
    subroutine ReadPopsHeadv4(iunithead,tPop64Bit,tPopHPHF,tPopLz,iPopLenof_Sign,iPopNel, &
                iPopAllTotWalkers,PopDiagSft,PopDiagSft2,PopSumNoatHF_out,PopAllSumENum,iPopIter,   &
                PopNIfD,PopNIfY,PopNIfSgn,Popinum_runs,PopNIfFlag,PopNIfTot,read_tau, &
                PopBlockingIter, PopRandomHash, read_psingles, read_pparallel, &
                read_nnodes, read_walkers_on_nodes)
        integer , intent(in) :: iunithead
        logical , intent(out) :: tPop64Bit,tPopHPHF,tPopLz
        integer, intent(out) :: iPopLenof_sign, iPopNel, iPopIter, PopNIfD
        integer, intent(out) :: PopNIfY, PopNIfSgn, PopNIfFlag, PopNIfTot
        integer, intent(out) :: PopBlockingIter, read_nnodes, Popinum_runs
        integer, intent(out) :: PopRandomHash(1024)
        integer(int64), intent(out) :: read_walkers_on_nodes(0:nProcessors-1)
        integer(int64) , intent(out) :: iPopAllTotWalkers
        real(dp) , intent(out) :: PopDiagSft,PopDiagSft2, read_tau, read_psingles
        real(dp), intent(out) :: PopSumNoatHF_out(lenof_sign/inum_runs)
        real(dp) :: PopSumNoatHF(1024)
        real(dp), intent(out) :: read_pparallel
        HElement_t , intent(out) :: PopAllSumENum
        integer :: PopsVersion
        !Variables for the namelist
        logical :: Pop64Bit,PopHPHF,PopLz
        integer :: PopLensign,PopNEl,PopCyc,PopiBlockingIter
        integer, parameter :: max_nodes = 30000
        integer(int64) :: PopTotwalk, PopWalkersOnNodes(max_nodes)
        integer :: PopNNodes
        real(dp) :: PopSft, PopTau, PopPSingles, PopPParallel, PopGammaSing
        real(dp) :: PopGammaDoub, PopGammaOpp, PopGammaPar, PopMaxDeathCpt
        real(dp) :: PopTotImagTime, PopSft2, PopParBias
        character(*), parameter :: t_r = 'ReadPopsHeadv4'
        HElement_t :: PopSumENum
        namelist /POPSHEAD/ Pop64Bit,PopHPHF,PopLz,PopLensign,PopNEl, &
                    PopTotwalk,PopSft,PopSft2,PopSumNoatHF,PopSumENum, &
                    PopCyc,PopNIfD,PopNIfY,PopNIfSgn,PopNIfFlag,PopNIfTot, &
                    PopTau,PopiBlockingIter,PopRandomHash,PopPSingles, &
                    PopPParallel, PopNNodes, PopWalkersOnNodes, PopGammaSing, &
                    PopGammaDoub, PopGammaOpp, PopGammaPar, PopMaxDeathCpt, &
                    PopTotImagTime, Popinum_runs, PopParBias

        PopsVersion=FindPopsfileVersion(iunithead)
        if(PopsVersion.ne.4) call stop_all("ReadPopsfileHeadv4","Wrong popsfile version for this routine.")

        PopParBias = 0.0_dp
        PopPParallel = 0.0_dp
        if(iProcIndex.eq.root) then
            read(iunithead,POPSHEAD)
        endif
        !Broadcast the read in values from the header to all nodes.
        call MPIBCast(Pop64Bit)
        call MPIBCast(PopHPHF)
        call MPIBCast(PopLz)
        call MPIBCast(PopLensign)
        call MPIBCast(PopNEl)
        call MPIBCast(PopTotwalk)
        call MPIBCast(PopSft)
        call MPIBCast(PopSft2)
        PopSumNoatHF_out = PopSumNoatHF(1:lenof_sign/inum_runs)
        call MPIBCast(PopSumNoatHF_out)
        call MPIBCast(PopSumENum)
        call MPIBCast(PopCyc)
        call MPIBCast(PopNIfD)
        call MPIBCast(PopNIfY)
        call MPIBCast(PopNIfSgn)
        call MPIBCast(PopNIfFlag)
        call MPIBCast(PopNIfTot)
        call MPIBCast(Popinum_runs)
        call MPIBCast(PopTau)
        call MPIBCast(PopiBlockingIter)
        call MPIBCast(PopPSingles)
        call MPIBCast(PopPParallel)
        call MPIBCast(PopParBias)
        call MPIBCast(PopNNodes)
        call MPIBcast(PopTotImagTime)
        if (PopNNodes == nProcessors) then
            ! What is the maximum number of nodes currently supported. We might
            ! need to update this...
            if (PopNNodes > max_nodes) &
                call stop_all(t_r, "Too many processors in POPSFILE. Update &
                                   &max_nodes")

            call MPIBCast(PopWalkersOnNodes(1:PopNNodes))
            read_walkers_on_nodes(0:PopNNodes-1) = &
                PopWalkersOnNodes(1:PopNNodes)
        end if
        call MPIBcast(PopGammaSing)
        call MPIBcast(PopGammaDoub)
        call MPIBcast(PopGammaOpp)
        call MPIBCast(PopGammaPar)
        call MPIBcast(PopMaxDeathCpt)
        call MPIBcast(PopRandomHash)
        tPop64Bit=Pop64Bit
        tPopHPHF=PopHPHF
        tPopLz=PopLz
        iPopLenof_sign=PopLensign
        iPopNel=PopNel
        iPopAllTotWalkers=PopTotwalk
        PopDiagSft=PopSft
        PopDiagSft2=PopSft2
        PopAllSumENum=PopSumENum
        iPopIter=PopCyc
        read_tau=PopTau 
        PopBlockingIter=PopiBlockingIter
        read_psingles = PopPSingles
        if (PopParBias /= 0.0_dp .and. PopPParallel == 0.0_dp) then
            popPParallel = (PopParBias * par_elec_pairs) &
                         / (PopParBias * par_elec_pairs + AB_elec_pairs)
        end if
        read_pParallel = PopPParallel
        read_nnodes = PopNNodes
        TotImagTime = PopTotImagTime

        ! Fill the tau-searching accumulators, to avoid blips in tau etc.
        gamma_sing = PopGammaSing
        gamma_doub = PopGammaDoub
        gamma_opp = PopGammaOpp
        gamma_par = PopGammaPar
        max_death_cpt = PopMaxDeathCpt

    end subroutine ReadPopsHeadv4
    
    !NOTE: This should only be used for the v3 POPSFILEs, since we only open the POPSFILE on the head node.
    subroutine open_pops_head(iunithead,formpops,binpops)
        integer , intent(out) :: iunithead
        logical , intent(out) :: formpops,binpops
        character(255) :: popsfile

        if(iProcIndex.eq.root) then
            iunithead=get_free_unit()
            call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            inquire(file=popsfile,exist=formpops)
            if(formpops) then
                open(iunithead,file=popsfile,status='old')
                binpops=.false.
            else
                ! If we are using split popsfiles, the filenames are a bit
                ! different.
                if (tSplitPops) then
                    call get_unique_filename ('POPSFILEBIN-0', tIncrementPops, &
                                              .false., iPopsFileNoRead, &
                                              popsfile)
                else
                    call get_unique_filename ('POPSFILEBIN', tIncrementPops, &
                                              .false., iPopsFileNoRead, &
                                              popsfile)
                end if
                inquire(file=popsfile,exist=binpops)
                if(binpops) then
                    call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                    open(iunithead,file=popsfile,status='old')
                else 
                    call stop_all("open_pops_head","No POPSFILEs detected...")
                endif
            endif
            rewind(iunithead)
        endif
        call MPIBCast(binpops) 
        call MPIBCast(formpops) 

    end subroutine open_pops_head

!Return the version number of the popsfile
    integer function FindPopsfileVersion(iunithead)
        integer, intent(in) :: iunithead
        logical :: formpops, binpops, tEOF
        integer :: stat
        character(255) :: FirstLine

        ! Default value
        tRealPOPSfile = .true.

        if(iProcIndex.eq.root) then
            rewind(iunithead)
            read(iunithead,'(a255)') FirstLine

            if(index(FirstLine,'VERSION').eq.0) then
                FindPopsfileVersion=1
            else
                rewind(iunithead)
                read(iunithead,*) FirstLine,FirstLine,FirstLine,FindPopsfileVersion
            endif

            ! We need to be able to deal with popsfiles created with the 
            ! (old) integer version of the code
            ! --> No direct option was included for the output files to
            !     indicate if they used integers or real coefficients
            ! --> We need to take a (slight) guess...
            if (FindPopsfileVersion == 4 .and. tBinPops) then
                do while (.true.)
                    ! Read until the end of the file.
                    read(iunithead, '(a255)', iostat=stat) FirstLine
                    if (stat < 0) exit

                    ! If we have the line with PopSumNoatHF in it, check if
                    ! the reported number is a real number. If it is, we know
                    ! the popsfile was created by the realcoeff branch.
                    if (index(FirstLine, 'PopSumNoatHF=') /= 0) then
                        if (index(FirstLine, '.') /= 0) then
                            tRealPOPSfile = .true.
                        else
                            tRealPOPSfile = .false.
                        end if
                        exit
                    end if
                end do

                ! Rewind to allow normal reading of the header.
                rewind(iunithead)
                read(iunithead,*) FirstLine,FirstLine,FirstLine,FindPopsfileVersion
            end if

        endif
        call MPIBCast(FindPopsfileVersion)
        call MPIBcast(tRealPOPSfile)

    end function FindPopsfileVersion


!This routine is the same as WriteToPopsfilePar, but does not require two 
! main arrays to hold the data.
!The root processors data will be stored in a temporary array while it 
! recieves the data from the other processors.
!This routine will write out to a popsfile. It transfers all walkers to the 
! head node sequentially, so does not want to be called too often
    SUBROUTINE WriteToPopsfileParOneArr(Dets,nDets)
        use constants, only: size_n_int,n_int
        use CalcData, only: iPopsFileNoWrite, InitiatorWalkNo
        use MemoryManager, only: TagIntType
        integer(int64),intent(in) :: nDets !The number of occupied entries in Dets
        integer(kind=n_int),intent(in) :: Dets(0:nIfTot,1:nDets)
        INTEGER :: error
        integer(int64) :: WalkersonNodes(0:nNodes-1),writeoutdet
        INTEGER :: Tag, Tag2
        INTEGER :: Total,i,j,k
        INTEGER(KIND=n_int), ALLOCATABLE :: Parts(:,:)
        INTEGER(TagIntType) :: PartsTag=0
        integer :: nMaxDets, TempDet(0:NIfTot), TempFlags
        integer :: iunit, iunit_2, iunit_3, Initiator_Count, nwrite
        integer(int64) :: write_count, write_count_sum
        CHARACTER(len=*) , PARAMETER :: this_routine='WriteToPopsfileParOneArr'
        character(255) :: popsfile
        real(dp) :: TempSign(lenof_sign)
        integer :: excit_lev
        character(1024) :: out_tmp
        character(12) :: num_tmp
        type(timer), save :: write_timer

        ! Tme the overall popsfile read in
        write_timer%timer_name = 'POPS-write'
        call set_timer(write_timer)

        CALL MPIBarrier(error)  !sync
!        WRITE(6,*) "Get Here",nDets
!        CALL neci_flush(6)

!First, make sure we have up-to-date information - again collect AllTotWalkers
! ,AllSumNoatHF and AllSumENum...
!Calculate the energy by summing all on HF and doubles - convert number at HF
!  to a real since no int*8 MPI data type
        CALL MPISum(SumNoatHF,1,AllSumNoatHF)
        CALL MPISum(SumENum,1,AllSumENum)

!We also need to tell the root processor how many particles to expect from 
!each node - these are gathered into WalkersonNodes
        if(bNodeRoot) CALL MPIAllGather(nDets,WalkersonNodes,error,Roots)

        Tag=125
        Tag2=126

!We have to make the distinction here between the number of entries to expect,
!and the number of determinants we are writing out. Since the list is not
!necessarily contiguous any more, we have to calculate Alltotwalkers seperately.
        if(tHashwalkerlist) then
            Writeoutdet=0
            do i=1,int(nDets,sizeof_int)
                call extract_sign(Dets(:,i),TempSign)
                if(.not.IsUnoccDet(TempSign)) then
                    !Count this det in AllTotWalkers
                    Writeoutdet=Writeoutdet+1
                endif
            enddo
            writeoutdet=int(writeoutdet/iPopsPartEvery)
            call mpisum(writeoutdet,1,AllTotWalkers)
        else
            if(iProcIndex.eq.Root) then
                AllTotWalkers=0
                do i=0,nNodes-1
                    AllTotWalkers=AllTotWalkers+INT(WalkersonNodes(i)/iPopsPartEvery)
                enddo
            endif
        endif

        ! We only want to be using space/speed saving devices if we are doing
        ! them to the best of our ability.
        if (tSplitPops .and. .not. tBinPops) then
            call stop_all (this_routine, "Split pops files only works with &
                                         &binary pops files")
        end if

        if (iProcIndex == root) then

            ! Construct an output string to give feedback about what is 
            ! being done.
            write(6,*)
            write(6,'("*********************************")')
            write(out_tmp, '("Writing a ", i2, "-bit")') bits_n_int
            if (iPopsPartEvery /= 1) out_tmp = trim(out_tmp) // " reduced"
            out_tmp = trim(out_tmp) // " POPSFILE"
            if (tBinPops) out_tmp = trim(out_tmp) // "BIN"
            out_tmp = trim(out_tmp) // "..."
            write(6, '(a)') trim(adjustl(out_tmp))

            ! If we aren't outputting all of the walkers, indicate how many
            ! we are going to output.
            !if (iPopsPartEvery /= 1) then
                write(num_tmp, '(i12)') AllTotWalkers
                write(6, '("Writing a total of ",a," determinants.")') &
                    trim(adjustl(num_tmp))
            !end if

            ! If we are only outputting determinants with a certain weight
            if (tBinPops .and. binarypops_min_weight /= 0) then
                write(num_tmp, '(f12.3)') binarypops_min_weight
                write(6, '("Only outputting determinants holding >= ",a,&
                           &" walkers")') trim(adjustl(num_tmp))
            end if
            write(6,'("*********************************")')
            write(6,*)

            ! With a normal popsfile, the header is written at the start.
            ! Thus we need to do that now.
            if (.not. tBinPops) then
                call get_unique_filename('POPSFILE', tIncrementPops, .true., &
                                         iPopsFileNoWrite, popsfile)
                iunit = get_free_unit()
                ! We set recl=50000, which allows the line length written to be
                ! up to 50000 characters long. This allows popsfiles to be
                ! written in jobs with up to around 2940 MPI processes (the
                ! default value caused crashes when using over 2000 processes).
                open(iunit, file=popsfile, status='replace', recl=50000)
                call write_popsfile_header (iunit, AllTotWalkers, &
                                            WalkersonNodes)
            end if

        end if

        write_count = 0
        write_count_sum = 0
        if ((tSplitPops .and. bNodeRoot) .or. iProcIndex == root) then

            ! For a binary popsfile, the actual data is stored more
            ! compactly.
            if (tBinPops) then
                if (tSplitPops) then
                    write(num_tmp, '(i12)') iProcIndex
                    out_tmp = 'POPSFILEBIN-' // trim(adjustl(num_tmp))
                else
                    out_tmp = 'POPSFILEBIN'
                end if

                ! We need to get another unit, as not all nodes will go
                ! through the header-output section.
                iunit = get_free_unit()
                call get_unique_filename(trim(out_tmp), tIncrementPops, &
                                         .true., iPopsFileNoWrite, popsfile)
                open(iunit, file=popsfile, status='replace', &
                     form='unformatted')
            end if

            ! Are we outputting initiators as well?
            if (tPrintInitiators) then
                iunit_2 = get_free_unit()
                if (tSplitPops) then
                    write(num_tmp, '(i12)') iProcIndex
                    out_tmp = 'INITIATORS-' // trim(adjustl(num_tmp))
                else
                    out_tmp = 'INITIATORS'
                end if
                open(iunit_2, file=trim(out_tmp), status='UNKNOWN')
                Initiator_Count = 0
            end if

            !if(tRDMonFly.and.(.not.tExplicitAllRDM)) then
            !iunit_3 = get_free_unit()
            ! 
            !if (tSplitPops) then
            !        write(num_tmp, '(i12)') iProcIndex
            !        out_tmp = 'RDM_AV_POP-' // trim(adjustl(num_tmp))
            !    else
            !        out_tmp = 'RDM_AV_POP'
            !    end if
            !    open(iunit_3, file=trim(out_tmp), status='UNKNOWN', form='unformatted')
            !endif

            ! Write out the dets from this node (the head node, unless we
            ! are in split-pops mode.
            do j = 1, int(ndets, sizeof_int)
                ! Count the number of written particles
                if (write_pops_det (iunit, iunit_2, Dets(:,j), j)) then
                    !if(tRDMonFly.and.(.not.tExplicitAllRDM)) then
                    !    write(iunit_3) CurrentH(1:1+2*lenof_sign,j)
                    !endif
                    write_count = write_count + 1
                else
                    ! Zero determinants don't get written to binary popsfiles
                    ! --> Adjust the counts to deal with this.
                    WalkersOnNodes(0) = WalkersOnNodes(0) - 1
                end if
            end do

            if (.not. tSplitPops) then
                ! Allocate a temporary array to store the data being received
                ! from the other nodes.
                ! BEWARE - this is a potential crash point with a big 
                !          calculation.
                ! TODO: If this is the end of the run (say on a big machine),
                !       should we have deallocated spawnedparts by here to 
                !       ensure we have room, or does the deallocated space from
                !       dealing with freezing give us plenty of room?
                nMaxDets = int(maxval(WalkersOnNodes), sizeof_int)
                allocate(Parts(0:NIfTot, nMaxDets), stat=error)

                call LogMemAlloc ('Parts', int(nMaxDets,int32)*(NIfTot+1), &
                                  size_n_int, this_routine, PartsTag, error)
                
                ! Loop over the other nodes in the system sequentially, receive
                ! their walker lists, and output to the popsfiles.
                do i = 1, nNodes - 1

                    ! Calculate the amount of data to receive, and get it.
                    j = int(WalkersonNodes(i), sizeof_int) * (NIfTot+1)
                    call MPIRecv (Parts(:, 1:WalkersonNodes(i)), j, &
                                  NodeRoots(i), Tag, error)
                    
                    ! SDS: Catherine seems to have disabled writing these out
                    !      so no need to communicate them.
                    !!!if(tRDMonFly.and.(.not.tExplicitAllRDM)) then
                    !!!    ! And now for CurrentH
                    !!!    j = int(WalkersonNodes(i), sizeof_int) * (1+2*lenof_sign)
                    !!!    call MPIRecv (AllCurrentH(:, 1:WalkersonNodes(i)), j, &
                    !!!                  NodeRoots(i), Tag2, error)
                    !!!endif

                    ! Then write it out in the same way as above.
                    nwrite = int(WalkersOnNodes(i), sizeof_int)
                    do j = 1, nwrite
                        if (write_pops_det(iunit, iunit_2, Parts(:,j), j)) then
                            !if(tRDMonFly.and.(.not.tExplicitAllRDM)) then
                            !    write(iunit_3) AllCurrentH(1:1+2*lenof_sign,j)
                            !endif
                            write_count = write_count + 1
                        else
                            ! We have found a zero determinant. This doesn't
                            ! get added to the binary popsfile, so adjust
                            ! the count here.
                            WalkersOnNodes(i) = WalkersOnNodes(i) - 1
                        end if
                    end do

                end do

            end if

            ! Close the output files
            close(iunit)
            if (tPrintInitiators) close(iunit_2)

        else if (bNodeRoot) then

            ! Send the data on this processor to the root node for outputting
            ! in a combined popsfile
            ASSERT(.not. tSplitPops)
            j = int(nDets, sizeof_int) * (NIfTot + 1)
            call MPISend (Dets(0:NIfTot, 1:nDets), j, root, Tag, error)
            
            !!!if(tRDMonFly.and.(.not.tExplicitAllRDM)) then
            !!!    j = int(nDets, sizeof_int) * (1+2*lenof_sign)
            !!!    call MPISend (CurrentH(1:1+2*lenof_sign, 1:nDets), j, root, Tag2, error)
            !!!endif

        end if



        ! With binary popsfiles, the header is written once the popsfiles
        ! have been created, so that we can store the number of particles
        ! actually written, rather than the total within the system.
        ! Note that for non-binary popsfiles, all particles are written.
        if (tBinPops) then

            ! Get a count of the number of particles written
            call MPISum(write_count, write_count_sum)
            if (binarypops_min_weight == 0 .and. &
                    write_count_sum /= AllTotwalkers) then
                write(6,*) 'WARNING: Number of particles written (', &
                    write_count_sum, ') does not equal AllTotWalkers (', &
                    AllTotWalkers, ')'
            end if

            ! Now we need to actually write the header
            if (iProcIndex == root) then
                call get_unique_filename('POPSFILEHEAD', tIncrementPops, &
                                         .true., iPopsFileNoWrite, popsfile)
                iunit = get_free_unit()
                ! We set recl=50000, which allows the line length written to be
                ! up to 50000 characters long. This allows popsfiles to be
                ! written in jobs with up to around 2940 MPI processes (the
                ! default value caused crashes when using over 2000 processes).
                open(iunit, file=popsfile, status='replace', recl=50000)
                call write_popsfile_header (iunit, write_count_sum, &
                                            WalkersonNodes)
                close(iunit)
            end if
        end if

        ! And stop timing
        call halt_timer(write_timer)

        ! Reset some globals
        AllSumNoatHF = 0
        AllSumENum = 0
        AllTotWalkers = 0

    end subroutine WriteToPopsfileParOneArr

    subroutine write_popsfile_header (iunit, num_walkers, WalkersonNodes)

        ! Write the popsfile header into the file specified.

        integer, intent(in) :: iunit
        integer(int64), intent(in) :: num_walkers
        integer :: pops_niftot, pops_nifflag, i 
        integer(int64), intent(in) :: WalkersonNodes(:)

        ! If the popsfile uses flags, but we have combined the
        ! representation of flags in memory, then the representation in
        ! the popsfile is one unit longer than in memory.
        pops_niftot = NIfTot
        pops_nifflag = NIfFlag
        if (tUseFlags .and. NIfFlag == 0) then
            pops_niftot = pops_niftot + 1
            pops_nifflag = 1
        end if

        write(iunit, '(a)') '# POPSFILE VERSION 4'
        write(iunit, '("&POPSHEAD Pop64Bit=",l1)') build_64bit
        write(iunit, '(a,l1,a,l1,a,i2,a,i3,a)') &
            'PopHPHF=', tHPHF, ',PopLz=', tFixLz, ',PopLensign=', &
            lenof_sign, ',PopNEl=', NEl, ','
        write(iunit, '(a,i15,a,f18.12,a)') &
            'PopTotwalk=', num_walkers, ',PopSft=', DiagSft(1), ','
        write(iunit, *) 'PopSumNoatHF=', AllSumNoatHF(1), ','
        write(iunit, *) 'PopSumENum=', AllSumENum(1), ','
        write(iunit, '(a,i16,a,i2,a,i2,a,i2,a)') &
            'PopCyc=', Iter+PreviousCycles, ',PopNIfD=', NIfD, &
            ',PopNIfY=', NIfY, ',PopNIfSgn=' ,NIfSgn, ','
        write(iunit, '(a,i2,a,i2,a,f18.12,a)') &
            'PopNIfFlag=', pops_nifflag, ',PopNIfTot=', pops_niftot, &
            ',PopTau=', Tau, ','
        write(iunit, '(a,i16)') 'PopiBlockingIter=', iBlockingIter(1)
        write(iunit, '(a,f18.12,a,f18.12)') 'PopPSingles=', pSingles, &
            ',PopPParallel=', pParallel

        ! What is the current total imaginary time? Should continue from where
        ! we left off, so that plots work correctly.
        write(iunit, '(a,f18.12)') 'PopTotImagTime=', TotImagTime

        ! Write out accumulated data used for tau searching, to ensure there
        ! are no blips in particle growth, tau, etc.
        write(iunit, '(5(a,f18.12))') 'PopGammaSing=', gamma_sing, &
                                      ',PopGammaDoub=', gamma_doub, &
                                      ',PopGammaOpp=', gamma_opp, &
                                      ',PopGammaPar=', gamma_par, &
                                      ',PopMaxDeathCpt=', max_death_cpt

        if (.not. tSplitPops) then
            ! Write out the number of particles on each processor.
            write(iunit, '(a,i6)') 'PopNNodes=', nNodes
            write(iunit, '(a)', advance='no') "PopWalkersOnNodes="
            do i = lbound(WalkersonNodes, 1), ubound(WalkersonNodes, 1)
                write(iunit, '(i16,",")', advance='no') WalkersOnNodes(i)
            end do
            write(iunit, *)
        end if

        ! Store the random hash in the header to allow later processing
        write(iunit, '(a)', advance='no') "PopRandomHash= "
        do i = 1, nbasis
            write(iunit,'(i12,",")', advance='no') RandomHash(i)
        end do
        write(iunit, *)

        write(iunit, '(a5)') '&END'

    end subroutine


    function write_pops_det (iunit, iunit_2, det, j) result(bWritten)

        ! Output a particle to a popsfile in format acceptable for popsfile v4

        integer, intent(in) :: iunit, iunit_2
        integer(n_int), intent(in) :: det(0:NIfTot)
        real(dp) :: real_sgn(lenof_sign), detenergy
        integer :: flg, j, k, ex_level, nopen, nI(nel) 
        logical :: bWritten

        bWritten = .false.

        call extract_sign(det, real_sgn)

        if (tBinPops) then
            ! We don't want to bother outputting empty particles, or those
            ! with a weight which is lower than specified as the cutoff
            if (sum(abs(real_sgn)) > binarypops_min_weight) then
                if (mod(j, iPopsPartEvery) == 0) then 
                    bWritten = .true.
                end if
            end if
        else
            ! If not using a binary popsfile, the header has already been written,
            ! and the total number of determinants was written out. Therefore, all
            ! determinants must be printed out, else there will be an error reading
            ! them in again.
            bWritten = .true.
        end if

        if (bWritten) then

            ! Write output in the desired format. If __INT64, we are 
            ! including the flag information with the signs in storage in
            ! memory --> need to extract these before outputting them.
            flg = extract_flags(det)
            if (tBinPops) then
                ! All write statements MUST be on the same line, or we end
                ! up with multiple records.
                ! TODO: For POPSFILE V5 --> stream output.
                if (tUseFlags) then
                    write(iunit) det(0:NIfD), real_sgn, int(flg, n_int)
                else
                    write(iunit) det(0:NIfD), real_sgn
                end if
            else
                do k = 0, NIfDBO
                    write(iunit, '(i24)', advance='no') det(k)
                end do
                do k = 1, lenof_sign
                    write(iunit, '(f30.8)', advance='no') real_sgn(k)
                end do
                if (tUseFlags) write(iunit, '(i24)', advance='no') flg
                write(iunit, *)
            end if

            if (tPrintInitiators .and. &
                abs(real_sgn(1)) > InitiatorWalkNo) then
                ! Testing using the sign now, because after annihilation
                ! the current flag will not necessarily be correct.
                ex_level = FindBitExcitLevel(ilutRef, det, nel)
                nopen = count_open_orbs(det)
                call decode_bit_det(nI, det)
                if(tHPHF)then
                    detenergy = hphf_diag_helement(nI, det)
                else
                    detenergy = get_helement(nI, nI, 0)
                endif
                write(iunit_2, '(f20.10,a20)', advance='no') &
                    abs(real_sgn(1)), ''
                call writebitdet (iunit_2, det, .false.)
                write(iunit_2,'(i30,i30,f20.10)') ex_level, nopen, detenergy

            end if

        end if

    end function write_pops_det

!This routine reads in particle configurations from a POPSFILE.
    SUBROUTINE ReadFromPopsfilePar()
        LOGICAL :: exists,tBinRead
        Real(dp) :: NodeSumNoatHF(nProcessors)
        INTEGER :: WalkerstoReceive(nProcessors)
        integer(int64) :: AvWalkers
        integer(int64) :: TempTotParts(lenof_sign),TempCurrWalkers
        INTEGER :: TempInitWalkers,error,i,j,l,total,ierr,MemoryAlloc,Tag,Proc,CurrWalkers,ii
        INTEGER , DIMENSION(lenof_sign) :: TempSign
        real(dp) :: RealTempSign(lenof_sign)
        integer(int64) :: iLutTemp64(0:nBasis/64+1)
        INTEGER :: iLutTemp32(0:nBasis/32+1)
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        INTEGER :: AvSumNoatHF,IntegerPart,TempnI(NEl),ExcitLevel
        INTEGER :: NIfWriteOut,pos,orb,PopsVersion, iunit, iunit_3
        real(dp) :: r, FracPart, Gap, DiagSftTemp, tmp_dp
        HElement_t :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        character(255) :: popsfile,FirstLine
        character(len=24) :: junk,junk2,junk3,junk4
        LOGICAL :: tPop64BitDets,tPopHPHF,tPopLz,tPopInitiator
        integer(n_int) :: ilut_largest(0:NIfTot)
        real(dp) :: sign_largest

        WRITE(6,*) "THIS IS THE POPSFILE ROUTINE WE'RE USING"

        if (lenof_sign /= 1) &
            call Stop_All(this_routine, "Popsfile V.2 does not work with &
                                        &complex walkers")

        PreviousCycles=0    !Zero previous cycles
        SumENum=0.0_dp
        TotParts=0
        SumNoatHF=0
        DiagSft=0.0_dp
        Tag=124             !Set Tag

        call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
        iunit = get_free_unit()
        INQUIRE(FILE=popsfile,EXIST=exists)
        IF(exists) THEN
            OPEN(iunit,FILE=popsfile,Status='old')
            tBinRead=.false.
        ELSE
            tBinRead=.true.
            call stop_all("READ PAR", "BINARY")
            call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            INQUIRE(FILE=popsfile,EXIST=exists)
            IF(.not.exists) THEN
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"No POPSFILEs of any kind found.")
                ELSE
                    CALL Stop_All(this_routine,"POPSFILEBIN(.x) found, but POPSFILEHEAD(.x) also needed for header information")
                ENDIF
            ELSE
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"POPSFILEHEAD(.x) found, but no POPSFILEBIN(.x) " &
                    & //"for particle information - this is also needed")
                ELSE
                    call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                    OPEN(iunit,FILE=popsfile,Status='old')
                ENDIF
            ENDIF
        ENDIF

        READ(iunit,'(a255)') FirstLine

        IF(INDEX(FirstLine,'VERSION').eq.0) THEN
!No version number to be found
            PopsVersion=1
            REWIND(iunit)
        ELSE
            !Found version - which number is it?
            REWIND(iunit)
            READ(iunit,*) FirstLine,FirstLine,FirstLine,PopsVersion
        ENDIF
        WRITE(6,"(A,I5,A)") "Version",PopsVersion," POPSFILE detected"

!Read in initial data on processors which have a popsfile
        IF(PopsVersion.eq.2) THEN
            READ(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') junk,tPop64BitDets,junk2,tPopHPHF,junk3,tPopLz,junk4,tPopInitiator
        ELSE
            WRITE(6,'(A)') "Reading in from depreciated POPSFILE - assuming that parameters " &
            & //"are the same as when POPSFILE was written"
        ENDIF
        READ(iunit,*) tmp_dp
        AllTotWalkers = int(tmp_dp,int64)
        READ(iunit,*) DiagSftTemp
        READ(iunit,*) AllSumNoatHF
        READ(iunit,*) AllSumENum
        READ(iunit,*) PreviousCycles

        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,I4,A)") "Although ",AllTotWalkers, &
                 " configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            ENDIF
        ENDIF

        IF(.not.tWalkContGrow) THEN
!If we want the walker number to continue growing, then take the diagonal 
! shift from the input, rather than the POPSFILE.
            DiagSft=DiagSftTemp
        ENDIF

        IF(DiagSftTemp.eq.0.0_dp) THEN
            tWalkContGrow=.true.
            DiagSft=DiagSftTemp
        ENDIF

        IF(tBinRead) THEN
!Test for the end of the file.
!If this is not the end of the file, there is one more keyword that tells us 
! the calculation had not entered variable shift mode yet.
!Want to put this test at the end of the non-binary file too.
            CLOSE(iunit)
            call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            OPEN(iunit,FILE=popsfile,Status='old',form='unformatted')
        ENDIF

        IF(iProcIndex.eq.Root) THEN

            WRITE(6,*) "Number of cycles in previous simulation: ",PreviousCycles
            IF(NEquilSteps.gt.0) THEN
                WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
                NEquilSteps=0
            ENDIF
            IF(TZeroProjE) THEN
!Reset energy estimator
                WRITE(6,*) "Resetting projected energy counters to zero..."
                AllSumENum=0.0_dp
                AllSumNoatHF = 0
            ENDIF

!Need to calculate the number of walkers each node will receive...
            AvWalkers=NINT(AllTotWalkers/real(nProcessors,dp),int64)

!Divide up the walkers to receive for each node
            do i=1,nProcessors-1
                WalkerstoReceive(i)=int(AvWalkers,sizeof_int)
            enddo
!The last processor takes the 'remainder'
            WalkerstoReceive(nProcessors)=int(AllTotWalkers-(AvWalkers*(nProcessors-1)),sizeof_int)

!Quick check to ensure we have all walkers accounted for
            total=0
            do i=1,nProcessors
                total=total+WalkerstoReceive(i)
            enddo
            if (total /= AllTotWalkers) then
                CALL Stop_All("ReadFromPopsfilePar","All Walkers not accounted for when reading in from POPSFILE")
            endif

!InitWalkers needs to be reset for the culling criteria
            IF(.not.tWalkContGrow) THEN
!Now, let the total space allocated for storing walkers which have been read in to be equal to the initwalkers from the input file.
!                InitWalkers=AvWalkers
            ELSE
                TSinglePartPhase=.true.
            ENDIF
            SumENum=AllSumENum/REAL(nProcessors,dp)     !Divide up the SumENum over all processors
            AvSumNoatHF = int(AllSumNoatHF(1)/nProcessors,sizeof_int) !This is the average Sumnoathf
            do i=1,nProcessors-1
                NodeSumNoatHF(i)=real(INT(AvSumNoatHF,int64),dp)
            enddo
            NodeSumNoatHF(nProcessors)=real(AllSumNoatHF(1)-INT((AvSumNoatHF*(nProcessors-1)),int64),dp)

            ProjectionE=AllSumENum/real(AllSumNoatHF(1),dp)

!Reset the global variables
            AllSumENum=0.0_dp
            AllSumNoatHF = 0

        ENDIF

        CALL MPIBarrier(error)  !Sync

!Now we need to scatter the WalkerstoReceive to each node, and allocate the desired memory to each node...
!Broadcast info which needs to go to all processors
        CALL MPIBCast(DiagSft)
        CALL MPIBCast(SumENum)
        CALL MPIBCast(InitWalkers)
        CALL MPIBCast(NEquilSteps)
        CALL MPIBCast(NShiftEquilSteps)
        CALL MPIBCast(TSinglePartPhase)
!        CALL MPI_BCast(tChangenProcessors,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
!Scatter the number of walkers each node will receive to TempInitWalkers, and 
! the SumNoatHF for each node which is distributed approximatly equally
        CALL MPIScatter(WalkerstoReceive,TempInitWalkers,error)
        CALL MPIScatter(NodeSumNoatHF,SumNoatHF(1),error)

        IF(MemoryFacPart.le.1.0_dp) THEN
            WRITE(6,*) 'MemoryFacPart must be larger than 1.0 when reading in a POPSFILE - increasing it to 1.50.'
            MemoryFacPart=1.50
        ENDIF

!Now we want to allocate memory on all nodes.
!InitWalkers here is simply the average number of walkers per node, not actual
        MaxWalkersPart=NINT(MemoryFacPart*(NINT(InitWalkers*ScaleWalkers)))
        MaxSpawned=NINT(MemoryFacSpawn*(NINT(InitWalkers*ScaleWalkers*inum_runs)))

        Gap=REAL(MaxSpawned,dp)/REAL(nProcessors,dp)
        do i=0,nProcessors-1
            InitialSpawnedSlots(i)=NINT(Gap*i)+1
        enddo
!ValidSpawndList now holds the next free position in the newly-spawned list, 
! but for each processor.
        ValidSpawnedList(:)=InitialSpawnedSlots(:)

        CALL MPIBarrier(error)  !Sync

!Allocate memory to hold walkers at least temporarily
        ALLOCATE(WalkVecDets(0:NIfTot,MaxWalkersPart),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All(this_routine,'Problem allocating WalkVecDets array.')
        CALL LogMemAlloc('WalkVecDets',MaxWalkersPart*(NIfTot+1),size_n_int,this_routine,WalkVecDetsTag,ierr)
        WalkVecDets(:,:)=0
        MemoryAlloc=(NIfTot+1)*MaxWalkersPart*size_n_int    !Memory Allocated in bytes

!Just allocating this here, so that the SpawnParts arrays can be used for sorting the determinants when using direct annihilation.
        WRITE(6,"(A,I12,A)") " Spawning vectors allowing for a total of ",MaxSpawned, &
         &" particles to be spawned in any one iteration."
        ALLOCATE(SpawnVec(0:NIfTot,MaxSpawned),stat=ierr)
        CALL LogMemAlloc('SpawnVec',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVecTag,ierr)
        SpawnVec(:,:)=0
        ALLOCATE(SpawnVec2(0:NIfTot,MaxSpawned),stat=ierr)
        CALL LogMemAlloc('SpawnVec2',MaxSpawned*(NIfTot+1),size_n_int,this_routine,SpawnVec2Tag,ierr)
        SpawnVec2(:,:)=0
        MemoryAlloc=MemoryAlloc+MaxSpawned*(1+NIfTot)*2*size_n_int

!Point at correct spawning arrays
        SpawnedParts=>SpawnVec
        SpawnedParts2=>SpawnVec2
!Allocate pointer to the correct walker array...
        CurrentDets=>WalkVecDets

        ! Allocate storage for persistent data to be stored alongside
        ! the current determinant list (particularly diagonal matrix
        ! elements, i.e. CurrentH; now global_determinant_data).
        call init_global_det_data()

        if(tRDMonFly.and.(.not.tExplicitAllRDM)) then
!Allocate memory to hold walkers spawned from one determinant at a time.
!Walkers are temporarily stored here, so we can check if we're spawning onto the same Dj multiple times.
            ALLOCATE(TempSpawnedParts(0:NIfDBO,20000),stat=ierr)
            CALL LogMemAlloc('TempSpawnedParts',20000*(NIfDBO+1),size_n_int,this_routine,TempSpawnedPartsTag,ierr)
            TempSpawnedParts(0:NIfDBO,1:20000)=0
            MemoryAlloc=MemoryAlloc + (NIfDBO+1)*20000*size_n_int    !Memory Allocated in bytes
            WRITE(6,"(A)") " Allocating temporary array for walkers spawned from a particular Di."
            WRITE(6,"(A,F14.6,A)") " This requires ", REAL(((NIfDBO+1)*20000*size_n_int),dp)/1048576.0_dp," Mb/Processor"
        endif

! The hashing will be different in the new calculation from the one where the
!  POPSFILE was produced, this means we must recalculate the processor each 
! determinant wants to go to.
! This is done by reading in all walkers to the root and then distributing 
! them in the same way as the spawning steps are done - by finding the
!  determinant and sending it there.
        IF((PopsVersion.ne.1).and.tHPHF.and.(.not.tPopHPHF)) THEN
            CALL Stop_All(this_routine,"HPHF on, but HPHF was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.tFixLz.and.(.not.tPopLz)) THEN
            CALL Stop_All(this_routine,"Lz on, but Lz was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tHPHF).and.tPopHPHF) THEN
            CALL Stop_All(this_routine,"HPHF off, but HPHF was used for creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tFixLz).and.tPopLz) THEN
            CALL Stop_All(this_routine,"Lz off, but Lz was used for creation of the POPSFILE")
        ENDIF
        ! TODO: Add tests for CSFs here.
        IF(PopsVersion.eq.1) THEN
            tPop64BitDets=.false.
            NIfWriteOut=nBasis/32
            IF(tCSF) NIfWriteOut=NIfWriteOut+1
        ELSE
            IF(.not.tPop64BitDets) THEN
                NIfWriteOut=nBasis/32
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ELSE
                NIfWriteOut=nBasis/64
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ENDIF
        ENDIF

        CurrWalkers=0
        sign_largest = 0
        ilut_largest = 0
        do i=1,int(AllTotWalkers,sizeof_int)
            iLutTemp(:)=0
            IF(PopsVersion.ne.1) THEN
                IF(tBinRead) THEN
                    IF(tPop64BitDets) THEN
                        READ(iunit) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ELSE
                    IF(tPop64BitDets) THEN
                        READ(iunit,*) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ENDIF
            ELSE
                !POPSFILE v. 1 only printed out 32 bit determinant strings.
                IF(tBinRead) THEN
                    READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                ELSE
                    READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                ENDIF
            ENDIF
            do j=1,lenof_sign
                RealTempSign(j) = transfer(TempSign(j), RealTempSign(j))
            enddo

#ifdef __INT64
            if (.not.tPop64BitDets) then
                ! If we are using 64 bit integers, but have read in 32 bit 
                ! integers, then we need to convert them.
                do ii=0,nBasis/32
                    do j=0,31
                        if(btest(iLutTemp32(ii),j)) then
                            orb=(ii*32)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp64(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp64(0:NIfDBO)
            endif

#else
            ! If we are using 32 bit integers, but have read in 64 bit 
            ! integers, then we need to convert them.
            if (tPop64BitDets) then
                do ii=0,nBasis/64
                    do j=0,63
                        if(btest(iLutTemp64(ii),j)) then
                            orb=(ii*64)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp32(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp32(0:NIfDBO)
            endif

#endif
            call decode_bit_det (TempnI, iLutTemp)
            Proc = DetermineDetNode(nel,TempnI,0)
            IF((Proc.eq.iNodeIndex).and.(abs(RealTempSign(1)).ge.iWeightPopRead)) THEN
                CurrWalkers=CurrWalkers+1
                !Do not need to send a flag here...
                call encode_bit_rep(CurrentDets(:,CurrWalkers),iLutTemp(0:NIfDBO),RealTempSign,0) 
                !TODO: Add flag for complex walkers to read in both
                
!>>>"                CurrentH(1:1+2*lenof_sign,CurrWalkers)=CurrentHEntry(1:1+2*lenof_sign)
            ENDIF

            ! Keep track of what the most highly weighted determinant is
            if (abs(RealTempSign(1)) > sign_largest) then
                sign_largest = abs(RealTempSign(1))
                ilut_largest = iLutTemp
            endif
        enddo
        CLOSE(iunit)
        TempCurrWalkers=int(CurrWalkers,int64)

        ! Sort the lists so that they are in order if we change the number
        ! of processors.
        call sort (currentdets(:,1:CurrWalkers), &
                   global_determinant_data(:,1:CurrWalkers))

        ! Check that the bit-det comparisons agree that it is in order.
        do i=2,currwalkers
            if(DetBitLT(CurrentDets(:,i),CurrentDets(:,i-1),NIfDBO) == 1) then
                print*, 'Walkers: ', i-1, i
                print*, 'bit reps: '
                print*, currentdets(:, i-1)
                print*, currentdets(:, i)
                call stop_all (this_routine, 'Main list out of order')
            endif
        enddo

        CALL MPIBarrier(error)  !Sync
        CALL MPIAllReduce(TempCurrWalkers,MPI_SUM,AllTotWalkers)

        IF(iProcIndex.eq.root) WRITE(6,'(I10,A)') INT(AllTotWalkers,int64)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers by a factor of: ",ScaleWalkers

! CurrWalkers is the number of determinants on a particular node, AllTotWalkers is the total over all nodes.
            IntegerPart=INT(ScaleWalkers)
            FracPart=ScaleWalkers-REAL(IntegerPart,dp)

            do l=1,CurrWalkers
                call extract_sign(CurrentDets(:,l),RealTempSign)
                RealTempSign=RealTempSign*IntegerPart
                r = genrand_real2_dSFMT() 
                IF(r.lt.FracPart) THEN
!Stochastically create another particle
                    IF(RealTempSign(1).lt.0) THEN
                        RealTempSign(1)=RealTempSign(1)-1
                    ELSE
                        RealTempSign(1)=RealTempSign(1)+1
                    ENDIF
                ENDIF
                call encode_sign(CurrentDets(:,l),RealTempSign)
            enddo

            InitWalkers=NINT(InitWalkers*ScaleWalkers)  !New (average) number of initial particles for culling criteria
!Other parameters don't change (I think) because the number of determinants isn't changing.                
            TotWalkers=CurrWalkers
            TotWalkersOld=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                AllTotWalkersOld=AllTotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ELSE
!We are not scaling the number of walkers...

            TotWalkers=CurrWalkers
            TotWalkersOld=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                AllTotWalkersOld=AllTotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ENDIF

        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft(1)
        WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,dp)/1048576.0_dp," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        WRITE(6,*) "Excitgens will be regenerated when they are needed..."
        CALL neci_flush(6)

        ! If we are changing the reference determinant to the largest
        ! weighted one in the file, do it here
        if (tReadPopsChangeRef .or. tReadPopsRestart) then
            if (.not. DetBitEq(ilut_largest, iLutRef, NIfDBO)) then

                ! Set new reference
                iLutRef = ilut_largest
                call decode_bit_det (ProjEDet, iLutRef)
                tNoBrillouin = .true.

                ! Recalculate the reference E
                if (tHPHF) then
                    HElemTemp = hphf_diag_helement (ProjEDet, iLutRef)
                else
                    HElemTemp = get_helement (ProjEDet, ProjEDet, 0)
                endif
                Hii = real(HElemTemp, dp)

                ! Output info on root node.
                if (iProcIndex == root) then
                    write(6, '(a)', advance='no') &
                        "Changing projected energy reference determinant to: "
                    call write_det (6, ProjEDet, .true.)
                    write (6, '(a)') &
                        "Ensuring that Brillouin's theorem is no longer used."
                    write (6, '(a,g25.15)') &
                        "Reference energy now set to: ", Hii
                endif
            endif
        endif

!Now find out the data needed for the particles which have been read in...
        TotParts=0
        do j=1,int(TotWalkers,sizeof_int)
            call decode_bit_det (TempnI, currentDets(:,j))
            Excitlevel = FindBitExcitLevel(iLutHF, CurrentDets(:,j), 2)
            IF(Excitlevel.eq.0) THEN
                call set_det_diagH(j, 0.0_dp)
            ELSE
                if (tHPHF) then
                    HElemTemp = hphf_diag_helement (TempnI, &
                                                    CurrentDets(:,j))
                else
                    HElemTemp = get_helement (TempnI, TempnI, 0)
                endif
                call set_det_diagH(j, real(HElemTemp, dp) - Hii)
            ENDIF
            call extract_sign(CurrentDets(:,j),RealTempSign)
            TotParts(1)=TotParts(1)+abs(RealTempSign(1))
            TotParts(inum_runs)=TotParts(1)
        enddo

        CALL MPIBarrier(error)  !Sync
        CALL MPIReduce(TotParts,MPI_SUM,AllTotParts)

        IF(iProcIndex.eq.root) then
            AllTotPartsOld=AllTotParts
            iter_data_fciqmc%tot_parts_old = AllTotParts
        endif
        write(6,'(A,i20)') ' The total number of particles read from the POPSFILE is: ',AllTotParts(1)

        if (tReadPopsRestart) then
            tPopsAlreadyRead = .true.
            call ChangeRefDet (ProjEDet)
            tPopsAlreadyRead = .false.
        endif

    END SUBROUTINE ReadFromPopsfilePar

!This routine reads in particle configurations from a POPSFILE.
! It's a bastardisation of ReadFromPopsfilePar but only reads a popsfile into an array and leaves the rest until later.
! nDets comes in with the max number of amplitudes possible in Dets.
! Dets is then filled with all the amplitudes in the POPSFILE and nDets is returned with that number.

! This will very likely still be tweaking the inner heart-strings of FciMCPar.  Caveatis stulti!  AJWT
! GHB says he will incorporate this functionality into a rewrite of ReadFromPopsfilePar. 19/8/2010
    SUBROUTINE ReadFromPopsfileOnly(Dets,nDets)
        use CalcData , only : MemoryFacPart,MemoryFacSpawn,iWeightPopRead
        use LoggingData, only: tZeroProjE
        use constants, only: size_n_int,bits_n_int
        integer(int64),intent(inout) :: nDets !The number of occupied entries in Dets
        integer(kind=n_int),intent(out) :: Dets(0:nIfTot,1:nDets)
        LOGICAL :: exists,tBinRead
        integer(int64) :: TempCurrWalkers
        REAL(dp) :: TempTotParts(lenof_sign)
        INTEGER :: TempInitWalkers,error,i,j,l,total,ierr,MemoryAlloc,Tag,Proc,CurrWalkers,ii
        INTEGER , DIMENSION(lenof_sign) :: TempSign
        REAL(dp) , DIMENSION(lenof_sign) :: RealTempSign
        integer(int64) :: iLutTemp64(0:nBasis/64+1)
        INTEGER :: iLutTemp32(0:nBasis/32+1)
        INTEGER(KIND=n_int) :: iLutTemp(0:NIfTot)
        INTEGER :: AvSumNoatHF,IntegerPart,TempnI(NEl),ExcitLevel
        INTEGER :: NIfWriteOut,pos,orb,PopsVersion, iunit
        real(dp) :: r,FracPart,Gap,DiagSftTemp
        HElement_t :: HElemTemp
        CHARACTER(len=*), PARAMETER :: this_routine='ReadFromPopsfilePar'
        character(255) :: popsfile,FirstLine
        character(len=24) :: junk,junk2,junk3,junk4
        LOGICAL :: tPop64BitDets,tPopHPHF,tPopLz,tPopInitiator
        integer(n_int) :: ilut_largest(0:NIfTot)
        real(dp) :: sign_largest

        MemoryAlloc = 0

        IF(lenof_sign.ne.1) CALL Stop_All("ReadFromPopsfilePar","Popsfile does not work with complex walkers")

        PreviousCycles=0    !Zero previous cycles
        SumENum=0.0_dp
        TotParts=0
        SumNoatHF=0
        DiagSft=0.0_dp
        Tag=124             !Set Tag

        call get_unique_filename('POPSFILE',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
        iunit = get_free_unit()
        INQUIRE(FILE=popsfile,EXIST=exists)
        IF(exists) THEN
            OPEN(iunit,FILE=popsfile,Status='old')
            tBinRead=.false.
        ELSE
            tBinRead=.true.
            call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            INQUIRE(FILE=popsfile,EXIST=exists)
            IF(.not.exists) THEN
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"No POPSFILEs of any kind found.")
                ELSE
                    CALL Stop_All(this_routine,"POPSFILEBIN(.x) found, but POPSFILEHEAD(.x) also needed for header information")
                ENDIF
            ELSE
                call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                INQUIRE(FILE=popsfile,EXIST=exists)
                IF(.not.exists) THEN
                    CALL Stop_All(this_routine,"POPSFILEHEAD(.x) found, but no POPSFILEBIN(.x) for " &
                    & //"particle information - this is also needed")
                ELSE
                    call get_unique_filename('POPSFILEHEAD',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
                    OPEN(iunit,FILE=popsfile,Status='old')
                ENDIF
            ENDIF
        ENDIF

        READ(iunit,'(a255)') FirstLine

        IF(INDEX(FirstLine,'VERSION').eq.0) THEN
!No version number to be found
            PopsVersion=1
            REWIND(iunit)
        ELSE
            !Found version - which number is it?
            REWIND(iunit)
            READ(iunit,*) FirstLine,FirstLine,FirstLine,PopsVersion
        ENDIF
        WRITE(6,"(A,I5,A)") "Version",PopsVersion," POPSFILE detected"



!Read in initial data on processors which have a popsfile
        IF(PopsVersion.eq.2) THEN
            READ(iunit,'(A12,L5,A8,L5,A8,L5,A12,L5)') junk,tPop64BitDets,junk2,tPopHPHF,junk3,tPopLz,junk4,tPopInitiator
        ELSE
            WRITE(6,'(A)') "Reading in from depreciated POPSFILE - assuming that parameters are same as when POPSFILE was written"
        ENDIF
        READ(iunit,*) AllTotWalkers
        READ(iunit,*) DiagSftTemp
        READ(iunit,*) AllSumNoatHF
        READ(iunit,*) AllSumENum
        READ(iunit,*) PreviousCycles

        IF(iProcIndex.eq.Root) THEN
            IF(iWeightPopRead.ne.0) THEN
                WRITE(6,"(A,I15,A,I4,A)") "Although ",AllTotWalkers,    &
                " configurations will be read in, only determinants with a weight of over ",iWeightPopRead," will be stored."
            ENDIF
        ENDIF

        IF(.not.tWalkContGrow) THEN
!If we want the walker number to continue growing, then take the diagonal shift from the input, rather than the POPSFILE.
            DiagSft=DiagSftTemp
        ENDIF

        IF(DiagSftTemp.eq.0.0_dp) THEN
            tWalkContGrow=.true.
            DiagSft=DiagSftTemp
        ENDIF

        IF(tBinRead) THEN
!Test for the end of the file.
!If this is not the end of the file, there is one more keyword that tells us 
! the calculation had not entered variable shift mode yet.
!Want to put this test at the end of the non-binary file too.
            CLOSE(iunit)
            call get_unique_filename('POPSFILEBIN',tIncrementPops,.false.,iPopsFileNoRead,popsfile)
            OPEN(iunit,FILE=popsfile,Status='old',form='unformatted')
        ENDIF

        IF(iProcIndex.eq.Root) THEN

            WRITE(6,*) "Number of cycles in previous simulation: ",PreviousCycles
            IF(NEquilSteps.gt.0) THEN
                WRITE(6,*) "Removing equilibration steps since reading in from POPSFILE."
                NEquilSteps=0
            ENDIF
            IF(TZeroProjE) THEN
!Reset energy estimator
                WRITE(6,*) "Resetting projected energy counters to zero..."
                AllSumENum=0.0_dp
                AllSumNoatHF = 0
            ENDIF

!InitWalkers needs to be reset for the culling criteria
            IF(.not.tWalkContGrow) THEN
!Now, let the total space allocated for storing walkers which have been read in to be equal to the initwalkers from the input file.
!                InitWalkers=AvWalkers
            ELSE
                TSinglePartPhase=.true.
            ENDIF
            SumENum=AllSumENum/REAL(nProcessors,dp)     !Divide up the SumENum over all processors
            AvSumNoatHF = int(AllSumNoatHF(1),sizeof_int)/nProcessors !This is the average Sumnoathf

            ProjectionE=AllSumENum/real(AllSumNoatHF(1),dp)

!Reset the global variables
            AllSumENum=0.0_dp
            AllSumNoatHF = 0

        ENDIF

        CALL MPIBarrier(error)  !Sync

!Now we need to scatter the WalkerstoReceive to each node, and allocate the desired memory to each node...
!Broadcast info which needs to go to all processors
        CALL MPIBCast(DiagSft)
        CALL MPIBCast(SumENum)
        CALL MPIBCast(InitWalkers)
        CALL MPIBCast(NEquilSteps)
        CALL MPIBCast(NShiftEquilSteps)
        CALL MPIBCast(TSinglePartPhase)
!        CALL MPI_BCast(tChangenProcessors,1,MPI_LOGICAL,root,MPI_COMM_WORLD,error)
!Scatter the number of walkers each node will receive to TempInitWalkers, 
! and the SumNoatHF for each node which is distributed approximatly equally

        IF(MemoryFacPart.le.1.0_dp) THEN
            WRITE(6,*) 'MemoryFacPart must be larger than 1.0 when reading in a POPSFILE - increasing it to 1.50.'
            MemoryFacPart=1.50
        ENDIF

        CALL MPIBarrier(error)  !Sync

        if(AllTotWalkers>nDets) CALL Stop_All(this_routine,'Not enough memory to read in POPSFILE.')

! The hashing will be different in the new calculation from the one where the
! POPSFILE was produced, this means we must recalculate the processor each 
! determinant wants to go to.                
! This is done by reading in all walkers to the root and then distributing 
! them in the same way as the spawning steps are done - by finding the 
! determinant and sending it there.
        IF((PopsVersion.ne.1).and.tHPHF.and.(.not.tPopHPHF)) THEN
            CALL Stop_All(this_routine,"HPHF on, but HPHF was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.tFixLz.and.(.not.tPopLz)) THEN
            CALL Stop_All(this_routine,"Lz on, but Lz was not used in creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tHPHF).and.tPopHPHF) THEN
            CALL Stop_All(this_routine,"HPHF off, but HPHF was used for creation of the POPSFILE")
        ENDIF
        IF((PopsVersion.ne.1).and.(.not.tFixLz).and.tPopLz) THEN
            CALL Stop_All(this_routine,"Lz off, but Lz was used for creation of the POPSFILE")
        ENDIF
        ! TODO: Add tests for CSFs here.
        IF(PopsVersion.eq.1) THEN
            tPop64BitDets=.false.
            NIfWriteOut=nBasis/32
            IF(tCSF) NIfWriteOut=NIfWriteOut+1
        ELSE
            IF(.not.tPop64BitDets) THEN
                NIfWriteOut=nBasis/32
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ELSE
                NIfWriteOut=nBasis/64
                IF(tCSF) NIfWriteOut=NIfWriteOut+1
            ENDIF
        ENDIF

        CurrWalkers=0
        sign_largest = 0
        ilut_largest = 0
        write(6,*) "Reading in ", AllTotWalkers, " walkers"
        do i=1,int(AllTotWalkers,sizeof_int)
            iLutTemp(:)=0
            IF(PopsVersion.ne.1) THEN
                IF(tBinRead) THEN
                    IF(tPop64BitDets) THEN
                        READ(iunit) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ELSE
                    IF(tPop64BitDets) THEN
                        READ(iunit,*) iLutTemp64(0:NIfWriteOut),TempSign
                    ELSE
                        READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                    ENDIF
                ENDIF
            ELSE
                !POPSFILE v. 1 only printed out 32 bit determinant strings.
                IF(tBinRead) THEN
                    READ(iunit) iLutTemp32(0:NIfWriteOut),TempSign
                ELSE
                    READ(iunit,*) iLutTemp32(0:NIfWriteOut),TempSign
                ENDIF
            ENDIF
            do j=1,lenof_sign
                RealTempSign(j) = transfer(TempSign(j), RealTempSign(j))
            enddo

#ifdef __INT64
            if (.not.tPop64BitDets) then
                ! If we are using 64 bit integers, but have read in 32 bit 
                ! integers, then we need to convert them.
                do ii=0,nBasis/32
                    do j=0,31
                        if(btest(iLutTemp32(ii),j)) then
                            orb=(ii*32)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp64(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp64(0:NIfDBO)
            endif

#else
            ! If we are using 32 bit integers, but have read in 64 bit 
            ! integers, then we need to convert them.
            if (tPop64BitDets) then
                do ii=0,nBasis/64
                    do j=0,63
                        if(btest(iLutTemp64(ii),j)) then
                            orb=(ii*64)+j+1
                            pos=(orb-1)/bits_n_int
                            iLutTemp(pos)=ibset(iLutTemp(pos),mod(orb-1,bits_n_int))
                        endif
                    enddo
                enddo
                iLutTemp(NIfD+1:NIfDBO) = iLutTemp32(ii:NIfWriteOut)
            else
                iLutTemp(0:NIfDBO)=iLutTemp32(0:NIfDBO)
            endif

#endif
            call decode_bit_det (TempnI, iLutTemp)
            Proc=0  !DetermineDetNode(TempnI)   !This wants to return a value between 0 -> nProcessors-1
            IF((Proc.eq.iProcIndex).and.(abs(RealTempSign(1)).ge.iWeightPopRead)) THEN
                CurrWalkers=CurrWalkers+1
                !Do not need to send a flag here...
                call encode_bit_rep(Dets(:,CurrWalkers),iLutTemp(0:NIfDBO),RealTempSign,0)
                !TODO: Add flag for complex walkers to read in both
            ENDIF

            ! Keep track of what the most highly weighted determinant is
            if (abs(RealTempSign(1)) > sign_largest) then
                sign_largest = abs(RealTempSign(1))
                ilut_largest = iLutTemp
            endif
        enddo
        CLOSE(iunit)
        TempCurrWalkers=int(CurrWalkers,int64)

        ! Sort the lists so that they are in order if we change the number
        ! of processors.
        call sort(Dets(:,1:CurrWalkers))

        ! Check that the bit-det comparisons agree that it is in order.
        do i=2,currwalkers
            if(DetBitLT(Dets(:,i),Dets(:,i-1),NIfDBO) == 1) then
                print*, 'Walkers: ', i-1, i
                print*, 'bit reps: '
                print*, dets(:, i-1)
                print*, dets(:, i)
                call stop_all (this_routine, 'Main list out of order')
            endif
        enddo

        CALL MPIBarrier(error)  !Sync
        CALL MPIAllReduce(TempCurrWalkers,MPI_SUM,AllTotWalkers)

        IF(iProcIndex.eq.root) WRITE(6,'(I10,A)') INT(AllTotWalkers,int64)," configurations read in from POPSFILE and distributed."

        IF(ScaleWalkers.ne.1) THEN

            WRITE(6,*) "Rescaling walkers by a factor of: ",ScaleWalkers

! CurrWalkers is the number of determinants on a particular node, AllTotWalkers is the total over all nodes.
            IntegerPart=INT(ScaleWalkers)
            FracPart=ScaleWalkers-REAL(IntegerPart,dp)

            do l=1,CurrWalkers
                call extract_sign(Dets(:,l),RealTempSign)
                RealTempSign=TempSign*IntegerPart
                r = genrand_real2_dSFMT() 
                IF(r.lt.FracPart) THEN
!Stochastically create another particle
                    IF(RealTempSign(1).lt.0) THEN
                        RealTempSign(1)=RealTempSign(1)-1
                    ELSE
                        RealTempSign(1)=RealTempSign(1)+1
                    ENDIF
                ENDIF
                call encode_sign(Dets(:,l),RealTempSign)
            enddo

            InitWalkers=NINT(InitWalkers*ScaleWalkers)  !New (average) number of initial particles for culling criteria
!Other parameters don't change (I think) because the number of determinants isn't changing.                
            nDets=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ELSE
!We are not scaling the number of walkers...

            nDets=CurrWalkers
            IF(iProcIndex.eq.root) THEN
!                AllTotWalkers=TotWalkers
                WRITE(6,'(A,I10)') " Number of initial walkers on this processor is now: ",INT(TotWalkers,int64)
            ENDIF

        ENDIF

        WRITE(6,*) "Initial Diagonal Shift (ECorr guess) is now: ",DiagSft
        WRITE(6,"(A,F14.6,A)") " Initial memory (without excitgens) consists of : ",REAL(MemoryAlloc,dp)/1048576.0_dp," Mb"
        WRITE(6,*) "Initial memory allocation successful..."
        WRITE(6,*) "Excitgens will be regenerated when they are needed..."
        CALL neci_flush(6)

!Now find out the data needed for the particles which have been read in...
        TotParts=0
        do j=1,int(nDets,sizeof_int)
            call extract_sign(Dets(:,j),RealTempSign)
            TotParts=TotParts+abs(RealTempSign(1))
        enddo

        TempTotParts=TotParts

        CALL MPIBarrier(error)  !Sync
        CALL MPIReduce(TempTotParts,MPI_SUM,AllTotParts)

        IF(iProcIndex.eq.root) AllTotPartsOld=AllTotParts
        write(6,'(A,i20)') ' The total number of particles read from the POPSFILE is: ',AllTotParts(1)

        if (tReadPopsRestart) then
            tPopsAlreadyRead = .true.
            call ChangeRefDet (ProjEDet)
            tPopsAlreadyRead = .false.
        endif

    END SUBROUTINE ReadFromPopsfileOnly

    subroutine add_pops_norm_contrib(ilut)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        real(dp) :: real_sign(lenof_sign)

        call extract_sign(ilut, real_sign)

#ifdef __DOUBLERUN
        pops_norm = pops_norm + real_sign(1)*real_sign(2)
#elif __CMPLX
        pops_norm = pops_norm + real_sign(1)**2 + real_sign(2)**2
#else
        pops_norm = pops_norm + real_sign(1)*real_sign(1)
#endif

    end subroutine add_pops_norm_contrib

    subroutine write_pops_norm()

        use CalcData, only: pops_norm_unit

        integer :: i

        if (iProcIndex /= root) return

        ! When calling for the first time, clear any current POPS_NORM file.
        if (pops_norm_unit == 0) then
            pops_norm_unit = get_free_unit()
            open(pops_norm_unit, file='POPS_NORM', status='replace')
        end if

        write(pops_norm_unit,'(1x,es19.12)') sqrt(pops_norm)

    end subroutine write_pops_norm

END MODULE PopsfileMod
