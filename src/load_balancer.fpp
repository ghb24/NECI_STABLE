#include "macros.h"
#:include "macros.fpph"

module load_balance

    use CalcData, only: tUniqueHFNode, tSemiStochastic, &
                        tCheckHighestPop, OccupiedThresh, &
                        tContTimeFCIMC, t_prone_walkers, &
                        tContTimeFull, tTrialWavefunction, &
                        tPairedReplicas, tau, tSeniorInitiators, &
                        t_activate_decay, tAutoAdaptiveShift, tMoveGlobalDetData
    use global_det_data, only: global_determinant_data, &
                               global_determinant_data_tmp, &
                               set_det_diagH, set_spawn_rate, &
                               set_all_spawn_pops, reset_all_tau_ints, &
                               reset_all_shift_ints, det_diagH, store_decoding, &
                               reset_all_tot_spawns, reset_all_acc_spawns
    use bit_rep_data, only: flag_initiator, NIfDBO, &
                            flag_connected, flag_trial, flag_prone, flag_removed
    use bit_reps, only: set_flag, nullify_ilut_part, &
                        encode_part_sign, nullify_ilut, clr_flag
    use FciMCData, only: HashIndex, FreeSlot, CurrentDets, iter_data_fciqmc, &
                         tFillingStochRDMOnFly, full_determ_vecs, ntrial_excits, &
                         con_space_size, NConEntry, con_send_buf, sFAlpha, sFBeta, &
                         n_prone_dets
    use SystemData, only: tHPHF
    use procedure_pointers, only: scaleFunction
    use searching, only: hash_search_trial, bin_search_trial
    use determinants, only: get_helement, write_det
    use hphf_integrals, only: hphf_diag_helement
    use LoggingData, only: tOutputLoadDistribution, tAccumPopsActive
    use cont_time_rates, only: spawn_rate_full
    use DetBitOps, only: DetBitEq, tAccumEmptyDet
    use sparse_arrays, only: con_ht, trial_ht, trial_hashtable
    use trial_ht_procs, only: buffer_trial_ht_entries, add_trial_ht_entries
    use load_balance_calcnodes
    use Parallel_neci
    use constants
    use util_mod
    use hash

    implicit none

    ! TODO:
    ! - Integrate with POPSFILES. Need to output the mapping for restarts.
    ! - Move global data around
    ! - Note if we move any of the reference sites around!!!!!!

contains

    subroutine init_load_balance()


        ! Initialise the load balancing.
        !
        ! n.b. The initialisation of RandomOrbIndex remains in SetupParameters
        !      to preserve sequencing, which maintains testcode results.

        integer :: oversample_factor, ierr, i
        character(*), parameter :: this_routine = 'init_load_balance'

        !
        ! Initialise the mapping of balancing blocks to nodes. By default
        ! this is just a uniform mapping.
        ASSERT(.not. (tLoadBalanceBlocks .and. tUniqueHFNode))
        if (tLoadBalanceBlocks) then
            oversample_factor = 100
        else
            oversample_factor = 1
        end if

        balance_blocks = nNodes * oversample_factor
        allocate(LoadBalanceMapping(balance_blocks), stat=ierr)
        @:log_alloc(LoadBalanceMapping, lb_tag, ierr)

        ! Generate a uniform mapping(by default)
        do i = 1, balance_blocks
            LoadBalanceMapping(i) = int((i - 1) / oversample_factor)
        end do

    end subroutine

    subroutine pops_init_balance_blocks(pops_blocks)

        use LoggingData, only: tHDF5PopsRead

        integer, intent(in) :: pops_blocks
        integer :: det(nel), block, j
        integer :: mapping_tmp(balance_blocks)
        integer :: mapping_test(balance_blocks)
        integer :: mapping_test_all(balance_blocks)
        real(dp) :: sgn(lenof_sign)
        character(*), parameter :: this_routine = 'pops_init_balance_blocks'

        ! If the particles have been read in using the 'general' routine, then
        ! they will have been distributed according to the already-initialised
        ! blocking structure
        ! --> all we need to do is rebalance things
        ! --> This also applies if we are using the HDF5 popsfile routines,
        !     which distribute the particles at read time
        ASSERT(allocated(LoadBalanceMapping))
        if (pops_blocks == -1 .or. tHDF5PopsRead) then
            call adjust_load_balance(iter_data_fciqmc)
            return
        end if

        write(6,"('Initialising load balancing blocks from data in POPSFILE')")

        ! We can only initialise blocking in this manner if the blocks match
        ! the number of blocks in the popsfile
        if (pops_blocks /= balance_blocks) &
            call stop_all(this_routine, 'Incorrect number of load balancing &
                                        &blocks.')

        ! Determine which blocks are found on which processors. Do this by
        ! setting all the blocks to zero, seeing what is on each processor and
        ! then summing.
        mapping_tmp = 0
        mapping_test = 0
        do j = 1, int(TotWalkers)
            call extract_sign(CurrentDets(:,j), sgn)
            if (IsUnoccDet(sgn) .and. .not. tAccumEmptyDet(CurrentDets(:,j))) cycle

            ! Use ceiling as part-integer particles involve the same
            ! computational cost...
            call decode_bit_det(det, CurrentDets(:,j))
            block = get_det_block(nel, det, 0)

            mapping_tmp(block) = iProcIndex
            mapping_test(block) = 1
        end do
        call MPISumAll(mapping_tmp, LoadBalanceMapping)

        ! Ensure that all of the mappings are set on one, and only one
        ! processor.
        call MPISumAll(mapping_test, mapping_test_all)
        if (.not. all(mapping_test_all == 1 .or. mapping_test_all == 0)) then
            call stop_all(this_routine, "Multi-processor mapping not &
                         &correctly determined")
        end if

        ! And do a load balancing before anything else happens
        call adjust_load_balance(iter_data_fciqmc)

    end subroutine

    subroutine clean_load_balance()

        character(*), parameter :: this_routine = 'clean_load_balance'

        if (allocated(LoadBalanceMapping)) then
            deallocate(LoadBalanceMapping)
            log_dealloc(lb_tag)
        end if

    end subroutine

    subroutine adjust_load_balance(iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer(int64) :: block_parts(balance_blocks)
        integer(int64) :: block_parts_all(balance_blocks)
        integer(int64) :: proc_parts(0:nProcessors-1)
        integer(int64) :: smallest_size
        integer :: j, proc, det(nel), block, TotWalkersTmp
        integer :: min_parts, max_parts, min_proc, max_proc
        integer :: smallest_block,iBlockMoves
        real(dp) :: sgn(lenof_sign), avg_parts
        logical :: unbalanced
        character(*), parameter :: this_routine = 'adjust_load_balance'

        ! TODO: Need to ensure we don't move around the semi-stochastic sites,
        ! or we need to be a bit more clever!!!
        ! TODO: What happens if we move reference sites around?
        ! Actually, we can call this with tSemiStochastic=.true. if we haven't yet
        ! set up the deterministic space.
        if(allocated(full_determ_vecs)) then
            call stop_all(this_routine, &
                'Should not be dynamically load-balancing with fixed deterministic space')
        endif
        if(tFillingStochRDMOnFly) then
            call stop_all(this_routine, &
                'Should not be dynamically load-balancing while sampling RDMs')
        endif

        ! Count the number of particles inside each of the blocks
        block_parts = 0
        HolesInList = 0
        if(tAccumPopsActive)then
            FreeSlot(1:iEndFreeSlot)=0
            iStartFreeSlot=1
            iEndFreeSlot=0
        endif
        do j = 1, int(TotWalkers, sizeof_int)

            call extract_sign(CurrentDets(:,j), sgn)
            if (IsUnoccDet(sgn) .and. .not. tAccumEmptyDet(CurrentDets(:,j))) then
                HolesInList = HolesInList + 1
                if(tAccumPopsActive)then
                    iEndFreeSlot = iEndFreeSlot + 1
                    FreeSlot(iEndFreeSlot) = j
                endif
                cycle
            end if

            ! Use ceiling as part-integer particles involve the same
            ! computational cost...
            call decode_bit_det(det, CurrentDets(:,j))
            block = get_det_block(nel, det, 0)
            block_parts(block) = block_parts(block) + sum(ceiling(abs(sgn)))
        end do

        ! Accumulate the data from all of the processors (on root)
        call MPISum(block_parts, block_parts_all)

        iBlockMoves = 0
        do while (.true.)

            ! n.b. the required data is only available on the root node.
            if (iProcIndex == root) then

                ! How many particles are on each processor?
                proc_parts = 0
                do block = 1, balance_blocks
                    proc = LoadBalanceMapping(block)
                    proc_parts(proc) = proc_parts(proc) + block_parts_all(block)
                end do

                ! Where are the minimal and maximal values found?
                ! n.b. min/maxloc treat all arrays as starting at index 1. sigh
                avg_parts = real(sum(proc_parts), dp) / real(nProcessors, dp)
                min_proc = minloc(proc_parts, dim=1) - 1
                min_parts = int(proc_parts(min_proc))
                max_proc = maxloc(proc_parts, dim=1) - 1
                max_parts = int(proc_parts(max_proc))
                ASSERT(max_proc <= ubound(proc_parts, 1))
                ASSERT(min_proc <= ubound(proc_parts, 1))
                ASSERT(min_proc >= 0)
                ASSERT(max_proc >= 0)

                if (min_proc > nProcessors-1 .or. max_proc > nProcessors-1)&
                    call stop_all(this_routine, 'invalid value')

                ! Create a list of the blocks associated with the most
                ! heavily utilised processor in increasing size order.
                smallest_block = 0
                smallest_size = -1
                do block = 1, balance_blocks
                    if (LoadBalanceMapping(block) == max_proc) then
                        if (block_parts_all(block) > 0 .and. &
                            (block_parts_all(block) < smallest_size .or. &
                                smallest_size == -1)) then
                            smallest_block = block
                            smallest_size = block_parts_all(block)
                        end if
                    end if
                end do

                ! If moving a block of the smallest size between the largest
                ! and the smallest is a helpful thing to do, then move it!
                if (smallest_block /= 0) then
                    if ((abs(min_parts + smallest_size - avg_parts) < abs(min_parts - avg_parts)) .and. &
                        (abs(max_parts - smallest_size - avg_parts) < abs(max_parts - avg_parts))) then
                        unbalanced = .true.
                    else
                        unbalanced = .false.
                    end if
                else
                    unbalanced = .false.
                endif
            end if
            call MPIBcast(unbalanced)

            ! If this is sufficiently balanced, then we make no (further)
            ! changes.
            if (.not. unbalanced) then
                exit
            else
                iBlockMoves = iBlockMoves + 1
            endif

            ! Broadcast the parameters for the change!
            call MPIBCast(min_proc)
            call MPIBCast(max_proc)
            call MPIBcast(smallest_block)

            ! Move the block from where it is to the currently least worked
            ! processor
            call move_block(smallest_block, min_proc)

        end do

        if (iProcIndex == root .and. tOutputLoadDistribution) then
            write(6, '("Load balancing distribution:")')
            write(6, '("node #, particles")')
            do j = 0, nNodes - 1
                write(6,'(i8,i10)') j, proc_parts(j)
            end do
            write(6,*) '--'
        end if

        if(iBlockMoves.gt.0) then
            !Only redo hash table if blocks have been moved around
            !Not only is this an optimization, but it also seems to hide the fact
            !that for gfortran and openmpi 1.6, going through this code at
            !the initialization step of the calculation starting from a single
            !walker seems to cause a race condition later on in the calculation,
            !even though this code has no mpi calls, and is not doing anything with
            !a single determinant.
            !Very confusing...
            TotWalkersTmp = int(TotWalkers, sizeof_int)
            call CalcHashTableStats(TotWalkersTmp, iter_data)
            TotWalkers = int(TotWalkersTmp, int64)
        endif
        !   -- Test if sufficiently uniform
        !   -- If not, pick largest, and smallest, sites
        !   -- Transfer the largest block that will not take either of the
        !      processors past the average number.
        !   -- If there is no such block, use the smallest, ensuring that
        !      we still make improvements.
        !   -- If there is no improving swap, do nothing!

    end subroutine


    subroutine move_block(block, tgt_proc)
        implicit none
        integer, intent(in) :: block, tgt_proc

        integer :: src_proc, ierr, nsend, nelem, j, k, det_block, hash_val, PartInd
        integer :: det(nel), TotWalkersTmp, nconsend, clashes, ntrial, ncon, err
        integer(n_int) :: con_state(0:NConEntry)
        real(dp) :: sgn(lenof_sign)
        real(dp) :: HDiag

        ! A tag is used to identify this send/recv pair over any others
        integer, parameter :: mpi_tag_nsend = 223456
        integer, parameter :: mpi_tag_dets = 223457
        integer, parameter :: mpi_tag_nconsend = 223458
        integer, parameter :: mpi_tag_con = 223459
        integer, parameter :: mpi_tag_ntrialsend = 223460
        integer, parameter :: mpi_tag_trial = 223461
        integer, parameter :: mpi_tag_glob = 223462

        src_proc = LoadBalanceMapping(block)

        ! Provide some feedback to the user.
        if (iProcIndex == root) then
            write(6,'(a,i9,a,i6,a,i6)') 'Moving load balancing block ', &
                     block, ' from processor ', src_proc, ' to ', tgt_proc
        end if

        if (iProcIndex == src_proc) then

            ! Loop over the available walkers, and broadcast them to the
            ! target processor. Use the SpawnedParts array as a buffer.
            nsend = 0
            nconsend = 0
            do j = 1, int(TotWalkers, sizeof_int)

                ! Skip unoccupied sites (non-contiguous)
                call extract_sign(CurrentDets(:,j), sgn)
                if (IsUnoccDet(sgn) .and. .not. tAccumEmptyDet(CurrentDets(:,j))) cycle

                call decode_bit_det(det, CurrentDets(:,j))
                det_block = get_det_block(nel, det, 0)
                if (det_block == block) then
                    nsend = nsend + 1
                    SpawnedParts(0:NIfTot,nsend) = CurrentDets(:,j)
                    if(tMoveGlobalDetData) then
                        global_determinant_data_tmp(:,nsend) = global_determinant_data(:,j)
                    endif

                    ! Remove the det from the main list.
                    call nullify_ilut(CurrentDets(:,j))
                    call RemoveHashDet(HashIndex, det, j)
                end if
            end do

            ! And send the data to the relevant (target) processor
            nelem = nsend * (1 + NIfTot)
            call MPISend(nsend, 1, tgt_proc, mpi_tag_nsend, ierr)
            call MPISend(SpawnedParts(0:NIfTot, 1:nsend), nelem, tgt_proc, &
                         mpi_tag_dets, ierr)
            if(tMoveGlobalDetData) then
                nelem = nsend * SIZE(global_determinant_data_tmp, 1)
                call MPISend(global_determinant_data_tmp(:, 1:nsend), nelem, &
                             tgt_proc, mpi_tag_glob, ierr)
            endif

            ! we only communicate the trial hashtable
            if(tTrialWavefunction .and. tTrialHash) then
               ! now get those connected determinants that need to be
               ! communicated (they might not be in currentdets)
               nconsend = buffer_trial_ht_entries(block, con_ht, con_space_size)
               ! And send the trial wavefunction connection information
               nelem = nconsend * (1 + NConEntry)
               call MPISend(nconsend,1,tgt_proc,mpi_tag_nconsend, ierr)
               if(nelem > 0) then
                  call MPISend(con_send_buf(0:NConEntry,1:nconsend),nelem,tgt_proc, &
                       mpi_tag_con, ierr)
               endif
               ! Do the same with the trial wavefunction itself
               nconsend = buffer_trial_ht_entries(block, trial_ht, trial_space_size)
               nelem = nconsend * (1 + NConEntry)
               call MPISend(nconsend,1,tgt_proc,mpi_tag_ntrialsend, ierr)

               if(nelem > 0) then
                    call MPISend(con_send_buf(0:NConEntry,1:nconsend),nelem,tgt_proc,&
                    mpi_tag_trial, ierr)
                 endif
            end if

            ! We have now created lots of holes in the main list
            HolesInList = HolesInList + nsend

        else if (iProcIndex == tgt_proc) then

            ! Receive walkers!
            call MPIRecv(nsend, 1, src_proc, mpi_tag_nsend, ierr)
            nelem = nsend * (1 + NIfTot)
            call MPIRecv(SpawnedParts(0:NIfTot, 1:nsend), nelem, src_proc, mpi_tag_dets, ierr)
            if(tMoveGlobalDetData) then
                nelem = nsend * SIZE(global_determinant_data_tmp, 1)
                call MPIRecv(global_determinant_data_tmp(:, 1:nsend), nelem, &
                             src_proc, mpi_tag_glob, ierr)
            endif

            do j = 1, nsend
                call decode_bit_det(det, SpawnedParts(:,j))
                call extract_sign(SpawnedParts(:,j), sgn)

                hash_val = FindWalkerHash(det, size(HashIndex))

                ! n.b. Ensure that Totwalkers passed in always has the correct
                !      type even on 32-bit machines
                TotWalkersTmp = int(TotWalkers)

                ! Calculate the diagonal hamiltonian matrix element for the new particle to be merged.
                HDiag = get_diagonal_matel(det, SpawnedParts(:,j))
                call AddNewHashDet(TotWalkersTmp, SpawnedParts(:, j), &
                                   hash_val, det, HDiag, PartInd, err)

                if(tMoveGlobalDetData) then
                    global_determinant_data(:,PartInd) = global_determinant_data_tmp(:,j)
                endif
                TotWalkers = TotWalkersTmp
            end do

            ! We have filled in some of the holes in the list (possibly all)
            ! and possibly extended the list
            HolesInList = max(0, HolesInList - nsend)

            ! Recieve information on the trial + connected determinants
            ! only if trial wavefunction is enabled, of course
            if(tTrialWavefunction .and. tTrialHash) then

               ! first, we get the connected ones
               call MPIRecv(nconsend, 1, src_proc, mpi_tag_nconsend, ierr)
               nelem = nconsend * (1 + NConEntry)
               if(nelem > 0) then
                  ! get the connected states themselves
                  call MPIRecv(con_send_buf(0:NConEntry, 1:nconsend), nelem, src_proc, mpi_tag_con, ierr)
                  ! add the recieved connected dets to the hashtable
                  call add_trial_ht_entries(con_send_buf(:,1:nconsend), nconsend, &
                       con_ht, con_space_size)
               endif
               ! Recieve the information on the trial wave function
               call MPIRecv(nconsend, 1, src_proc, mpi_tag_ntrialsend, ierr)
               nelem = nconsend * (1 + NConEntry)
               if(nelem > 0) then
                  ! get the states
                  call MPIRecv(con_send_buf(0:NConEntry, 1:nconsend), nelem, src_proc, mpi_tag_trial, ierr)
                  ! add them to the hashtable
                  call add_trial_ht_entries(con_send_buf(:,1:nconsend), nconsend, &
                       trial_ht, trial_space_size)
               endif
            endif

        end if

        ! Adjust the load balancing mapping
        LoadBalanceMapping(block) = tgt_proc

        ! And synchronise when everything is done
        call MPIBarrier(ierr)

    end subroutine

    subroutine AddNewHashDet(TotWalkersNew, iLutCurr, DetHash, nJ, HDiag, DetPosition, err)
        ! Add a new determinant to the main list. This involves updating the
        ! list length, copying it across, updating its flag, adding its diagonal
        ! helement (if neccessary). We also need to update the hash table to
        ! point at it correctly.
        integer, intent(inout) :: TotWalkersNew
        integer(n_int), intent(inout) :: iLutCurr(0:NIfTot)
        integer, intent(in) :: DetHash, nJ(nel)
        integer, intent(out) :: DetPosition
        real(dp), intent(in) :: HDiag
        integer, intent(out) :: err
        HElement_t(dp) :: trial_amps(ntrial_excits)
        logical :: tTrial, tCon
        real(dp), dimension(lenof_sign) :: SignCurr
        character(len=*), parameter :: t_r = "AddNewHashDet"

        err = 0
        if (iStartFreeSlot <= iEndFreeSlot) then
            ! We can slot it into a free slot in the main list, rather than increase its length
            DetPosition = FreeSlot(iStartFreeSlot)
            CurrentDets(:, DetPosition) = iLutCurr(:)
            iStartFreeSlot = iStartFreeSlot + 1
        else
            ! We must increase the length of the main list to slot the new walker in
            TotWalkersNew = TotWalkersNew + 1
            DetPosition = TotWalkersNew
            if (TotWalkersNew >= MaxWalkersPart) then
               ! return with an error
               err = 1
               return
            end if
            CurrentDets(:,DetPosition) = iLutCurr(:)

            ! if the list is almost full, activate the walker decay
            if(t_prone_walkers .and. TotWalkersNew > 0.95_dp * real(MaxWalkersPart,dp)) then
               t_activate_decay = .true.
               write(iout,*) "Warning: Starting to randomly kill singly-spawned walkers"
            endif
        end if
        CurrentDets(:,DetPosition) = iLutCurr(:)

        ! For the RDM code we need to set all of the elements of CurrentH to 0,
        ! except the first one, holding the diagonal Hamiltonian element.
        global_determinant_data(:,DetPosition) = 0.0_dp
        call set_det_diagH(DetPosition, real(HDiag,dp) - Hii)

        ! we reset the death timer, so this determinant can linger again if
        ! it died before

        ! we add the determinant to the cache
        call store_decoding(DetPosition, nJ)

        if(tSeniorInitiators) then
            call extract_sign (ilutCurr, SignCurr)
            call set_all_spawn_pops(DetPosition, SignCurr)
            call reset_all_tau_ints(DetPosition)
            call reset_all_shift_ints(DetPosition)
        end if

        if(tAutoAdaptiveShift) then
            call reset_all_tot_spawns(DetPosition)
            call reset_all_acc_spawns(DetPosition)
        end if

        ! If using a trial wavefunction, search to see if this state is in
        ! either the trial or connected space. If so, *_search_trial returns
        ! the corresponding amplitude, which is stored.
        !
        ! n.b. if this routine is called from load balancing whilst loading a popsfile, the
        !      trial wavefunction code will not be enabled yet (despite being enabled).
        !      Therefore, skip this functionality
        if (tTrialWavefunction .and. allocated(current_trial_amps)) then
            ! Search to see if this is a trial or connected state, and
            ! retreive the corresponding amplitude (zero if neither a trial or
            ! connected state).
            if (tTrialHash) then
                call hash_search_trial(CurrentDets(:,DetPosition), nJ, trial_amps, tTrial, tCon)
            else
                call bin_search_trial(CurrentDets(:,DetPosition), trial_amps, tTrial, tCon)
            end if

            ! Set the appropraite flag (if any). Unset flags which aren't
            ! appropriate, just in case.
            if (tTrial) then
                call set_flag(CurrentDets(:,DetPosition), flag_trial, .true.)
                call set_flag(CurrentDets(:,DetPosition), flag_connected, .false.)
             else if (tCon) then
                call set_flag(CurrentDets(:,DetPosition), flag_trial, .false.)
                call set_flag(CurrentDets(:,DetPosition), flag_connected, .true.)
            else
                call set_flag(CurrentDets(:,DetPosition), flag_trial, .false.)
                call set_flag(CurrentDets(:,DetPosition), flag_connected, .false.)
            end if

            ! Set the amplitude (which may be zero).
            current_trial_amps(:,DetPosition) = trial_amps
        end if

        ! If we are storing spawning rates for continuous time propagation, do
        ! it here
        if (tContTimeFCIMC .and. tContTimeFull) then
            call set_spawn_rate(DetPosition, spawn_rate_full(nJ, ilutCurr))
        end if

        !In case we are filling a hole, clear the removed flag
        call set_flag(CurrentDets(:,DetPosition), flag_removed, .false.)

        ! Add the new determinant to the hash table.
        call add_hash_table_entry(HashIndex, DetPosition, DetHash)

    end subroutine AddNewHashDet

    subroutine RemoveHashDet(HashIndex, nJ, partInd)
      implicit none
      type(ll_node), pointer, intent(inout) :: HashIndex(:)
      integer, intent(in) :: nJ(nel), partInd

      ! remove a determinant from the hashtable
      call remove_hash_table_entry(HashIndex, nJ, PartInd)
      ! Add to "freeslot" list so it can be filled in.
      iEndFreeSlot = iEndFreeSlot + 1
      FreeSlot(iEndFreeSlot) = PartInd
      ! Mark it as removed
      call set_flag(CurrentDets(:, PartInd), flag_removed, .true.)

    end subroutine RemoveHashDet

    function get_diagonal_matel(nI, ilut) result(diagH)
      ! Get the diagonal element for a determinant nI with ilut representation ilut

      ! In:  nI        - The determinant to evaluate
      !      ilut      - Bit representation (only used with HPHF
      ! Ret: diagH     - The diagonal matrix element
      implicit none
      integer, intent(in) :: nI(nel)
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      real(dp) :: diagH

      if(tHPHF) then
         diagH = hphf_diag_helement(nI, ilut)
      else
         diagH = get_helement(nI,nI,0)
      endif

    end function get_diagonal_matel

    subroutine CalcHashTableStats(TotWalkersNew, iter_data)

        use DetBitOps, only: FindBitExcitLevel
        use hphf_integrals, only: hphf_off_diag_helement
        use FciMCData, only: ProjEDet, CurrentDets, n_prone_dets
        use LoggingData, only: FCIMCDebug
        use bit_rep_data, only: NOffSgn

        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data

        integer :: i, j, AnnihilatedDet, lbnd, ubnd, part_type
        real(dp) :: CurrentSign(lenof_sign)
        real(dp) :: pRemove, r
        integer :: nI(nel), run, ic
        logical :: tIsStateDeterm
        real(dp) :: hij, scaledOccupiedThresh
        character(*), parameter :: t_r = 'CalcHashTableStats'

        if (.not. bNodeRoot) return

        TotParts = 0.0_dp
        norm_psi_squared = 0.0_dp
        norm_semistoch_squared = 0.0_dp
        iHighestPop = 0
        AnnihilatedDet = 0
        tIsStateDeterm = .false.
        InstNoAtHf = 0.0_dp
        n_prone_dets = 0

        if (TotWalkersNew > 0) then
            do i=1,TotWalkersNew

                call extract_sign(CurrentDets(:,i),CurrentSign)
                if (tSemiStochastic) tIsStateDeterm = test_flag(CurrentDets(:,i), flag_deterministic)

                if (IsUnoccDet(CurrentSign) .and. (.not. tIsStateDeterm)) then
                    if(.not. tAccumEmptyDet(CurrentDets(:,i))) AnnihilatedDet = AnnihilatedDet + 1
                else

                   ! count the number of walkers that are single-spawns at the threshold
                   if(t_prone_walkers) then
                      if(test_flag(CurrentDets(:,i), flag_prone)) n_prone_dets = n_prone_dets + 1
                   endif

                   if(tEScaleWalkers) then
                      scaledOccupiedThresh = OccupiedThresh * scaleFunction(det_diagH(i))
                   else
                      scaledOccupiedThresh = OccupiedThresh
                   endif
                    do j=1, lenof_sign
                        run = part_type_to_run(j)
                        if (.not. tIsStateDeterm) then
                            if ((abs(CurrentSign(j)) > 1.e-12_dp) .and. (abs(CurrentSign(j)) < scaledOccupiedThresh)) then
                                !We remove this walker with probability 1-RealSignTemp
                                pRemove=(scaledOccupiedThresh-abs(CurrentSign(j)))/scaledOccupiedThresh
                                r = genrand_real2_dSFMT ()
                                if (pRemove  >  r) then
                                   !Remove this walker
                                   NoRemoved(run) = NoRemoved(run) + abs(CurrentSign(j))
                                   iter_data%nremoved(j) = iter_data%nremoved(j) &
                                        + abs(CurrentSign(j))
                                   CurrentSign(j) = 0.0_dp
                                   call nullify_ilut_part(CurrentDets(:,i), j)
                                   call decode_bit_det(nI, CurrentDets(:,i))
                                   if (IsUnoccDet(CurrentSign) .and. .not. tAccumEmptyDet(CurrentDets(:,i))) then
                                      call RemoveHashDet(HashIndex, nI, i)
                                      ! also update both the number of annihilated dets
                                      AnnihilatedDet = AnnihilatedDet + 1
                                      ! and the number of holes
                                      HolesInList = HolesInList + 1
                                   end if
                                else
                                   NoBorn(run) = NoBorn(run) + scaledOccupiedThresh - abs(CurrentSign(j))
                                   iter_data%nborn(j) = iter_data%nborn(j) &
                                        + scaledOccupiedThresh - abs(CurrentSign(j))
                                   CurrentSign(j) = sign(scaledOccupiedThresh, CurrentSign(j))
                                   call encode_part_sign (CurrentDets(:,i), CurrentSign(j), j)
                                end if
                             end if
                          end if
                       end do

                    TotParts = TotParts + abs(CurrentSign)

                    call addNormContribution(CurrentSign, tIsStateDeterm)

                    if (tCheckHighestPop) then
                        ! If this option is on, then we want to compare the
                        ! weight on each determinant to the weight at the HF
                        ! determinant.
                        !
                        ! Record the highest weighted determinant on each
                        ! processor. If double run, only consider set 1 to keep things simple.
                        do run = 1, inum_runs
                            lbnd = min_part_type(run)
                            ubnd = max_part_type(run)
                            if (abs_sign(CurrentSign(lbnd:ubnd)) > iHighestPop(run)) then
                                iHighestPop(run) = int(abs_sign(CurrentSign(lbnd:ubnd)))
                                HighestPopDet(:,run)=CurrentDets(:,i)
                            end if
                        end do
                    end if
                end if

                if (tFillingStochRDMonFly) then
                    if (IsUnoccDet(CurrentSign) .and. (.not. tIsStateDeterm)) then
                        if (DetBitEQ(CurrentDets(:,i), iLutHF_True, NIfDBO)) then
                            ! We have to do this such that AvNoAtHF matches up with AvSign.
                            ! AvSign is extracted from CurrentH, and if the HFDet is unoccupied
                            ! at this moment during annihilation, it's CurrentH entry is removed
                            ! and the averaging information in it is lost.
                            ! In some cases (a successful spawning event) a CurrentH entry will
                            ! be recreated, but with AvSign 0, so we must match this here.
                            AvNoAtHF = 0.0_dp
                            IterRDM_HF = Iter + 1
                        end if
                    end if
                end if

                ! This InstNoAtHF call must be placed at the END of the routine
                ! as the value of CurrentSign can change during it!
                if (DetBitEQ(CurrentDets(:,i), iLutHF_True, NIfDBO)) then
                    InstNoAtHF=CurrentSign
                end if

            end do
        end if

        IFDEBUGTHEN(FCIMCDebug,6)
            write(6,*) "After annihilation: "
            write(6,*) "TotWalkersNew: ", TotWalkersNew
            write(6,*) "AnnihilatedDet: ", AnnihilatedDet
            write(6,*) "HolesInList: ", HolesInList
            write(iout,"(A,I12)") "Walker list length: ",TotWalkersNew
            write(iout,"(A)") "TW: Walker  Det"
            do j = 1, int(TotWalkersNew,sizeof_int)
                CurrentSign = transfer(CurrentDets(NOffSgn:NOffSgn+lenof_sign-1,j),CurrentSign)
                write(iout, "(A,I10,a)", advance='no') 'TW:', j, '['
                do part_type = 1, lenof_sign
                    write(iout, "(f16.3)", advance='no') CurrentSign(part_type)
                end do
                call WriteBitDet(iout,CurrentDets(:,j),.true.)
                call neci_flush(iout)
            enddo
        ENDIFDEBUG

        ! RT_M_Merge: i have to ask werner why this check makes sense
        ! AnnihilatedDet is only affected by empty dets and emptying a det increses HolesInList
        ! But adding a new det decreases HolesInList and does not affect AnnihilatedDet ->?
        if (AnnihilatedDet /= HolesInList) then
            write(6,*) "TotWalkersNew: ", TotWalkersNew
            write(6,*) "AnnihilatedDet: ", AnnihilatedDet
            write(6,*) "HolesInList: ", HolesInList
            write(6,*) "iStartFreeSlot, iEndFreeSlot:", iStartFreeSlot, iEndFreeSlot
            write(6,*) "TotParts: ", TotParts
            call neci_flush(6)
            call stop_all(t_r, "Error in determining annihilated determinants")
        end if
    end subroutine CalcHashTableStats

    subroutine addNormContribution(CurrentSign, tIsStateDeterm)
      implicit none
      real(dp), intent(in) :: CurrentSign(lenof_sign)
      logical, intent(in) :: tIsStateDeterm
      integer :: run

#if defined(CMPLX_)
      do run = 1, inum_runs
         norm_psi_squared(run) = norm_psi_squared(run) + sum(CurrentSign(min_part_type(run):max_part_type(run))**2)
         if (tIsStateDeterm) then
            norm_semistoch_squared(run) = norm_semistoch_squared(run) &
                 + sum(CurrentSign(min_part_type(run):max_part_type(run))**2)
         endif
      enddo

#else
      norm_psi_squared = norm_psi_squared + CurrentSign**2
      if (tIsStateDeterm) norm_semistoch_squared = norm_semistoch_squared + CurrentSign**2
#endif
  end subroutine addNormContribution

!------------------------------------------------------------------------------------------!

  !> Gauge if a load balancing step shall be taken given the current load-imbalance
  !! measure lt_imb
  !> @param[in] lt_imb  current load imbalance measure: Time lost due to load imbalance
  !!                    during the last 100 iterations divided by the total time taken for these
  !> @result t_lb  true if a load balancing step is justified
  function need_load_balancing(lt_imb) result(t_lb)
      real(dp), intent(in) :: lt_imb
      logical :: t_lb
      
      real(dp), save :: last_imb = 0.0_dp
      logical, save :: last_t_lb = .false.

      ! In the cycle immediately after a load balancing step, we do not load balance
      ! again, but instead log the imbalance measure for comparison
      if(last_t_lb) then
          last_imb = lt_imb
          t_lb = .false.
          last_t_lb = .false.
      else
          ! Load balance if the measure is sufficiently high (both absolute and relative
          ! to what we had after the last load balancing)
          t_lb = lt_imb > max(0.1_dp, 2*last_imb)
          last_t_lb = t_lb
      end if
  end function need_load_balancing

end module
