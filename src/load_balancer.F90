#include "macros.h"
module load_balance

    use CalcData, only: tUniqueHFNode, tSemiStochastic, tTruncInitiator, &
                        tCheckHighestPop, tEnhanceRemainder, OccupiedThresh, &
                        InitiatorOccupiedThresh, tContTimeFCIMC, &
                        tContTimeFull, tTrialWavefunction, tInitOccThresh
    use global_det_data, only: global_determinant_data, get_iter_occ, &
                               set_det_diagH, set_part_init_time, &
                               inc_spawn_count, set_spawn_rate
    use bit_rep_data, only: flag_initiator, NIfDBO, flag_has_been_initiator, &
                            flag_connected, flag_trial
    use bit_reps, only: set_flag, nullify_ilut_part, clear_has_been_initiator,&
                        encode_part_sign, nullify_ilut
    use FciMCData, only: HashIndex, FreeSlot, CurrentDets, iter_data_fciqmc, &
                         tFillingStochRDMOnFly
    use searching, only: hash_search_trial, bin_search_trial
    use Determinants, only: get_helement, write_det
    use rdm_filling, only: det_removed_fill_diag_rdm
    use hphf_integrals, only: hphf_diag_helement
    use cont_time_rates, only: spawn_rate_full
    use SystemData, only: nel, tHPHF
    use DetBitOps, only: DetBitEq
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
        character(*), parameter :: t_r = this_routine

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
        log_alloc(LoadBalanceMapping, lb_tag, ierr)

        ! Generate a uniform mapping(by default)
        do i = 1, balance_blocks
            LoadBalanceMapping(i) = int((i - 1) / oversample_factor)
        end do

    end subroutine

    subroutine pops_init_balance_blocks(pops_blocks)

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
        ASSERT(allocated(LoadBalanceMapping))
        if (pops_blocks == -1) then
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
            if (IsUnoccDet(sgn)) cycle

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
        if (.not. all(mapping_test_all == 1)) then
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
        integer :: j, proc, nblocks, det(nel), block, TotWalkersTmp
        integer :: min_parts, max_parts, min_proc, max_proc
        integer :: smallest_block
        real(dp) :: sgn(lenof_sign), avg_parts
        logical :: unbalanced
        character(*), parameter :: this_routine = 'adjust_load_balance'

        ! TODO: Need to ensure we don't move around the semi-stochastic sites,
        ! or we need to be a bit more clever!!!
        ! TODO: What happens if we move reference sites around?
        ASSERT(.not. tSemiStochastic)
        ASSERT(.not. tFillingStochRDMOnFly)

        ! Count the number of particles inside each of the blocks
        block_parts = 0
        do j = 1, int(TotWalkers, sizeof_int)

            call extract_sign(CurrentDets(:,j), sgn)
            if (IsUnoccDet(sgn)) then
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
                min_parts = proc_parts(min_proc)
                max_proc = maxloc(proc_parts, dim=1) - 1
                max_parts = proc_parts(max_proc)
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
            if (.not. unbalanced) exit

            ! Broadcast the parameters for the change!
            call MPIBCast(min_proc)
            call MPIBCast(max_proc)
            call MPIBcast(smallest_block)

            ! Move the block from where it is to the currently least worked
            ! processor
            call move_block(smallest_block, min_proc)

        end do

        write(6, '("Load balancing distribution:")')
        write(6, '("node #, particles")')
        do j = 0, nNodes - 1
            write(6,'(i7,i9)') j, proc_parts(j)
        end do
        write(6,*) '--'

        ! TODO: Only call this if we have made changes!
        TotWalkersTmp = int(TotWalkers, sizeof_int)
        call CalcHashTableStats(TotWalkersTmp, iter_data)
        TotWalkers = TotWalkersTmp

        !   -- Test if sufficiently uniform
        !   -- If not, pick largest, and smallest, sites
        !   -- Transfer the largest block that will not take either of the
        !      processors past the average number.
        !   -- If there is no such block, use the smallest, ensuring that
        !      we still make improvements.
        !   -- If there is no improving swap, do nothing!

    end subroutine

        
    subroutine move_block(block, tgt_proc)

        integer, intent(in) :: block, tgt_proc
        integer :: src_proc, ierr, nsend, nelem, j, det_block, hash_val
        integer :: det(nel), TotWalkersTmp
        real(dp) :: sgn(lenof_sign)
        
        ! A tag is used to identify this send/recv pair over any others
        integer, parameter :: mpi_tag_nsend = 223456
        integer, parameter :: mpi_tag_dets = 223457

        src_proc = LoadBalanceMapping(block)

        ! Provide some feedback to the user.
        if (iProcIndex == root) then
            write(6,'(a,i6,a,i6,a,i6)') 'Moving load balancing block ', &
                     block, ' from processor ', src_proc, ' to ', tgt_proc
        end if

        if (iProcIndex == src_proc) then

            ! Loop over the available walkers, and broadcast them to the
            ! target processor. Use the SpawnedParts array as a buffer.
            nsend = 0
            do j = 1, int(TotWalkers, sizeof_int)

                ! Skip unoccupied sites (non-contiguous)
                call extract_sign(CurrentDets(:,j), sgn)
                if (IsUnoccDet(sgn)) cycle

                call decode_bit_det(det, CurrentDets(:,j))
                det_block = get_det_block(nel, det, 0)
                if (det_block == block) then
                    nsend = nsend + 1
                    SpawnedParts(:,nsend) = CurrentDets(:,j)

                    ! Remove the det from the main list.
                    call nullify_ilut(CurrentDets(:,j))
                    call remove_hash_table_entry(HashIndex, det, j)
                    iEndFreeSlot = iEndFreeSlot + 1
                    FreeSlot(iEndFreeSlot) = j
                end if
            end do

            ! And send the data to the relevant (target) processor
            nelem = nsend * (1 + NIfTot)
            call MPISend(nsend, 1, tgt_proc, mpi_tag_nsend, ierr)
            call MPISend(SpawnedParts(:, 1:nsend), nelem, tgt_proc, &
                         mpi_tag_dets, ierr)

            ! We have now created lots of holes in the main list
            HolesInList = HolesInList + nsend

        else if (iProcIndex == tgt_proc) then

            ! Receive walkers!
            call MPIRecv(nsend, 1, src_proc, mpi_tag_nsend, ierr)
            nelem = nsend * (1 + NIfTot)
            call MPIRecv(SpawnedParts, nelem, src_proc, mpi_tag_dets, ierr)

            do j = 1, nsend
                call decode_bit_det(det, SpawnedParts(:,j))
                call extract_sign(SpawnedParts(:,j), sgn)
                hash_val = FindWalkerHash(det, size(HashIndex))

                ! n.b. Ensure that Totwalkers passed in always has the correct
                !      type even on 32-bit machines
                TotWalkersTmp = TotWalkers
                call AddNewHashDet(TotWalkersTmp, SpawnedParts(:, j), &
                                   hash_val, det)
                TotWalkers = TotWalkersTmp
            end do

            ! We have filled in some of the holes in the list (possibly all)
            ! and possibly extended the list
            HolesInList = max(0, HolesInList - nsend)

        end if

        ! Adjust the load balancing mapping
        LoadBalanceMapping(block) = tgt_proc

        ! And synchronise when everything is done
        call MPIBarrier(ierr)

    end subroutine


    subroutine AddNewHashDet(TotWalkersNew, iLutCurr, DetHash, nJ)

        ! Add a new determinant to the main list. This involves updating the
        ! list length, copying it across, updating its flag, adding its diagonal
        ! helement (if neccessary). We also need to update the hash table to
        ! point at it correctly.

        integer, intent(inout) :: TotWalkersNew 
        integer(n_int), intent(inout) :: iLutCurr(0:NIfTot)
        integer, intent(in) :: DetHash, nJ(nel)
        integer :: DetPosition
        HElement_t :: HDiag
        real(dp) :: trial_amps(ntrial_excits)
        logical :: tTrial, tCon
        character(len=*), parameter :: t_r = "AddNewHashDet"

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
                call stop_all(t_r, "Not enough memory to merge walkers into main list. Increase MemoryFacPart")
            end if
            CurrentDets(:,DetPosition) = iLutCurr(:)
        end if

        ! Calculate the diagonal hamiltonian matrix element for the new particle to be merged.
        if (tHPHF) then
            HDiag = hphf_diag_helement (nJ,CurrentDets(:,DetPosition))
        else
            HDiag = get_helement (nJ, nJ, 0)
        end if

        ! For the RDM code we need to set all of the elements of CurrentH to 0,
        ! except the first one, holding the diagonal Hamiltonian element.
        global_determinant_data(:,DetPosition) = 0.0_dp
        call set_det_diagH(DetPosition, real(HDiag,dp) - Hii)

        ! Store the iteration, as this is the iteration on which the particle
        ! is created
        call set_part_init_time(DetPosition, TotImagTime)

        ! There is at least one spawning count here
        call inc_spawn_count(DetPosition)

        ! If using a trial wavefunction, search to see if this state is in
        ! either the trial or connected space. If so, *_search_trial returns
        ! the corresponding amplitude, which is stored.
        if (tTrialWavefunction) then
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

        ! Add the new determinant to the hash table.
        call add_hash_table_entry(HashIndex, DetPosition, DetHash)

    end subroutine AddNewHashDet

    subroutine CalcHashTableStats(TotWalkersNew, iter_data)

        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: i, j, AnnihilatedDet, lbnd, ubnd
        real(dp) :: CurrentSign(lenof_sign), SpawnedSign(lenof_sign)
        real(dp) :: pRemove, r
        integer :: nI(nel), run
        logical :: tIsStateDeterm
        character(*), parameter :: t_r = 'CalcHashTableStats'

        if (.not. bNodeRoot) return

        TotParts = 0.0_dp
        norm_psi_squared = 0.0_dp
        norm_semistoch_squared = 0.0_dp
        iHighestPop = 0
        AnnihilatedDet = 0
        tIsStateDeterm = .false.
        InstNoAtHf = 0.0_dp

        if (TotWalkersNew > 0) then
            do i=1,TotWalkersNew
                call extract_sign(CurrentDets(:,i),CurrentSign)
                if (tSemiStochastic) tIsStateDeterm = test_flag(CurrentDets(:,i), flag_deterministic)

                if (IsUnoccDet(CurrentSign) .and. (.not. tIsStateDeterm)) then
                    AnnihilatedDet = AnnihilatedDet + 1 
                else
                    do j=1, lenof_sign
                        run = part_type_to_run(j)
                        if (.not. tIsStateDeterm) then
                            if (tInitOccThresh.and.test_flag(CurrentDets(:,i), flag_has_been_initiator(1)))then
                                if ((abs(CurrentSign(j)) > 0.0) .and. (abs(CurrentSign(j)) < InitiatorOccupiedThresh)) then
                                    ! We remove this walker with probability 1-RealSignTemp.
                                    pRemove = (InitiatorOccupiedThresh-abs(CurrentSign(j)))/InitiatorOccupiedThresh
                                    r = genrand_real2_dSFMT ()
                                    if (pRemove > r) then
                                        ! Remove this walker.
                                        NoRemoved(run) = NoRemoved(run) + abs(CurrentSign(j))
                                        iter_data%nremoved(j) = iter_data%nremoved(j) &
                                                              + abs(CurrentSign(j))
                                        CurrentSign(j) = 0.0_dp
                                        call nullify_ilut_part(CurrentDets(:,i), j)
                                        call decode_bit_det(nI, CurrentDets(:,i))
                                        call clear_has_been_initiator(CurrentDets(:,i),flag_has_been_initiator(1))
                                        if (IsUnoccDet(CurrentSign)) then
                                            call remove_hash_table_entry(HashIndex, nI, i)
                                            iEndFreeSlot=iEndFreeSlot+1
                                            FreeSlot(iEndFreeSlot)=i
                                        end if
                                    else if (tEnhanceRemainder) then
                                        NoBorn(run) = NoBorn(run) + InitiatorOccupiedThresh - abs(CurrentSign(j))
                                        iter_data%nborn(j) = iter_data%nborn(j) &
                                             + InitiatorOccupiedThresh - abs(CurrentSign(j))
                                        CurrentSign(j) = sign(InitiatorOccupiedThresh, CurrentSign(j))
                                        call encode_part_sign (CurrentDets(:,i), CurrentSign(j), j)
                                    end if
                                end if
                            else
                                if ((abs(CurrentSign(j)) > 0.0) .and. (abs(CurrentSign(j)) < OccupiedThresh)) then
                                !We remove this walker with probability 1-RealSignTemp
                                pRemove=(OccupiedThresh-abs(CurrentSign(j)))/OccupiedThresh
                                r = genrand_real2_dSFMT ()
                                if (pRemove  >  r) then
                                    !Remove this walker
                                    NoRemoved(run) = NoRemoved(run) + abs(CurrentSign(j))
                                    iter_data%nremoved(j) = iter_data%nremoved(j) &
                                                          + abs(CurrentSign(j))
                                    CurrentSign(j) = 0.0_dp
                                    call nullify_ilut_part(CurrentDets(:,i), j)
                                    call decode_bit_det(nI, CurrentDets(:,i))
                                    if (IsUnoccDet(CurrentSign)) then
                                        call remove_hash_table_entry(HashIndex, nI, i)
                                        iEndFreeSlot=iEndFreeSlot+1
                                        FreeSlot(iEndFreeSlot)=i
                                    end if
                                else if (tEnhanceRemainder) then
                                    NoBorn(run) = NoBorn(run) + OccupiedThresh - abs(CurrentSign(j))
                                    iter_data%nborn(j) = iter_data%nborn(j) &
                                         + OccupiedThresh - abs(CurrentSign(j))
                                    CurrentSign(j) = sign(OccupiedThresh, CurrentSign(j))
                                    call encode_part_sign (CurrentDets(:,i), CurrentSign(j), j)
                                end if
                            end if
                            end if
                            !!!!
                        end if
                    end do

                    TotParts = TotParts + abs(CurrentSign)
#if defined(__CMPLX)
                    norm_psi_squared = norm_psi_squared + sum(CurrentSign**2)
                    if (tIsStateDeterm) norm_semistoch_squared = norm_semistoch_squared + sum(CurrentSign**2)
#else
                    norm_psi_squared = norm_psi_squared + CurrentSign**2
                    if (tIsStateDeterm) norm_semistoch_squared = norm_semistoch_squared + CurrentSign**2
#endif
                    
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

                if (tFillingStochRDMonFly .and. (.not. tIsStateDeterm)) then
                    if (inum_runs == 2) then

                        if ((CurrentSign(1) == 0 .and. get_iter_occ(i, 1) /= 0) .or. &
                            (CurrentSign(inum_runs) == 0 .and. get_iter_occ(i, 2) /= 0) .or. &
                            (CurrentSign(1) /= 0 .and. get_iter_occ(i, 1) == 0) .or. &
                            (CurrentSign(inum_runs) /= 0 .and. get_iter_occ(i, 2) == 0)) then
                               
                            ! At least one of the signs has just gone to zero or just become reoccupied
                            ! so we need to consider adding in diagonal elements and connections to HF
                            ! The block that's just ended was occupied in at least one population.
                            call det_removed_fill_diag_rdm(CurrentDets(:,i), i)
                        end if
                    else
                        if (IsUnoccDet(CurrentSign)) then
                            call det_removed_fill_diag_rdm(CurrentDets(:,i), i)
                        end if
                    end if
                end if

                if (IsUnoccDet(CurrentSign) .and. (.not. tIsStateDeterm) .and. tTruncInitiator) then
                    do j=1,lenof_sign
                        if (test_flag(CurrentDets(:,i),flag_initiator(j))) then
                            !determinant was an initiator...it obviously isn't any more...
                            NoAddedInitiators(j)=NoAddedInitiators(j)-1
                        end if
                    end do
                end if

                ! This InstNoAtHF call must be placed at the END of the routine
                ! as the value of CurrentSign can change during it!
                if (DetBitEQ(CurrentDets(:,i), iLutHF_True, NIfDBO)) then
                    InstNoAtHF=CurrentSign
                end if

            end do
        end if

        if (AnnihilatedDet /= HolesInList) then
            write(6,*) "TotWalkersNew: ", TotWalkersNew
            write(6,*) "AnnihilatedDet: ", AnnihilatedDet
            write(6,*) "HolesInList: ", HolesInList
            call stop_all(t_r, "Error in determining annihilated determinants")
        end if

    end subroutine CalcHashTableStats
    

end module
