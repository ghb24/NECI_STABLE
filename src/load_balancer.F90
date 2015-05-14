#include "macros.h"
module load_balance

    use global_det_data, only: set_det_diagH, get_iter_occ, set_spawn_rate, &
                               global_determinant_data, set_part_init_time, &
                               inc_spawn_count, get_spawn_count, pos_spawn_cnt
    use bit_rep_data, only: flag_trial, flag_connected
    use CalcData, only: tUniqueHFNode, tSemiStochastic
    use FciMCData, only: HFDet, hash_iter, hash_shift, TotImagTime, &
                         current_trial_amps, HashIndex, Hii, &
                         FreeSlot, CurrentDets, MaxWalkersPart, tTrialHash, &
                         ntrial_excits, iStartFreeSlot, iEndFreeSlot
    use searching, only: hash_search_trial, bin_search_trial
    use SystemData, only: tCSF, nBasis, nel, tHPHF
    use hphf_integrals, only: hphf_diag_helement
    use cont_time_rates, only: spawn_rate_full
    use Determinants, only: get_helement
    use csf_data, only: csf_orbital_mask
    use bit_reps, only: set_flag
    use Parallel_neci
    use constants
    use util_mod
    use hash

    implicit none

    integer, allocatable :: RandomOrbIndex(:), LoadBalanceMapping(:)
    integer(TagIntType) :: lb_tag
    integer :: balance_blocks
    logical :: tLoadBalanceBlocks

    ! TODO:
    ! - Initialise mapping
    ! - Modify DetermineDetNode to use the mapping
    ! - Add mechanism to shuffle walkers around
    ! - Add redistribution of walkers
    ! - Integrate with POPSFILES. Need to output the mapping for restarts.
    ! - Ensure that we never re-load balance once using SemiStochastic
    ! - Consider if we want to load balance based on sync-time rather than
    !   particles?

contains

    subroutine init_load_balance()

        ! Initialise the load balancing.
        !
        ! n.b. The initialisation of RandomOrbIndex remains in SetupParameters
        !      to preserve sequencing, which maintains testcode results.

        integer :: oversample_factor, ierr, i
        character(*), parameter :: t_r = 'init_load_balance'

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

    subroutine clean_load_balance()

        character(*), parameter :: this_routine = 'clean_load_balance'

        if (allocated(LoadBalanceMapping)) then
            deallocate(LoadBalanceMapping)
            log_dealloc(lb_tag)
        end if

    end subroutine

    subroutine adjust_load_balance()

        integer(int64) :: block_parts(balance_blocks)
        integer(int64) :: block_parts_all(balance_blocks)
        integer(int64) :: proc_parts(0:nProcessors-1)
        integer(int64) :: smallest_size
        integer :: j, proc, nblocks, det(nel), block
        integer :: min_parts, max_parts, min_proc, max_proc
        integer :: smallest_block
        real(dp) :: sgn(lenof_sign), avg_parts
        logical :: unbalanced
        character(*), parameter :: this_routine = 'adjust_load_balance'

        ! TODO: Need to ensure we don't move around the semi-stochastic sites,
        ! or we need to be a bit more clever!!!
        ! TODO: What happens if we move reference sites around?
        if (tSemiStochastic) &
            call stop_all(this_routine, "Not yet implemented")

        ! Count the number of particles inside each of the blocks
        block_parts = 0
        do j = 1, int(TotWalkers, sizeof_int)

            call extract_sign(CurrentDets(:,j), sgn)
            if (IsUnoccDet(sgn)) cycle

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
                avg_parts = real(sum(proc_parts), dp) / real(nProcessors, dp)
                min_proc = minloc(proc_parts, dim=1)
                min_parts = proc_parts(min_proc)
                max_proc = maxloc(proc_parts, dim=1)
                max_parts = proc_parts(max_proc)

                ! Create a list of the blocks associated with the most
                ! heavily utilised processor in increasing size order.
                smallest_block = 0
                smallest_size = -1
                do block = 1, balance_blocks
                    if (LoadBalanceMapping(block) == max_proc) then
                        if (block_parts(block) < smallest_size .or. &\
                                smallest_size == -1) then
                            smallest_block = block
                            smallest_size = block_parts(block)
                        end if
                    end if
                end do

                ! If moving a block of the smallest size between the largest
                ! and the smallest is a helpful thing to do, then move it!
                if ((abs(min_parts + smallest_size - avg_parts) < abs(min_parts - avg_parts)) .and. &
                    (abs(max_parts - smallest_size - avg_parts) < abs(max_parts - avg_parts))) then
                    unbalanced = .true.
                else
                    unbalanced = .false.
                end if
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
        integer :: det(nel)
        real(dp) :: sgn(lenof_sign)
        
        ! A tag is used to identify this send/recv pair over any others
        integer, parameter :: mpi_tag_nsend = 223456
        integer, parameter :: mpi_tag_dets = 223457

        src_proc = LoadBalanceMapping(block)

        ! Provide some feedback to the user.
        if (iProcIndex == root) then
            write(6,*) 'Moving load balancing block ', block, &
                       'from processor ', src_proc, ' to ', tgt_proc
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

            ! Adjust local counters
            TotWalkers = TotWalkers - nsend

        else if (iProcIndex == tgt_proc) then

            ! Receive walkers!
            call MPIRecv(nsend, 1, src_proc, mpi_tag_nsend, ierr)
            nelem = nsend * (1 + NIfTot)
            call MPIRecv(SpawnedParts, nelem, src_proc, mpi_tag_dets, ierr)

            do j = 1, nsend
                call decode_bit_det(det, CurrentDets(:,j))
                hash_val = FindWalkerHash(det, MaxWalkersPart)
                call AddNewHashDet(TotWalkers, SpawnedParts(:, nsend), &
                                   hash_val, det)
            end do

            ! Todo: remember to regenerate global stored data!

            ! Adjust local counters
            TotWalkers = TotWalkers + nsend

        end if

        ! And synchronise when everything is done
        call MPIBarrier(ierr)

    end subroutine


    pure function DetermineDetNode (nel_loc, nI, iIterOffset) result(node)

        ! Depending on the Hash, determine which node determinant nI
        ! belongs to in the DirectAnnihilation scheme. NB FCIMC has each
        ! processor as a separate logical node.

        ! In:  nI   - Integer ordered list for the determinant
        ! In:  iIterOffset - Offset this iteration by this amount
        ! Out: proc - The (0-based) processor index.

        ! --> This function takes the calculated block (from get_det_block)
        !     and converts it via a simple lookup into the required node


        integer, intent(in) :: nel_loc
        integer, intent(in) :: nI(nel_loc)
        integer, intent(in) :: iIterOffset
        integer :: node
        
        integer :: block

        block = get_det_block(nel_loc, nI, iIterOffset)

        ! Look up the relevant node in the block-mapping.
        node = LoadBalanceMapping(block)

    end function


    pure function get_det_block(nel_loc, nI, iIterOffset) result(block)

        ! Depending on the Hash, determine which node determinant nI
        ! belongs to in the DirectAnnihilation scheme. NB FCIMC has each
        ! processor as a separate logical node.

        ! In:  nI   - Integer ordered list for the determinant
        ! In:  iIterOffset - Offset this iteration by this amount
        ! Out: proc - The (0-based) processor index.

        ! --> This function calculates the hash, and takes it mod the number
        !     of blocks in use

        integer, intent(in) :: nel_loc
        integer, intent(in) :: nI(nel_loc)
        integer, intent(in) :: iIterOffset
        integer :: block
        
        integer :: i
        integer(int64) :: acc
        integer(int64) :: offset
        integer(int64), parameter :: large_prime = 1099511628211_int64

        ! If we are assigning the HF to a unique processor, then put it 
        ! on the last available processor.
        if (size(nI) == size(HFDet)) then
            if (tUniqueHFNode .and. all(nI == HFDet)) then
                block = nNodes-1
                return
            end if
        end if

        ! sum(nI) ensures that a random number is generated for each different
        ! nI, which is then added to the iteration, and the result shifted.
        ! Consequently, the less significant bits (but not the least, as these
        ! have been shifted away) for each nI will change on average once per
        ! 2^hash_shift iterations, but the change spread throughout the
        ! different iters.
        ! Generate a hash to work out an offset.  Probably very inefficient.
        if (hash_iter>0) then  
           acc = 0
           do i = 1, nel_loc
               acc = (large_prime * acc) + &
                       (RandomOrbIndex(mod(iand(nI(i), csf_orbital_mask)-1,nBasis)+1) * i)
           enddo
           offset=ishft(abs(ishft(acc+hash_iter+iIterOffset, -hash_shift) ),-4)
        else
           offset=0
        endif
        acc = 0
        if(tCSF) then
            do i = 1, nel_loc
                acc = (large_prime * acc) + &
                        (RandomOrbIndex(mod(iand(nI(i), csf_orbital_mask)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        else
            do i = 1, nel_loc
                acc = (large_prime * acc) + &
                        (RandomOrbIndex(mod(nI(i)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        endif

        ! If the last available processor is being used for the HF det, then
        ! we can only use the the remaining (nNodes-1) processors to
        ! distribute the rest ofthe particles
        if (tUniqueHFNode) then
            block = int(abs(mod(acc, int(balance_blocks-1, int64))),sizeof_int)
        else
            block = int(abs(mod(acc, int(balance_blocks, int64))),sizeof_int)
        end if

    end function

    !
    ! --- This following function really wants to be in Annihilation, but
    !     putting it here is easier for circular dependencie resolution.

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

end module
