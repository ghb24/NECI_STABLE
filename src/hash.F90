#include "macros.h"

module hash

    use bit_rep_data, only: NIfTot, extract_sign, test_flag, flag_deterministic
    use bit_reps, only: set_flag, decode_bit_det
    use constants
    use FciMCData , only : hash_iter, hash_shift, RandomHash, RandomHash2, HFDet, ll_node
    use Parallel_neci , only : nNodes
    use constants , only : int64,sizeof_int
    use csf_data, only: csf_orbital_mask
    use Systemdata, only: nel, tCSF, nBasis
    use CalcData, only: tUniqueHFNode, tSemiStochastic

    implicit none

    contains

    pure function DetermineDetNode (nel_loc, nI, iIterOffset) result(node)

        ! Depending on the Hash, determine which node determinant nI
        ! belongs to in the DirectAnnihilation scheme. NB FCIMC has each processor as a separate logical node.
        !
        ! In:  nI   - Integer ordered list for the determinant
        ! In:  iIterOffset - Offset this iteration by this amount
        ! Out: proc - The (0-based) processor index.

        implicit none
        integer, intent(in) :: nel_loc
        integer, intent(in) :: nI(nel_loc)
        integer, intent(in) :: iIterOffset
        integer :: node
        
        integer :: i
        integer(int64) :: acc
        integer(int64) :: offset
        integer(int64), parameter :: large_prime = 1099511628211_int64

        ! If we are assigning the HF to a unique processor, then put it 
        ! on the last available processor.
        if (size(nI) == size(HFDet)) then
            if (tUniqueHFNode .and. all(nI == HFDet)) then
                node = nNodes-1
                return
            end if
        end if

!sum(nI) ensures that a random number is generated for each different nI, which is then added to the iteration,
!  and the result shifted.  Consequently, the less significant bits (but not the least, as these have been shifted away)
!  for each nI will change on average once per 2^hash_shift iterations, but the change spread throughout the different iters.
        if (hash_iter>0) then  !generate a hash to work out an offset.  Probably very inefficient
           acc = 0
           do i = 1, nel_loc
               acc = (large_prime * acc) + &
                       (RandomHash(mod(iand(nI(i), csf_orbital_mask)-1,nBasis)+1) * i)
           enddo
           offset=ishft(abs(ishft(acc+hash_iter+iIterOffset, -hash_shift) ),-4)
        else
           offset=0
        endif
        acc = 0
        if(tCSF) then
            do i = 1, nel_loc
                acc = (large_prime * acc) + &
                        (RandomHash(mod(iand(nI(i), csf_orbital_mask)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        else
            do i = 1, nel_loc
                acc = (large_prime * acc) + &
                        (RandomHash(mod(nI(i)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        endif

        ! If the last available processor is being used for the HF det, then
        ! we can only use the the remaining (nNodes-1) processors to
        ! distribute the rest ofthe particles
        if (tUniqueHFNode) then
            node = int(abs(mod(acc, int(nNodes-1, int64))),sizeof_int)
        else
            node = int(abs(mod(acc, int(nNodes, int64))),sizeof_int)
        end if

    end function DetermineDetNode

    ! Routin to find the correct position in the hash table
    pure function FindWalkerHash(nJ, HashIndexLength) result(hashInd)
        implicit none
        integer, intent(in) :: nJ(nel)
        integer, intent(in) :: HashIndexLength
        integer :: hashInd
        integer :: i
        integer(int64) :: hash
        hash = 0
!        write(6,*) "nJ: ",nJ(:)
        if(tCSF) then
            do i = 1, nel
                hash = (1099511628211_int64 * hash) + &
                        int(RandomHash2(mod(iand(nJ(i), csf_orbital_mask)-1,nBasis)+1) * i,int64)
            enddo
        else
            do i = 1, nel
                hash = (1099511628211_int64 * hash) + &
                        int(RandomHash2(nJ(i))*i,int64)

!                        (RandomHash(mod(nI(i)+offset-1,int(nBasis,int64))+1) * i)
            enddo
        endif
        hashInd = int(abs(mod(hash, int(HashIndexLength, int64))),sizeof_int)+1

    end function FindWalkerHash

    FUNCTION CreateHash(DetCurr)
        INTEGER :: DetCurr(NEl),i
        INTEGER(KIND=int64) :: CreateHash

        CreateHash=0
        do i=1,NEl
!            CreateHash=13*CreateHash+i*DetCurr(i)
            CreateHash=(1099511628211_int64*CreateHash)+i*DetCurr(i)
            
!            CreateHash=mod(1099511628211*CreateHash,2**64)
!            CreateHash=XOR(CreateHash,DetCurr(i))
        enddo
!        WRITE(6,*) CreateHash
        RETURN

    END FUNCTION CreateHash

    pure subroutine init_hash_table(hash_table)

        type(ll_node), pointer, intent(inout) :: hash_table(:)
        integer :: i

        do i = 1, size(hash_table)
            hash_table(i)%ind = 0
            nullify(hash_table(i)%next)
        end do

    end subroutine init_hash_table

    pure subroutine reset_hash_table(hash_table)

        type(ll_node), pointer, intent(inout) :: hash_table(:)
        type(ll_node), pointer :: curr, prev
        integer :: i

        ! Reset the hash index array.
        do i = 1, size(hash_table)
            curr => hash_table(i)%next
            prev => hash_table(i)
            prev%ind = 0
            nullify(prev%next)
            do while (associated(curr))
                prev => curr
                curr => curr%next
                deallocate(prev)
            end do
        end do
        nullify(curr)
        nullify(prev)

    end subroutine reset_hash_table

    subroutine fill_in_hash_table(hash_table, table_length, walker_list, list_length, ignore_unocc)

        ! This assumes that the input hash table is clear (use reset_hash_table)
        ! and that there are no repeats in walker_list.

        ! If ignore_unocc dets is input as .true. then unoccupied determinants
        ! will not be included in the hash table, unless they are core
        ! determinants.

        integer, intent(in) :: table_length, list_length
        type(ll_node), pointer, intent(inout) :: hash_table(:)
        integer(n_int), intent(in) :: walker_list(0:,:)
        logical, intent(in) :: ignore_unocc

        type(ll_node), pointer :: temp_node
        integer :: i, hash_val, nI(nel)
        real(dp) :: real_sign(lenof_sign)
        logical :: tCoreDet

        tCoreDet = .false.

        do i = 1, list_length
            ! If the ignore_unocc option is true then we don't want to include
            ! unoccupied determinants in the hash table, unless they're
            ! deterministic states.
            if (ignore_unocc) then
                if (tSemiStochastic) tCoreDet = test_flag(walker_list(:,i), flag_deterministic)
                call extract_sign(walker_list(:,i), real_sign)
                if (IsUnoccDet(real_sign) .and. (.not. tCoreDet)) cycle
            end if

            call decode_bit_det(nI, walker_list(:,i))
            ! Find the hash value corresponding to this determinant.
            hash_val = FindWalkerHash(nI, table_length)
            ! Point to the start of the linked list corresponding to hash_val.
            temp_node => hash_table(hash_val)
            if (temp_node%ind == 0) then
                ! If we get here then this linked list is currently empty.
                temp_node%ind = i
            else
                do while (associated(temp_node%next))
                    temp_node => temp_node%next
                end do
                allocate(temp_node%next)
                nullify(temp_node%next%next)
                temp_node%next%ind = i
            end if
            nullify(temp_node)
        end do

    end subroutine fill_in_hash_table

    subroutine remove_hash_table_entry(hash_table, nI, ind)

        ! Find and remove the entry in hash_table corresponding to nI, which
        ! must have index ind in the hash table. If not found then an error
        ! will be thrown.

        type(ll_node), pointer, intent(inout) :: hash_table(:)
        integer, intent(in) :: nI(:)
        integer, intent(in) :: ind

        integer :: hash_val
        type(ll_node), pointer :: prev, curr
        logical :: found

        found = .false.

        ! Find the hash value corresponding to this determinant.
        hash_val = FindWalkerHash(nI, size(hash_table))
        ! Point at the start of the linked list for this hash value.
        curr => hash_table(hash_val)
        prev => null()
        ! Loop over all entries in the linked list until we find the one equal
        ! to ind, the entry that we want to remove.
        do while (associated(curr))
            if (curr%ind == ind) then
                ! If this is the state to be removed.
                found = .true.
                call remove_node(prev, curr)
                exit
            end if
            prev => curr
            curr => curr%next
        end do

        ASSERT(found)

    end subroutine remove_hash_table_entry

    pure subroutine remove_node(prev, curr)

        ! On input, prev should point to the the node before curr in the linked list,
        ! or should point to null if curr is the first node in the list. curr should
        ! point to the node which is to be removed.

        ! On output, both prev and curr will be nullified.

        type(ll_node), pointer, intent(inout) :: prev, curr
        type(ll_node), pointer :: temp_node

        if (associated(prev)) then
            ! If not the first state in the list.
            prev%next => curr%next
            deallocate(curr)
        else
            ! If the first state in the list.
            if (associated(curr%next)) then
                ! If the first but not the only state in the list.
                ! Move the details of the second entry in the list to the
                ! first entry, and then deallocate the second entry.
                curr%ind = curr%next%ind
                temp_node => curr%next
                curr%next => curr%next%next
                deallocate(temp_node)
            else
                ! If the first and only state in the list.
                curr%ind = 0
                curr%next => null()
            end if
        end if

        prev => null()
        curr => null()
        temp_node => null()

    end subroutine remove_node

    subroutine add_hash_table_entry(hash_table, ind, hash_val)

        ! Add an entry of ind into hash_table at an index specified by hash_val.

        type(ll_node), pointer, intent(inout) :: hash_table(:)
        integer, intent(in) :: ind
        integer, intent(in) :: hash_val

        type(ll_node), pointer :: temp_node
        
        ! Point to the start of the linked list corresponding to hash_val.
        temp_node => hash_table(hash_val)
        if (temp_node%ind == 0) then
            ! If here then this linked list is empty.
            ! Just need to add the index to the hash table and exit.
            temp_node%ind = ind
        else
            ! If here then there is at least one entry in this linked list.
            ! Cycle to the end of the linked list, and add this new entry on
            ! the end.
            do while (associated(temp_node%next))
                temp_node => temp_node%next
            end do
            allocate(temp_node%next)
            nullify(temp_node%next%next)
            temp_node%next%ind = ind
        end if

        nullify(temp_node)

    end subroutine add_hash_table_entry
      
end module hash
