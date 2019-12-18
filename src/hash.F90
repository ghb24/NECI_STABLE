#include "macros.h"

module hash

    use FciMCData, only: hash_iter, hash_shift, RandomHash2, HFDet, ll_node
    use bit_rep_data, only: flag_deterministic, test_flag
    use bit_reps, only: extract_sign, decode_bit_det
    use Systemdata, only: nel, tCSF, nBasis
    use csf_data, only: csf_orbital_mask
    use CalcData, only: tSemiStochastic
    use constants

    implicit none

    interface hash_table_lookup
        module procedure hash_table_lookup_int_32
        module procedure hash_table_lookup_int_64
    end interface

    contains

    ! Routine to find the correct position in the hash table.
    pure function FindWalkerHash(orb_array, HashIndexLength) result(hashInd)

        integer, intent(in) :: orb_array(:)
        integer, intent(in) :: HashIndexLength
        integer :: hashInd
        integer :: i
        integer(int64) :: hash
        hash = 0
        if(tCSF) then
            do i = 1, size(orb_array)
                hash = (1099511628211_int64 * hash) + &
                        int(RandomHash2(mod(iand(orb_array(i), csf_orbital_mask)-1,nBasis)+1) * i,int64)
            enddo
        else
            do i = 1, size(orb_array)
                hash = (1099511628211_int64 * hash) + &
                        int(RandomHash2(orb_array(i))*i,int64)
            enddo
        endif
        hashInd = int(abs(mod(hash, int(HashIndexLength, int64))),sizeof_int)+1

    end function FindWalkerHash

    ! -- Hash table routines applicable to hash tables that reference lists of
    ! *any* type of data -----------------------------------------------------

    pure subroutine init_hash_table(hash_table)

        ! Take a just-allocated hash table, which must be empty, and
        ! initialise it by nullifying all pointers and setting all entries to
        ! zero.

        type(ll_node), pointer, intent(inout) :: hash_table(:)
        integer :: i

        do i = 1, size(hash_table)
            hash_table(i)%ind = 0
            nullify(hash_table(i)%next)
        end do

    end subroutine init_hash_table

    pure subroutine clear_hash_table(hash_table)

        ! Take hash_table and clear it. This is done by nullifying all pointers
        ! in all the linked lists that form the hash table, and setting the
        ! first index to zero.

        type(ll_node), pointer, intent(inout) :: hash_table(:)
        type(ll_node), pointer :: curr, prev
        integer :: i

        ! Loop over all entries corresponding to different hash values.
        do i = 1, size(hash_table)
            ! Point to the second entry in this linked list.
            curr => hash_table(i)%next
            ! Point to the first entry in this linked list.
            prev => hash_table(i)
            ! Set the first index to zero.
            prev%ind = 0
            nullify(prev%next)
            ! Loop over the whole linked list and deallocate all pointers.
            do while (associated(curr))
                prev => curr
                curr => curr%next
                deallocate(prev)
            end do
        end do

        nullify(curr)
        nullify(prev)

    end subroutine clear_hash_table

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
#ifdef DEBUG_
        character(len=*), parameter :: this_routine = "remove_hash_table_entry"
#endif

        ASSERT(all(nI <= nBasis))
        ASSERT(all(nI > 0))

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

    subroutine update_hash_table_ind(hash_table, nI, ind_old, ind_new)

        ! Find and remove the entry in hash_table corresponding to nI, which
        ! must have index ind in the hash table. If not found then an error
        ! will be thrown.

        type(ll_node), pointer, intent(inout) :: hash_table(:)
        integer, intent(in) :: nI(:)
        integer, intent(in) :: ind_old, ind_new

        integer :: hash_val
        type(ll_node), pointer :: prev, curr
        logical :: found
#ifdef DEBUG_
        character(len=*), parameter :: this_routine = "update_hash_table_ind"
#endif

        ASSERT(all(nI <= nBasis))
        ASSERT(all(nI > 0))

        found = .false.

        ! Find the hash value corresponding to this determinant.
        hash_val = FindWalkerHash(nI, size(hash_table))
        ! Point at the start of the linked list for this hash value.
        curr => hash_table(hash_val)
        prev => null()
        ! Loop over all entries in the linked list until we find the one equal
        ! to ind, the entry that we want to remove.
        do while (associated(curr))
            if (curr%ind == ind_old) then
                ! If this is the state to be removed.
                found = .true.
                curr%ind = ind_new
                exit
            end if
            prev => curr
            curr => curr%next
        end do

        ASSERT(found)

    end subroutine update_hash_table_ind

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

    pure subroutine add_hash_table_entry(hash_table, ind, hash_val)

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

    pure subroutine hash_table_lookup_int_32(orbs, targ, max_elem, hash_table, targ_array, ind, hash_val, found)

        ! Perform a search of targ_array for targ, by using hash_table.
        ! Only elements 0:max_elem will be used in comparing targ.
        ! If targ is found then found will be returned as .true. and ind will
        ! contain the index of targ in targ_array. Else, found will be
        ! returned as .false. and ind will be unset.

        integer, intent(in) :: orbs(:)
        integer(int32), intent(in) :: targ(0:)
        integer, intent(in) :: max_elem
        ! Note that hash_table won't actually change, but we need the inout
        ! label to make this routine pure.
        type(ll_node), pointer, intent(inout) :: hash_table(:)
        integer(int32), intent(in) :: targ_array(0:,:)
        integer, intent(out) :: ind
        integer, intent(out) :: hash_val
        logical, intent(out) :: found

        type(ll_node), pointer :: temp_node

        found = .false.

        hash_val = FindWalkerHash(orbs, size(hash_table))
        ! Point to the first node in the linked list for this hash value.
        temp_node => hash_table(hash_val)
        ! If temp_node%ind = 0 then there are no entries in targ_array
        ! with this hash, so the target is not in targ_array, so return.
        if (temp_node%ind /= 0) then
            do while (associated(temp_node))
                ! Check using the index stored in temp_node to see if we have
                ! found the searched-for determinant.
                if (all(targ(0:max_elem) == targ_array(0:max_elem, temp_node%ind))) then
                    found = .true.
                    ind = temp_node%ind
                    exit
                end if
                ! Move on to the next determinant with this hash value.
                temp_node => temp_node%next
            end do
        end if

        nullify(temp_node)

    end subroutine hash_table_lookup_int_32

    pure subroutine hash_table_lookup_int_64(orbs, targ, max_elem, hash_table, targ_array, ind, hash_val, found)

        ! Perform a search of targ_array for targ, by using hash_table.
        ! Only elements 0:max_elem will be used in comparing targ.
        ! If targ is found then found will be returned as .true. and ind will
        ! contain the index of targ in targ_array. Else, found will be
        ! returned as .false. and ind will be unset.

        integer, intent(in) :: orbs(:)
        integer(int64), intent(in) :: targ(0:)
        integer, intent(in) :: max_elem
        ! Note that hash_table won't actually change, but we need the inout
        ! label to make this routine pure.
        type(ll_node), pointer, intent(inout) :: hash_table(:)
        integer(int64), intent(in) :: targ_array(0:,:)
        integer, intent(out) :: ind
        integer, intent(out) :: hash_val
        logical, intent(out) :: found

        type(ll_node), pointer :: temp_node

        found = .false.

        hash_val = FindWalkerHash(orbs, size(hash_table))
        ! Point to the first node in the linked list for this hash value.
        temp_node => hash_table(hash_val)
        ! If temp_node%ind = 0 then there are no entries in targ_array
        ! with this hash, so the target is not in targ_array, so return.
        if (temp_node%ind /= 0) then
            do while (associated(temp_node))
                ! Check using the index stored in temp_node to see if we have
                ! found the searched-for determinant.
                if (all(targ(0:max_elem) == targ_array(0:max_elem, temp_node%ind))) then
                    found = .true.
                    ind = temp_node%ind
                    exit
                end if
                ! Move on to the next determinant with this hash value.
                temp_node => temp_node%next
            end do
        end if

        nullify(temp_node)

    end subroutine hash_table_lookup_int_64

    ! -- Hash table routines applicable to hash tables that reference lists of
    ! determinants which are stored in the standard ilut form ----------------

    subroutine fill_in_hash_table(hash_table, table_length, walker_list, list_length, ignore_unocc)

        ! This assumes that the input hash table is clear (use clear_hash_table)
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
#ifdef DEBUG_
        character(*), parameter :: this_routine = "fill_in_hash_table"
#endif

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
            ASSERT(all(nI <= nBasis))
            ASSERT(all(nI > 0))

            hash_val = FindWalkerHash(nI, table_length)

            call add_hash_table_entry(hash_table, i, hash_val)
        end do

    end subroutine fill_in_hash_table

    subroutine rm_unocc_dets_from_hash_table(hash_table, walker_list, list_length)

        ! This routine loops through all determinants in walker_list removes
        ! all entries from hash_table for determinants which are both
        ! unoccupied and not core determinants.

        type(ll_node), pointer, intent(inout) :: hash_table(:)
        integer(n_int), intent(in) :: walker_list(0:,:)
        integer, intent(in) :: list_length

        integer :: i, hash_val, nI(nel)
        real(dp) :: real_sign(lenof_sign)
        logical :: found, tCoreDet
        type(ll_node), pointer :: temp_node, prev
#ifdef DEBUG_
        character(len=*), parameter :: this_routine = "rm_unocc_dets_from_hash_table"
#endif

        do i = 1, list_length
            call extract_sign(walker_list(:,i), real_sign)
            tCoreDet = .false.
            if (tSemiStochastic) tCoreDet = test_flag(walker_list(:,i), flag_deterministic)
            if ((.not. IsUnoccDet(real_sign)) .or. tCoreDet) cycle
            found = .false.
            call decode_bit_det(nI, walker_list(:, i))

            ASSERT(all(nI <= nBasis))
            ASSERT(all(nI > 0))

            hash_val = FindWalkerHash(nI, size(hash_table))
            temp_node => hash_table(hash_val)
            prev => null()
            if (.not. temp_node%ind == 0) then
                ! Loop over all entries with this hash value.
                do while (associated(temp_node))
                    if (temp_node%ind == i) then
                        found = .true.
                        call remove_node(prev, temp_node)
                        exit
                    end if
                    ! Move on to the next determinant with this hash value.
                    prev => temp_node
                    temp_node => temp_node%next
                end do
            end if
            ASSERT(found)
        end do

    end subroutine rm_unocc_dets_from_hash_table

end module hash
