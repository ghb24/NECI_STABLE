#include "macros.h"
module trial_ht_procs

  use FciMCData, only: NConEntry, CurrentDets, TotWalkers, con_send_buf
  use load_balance_calcnodes
  use searching, only: hash_search_trial
  use sparse_arrays, only: trial_hashtable
  use hash
  use constants

  implicit none

contains

!------------------------------------------------------------------------------------------!

    function buffer_trial_ht_entries(block, source_ht, source_ht_size) result(nsend)
      integer, intent(in) :: block
      type(trial_hashtable), intent(inout) :: source_ht(:)
      integer, intent(in) :: source_ht_size
      integer :: nsend
      integer :: clashes, j, k, det_block, det(nel)
      integer(n_int) :: source_state(0:NConEntry)
      ! get all entries from source_ht that belong to block and move them
      ! to con_send_buf, deleting them from source_ht in the process
      nsend = 0
      do j = 1, source_ht_size
         clashes = source_ht(j)%nclash
         if(clashes > 0) then
            k = 0
            do
               k = k + 1
               call decode_bit_det(det,source_ht(j)%states(:,k))
               det_block = get_det_block(nel, det, 0)
               if(det_block == block) then
                  call extract_trial_ht_entry(j,k,source_state,source_ht)
                  nsend = nsend + 1
                  con_send_buf(:,nsend) = source_state
                  clashes = clashes - 1
                  k = k - 1
               endif
               if(k==clashes) exit
            enddo
         endif
      end do
    end function buffer_trial_ht_entries

!------------------------------------------------------------------------------------------!

    subroutine extract_trial_ht_entry(hash_val, i, ht_entry,source_ht)
      implicit none
      integer(n_int), intent(out) :: ht_entry(0:NConEntry)
      integer, intent(in) :: hash_val, i
      type(trial_hashtable), intent(inout) :: source_ht(:)
      integer :: clashes
      character(*), parameter :: this_routine = "extract_trial_ht_entry"

      ! get the stores state
      ht_entry = source_ht(hash_val)%states(:,i)
      ! then remove it from the table
      clashes = source_ht(hash_val)%nclash
      call remove_trial_ht_entry(hash_val,i,clashes,source_ht)
    end subroutine extract_trial_ht_entry

!------------------------------------------------------------------------------------------!

    subroutine remove_trial_ht_entry(hash_val, index, clashes, source_ht)
      implicit none
      integer, intent(in) :: hash_val, index, clashes
      type(trial_hashtable), intent(inout) :: source_ht(:)
      integer(n_int), allocatable :: tmp(:,:)
      integer :: i, ierr
      character(*), parameter :: this_routine = "remove_trial_ht_entry"

      ! first, copy the contnet of the source_ht entry to a temporary
      ! if there is any to be left
      if(clashes-1 > 0) then
         allocate(tmp(0:NConEntry,clashes-1), stat = ierr)
         if(ierr .ne. 0) call stop_all(this_routine, "Failed allocation")
         do i = 1, index - 1
            tmp(:,i) = source_ht(hash_val)%states(:,i)
         end do
         ! omitting the element to remove
         do i = index + 1, clashes
            tmp(:,i-1) = source_ht(hash_val)%states(:,i)
         end do

         ! then, reallocate the source_ht entry (if required)
         deallocate(source_ht(hash_val)%states, stat = ierr)
         if(ierr .ne. 0) call stop_all(this_routine, "Failed deallocation")
         allocate(source_ht(hash_val)%states(0:NConEntry,clashes-1), stat = ierr)
         if(ierr .ne. 0) call stop_all(this_routine, "Failed allocation")
         ! and copy the temporary back (if it is non-empty)
         source_ht(hash_val)%states(0:NConEntry,:) = tmp(0:NConEntry,:)
         deallocate(tmp,stat = ierr)
         if(ierr .ne. 0) call stop_all(this_routine, "Failed deallocation")
      else
         ! just to be sure, allocate with size 0
         deallocate(source_ht(hash_val)%states)
         allocate(source_ht(hash_val)%states(0:NConEntry,0))
      endif

      ! finally, update the nclashes information
      source_ht(hash_val)%nclash = clashes - 1
    end subroutine remove_trial_ht_entry

!------------------------------------------------------------------------------------------!

    subroutine add_trial_ht_entries(entries, n_entries, source_ht, source_ht_size)
      ! this adds n_entries entries to the source_ht hashtable
      implicit none
      integer, intent(in) :: n_entries
      integer(n_int), intent(in) :: entries(0:NConEntry,n_entries)
      type(trial_hashtable), allocatable, intent(inout) :: source_ht(:)
      ! in principle, the source_ht can be resized if really required,
      ! thus changing source_ht_size
      integer, intent(inout) :: source_ht_size
      integer :: i, hash_val, nI(nel), clashes
      ! we need to be careful: if the source_ht happens to be empty,
      ! it needs to be resized, as lookups are not possible on
      ! empty hashtables and will throw an error
      if(source_ht_size == 0) then
         call resize_trial_ht(source_ht,source_ht_size,1)
      endif

      do i = 1, n_entries
         call decode_bit_det(nI,entries(:,i))
         hash_val = FindWalkerHash(nI, source_ht_size)
         ! just add them one by one
         call add_single_trial_ht_entry(entries(:,i),hash_val, source_ht)
      enddo
    end subroutine add_trial_ht_entries

!------------------------------------------------------------------------------------------!

    subroutine add_single_trial_ht_entry(ht_entry, hash_val, source_ht)
      implicit none
      integer(n_int), intent(in) :: ht_entry(0:NConEntry)
      integer, intent(in) :: hash_val
      type(trial_hashtable), intent(inout) :: source_ht(:)
      integer :: clashes, ntrial ,ncon
      integer(n_int), allocatable :: tmp(:,:)

      ! add a single entry to trial_ht with hash_val
      clashes = source_ht(hash_val)%nclash
      ! store the current entries in a temporary
      allocate(tmp(0:NConEntry,clashes+1))
      ! if there are any, copy them now
      if(allocated(source_ht(hash_val)%states)) then
         tmp(:,:clashes) = source_ht(hash_val)%states(:,:)
         ! then deallocate
         deallocate(source_ht(hash_val)%states)
      endif
      ! add the new entry
      tmp(:,clashes+1) = ht_entry

      ! and allocoate the new entry
      allocate(source_ht(hash_val)%states(0:NConEntry,clashes+1))
      ! fill it
      source_ht(hash_val)%states = tmp
      deallocate(tmp)
      ! and update the nclashes info
      source_ht(hash_val)%nclash = clashes + 1

    end subroutine add_single_trial_ht_entry

!------------------------------------------------------------------------------------------!

    subroutine resize_trial_ht(source_ht, source_ht_size, new_size)
      ! take a trial_hashtable of size source_ht_size and resize it to new_size
      implicit none
      type(trial_hashtable), allocatable, intent(inout) :: source_ht(:)
      integer, intent(inout) :: source_ht_size
      integer, intent(in) :: new_size

      type(trial_hashtable), allocatable :: tmp(:)
      integer :: ierr, i

      ! no support for shrinking down hashtables
      if(source_ht_size > new_size) return

      if(allocated(source_ht)) then
         ! if the source_ht contains data, save it
         if(source_ht_size > 0) then
            allocate(tmp(source_ht_size), stat = ierr)
            tmp = source_ht
         endif
         ! if source_ht is already there, deallocate it
         deallocate(source_ht)
      endif

      ! now, create the new source_ht
      allocate(source_ht(new_size), stat = ierr)
      ! if we saved data, move it now
      if(allocated(tmp)) then
         source_ht(1:source_ht_size) = tmp(1:source_ht_size)
         deallocate(tmp)
      endif

      ! then, fill up the newly created entries with empty lists
      do i = source_ht_size+1, new_size
         source_ht(i)%nclash = 0
      end do

      ! finally, declare that source_ht is now of size new_size
      source_ht_size = new_size
    end subroutine resize_trial_ht

!------------------------------------------------------------------------------------------!

    subroutine count_trial()
      use Parallel_neci, only: MPISumAll
      implicit none
      integer ::  ntrialtot, ncontot, ntrial, ncon

      call count_trial_this_proc(ntrial, ncon)
      call MPISumAll(ntrial,ntrialtot)
      call MPISumAll(ncon,ncontot)
      write(iout,*) "Trial states ", ntrialtot
      write(iout,*) "Connected states ", ncontot
    end subroutine count_trial

!------------------------------------------------------------------------------------------!

    subroutine count_trial_this_proc(ntrial, ncon)
      use searching, only: hash_search_trial
      use FciMCData, only: ntrial_excits
      implicit none
      integer, intent(out) :: ntrial, ncon
!       integer :: i, nI(nel)
      integer(int64) :: i
      integer :: nI(nel)
      real(dp) :: sgn(lenof_sign)
      logical :: tTrial, tCon
      HElement_t(dp) :: amp(ntrial_excits)

      ntrial = 0
      ncon = 0
      do i = 1, TotWalkers
         call decode_bit_det(nI, CurrentDets(:,i))
         call extract_sign(CurrentDets(:,i),sgn)
         if(IsUnoccDet(sgn)) cycle
         call hash_search_trial(CurrentDets(:,i),nI,amp,tTrial,tCon)
         if(tTrial) ntrial = ntrial + 1
         if(tCon) ncon = ncon + 1
      end do

    end subroutine count_trial_this_proc

end module trial_ht_procs
