#include "macros.h"

module real_time_aux
  use real_time_data, only: overlap_states, gf_count
  use bit_rep_data, only: extract_sign, niftot, nifdbo
  use bit_reps, only: decode_bit_det, flag_deterministic
  use FciMCData, only: SpawnedParts, ValidSpawnedList, MaxSpawned, InitialSpawnedSlots, &
       core_space, determ_sizes, determ_space_size, CurrentDets, hashindex
  use SystemData, only: nel
  use constants, only: dp, lenof_sign, n_int, iout
  use ParallelHelper, only: iProcIndex, nNodes, MPIBarrier
  use Parallel_neci

  contains

    subroutine write_overlap_state_serial(state,length,index)
      implicit none
      integer(n_int), intent(in) :: state(0:nIfTot,length)
      integer, intent(in) :: index, length
      integer :: iProc, ierr
      do iProc = 0, nProcessors-1
         ! sequentialize to overcome bottlenecks for shared memory systems (this is not performance critical)
         if(iProc .eq. iProcIndex) then
            call write_overlap_state(state,length,index)
         endif
         call MPIBarrier(ierr)
      enddo

    end subroutine write_overlap_state_serial

    subroutine write_overlap_state(state, length, index)
      implicit none
      integer(n_int), intent(in) :: state(0:nIfTot,length)
      integer, intent(in) :: index, length
      integer :: nOccDets, i, ierr
      integer(n_int), allocatable :: buffer(:,:)
      real(dp) :: tmp_sign(lenof_sign)
      character(*), parameter :: this_routine = "write_overlap_state"

      ! check how many determinants are stored for this state on this core
      allocate(buffer(0:niftot,length),stat=ierr)
      nOccDets = 0
      do i=1, length
         call extract_sign(state(:,i), tmp_sign)
         if(IsUnoccDet(tmp_sign)) then
            cycle
         endif
         nOccDets = nOccDets + 1
         buffer(:,nOccDets) = state(:,i)
      enddo

      ! copy them to overlap_states
      if(allocated(overlap_states(index)%dets)) deallocate(overlap_states(index)%dets)
      allocate(overlap_states(index)%dets(0:nIfTot,nOccDets))
      overlap_states(index)%dets(:,1:nOccDets) = buffer(:,1:nOccDets)
      overlap_states(index)%nDets = nOccDets
      deallocate(buffer)
    end subroutine write_overlap_state

!------------------------------------------------------------------------------------------!

    ! for load balancing, the load balancing blocks of the overlap states have to be
    ! communicated alongside the currentdets.
    subroutine move_overlap_block(block,tgt_proc)
      use load_balance_calcnodes, only: LoadBalanceMapping, get_det_block
      use bit_reps, only: nullify_ilut
      implicit none
      integer, intent(in) :: block, tgt_proc
      integer(n_int), allocatable :: perturbed_buf(:,:)
      integer :: tmp_totwalkers, iGf, ierr, src_proc, nsend, idet, &
           det_block, offset, nelem
      integer :: det(nel)
      real(dp) :: sgn(lenof_sign)
      integer, parameter :: mpi_tag_nsend = 223458
      integer, parameter :: mpi_tag_dets = 223459

      if(allocated(overlap_states)) then
         call MPIBarrier(ierr)

         if(iProcIndex == root) print *, "Moving overlap states"

         src_proc = LoadBalanceMapping(block)
         ! send the corresponding block
         do iGf = 1, gf_count
            if (iProcIndex == src_proc) then
               nsend = 0

               do idet = 1, overlap_states(iGf)%nDets
                  call extract_sign(overlap_states(iGf)%dets(:,idet), sgn)
                  if(IsUnoccDet(sgn)) cycle
                  call decode_bit_det(det, overlap_states(iGf)%dets(:,idet))
                  det_block = get_det_block(nel,det,0)
                  if(det_block == block) then
                     nsend = nsend + 1
                     SpawnedParts(:,nsend) = overlap_states(iGf)%dets(:,idet)
                     call nullify_ilut(overlap_states(iGf)%dets(:,idet))
                  endif
               end do
               nelem = nsend * (1 + niftot)

               call MPISend(nsend,1,tgt_proc,mpi_tag_nsend,ierr)
               call MPISend(SpawnedParts(:,1:nsend),nelem,tgt_proc,mpi_tag_dets,ierr)

               ! compress the corresponding overlap state
               call write_overlap_state(overlap_states(iGf)%dets,overlap_states(iGf)%nDets,iGf)
            else if (iProcIndex == tgt_proc) then
               ! here we have to expand the overlap state array, so we need a buffer

               call MPIRecv(nsend,1,src_proc,mpi_tag_nsend,ierr)
               nelem = nsend * (1 + niftot)
               call MPIRecv(SpawnedParts,nelem,src_proc,mpi_tag_dets,ierr)

               offset = overlap_states(iGf)%nDets
               tmp_totwalkers = offset + nsend
               allocate(perturbed_buf(0:niftot,tmp_totwalkers),stat=ierr)
               perturbed_buf(:,1:offset) = overlap_states(iGf)%dets(:,1:offset)

               do idet = 1, nsend
                  perturbed_buf(:,offset+idet) = SpawnedParts(:,idet)
               enddo

               deallocate(overlap_states(iGf)%dets)
               allocate(overlap_states(iGf)%dets(0:niftot,tmp_totwalkers),stat=ierr)
               if(ierr /= 0) print *, "WARNING: Allocation failed in overlap_states load balancing"
               overlap_states(iGf)%dets(:,1:tmp_totwalkers) = perturbed_buf(:,1:tmp_totwalkers)
               deallocate(perturbed_buf)
               overlap_states(iGf)%nDets = tmp_totwalkers
            endif
            call MPIBarrier(ierr)
         enddo
      endif

    end subroutine move_overlap_block

!------------------------------------------------------------------------------------------!

    subroutine add_semistochastic_state(ilut_list, list_size, ssht, ilut)
      use FciMCData, only: ll_node, MaxWalkersPart
      use hash, only: hash_table_lookup, add_hash_table_entry
      use bit_reps, only: encode_sign
      implicit none
      type(ll_node), pointer, intent(inout) :: ssht(:)
      integer(n_int), intent(inout) :: ilut_list(0:,1:)
      integer, intent(inout) :: list_size
      integer(n_int), intent(in) :: ilut(0:nifTot)
      integer :: index, hash_val, nI(nel)
      logical :: tSuccess
      real(dp) :: new_sgn(lenof_sign), old_sgn(lenof_sign)
      character(*), parameter :: this_routine = "add_semistochastic_state"

      call decode_bit_det(nI,ilut)
      call hash_table_lookup(nI,ilut,nifDBO,ssht,ilut_list,index,hash_val,tSuccess)
      ! If it is already in corespace, we check which population is higher
      if(tSuccess) then
         call extract_sign(ilut_list(:,index),old_sgn)
         call extract_sign(ilut,new_sgn)
         if(sum(abs(new_sgn)) > sum(abs(old_sgn))) call encode_sign(ilut_list(:,index),new_sgn)
      else
         if(list_size > MaxWalkersPart) call stop_all(this_routine, &
              "Out of memory for corespace construction.")
         list_size = list_size + 1
         ilut_list(:,list_size) = ilut
         call add_hash_table_entry(ssht,list_size,hash_val)
      endif
    end subroutine add_semistochastic_state

!------------------------------------------------------------------------------------------!

end module real_time_aux
