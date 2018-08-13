#include "macros.h"

module sign_coherence_check
  use SystemData, only: tHPHF
  use hash, only: clear_hash_table, add_hash_table_entry, hash_table_lookup, &
       remove_hash_table_entry, FindWalkerHash
  use FciMCData, only: clistIndex, AllConflictingDets, ConflictingDets, &
       tStoreConflicts, maxConflictExLvl, avSigns, fcimc_iter_data, conflictHash, &
       conflictExLvl, nCDetsStore, max_calc_ex_level, HashIndex, CurrentDets, &
       cAccIter, cHashSize, correctionInterval, TotParts, sign_correction_time, &
       tSinglePartPhase
  use CalcData, only: tAvReps, DiagSft
  use SystemData, only: nel
  use global_det_data, only: det_diagH
  use constants
  use Parallel_neci
  use util_mod
  use global_det_data, only: increase_conflict_counter, get_conflict_counter
  use DetBitOps, only: FindBitExcitLevel
  use Determinants, only: get_helement
  use bit_reps, only: NIfTot, encode_sign, extract_sign, decode_bit_det, NIfDBO
  use hphf_integrals, only: hphf_off_diag_helement

contains

!------------------------------------------------------------------------------------------!

    subroutine communicate_conflict_buffers(ilut_list, all_ilut_list,this_size,tot_size)
      ! send the full list of sign conflicting dets to all other procs
      implicit none
      integer(n_int), intent(in) :: ilut_list(:,:)
      integer(n_int), intent(inout) :: all_ilut_list(:,:)
      integer, intent(in) :: this_size
      integer, intent(out) :: tot_size
      integer(MPIArg), dimension(nProcessors) :: sendcounts, disps, recvcounts, recvdisps
      integer(MPIArg) :: entrySize
      integer :: i, error
      real(dp) :: zero_sign(lenof_sign)

      entrySize = int(size(ilut_list,1),MPIArg)

      call set_timer(sign_correction_time)

      ! first, get the number to send
      sendcounts = this_size   ! we just send everything, the number of elements is clistIndex
      call MPIAllToAll(sendcounts,1,recvcounts,1,error)

      ! get the displacements
      recvdisps(1) = 0
      do i = 2, nProcessors
         recvdisps(i) = recvcounts(i-1) + recvdisps(i-1)
      end do   
      
      disps = 0
      ! get the total number of entries
      tot_size = (recvdisps(nProcessors)+recvcounts(nProcessors))
      ! scale with respect to entrySize
      do i = 1, nProcessors
         recvdisps(i) = recvdisps(i)*entrySize
         recvcounts(i) = recvcounts(i)*entrySize
         sendcounts(i) = sendcounts(i)*entrySize
         disps(i) = disps(i)*entrySize
      end do

      all_ilut_list = 0
      ! now send the stuff
      call MPIAllToAllv(ilut_list,sendcounts,disps,all_ilut_list,recvcounts,&
           recvdisps,error)

      ! and now: set all the signs stored in AllConflictingDets to 0 so we can
      ! accumulate corrections
      zero_sign = 0.0_dp
      do i = 1, tot_size
         call encode_sign(AllConflictingDets(:,i),zero_sign)
      end do

      ! omit the barrier and know that the time indicated is the time
      ! taken on some random proc (for efficiency estimation only)
      call halt_timer(sign_correction_time)

    end subroutine communicate_conflict_buffers

!------------------------------------------------------------------------------------------!

    subroutine apply_conflict_correction(iter_data, ilut_list, all_ilut_list, this_list_size)
      implicit none
      type(fcimc_iter_data), intent(inout) :: iter_data
      integer(n_int), intent(inout) :: ilut_list(:,:)
      integer(n_int), intent(in) :: all_ilut_list(:,:)
      integer, intent(in) :: this_list_size

      integer :: i, nI(nel), PartInd, DetHash, part
      real(dp) :: correctedSign(lenof_sign), oldSign(lenof_sign)
      logical :: tSuccess
      real(sp), save :: usedTime = 0.0_sp

      call set_timer(sign_correction_time)

      ! first, gather the correction contribution from all processors
      call gather_conflict_correction(ilut_list, all_ilut_list, int(this_list_size,MPIArg))
      ! then, loop over all determinants in ilut_list
      do i = 1, this_list_size
         call decode_bit_det(nI,ilut_list(:,i))
         
         call hash_table_lookup(nI,ilut_list(:,i),NIfDBO,HashIndex,CurrentDets,PartInd,&
              DetHash,tSuccess)
         if(tSuccess) then
            ! correct the sign
            call extract_sign(ilut_list(:,i),correctedSign)
            ! only for logging: old sign
            call extract_sign(CurrentDets(:,PartInd),oldSign)
            ! divide by Hii - S
            do part = 1, lenof_sign
               correctedSign(part)  =  &
                    - correctedSign(part) / &
                    ( cAccIter * (det_diagH(PartInd) - DiagSft(part)))
               ! additional normalziation 1/cAccIter to average over a number of
               ! iterations
               iter_data%nremoved(part) = iter_data%nremoved(part) &
                    - abs(correctedSign(part)) &
                    + abs(oldSign(part))
            end do
            ! and update the sign
            TotParts = TotParts - abs(oldSign) + abs(correctedSign)
            call encode_sign(CurrentDets(:,PartInd),correctedSign)

         endif
      end do      

      call halt_timer(sign_correction_time)
      write(iout,*) "Time (seconds) taken for sign correction across replicas", &
           get_total_time(sign_correction_time) - usedTime
      usedTime = get_total_time(sign_correction_time) 
    end subroutine apply_conflict_correction

!------------------------------------------------------------------------------------------!

    subroutine gather_conflict_correction(ilut_list, all_ilut_list, this_list_size)
      implicit none
      integer(n_int), intent(inout) :: ilut_list(:,:)
      integer(n_int), intent(in) :: all_ilut_list(:,:)
      integer(MPIArg), intent(in) :: this_list_size
      
      integer(MPIArg) :: list_sizes(nProcessors)
      integer :: i, j, error, disp(nProcessors), maxListSize
      real(dp), allocatable :: accumulatedSign(:,:)
      real(dp), allocatable :: sendSign(:,:)
      character(*), parameter :: t_r = "gather_conflict_correction"

      call set_timer(sign_correction_time)

      call MPIAllGather(this_list_size,list_sizes,error)

      maxListSize = maxVal(list_sizes)
      allocate(accumulatedSign(lenof_sign,maxListSize))
      allocate(sendSign(lenof_sign,maxListSize))

      ! get the positions where this procs sign for proc i is stored
      disp(1) = 0
      do i = 2, nProcessors
         disp(i) = disp(i-1) + list_sizes(i-1)
      end do

      do i = 1, nProcessors
         do j = 1, list_sizes(i)
            call extract_sign(all_ilut_list(:,j+disp(i)),sendSign(:,j))
         end do

         ! reduce the sign (use MPI_Reduce to be able to set the root proc)
         call MPI_Reduce(sendSign,accumulatedSign,lenof_sign*list_sizes(i),&
              MPI_DOUBLE_PRECISION,MPI_SUM,i-1,CommGlobal,error)

         if(error .ne. 0) call stop_all(t_r, "Error in MPI_Reduce")
      end do

      ! unpack the sign
      do i = 1, this_list_size
         call encode_sign(ilut_list(:,i),accumulatedSign(:,i))
      end do

      deallocate(accumulatedSign)
      deallocate(sendSign)

      call halt_timer(sign_correction_time)
      
    end subroutine gather_conflict_correction

!------------------------------------------------------------------------------------------!

    subroutine accumulate_conflict_correction(nI,ilut,sgn)
      use FciMCData, only: AllCListSize ! total number of conflicting dets
      implicit none
      integer, intent(in) :: nI(nel)
      integer(n_int), intent(in) :: ilut(0:NIfTot)
      real(dp), intent(in) :: sgn(lenof_sign)

      integer :: i, exLvl, cDet(nel)
      real(dp) :: accumulatedSign(lenof_sign)
      integer :: PartInd, DetHash
      logical :: tSuccess, tParity
      HElement_t(dp), parameter :: helZero = 0.0_dp
      integer :: exMat(2,2)
      
      call set_timer(sign_correction_time)

      do i = 1, AllCListSize
         ! check if it is a double/single
         exLvl = FindBitExcitLevel(ilut,AllConflictingDets(:,i),max_calc_ex_level)
         if(exLvl == 1 .or. exLvl == 2) then
            !if yes, get the matrix element


            call decode_bit_det(cDet,AllConflictingDets(:,i))
            ! get the stuff requird for computing the matrix element
            exMat = 0
            exMat(1,1) = 2
            call GetExcitation(nI,cDet,nel,exMat,tParity)

            ! and add it to the accumulated sign correction
            call extract_sign(AllConflictingDets(:,i),accumulatedSign)
            if(tHPHF) then 
               accumulatedSign = accumulatedSign + sgn*hphf_off_diag_helement(&
                    nI,cDet,ilut,AllConflictingDets(:,i),exLvl,exMat,tParity,helZero)
            else
               accumulatedSign = accumulatedSign + sgn*get_helement(nI,cDet,exLvl,&
                    ilut,AllConflictingDets(:,i))
            endif
            ! store the updated sign in the AllConflictingDets
            call encode_sign(AllConflictingDets(:,i),accumulatedSign)
         endif
      end do

      call halt_timer(sign_correction_time)
            
    end subroutine accumulate_conflict_correction

!------------------------------------------------------------------------------------------!

    subroutine reset_conflict_storage()
      implicit none
      ! reset the stored sign conflicts and clear the cache
      clistIndex = 0
      call clear_hash_table(conflictHash)
      AllConflictingDets = 0
      ConflictingDets = 0

    end subroutine reset_conflict_storage

!------------------------------------------------------------------------------------------!

    subroutine replica_coherence_check(iter,iter_data,ilut,j,sgn,exLvl)

      ! do a check if a determinant has sign consistent walkers across replicas
      ! input: 
      ! iter = iteration
      ! iter_data = fcimc_iter_data object for storing population changes
      ! ilut = determinant + population in ilut-format
      ! j = global index of ilut in CurrentDets
      ! sgn = sign of ilut (i.e. #walkers)
      ! exLvl = excitation level of ilut for logging purposes
      ! depending on the chosen options, the sign is adjusted
      ! in any case, statistics are gahtered on sign conflicts
      implicit none
      integer, intent(in) :: iter,j ! iter = iteration, j = global index of ilut
      type(fcimc_iter_data), intent(inout) :: iter_data
      integer(n_int), intent(inout) :: ilut(0:NIfTot)
      real(dp), intent(inout) :: sgn(lenof_sign)
      integer, intent(in) :: exLvl
      
      integer :: run
      real(dp) :: avSign

#ifdef __CMPLX
#else
      ! check if there are any sign changes within sgn
      if(any(sgn*sgn(1) < 0)) then
         avSign = sum(sgn)/inum_runs
         ! log the change in population
         do run = 1, inum_runs
            ! are we looking at an erroneous sign?
            if(sgn(run)*avSign < 0) then
               ! log the conflict
               avSigns = avSigns + abs(sgn(run))
               if(exLvl < maxConflictExLvl) ConflictExLvl(exLvl) = ConflictExLvl(exLvl) + 1

               if(tStoreConflicts .and. all(.not. tSinglePartPhase)) then
                  call store_conflict(iter,ilut,j,sgn)
               endif

               if(tAvReps) then
                  ! log the sign change
                  iter_data%nremoved(run) = iter_data%nremoved(run) &
                       + (abs(sgn(run)))

                  ! set the sign on the conflicting run to 0
                  sgn(run) = 0.0
               endif
            endif
         end do
         call encode_sign(ilut,sgn)
      endif
#endif
    end subroutine replica_coherence_check

!------------------------------------------------------------------------------------------!

    subroutine store_conflict(iter,ilut,j,sgn)
      implicit none
      integer, intent(in) :: iter,j
      integer(n_int), intent(in) :: ilut(0:NifTot)
      real(dp), intent(in) :: sgn(lenof_sign)
      integer :: nI(nel)
      integer :: PartInd, DetHash, i
      integer :: conflictPersistence, thisPersistence, iMin
      real(dp) :: av_sgn(lenof_sign), r
      logical :: tSuccess

      call decode_bit_det(nI,ilut)
      call hash_table_lookup(nI,ilut,NIfDBO, conflictHash, conflictingDets, PartInd, &
           DetHash, tSuccess)

      ! increase the conflict counter in global_det_data
      call increase_conflict_counter(j)

      ! if there already was a conflict here, do three things:
      ! 1. update the sign of the stored determinant
      ! 2. increase the conflict-counter by 1
      ! 3. set the iteration with the latest conflict to this one
      if(tSuccess) then
         ! update the sign
         call extract_sign(conflictingDets(:,PartInd),av_sgn)
         av_sgn = av_sgn + sgn

         ! increase the conflict-counter by 1
         conflictingDets(NIfTot,PartInd) = conflictingDets(NIfTot,PartInd) + 1
         
         ! set the iteration
         conflictingDets(NIfTot+1,PartInd) = iter
      else
         ! if the conflict is new, first check if the list is already full
         if(clistIndex == nCDetsStore) then
            ! if yes, throw out the least important one
            iMin = 1
            conflictPersistence = persistence(conflictingDets(:,1))
            do i = 2, nCDetsStore
               thisPersistence = persistence(conflictingDets(:,i))
               if(conflictPersistence > thisPersistence) then
                  conflictPersistence = thisPersistence
                  iMin = i
               endif
            end do
            ! if there is a conflict which is less important than the current
            ! one, replace that one
            if(conflictPersistence <= get_conflict_counter(j)) then
            ! remvoe the outdated hash-table entry
               call decode_bit_det(nI, conflictingDets(0:NifTot,iMin))
               call remove_hash_table_entry(conflictHash,nI,iMin)
            
            ! assign the new determinant in the conflictingDets list
               call assign_conflicting_det(iter,ilut,iMin,j,DetHash)
            endif


         else
            ! else, add a new entry
            clistIndex = clistIndex + 1
            call assign_conflicting_det(iter,ilut,clistIndex,j,DetHash)
         endif
      endif
      
    end subroutine store_conflict

!------------------------------------------------------------------------------------------!

    subroutine assign_conflicting_det(iter,ilut,index,j,DetHash)
      implicit none

      integer, intent(in) :: iter, index,j, DetHash
      integer(n_int), intent(in) :: ilut(0:NifTot)

      conflictingDets(0:NifTot,index) = ilut(0:NIfTot)
      conflictingDets(NifTot,index) = get_conflict_counter(j)
      conflictingDets(NifTot+1,index) = iter

      call add_hash_table_entry(conflictHash, index, DetHash)
    end subroutine assign_conflicting_det

!------------------------------------------------------------------------------------------!
    
    function persistence(cilut) result(p)
      implicit none
      integer(n_int), intent(in) :: cilut(0:NIfTot+1)
      integer :: p

      ! a measure of how persistent a conflict is if it occurred nOccur times 
      ! and the last occurrence was at iteration iter
      ! (the larger the more persistent)
      p = cilut(NIfTot)
    end function persistence

end module sign_coherence_check
