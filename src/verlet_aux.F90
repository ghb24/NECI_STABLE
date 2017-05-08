#include "macros.h"

! module for auxiliary routines used in the third order verlet algorithm

module verlet_aux
  use real_time_data, only: spawnBuf, spawnBufSize
  
  contains

    subroutine init_verlet_iteration()
      implicit none         

      FreeSlot(1:iEndFreeSlot) = 0
      iStartFreeSlot = 1
      iEndFreeSlot = 0

      ! verlet always uses a spawn hashtable
      call clear_hash_table(spawn_ht)      

    end subroutine init_verlet_iteration

    ! applies the hamiltonian twice to the current population and stores the result
    ! in spawnedParts
    subroutine obtain_h2_psi()
      implicit none

      ! apply H once, we now have the spawnedParts from a single iteration
      call apply_hamiltonian(CurrentDets,TotWalkers,.true.,tTruncInitiator)

      ! communicate the spawns between processors and store the compressed spawns into a buffer
      call generate_spawn_buf()

      ! apply H to the buffer to get H^2 on the original population. The result is stored
      ! in spawnedParts (uncommunicated)
      call apply_hamiltonian(spawnBuf,spawnBufSize,.false.,.false.)
      
      ! communicate the result and compress the population (such that each determinant
      ! only occurs once)
      call SendProcNewParts(spawnBufSize,.false.)
      call CompressSpawnedList(spawnBufSize,iter_data_fciqmc)
    end subroutine obtain_h2_psi
    
    subroutine apply_hamiltonian(population,popsize,tGetFreeSlots,tGetInitFlags)
      ! this subroutine performs spawning from population to spawnVec by applying
      ! delta t H once. No annihilation is performed and no other steps performed
      ! keep it minimalistic and stick to the SRP principle

      ! by default writes into SpawnedParts
      implicit none
      integer, intent(in) :: population(0:,1:), popsize
      ! store the free slots in population
      integer, intent(out) :: spawnVec(0:,1:)
      logical, intent(in) :: tGetInitflags, tGetFreeSlots
      integer :: idet, nI(nel), determ_index
      real(dp) :: parent_sign(lenof_sign), unused_diagH
      logical :: tEmptyDet, tCoreDet
      
      ! where to put this?
      attempt_create => attempt_create_normal

      ValidSpawnedList = InitialSpawnedSlots

      do idet = 1, popsize
         ! apply spawn and death for each walker
         fcimc_excit_gen_store%tFilled = .false.
         parent_flags = 0

         call extract_bit_rep(population(:,idet)), nI, parent_sign, unused_flags, &
              fcimc_excit_gen_store)

         tEmptyDet = IsUnoccDet(parent_sign)
         ! we do not collect freeSlots here, this is done before
         if(tEmptyDet) then
            if(tGetFreeSlots) then
               iEndFreeSlot = iEndFreeSlot + 1
               FreeSlot(iEndFreeSlot) = idet
            endif
            cycle
         endif

         if(tGetInitFlags) call CalcParentFlag(idet, unused_flags, unused_diagH)
         
         tCoreDet = check_determ_flat(population(:,idet))

         if(tCoreDet) then
            indices_of_determ_states(determ_index) = idet
            partial_determ_vecs(:,determ_index) = parent_sign
            denterm_index = determ_index + 1
            if(IsUnoccDet(parent_sign)) cycle
         endif

         ! the initiator flags are set upon the first iteration
         
         ! TODO: implement semi-stochastic approach
         call perform_spawn(idet,parent_sign,nI_parent,fcimc_excit_gen_store,tCoreDet)
    end subroutine apply_hamiltonian

    subroutine perform_spawn(nI,iLut_parent,det_sign)
      implicit none
      real(dp), intent(in) :: parent_sign(lenof_sign)
      integer, intent(in) :: nI(nel)
      integer(n_int), intent(in) :: ilut_parent(0:niftot)
      integer :: part, nspawn, ispawn, nI_child(nel), ic, ex(2,2), unused_ex_level
      integer(n_int) :: ilut_child(0:niftot)
      real(dp) :: prob, child_sign(lenof_sign), hdiag
      logical :: tParity
      HElement_t(dp) :: HElGen

      hdiag = get_helement(nI,nI,0)
      
      unused_ex_leve = 0
      do part = 1, lenof_sign
         call decide_num_to_spawn(parent_sign(part),AvMCExcits,nspawn)
         do ispawn = 1, nspawn
            ilut_child = 0_n_int
            
            call generate_excitation(nI,iLut_parent,nI_child,ilut_child,exFlag,ic,&
                 ex,tParity,prob,HElGen,fcimc_excit_gen_store)
            
            if(.not. IsNullDet(nI_child)) then
               child_sign = attempt_create(nI,iLut_parent,parent_sign,nI_child,iLut_child, &
                    prob, HElGen, ic, ex, tParity, unused_ex_level)
            else
               child_sign = 0.0_dp
            endif
            
            if ((any(child_sign /= 0)) .and. (ic /= 0) .and. (ic <= 2)) then               
               call create_particle_with_hash_table (nI_child, ilut_child, child_sign, &
                    part, population(:,idet), iter_data_fciqmc)
            end if ! If a child was spawned.

         end do ! Over mulitple particles on same determinant.
         
         diag_sign = -attempt_die_normal(nI,hdiag,parent_sign,unused_excit_level)
         if(any(abs(diag_sign)) > EPS) &
            call create_diagonal_as_spawn(nI_parent, iLut_parent, diag_sing, iter_data_fciqmc)
      end do

    end subroutine perform_spawn

    subroutine generate_spawn_buf()
      implicit none
      character(*), parameter :: this_routine = "generate_spawn_buf"
      
      call SendProcNewParts(spawnBufSize,.false.)
      call CompressSpawnedList(spawnBufSize,iter_data_fciqmc)
      spawnBuf(:,1:MaxIndex) = SpawnedParts(:,1:MaxIndex)
      spawnBufSize = MaxIndex
    end subroutine generate_spawn_buf

    subroutine merge_ilut_lists(listA, listB, hashTable, sizeA, sizeB, maxSizeA)
      implicit none
      integer(n_int), intent(in) :: listA(0:,1:), listB(0:,1:)
      integer, intent(in) :: sizeB, maxSizeA
      integer, intent(inout) :: sizeA
      type(ll_node), pointer, intent(inout) :: hashTable
      integer :: nI(nel), nJ(nel), hashIndex, ilutIndex, hashValue, i, j, run
      integer :: signA(lenof_sign), signB(lenof_sign), insertPos
      logical :: tSuccess
      
      ! this merges listB into listA. In the current form, empty slots in listA are 
      ! not exploited, because it should not make a difference currently (optimization follows)
      do i = 1, sizeB
         call decode_bit_det(nJ, listB(:,i))
         
         call hash_table_lookup(nJ, listB(:,i), nifdbo, hashTable, listA, ilutIndex, &
              hashValue, tSuccess)
         
         if(tSuccess) then 
            ! the i-th determinant in listB is already present in listA
            ! -> add up the signs
            call extract_sign(listA(:,ilutIndex), signA)
            call extract_sign(listB(:,i), signB)
            
            ! we do not fill up empty slots, so we do not care if signA == 0
            ! this differs from AnnihilateSpawnedParts
            call encode_sign(listA(:,ilutIndex),signA+signB)
            call encode_sign(listB(:,i),null_part)            

            ! the initiator criterium is checked upon annihilation, no need to do so here
         else
            ! if the entry in listB is not in listA, add it as the last entry
            sizeA = sizeA + 1
            insertPos = sizeA
            if( sizeA > maxSizeA) &
               call stop_all(this_routine, "Out of memory for merging ilut lists")
            listA(:,insertPos) = listB(:,i)
            call add_hash_table_entry(hashTable,insertPos,hashValue)
         end if
    end subroutine merge_ilut_lists

    subroutine update_delta_psi()
      ! here, we add H^2 psi to delta_psi to generate the new delta_psi
      ! this might cause trouble as stochastic error is not averaged out, 
      ! but we just carry over the error from the last iteration
      implicit none
      
      ! add up delta_psi from the last iteration and spawnedParts (i.e. H^2 psi)
      call merge_ilut_list(spawnedParts, dpsi_cache, spawn_ht, spawnBufSize, &
           dpsi_size, maxSpawned)
      ! cache delta_psi for the next iteration
      dpsi_cache(:,spawnBufSize) = spawnedParts(:,spawnBufSize)
      dpsi_size = spawnBufSize
      ! for now, consider all entries in delta_psi as safe spawns for the next iteration
      call set_initiator_flags_array(dpsi_cache,dpsi_size)
    end subroutine update_delta_psi

    subroutine set_initator_flags_array(list, listSize)
      implicit none
      integer(n_int), intent(inout) :: list(0:,1:)
      integer, intent(in) :: listSize
      integer :: i, run
      
      do i = 1, listSize
         do run = 1, inum_runs
            call set_flag(list(:,i),get_initiator_flag_by_run(run))
         end do
      end do
    end subroutine set_initator_flags
end module verlet_aux
