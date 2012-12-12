#include "macros.h"
!This module is to be used for various types of walker MC annihilation in serial and parallel.
MODULE AnnihilationMod
    use SystemData , only : NEl, tHPHF, nBasis, tCSF, tSemiStochastic
    use CalcData , only : TRegenExcitgens,tRegenDiagHEls, tEnhanceRemainder
    USE DetCalcData , only : Det,FCIDetIndex
    USE Parallel_neci
    USE dSFMT_interface , only : genrand_real2_dSFMT
    USE FciMCData
    use DetBitOps, only: DetBitEQ, DetBitLT, FindBitExcitLevel, ilut_lt, &
                         ilut_gt, DetBitZero
    use spatial_initiator, only: add_initiator_list, rm_initiator_list, &
                                 is_spatial_init
    use CalcData , only : tTruncInitiator, tSpawnSpatialInit, tSortDetermToTop
    use Determinants, only: get_helement
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use sort_mod
    use constants, only: n_int,lenof_sign,null_part,sizeof_int
    use bit_rep_data
    use bit_reps, only: decode_bit_det, extract_sign, extract_flags, &
                        encode_sign, encode_flags, test_flag, set_flag, &
                        clr_flag, flag_parent_initiator, encode_part_sign, &
                        extract_part_sign, copy_flag, nullify_ilut, &
                        nullify_ilut_part
    use csf_data, only: csf_orbital_mask
    use hist_data, only: tHistSpawn, HistMinInd2
    use Logging , only : tHF_Ref_Explicit
    IMPLICIT NONE

    contains

    !TODO:
    !   H Elements - send through logical to decide whether to create or not.
    !   Parallel spawned parts - create the ValidSpawnedList itself.
    !   Going to have to sort this out for the new packaged walkers - will have to package them up in this interface.
    subroutine AnnihilationInterface (TotDets, MainParts, MaxMainInd, &
                                      SpawnDets, SpawnParts, MaxSpawnInd, &
                                      iter_data)
        use constants, only: size_n_int
        use shared_alloc, only: shared_allocate_iluts, shared_deallocate
!This is an interface routine to the Direct Annihilation routines.
!It is not quite as fast as the main annihilation routines since there is a small degree of initialisation required
!which can be achieved on-the-fly if increased performance is required.

!       MainParts(:,:)      This is the main list of particles as determinants. It must be ordered, sign-coherent,
!                           (i.e. annihilation-free), and each determinant must only be specified once.
!                           The number of particles (with sign) on each determinant should be stored in MainSign,
!                           and there should not be a determinant entry with 'zero' particles associated with it.
!                           The fastest-moving index is associated with the bit-representation, i.e. 0 -> NIfTot
!                           This is returned as a list of all particles fully annihilated and merged maintaining order
!                           with the spawned list.
!       MaxMainInd          This is the size of the 'Main' lists (same on all processes).
!       TotDets             in: This is number of determinants specified in MainParticles on each process.
!                           out: This is the new number of determinants, having been annihilated and merged with the
!                                spawned list.
!       SpawnParts(:,:)     This is the list of particles to attempt to annihilate. Unlike the Main list, this list
!                           does *not* need to be ordered or sign coherent, and can also contain 'zero' sign particles.  
!                           Each particle contains its own sign
!       MaxSpawnInd         This is the size of the SpawnParts array.
!       SpawnDets           This is the number of spawned particles in SpawnParts.
!                           entry in the SpawnParts array.

        INTEGER, INTENT(IN) :: MaxMainInd,MaxSpawnInd
        INTEGER, INTENT(INOUT) :: TotDets
        type(fcimc_iter_data), intent(inout) :: iter_data
        INTEGER(KIND=n_int), INTENT(INOUT) , TARGET :: MainParts(0:NIfTot,MaxMainInd),SpawnParts(0:NIfTot,MaxSpawnInd)
!        INTEGER, INTENT(INOUT) , TARGET :: MainSign(MaxMainInd)
        INTEGER, INTENT(INOUT) :: SpawnDets
        INTEGER :: ierr,i
        CHARACTER(len=*) , PARAMETER :: this_routine='AnnihilationInterface'
!        INTEGER, DIMENSION(lenof_sign) :: TempSign
        TYPE(timer),save :: Annihil_time
        TYPE(timer),save :: Sync_time
        integer(kind=n_int), pointer,save :: SpawnVecLocal(:,:)
        Sync_time%timer_name='AnnihSync   innterface'
        call set_timer(Sync_time,20)

        if(tHashWalkerList) call stop_all(this_routine,"Cannot use annihilation interface with Hash table of particles")

        IF(.not.(ALLOCATED(ValidSpawnedList))) THEN
!This needs to be filled correctly before annihilation can take place.
            ALLOCATE(ValidSpawnedList(0:nNodes-1),stat=ierr)
        ENDIF
        call MPIBarrier(ierr)
        call halt_timer(Sync_time)
        Annihil_time%timer_name='Annihilation interface'
        call set_timer(Annihil_time,20)
        IF(.not.(ASSOCIATED(SpawnVecLocal))) THEN
!This is required scratch space of the size of the spawned arrays
            call shared_allocate_iluts("SpawnVecLocal",SpawnVecLocal,(/NIfTot,MaxSpawnInd/),iNodeIndex)
            ierr=0
            CALL LogMemAlloc('SpawnVecLocal',MaxSpawnInd*(NIfTot+1),size_n_int,this_routine,SpawnVec2Tag,ierr)
            call MPIBarrier(ierr)
!            SpawnVecLocal(:,:)=0
!            ALLOCATE(SpawnSignVec2(0:MaxSpawnInd),stat=ierr)
!            CALL LogMemAlloc('SpawnSignVec2',MaxSpawnInd+1,size_n_int,this_routine,SpawnSignVec2Tag,ierr)
!            SpawnSignVec2(:)=0
        ENDIF
!ValidSpawnedList indicates the next free index for each processor.
!For CCMC using shared memory, we will have a single processor on each node handling this.
!Say there are 2 nodes with 4 processors on each.  ValidSpawnedList will contain
! (#Det on Node 1)
! (#End of Node 1's spawned list)
! (#End of Node 1's spawned list)
! (#End of Node 1's spawned list)
! (#Det on Node 2)
! (#End of Node 2's spawned list)
! (#End of Node 2's spawned list)
! (#End of Node 2's spawned list)

! Since we've no way of knowing about nodes as yet, we just assume all processors are on the same node.  AJWT  TODO Multinodes.

!        ValidSpawnedList(0)=SpawnDets+1   !Add one since it always indicates the next free slot.
!        do i=1,nNodes-1
!            ValidSpawnedList(i)=MaxSpawnInd
!        enddo 

!        TempSign=0
!        do i=1,TotDets
!            TempSign(1)=MainSign(i)
!            call encode_sign(MainParts(:,i),TempSign)
!        enddo
! The SpawnParts already have their signs inside them

        MaxWalkersPart=MaxMainInd
!Point at correct arrays... will need to sort out how these are swapped in the main routine.
        CurrentDets => MainParts
        SpawnedParts => SpawnParts
!These point to the scratch space
        SpawnedParts2 => SpawnVecLocal

!          CALL DirectAnnihilation(TotDets, iter_data,.true.) !.true. for single processor annihilation
        CALL DirectAnnihilation(TotDets, iter_data,.false.) !.true. for single processor annihilation
!       W
        !if(iProcIndex==root) then        
!Signs put back again into seperate array
!           do i=1,TotDets
!               call extract_sign(CurrentDets(:,i),TempSign)
!               MainSign(i)=TempSign(1)
!           enddo
!         endif
         call MPIBarrier(ierr)

        call halt_timer(Annihil_time)
    END SUBROUTINE AnnihilationInterface


!This is a new annihilation algorithm. In this, determinants are kept on predefined processors, 
!and newlyspawned particles are sent here so that all the annihilations are
!done on a predetermined processor, and not rotated around all of them.
    SUBROUTINE DirectAnnihilation(TotWalkersNew, iter_data, tSingleProc)
        use bit_reps, only: test_flag
        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        INTEGER :: MaxIndex,ierr
        INTEGER(Kind=n_int) , POINTER :: PointTemp(:,:)
        logical, intent(in) :: tSingleProc
        TYPE(timer),save :: Compress_time
        integer :: i

!This routine will send all the newly-spawned particles to their correct processor. 
!MaxIndex is returned as the new number of newly-spawned particles on the processor. May have duplicates.
!The particles are now stored in SpawnedParts2/SpawnedSign2.
!        call WriteExcitorListP2(6,SpawnedParts,InitialSpawnedSlots,ValidSpawnedList,0,"Local")
!        if(bNodeRoot) 
        CALL SendProcNewParts(MaxIndex,tSingleProc)

!CompressSpawnedList works on SpawnedParts arrays, so swap the pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp

!Now we want to order and compress the spawned list of particles. 
!This will also annihilate the newly spawned particles amongst themselves.
!MaxIndex will change to reflect the final number of unique determinants in the newly-spawned list, 
!and the particles will end up in the spawnedSign/SpawnedParts lists.
!        WRITE(6,*) "Transferred",MaxIndex

        Compress_time%timer_name='Compression interface'
        call set_timer(Compress_time,20)

        CALL CompressSpawnedList(MaxIndex, iter_data)  

!        if(bNodeRoot) call sort(SpawnedParts(:,1:MaxIndex), ilut_lt, ilut_gt)
        call halt_timer(Compress_time)

!        WRITE(6,*) "List compressed",MaxIndex,TotWalkersNew

         ! If the semi-stochastic approach is being used then the following routine performs the
         ! annihilation of the deterministic states. These states are subsequently skipped in the
         ! AnnihilateSpawnedParts routine.
         if (tSemiStochastic) call deterministic_annihilation(iter_data)

!Binary search the main list and copy accross/annihilate determinants which are found.
!This will also remove the found determinants from the spawnedparts lists.
        CALL AnnihilateSpawnedParts(MaxIndex,TotWalkersNew, iter_data)  
        
        CALL set_timer(Sort_Time,30)
        if(tHashWalkerList) then
            call CalcHashTableStats(TotWalkersNew) 
        else
!Put the surviving particles in the main list, maintaining order of the main list (unless tHashWalkerList specified).
!Now we insert the remaining newly-spawned particles back into the original list (keeping it sorted), 
!and remove the annihilated particles from the main list.
            CALL InsertRemoveParts(MaxIndex, TotWalkersNew, iter_data)
        endif
        CALL halt_timer(Sort_Time)
        
    END SUBROUTINE DirectAnnihilation

!This routine is used for sending the determinants to the correct processors. 
    SUBROUTINE SendProcNewParts(MaxIndex,tSingleProc)
        REAL :: Gap
        integer :: i, j, error
        integer(MPIArg), dimension(nProcessors) :: sendcounts, disps, &
                                                   recvcounts, recvdisps
        INTEGER :: MaxSendIndex
        INTEGER, INTENT(OUT) :: MaxIndex
        LOGICAL, intent(in) :: tSingleProc

        if(tSingleProc) then
!Put all particles and gap on one proc.

!ValidSpawnedList(0:nNodes-1) indicates the next free index for each processor (for spawnees from this processor)
!  i.e. the list of spawned particles has already been arranged so that newly spawned particles are 
!grouped according to the processor they go to.

! sendcounts(1:) indicates the number of spawnees to send to each processor
! disps(1:) is the index into the spawned list of the beginning of the list to send to each processor (0-based)
           sendcounts(1)=int(ValidSpawnedList(0)-1,MPIArg)
           disps(1)=0
           if(nNodes>1) then
              sendcounts(2:nNodes)=0
              disps(2:nNodes)=int(ValidSpawnedList(1),MPIArg)
           endif
                                                                       
        else
!Distribute the gaps on all procs
           do i=0,nProcessors-1
               if(NodeRoots(ProcNode(i))==i) then  !This is a root 
                  sendcounts(i+1)=int(ValidSpawnedList(ProcNode(i))-    &
                        InitialSpawnedSlots(ProcNode(i)),MPIArg)
! disps is zero-based, but InitialSpawnedSlots is 1-based
                  disps(i+1)=int(InitialSpawnedSlots(ProcNode(i))-1,MPIArg)
               else
                  sendcounts(i+1)=0
                  disps(i+1)=disps(i)
               endif
           enddo
        endif


!        WRITE(6,*) 'sendcounts',sendcounts
!        WRITE(6,*) 'disps',disps

        MaxSendIndex=ValidSpawnedList(nNodes-1)-1

!We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
        recvcounts(1:nProcessors)=0
        
        call MPIBarrier(error)
        CALL set_timer(Comms_Time,30)

        CALL MPIAlltoAll(sendcounts,1,recvcounts,1,error)

!We can now get recvdisps from recvcounts, since we want the data to be contiguous after the move.
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo
        MaxIndex=recvdisps(nProcessors)+recvcounts(nProcessors)
        IF(tFillingStochRDMonFly.and.(.not.tHF_Ref_Explicit)) THEN            
            ! When we are filling the RDM, the SpawnedParts array contains 
            ! | Dj (0:NIfTot) | Di (0:NIfDBO) | Ci (1) | 
            ! All this needs to be passed around to the processor where Dj will be stored if 
            ! already occupied.  We then need to search the CurrentDets of that processor 
            ! to find Cj - while remembering the Di (and Ci) it goes with.
            do i=1,nProcessors
                recvdisps(i)=recvdisps(i)*int(NIfTot+NIfDBO+3,MPIArg)
                recvcounts(i)=recvcounts(i)*int(NIfTot+NIfDBO+3,MPIArg)
                sendcounts(i)=sendcounts(i)*int(NIfTot+NIfDBO+3,MPIArg)
                disps(i)=disps(i)*int(NIfTot+NIfDBO+3,MPIArg)
            enddo
        ELSE
            do i=1,nProcessors
                recvdisps(i)=recvdisps(i)*int(NIfTot+1,MPIArg)
                recvcounts(i)=recvcounts(i)*int(NIfTot+1,MPIArg)
                sendcounts(i)=sendcounts(i)*int(NIfTot+1,MPIArg)
                disps(i)=disps(i)*int(NIfTot+1,MPIArg)
            enddo
        ENDIF

!Max index is the largest occupied index in the array of hashes to be ordered in each processor 
        IF(MaxIndex.gt.(0.9*MaxSpawned)) THEN
            write(6,*) MaxIndex,MaxSpawned
            CALL Warning_neci("SendProcNewParts","Maximum index of newly-spawned array is " &
            & //"close to maximum length after annihilation send. Increase MemoryFacSpawn")
        ENDIF

!        WRITE(6,*) "sendcounts: ",sendcounts(:)
!        WRITE(6,*) "disps: ",disps(:)
!        WRITE(6,*) "recvcounts: ", recvcounts(:)
!        WRITE(6,*) "recvdisps: ",recvdisps(:)
!        WRITE(6,*) "Sent Sign: ", NINT(Gap),sendcounts(2)
!        do i=NINT(Gap),NINT(Gap)+sendcounts(2)
!            WRITE(6,*) i,"***",SpawnedSign(i)
!        enddo
        
!This is the main send of newly-spawned particles and signs to each determinants correct processor.
!        CALL MPIAlltoAllvI(SpawnedSign(1:MaxSendIndex),sendcounts,disps,SpawnedSign2(1:MaxIndex),recvcounts,recvdisps,error)
        
!        WRITE(6,*) MaxIndex, "Recieved signs: "
!        do i=1,MaxIndex
!            WRITE(6,*) SpawnedSign2(i)
!        enddo

!Update the number of integers we need to send.
!        do i=1,nNodes
!            sendcounts(i)=sendcounts(i)*(NIfTot+1)
!            disps(i)=disps(i)*(NIfTot+1)
!            recvcounts(i)=recvcounts(i)*(NIfTot+1)
!            recvdisps(i)=recvdisps(i)*(NIfTot+1)
!        enddo

!        WRITE(6,*) "Sent Particles: ", NINT(Gap),sendcounts(2)
!        do i=NINT(Gap)+1,NINT(Gap)+sendcounts(2)
!            write(6,*) i, '***', CountBits(spawnedparts(:,i), nifd)
!            WRITE(6,*) i,"***",SpawnedParts(:,i)
!        enddo

!        write(6,*) "Second all to all"
!        call neci_flush(6)
        CALL MPIAlltoAllv(SpawnedParts,sendcounts,disps,SpawnedParts2,recvcounts,recvdisps,error)
!        write(6,*) "Second all to all finish"
!        call neci_flush(6)

!        WRITE(6,*) MaxIndex, "Recieved particles: "
!        do i=1,MaxIndex
!            WRITE(6,*) SpawnedParts2(:,i)
!        enddo
        
        CALL halt_timer(Comms_Time)

    END SUBROUTINE SendProcNewParts

!This sorts and compresses the spawned list to make it easier for the rest of the annihilation process.
!This is not essential, but should proove worthwhile
    SUBROUTINE CompressSpawnedList(ValidSpawned, iter_data)
        type(fcimc_iter_data), intent(inout) :: iter_data
        INTEGER :: VecInd,ValidSpawned,DetsMerged,i,BeginningBlockDet,FirstInitIndex,CurrentBlockDet
        real(dp) :: SpawnedSign(lenof_sign), Temp_Sign(lenof_sign)
        integer :: EndBlockDet, part_type, StartCycleInit, cum_count, j, Parent_Array_Ind
        integer :: No_Spawned_Parents
        LOGICAL :: tSuc, tInc
        INTEGER(Kind=n_int) , POINTER :: PointTemp(:,:)
        integer(n_int) :: cum_det (0:niftot)
        CHARACTER(len=*), parameter :: this_routine='CompressSpawnedList'
        TYPE(timer),save :: Sort_time

        !print *, "SpawnedParts before:"
        !do j = 1, ValidSpawned
        !    print *, SpawnedParts(:,j)
        !end do
        !print *,

!We want to sort the list of newly spawned particles, in order for quicker binary searching later on. 
!(this is not essential, but should proove faster)
!They should remain sorted after annihilation between spawned
        if(.not.bNodeRoot) return
!        call WriteExcitorListP(6,SpawnedParts,0,ValidSpawned,0,"UnOrdered")

        index_of_first_non_determ = 1

        Sort_time%timer_name='Compress Sort interface'
        call set_timer(Sort_time,20)

        ! Note: If tSortDetermToTop is true then the sorting routine moves all deterministic states
        ! to the top of the list, for optimisation elsewhere. If tSortDetermToTop is false but we are
        ! still performing a semi-stochastic calculation, we want to sort the spawned list so that psips
        ! which are spawned from deterministic parents are separated from psips which aren't, so that they
        ! can be considered separately later (if a psip is spawned from a deterministic parent into the
        ! deterministic space then it should be aborted - we only know if the state it occupies is
        ! deterministic later on, if tSortDetermToTop is false).
        ! If not performing a semi-stochastic simulation then these conditions are just ignored.
        call sort(SpawnedParts(:,1:ValidSpawned), ilut_lt, ilut_gt)
        
        CALL halt_timer(Sort_time)

        IF(tHistSpawn) HistMinInd2(1:NEl)=FCIDetIndex(1:NEl)

!        call WriteExcitorListP(6,SpawnedParts,0,ValidSpawned,0,"Ordered")


!!        WRITE(6,*) "************ - Ordered",ValidSpawned,NIfTot,Iter
!!        do i=1,ValidSpawned
!!            WRITE(6,*) i,SpawnedParts(0:NIfTot-1,i),SpawnedParts(NIfTot,i)-2
!!        enddo

!First, we compress the list of spawned particles, so that they are only specified at most once in each processors list.
!During this, we transfer the particles from SpawnedParts to SpawnedParts2
!If we are working with complex walkers, we essentially do the same thing twice, annihilating real and imaginary
!particles seperately.
        VecInd=1    !This is the index in the SpawnedParts2 array to copy the compressed walkers into
        !BeginningBlockDet will indicate the index of the first entry for a given determinant in SpawnedParts
        BeginningBlockDet=1         
        DetsMerged=0
        Parent_Array_Ind = 1
        Spawned_Parts_Zero = 0

        do while(BeginningBlockDet.le.ValidSpawned)
            !loop in blocks of the same determinant to the end of the list of walkers

            FirstInitIndex=0
            CurrentBlockDet=BeginningBlockDet+1

            do while(CurrentBlockDet.le.ValidSpawned)
                if(.not.(DetBitEQ(SpawnedParts(0:NIfTot,BeginningBlockDet),SpawnedParts(0:NIfTot,CurrentBlockDet),NIfDBO))) exit
                if (tSemiStochastic .and. (.not. tSortDetermToTop)) then
                    if (.not. (test_flag(SpawnedParts(0:NIfTot,BeginningBlockDet), flag_determ_parent) .eqv. &
                               test_flag(SpawnedParts(0:NIfTot,CurrentBlockDet), flag_determ_parent))) exit
                end if
                ! Loop over walkers on the same determinant in SpawnedParts.
                ! Also, seperate out states which are and aren't spawned from the deterministic space.
                ! (Ignore the options relating to these deterministic flags if not concerned with the
                ! deterministic code).
                CurrentBlockDet=CurrentBlockDet+1
            enddo

            EndBlockDet=CurrentBlockDet-1 !EndBlockDet indicates that we have reached the end of the block of similar dets
!            WRITE(6,*) "Found Block: ",BeginningBlockDet," -> ",EndBlockDet

            if(EndBlockDet.eq.BeginningBlockDet) then
                !Optimisation: This block only consists of one entry. Simply copy it across rather than 
                !               explicitly searching the list.
                
                SpawnedParts2(:,VecInd)=SpawnedParts(:,BeginningBlockDet)   !Transfer all info to the other array

                IF(tFillingStochRDMonFly.and.(.not.tHF_Ref_Explicit)) THEN
                    ! SpawnedParts contains the determinants spawned on (Dj), and it's parent (Di) plus it's sign (Cj).
                    ! As in | Dj | Di | Ci |
                    ! We then compress multiple occurances of Dj, but these may have come from different parents, and 
                    ! we want to keep track of all Di's.
                    ! As we compress SpawnedParts, we therefore move all the parents (Di's) into Spawned_Parents.
                    ! If the compressed Dj is at position VecInd in SpawnedParts, then Spawned_Parents_Index(1,VecInd) 
                    ! is the starting point of it's parents (Di) in Spawned_Parents, and there are 
                    ! Spawned_Parents_Index(2,VecInd) entries corresponding to this Dj.
                    if(.not.(DetBitZero(SpawnedParts(NIfTot+1:NIfTot+NIfDBO+1,BeginningBlockDet),NIfDBO))) then
                        ! If the parent determinant is null, the contribution to the RDM is zero.  
                        ! No point in doing anything more with it.

                        Spawned_Parents(0:NIfDBO+1,Parent_Array_Ind) = SpawnedParts(NIfTot+1:NIfTot+NIfDBO+2,BeginningBlockDet)
                        ! The first NIfDBO of the Spawned_Parents entry is the parent determinant, the NIfDBO + 1 entry 
                        ! is the biased Ci. Parent_Array_Ind keeps track of the position in Spawned_Parents.
                        Spawned_Parents_Index(1,VecInd) = Parent_Array_Ind
                        Spawned_Parents_Index(2,VecInd) = 1
                        ! In this case there is only one instance of Dj - so therefore only 1 parent Di.
                        Parent_Array_Ind = Parent_Array_Ind + 1
                    else
                        Spawned_Parents_Index(1,VecInd) = Parent_Array_Ind
                        Spawned_Parents_Index(2,VecInd) = 0
                    endif
                    call extract_sign (SpawnedParts(:,BeginningBlockDet), temp_sign)
                    if(IsUnoccDet(temp_sign)) then
                        Spawned_Parts_Zero = Spawned_Parts_Zero + 1
                    endif
                ENDIF

                ! If this is a deterministic state then update the value of index_of_first_non_determ, which we eventually
                ! want to store the index of the first non-deterministic state in the space.
                ! Note, this is only used if tSortDetermToTop is true. If not true then the deterministic flag won't have
                ! been set.
                if (tSemiStochastic) then
                    if (test_flag(SpawnedParts(:,BeginningBlockDet), flag_deterministic)) index_of_first_non_determ = VecInd + 1
                end if

                VecInd=VecInd+1
                BeginningBlockDet=CurrentBlockDet           !Move onto the next block of determinants
                CYCLE   !Skip the rest of this block
            endif

            ! Reset the cumulative determinant
            cum_det = 0
            cum_det (0:nifdbo) = SpawnedParts(0:nifdbo, BeginningBlockDet)

            IF(tFillingStochRDMonFly.and.(.not.tHF_Ref_Explicit)) THEN
                ! This is the first Dj determinant - set the index for the beginning of where 
                ! the parents for this Dj can be found in Spawned_Parents.
                Spawned_Parents_Index(1,VecInd) = Parent_Array_Ind
 
                ! In this case, multiple Dj's must be compressed, and therefore the Di's dealt with as 
                ! described above. We first just initialise the position in the Spawned_Parents array to enter the Di's.
                Spawned_Parents_Index(2,VecInd) = 0
            ENDIF
            
            do part_type=1,lenof_sign   !Annihilate in this block seperately for real and imag walkers

                ! How many of either real/imaginary spawns are there onto each det
                cum_count = 0
                
!                WRITE(6,*) "Testing particle types: ",part_type
                    
                if(tTruncInitiator) then
                    !Need to find if there are any initiators in this block
                    do i=BeginningBlockDet,EndBlockDet  !Loop over the block
                        if (test_flag(SpawnedParts(:,i), flag_parent_initiator(part_type))) then
                            !Found an initiator walker
!                            WRITE(6,*) "Found another initiator: ",i
                            if (tHistSpawn) then
                                call extract_sign (SpawnedParts(:,i), SpawnedSign)
                                call extract_sign (SpawnedParts(:,BeginningBlockDet), temp_sign)
                                call HistAnnihilEvent(SpawnedParts, SpawnedSign, temp_sign, part_type)
                            endif
                            call FindResidualParticle (cum_det, SpawnedParts(:,i), cum_count, part_type, iter_data, &
                                                            VecInd, Parent_Array_Ind)
                        endif
                    enddo
                endif

!                WRITE(6,*) "After initiators: ",Cum_Sign(part_type),Cum_Flag

                !Now loop over the same block again, but this time calculating the contribution from non-initiators
                !We want to loop over the whole block.
                do i=BeginningBlockDet,EndBlockDet
                    tInc = .true.
                    if (tTruncInitiator) then
                        if (test_flag (SpawnedParts(:,i), flag_parent_initiator(part_type))) tInc = .false.
                    endif
                    if (tInc) then
                        ! If truncinitiator, only consider the non-initiators here 
                        ! (the initiators have already been dealt with).
                        if (tHistSpawn) then
                            call extract_sign (SpawnedParts(:,i), SpawnedSign)
                            call extract_sign (SpawnedParts(:,BeginningBlockDet), temp_sign)
                            call HistAnnihilEvent (SpawnedParts, SpawnedSign, temp_sign, part_type)
                        endif
                        call FindResidualParticle (cum_det, SpawnedParts(:,i), cum_count, part_type, iter_data, &
                                                        VecInd, Parent_Array_Ind)
                    endif
                enddo

            enddo ! End loop over particle type

            ! Copy details into the final array
            call extract_sign (cum_det, temp_sign)

            if ((sum(abs(temp_sign)) > 0.0_dp).or.(tFillingStochRDMonFly.and.(.not.tHF_Ref_Explicit))) then
                ! Transfer all info into the other array.
                ! Usually this is only done if the final sign on the compressed Dj is not equal to zero.
                ! But in the case of the stochastic RDM, we are concerned with the sign of Dj in the CurrentDets 
                ! array - not the newly spawned sign.  We still want to check if Dj has a non-zero Cj in Current Dets - 
                ! so we need to carry this Dj through to the stage of checking CurrentDets regardless of the sign here. 
                ! Also getting rid of them here would make the biased sign of Ci slightly wrong.

                ! If this is a deterministic state then update the value of index_of_first_non_determ, which we eventually
                ! want to store the index of the first non-deterministic state in the space.
                ! Note, this is only used if tSortDetermToTop is true. If not true then the deterministic flag won't have
                ! been set.
                if (tSemiStochastic) then
                    if (test_flag(SpawnedParts(:,BeginningBlockDet), flag_deterministic)) index_of_first_non_determ = VecInd + 1
                end if

                SpawnedParts2(0:NIfTot,VecInd) = cum_det(0:NIfTot)
                VecInd = VecInd + 1
                DetsMerged = DetsMerged + EndBlockDet - BeginningBlockDet

                ! Spawned_Parts_Zero is the number of spawned parts that are zero after compression 
                ! of the spawned_parts list - and should have been 
                ! removed from SpawnedParts if we weren't calculating the RDM. - need this for a check later.
                if(IsUnoccDet(temp_sign)) Spawned_Parts_Zero = Spawned_Parts_Zero + 1
            else
                ! All particles from block have been annihilated.
                DetsMerged = DetsMerged + EndBlockDet - BeginningBlockDet + 1
                ! This spawned entry will be removed - don't want to store any parents.
                ! Reset Parent_Array_Ind so that the parents will be written over.
            endif

            BeginningBlockDet=CurrentBlockDet           !Move onto the next block of determinants

        enddo   

        IF(tFillingStochRDMonFly.and.(.not.tHF_Ref_Explicit)) No_Spawned_Parents = Parent_Array_Ind - 1
        ValidSpawned=ValidSpawned-DetsMerged    !This is the new number of unique spawned determinants on the processor
        IF(ValidSpawned.ne.(VecInd-1)) THEN
            CALL Stop_All(this_routine,"Error in compression of spawned list")
        ENDIF

!Want the compressed list in spawnedparts at the end of it - swap pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp

!        WRITE(6,*) 'Spawned after compress'
!        do i = 1, ValidSpawned
!            WRITE(6,*) SpawnedParts(:,i)
!        enddo
!

!        WRITE(6,*) 'Spawned Parents'
!        do i = 1, No_Spawned_Parents
!            WRITE(6,*) Spawned_Parents(:,i)
!        enddo
!
!        WRITE(6,*) 'Spawned Parents Index'
!        do i = 1, ValidSpawned
!            WRITE(6,*) Spawned_Parents_Index(:,i)
!        enddo

!        IF(tFillingRDMonFly) THEN
!            WRITE(6,*) 'Parents'
!            do i = 1, No_Spawned_Parents
!                WRITE(6,*) Spawned_Parents(:,i)
!            enddo
!
!            WRITE(6,*) 'Index'
!            do i = 1, No_Spawned_Parents
!                WRITE(6,*) Spawned_Parents_Index(:,i)
!            enddo
!        ENDIF

!        WRITE(6,*) "************************"
!        call WriteExcitorListP(6,SpawnedParts,0,ValidSpawned,0,"Compressed Spawned")

!!        WRITE(6,*) "************************"
!!        WRITE(6,*) "Compressed List: ",ValidSpawned,DetsMerged
!!        do i=1,ValidSpawned
!!            WRITE(6,*) SpawnedParts(0:NIfTot-1,i),SpawnedParts(NIfTot,i)-2
!!        enddo
!!        WRITE(6,*) "***","iter_data%naborted: ",iter_data%naborted
!!        CALL neci_flush(6)
        
    END SUBROUTINE CompressSpawnedList

!Histogram a possible annihilation event
    subroutine HistAnnihilEvent(iLut,Sign1,Sign2,part_type)
        implicit none
        integer(kind=n_int), intent(in) :: iLut(0:NIfTot)
        real(dp), dimension(lenof_sign), intent(in) :: Sign1,Sign2
        integer, intent(in) :: part_type
        integer :: ExcitLevel,PartIndex
        logical :: tSuc

!We want to histogram where the particle annihilations are taking place.
        if((Sign1(part_type)*Sign2(part_type)).ge.0.0) return   !No annihilation occuring - particles same sign

        ExcitLevel = FindBitExcitLevel(iLut,iLutHF, nel)
        IF(ExcitLevel.eq.NEl) THEN
            CALL BinSearchParts2(iLut(:),HistMinInd2(ExcitLevel),Det,PartIndex,tSuc)
            HistMinInd2(ExcitLevel)=PartIndex
        ELSEIF(ExcitLevel.eq.0) THEN
            PartIndex=1
            tSuc=.true.
        ELSE
            CALL BinSearchParts2(iLut(:),HistMinInd2(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartIndex,tSuc)
            HistMinInd2(ExcitLevel)=PartIndex
        ENDIF
        IF(tSuc) THEN
            AvAnnihil(part_type,PartIndex)=AvAnnihil(part_type,PartIndex)+ &
                2*(MIN(abs(Sign1(part_type)),abs(Sign2(part_type))))
            InstAnnihil(part_type,PartIndex)=InstAnnihil(part_type,PartIndex)+ &
                2*(MIN(abs(Sign1(part_type)),abs(Sign2(part_type))))
        ELSE
            CALL Stop_All("CompressSpawnedList","Cannot find corresponding FCI determinant when histogramming")
        ENDIF

    end subroutine HistAnnihilEvent

    !This routine is called when compressing the spawned list.
    !It takes the sign and flags from two particles on the same determinant, 
    !and calculates what the residual particles and signs should be from them.
    !This deals with real and imaginary signs seperately, and so the 'signs' are integers.

    subroutine FindResidualParticle (cum_det, new_det, cum_count, part_type, &
                                     iter_data, Spawned_No, Parent_Array_Ind)

        ! This routine is called whilst compressing the spawned list during
        ! annihilation. It considers the sign and flags from two particles
        ! on the same determinant, and calculates the correct sign/flags
        ! for the compressed particle.
        !
        ! --> The information is stored within the first particle in a block
        ! --> Should be called for real/imaginary particles seperately

        integer(n_int), intent(inout) :: cum_det(0:nIfTot)
        integer(n_int), intent(in) :: new_det(0:niftot+nifdbo+2)
        integer, intent(inout) :: cum_count
        integer, intent(in) :: part_type, Spawned_No 
        integer, intent(inout) :: Parent_Array_Ind
        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp) :: new_sgn, cum_sgn, updated_sign, sgn_prod

        ! Obtain the signs and sign product. Ignore new particle if zero.
        new_sgn = extract_part_sign (new_det, part_type)

        if (new_sgn == 0) then
            ! New sign is just an entry from SpawnedParts - this should only ever be zero
            ! in the complex case. 
            ! If it is 0 and we're not filling the RDM (and therefore filling up the 
            ! Spawned_Parents array), can just ignore the zero entry.
            if(.not.tFillingStochRDMonFly) return
        endif
        cum_sgn = extract_part_sign (cum_det, part_type)
        
        sgn_prod = cum_sgn * new_sgn

        ! If we are including this det, then increment the count
        cum_count = cum_count + 1

        if (tTruncInitiator) then
            !The rules are as follows:
            !If particles are of the same sign:
                !Init + Init = Init
                !Init + Non = Init
                !Non + Non = Init

            !If particles are of opposite sign:
                !Init + Init = Init
                !Non + Init = Whichever is largest
                !Non + Non = Non

            ! If flag_make_initiator is set, it always becomes an initiator.
            ! n.b. we are using flag_is_initiator and flag_parent_initiator
            !      somwhat interchangably here...
            if (sgn_prod > 0.0_dp) then
                ! Signs are the same. This must be an initiator.
                ! Equivalent to (deprecated) tKeepDoubSpawns
                call set_flag (cum_det, flag_is_initiator(part_type))
            elseif (sgn_prod < 0.0_dp) then
                ! Annihilating (serially)
                ! --> Retain the initiator flag of the largest term.
                if (abs(new_sgn) > abs(cum_sgn)) then
                    call copy_flag (new_det, cum_det, &
                                    flag_is_initiator(part_type))
                endif
            else
                ! cum_sgn == 0, new_sgn /= 0, therefore just take the flags
                ! from the new particles.
                call encode_flags (cum_det, extract_flags (new_det))
            endif

            ! If we have set the make_initiator flag, then the target
            ! particle must become an initiator.
            if (test_flag (new_det, flag_make_initiator(part_type))) then
                ! If make_initiator is set, then the particle must become
                ! an initiator.
                call set_flag (cum_det, flag_make_initiator(part_type))
                call set_flag (cum_det, flag_is_initiator(part_type))
            endif

            ! This could be carried via more flags, but this is cleaner.
            if (cum_count > 2) then
                call set_flag (cum_det, flag_parent_initiator(part_type))
            endif
        endif

        if (tSemiStochastic) then
            if (test_flag(new_det, flag_deterministic)) call set_flag(cum_det, flag_deterministic)
            if (test_flag(new_det, flag_determ_parent)) call set_flag(cum_det, flag_determ_parent)
        end if

        ! Update annihilation statistics (is this really necessary?)
        if (sgn_prod < 0.0_dp) then
            Annihilated = Annihilated + 2*min(abs(cum_sgn), abs(new_sgn))
            iter_data%nannihil(part_type) = iter_data%nannihil(part_type)&
                + 2 * min(abs(cum_sgn), abs(new_sgn))
        endif

        ! Update the cumulative sign count
        updated_sign = cum_sgn + new_sgn
        call encode_part_sign (cum_det, updated_sign, part_type)

        ! Obviously only add the parent determinant into the parent array if it is 
        ! actually being stored - and is therefore not zero.
        if((tFillingStochRDMonFly.and.(.not.tHF_Ref_Explicit)).and.&
            (.not.DetBitZero(new_det(NIfTot+1:NIfTot+NIfDBO+1),NIfDBO))) then
            ! No matter what the final sign is, always want to add any Di stored in 
            ! SpawnedParts to the parent array.
            Spawned_Parents(0:NIfDBO+1,Parent_Array_Ind) = new_det(NIfTot+1:NIfTot+NIfDBO+2)
            Parent_Array_Ind = Parent_Array_Ind + 1
            Spawned_Parents_Index(2,Spawned_No) = Spawned_Parents_Index(2,Spawned_No) + 1
        endif

    end subroutine FindResidualParticle

    subroutine deterministic_annihilation(iter_data)

        ! A general explanation of this routine, and also the changes made in other annihilation
        ! routines for the semi-stochastic code:

        ! We have two sets of psips spawned onto the deterministic states. Firstly, those which
        ! are stored in the vector partial_determ_vector, which comes from the deterministic projection.
        ! Secondly, some of the states in the SpawnedParts, spawned from outside of the deterministic
        ! space, may also be in the deterministic space. *All* of the psips spawned onto the
        ! deterministic space are kept, even if they are not spawned from initiators. We therefore
        ! can simply move all of the states from both lists into the main list.

        ! partial_determ_vector is already in the same order as the the states in the main list (note
        ! that *all* deterministic states are always kept in the main list, even if they have zero
        ! weight). It is straightforward therefore to loop through these states and add in the
        ! correct weights.

        ! There are two ways that the annihilation of the SpawnedParts psips can occur, depending on
        ! the deterministic states being used. If it can be checked quickly if a state is in the
        ! deterministic space then the deterministic flag is already set in SpawnedParts and the
        ! ordering of the states is defined so that all the deterministic states are at the top of
        ! the list (tSortDetermToTop is true). In this case, index_of_first_non_determ denotes the index
        ! in SpawnedParts of the first non-deterministic state, and therefore index_of_first_non_determ-1
        ! is the index of the last deterministic state. Hence, we only need to perform a binary search
        ! over these states in the main list, and can skip these states out in AnnihilateSpawnedParts
        ! as an optimisation, as in that routine we the only deal with non-deterministic states.

        ! If the deterministic states are not sorted to the top (if tSortDetermToTop is false), then
        ! this optimisation is not possible, so we simply check in AnnihilateSpawnedParts if each
        ! state is deterministic. Because it is not possible to check if a state is deterministic
        ! earlier, it is also not possible to abort it ealier if it is spawned from the deterministic
        ! space (this type of spawning is done deterministically). Hence, in CompressSpawned list, we
        ! separate out states which do and do not have flag_determ_parent set. Then, in
        ! AnnihilateSpawnedParts, once the state is found in CurrentDets, we check if the state is
        ! deterministic and can then cancel the spawning. Also, because the same state can be
        ! specified twice in SpawnedParts, once with the determ_parent flag set, and once with it
        ! not set, when not spawning to the deterministic space, we merge these two states together.

        ! Also, regarding the initiator approximation, all spawning attempts to the deterministic
        ! space are allowed, and if any of the psips spawned on a site are from the deterministic space
        ! then the spawning is always allowed.

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: i, MinInd, MaxInd, PartInd
        real(dp), dimension(lenof_sign) :: SpawnedSign, CurrentSign, SignProd
        logical :: tSuccess

        ! First, copy across the weights from partial_determ_vector:
        do i = 1, deterministic_proc_sizes(iProcIndex)
            call extract_sign(CurrentDets(:, indices_of_determ_states(i)), CurrentSign)
            SpawnedSign = partial_determ_vector(i)
            call encode_sign(CurrentDets(:, indices_of_determ_states(i)), SpawnedSign + CurrentSign)

            ! Update stats:
            ! Number born:
            iter_data%nborn(1) = iter_data%nborn(1) + abs(SpawnedSign(1))
            ! Number annihilated:
            SignProd = CurrentSign*SpawnedSign
            if (SignProd(1) < 0.0_dp) iter_data%nannihil(1) = iter_data%nannihil(1) + &
                2*(min(abs(CurrentSign(1)), abs(SpawnedSign(1))))
        end do

        ! Second, *if* the deterministic states are sorted to the top, copy across the weights from
        ! SpawnedParts to CurrentDets here:

        if (tSortDetermToTop) then
            ! We only need to search the deterministic part of CurrentDets.
            MinInd = 1
            MaxInd = deterministic_proc_sizes(iProcIndex)
            do i = 1, index_of_first_non_determ-1
                ! Search the CurrentDets list for this state.
                call BinSearchParts(SpawnedParts(:, i), MinInd, MaxInd, PartInd, tSuccess)
                if (tSuccess) then
                    call extract_sign(CurrentDets(:, PartInd), CurrentSign)
                    call extract_sign(SpawnedParts(:, i), SpawnedSign)
                    call encode_sign(CurrentDets(:, PartInd), SpawnedSign + CurrentSign)
                    call encode_sign(SpawnedParts(:, i), null_part)

                    ! Update stats:
                    SignProd = CurrentSign*SpawnedSign
                    if (SignProd(1) < 0.0_dp) iter_data%nannihil(1) = iter_data%nannihil(1) + &
                        2*(min(abs(CurrentSign(1)), abs(SpawnedSign(1))))
                else
                    ! All deterministic states should be kept in CurrentDets always, so if this search
                    ! does not return a state then something has gone badly wrong!
                    call stop_all("deterministic_annihilation","Deterministic state not found in main list.")
                end if
            end do
        end if

    end subroutine
    
!In this routine, we want to search through the list of spawned particles. For each spawned particle, 
!we binary search the list of particles on the processor
!to see if an annihilation event can occur. The annihilated particles are then removed from the spawned list
!to the whole list of spawned particles at the end of the routine.
!In the main list, we change the 'sign' element of the array to zero. 
!These will be deleted at the end of the total annihilation step.
    SUBROUTINE AnnihilateSpawnedParts(ValidSpawned,TotWalkersNew, iter_data)
        use nElRDMMod , only : check_fillRDM_DiDj
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, intent(inout) :: TotWalkersNew
        integer, intent(inout) :: ValidSpawned 
        INTEGER :: MinInd,PartInd,i,j,ToRemove,DetsMerged,PartIndex
        real(dp), dimension(lenof_sign) :: CurrentSign, SpawnedSign, SignTemp
        real(dp), dimension(lenof_sign) :: SignProd, NewSignTemp
        REAL(dp) :: pRemove, r
        INTEGER :: ExcitLevel, nJ(NEl),DetHash,FinalVal,clash,walkExcitLevel, dettemp(NEl)
        INTEGER(KIND=n_int) , POINTER :: PointTemp(:,:)
        LOGICAL :: tSuccess,tSuc,tPrevOcc, tSkip, tAlwaysAllow
        character(len=*), parameter :: this_routine="AnnihilateSpawnedParts"
        integer :: comp

        !print *, "SpawnedParts after:"
        !do j = 1, ValidSpawned
        !    print *, SpawnedParts(:,j)
        !end do
        !print *,

        if(.not.bNodeRoot) return  !Only node roots to do this.

        CALL set_timer(AnnMain_time,30)

!        write(6,*) "Entering AnnihilateSpawnedParts"
!        call flush(6)

        ToRemove=0  !The number of particles to annihilate

        ! Logicals for semi-stochastic code.
        tSkip = .false.
        tAlwaysAllow = .false.

        ! If we are performing a semi-stochastic simulation *where the deterministic states
        ! are sorted to the top*, then annihilation in the deterministic space has already
        ! been taken care of. Therefore we should only loop over the non-deterministic states.
        ! In this case, index_of_first_non_determ gives the index of the first state in SpawnedParts
        ! which is not in the deterministic space, and therefore the first state to loop over. If
        ! not performing a semi-stochastic simulation then we want to loop over all states, so
        ! simply set this to 1.
        if (.not. tSortDetermToTop) then
            index_of_first_non_determ = 1
        else
            ToRemove = index_of_first_non_determ - 1
        end if

!MinInd indicates the minimum bound of the main array in which the particle can be found.
!Since the spawnedparts arrays are ordered in the same fashion as the main array, 
!we can find the particle position in the main array by only searching a subset.

        ! As noted above, if we are performing a semi-stochastic simulation with tSortDetermToTop
        ! true, then the deterministic states have been taken care of. So, we only need to perform the
        ! binary search over the non-deterministic states, so set MinInd to ignore deterministic states.
        if (tSortDetermToTop) then
            MinInd = deterministic_proc_sizes(iProcIndex) + 1
        else
            MinInd = 1
        end if

        IF(tHistSpawn) HistMinInd2(1:NEl)=FCIDetIndex(1:NEl)

!        WRITE(6,*) "Annihilating between ",ValidSpawned, " spawned particles and ",TotWalkersNew," original particles..."
!        WRITE(6,*) "SpawnedParts: "
!        do i=1,ValidSpawned
!            WRITE(6,*) SpawnedParts(:,i)
!        enddo
!        WRITE(6,*) "Original Parts: "
!        do i=1,TotWalkersNew
!            WRITE(6,*) CurrentDets(:,i)
!        enddo
        

!          write(6,*) "CurrentDets",TotWalkersNew
!          do i=1,TotWalkersNew
!             write(6,"(I5)",advance='no') i
!             call WriteBitDet(6,CurrentDets(:,i),.true.)
!          enddo
!             write(6,*) "ValidSpawned"
!             do i=1,ValidSpawned
!                write(6,"(I5)",advance='no') i
!                call WriteBitDet(6,SpawnedParts(:,i),.true.)
!             enddo

        CALL set_timer(BinSearch_time,45)

        do i = index_of_first_non_determ, ValidSpawned

            if (tSkip) then
                tSkip = .false.
                cycle
            end if

!This will binary search the CurrentDets array to find the desired particle. 
!tSuccess will determine whether the particle has been found or not.
!It will also return the index of the position one below where the particle would be found if was in the list.
!            CALL LinSearchParts(CurrentDets(:,1:TotWalkersNew),SpawnedParts(0:NIfD,i),MinInd,TotWalkersNew,PartInd,tSuccess)
!            WRITE(6,*) "MinInd",MinInd
!            do j=MinInd,min(MinInd+10,TotWalkersNew)
!               write(6,"(I5)",advance='no') j 
!               call WriteBitDet(6,CurrentDets(:,j),.true.)
!            enddo

            if(tHashWalkerList) then
                !Do not need to binary search list. It is not sorted, but there is a hash table to it.
                tSuccess=.false.
                call decode_bit_det (nJ, SpawnedParts(:,i))              
!                write(6,*) "Sending to hash func 1: ",nJ(:)
                DetHash=FindWalkerHash(nJ)
!                write(6,*) "DetHash: ",DetHash
                FinalVal=HashIndex(0,DetHash)-1
!                write(6,*) "FinalVal: ",FinalVal
                do clash=1,FinalVal
                    ASSERT(HashIndex(clash,DetHash).le.TotWalkersNew)
                    if(DetBitEQ(SpawnedParts(:,i),CurrentDets(:,HashIndex(clash,DetHash)),NIfDBO)) then
                        !We have found the matching determinant
                        tSuccess=.true.
                        PartInd=HashIndex(clash,DetHash)
!                        write(6,*) "Found it at: ",PartInd
                        exit
                    endif
                enddo
            else
                CALL BinSearchParts(SpawnedParts(:,i),MinInd,TotWalkersNew,PartInd,tSuccess)
            endif

!            WRITE(6,"(A,2I6,L)",advance="no") "Binary search complete: ",i,PartInd,tSuccess
!            call WriteBitDet(6,SpawnedParts(:,i),.true.)
!            call WriteBitDet(6,CurrentDets(:,PartInd),.true.)

!            WRITE(6,*) 'i,SpawnedParts(:,i)',i,SpawnedParts(:,i)

            ! Abort spawning from the deterministic space to the deterministic space, and also
            ! merge this state and next state in SpawnedParts if they are the same states but
            ! with the determ_parent flag set in the first state and not set in the second state.
            if (tSemiStochastic .and. (.not. tSortDetermToTop)) then
                if (test_flag(SpawnedParts(:,i), flag_determ_parent)) then
                    ! If psips spawned onto this state from within the deterministic space.
                    if (i < ValidSpawned) then
                        ! If not on the last state.
                        ! Is this state the same as the next state in SpawnedParts?
                        ! (comp is returned as 0 if true)
                        comp = DetBitLT(SpawnedParts(:, i), SpawnedParts(:, i+1), NIfDBO, .false.)
                    else
                        ! Set comp to anything but 0 in this case.
                        comp = 1
                    end if

                    if (test_flag(CurrentDets(:,PartInd), flag_deterministic) .and. tSuccess) then
                        ! If this state is in the deterministic space then abort it, as it was also spawned
                        ! from the deterministic space. Then simply move onto the next state in SpawnedParts,
                        ! even if it is the same state.
                        call extract_sign(SpawnedParts(:, i), SpawnedSign)
                        call encode_sign(SpawnedParts(:, i), null_part)
                        ToRemove = ToRemove + 1
                        MinInd = PartInd

                        ! Update stats:
                        iter_data%naborted(1) = iter_data%naborted(1) + abs(SpawnedSign(1))

                        cycle
                    else
                        if (comp == 0) then
                            ! If the next state is the same, and this state is not in the deterministic space
                            ! then merge the two amplitudes onto this one state, and then skip the next state
                            ! in SpawnedParts on the next loop.
                            call extract_sign(SpawnedParts(:, i), SpawnedSign)
                            call extract_sign(SpawnedParts(:, i+1), SignTemp)
                            call encode_sign(SpawnedParts(:, i), SpawnedSign + SignTemp)
                            call encode_sign(SpawnedParts(:, i+1), null_part)
                            ToRemove = ToRemove + 1
                            tSkip = .true.

                            ! Update stats:
                            SignProd = SpawnedSign*SignTemp
                            if (SignProd(1) < 0.0_dp) iter_data%nannihil(1) = iter_data%nannihil(1) + &
                                2*(min(abs(SpawnedSign(1)), abs(SignTemp(1))))

                            if (IsUnoccDet(SpawnedSign + SignTemp)) then
                                ! If no amplitude remaining after the states are merged.
                                MinInd = PartInd
                                ToRemove = ToRemove + 1
                                cycle
                            end if
                        end if
                    end if
                end if
            end if

            IF(tSuccess) THEN

                 !Our SpawnedParts determinant is found in CurrentDets
                 !If we're using real coefficients, the "removal" step for these
                 !Determinants (and the others in currentdets) will occur later
                 !in InsertRemoveParts

!                SearchInd=PartInd   !This can actually be min(1,PartInd-1) once we know that the binary search is working, 
!                                   as we know that PartInd is the same particle.
!                MinInd=PartInd      !Make sure we only have a smaller list to search next time since the next 
!                                   particle will not be at an index smaller than PartInd
!                AnnihilateInd=0    !AnnihilateInd indicates the index in CurrentDets of the particle we want to annihilate. 
!                                   It will remain 0 if we find not complimentary particle.
!                WRITE(6,'(3I20,A,3I20)') SpawnedParts(:,i),' equals ',CurrentDets(:,PartInd)

                ! If spawning has occured from the deterministic (even partially) or to the deterministic space,
                ! then we always allow the spawning to occur. Then simply cycle.
                if (tSemiStochastic) then
                    if (test_flag(CurrentDets(:,PartInd), flag_deterministic) .or. &
                        test_flag(SpawnedParts(:, i), flag_determ_parent)) then
                        call extract_sign(CurrentDets(:, PartInd), CurrentSign)
                        call extract_sign(SpawnedParts(:, i), SpawnedSign)
                        call encode_sign(CurrentDets(:, PartInd), SpawnedSign + CurrentSign)
                        call encode_sign(SpawnedParts(:,i),null_part)
                        ToRemove = ToRemove + 1
                        MinInd = PartInd

                        ! Update stats:
                        SignProd = CurrentSign*SpawnedSign
                        if (SignProd(1) < 0.0_dp) iter_data%nannihil(1) = iter_data%nannihil(1) + &
                            2*(min(abs(CurrentSign(1)), abs(SpawnedSign(1))))

                        cycle
                    end if
                end if

                call extract_sign(CurrentDets(:,PartInd),CurrentSign)
                call extract_sign(SpawnedParts(:,i),SpawnedSign)

                SignProd=CurrentSign*SpawnedSign

!                WRITE(6,*) 'DET FOUND in list'

                ! The spawned parts contain the Dj's spawned by the Di's in CurrentDets.
                ! If the SpawnedPart is found in the CurrentDets list, it means that the Dj has a non-zero 
                ! cj - and therefore the Di.Dj pair will have a non-zero ci.cj to contribute to the RDM.
                ! The index i tells us where to look in the parent array, for the Di's to go with this Dj.
                if(tFillingStochRDMonFly.and.(.not.tHF_Ref_Explicit)) &
                    call check_fillRDM_DiDj(i,CurrentDets(:,PartInd),CurrentH(2,PartInd))

                if(sum(abs(CurrentSign)) .ne. 0) then
                    !Transfer across
                    call encode_sign(CurrentDets(:,PartInd),SpawnedSign+CurrentSign)
                    call encode_sign(SpawnedParts(:,i),null_part)

                    ! The only way SpawnedSign can be zero is if we are calculating the RDM.
                    ! If this is the case, we would have already added the SpawnedDet to Spawned_Parts_Zero
                    ! when it was compressed and all walkers were annihilated.
                    ! This only counts the walkers where the SpawnedSign has newly become zero, by merging with CurrentDets.
                    if(sum(abs(SpawnedSign)) .ne. 0) ToRemove=ToRemove+1

                    do j=1,lenof_sign   !Run over real (& imag ) components

                        if (SignProd(j) < 0) then
                            ! This indicates that the particle has found the same particle of 
                            ! opposite sign to annihilate with
                            Annihilated=Annihilated+2*(min(abs(CurrentSign(j)),abs(SpawnedSign(j))))
                            iter_data%nannihil(j) = iter_data%nannihil(j) + &
                                               2*(min(abs(CurrentSign(j)), abs(SpawnedSign(j))))

                            IF(tTruncInitiator) THEN
                                ! If we are doing an initiator calculation - then if the walkers that
                                ! are left after annihilation came from the SpawnedParts array, and 
                                ! had spawned from determinants outside the active space, then it is 
                                ! like these have been spawned on an unoccupied determinant and they 
                                ! are killed.
                                ! 
                                ! If flag_make_initiator is set, it is allowed to survive anyway, and 
                                ! is actually made into an initiator.
                                if (abs(SpawnedSign(j)) > abs(CurrentSign(j))) then
                                    if (test_flag (SpawnedParts(:,i), flag_make_initiator(j))) then
                                        call set_flag (CurrentDets(:,PartInd), flag_is_initiator(j))
                                        call set_flag (CurrentDets(:,PartInd), flag_make_initiator(j))
                                        NoAddedInitiators = NoAddedInitiators + 1
                                        if (tSpawnSpatialInit) &
                                            call add_initiator_list (CurrentDets(:,PartInd))
                                    else
                                        ! If the residual particles were spawned from non-initiator 
                                        ! particles, abort them. Encode only the correct 'type'
                                        ! of sign.
                                        if (.not. test_flag (SpawnedParts(:,i), flag_parent_initiator(j))) then
                                            NoAborted = NoAborted + abs(SpawnedSign(j)) - abs(CurrentSign(j))
                                            iter_data%naborted(j) = iter_data%naborted(j) + &
                                                                    abs(SpawnedSign(j)) - abs(CurrentSign(j))
                                            call encode_part_sign (CurrentDets(:,PartInd), 0.0_dp , j)
                                        endif
                                    endif
                                endif
                            ENDIF

                            if(tHashWalkerList) then
                                call extract_sign (CurrentDets(:,PartInd), SignTemp)
                                if (IsUnoccDet(SignTemp)) then
                                    !All walkers in this main list have been annihilated away
                                    !Remove it from the hash index array so that no others find it (it is impossible to have
                                    !another spawned walker yet to find this determinant)
                                    call RemoveDetHashIndex(nJ,PartInd)
                                    !Add to "freeslot" list so it can be filled in
                                    iEndFreeSlot=iEndFreeSlot+1
                                    FreeSlot(iEndFreeSlot)=PartInd
                                endif
                            endif
     
                            IF(tHistSpawn) THEN
!We want to histogram where the particle annihilations are taking place.
                                ExcitLevel = FindBitExcitLevel(SpawnedParts(:,i),&
                                                               iLutHF, nel)
                                IF(ExcitLevel.eq.NEl) THEN
                                    CALL BinSearchParts2(SpawnedParts(:,i),HistMinInd2(ExcitLevel),Det,PartIndex,tSuc)
                                ELSEIF(ExcitLevel.eq.0) THEN
                                    PartIndex=1
                                    tSuc=.true.
                                ELSE
                                    CALL BinSearchParts2(SpawnedParts(:,i),HistMinInd2(ExcitLevel), &
                                            FCIDetIndex(ExcitLevel+1)-1,PartIndex,tSuc)
                                ENDIF
                                HistMinInd2(ExcitLevel)=PartIndex
                                IF(tSuc) THEN
                                    AvAnnihil(j,PartIndex)=AvAnnihil(j,PartIndex)+ &
                                    REAL(2*(min(abs(CurrentSign(j)),abs(SpawnedSign(j)))))
                                    InstAnnihil(j,PartIndex)=InstAnnihil(j,PartIndex)+ &
                                    REAL(2*(min(abs(CurrentSign(j)),abs(SpawnedSign(j)))))
                                ELSE
                                    WRITE(6,*) "***",SpawnedParts(0:NIftot,i)
                                    Call WriteBitDet(6,SpawnedParts(0:NIfTot,i),.true.)
                                    CALL Stop_All("AnnihilateSpawnedParts","Cannot find corresponding FCI "&
                                        & //"determinant when histogramming")
                                ENDIF
                            ENDIF

                        else
                            ! Spawning onto an existing determinant with the same sign
                            if (tTruncInitiator) then
                                if (test_flag(SpawnedParts(:,i), flag_make_initiator(j)) .and. &
                                    .not. test_flag(CurrentDets(:,PartInd), flag_is_initiator(j))) then
                                    call set_flag (CurrentDets(:,PartInd), flag_is_initiator(j))
                                    call set_flag (CurrentDets(:,PartInd), flag_make_initiator(j))
                                    NoAddedInitiators = NoAddedInitiators + 1
                                    if (tSpawnSpatialInit) &
                                        call add_initiator_list (CurrentDets(:,PartInd))
                                endif
                            endif
                        ENDIF

                    enddo   !Finish running over components of signs
                endif
            endif
                
            if((.not.tSuccess).or.(sum(abs(CurrentSign)) .eq. 0)) then

                if (tSemiStochastic) then
                    ! If performing a semi-stochastic simulation and spawning from the deterministic
                    ! space always allow the spawning.
                    if (test_flag(SpawnedParts(:,i), flag_determ_parent)) then
                        tAlwaysAllow = .true.
                    else
                        tAlwaysAllow = .false.
                    end if
                end if

                if(tTruncInitiator .and. (.not. tAlwaysAllow)) then
                    ! Determinant in newly spawned list is not found in currentdets - usually this 
                    ! would mean the walkers just stay in this list and get merged later - but in 
                    ! this case we want to check where the walkers came from - because if the newly 
                    ! spawned walkers are from a parent outside the active space they should be 
                    ! killed - as they have been spawned on an unoccupied determinant.
                    !
                    ! If flag_make_initiator is set, then obviously these are allowed to survive
                    call extract_sign (SpawnedParts(:,i), SignTemp)

                    !If we abort these particles, we'll still need to add them to ToRemove
                    tPrevOcc=.false.
                    if (.not. IsUnoccDet(SignTemp)) tPrevOcc=.true.   
                        
                    do j = 1, lenof_sign
                        if (.not. test_flag (SpawnedParts(:,i), flag_parent_initiator(j)) .and. &
                            .not. test_flag (SpawnedParts(:,i), flag_make_initiator(j))) then
                            ! Are we allowing particles to survive if there is an
                            ! initiator with the same spatial structure?
                            ! TODO: optimise this. Only call it once?

                            ! TODO: Surely this doesn't work? Need to avoid aborting
                            !       the particle?
                            if (tSpawnSpatialInit) then
                                if (is_spatial_init(SpawnedParts(:,i))) then
                                    call set_flag (SpawnedParts(:,i), &
                                                   flag_parent_initiator(j))
                                endif
                            endif

                            ! Walkers came from outside initiator space.
                            NoAborted = NoAborted + abs(SignTemp(j))
                            iter_data%naborted(j) = iter_data%naborted(j) + abs(SignTemp(j))
                            ! We've already counted the walkers where SpawnedSign become zero in the compress,
                            ! and in the merge, all that's left is those which get aborted which are counted here
                            ! only if the sign was not already zero (when it already would have been counted).
!                            if(SignTemp(j).ne.0) ToRemove = ToRemove + 1
                            SignTemp(j) = 0.0_dp
                            call encode_part_sign (SpawnedParts(:,i), SignTemp(j), j)

                            if(tHashWalkerList.and.(sum(abs(CurrentSign)) .eq. 0)) then
                                !All walkers in this main list have died, and none have been spawned onto it.
                                !Remove it from the hash index array so that no others find it (it is impossible to have
                                !another spawned walker yet to find this determinant)
                                call RemoveDetHashIndex(nJ,PartInd)
                                !Add to "freeslot" list so it can be filled in
                                iEndFreeSlot=iEndFreeSlot+1
                                FreeSlot(iEndFreeSlot)=PartInd
                            endif
                        endif
                        if ((abs(SignTemp(j)).gt.0.0) .and. (abs(SignTemp(j)).lt.1.0)) then
                            !We remove this walker with probability 1-RealSignTemp
                            pRemove=1-abs(SignTemp(j))
                            r = genrand_real2_dSFMT ()
                            if (pRemove .gt. r) then
                                !Remove this walker
                                NoRemoved = NoRemoved + abs(SignTemp(j))
                                Annihilated = Annihilated + abs(SignTemp(j))
                                iter_data%nannihil = iter_data%nannihil + abs(SignTemp(j))
                                SignTemp(j) = 0
                                call nullify_ilut_part (SpawnedParts(:,i), j)
                            elseif (tEnhanceRemainder) then
                                NoBorn = NoBorn + 1.0_dp - abs(SignTemp(j))
                                iter_data%nborn = iter_data%nborn + 1.0_dp - abs(SignTemp(j))
                                SignTemp(j) = sign(1.0_dp, SignTemp(j))
                                call encode_part_sign (SpawnedParts(:,i), SignTemp(j), j)
                            endif
                        endif
                    enddo

                    if (IsUnoccDet(SignTemp) .and. tPrevOcc) then
                        ! All particle 'types' have been aborted
!                        The zero sign has already been taken into account in Spawned_Parts_Zero, if it was zero
!                        directly after the compress. Only add in here if not already taken care of there.
                         ToRemove = ToRemove + 1
                    elseif(tHashWalkerList) then
                        !Walkers have not been aborted, and so we should copy the determinant straight over to the main list
                        !We do not need to recompute the hash, since this should be the same one as was generated at the
                        !beginning of the loop
                        call AddNewHashDet(TotWalkersNew,SpawnedParts(:,i),DetHash,nJ)
                    endif
                else
                    ! Running Full-Scheme calc
                    ! Determinant in newly spawned list is not found in currentdets
                    ! If coeff <1, apply removal criterion
                    call extract_sign (SpawnedParts(:,i), SignTemp)
                    
                    !If we abort these particles, we'll still need to add them to ToRemove
                    tPrevOcc=.false.
                    if (.not. IsUnoccDet(SignTemp)) tPrevOcc=.true.   
                    
                    do j = 1, lenof_sign
                        if ((abs(SignTemp(j)).gt.0.0) .and. (abs(SignTemp(j)).lt.1.0)) then
                            !We remove this walker with probability 1-RealSignTemp
                            pRemove=1-abs(SignTemp(j))
                            r = genrand_real2_dSFMT ()
                            if (pRemove .gt. r) then
                                !Remove this walker
                                NoRemoved = NoRemoved + abs(SignTemp(j))
                                Annihilated = Annihilated + abs(SignTemp(j))
                                iter_data%nannihil = iter_data%nannihil + abs(SignTemp(j))
                                SignTemp(j) = 0
                                call nullify_ilut_part (SpawnedParts(:,i), j)
                            elseif (tEnhanceRemainder) then
                                NoBorn = NoBorn + 1.0_dp - abs(SignTemp(j))
                                iter_data%nborn = iter_data%nborn + 1.0_dp - abs(SignTemp(j))
                                SignTemp(j) = sign(1.0_dp, SignTemp(j))
                                call encode_part_sign (SpawnedParts(:,i), SignTemp(j), j)
                            endif
                        endif
                    enddo
                    
                    if (IsUnoccDet(SignTemp) .and. tPrevOcc) then
                        ! All particle 'types' have been aborted
!                        The zero sign has already been taken into account in Spawned_Parts_Zero, if it was zero
!                        directly after the compress. Only add in here if not already taken care of there.
                         ToRemove = ToRemove + 1
                    elseif(tHashWalkerList) then
                        !Walkers have not been aborted, and so we should copy the determinant straight over to the main list
                        !We do not need to recompute the hash, since this should be the same one as was generated at the
                        !beginning of the loop
                        call AddNewHashDet(TotWalkersNew,SpawnedParts(:,i),DetHash,nJ)
                    endif
                endif
            endif

            ! Even if a corresponding particle wasn't found, we can still
            ! search a smaller list next time....so not all bad news then...
            MinInd=PartInd

            ! We don't want the determ_parent flag to ever be set in CurrentDets, so make sure it is unset now.
            if (tSemiStochastic) call clr_flag(SpawnedParts(:,i), flag_determ_parent)

        enddo

        CALL halt_timer(BinSearch_time)

!        WRITE(6,*) "Leftover Parts..."
!        do i=1,ValidSpawned
!            WRITE(6,*) SpawnedParts(0:NIfTot-1,i),SpawnedParts(NIfTot,i)-2
!        enddo

        if(.not.tHashWalkerList) then

!Now we have to remove the annihilated particles from the spawned list. 
!They will be removed from the main list at the end of the annihilation process.
!It may actually be easier to just move the annihilated particles to the end of the list and resort the list?
!Or, the removed indices could be found on the fly? This may have little benefit though if the memory isn't needed.
            IF((ToRemove+Spawned_Parts_Zero).gt.0) THEN

!Since reading and writing from the same array is slow, copy the information
!accross to the other spawned array, and just swap the pointers around after.
                DetsMerged=0
                do i=1,ValidSpawned
!We want to move all the elements above this point down to 'fill in' the annihilated determinant.
                    call extract_sign(SpawnedParts(:,i),SignTemp)
                    IF(IsUnoccDet(SignTemp)) THEN
                        DetsMerged=DetsMerged+1
                    ELSE
                        SpawnedParts2(0:NIfTot,i-DetsMerged)=SpawnedParts(0:NIfTot,i)
                    ENDIF
                enddo
                ValidSpawned=ValidSpawned-DetsMerged
                IF(DetsMerged.ne.(ToRemove+Spawned_Parts_Zero)) THEN
                    WRITE(6,*) "***", Iter, DetsMerged, ToRemove, Spawned_Parts_Zero
                    CALL Stop_All("AnnihilateSpawnedParts","Incorrect number of particles removed from spawned list")
                ENDIF
!We always want to annihilate from the SpawedParts and SpawnedSign arrays, so swap them around.
                PointTemp => SpawnedParts2
                SpawnedParts2 => SpawnedParts
                SpawnedParts => PointTemp

            ENDIF
        else

            !Update remaining number of holes in list for walkers stats.
            if(iStartFreeSlot.gt.iEndFreeSlot) then
                !All slots filled
                HolesInList=0
!                write(6,*) "No holes in list"
            else
                HolesInList=iEndFreeSlot-(iStartFreeSlot-1)
!                write(6,*) "Holes in main list: ",HolesInList
            endif

        endif
!        call WriteExcitorListP(6,SpawnedParts,0,ValidSpawned,0,"After zero-removal")

!        WRITE(6,*) "After removal of zeros: "
!        do i=1,ValidSpawned
!            WRITE(6,*) SpawnedParts(:,i)
!        enddo

        CALL halt_timer(AnnMain_time)

    END SUBROUTINE AnnihilateSpawnedParts

    !Add a new determinant to the main list when tHashWalkerList is true
    !This involves updating the list length, copying it across, updating its flag, adding its diagonal helement(if neccessary)
    !Also need to update the hash table to point at it correctly
    subroutine AddNewHashDet(TotWalkersNew,iLutCurr,DetHash,nJ)
        implicit none
        integer, intent(inout) :: TotWalkersNew 
        integer(n_int), intent(inout) :: iLutCurr(0:NIfTot)
        integer, intent(in) :: DetHash, nJ(nel)
        integer :: DetPosition
        HElement_t :: HDiag
        character(len=*), parameter :: t_r="AddNewHashDet"

        !update its flag
        if(tTruncInitiator) call FlagifDetisInitiator(iLutCurr)

        if(iStartFreeSlot.le.iEndFreeSlot) then
            !We can slot it into a free slot in the main list, rather than increase its length
            DetPosition =  FreeSlot(iStartFreeSlot)
            CurrentDets(:,DetPosition)=iLutCurr(:)
            iStartFreeSlot=iStartFreeSlot+1
        else
            !We must increase the length of the main list to slot the new walker in
            TotWalkersNew=TotWalkersNew+1
            DetPosition=TotWalkersNew
            if(TotWalkersNew.ge.MaxWalkersPart) then
                call stop_all(t_r,"Not enough memory to merge walkers into main list. Increase MemoryFacPart")
            endif
            CurrentDets(:,DetPosition)=iLutCurr(:)
        endif
!        write(6,*) "Putting into position: ",DetPosition
        
        if(.not.tRegenDiagHEls) then
!We want to calculate the diagonal hamiltonian matrix element for the new particle to be merged.
!            call decode_bit_det (nJ, CurrentDets(:,CurrDetsIndex))
            if (tHPHF) then
                HDiag = hphf_diag_helement (nJ,CurrentDets(:,DetPosition))
            else
                HDiag = get_helement (nJ, nJ, 0)
            endif
            CurrentH(1,DetPosition)=real(HDiag,dp)-Hii
        endif

        HashIndex(HashIndex(0,DetHash),DetHash)=DetPosition
!        write(6,*) "Adding index to hash position: ",HashIndex(0,DetHash)
        HashIndex(0,DetHash)=HashIndex(0,DetHash)+1
!        write(6,*) "DetHash: ",DetHash
!        write(6,*) "HashIndex(0,DetHash): ",HashIndex(0,DetHash)
!        write(6,*) "nClashMax: ",nClashMax
        if(HashIndex(0,DetHash).gt.nClashMax) then
            call EnlargeHashTable()
        endif

    end subroutine AddNewHashDet

    subroutine EnlargeHashTable()
        implicit none
        integer :: ierr,i,FinalVal
        character(len=*), parameter :: t_r="EnlargeHashTable"

        nClashMax=nClashMax+1
        write(6,"(A)") "Enlarging hash table since it is not optimally load balanced."
        write(6,"(A,I5,A,I14,A)") "Allowing for a maximum of ",nClashMax," hash clashes, and potential storage for "&
                            ,nClashMax*nWalkerHashes," occupied determinants"

        if(allocated(HashIndexArr1)) then
!            write(6,*) "HashIndexArr1 allocated"
            allocate(HashIndexArr2(0:nClashMax,nWalkerHashes),stat=ierr)
!            write(6,*) "Ierr: ",ierr
            if(ierr.ne.0) then
                call stop_all(t_r,"Error reallocating HashIndexArr2")
            endif
!            write(6,*) "Zeroing array"
!            call flush(6)
            HashIndexArr2(:,:)=0
!            write(6,*) "Filling new hash array"
            do i=1,nWalkerHashes
!                write(6,*) i
!                call flush(6)
                FinalVal=HashIndex(0,i)-1
                HashIndexArr2(0:FinalVal,i) = HashIndex(0:FinalVal,i)
            enddo
!            write(6,*) "Deallocating old hash array"
            deallocate(HashIndexArr1)
            nullify(HashIndex)
            HashIndex => HashIndexArr2
        else
!            write(6,*) "HashIndexArr2 allocated"
            if(.not.allocated(HashIndexArr2)) then
                call stop_all(t_r,"HashIndexArr2 should be allocated")
            endif
            if(allocated(HashIndexArr1)) then
                call stop_all(t_r,"HashIndexArr1 should not be allocated")
            endif
            allocate(HashIndexArr1(0:nClashMax,nWalkerHashes),stat=ierr)
            if(ierr.ne.0) then
                call stop_all(t_r,"Error reallocating HashIndexArr2")
            endif
            HashIndexArr1(:,:)=0
            do i=1,nWalkerHashes
                FinalVal=HashIndex(0,i)-1
                HashIndexArr1(0:FinalVal,i)=HashIndex(0:FinalVal,i)
            enddo
            deallocate(HashIndexArr2)
            nullify(HashIndex)
            HashIndex => HashIndexArr1
        endif

    end subroutine EnlargeHashTable

    SUBROUTINE CalcHashTableStats(TotWalkersNew)
        use util_mod, only: abs_sign
        use CalcData , only : tCheckHighestPop
        integer, intent(in) :: TotWalkersNew
        INTEGER :: i,part_type,AnnihilatedDet
        real(dp) :: CurrentSign(lenof_sign), SpawnedSign(lenof_sign)
        character(*), parameter :: t_r = 'CalcHashTableStats'

        if(.not.bNodeRoot) return
        TotParts=0.0
        norm_psi_squared = 0.0
        iHighestPop=0
        AnnihilatedDet=0
        IF(TotWalkersNew.gt.0) THEN
            do i=1,TotWalkersNew
                call extract_sign(CurrentDets(:,i),CurrentSign)
                IF(IsUnoccDet(CurrentSign)) then
                    AnnihilatedDet=AnnihilatedDet+1 
                    IF(tTruncInitiator) THEN
                        do part_type=1,lenof_sign
                            if (test_flag(CurrentDets(:,i),flag_parent_initiator(part_type))) then
                                !determinant was an initiator...it obviously isn't any more...
                                NoAddedInitiators=NoAddedInitiators-1
                                if (tSpawnSpatialInit) &
                                    call rm_initiator_list (CurrentDets(:,i))
                            endif
                        enddo
                    ENDIF
                ELSE
!We want to move all the elements above this point down to 'fill in' the annihilated determinant.
                    TotParts=TotParts+abs(CurrentSign)
                    norm_psi_squared = norm_psi_squared + sum(CurrentSign**2)
                    IF(tCheckHighestPop) THEN
                        ! If this option is on, then we want to compare the 
                        ! weight on each determinant to the weight at the HF 
                        ! determinant.
                        !
                        ! Record the highest weighted determinant on each 
                        ! processor.
                        if (abs_sign(ceiling(CurrentSign)) > iHighestPop) then
                            iHighestPop = abs_sign(ceiling(CurrentSign))
                            HighestPopDet(:)=CurrentDets(:,i)
                        end if
                    ENDIF
                ENDIF
            enddo
        ENDIF
        if(AnnihilatedDet.ne.HolesInList) then
            write(6,*) "TotWalkersNew: ",TotWalkersNew
            write(6,*) "AnnihilatedDet: ",AnnihilatedDet
            write(6,*) "HolesInList: ",HolesInList
            call stop_all(t_r,"Error in determining annihilated determinants")
        endif

    END SUBROUTINE CalcHashTableStats
    
    !routine to find and remove the index to a determinant from the HashIndex array
    !This could potentially be speeded up by an ordered HashIndex array, and binary searching
    subroutine RemoveDetHashIndex(nI,DetPosition)
        implicit none
!        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: nI(nel) 
        integer, intent(in) :: DetPosition
        integer :: DetHash,FinalVal,i
        character(*), parameter :: this_routine="RemoveDetHashIndex"

!        write(6,*) "Sending to hash func 2: ",nI(:)
        DetHash=FindWalkerHash(nI)
        FinalVal=HashIndex(0,DetHash)-1
        do i=1,FinalVal !Run through all indices
!            write(6,*) "i, HashIndex(i,DetHash): ",i, HashIndex(i,DetHash)
            if(HashIndex(i,DetHash).eq.DetPosition) exit
        enddo
!        write(6,*) "DetHash: ",DetHash
!        write(6,*) "FinalVal: ",FinalVal
!        write(6,*) "DetPosition: ",DetPosition
        ASSERT(i.le.FinalVal)

        !i is now the position of the determinant to remove
        !Indicate that there is now one less determinant at this hash value
        HashIndex(0,DetHash)=HashIndex(0,DetHash)-1
        !Move everything from i+1 -> HashIndex(0,DetHash) down one
        if(i.ne.FinalVal) then
            HashIndex(i:FinalVal-1,DetHash)=HashIndex(i+1:FinalVal,DetHash)
        endif

    end subroutine RemoveDetHashIndex

    
!This routine will run through the total list of particles (TotWalkersNew in CurrentDets with sign CurrentSign) and the list of 
!newly-spawned but non annihilated particles (ValidSpawned in SpawnedParts and SpawnedSign) and move the new particles into the
!correct place in the new list, while removing the particles with sign = 0 from CurrentDets. 
!Binary searching can be used to speed up this transfer substantially.
!The key feature which makes this work, is that it is impossible for the same determinant to be specified in both the spawned 
!and main list at the end of the annihilation process. Therefore we won't multiply specify determinants when we merge the lists.
    SUBROUTINE InsertRemoveParts(ValidSpawned, TotWalkersNew, iter_data)
        use util_mod, only: abs_sign
        use SystemData, only: tHPHF, tRef_Not_HF
        use bit_reps, only: NIfD
        use CalcData , only : tCheckHighestPop, NMCyc, InitiatorWalkNo
        use Logging , only : tRDMonFly, tExplicitAllRDM, tInitiatorRDM
        use nElRDMMod , only : det_removed_fill_diag_rdm 
        INTEGER, intent(in) :: ValidSpawned
        integer, intent(inout) :: TotWalkersNew
        real(dp) :: CurrentSign(lenof_sign), SpawnedSign(lenof_sign)
        real(dp) :: HDiag, pRemove, r
        INTEGER :: i,DetsMerged,nJ(NEl),part_type, ExcitLevelCurr, j
        LOGICAL :: TestClosedShellDet, tIsStateDeterm
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: this_routine = 'InsertRemoveParts'
        HElement_t :: HDiagTemp

!It appears that the rest of this routine isn't thread-safe if ValidSpawned is zero.
        
    ! This logical is only used for the semi-stochastic code. If true then we don't
    ! try to remove the state.
    tIsStateDeterm = .false.

    if(.not.bNodeRoot) return
!Annihilated determinants first are removed from the main array (zero sign). 
!Surely we only need to perform this loop if the number of annihilated particles > 0?
        TotParts=0.0
        norm_psi_squared = 0.0
        DetsMerged=0
        iHighestPop=0
        InstNoatHF = 0.0
        IF(TotWalkersNew.gt.0) THEN
            do i=1,TotWalkersNew
                call extract_sign(CurrentDets(:,i),CurrentSign)
                if (tSemiStochastic) tIsStateDeterm = test_flag(CurrentDets(:,i), flag_deterministic)
                do j=1, lenof_sign
                    if (.not. tIsStateDeterm) then
                        if ((abs(CurrentSign(j)).gt.0.0) .and. (abs(CurrentSign(j)).lt.1.0)) then
                            !We remove this walker with probability 1-RealSignTemp
                            pRemove=1-abs(CurrentSign(j))
                            r = genrand_real2_dSFMT ()
                            if (pRemove .gt. r) then
                                !Remove this walker
                                NoRemoved = NoRemoved + abs(CurrentSign(j))
                                Annihilated = Annihilated + abs(CurrentSign(j))
                                iter_data%nannihil = iter_data%nannihil + abs(CurrentSign(j))
                                CurrentSign(j) = 0
                                call nullify_ilut_part (CurrentDets(:,i), j)
                            elseif (tEnhanceRemainder) then
                                ! SDS: TODO: Account for the TotParts Changes
                                ! Should we always do this here? Probably. Should

                                NoBorn = NoBorn + 1.0_dp - abs(CurrentSign(j))
                                iter_data%nborn = iter_data%nborn + 1.0_dp - abs(CurrentSign(j))
                                CurrentSign(j) = sign(1.0_dp, CurrentSign(j))
                                call encode_part_sign (CurrentDets(:,i), CurrentSign(j), j)
                            endif
                        endif
                    endif
                enddo
                
                if(i.eq.HFInd)  InstNoatHF = CurrentSign

                IF(IsUnoccDet(CurrentSign) .and. (.not. tIsStateDeterm)) THEN
                    DetsMerged=DetsMerged+1
                    if(tFillingStochRDMonFly) &
                        call det_removed_fill_diag_rdm(CurrentDets(:,i), CurrentH(1:NCurrH,i))
                    if(i.eq.HFInd) then
                        !We have to do this such that AvNoAtHF matches up with AvSign.
                        !AvSign is extracted from CurrentH, and if the HFDet is unoccupied
                        !at this moment during annihilation, it's CurrentH entry is removed
                        !and the averaging information in it is lost.
                        !In some cases (a successful spawning event) a CurrentH entry will
                        !be recreated, but with AvSign 0, so we must match this here.
                        AvNoAtHF=0.0_dp 
                        IterRDM_HF = Iter + 1 
                    endif
                    IF(tTruncInitiator) THEN
                        do part_type=1,lenof_sign
                            if (test_flag(CurrentDets(:,i),flag_parent_initiator(part_type))) then
                                !determinant was an initiator...it obviously isn't any more...
                                NoAddedInitiators=NoAddedInitiators-1
                                if (tSpawnSpatialInit) &
                                    call rm_initiator_list (CurrentDets(:,i))
                            endif
                        enddo
                    ENDIF
                ELSE
!We want to move all the elements above this point down to 'fill in' the annihilated determinant.
                    IF(DetsMerged.ne.0) THEN
                        CurrentDets(0:NIfTot,i-DetsMerged)=CurrentDets(0:NIfTot,i)
                        IF(.not.tRegenDiagHEls) &
                            CurrentH(:,i-DetsMerged)=CurrentH(:,i)
                    ENDIF
                    TotParts=TotParts+abs(CurrentSign)
                    ! CMO
                    norm_psi_squared = norm_psi_squared + sum(CurrentSign**2)
                    IF(tCheckHighestPop) THEN
!If this option is on, then we want to compare the weight on each determinant to the weight at the HF determinant.
!Record the highest weighted determinant on each processor.
                        if (abs_sign(ceiling(CurrentSign)) > iHighestPop) then
                            iHighestPop = abs_sign(ceiling(CurrentSign))
                            HighestPopDet(:)=CurrentDets(:,i)
                        endif
                    ENDIF
                ENDIF
            enddo
            TotWalkersNew=TotWalkersNew-DetsMerged
        ENDIF

!        do i=1,TotWalkersNew
!            IF(CurrentSign(i).eq.0) THEN
!                CALL Stop_All("InsertRemoveParts","Particles not removed correctly")
!            ENDIF
!        enddo
!        CALL CheckOrdering(CurrentDets,CurrentSign(1:TotWalkersNew),TotWalkersNew,.true.)

!We now calculate the contribution to the total number of particles from the spawned lists.
!The list has previously been compressed.
        NumMerged=ValidSpawned        

        IF(ValidSpawned.gt.0) THEN
            call extract_sign(SpawnedParts(:,1),SpawnedSign)
            TotParts=TotParts+abs(SpawnedSign)
            norm_psi_squared = norm_psi_squared + sum(SpawnedSign**2)
        ENDIF
        do i=2,ValidSpawned
            call extract_sign(SpawnedParts(:,i),SpawnedSign)
            TotParts=TotParts+abs(SpawnedSign)
            norm_psi_squared = norm_psi_squared + sum(SpawnedSign**2)
        enddo

!        CALL CheckOrdering(SpawnedParts,SpawnedSign(1:ValidSpawned),ValidSpawned,.true.)
!
!        WRITE(6,*) "Before merging...",TotWalkersNew,Iter
!        do i=1,TotWalkersNew
!            WRITE(6,*) i,CurrentDets(:,i),CurrentSign(i)
!        enddo
!        WRITE(6,*) "***"
!        do i=1,ValidSpawned
!            WRITE(6,*) i,SpawnedParts(:,i),SpawnedSign(i)
!        enddo
!        WRITE(6,*) "***"
!        CALL neci_flush(6)

!TotWalkersNew is now the number of determinants in the main list left.
!We now want to merge the main list with the spawned list of non-annihilated spawned particles.
!The final list will be of length TotWalkersNew+ValidSpawned. This will be returned in the first element of MergeLists updated.
        if(TotWalkersNew+ValidSpawned>MaxWalkersPart) then
            WRITE(6,*) "Non-annihilated old walkers:",TotWalkersNew
            WRITE(6,*) "Non-annihilated spawned:",ValidSpawned
            WRITE(6,*) "Total walkers to remain:", TotWalkersNew+ValidSpawned
            WRITE(6,*) "Size of Particle List:", MaxWalkersPart
            call stop_all(this_routine, "Not enough space in particle list for merge.") 
        endif
        IF(tRegenDiagHEls) THEN
            IF(TotWalkersNew.eq.0) THEN
!Merging algorithm will not work with no determinants in the main list.
                TotWalkersNew=ValidSpawned
                do i=1,ValidSpawned
                    CurrentDets(:,i)=SpawnedParts(:,i)
                    IF(tTruncInitiator) CALL FlagifDetisInitiator(CurrentDets(0:NIfTot,i))
                enddo
            ELSE
                CALL MergeLists(TotWalkersNew,ValidSpawned,SpawnedParts(0:NIfTot,1:ValidSpawned))
            ENDIF
        ELSE
            IF(TotWalkersNew.eq.0) THEN
!Merging algorithm will not work with no determinants in the main list.
                TotWalkersNew=ValidSpawned
                do i=1,ValidSpawned
                    CurrentDets(:,i)=SpawnedParts(:,i)
                    IF(tTruncInitiator) CALL FlagifDetisInitiator(CurrentDets(0:NIfTot,i))
!We want to calculate the diagonal hamiltonian matrix element for the new particle to be merged.
                    if (DetBitEQ(CurrentDets(:,i), iLutRef, NIfDBO)) then
!We know we are at HF - HDiag=0
                        HDiag=0.D0
                        call extract_sign(CurrentDets(:,i),InstNoatHF)
                    else
                        call decode_bit_det (nJ, CurrentDets(:,i))
                        if (tHPHF) then
                            HDiagTemp = hphf_diag_helement (nJ, &
                                                            CurrentDets(:,i))
                        else
                            HDiagTemp = get_helement (nJ, nJ, 0)
                        endif
                        HDiag=(REAL(HDiagTemp,dp))-Hii
                    endif
                    CurrentH(1,i)=HDiag
                enddo
            ELSE
                CALL MergeListswH(TotWalkersNew,ValidSpawned,SpawnedParts(0:NIfTot,1:ValidSpawned))
            ENDIF

        ENDIF
!        CALL CheckOrdering(CurrentDets,CurrentSign(1:TotWalkers),TotWalkers,.true.)

    END SUBROUTINE InsertRemoveParts

    !Another hashing routine - this time to find the correct position in the hash table
    pure function FindWalkerHash(nJ) result(hashInd)
        implicit none
        integer, intent(in) :: nJ(nel)
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
        hashInd = int(abs(mod(hash, int(nWalkerHashes, int64))),sizeof_int)+1
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

    SUBROUTINE LinSearchParts(DetArray,iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: MinInd,MaxInd,PartInd
        INTEGER(KIND=n_int) :: iLut(0:NIfTot),DetArray(0:NIfTot,1:MaxInd)
        INTEGER :: N,Comp
        LOGICAL :: tSuccess

        N=MinInd
        do while(N.le.MaxInd)
            Comp=DetBitLT(DetArray(:,N),iLut(:),NIfDBO)
            IF(Comp.eq.1) THEN
                N=N+1
            ELSEIF(Comp.eq.-1) THEN
                PartInd=N-1
                tSuccess=.false.
                RETURN
            ELSE
                tSuccess=.true.
                PartInd=N
                RETURN
            ENDIF
        enddo
        tSuccess=.false.
        PartInd=MaxInd-1

    END SUBROUTINE LinSearchParts

!Do a binary search in CurrentDets, between the indices of MinInd and MaxInd. If successful, tSuccess will be true and 
!PartInd will be a coincident determinant. If there are multiple values, the chosen one may be any of them...
!If failure, then the index will be one less than the index that the particle would be in if it was present in the list.
!(or close enough!)
    SUBROUTINE BinSearchParts(iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER(KIND=n_int) :: iLut(0:NIfTot)
        INTEGER :: MinInd,MaxInd,PartInd
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(CurrentDets(:,MaxInd),iLut(:),NIfDBO,.false.)
            IF(Comp.eq.0) THEN
                tSuccess=.true.
                PartInd=MaxInd
                RETURN
            ELSE
                tSuccess=.false.
                PartInd=MinInd
            ENDIF
        ENDIF

        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(CurrentDets(:,N),iLut(:),NIfDBO,.false.)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower
!bound of the search up to N.  However, if the lower bound is already equal to N then the two bounds are 
!consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry 
!of the list and leave.  We need to check the last index.
                    Comp=DetBitLT(CurrentDets(:,i+1),iLut(:),NIfDBO,.false.)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN

                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchParts
    
END MODULE AnnihilationMod


