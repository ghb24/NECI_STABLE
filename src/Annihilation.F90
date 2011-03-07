!This module is to be used for various types of walker MC annihilation in serial and parallel.
MODULE AnnihilationMod
    use SystemData , only : NEl, tHPHF
    use CalcData , only : TRegenExcitgens,tRegenDiagHEls
    USE DetCalcData , only : Det,FCIDetIndex
    USE Parallel
    USE dSFMT_interface , only : genrand_real2_dSFMT
    USE FciMCData
    use DetBitOps, only: DetBitEQ, DetBitLT, FindBitExcitLevel, ilut_lt, &
                         ilut_gt
    use spatial_initiator, only: add_initiator_list, rm_initiator_list, &
                                 is_spatial_init
    use CalcData , only : tTruncInitiator, tSpawnSpatialInit
    use Determinants, only: get_helement
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use sort_mod
    use constants, only: n_int,lenof_sign,null_part
    use bit_rep_data
    use bit_reps, only: decode_bit_det, extract_sign, extract_flags, &
                        encode_sign, encode_flags, test_flag, set_flag, &
                        clr_flag, flag_parent_initiator, encode_part_sign, &
                        extract_part_sign, copy_flag
    use csf_data, only: csf_orbital_mask
    use hist_data, only: tHistSpawn, HistMinInd2
    IMPLICIT NONE

    contains

    !TODO:
    !   H Elements - send through logical to decide whether to create or not.
    !   Parallel spawned parts - create the ValidSpawnedList itself.
    !   Going to have to sort this out for the new packaged walkers - will have to package them up in this interface.
    SUBROUTINE AnnihilationInterface(TotDets,MainParts,MaxMainInd,SpawnDets,SpawnParts,MaxSpawnInd,iter_data)
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
!                           does *not* need to be ordered or sign coherent, and can also contain 'zero' sign particles.  Each particle contains its own sign
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
        integer(kind=n_int), pointer,save :: SpawnVecLocal(:,:)
        Annihil_time%timer_name='Annihilation interface'
        call set_timer(Annihil_time,20)

        IF(.not.(ALLOCATED(ValidSpawnedList))) THEN
!This needs to be filled correctly before annihilation can take place.
            ALLOCATE(ValidSpawnedList(0:nNodes-1),stat=ierr)
        ENDIF
        call MPIBarrier(ierr)
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
!        if(iProcIndex==root) then        
!Signs put back again into seperate array
!           do i=1,TotDets
!               call extract_sign(CurrentDets(:,i),TempSign)
!               MainSign(i)=TempSign(1)
!           enddo
!         endif
         call MPIBarrier(ierr)

        call halt_timer(Annihil_time)
    END SUBROUTINE AnnihilationInterface


!This is a new annihilation algorithm. In this, determinants are kept on predefined processors, and newlyspawned particles are sent here so that all the annihilations are
!done on a predetermined processor, and not rotated around all of them.
    SUBROUTINE DirectAnnihilation(TotWalkersNew, iter_data, tSingleProc)
        use bit_reps, only: test_flag
        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        INTEGER :: MaxIndex,ierr
        INTEGER(Kind=n_int) , POINTER :: PointTemp(:,:)
        logical, intent(in) :: tSingleProc
         

!        WRITE(6,*) "Direct annihilation"
!        CALL FLUSH(6)

!This routine will send all the newly-spawned particles to their correct processor. MaxIndex is returned as the new number of newly-spawned particles on the processor. May have duplicates.
!The particles are now stored in SpawnedParts2/SpawnedSign2.
!        call WriteExcitorListP2(6,SpawnedParts,InitialSpawnedSlots,ValidSpawnedList,0,"Local")
!        if(bNodeRoot) 
        CALL SendProcNewParts(MaxIndex,tSingleProc)

!        WRITE(6,*) "Sent particles"
!        WRITE(6,*) 'MaxIndex',MaxIndex
!        CALL FLUSH(6)

!CompressSpawnedList works on SpawnedParts arrays, so swap the pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp

!Now we want to order and compress the spawned list of particles. This will also annihilate the newly spawned particles amongst themselves.
!MaxIndex will change to reflect the final number of unique determinants in the newly-spawned list, and the particles will end up in the spawnedSign/SpawnedParts lists.
!        WRITE(6,*) "Transferred",MaxIndex

        CALL CompressSpawnedList(MaxIndex, iter_data)  

!        WRITE(6,*) "List compressed",MaxIndex,TotWalkersNew

!Binary search the main list and copy accross/annihilate determinants which are found.
!This will also remove the found determinants from the spawnedparts lists.

        CALL AnnihilateSpawnedParts(MaxIndex,TotWalkersNew, iter_data)  

 !       WRITE(6,*) "Annihilation finished",MaxIndex,TotWalkersNew

!Put the surviving particles in the main list, maintaining order of the main list.
!Now we insert the remaining newly-spawned particles back into the original list (keeping it sorted), and remove the annihilated particles from the main list.
        CALL set_timer(Sort_Time,30)
        CALL InsertRemoveParts(MaxIndex,TotWalkersNew)

        CALL halt_timer(Sort_Time)
        call MPIBarrier(ierr)

    END SUBROUTINE DirectAnnihilation

!This routine is used for sending the determinants to the correct processors. 
    SUBROUTINE SendProcNewParts(MaxIndex,tSingleProc)
        use constants, only: MpiDetInt
        REAL :: Gap
        INTEGER :: i,sendcounts(nProcessors),disps(nProcessors),recvcounts(nProcessors),recvdisps(nProcessors),error
        INTEGER :: MaxSendIndex
        INTEGER, INTENT(OUT) :: MaxIndex
        LOGICAL, intent(in) :: tSingleProc
      

        if(tSingleProc) then
!Put all particles and gap on one proc.

!ValidSpawnedList(0:nNodes-1) indicates the next free index for each processor (for spawnees from this processor)
!  i.e. the list of spawned particles has already been arranged so that newly spawned particles are grouped according to the processor they go to.

! sendcounts(1:) indicates the number of spawnees to send to each processor
! disps(1:) is the index into the spawned list of the beginning of the list to send to each processor (0-based)
           sendcounts(1)=ValidSpawnedList(0)-1
           disps(1)=0
           if(nNodes>1) then
              sendcounts(2:nNodes)=0
              disps(2:nNodes)=ValidSpawnedList(1)
           endif
        else
!Distribute the gaps on all procs
           do i=0,nProcessors-1
               if(NodeRoots(ProcNode(i))==i) then  !This is a root 
                  sendcounts(i+1)=ValidSpawnedList(ProcNode(i))-InitialSpawnedSlots(ProcNode(i))
! disps is zero-based, but InitialSpawnedSlots is 1-based
                  disps(i+1)=InitialSpawnedSlots(ProcNode(i))-1
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
        
        CALL set_timer(Comms_Time,30)

        CALL MPIAlltoAll(sendcounts,1,recvcounts,1,error)

!We can now get recvdisps from recvcounts, since we want the data to be contiguous after the move.
        recvdisps(1)=0
        do i=2,nProcessors
            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
        enddo
        MaxIndex=recvdisps(nProcessors)+recvcounts(nProcessors)
        do i=1,nProcessors
            recvdisps(i)=recvdisps(i)*(NIfTot+1)
            recvcounts(i)=recvcounts(i)*(NIfTot+1)
            sendcounts(i)=sendcounts(i)*(NIfTot+1)
            disps(i)=disps(i)*(NIfTot+1)
        enddo

!Max index is the largest occupied index in the array of hashes to be ordered in each processor 
        IF(MaxIndex.gt.(0.9*MaxSpawned)) THEN
            write(6,*) MaxIndex,MaxSpawned
            CALL Warning("SendProcNewParts","Maximum index of newly-spawned array is " &
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
!        call flush(6)
        CALL MPIAlltoAllv(SpawnedParts,sendcounts,disps,SpawnedParts2,recvcounts,recvdisps,error)
!        write(6,*) "Second all to all finish"
!        call flush(6)

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
        integer :: EndBlockDet, part_type, StartCycleInit, cum_count
        INTEGER, DIMENSION(lenof_sign) :: SpawnedSign,Temp_Sign
        LOGICAL :: tSuc, tInc
        INTEGER(Kind=n_int) , POINTER :: PointTemp(:,:)
        integer(n_int) :: cum_det (0:niftot)
        CHARACTER(len=*), parameter :: this_routine='CompressSpawnedList'

!We want to sort the list of newly spawned particles, in order for quicker binary searching later on. (this is not essential, but should proove faster)
!They should remain sorted after annihilation between spawned
        if(.not.bNodeRoot) return
!        call WriteExcitorListP(6,SpawnedParts,0,ValidSpawned,0,"UnOrdered")
        call sort(SpawnedParts(:,1:ValidSpawned), ilut_lt, ilut_gt)
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
        do while(BeginningBlockDet.le.ValidSpawned)
            !loop in blocks of the same determinant to the end of the list of walkers

            FirstInitIndex=0
            CurrentBlockDet=BeginningBlockDet+1

            do while(CurrentBlockDet.le.ValidSpawned)
                if(.not.(DetBitEQ(SpawnedParts(:,BeginningBlockDet),SpawnedParts(:,CurrentBlockDet),NIfDBO))) exit
                !loop over walkers on the same determinant in SpawnedParts
                CurrentBlockDet=CurrentBlockDet+1
            enddo

            EndBlockDet=CurrentBlockDet-1 !EndBlockDet indicates that we have reached the end of the block of similar dets
!            WRITE(6,*) "Found Block: ",BeginningBlockDet," -> ",EndBlockDet

            if(EndBlockDet.eq.BeginningBlockDet) then
                !Optimisation: This block only consists of one entry. Simply copy it across rather than 
                !               explicitly searching the list.
                SpawnedParts2(:,VecInd)=SpawnedParts(:,BeginningBlockDet)   !Transfer all info to the other array
                VecInd=VecInd+1
                BeginningBlockDet=CurrentBlockDet           !Move onto the next block of determinants
                CYCLE   !Skip the rest of this block
            endif

            ! Reset the cumulative determinant
            cum_det = 0
            cum_det (0:nifdbo) = SpawnedParts(0:nifdbo, BeginningBlockDet)
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
                            call FindResidualParticle (cum_det, SpawnedParts(:,i), cum_count, part_type, iter_data)
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
                        call FindResidualParticle (cum_det, SpawnedParts(:,i), cum_count, part_type, iter_data)
                    endif
                enddo

            enddo ! End loop over particle type

            ! Copy details into the final array
            call extract_sign (cum_det, temp_sign)
            if (sum(abs(temp_sign)) > 0) then
                ! Transfer all ino into the other array.
                SpawnedParts2(:,VecInd) = cum_det
                VecInd = VecInd + 1
                DetsMerged = DetsMerged + EndBlockDet - BeginningBlockDet
            else
                ! All particles from block have been annihilated.
                DetsMerged = DetsMerged + EndBlockDet - BeginningBlockDet + 1
            endif

            BeginningBlockDet=CurrentBlockDet           !Move onto the next block of determinants

        enddo   

        ValidSpawned=ValidSpawned-DetsMerged    !This is the new number of unique spawned determinants on the processor
        IF(ValidSpawned.ne.(VecInd-1)) THEN
            CALL Stop_All(this_routine,"Error in compression of spawned list")
        ENDIF

!Want the compressed list in spawnedparts at the end of it - swap pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp

!        WRITE(6,*) "************************"
!        call WriteExcitorListP(6,SpawnedParts,0,ValidSpawned,0,"Compressed Spawned")


!!        WRITE(6,*) "Compressed List: ",ValidSpawned,DetsMerged
!!        do i=1,ValidSpawned
!!            WRITE(6,*) SpawnedParts(0:NIfTot-1,i),SpawnedParts(NIfTot,i)-2
!!        enddo
!!        WRITE(6,*) "************************"
!!        WRITE(6,*) "Compressed List: ",ValidSpawned,DetsMerged
!!        do i=1,ValidSpawned
!!            WRITE(6,*) SpawnedParts(0:NIfTot-1,i),SpawnedParts(NIfTot,i)-2
!!        enddo
!!        WRITE(6,*) "***","iter_data%naborted: ",iter_data%naborted
!!        CALL FLUSH(6)
        
    END SUBROUTINE CompressSpawnedList

!Histogram a possible annihilation event
    subroutine HistAnnihilEvent(iLut,Sign1,Sign2,part_type)
        implicit none
        integer(kind=n_int), intent(in) :: iLut(0:NIfTot)
        integer, dimension(lenof_sign), intent(in) :: Sign1,Sign2
        integer, intent(in) :: part_type
        integer :: ExcitLevel,PartIndex
        logical :: tSuc

!We want to histogram where the particle annihilations are taking place.
        if((Sign1(part_type)*Sign2(part_type)).ge.0) return   !No annihilation occuring - particles same sign

        ExcitLevel = FindBitExcitLevel(iLut,iLutHF, nel)
        IF(ExcitLevel.eq.NEl) THEN
            CALL BinSearchParts2(iLut(:),HistMinInd2(ExcitLevel),Det,PartIndex,tSuc)
        ELSEIF(ExcitLevel.eq.0) THEN
            PartIndex=1
            tSuc=.true.
        ELSE
            CALL BinSearchParts2(iLut(:),HistMinInd2(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartIndex,tSuc)
        ENDIF
        HistMinInd2(ExcitLevel)=PartIndex
        IF(tSuc) THEN
            AvAnnihil(part_type,PartIndex)=AvAnnihil(part_type,PartIndex)+ &
                REAL(2*(MIN(abs(Sign1(part_type)),abs(Sign2(part_type)))),dp)
            InstAnnihil(part_type,PartIndex)=InstAnnihil(part_type,PartIndex)+ &
                REAL(2*(MIN(abs(Sign1(part_type)),abs(Sign2(part_type)))),dp)
        ELSE
            CALL Stop_All("CompressSpawnedList","Cannot find corresponding FCI determinant when histogramming")
        ENDIF

    end subroutine HistAnnihilEvent

    !This routine is called when compressing the spawned list.
    !It takes the sign and flags from two particles on the same determinant, 
    !and calculates what the residual particles and signs should be from them.
    !This deals with real and imaginary signs seperately, and so the 'signs' are integers.

    subroutine FindResidualParticle (cum_det, new_det, cum_count, part_type, &
                                     iter_data)

        ! This routine is called whilst compressing the spawned list during
        ! annihilation. It considers the sign and flags from two particles
        ! on the same determinant, and calculates the correct sign/flags
        ! for the compressed particle.
        !
        ! --> The information is stored within the first particle in a block
        ! --> Should be called for real/imaginary particles seperately

        integer(n_int), intent(inout) :: cum_det(0:nIfTot)
        integer(n_int), intent(in) :: new_det(0:niftot)
        integer, intent(inout) :: cum_count
        integer, intent(in) :: part_type
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: new_sgn, cum_sgn, sgn_prod

        ! Obtain the signs and sign product. Ignore new particel if zero.
        new_sgn = extract_part_sign (new_det, part_type)
        if (new_sgn == 0) return
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
            if (sgn_prod > 0) then
                ! Signs are the same. This must be an initiator.
                ! Equivalent to (deprecated) tKeepDoubSpawns
                call set_flag (cum_det, flag_is_initiator(part_type))
            elseif (sgn_prod < 0) then
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

        ! Update annihilation statistics (is this really necessary?)
        if (sgn_prod < 0) then
            Annihilated = Annihilated + 2*min(abs(cum_sgn), abs(new_sgn))
            iter_data%nannihil(part_type) = iter_data%nannihil(part_type)&
                + 2 * min(abs(cum_sgn), abs(new_sgn))
        endif

        ! Update the cumulative sign count
        call encode_part_sign (cum_det, cum_sgn + new_sgn, part_type)

    end subroutine FindResidualParticle

    
!In this routine, we want to search through the list of spawned particles. For each spawned particle, 
!we binary search the list of particles on the processor
!to see if an annihilation event can occur. The annihilated particles are then removed from the spawned list
!to the whole list of spawned particles at the end of the routine.
!In the main list, we change the 'sign' element of the array to zero. These will be deleted at the end of the total annihilation step.
    SUBROUTINE AnnihilateSpawnedParts(ValidSpawned,TotWalkersNew, iter_data)
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, intent(in) :: TotWalkersNew
        integer, intent(inout) :: ValidSpawned 
        INTEGER :: MinInd,PartInd,i,j,ToRemove,DetsMerged,PartIndex
        INTEGER, DIMENSION(lenof_sign) :: SignProd,CurrentSign,SpawnedSign,SignTemp
        INTEGER :: ExcitLevel
        INTEGER(KIND=n_int) , POINTER :: PointTemp(:,:)
        LOGICAL :: tSuccess,tSuc

        if(.not.bNodeRoot) return  !Only node roots to do this.

        CALL set_timer(AnnMain_time,30)

!MinInd indicates the minimum bound of the main array in which the particle can be found.
!Since the spawnedparts arrays are ordered in the same fashion as the main array, we can find the particle position in the main array by only searching a subset.
        MinInd=1
        IF(tHistSpawn) HistMinInd2(1:NEl)=FCIDetIndex(1:NEl)
        ToRemove=0  !The number of particles to annihilate
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

        do i=1,ValidSpawned

!This will binary search the CurrentDets array to find the desired particle. tSuccess will determine whether the particle has been found or not.
!It will also return the index of the position one below where the particle would be found if was in the list.
!            CALL LinSearchParts(CurrentDets(:,1:TotWalkersNew),SpawnedParts(0:NIfD,i),MinInd,TotWalkersNew,PartInd,tSuccess)
!            WRITE(6,*) "MinInd",MinInd
!            do j=MinInd,min(MinInd+10,TotWalkersNew)
!               write(6,"(I5)",advance='no') j 
!               call WriteBitDet(6,CurrentDets(:,j),.true.)
!            enddo
            CALL BinSearchParts(SpawnedParts(:,i),MinInd,TotWalkersNew,PartInd,tSuccess)
!            WRITE(6,"(A,2I6,L)",advance="no") "Binary search complete: ",i,PartInd,tSuccess
!            call WriteBitDet(6,SpawnedParts(:,i),.true.)
!            call WriteBitDet(6,CurrentDets(:,PartInd),.true.)

            IF(tSuccess) THEN
!                SearchInd=PartInd   !This can actually be min(1,PartInd-1) once we know that the binary search is working, as we know that PartInd is the same particle.
!                MinInd=PartInd      !Make sure we only have a smaller list to search next time since the next particle will not be at an index smaller than PartInd
!                AnnihilateInd=0     !AnnihilateInd indicates the index in CurrentDets of the particle we want to annihilate. It will remain 0 if we find not complimentary particle.
!                WRITE(6,'(3I20,A,3I20)') SpawnedParts(:,i),' equals ',CurrentDets(:,PartInd)
                call extract_sign(CurrentDets(:,PartInd),CurrentSign)
                call extract_sign(SpawnedParts(:,i),SpawnedSign)
                SignProd=CurrentSign*SpawnedSign

                !Transfer across
                call encode_sign(CurrentDets(:,PartInd),SpawnedSign+CurrentSign)
                call encode_sign(SpawnedParts(:,i),null_part)
                ToRemove=ToRemove+1

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
                                        iter_data%naborted(j) = iter_data%naborted(j) + abs(SpawnedSign(j)) - abs(CurrentSign(j))
                                        call encode_part_sign (CurrentDets(:,PartInd), 0, j)
                                    endif
                                endif
                            endif
                        ENDIF
                        
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

                    endif

                enddo   !Finish running over components of signs
            
            elseif (tTruncInitiator) then
                ! Determinant in newly spawned list is not found in currentdets - usually this 
                ! would mean the walkers just stay in this list and get merged later - but in 
                ! this case we want to check where the walkers came from - because if the newly 
                ! spawned walkers are from a parent outside the active space they should be 
                ! killed - as they have been spawned on an unoccupied determinant.
                !
                ! If flag_make_initiator is set, then obviously these are allowed to survive
                call extract_sign (SpawnedParts(:,i), SignTemp)
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
                        SignTemp(j) = 0
                        call encode_part_sign (SpawnedParts(:,i), 0, j)
                    endif
                enddo
                if (IsUnoccDet(SignTemp)) then
                    ! All particle 'types' have been aborted
                    ToRemove = ToRemove + 1
                endif
            endif

            ! Even if a corresponding particle wasn't found, we can still
            ! search a smaller list next time....so not all bad news then...
            MinInd=PartInd

        enddo
        
        CALL halt_timer(BinSearch_time)

!        WRITE(6,*) "Leftover Parts..."
!        do i=1,ValidSpawned
!            WRITE(6,*) SpawnedParts(0:NIfTot-1,i),SpawnedParts(NIfTot,i)-2
!        enddo

!Now we have to remove the annihilated particles from the spawned list. They will be removed from the main list at the end of the annihilation process.
!It may actually be easier to just move the annihilated particles to the end of the list and resort the list?
!Or, the removed indices could be found on the fly? This may have little benefit though if the memory isn't needed.
        IF(ToRemove.gt.0) THEN

!Since reading and writing from the same array is slow, copy the information accross to the other spawned array, and just swap the pointers around after.
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
            IF(DetsMerged.ne.ToRemove) THEN
                WRITE(6,*) "***", Iter, DetsMerged, ToRemove
                CALL Stop_All("AnnihilateSpawnedParts","Incorrect number of particles removed from spawned list")
            ENDIF
!We always want to annihilate from the SpawedParts and SpawnedSign arrays, so swap them around.
            PointTemp => SpawnedParts2
            SpawnedParts2 => SpawnedParts
            SpawnedParts => PointTemp

        ENDIF
!        call WriteExcitorListP(6,SpawnedParts,0,ValidSpawned,0,"After zero-removal")

!        WRITE(6,*) "After removal of zeros: "
!        do i=1,ValidSpawned
!            WRITE(6,*) SpawnedParts(:,i)
!        enddo

        CALL halt_timer(AnnMain_time)

    END SUBROUTINE AnnihilateSpawnedParts

    PURE LOGICAL FUNCTION IsUnoccDet(CurrentSign)
        INTEGER, DIMENSION(lenof_sign), INTENT(IN) :: CurrentSign

        IF(lenof_sign.eq.1) THEN
            IsUnoccDet=CurrentSign(1).eq.0
        ELSE
            IF((CurrentSign(1).eq.0).and.(CurrentSign(lenof_sign).eq.0)) THEN
                IsUnoccDet=.true.
            ELSE
                IsUnoccDet=.false.
            ENDIF
        ENDIF
    END FUNCTION IsUnoccDet

    
!This routine will run through the total list of particles (TotWalkersNew in CurrentDets with sign CurrentSign) and the list of newly-spawned but
!non annihilated particles (ValidSpawned in SpawnedParts and SpawnedSign) and move the new particles into the correct place in the new list,
!while removing the particles with sign = 0 from CurrentDets. 
!Binary searching can be used to speed up this transfer substantially.
!The key feature which makes this work, is that it is impossible for the same determinant to be specified in both the spawned and main list at the end of
!the annihilation process. Therefore we will not multiply specify determinants when we merge the lists.
    SUBROUTINE InsertRemoveParts(ValidSpawned,TotWalkersNew)
        use SystemData, only: tHPHF
        use bit_reps, only: NIfD
        use CalcData , only : tCheckHighestPop
        INTEGER :: TotWalkersNew,ValidSpawned
        INTEGER :: i,DetsMerged,nJ(NEl),part_type
        INTEGER, DIMENSION(lenof_sign) :: CurrentSign,SpawnedSign
        REAL*8 :: HDiag
        LOGICAL :: TestClosedShellDet
        character(*), parameter :: this_routine = 'InsertRemoveParts'
        HElement_t :: HDiagTemp

!It appears that the rest of this routine isn't thread-safe if ValidSpawned is zero.
        if(.not.bNodeRoot) return
!Annihilated determinants first are removed from the main array (zero sign). 
!Surely we only need to perform this loop if the number of annihilated particles > 0?
        TotParts=0
        norm_psi_squared = 0
        DetsMerged=0
        iHighestPop=0
        IF(TotWalkersNew.gt.0) THEN
            do i=1,TotWalkersNew
                call extract_sign(CurrentDets(:,i),CurrentSign)
                IF(IsUnoccDet(CurrentSign)) THEN
                    DetsMerged=DetsMerged+1
                    IF(tTruncInitiator) THEN
                        do part_type=1,lenof_sign
                            if (test_flag(CurrentDets(:,i),flag_parent_initiator(part_type))) then
                                !determinant was an initiator...it obviously isn't any more...
                                NoAddedInitiators=NoAddedInitiators-1.D0
                                if (tSpawnSpatialInit) &
                                    call rm_initiator_list (CurrentDets(:,i))
                            endif
                        enddo
                    ENDIF
                ELSE
!We want to move all the elements above this point down to 'fill in' the annihilated determinant.
                    IF(DetsMerged.ne.0) THEN
                        CurrentDets(0:NIfTot,i-DetsMerged)=CurrentDets(0:NIfTot,i)
                        IF(.not.tRegenDiagHEls) THEN
                            CurrentH(i-DetsMerged)=CurrentH(i)
                        ENDIF
                    ENDIF
                    TotParts=TotParts+abs(CurrentSign)
                    norm_psi_squared = norm_psi_squared + sum(CurrentSign**2)
                    IF(tCheckHighestPop) THEN
!If this option is on, then we want to compare the weight on each determinant to the weight at the HF determinant.
!Record the highest weighted determinant on each processor.
!TODO: NOTE: THIS STILL ONLY WORKS EXPLICITLY FOR REAL WALKERS ONLY
                        IF((abs(CurrentSign(1))).gt.iHighestPop) THEN
                            iHighestPop=abs(CurrentSign(1))
                            HighestPopDet(:)=CurrentDets(:,i)
                        ENDIF
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
!        CALL FLUSH(6)

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
                    if (DetBitEQ(CurrentDets(:,i), iLutHF, NIfDBO)) then
!We know we are at HF - HDiag=0
                        HDiag=0.D0
                    else
                        call decode_bit_det (nJ, CurrentDets(:,i))
                        if (tHPHF) then
                            HDiagTemp = hphf_diag_helement (nJ, &
                                                            CurrentDets(:,i))
                        else
                            HDiagTemp = get_helement (nJ, nJ, 0)
                        endif
                        HDiag=(REAL(HDiagTemp,8))-Hii
                    endif
                    CurrentH(i)=HDiag
                enddo
            ELSE
                CALL MergeListswH(TotWalkersNew,ValidSpawned,SpawnedParts(0:NIfTot,1:ValidSpawned))
            ENDIF

        ENDIF
!        CALL CheckOrdering(CurrentDets,CurrentSign(1:TotWalkers),TotWalkers,.true.)

    END SUBROUTINE InsertRemoveParts

    pure function DetermineDetNode (nI, iIterOffset) result(node)
        use SystemData, only: nBasis
        ! Depending on the Hash, determine which node determinant nI
        ! belongs to in the DirectAnnihilation scheme. NB FCIMC has each processor as a separate logical node.
        !
        ! In:  nI   - Integer ordered list for the determinant
        ! In:  iIterOffset - Offset this iteration by this amount
        ! Out: proc - The (0-based) processor index.

        integer, intent(in) :: nI(nel)
        integer, intent(in) :: iIterOffset
        integer :: node
        
        integer :: i
        integer(int64) :: acc
        integer offset

        acc = 0
        offset=hash_iter+iIterOffset
        do i = 1, nel
            acc = (1099511628211_int64 * acc) + &
                    (ishft(RandomHash(mod(iand(nI(i), csf_orbital_mask)+offset-1,nBasis)+1),hash_shift) * i)
!            offset=0
        enddo
        node = abs(mod(acc, int(nNodes, 8)))

    end function
    
    FUNCTION CreateHash(DetCurr)
        INTEGER :: DetCurr(NEl),i
        INTEGER(KIND=int64) :: CreateHash

        CreateHash=0
        do i=1,NEl
!            CreateHash=13*CreateHash+i*DetCurr(i)
            CreateHash=(1099511628211_8*CreateHash)+i*DetCurr(i)
            
!            CreateHash=mod(1099511628211*CreateHash,2**64)
!            CreateHash=XOR(CreateHash,DetCurr(i))
        enddo
!        WRITE(6,*) CreateHash
        RETURN

    END FUNCTION CreateHash

    SUBROUTINE LinSearchParts(DetArray,iLut,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER(KIND=n_int) :: iLut(0:NIfTot),DetArray(0:NIfTot,1:MaxInd)
        INTEGER :: MinInd,MaxInd,PartInd
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
            Comp=DetBitLT(CurrentDets(:,MaxInd),iLut(:),NIfDBO)
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
            Comp=DetBitLT(CurrentDets(:,N),iLut(:),NIfDBO)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
                    Comp=DetBitLT(CurrentDets(:,i+1),iLut(:),NIfDBO)
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
