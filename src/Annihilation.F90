!This module is to be used for various types of walker MC annihilation in serial and parallel.
MODULE AnnihilationMod
    use SystemData , only : NEl, tHPHF
    use CalcData , only : TRegenExcitgens,tRegenDiagHEls
    USE DetCalcData , only : Det,FCIDetIndex
    USE Logging , only : tHistSpawn
    USE Parallel
    USE dSFMT_interface , only : genrand_real2_dSFMT
    USE FciMCData
    use DetBitOps, only: DetBitEQ, DetBitLT, FindBitExcitLevel, ilut_lt, &
                         ilut_gt
    use CalcData , only : tTruncInitiator
    use Determinants, only: get_helement
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use sort_mod
    use constants, only: n_int,lenof_sign,null_part
    use bit_rep_data
    use bit_reps, only: decode_bit_det, extract_sign, extract_flags, &
                        encode_sign, encode_flags, test_flag, set_flag, &
                        clr_flag, flag_parent_initiator, encode_part_sign, &
                        extract_part_sign, copy_flag
    IMPLICIT NONE

    contains

    !TODO:
    !   H Elements - send through logical to decide whether to create or not.
    !   Parallel spawned parts - create the ValidSpawnedList itself.
    !   Going to have to sort this out for the new packaged walkers - will have to package them up in this interface.
    SUBROUTINE AnnihilationInterface(TotDets,MainParts,MainSign,MaxMainInd,SpawnDets,SpawnParts,MaxSpawnInd,iter_data)
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
!       MainSign(:)         This is the signed number of particles on the determinants specified in the equivalent
!                           entry in the MainParticles array.
!       SpawnParts(:,:)     This is the list of particles to attempt to annihilate. Unlike the Main list, this list
!                           does *not* need to be ordered or sign coherent, and can also contain 'zero' sign particles.  Each particle contains its own sign
!       MaxSpawnInd         This is the size of the SpawnParts array.
!       SpawnDets           This is the number of spawned particles in SpawnParts.
!                           entry in the SpawnParts array.

        INTEGER, INTENT(IN) :: MaxMainInd,MaxSpawnInd
        INTEGER, INTENT(INOUT) :: TotDets
        type(fcimc_iter_data), intent(inout) :: iter_data
        INTEGER(KIND=n_int), INTENT(INOUT) , TARGET :: MainParts(0:NIfTot,MaxMainInd),SpawnParts(0:NIfTot,MaxSpawnInd)
        INTEGER, INTENT(INOUT) , TARGET :: MainSign(MaxMainInd)
        INTEGER, INTENT(INOUT) :: SpawnDets
        INTEGER :: ierr,i
        CHARACTER(len=*) , PARAMETER :: this_routine='AnnihilationInterface'
        INTEGER, DIMENSION(lenof_sign) :: TempSign
        TYPE(timer),save :: Annihil_time
        integer(kind=n_int), pointer,save :: SpawnVecLocal(:,:)
        Annihil_time%timer_name='Annihilation interface'
        call set_timer(Annihil_time,20)

        IF(.not.(ALLOCATED(ValidSpawnedList))) THEN
!This needs to be filled correctly before annihilation can take place.
            ALLOCATE(ValidSpawnedList(0:nProcessors-1),stat=ierr)
        ENDIF
        IF(.not.(ASSOCIATED(SpawnVecLocal))) THEN
!This is required scratch space of the size of the spawned arrays
            call shared_allocate_iluts("SpawnVecLocal",SpawnVecLocal,(/NIfTot,MaxSpawnInd/))
            ierr=0
            CALL LogMemAlloc('SpawnVecLocal',MaxSpawnInd*(NIfTot+1),size_n_int,this_routine,SpawnVec2Tag,ierr)
            SpawnVecLocal(:,:)=0
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

        ValidSpawnedList(0)=SpawnDets+1   !Add one since it always indicates the next free slot.
        do i=1,nProcessors-1
            ValidSpawnedList(i)=MaxSpawnInd
        enddo 

        TempSign=0
        do i=1,TotDets
            TempSign(1)=MainSign(i)
            call encode_sign(MainParts(:,i),TempSign)
        enddo
! The SpawnParts already have their signs inside them

        MaxWalkersPart=MaxMainInd
!Point at correct arrays... will need to sort out how these are swapped in the main routine.
        CurrentDets => MainParts
        SpawnedParts => SpawnParts
!These point to the scratch space
        SpawnedParts2 => SpawnVecLocal

        CALL DirectAnnihilation(TotDets, iter_data,.true.) !.true. for single processor annihilation
        if(iProcIndex==root) then        
!Signs put back again into seperate array
           do i=1,TotDets
               call extract_sign(CurrentDets(:,i),TempSign)
               MainSign(i)=TempSign(1)
           enddo
         endif
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

!        WRITE(6,*) "Annihilation finished",MaxIndex,TotWalkersNew

!Put the surviving particles in the main list, maintaining order of the main list.
!Now we insert the remaining newly-spawned particles back into the original list (keeping it sorted), and remove the annihilated particles from the main list.
        CALL set_timer(Sort_Time,30)
        CALL InsertRemoveParts(MaxIndex,TotWalkersNew,tSingleProc)

!       WRITE(6,*) "Surviving particles merged"

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

!        WRITE(6,*) "ValidSpawnedList ",ValidSpawnedList(:)
        if(tSingleProc) then
!Put all particles and gap on one proc.

!ValidSpawnedList(0:nProcessors-1) indicates the next free index for each processor (for spawnees from this processor)
!  i.e. the list of spawned particles has already been arranged so that newly spawned particles are grouped according to the processor they go to.

! sendcounts(1:) indicates the number of spawnees to send to each processor (1-based)
! disps(1:) is the index into the spawned list of the beginning of the list to send to each processor (1-based)
           sendcounts(1)=ValidSpawnedList(0)-1
           disps(1)=0
           if(nProcessors>1) then
              sendcounts(2:nProcessors)=0
              disps(2:nProcessors)=ValidSpawnedList(1)
           endif
        else
!Distribute the gaps on all procs
!TODO:  This should use information from InitialSpawnedSlots for consistency.
            Gap=REAL(MaxSpawned)/REAL(nProcessors)
   !        WRITE(6,*) "Gap: ",Gap

           do i=0,nProcessors-1
               sendcounts(i+1)=ValidSpawnedList(i)-(NINT(Gap*i)+1)
               disps(i+1)=NINT(Gap*i)
           enddo
        endif


!        WRITE(6,*) 'sendcounts',sendcounts
!        WRITE(6,*) 'disps',disps

        MaxSendIndex=ValidSpawnedList(nProcessors-1)-1

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
            CALL Warning("SendProcNewParts","Maximum index of newly-spawned array is close to maximum length after annihilation send. Increase MemoryFacSpawn")
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
!        do i=1,nProcessors
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
        CALL MPIAlltoAllv(SpawnedParts,sendcounts,disps,SpawnedParts2,recvcounts,recvdisps,error)

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
        LOGICAL :: tSuc
        INTEGER(Kind=n_int) , POINTER :: PointTemp(:,:)
        integer(n_int) :: cum_det (0:niftot)
        CHARACTER(len=*), parameter :: this_routine='CompressSpawnedList'

!We want to sort the list of newly spawned particles, in order for quicker binary searching later on. (this is not essential, but should proove faster)
!They should remain sorted after annihilation between spawned
        call sort(SpawnedParts(:,1:ValidSpawned), ilut_lt, ilut_gt)
        IF(tHistSpawn) HistMinInd2(1:NEl)=FCIDetIndex(1:NEl)

!        WRITE(6,*) "************ - Ordered",ValidSpawned,NIfTot,Iter
!        do i=1,ValidSpawned
!            WRITE(6,*) i,SpawnedParts(0:NIfTot-1,i),SpawnedParts(NIfTot,i)-2
!        enddo

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
                    if (.not. tTruncInitiator .or. &
                        .not. test_flag (SpawnedParts(:,i), flag_parent_initiator(part_type))) then
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
!        WRITE(6,*) "Compressed List: ",ValidSpawned,DetsMerged
!        do i=1,ValidSpawned
!            WRITE(6,*) SpawnedParts(0:NIfTot-1,i),SpawnedParts(NIfTot,i)-2
!        enddo
!        WRITE(6,*) "***","iter_data%naborted: ",iter_data%naborted
!        CALL FLUSH(6)
        
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
                                CALL BinSearchParts2(SpawnedParts(:,i),HistMinInd2(ExcitLevel),FCIDetIndex(ExcitLevel+1)-1,PartIndex,tSuc)
                            ENDIF
                            HistMinInd2(ExcitLevel)=PartIndex
                            IF(tSuc) THEN
                                AvAnnihil(j,PartIndex)=AvAnnihil(j,PartIndex)+REAL(2*(min(abs(CurrentSign(j)),abs(SpawnedSign(j)))))
                                InstAnnihil(j,PartIndex)=InstAnnihil(j,PartIndex)+REAL(2*(min(abs(CurrentSign(j)),abs(SpawnedSign(j)))))
                            ELSE
                                WRITE(6,*) "***",SpawnedParts(0:NIftot,i)
                                Call WriteBitDet(6,SpawnedParts(0:NIfTot,i),.true.)
                                CALL Stop_All("AnnihilateSpawnedParts","Cannot find corresponding FCI determinant when histogramming")
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
    SUBROUTINE InsertRemoveParts(ValidSpawned,TotWalkersNew,tSingleProc)
        use SystemData, only: tHPHF
        use bit_reps, only: NIfD
        use CalcData , only : tCheckHighestPop
        INTEGER :: TotWalkersNew,ValidSpawned
        INTEGER :: i,DetsMerged,nJ(NEl),part_type
        INTEGER, DIMENSION(lenof_sign) :: CurrentSign,SpawnedSign
        REAL*8 :: HDiag
        LOGICAL :: TestClosedShellDet
        LOGICAL, intent(in) ::tSingleProc
        HElement_t :: HDiagTemp

!It appears that the rest of this routine isn't thread-safe if ValidSpawned is zero.
        if(tSingleProc.and.iProcIndex/=root) return
!Annihilated determinants first are removed from the main array (zero sign). 
!Surely we only need to perform this loop if the number of annihilated particles > 0?
        TotParts=0
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
                    IF(tCheckHighestPop) THEN
!If this option is on, then we want to compare the weight on each determinant to the weight at the HF determinant.
!Record the highest weighted determinant on each processor.
!TODO: NOTE: THIS STILL ONLY WORKS EXPLICITLY FOR REAL WALKERS ONLY
                        IF((abs(CurrentSign(1))).gt.iHighestPop) THEN
                            IF(tHPHF) THEN
                                !For HPHF functions, we restrict ourselves to closed shell determinants for simplicity.
                                IF(TestClosedShellDet(CurrentDets(0:NIfDBO,i))) THEN
                                    !HPHF func is closed shell - we can move to this without problems.
                                    iHighestPop=abs(CurrentSign(1))
                                    HighestPopDet(:)=CurrentDets(:,i)
                                ENDIF
                            ELSE
                                iHighestPop=abs(CurrentSign(1))
                                HighestPopDet(:)=CurrentDets(:,i)
                            ENDIF
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
        ENDIF
        do i=2,ValidSpawned
            call extract_sign(SpawnedParts(:,i),SpawnedSign)
            TotParts=TotParts+abs(SpawnedSign)
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
    
    INTEGER FUNCTION DetermineDetProc(iLut)
        use bit_reps, only: NIfDBO
        use CalcData , only: tRandomiseHashOrbs
        use FciMCData, only: RandomHash 
        use constants, only: bits_n_int
        INTEGER(KIND=n_int) :: iLut(0:NIfTot)
        INTEGER :: i,j,Elecs!,TempDet(NEl),MurmurHash2Wrapper
        INTEGER(KIND=int64) :: Summ!,RangeofBins,NextBin

!        CALL DecodeBitDet(TempDet,iLut)
!        i=MurmurHash2Wrapper(TempDet,NEl,13)
!        write(6,*) i
        IF(tRandomiseHashOrbs) THEN
            Summ=0
            Elecs=0
            lp1: do i=0,NIfDBO
                do j=0,bits_n_int-1
                    IF(BTEST(iLut(i),j)) THEN
                        Elecs=Elecs+1
#ifdef __INT64
                        Summ=(1099511628211_8*Summ)+RandomHash((i*64)+(j+1))*Elecs
#else
                        Summ=(1099511628211_8*Summ)+RandomHash((i*32)+(j+1))*Elecs
#endif
                        IF(Elecs.eq.NEl) EXIT lp1
                    ENDIF
                enddo
            enddo lp1
            DetermineDetProc=abs(mod(Summ,INT(nProcessors,8)))

        ELSE
            Summ=0
            Elecs=0
            lp2: do i=0,NIfDBO
                do j=0,bits_n_int-1
                    IF(BTEST(iLut(i),j)) THEN
                        Elecs=Elecs+1
#ifdef __INT64
                        Summ=(1099511628211_8*Summ)+((i*64)+(j+1))*Elecs
#else
                        Summ=(1099511628211_8*Summ)+((i*32)+(j+1))*Elecs
#endif
                        IF(Elecs.eq.NEl) EXIT lp2
                    ENDIF
                enddo
            enddo lp2
            DetermineDetProc=abs(mod(Summ,INT(nProcessors,8)))
        ENDIF
!        WRITE(6,*) iLut(0:niftot),DetermineDetProc,Summ,nProcessors

    END FUNCTION DetermineDetProc

    
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
    
!!This is a new routine to totally annihilate all particles on the same determinant. This is not done using an all-to-all, but rather
!!by rotating the newly spawned particles around all determinants and annihilating with the particles on their processor.
!!Valid spawned is the number of newly-spawned particles. These 'particles' can be multiply specified on the same determinant.
!!Each rotation and annihilation step, the number corresponds to a different processors spawned particles.
!!TotWalkersNew indicates the number of particles in NewDets - the list of particles to compare for annihilation.
!!Improvements in AnnihilateBetweenSpawned:
!!Binary search for sendcounts and others, and only transfer all data when need to.
!!Memory improvements
!!Call as one array for All-to-alls
!!Make sure only sort what need to
!    SUBROUTINE RotoAnnihilation(ValidSpawned,TotWalkersNew)
!        INTEGER :: ValidSpawned,TotWalkersNew,i
!        INTEGER :: ierr,error!,SpawnedBeforeRoto
!        CHARACTER , ALLOCATABLE :: mpibuffer(:)
!
!!        InitialSpawned=TotSpawned     !Initial spawned will store the original number of spawned particles, so that we can compare afterwards.
!!        InitialSpawned=Annihilated
!        
!        CALL CompressSpawnedList(ValidSpawned)
!
!!        CALL SortBitDets(ValidSpawned,SpawnedParts(0:NIfD,1:ValidSpawned),NIfD,SpawnedSign(1:ValidSpawned))
!        CALL MPIBarrier(error)
!!        WRITE(6,*) "Entering rotoannilation: ",Iter,InitialSpawned,TotWalkersNew
!!        CALL FLUSH(6)
!
!!First, annihilate between newly spawned particles. Memory for this will be allocated dynamically.
!!This will be done in the usual fashion using the All-to-All communication and hashes.
!        CALL AnnihilateBetweenSpawned(ValidSpawned)
!!        CALL AnnihilateBetweenSpawnedOneProc(ValidSpawned)
!!        Annihilated=Annihilated+(InitialSpawned-TotSpawned)
!!        IF(Annihilated.ne.InitialSpawned) THEN
!!            WRITE(6,*) "Have annihilated between newly-spawned...",Annihilated-InitialSpawned,Iter
!!        ENDIF
!            
!!We want to sort the list of newly spawned particles, in order for quicker binary searching later on. (this is not essential, but should proove faster)
!!        CALL SortBitDets(ValidSpawned,SpawnedParts(0:NIfD,1:ValidSpawned),NIfD,SpawnedSign(1:ValidSpawned))
!!        CALL CheckOrdering(SpawnedParts(:,1:ValidSpawned),SpawnedSign(1:ValidSpawned),ValidSpawned,.true.)
!!        do i=1,ValidSpawned
!!            WRITE(6,*) 1,i,SpawnedParts(:,i),SpawnedSign(i),Iter
!!            CALL FLUSH(6)
!!        enddo
!
!!        SpawnedBeforeRoto=ValidSpawned
!!        WRITE(6,*) "SpawnedBeforeRoto: ",ValidSpawned
!
!!This RemoveInds is useful scratch space for the removal of particles from lists. It probably isn't essential, but keeps things simpler initially.
!!        ALLOCATE(RemoveInds(MaxSpawned),stat=ierr)
!        
!!This routine annihilates the processors set of newly-spawned particles, with the complete set of particles on the processor.
!        CALL AnnihilateSpawnedParts(ValidSpawned,TotWalkersNew)
!
!        CALL MPIBarrier(error)
!
!!Allocate a buffer here to hold particles when using a buffered send...
!!The buffer wants to be able to hold (MaxSpawned+1)x(NIfD+2) integers (*4 for in bytes). If we could work out the maximum ValidSpawned accross the determinants,
!!it could get reduced to this... 
!        IF(nProcessors.ne.1) THEN
!            ALLOCATE(mpibuffer(8*(MaxSpawned+1)*(2*NIfTot+2)),stat=ierr)
!            IF(ierr.ne.0) THEN
!                CALL Stop_All("RotoAnnihilation","Error allocating memory for transfer buffers...")
!            ENDIF
!#ifdef PARALLEL
!            CALL MPI_Buffer_attach(mpibuffer,8*(MaxSpawned+1)*(2*NIfTot+2),error)
!#endif
!            IF(error.ne.0) THEN
!                CALL Stop_All("RotoAnnihilation","Error allocating memory for transfer buffers...")
!            ENDIF
!        ENDIF
!
!        do i=1,nProcessors-1
!!Move newly-spawned particles which haven't been annihilated around the processors in sequence, annihilating locally each step.
!!This moves the set of newly-spawned particles on this processor one to the right, and recieves from the left.
!!This also updates the ValidSpawned variable so that it now refers to the new set of spawned-particles.
!            CALL RotateParticles(ValidSpawned)
!!            WRITE(6,*) "Rotating particles for the ",i," time...",Iter
!!            CALL FLUSH(6)
!
!!This routine annihilates the processors set of newly-spawned particles, with the complete set of particles on the processor.
!            CALL AnnihilateSpawnedParts(ValidSpawned,TotWalkersNew)
!!            CALL CheckOrdering(SpawnedParts(:,1:ValidSpawned),SpawnedSign(1:ValidSpawned),ValidSpawned,.true.)
!!            WRITE(6,*) "Annihilated locally....",i
!!            CALL FLUSH(6)
!
!        enddo
!
!!One final rotation means that the particles are all on their original processor.
!        IF(nProcessors.ne.1) THEN
!            CALL RotateParticles(ValidSpawned)
!
!#ifdef PARALLEL
!!Detach buffers
!            CALL MPI_Buffer_detach(mpibuffer,8*(MaxSpawned+1)*(2*NIfTot+2),error)
!#endif
!            DEALLOCATE(mpibuffer)
!        ENDIF
!        
!!Test that we have annihilated the correct number here (from each lists), and calculate Annihilated for each processor.
!!Now we insert the remaining newly-spawned particles back into the original list (keeping it sorted), and remove the annihilated particles from the main list.
!        CALL set_timer(Sort_Time,30)
!!        WRITE(6,*) "Entering insert/remove..."
!!        CALL FLUSH(6)
!        CALL InsertRemoveParts(ValidSpawned,TotWalkersNew)
!        CALL halt_timer(Sort_Time)
!
!!        DEALLOCATE(RemoveInds)
!
!
!    END SUBROUTINE RotoAnnihilation
!
!
!    
!    SUBROUTINE AnnihilateBetweenSpawnedOneProc(ValidSpawned)
!        INTEGER(KIND=n_int) :: DetCurr(0:NIfTot)
!        INTEGER :: ValidSpawned,i,j,k,LowBound,HighBound,WSign
!        INTEGER :: VecSlot,TotSign
!
!        call sort (SpawnedParts(:,1:ValidSpawned),SpawnedSign(1:ValidSpawned))
!
!        VecSlot=1
!        i=1
!        do while(i.le.ValidSpawned)
!            LowBound=i
!            DetCurr(0:NIfTot)=SpawnedParts(0:NIfTot,i)
!            i=i+1
!            do while(DetBitEQ(DetCurr(0:NIfTot),SpawnedParts(0:NIfTot,i),NIfDBO).and.(i.le.ValidSpawned))
!                i=i+1
!            enddo
!            HighBound=i-1
!
!!Now, run through the block of common particles again, counting the residual sign
!            TotSign=0
!            do j=LowBound,HighBound
!                TotSign=TotSign+SpawnedSign(j)
!            enddo
!
!!Now, fill up SpawnedSign2 and SpawnedParts2 with the residual particles
!            IF(TotSign.ne.0) THEN
!                WSign=INT(TotSign/abs(TotSign))
!                do k=1,abs(TotSign)
!                    SpawnedParts2(0:NIfTot,VecSlot)=DetCurr(0:NIfTot)
!                    SpawnedSign2(VecSlot)=WSign
!                    VecSlot=VecSlot+1
!                enddo
!            ENDIF
!
!        enddo
!
!        ValidSpawned=VecSlot-1
!
!        do i=1,ValidSpawned
!            SpawnedParts(0:NIfTot,i)=SpawnedParts2(0:NIfTot,i)
!            SpawnedSign(i)=SpawnedSign2(i)
!        enddo
!
!    END SUBROUTINE AnnihilateBetweenSpawnedOneProc
!
!
!!This routine wants to take the ValidSpawned particles in the SpawnedParts array and perform All-to-All communication so that 
!!we can annihilate all common particles with opposite signs.
!!Particles are fed in on the SpawnedParts and SpawnedSign array, and are returned in the same arrays.
!!It requires MaxSpawned*36 bytes of memory (on top of the memory of the arrays fed in...)
!!Might not need to send hashes in all-to-all - could just use them for determining where they go
!!Package up temp arrays?
!    SUBROUTINE AnnihilateBetweenSpawned(ValidSpawned)
!        use CalcData, only: tReadPops,tAnnihilatebyRange
!        integer(int64), allocatable :: HashArray1(:), HashArray2(:), HashArrayTmp(:)
!        INTEGER , ALLOCATABLE :: IndexTable1(:),IndexTable2(:),ProcessVec1(:),ProcessVec2(:),TempSign(:)
!        INTEGER :: i,j,k,ToAnnihilateIndex,ValidSpawned,ierr,error,sendcounts(nProcessors)
!        INTEGER :: TotWalkersDet,InitialBlockIndex,FinalBlockIndex,ToAnnihilateOnProc,VecSlot
!        INTEGER :: disps(nProcessors),recvcounts(nProcessors),recvdisps(nProcessors),nJ(NEl)
!        INTEGER :: Minsendcounts,Maxsendcounts,DebugIter,SubListInds(2,nProcessors),MinProc,MinInd
!        INTEGER(KIND=int64) :: HashCurr,MinBin,RangeofBins,NextBinBound,MinHash
!        CHARACTER(len=*), PARAMETER :: this_routine='AnnihilateBetweenSpawned'
!
!        CALL set_timer(AnnSpawned_time,30)
!
!!First, we need to allocate memory banks. Each array needs a hash value, a processor value, and an index value.
!!We also want to allocate a temporary sign value
!        ALLOCATE(TempSign(ValidSpawned),stat=ierr)
!
!!These arrays may as well be kept all the way through the simulation?
!        ALLOCATE(HashArray1(MaxSpawned),stat=ierr)
!        ALLOCATE(HashArray2(MaxSpawned),stat=ierr)
!        ALLOCATE(HashArrayTmp(MaxSpawned),stat=ierr)
!        ALLOCATE(IndexTable1(MaxSpawned),stat=ierr)
!        ALLOCATE(IndexTable2(MaxSpawned),stat=ierr)
!        ALLOCATE(ProcessVec1(MaxSpawned),stat=ierr)
!        ALLOCATE(ProcessVec2(MaxSpawned),stat=ierr)
!
!        IF(ierr.ne.0) THEN
!            CALL Stop_All("AnnihilateBetweenSpawned","Error in allocating initial data")
!        ENDIF
!
!        TempSign(1:ValidSpawned)=SpawnedSign(1:ValidSpawned)
!        ProcessVec1(1:ValidSpawned)=iProcIndex
!
!!        WRITE(6,*) "***************************************"
!        do i=1,ValidSpawned
!            IndexTable1(i)=i
!            call decode_bit_det (nJ, SpawnedParts(0:NIfTot,i))
!            HashArray1(i)=CreateHash(nJ)
!!            IF(Iter.eq.1346.and.(HashArray1(i).eq.2905380077198165348)) THEN
!!                WRITE(6,*) "Hash found, ",i,SpawnedSign(i),HashArray1(i),SpawnedParts(0:NIfTot,i)
!!            ENDIF
!        enddo
!
!!Next, order the hash array, taking the index, CPU and sign with it...
!        IF(.not.tAnnihilatebyRange) THEN
!!Order the array by abs(mod(Hash,nProcessors)). This will result in a more load-balanced system (no need to actually take ProcessVec1 - this will always be iProcIndex here.
!
!            HashArrayTmp(1:ValidSpawned) = &
!                        abs(mod(HashArray1(1:ValidSpawned), int(nProcessors,8)))
!            call sort (HashArrayTmp(1:ValidSpawned), &
!                       HashArray1(1:ValidSpawned), &
!                       IndexTable1(1:ValidSpawned), &
!                       ProcessVec1(1:ValidSpawned), &
!                       SpawnedSign(1:ValidSpawned))
!
!!Send counts is the size of each block of ordered dets which are going to each processor. This could be binary searched for extra speed
!            IF(ValidSpawned.gt.0) THEN
!                j=1
!                do i=0,nProcessors-1    !Search through all possible values of abs(mod(Hash,nProcessors))
!                    do while((abs(mod(HashArray1(j),INT(nProcessors,8))).eq.i).and.(j.le.ValidSpawned))
!                        j=j+1
!                    enddo
!                    sendcounts(i+1)=j-1
!                enddo
!            ELSE
!                sendcounts(1:nProcessors)=0
!            ENDIF
!
!        ELSE
!            ! We can try to sort the hashes by range, which may result in 
!            ! worse load-balancing, but will remove the need for a second 
!            ! sort of the hashes once they have been sent to the correct 
!            ! processor.
!            call sort (HashArray1(1:ValidSpawned), &
!                       IndexTable1(1:ValidSpawned), &
!                       ProcessVec1(1:ValidSpawned), &
!                       SpawnedSign(1:ValidSpawned))
!!We also need to know the ranges of the hashes to send to each processor. Each range should be the same.
!            IF(nProcessors.ne.1) THEN
!                Rangeofbins=INT(HUGE(Rangeofbins)/(nProcessors/2),8)
!                MinBin=-HUGE(MinBin)
!                NextBinBound=MinBin+Rangeofbins
!
!!We need to find the indices for each block of hashes which are to be sent to each processor.
!!Sendcounts is the size of each block of ordered dets which are going to each processors. This could be binary searched for extra speed.
!                j=1
!                do i=1,nProcessors    !Search through all possible values of the hashes
!                    do while((HashArray1(j).le.NextBinBound).and.(j.le.ValidSpawned))
!                        j=j+1
!                    enddo
!                    sendcounts(i)=j-1
!                    IF(i.eq.nProcessors-1) THEN
!!Make sure the final bin catches everything...
!                        NextBinBound=HUGE(NextBinBound)
!                    ELSE
!                        NextBinBound=NextBinBound+Rangeofbins
!                    ENDIF
!                enddo
!            ELSE
!                sendcounts(1)=ValidSpawned
!!                do j=1,ValidSpawned
!!                    WRITE(6,*) Iter,j,HashArray1(j),SpawnedSign(j)
!!                enddo
!                    
!            ENDIF
!        ENDIF
!
!        IF(sendcounts(nProcessors).ne.ValidSpawned) THEN
!            WRITE(6,*) "SENDCOUNTS is: ",sendcounts(:)
!            WRITE(6,*) "VALIDSPAWNED is: ",ValidSpawned
!            CALL FLUSH(6)
!            CALL Stop_All("AnnihilateBetweenSpawned","Incorrect calculation of sendcounts")
!        ENDIF
!
!!Oops, we have calculated them cumulativly - undo this
!        maxsendcounts=sendcounts(1)
!        minsendcounts=sendcounts(1)     !Find max & min sendcounts, so that load-balancing can be checked
!!        WRITE(6,*) maxsendcounts,minsendcounts
!        do i=2,nProcessors
!            do j=1,i-1
!                sendcounts(i)=sendcounts(i)-sendcounts(j)
!            enddo
!            IF(sendcounts(i).gt.maxsendcounts) THEN
!                maxsendcounts=sendcounts(i)
!            ELSEIF(sendcounts(i).lt.minsendcounts) THEN
!                minsendcounts=sendcounts(i)
!            ENDIF
!        enddo
!
!!The disps however do want to be cumulative - this is the array indexing the start of the data block
!        disps(1)=0      !Starting element is always the first element
!        do i=2,nProcessors
!            disps(i)=disps(i-1)+sendcounts(i-1)
!        enddo
!
!!We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
!        recvcounts(1:nProcessors)=0
!
!        CALL MPIAlltoAllI(sendcounts,1,recvcounts,1,error)
!
!!We can now get recvdisps from recvcounts in the same way we obtained disps from sendcounts
!        recvdisps(1)=0
!        do i=2,nProcessors
!            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
!        enddo
!
!        MaxIndex=recvdisps(nProcessors)+recvcounts(nProcessors)
!!Max index is the largest occupied index in the array of hashes to be ordered in each processor 
!        IF(MaxIndex.gt.(0.93*MaxSpawned)) THEN
!            CALL Warning("AnnihilateBetweenSpawned","Maximum index of annihilation array is close to maximum length. Increase MemoryFacSpawn")
!            IF(tReadPops) CALL Warning("AnnihilateBetweenSpawned","When reading in a POPSFILE, MemoryFacSpawn must be greater than 1.0")
!        ENDIF
!
!!Uncomment this if you want to write out load-balancing statistics.
!!        AnnihilPart(:)=0
!!        CALL MPI_Gather(MaxIndex,1,MPI_INTEGER,AnnihilPart,1,MPI_INTEGER,root,MPI_COMM_WORLD,error)
!!        IF(iProcIndex.eq.root) THEN
!!            WRITE(13,"(I10)",advance='no') Iter
!!            do i=1,nProcessors
!!                WRITE(13,"(I10)",advance='no') AnnihilPart(i)
!!            enddo
!!            WRITE(13,"(A)") ""
!!            CALL FLUSH(13)
!!        ENDIF
!
!!        IF(Iter.eq.DebugIter) THEN
!!            WRITE(6,*) "RECVCOUNTS: "
!!            WRITE(6,*) recvcounts(:)
!!            WRITE(6,*) "RECVDISPS: "
!!            WRITE(6,*) recvdisps(:),MaxIndex
!!            CALL FLUSH(6)
!!        ENDIF
!
!!Insert a load-balance check here...maybe find the s.d. of the sendcounts array - maybe just check the range first.
!!        IF(TotWalkersNew.gt.200) THEN
!!            IF((Maxsendcounts-Minsendcounts).gt.(TotWalkersNew/3)) THEN
!!                WRITE(6,"(A,I12)") "**WARNING** Parallel annihilation not optimally balanced on this node, for iter = ",Iter
!!                WRITE(6,*) "Sendcounts is: ",sendcounts(:)
!!                CALL FLUSH(6)
!!            ENDIF
!!        ENDIF
!
!!Now send the chunks of hashes to the corresponding processors
!        CALL MPIAlltoAllvI8(HashArray1(1:ValidSpawned),sendcounts,disps,HashArray2(1:MaxIndex),recvcounts,recvdisps,error)
!
!!The signs of the hashes, index and CPU also need to be taken with them.
!        CALL MPIAlltoAllvI(SpawnedSign(1:ValidSpawned),sendcounts,disps,SpawnedSign2(1:MaxIndex),recvcounts,recvdisps,error)
!        CALL MPIAlltoAllvI(IndexTable1(1:ValidSpawned),sendcounts,disps,IndexTable2,recvcounts,recvdisps,error)
!        CALL MPIAlltoAllvI(ProcessVec1(1:ValidSpawned),sendcounts,disps,ProcessVec2,recvcounts,recvdisps,error)
!
!        IF(.not.tAnnihilatebyrange) THEN
!            ! The hashes now need to be sorted again - this time by their 
!            ! number. This sorting would be redundant if we had initially 
!            ! sorted the hashes by range (ie tAnnihilatebyrange).
!            call sort (HashArray2(1:MaxIndex), &
!                       IndexTable2(1:MaxIndex), &
!                       ProcessVec2(1:MaxIndex), &
!                       SpawnedSign2(1:MaxIndex))
!        ELSE
!!Here, because we have ordered the hashes initially numerically, we have a set of ordered lists. It is therefore easier to sort them.
!!We have to work out how to run sequentially through the hashes, which are a set of nProc seperate ordered lists.
!!We would need to have 2*nProc indices, since we will have a set of nProc disjoint ordered sublists.
!!SubListInds(1,iProc)=index of current hash from processor iProc
!!SubListInds(2,iProc)=index of final hash from processor iProc
!!Indices can be obtained from recvcounts and recvdisps - recvcounts(iProc-1) is number of hashes from iProc
!!recvdisps(iProc-1) is the displacement to the start of the hashes from iProc
!            do i=1,nProcessors-1
!                SubListInds(1,i)=recvdisps(i)+1
!                SubListInds(2,i)=recvdisps(i+1)
!            enddo
!            SubListInds(1,nProcessors)=recvdisps(nProcessors)+1
!            SubListInds(2,nProcessors)=MaxIndex
!!            WRITE(6,*) "SubListInds(1,:) ", SubListInds(1,:)
!!            WRITE(6,*) "SubListInds(2,:) ", SubListInds(2,:)
!!            WRITE(6,*) "Original hash list is: "
!!Reorder the lists so that they are in numerical order.
!            j=1
!            do while(j.le.MaxIndex)
!                do i=1,nProcessors
!                    IF(SubListInds(1,i).le.SubListInds(2,i)) THEN
!!This block still has hashes which want to be sorted
!                        MinHash=HashArray2(SubListInds(1,i))
!                        MinProc=i
!                        MinInd=SubListInds(1,i)
!                        EXIT
!                    ENDIF
!!                    IF(i.eq.nProcessors) THEN
!!                        WRITE(6,*) "ERROR HERE!!"
!!                        CALL FLUSH(6)
!!                    ENDIF
!                enddo
!                IF(MinHash.ne.HashCurr) THEN
!                    do i=MinProc+1,nProcessors
!                        IF((SubListInds(1,i).le.SubListInds(2,i)).and.(HashArray2(SubListInds(1,i)).lt.MinHash)) THEN
!                            MinHash=HashArray2(SubListInds(1,i))
!                            MinProc=i
!                            MinInd=SubListInds(1,i)
!                            IF(MinHash.eq.HashCurr) THEN
!                                EXIT
!                            ENDIF
!                        ENDIF
!                    enddo
!                ENDIF
!!Next smallest hash is MinHash - move the ordered elements into the other array.
!                HashArray1(j)=MinHash
!                IndexTable1(j)=IndexTable2(MinInd)
!                ProcessVec1(j)=ProcessVec2(MinInd)
!                SpawnedSign(j)=SpawnedSign2(MinInd)
!                HashCurr=MinHash
!!Move through the block
!                j=j+1
!                SubListInds(1,MinProc)=SubListInds(1,MinProc)+1
!            enddo
!
!            IF((j-1).ne.MaxIndex) THEN
!                CALL Stop_All(this_routine,"Error here in the merge sort algorithm")
!            ENDIF
!
!!Need to copy the lists back to the original array to fit in with the rest of the code
!            do i=1,MaxIndex
!                IndexTable2(i)=IndexTable1(i)
!                ProcessVec2(i)=ProcessVec1(i)
!                SpawnedSign2(i)=SpawnedSign(i)
!                HashArray2(i)=HashArray1(i)
!            enddo
!
!        ENDIF
!
!!Work out the index of the particles which want to be annihilated
!        j=1
!        ToAnnihilateIndex=1
!        do while(j.le.MaxIndex)
!            TotWalkersDet=0
!            InitialBlockIndex=j
!            FinalBlockIndex=j-1         !Start at j-1 since we are increasing FinalBlockIndex even with the first det in the next loop
!            HashCurr=HashArray2(j)
!            do while((HashArray2(j).eq.HashCurr).and.(j.le.MaxIndex))
!!First loop counts walkers in the block - TotWalkersDet is then the residual sign of walkers on that determinant
!                TotWalkersDet=TotWalkersDet+SpawnedSign2(j)
!
!!                IF(SpawnedSign2(j).eq.1) THEN
!!                    TotWalkersDet=TotWalkersDet+1
!!                ELSE
!!                    TotWalkersDet=TotWalkersDet-1
!!                ENDIF
!                FinalBlockIndex=FinalBlockIndex+1
!                j=j+1
!            enddo
!
!!            IF((Iter.eq.1877)) THEN
!!                WRITE(6,*) "Common block of dets found from ",InitialBlockIndex," ==> ",FinalBlockIndex
!!                WRITE(6,*) "Sum of signs in block is: ",TotWalkersDet,HashCurr
!!                do k=InitialBlockIndex,FinalBlockIndex
!!                    WRITE(6,*) TotWalkersDet,ToAnnihilateIndex,IndexTable2(k),ProcessVec2(k),SpawnedSign2(k)
!!                enddo
!!                CALL FLUSH(6)
!!            ENDIF
!!We need to now run through the block, and count of the same number of surviving particles as given by TotWalkersDet
!    ! 1. If particles are of opposite sign, then annihilation
!    ! 2. If particles are of same sign, then count out until we have the required number and annihilate the rest.
!    ! Now, the sign has to be passed back. This will indicate the sign of the SURVIVING particles on that determinant.
!    ! ToAnnihilateIndex now indicates the number of particles who want their sign changed at all...
!
!            do k=InitialBlockIndex,FinalBlockIndex
!!Second run through the block of same determinants marks walkers for annihilation
!                IF(TotWalkersDet.eq.0) THEN
!!All walkers in block want to be annihilated from now on.
!                    IndexTable1(ToAnnihilateIndex)=IndexTable2(k)
!                    ProcessVec1(ToAnnihilateIndex)=ProcessVec2(k)
!                    SpawnedSign(ToAnnihilateIndex)=0   
!                    ToAnnihilateIndex=ToAnnihilateIndex+1
!                    Annihilated=Annihilated+abs(SpawnedSign2(k))
!!                    TotSpawned=TotSpawned-abs(SpawnedSign2(k))
!                ELSEIF((TotWalkersDet.lt.0).and.(SpawnedSign2(k).gt.0)) THEN
!!Annihilate if block has a net negative walker count, and current walker is positive
!                    IndexTable1(ToAnnihilateIndex)=IndexTable2(k)
!                    ProcessVec1(ToAnnihilateIndex)=ProcessVec2(k)
!                    SpawnedSign(ToAnnihilateIndex)=0
!                    ToAnnihilateIndex=ToAnnihilateIndex+1
!                    Annihilated=Annihilated+SpawnedSign2(k)
!!                    TotSpawned=TotSpawned-SpawnedSign2(k)
!                ELSEIF((TotWalkersDet.gt.0).and.(SpawnedSign2(k).lt.0)) THEN
!!Annihilate if block has a net positive walker count, and current walker is negative
!                    IndexTable1(ToAnnihilateIndex)=IndexTable2(k)
!                    ProcessVec1(ToAnnihilateIndex)=ProcessVec2(k)
!                    SpawnedSign(ToAnnihilateIndex)=0
!                    ToAnnihilateIndex=ToAnnihilateIndex+1
!                    Annihilated=Annihilated-SpawnedSign2(k)
!!                    TotSpawned=TotSpawned+SpawnedSign2(k)
!                ELSE
!!If net walkers is positive, and we have a positive walkers, then remove one from the net positive walkers and continue through the block
!!Now, we have a particle which is the same sign as the residual sign we want to pass through.
!!If the sign on the particle is equal to, or less than the residual sign, then we want to let all particles live.
!!Otherwise, we want to annihilate a fraction of them...
!                    IF((abs(TotWalkersDet)).ge.(abs(SpawnedSign2(k)))) THEN
!!All these particles are ok to be transferred accross...Increase (SpawnedSign2(k) < 0) the sign on totwalkersdet)
!                        TotWalkersDet=TotWalkersDet-SpawnedSign2(k)
!                    ELSE
!!There is a greater number of particles in this entry than the total residual sign. Therefore, this entry want to be PARTIALLY annihilated.
!!SpawnedSign will indicate the number of particles we want to remain on this entry.
!                        IndexTable1(ToAnnihilateIndex)=IndexTable2(k)
!                        ProcessVec1(ToAnnihilateIndex)=ProcessVec2(k)
!                        Annihilated=Annihilated+(abs(SpawnedSign2(k))-abs(TotWalkersDet))
!!                        TotSpawned=TotSpawned-(abs(SpawnedSign2(k))-abs(TotWalkersDet))
!                        SpawnedSign(ToAnnihilateIndex)=TotWalkersDet    !The number of particles that we want left to copy accross is simply the remaining residual sign
!                        ToAnnihilateIndex=ToAnnihilateIndex+1
!                        TotWalkersDet=0     !All the residual sign has now been compensated for.
!
!                    ENDIF
!
!                ENDIF
!            enddo
!            IF(TotWalkersDet.ne.0) THEN
!                CALL Stop_All("AnnihilateBetweenSpawned","Problem counting residual sign...")
!            ENDIF
!
!        enddo
!
!        ToAnnihilateIndex=ToAnnihilateIndex-1   !ToAnnihilateIndex now tells us the total number of particles to annihilate from the list on this processor
!
!!The annihilation is complete - particles to be annihilated are stored in IndexTable and need to be sent back to their original processor
!!To know which processor that is, we need to order the particles to be annihilated in terms of their CPU, i.e. ProcessVec(1:ToAnnihilateIndex)
!!Is the list already ordered according to CPU? Is this further sort even necessary?
!
!        IF(ToAnnihilateIndex.gt.1) THEN
!            ! Do not actually have to take indextable, hash2array or newsign 
!            ! with it...
!            call sort (ProcessVec1(1:ToAnnihilateIndex), &
!                       IndexTable1(1:ToAnnihilateIndex), &
!                       HashArray1(1:ToAnnihilateIndex), &
!                       SpawnedSign(1:ToAnnihilateIndex))
!        ENDIF
!
!!We now need to regenerate sendcounts and disps
!        sendcounts(1:nProcessors)=0
!        do i=1,ToAnnihilateIndex
!            IF(ProcessVec1(i).gt.(nProcessors-1)) THEN
!                WRITE(6,*) i,ToAnnihilateIndex
!                WRITE(6,*) "***"
!                WRITE(6,*) ProcessVec1(1:ToAnnihilateIndex)
!                WRITE(6,*) "***"
!                WRITE(6,*) sendcounts(:)
!                WRITE(6,*) "***"
!                WRITE(6,*) HashArray1(1:ToAnnihilateIndex)
!                WRITE(6,*) "***"
!                WRITE(6,*) IndexTable1(1:ToAnnihilateIndex)
!
!                CALL Stop_All("AnnihilateBetweenSpawned","Annihilation error")
!            ENDIF
!            sendcounts(ProcessVec1(i)+1)=sendcounts(ProcessVec1(i)+1)+1
!        enddo
!!The disps however do want to be cumulative
!        disps(1)=0      !Starting element is always the first element
!        do i=2,nProcessors
!            disps(i)=disps(i-1)+sendcounts(i-1)
!        enddo
!
!!We now need to calculate the recvcounts and recvdisps - this is a job for AlltoAll
!        recvcounts(1:nProcessors)=0
!
!        CALL MPIAlltoAllI(sendcounts,1,recvcounts,1,error)
!
!!We can now get recvdisps from recvcounts in the same way we obtained disps from sendcounts
!        recvdisps(1)=0
!        do i=2,nProcessors
!            recvdisps(i)=recvdisps(i-1)+recvcounts(i-1)
!        enddo
!
!        ToAnnihilateonProc=recvdisps(nProcessors)+recvcounts(nProcessors)
!
!        CALL MPIAlltoAllvI(IndexTable1(1:ToAnnihilateonProc),sendcounts,disps,IndexTable2,recvcounts,recvdisps,error)
!        CALL MPIAlltoAllvI(SpawnedSign(1:ToAnnihilateonProc),sendcounts,disps,SpawnedSign2(1:),recvcounts,recvdisps,error)
!
!        ! We now need to take with the index, the sign to remain on the entry,
!        ! as it does not necessarily want to be totally annihilated
!        call sort (IndexTable2(1:ToAnnihilateonProc), &
!                   SpawnedSign2(1:ToAnnihilateonProc))
!
!        IF(ToAnnihilateonProc.ne.0) THEN
!!Copy across the data, apart from ones which have an index given by the indicies in Index2Table(1:ToAnnihilateonProc)
!            VecSlot=1       !VecSlot is the index in the final array of TotWalkers
!            i=1             !i is the index in the original array of TotWalkersNew
!            do j=1,ToAnnihilateonProc
!!Loop over all particles to be annihilated
!                do while(i.lt.IndexTable2(j))
!!Copy accross all particles less than this number
!                    SpawnedParts2(:,VecSlot)=SpawnedParts(:,i)
!                    SpawnedSign(VecSlot)=TempSign(i)
!!                    IF(SpawnedSign(VecSlot).eq.0) THEN
!!                        CALL Stop_All("AnnihilateBetweenSpawned","Should have non-zero number of particles in this entry")
!!                    ENDIF
!                    i=i+1
!                    VecSlot=VecSlot+1
!                enddo
!                IF(SpawnedSign2(j).ne.0) THEN
!!We want the entry to be partially annihilated. Keep the particle, but change its value to be that given by SpawnedSign2
!                    SpawnedParts2(:,VecSlot)=SpawnedParts(:,i)
!                    IF(abs(SpawnedSign2(j)).ge.abs(TempSign(i))) THEN
!                        WRITE(6,*) "***",Iter,ToannihilateonProc,ValidSpawned
!                        CALL Stop_All("AnnihilateBetweenSpawned","Incorrect annihilating here...")
!                    ENDIF
!                    SpawnedSign(VecSlot)=SpawnedSign2(j)
!                    IF(SpawnedSign(VecSlot).eq.0) THEN
!                        CALL Stop_All("AnnihilateBetweenSpawned","Should have non-zero number of particles in this entry")
!                    ENDIF
!                    VecSlot=VecSlot+1
!                ENDIF
!                i=i+1
!            enddo
!
!!Now need to copy accross the residual - from Index2Table(ToAnnihilateonProc) to TotWalkersNew
!            do i=IndexTable2(ToAnnihilateonProc)+1,ValidSpawned
!                SpawnedParts2(:,VecSlot)=SpawnedParts(:,i)
!                SpawnedSign(VecSlot)=TempSign(i)
!                VecSlot=VecSlot+1
!            enddo
!
!        ELSE
!!No particles annihilated
!            VecSlot=1
!            do i=1,ValidSpawned
!                SpawnedParts2(:,VecSlot)=SpawnedParts(:,i)
!                SpawnedSign(VecSlot)=TempSign(i)
!                VecSlot=VecSlot+1
!            enddo
!        ENDIF
!
!!Have to swap arrays around here, since the pointers must stay in sync with the arrays they're pointing at.
!        ValidSpawned=VecSlot-1
!        do i=1,ValidSpawned
!            SpawnedSign2(i)=SpawnedSign(i)
!        enddo
!
!!        IF((TotWalkersNew-TotWalkers).ne.ToAnnihilateonProc) THEN
!!            WRITE(6,*) TotWalkers,TotWalkersNew,ToAnnihilateonProc,Iter
!!            CALL FLUSH(6)
!!            CALL Stop_All("AnnihilatePartPar","Problem with numbers when annihilating")
!!        ENDIF
!
!!Deallocate temp arrays
!        DEALLOCATE(TempSign)
!        DEALLOCATE(HashArray1)
!        DEALLOCATE(HashArray2)
!        DEALLOCATE(HashArrayTmp)
!        DEALLOCATE(IndexTable1)
!        DEALLOCATE(IndexTable2)
!        DEALLOCATE(ProcessVec1)
!        DEALLOCATE(ProcessVec2)
!
!!We also need to swap round the pointers to the two arrays, since the next annihilation steps take place on SpawnedParts, not SpawnedParts2 
!        IF(associated(SpawnedParts2,target=SpawnVec2)) THEN
!            SpawnedParts2 => SpawnVec
!            SpawnedSign2 => SpawnSignVec
!            SpawnedParts => SpawnVec2
!            SpawnedSign => SpawnSignVec2
!        ELSE
!            SpawnedParts => SpawnVec
!            SpawnedSign => SpawnSignVec
!            SpawnedParts2 => SpawnVec2
!            SpawnedSign2 => SpawnSignVec2
!        ENDIF
!
!        CALL halt_timer(AnnSpawned_time)
!
!    END SUBROUTINE AnnihilateBetweenSpawned
!
!
!
!        
!!This rotates the spawned (and still alive particles) around the processors. Particles are sent to MOD(iProcIndex+1,nProcessors) and received from MOD(iProcIndex+nProcessors-1,nProcessors).
!!Issues here:
!!1) Want to avoid deadlock, but also want to avoid having to send data sequentially, therefore blocking is going to be necessary.
!!2) This will also mean we have to beware of buffer overflow. Do we need to attach a specific buffer for the particles?
!!3) Do we want one of two sets of data? If two, then we need to set up a pointer system. If not, then how do we know how many particles to recieve without
!!       extra communication?
!    SUBROUTINE RotateParticles(ValidSpawned)
!        use constants, only: MpiDetInt
!        INTEGER :: error,ValidSpawned
!#ifdef PARALLEL
!        INTEGER, DIMENSION(MPI_STATUS_SIZE) :: Stat 
!        
!        CALL set_timer(Comms_Time,30)
!
!!ValidSpawned is the number of particles spawned (and still alive) for this set of particles (index is iProcIndex-no.rotates)
!        SpawnedSign(0)=ValidSpawned
!
!!        WRITE(6,*) "Particles to send: ",ValidSpawned
!!        CALL FLUSH(6)
!
!!Send the signs of the particles (number sent is in the first element)
!        CALL MPI_BSend(SpawnedSign(0:ValidSpawned),ValidSpawned+1,MPI_INTEGER,MOD(iProcIndex+1,nProcessors),123,MPI_COMM_WORLD,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            CALL Stop_All("RotateParticles","Error in sending signs")
!        ENDIF
!!        WRITE(6,*) "Sent sign",ValidSpawned+1
!!        CALL FLUSH(6)
!
!!...and then send the particles themselves...
!        CALL MPI_BSend(SpawnedParts(0:NIfTot,1:ValidSpawned),ValidSpawned*(NIfTot+1),MpiDetInt,MOD(iProcIndex+1,nProcessors),456,MPI_COMM_WORLD,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            CALL Stop_All("RotateParticles","Error in sending particles")
!        ENDIF
!!        WRITE(6,*) "Sent particles",ValidSpawned
!!        CALL FLUSH(6)
!
!!Receive signs (let it receive the maximum possible (only the first ValidSpawned will be updated.))
!        CALL MPI_Recv(SpawnedSign2(0:MaxSpawned),MaxSpawned+1,MPI_INTEGER,MOD(iProcIndex+nProcessors-1,nProcessors),123,MPI_COMM_WORLD,Stat,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            CALL Stop_All("RotateParticles","Error in receiving signs")
!        ENDIF
!!        WRITE(6,*) "Recieved sign",MaxSpawned+1
!!        CALL FLUSH(6)
!
!!Update the ValidSpawned variable for this new set of data we are about to receive...
!        ValidSpawned=SpawnedSign2(0)
!
!        CALL MPI_Recv(SpawnedParts2(0:NIfTot,1:ValidSpawned),ValidSpawned*(NIfTot+1),MpiDetInt,MOD(iProcIndex+nProcessors-1,nProcessors),456,MPI_COMM_WORLD,Stat,error)
!        IF(error.ne.MPI_SUCCESS) THEN
!            CALL Stop_All("RotateParticles","Error in receiving particles")
!        ENDIF
!!        WRITE(6,*) "Recieving particles, ",ValidSpawned
!!        CALL FLUSH(6)
!
!!We now want to make sure that we are working on the correct array. We have now received particles in SpawnedParts2 - switch it so that we are pointing at the other array.
!!We always want to annihilate from the SpawedParts and SpawnedSign arrays.
!        IF(associated(SpawnedParts2,target=SpawnVec2)) THEN
!            SpawnedParts2 => SpawnVec
!            SpawnedSign2 => SpawnSignVec
!            SpawnedParts => SpawnVec2
!            SpawnedSign => SpawnSignVec2
!        ELSE
!            SpawnedParts => SpawnVec
!            SpawnedSign => SpawnSignVec
!            SpawnedParts2 => SpawnVec2
!            SpawnedSign2 => SpawnSignVec2
!        ENDIF
!!        WRITE(6,*) "Switched arrays around..."
!!        CALL FLUSH(6)
!
!        CALL halt_timer(Comms_Time)
!
!#endif
!
!    END SUBROUTINE RotateParticles
!
!
!
END MODULE AnnihilationMod
