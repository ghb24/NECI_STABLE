#include "macros.h"

module AnnihilationMod

    use SystemData, only: NEl, tHPHF, nBasis, tCSF
    use CalcData, only: TRegenExcitgens, tEnhanceRemainder, &
                          tTruncInitiator, OccupiedThresh, tSemiStochastic, &
                          tTrialWavefunction, tKP_FCIQMC, &
                          InitiatorOccupiedThresh, tInitOccThresh
    USE DetCalcData, only: Det,FCIDetIndex
    USE Parallel_neci
    USE dSFMT_interface, only: genrand_real2_dSFMT
    USE FciMCData
    use DetBitOps, only: DetBitEQ, DetBitLT, FindBitExcitLevel, ilut_lt, &
                         ilut_gt, DetBitZero
    use DeterminantData, only: write_det
    use Determinants, only: get_helement
    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement
    use sort_mod
    use constants, only: n_int,lenof_sign,null_part,sizeof_int
    use bit_rep_data
    use bit_reps, only: decode_bit_det, extract_flags, &
                        encode_sign, encode_flags, test_flag, set_flag, &
                        clr_flag, flag_initiator, encode_part_sign, &
                        extract_part_sign, copy_flag, nullify_ilut, &
                        nullify_ilut_part, clear_has_been_initiator, &
                        set_has_been_initiator, flag_has_been_initiator
    use csf_data, only: csf_orbital_mask
    use hist_data, only: tHistSpawn, HistMinInd2
    use LoggingData, only: tNoNewRDMContrib
    use util_mod, only: get_free_unit, binary_search_custom
    use sparse_arrays, only: trial_ht, con_ht
    use global_det_data, only: set_det_diagH, get_iter_occ, &
                               global_determinant_data, set_part_init_time
    use searching
    use hash

    implicit none

    contains

    !TODO:
    !   H Elements - send through logical to decide whether to create or not.
    !   Parallel spawned parts - create the ValidSpawnedList itself.
    !   Going to have to sort this out for the new packaged walkers - will have to package them up in this interface.

    subroutine AnnihilationInterface (TotDets, MainParts, MaxMainInd, &
                                      SpawnDets, SpawnParts, MaxSpawnInd, &
                                      iter_data)

        ! This is an interface routine to the Direct Annihilation routines.
        ! It is not quite as fast as the main annihilation routines since there is a small degree of initialisation required
        ! which can be achieved on-the-fly if increased performance is required.

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

        use constants, only: size_n_int
        use shared_alloc, only: shared_allocate_iluts, shared_deallocate

        integer, intent(in) :: MaxMainInd, MaxSpawnInd
        integer, intent(inout) :: TotDets
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer(kind=n_int), intent(inout), TARGET :: MainParts(0:NIfTot,MaxMainInd),SpawnParts(0:NIfTot,MaxSpawnInd)
        integer, intent(inout) :: SpawnDets
        integer :: ierr,i
        character(len=*), parameter :: t_r = 'AnnihilationInterface'
        type(timer), save :: Annihil_time
        type(timer), save :: Sync_time
        integer(kind=n_int), pointer, save :: SpawnVecLocal(:,:)
        Sync_time%timer_name='AnnihSync   innterface'
        call set_timer(Sync_time,20)

        if (tHashWalkerList) call stop_all(t_r,"Cannot use annihilation interface with Hash table of particles")

        if (.not. (allocated(ValidSpawnedList))) then
        ! This needs to be filled correctly before annihilation can take place.
            allocate(ValidSpawnedList(0:nNodes-1),stat=ierr)
        end if

        call MPIBarrier(ierr)
        call halt_timer(Sync_time)

        Annihil_time%timer_name='Annihilation interface'
        call set_timer(Annihil_time,20)

        if (.not.(associated(SpawnVecLocal))) then
            ! This is required scratch space of the size of the spawned arrays.
            call shared_allocate_iluts("SpawnVecLocal",SpawnVecLocal,(/NIfTot,MaxSpawnInd/),iNodeIndex)
            ierr=0
            call LogMemAlloc('SpawnVecLocal',MaxSpawnInd*(NIfTot+1),size_n_int,t_r,SpawnVec2Tag,ierr)
            call MPIBarrier(ierr)
        end if

        ! ValidSpawnedList indicates the next free index for each processor.
        ! For CCMC using shared memory, we will have a single processor on each node handling this.
        ! Say there are 2 nodes with 4 processors on each.  ValidSpawnedList will contain
        ! (#Det on Node 1)
        ! (#End of Node 1's spawned list)
        ! (#End of Node 1's spawned list)
        ! (#End of Node 1's spawned list)
        ! (#Det on Node 2)
        ! (#End of Node 2's spawned list)
        ! (#End of Node 2's spawned list)
        ! (#End of Node 2's spawned list)

        ! Since we've no way of knowing about nodes as yet, we just assume all
        ! processors are on the same node.

        ! The SpawnParts already have their signs inside them

        MaxWalkersPart=MaxMainInd
        ! Point at correct arrays... will need to sort out how these are
        ! swapped in the main routine.
        CurrentDets => MainParts
        SpawnedParts => SpawnParts
        ! These point to the scratch space.
        SpawnedParts2 => SpawnVecLocal

        ! .true. for single processor annihilation.
        call DirectAnnihilation(TotDets, iter_data,.false.) 
        call MPIBarrier(ierr)

        call halt_timer(Annihil_time)

    end subroutine AnnihilationInterface

    subroutine DirectAnnihilation(TotWalkersNew, iter_data, tSingleProc)

        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: MaxIndex,ierr
        integer(kind=n_int), pointer :: PointTemp(:,:)
        logical, intent(in) :: tSingleProc
        type(timer), save :: Compress_time
        integer :: i

        ! This routine will send all the newly-spawned particles to their
        ! correct processor. 
        call SendProcNewParts(MaxIndex,tSingleProc)

        ! CompressSpawnedList works on SpawnedParts arrays, so swap the pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp

        Compress_time%timer_name='Compression interface'
        call set_timer(Compress_time,20)
    

        ! Now we want to order and compress the spawned list of particles. 
        ! This will also annihilate the newly spawned particles amongst themselves.
        ! MaxIndex will change to reflect the final number of unique determinants in the newly-spawned list, 
        ! and the particles will end up in the spawnedSign/SpawnedParts lists.
        call CompressSpawnedList(MaxIndex, iter_data)  

        call halt_timer(Compress_time)

         ! If the semi-stochastic approach is being used then the following routine performs the
         ! annihilation of the deterministic states. These states are subsequently skipped in the
         ! AnnihilateSpawnedParts routine.
         if (tSemiStochastic) call deterministic_annihilation(iter_data)

        ! Binary search the main list and copy accross/annihilate determinants which are found.
        ! This will also remove the found determinants from the spawnedparts lists.
        call AnnihilateSpawnedParts(MaxIndex,TotWalkersNew, iter_data)  

        call set_timer(Sort_Time,30)
        if (tHashWalkerList) then
            call CalcHashTableStats(TotWalkersNew, iter_data) 
        else
            ! Put the surviving particles in the main list, maintaining order of
            ! the main list (unless tHashWalkerList specified). Now we insert
            ! the remaining newly-spawned particles back into the original list
            ! (keeping it sorted),  and remove the annihilated particles from
            ! the main list.
            call InsertRemoveParts(MaxIndex, TotWalkersNew, iter_data)
        end if

        call halt_timer(Sort_Time)
        
    end subroutine DirectAnnihilation

    subroutine SendProcNewParts(MaxIndex,tSingleProc)

        ! This routine is used for sending the determinants to the correct
        ! processors. 

        integer, intent(out) :: MaxIndex
        logical, intent(in) :: tSingleProc

        integer :: i, j, error
        integer(MPIArg), dimension(nProcessors) :: sendcounts, disps, &
                                                   recvcounts, recvdisps
        integer :: MaxSendIndex
        integer(MPIArg) :: SpawnedPartsWidth
        real :: Gap

        if (tSingleProc) then
            ! Put all particles and gap on one proc.

            ! ValidSpawnedList(0:nNodes-1) indicates the next free index for each
            ! processor (for spawnees from this processor) i.e. the list of spawned
            ! particles has already been arranged so that newly spawned particles are 
            ! grouped according to the processor they go to.

            ! sendcounts(1:) indicates the number of spawnees to send to each processor.
            ! disps(1:) is the index into the spawned list of the beginning of the list
            ! to send to each processor (0-based).
           sendcounts(1)=int(ValidSpawnedList(0)-1,MPIArg)
           disps(1)=0
           if (nNodes>1) then
              sendcounts(2:nNodes)=0
              disps(2:nNodes)=int(ValidSpawnedList(1),MPIArg)
           end if
                                                                       
        else
          ! Distribute the gaps on all procs.
           do i = 0 ,nProcessors-1
               if (NodeRoots(ProcNode(i)) == i) then
                  sendcounts(i+1) = int(ValidSpawnedList(ProcNode(i)) - &
                        InitialSpawnedSlots(ProcNode(i)),MPIArg)
                  ! disps is zero-based, but InitialSpawnedSlots is 1-based.
                  disps(i+1)=int(InitialSpawnedSlots(ProcNode(i))-1,MPIArg)
               else
                  sendcounts(i+1) = 0
                  disps(i+1) = disps(i)
               end if
           end do
        end if

        MaxSendIndex = ValidSpawnedList(nNodes-1) - 1

        ! We now need to calculate the recvcounts and recvdisps - this is a
        ! job for AlltoAll
        recvcounts(1:nProcessors) = 0
        
        call MPIBarrier(error)
        call set_timer(Comms_Time,30)

        call MPIAlltoAll(sendcounts,1,recvcounts,1,error)

        ! We can now get recvdisps from recvcounts, since we want the data to
        ! be contiguous after the move.
        recvdisps(1) = 0
        do i = 2, nProcessors
            recvdisps(i) = recvdisps(i-1) + recvcounts(i-1)
        end do
        MaxIndex = recvdisps(nProcessors) + recvcounts(nProcessors)

        SpawnedPartsWidth = int(size(SpawnedParts, 1), MPIArg)
        do i = 1, nProcessors
            recvdisps(i) = recvdisps(i)*SpawnedPartsWidth
            recvcounts(i) = recvcounts(i)*SpawnedPartsWidth
            sendcounts(i) = sendcounts(i)*SpawnedPartsWidth
            disps(i) = disps(i)*SpawnedPartsWidth
        end do

        ! Max index is the largest occupied index in the array of hashes to be
        ! ordered in each processor 
        if (MaxIndex > (0.9_dp*MaxSpawned)) then
            write(6,*) MaxIndex,MaxSpawned
            call Warning_neci("SendProcNewParts","Maximum index of newly-spawned array is " &
            & //"close to maximum length after annihilation send. Increase MemoryFacSpawn")
        end if

        call MPIAlltoAllv(SpawnedParts,sendcounts,disps,SpawnedParts2,recvcounts,recvdisps,error)

        call halt_timer(Comms_Time)

    end subroutine SendProcNewParts

    subroutine CompressSpawnedList(ValidSpawned, iter_data)

        ! This sorts and compresses the spawned list to make it easier for the
        ! rest of the annihilation process. This is not essential, but should
        ! prove worthwhile.

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: VecInd,ValidSpawned,DetsMerged,i,BeginningBlockDet,FirstInitIndex,CurrentBlockDet
        real(dp) :: SpawnedSign(lenof_sign), Temp_Sign(lenof_sign)
        integer :: EndBlockDet, part_type, StartCycleInit, j, Parent_Array_Ind
        integer :: No_Spawned_Parents
        logical :: tSuc, tInc
        integer(kind=n_int), pointer :: PointTemp(:,:)
        integer(n_int) :: cum_det(0:niftot), temp_det(0:niftot)
        character(len=*), parameter :: t_r = 'CompressSpawnedList'
        type(timer), save :: Sort_time

        ! We want to sort the list of newly spawned particles, in order for
        ! quicker binary searching later on. They should remain sorted after
        ! annihilation between spawned.

        if (.not. bNodeRoot) return

        Sort_time%timer_name='Compress Sort interface'
        call set_timer(Sort_time, 20)

        call sort(SpawnedParts(:,1:ValidSpawned), ilut_lt, ilut_gt)

        call halt_timer(Sort_time)

        if (tHistSpawn) HistMinInd2(1:NEl)=FCIDetIndex(1:NEl)

        ! First, we compress the list of spawned particles, so that they are
        ! only specified at most once in each processors list. During this, we
        ! transfer the particles from SpawnedParts to SpawnedParts2. If we are
        ! working with complex walkers, we essentially do the same thing twice,
        ! annihilating real and imaginary particles seperately.

        ! This is the index in the SpawnedParts2 array to copy the compressed
        ! walkers into.
        VecInd = 1
        ! BeginningBlockDet will indicate the index of the first entry for a
        ! given determinant in SpawnedParts.
        BeginningBlockDet = 1
        DetsMerged = 0
        Parent_Array_Ind = 1
        Spawned_Parts_Zero = 0

        do while (BeginningBlockDet <= ValidSpawned)

            ! Loop in blocks of the same determinant to the end of the list of
            ! walkers.

            FirstInitIndex = 0
            CurrentBlockDet = BeginningBlockDet + 1
    
            do while (CurrentBlockDet <= ValidSpawned)
                if (.not. (DetBitEQ(SpawnedParts(0:NIfTot,BeginningBlockDet), &
                                    SpawnedParts(0:NIfTot,CurrentBlockDet),NIfDBO))) exit
                ! Loop over walkers on the same determinant in SpawnedParts.
                CurrentBlockDet = CurrentBlockDet + 1
            end do

            ! EndBlockDet indicates that we have reached the end of the block
            ! of similar dets.
            EndBlockDet = CurrentBlockDet - 1 
            
            if (EndBlockDet == BeginningBlockDet) then
                ! Optimisation: This block only consists of one entry. Simply
                !               copy it across rather than explicitly searching
                !               the list.

                ! Transfer all info to the other array.
                SpawnedParts2(:,VecInd) = SpawnedParts(:, BeginningBlockDet)   

                if (tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib)) then
                    ! SpawnedParts contains the determinants spawned on (Dj),
                    ! and it's parent (Di) plus it's sign (Cj).
                    ! As in | Dj | Di | Ci |
                    ! We then compress multiple occurances of Dj, but these may
                    ! have come from different parents, and we want to keep
                    ! track of all Di's. As we compress SpawnedParts, we
                    ! therefore move all the parents (Di's) into Spawned_Parents.
                    ! If the compressed Dj is at position VecInd in SpawnedParts,
                    ! then Spawned_Parents_Index(1,VecInd)  is the starting point
                    ! of it's parents (Di) in Spawned_Parents, and there are 
                    ! Spawned_Parents_Index(2,VecInd) entries corresponding to
                    ! this Dj.
                    if (.not. (DetBitZero(SpawnedParts(NIfTot+1:NIfTot+NIfDBO+1,BeginningBlockDet),NIfDBO))) then
                        ! If the parent determinant is null, the contribution to
                        ! the RDM is zero. No point in doing anything more with it.

                        Spawned_Parents(0:NIfDBO+1,Parent_Array_Ind) = &
                            SpawnedParts(NIfTot+1:NIfTot+NIfDBO+2,BeginningBlockDet)

                        call extract_sign (SpawnedParts(:,BeginningBlockDet), temp_sign)
                        
                        if (temp_sign(1) /= 0) then
                            ! The child (and therefore parent) are from population 1.
                            Spawned_Parents(NIfDBO+2,Parent_Array_Ind) = 1
                        else if (temp_sign(lenof_sign) /= 0) then
                            ! The child (and therefore parent) are from population 2.
                            Spawned_Parents(NIfDBO+2,Parent_Array_Ind) = lenof_sign
                        end if
                        
                        ! The first NIfDBO of the Spawned_Parents entry is the
                        ! parent determinant, the NIfDBO + 1 entry is the biased Ci.
                        ! Parent_Array_Ind keeps track of the position in
                        ! Spawned_Parents.
                        Spawned_Parents_Index(1,VecInd) = Parent_Array_Ind
                        Spawned_Parents_Index(2,VecInd) = 1
                        ! In this case there is only one instance of Dj - so
                        ! therefore only 1 parent Di.
                        Parent_Array_Ind = Parent_Array_Ind + 1
                    else
                        Spawned_Parents_Index(1,VecInd) = Parent_Array_Ind
                        Spawned_Parents_Index(2,VecInd) = 0
                    end if
                    call extract_sign (SpawnedParts(:,BeginningBlockDet), temp_sign)
                    if (IsUnoccDet(temp_sign)) then
                        Spawned_Parts_Zero = Spawned_Parts_Zero + 1
                    end if
                end if

                VecInd = VecInd + 1
                ! Move onto the next block of determinants.
                BeginningBlockDet=CurrentBlockDet 
                cycle ! Skip the rest of this block.
            end if

            ! Reset the cumulative determinant
            cum_det = 0
            cum_det (0:nifdbo) = SpawnedParts(0:nifdbo, BeginningBlockDet)
        
            if (tFillingStochRDMonFly.and.(.not.tNoNewRDmContrib)) then
                ! This is the first Dj determinant - set the index for the
                ! beginning of where the parents for this Dj can be found in
                ! Spawned_Parents.
                Spawned_Parents_Index(1,VecInd) = Parent_Array_Ind
 
                ! In this case, multiple Dj's must be compressed, and therefore
                ! the Di's dealt with as  described above. We first just
                ! initialise the position in the Spawned_Parents array to enter
                ! the Di's.
                Spawned_Parents_Index(2,VecInd) = 0
            end if

            ! Annihilate in this block seperately for real and imag walkers.
            do part_type = 1, lenof_sign   

                do i = BeginningBlockDet, EndBlockDet
                    if (tHistSpawn) then
                        call extract_sign (SpawnedParts(:,i), SpawnedSign)
                        call extract_sign (SpawnedParts(:,BeginningBlockDet), temp_sign)
                        call HistAnnihilEvent (SpawnedParts, SpawnedSign, temp_sign, part_type)
                    end if
                    call FindResidualParticle (cum_det, SpawnedParts(:,i), part_type, iter_data, &
                                                    VecInd, Parent_Array_Ind)
                end do

            end do ! Loop over particle type.


            ! Copy details into the final array.
            call extract_sign (cum_det, temp_sign)

            if ((sum(abs(temp_sign)) > 1.e-12_dp) .or. (tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib))) then
                ! Transfer all info into the other array.
                ! Usually this is only done if the final sign on the compressed
                ! Dj is not equal to zero. But in the case of the stochastic RDM,
                ! we are concerned with the sign of Dj in the CurrentDets  array,
                ! not the newly spawned sign.  We still want to check if Dj has
                ! a non-zero Cj in Current Dets, so we need to carry this Dj
                ! through to the stage of checking CurrentDets regardless of
                ! the sign here.  Also getting rid of them here would make the
                ! biased sign of Ci slightly wrong.

                SpawnedParts2(0:NIfTot,VecInd) = cum_det(0:NIfTot)
                VecInd = VecInd + 1
                DetsMerged = DetsMerged + EndBlockDet - BeginningBlockDet

                ! Spawned_Parts_Zero is the number of spawned parts that are
                ! zero after compression of the spawned_parts list - and should
                ! have been removed from SpawnedParts if we weren't calculating
                ! the RDM. - need this for a check later.
                if (IsUnoccDet(temp_sign)) Spawned_Parts_Zero = Spawned_Parts_Zero + 1
            else
                ! All particles from block have been annihilated.
                DetsMerged = DetsMerged + EndBlockDet - BeginningBlockDet + 1
                ! This spawned entry will be removed - don't want to store any
                ! parents. Reset Parent_Array_Ind so that the parents will be
                ! written over.

                ! For semi-stochastic simulations, store this state so that we
                ! can set the flags of the next state to be the same, if it is
                ! the same state.
                temp_det(0:NIfTot) = cum_det(0:NIfTot)
            end if

            ! Move onto the next block of determinants.
            BeginningBlockDet = CurrentBlockDet

        end do   

        if (tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib)) No_Spawned_Parents = Parent_Array_Ind - 1
        ! This is the new number of unique spawned determinants on the processor.
        ValidSpawned = ValidSpawned - DetsMerged    
        if (ValidSpawned /= (VecInd-1)) then
            call stop_all(t_r, "Error in compression of spawned list")
        end if

        ! Want the compressed list in spawnedparts at the end of it - swap
        ! pointers around.
        PointTemp => SpawnedParts2
        SpawnedParts2 => SpawnedParts
        SpawnedParts => PointTemp

        ! For kp-fciqmc calculations, store the current state of the spawning
        ! array.
        if (tKP_FCIQMC) then
            max_spawned_ind = ValidSpawned
            do i = 1, ValidSpawned
                SpawnedPartsKP(0:NIfDBO+lenof_sign,i) = SpawnedParts(0:NIfDBO+lenof_sign,i)
            end do
        end if

    end subroutine CompressSpawnedList

    subroutine HistAnnihilEvent(iLut, Sign1, Sign2, part_type)

        ! Histogram a possible annihilation event.

        integer(kind=n_int), intent(in) :: iLut(0:NIfTot)
        real(dp), dimension(lenof_sign), intent(in) :: Sign1,Sign2
        integer, intent(in) :: part_type
        integer :: ExcitLevel,PartIndex
        logical :: tSuc

        ! We want to histogram where the particle annihilations are taking place.
        ! No annihilation occuring - particles have the same sign. 
        if ((Sign1(part_type)*Sign2(part_type)) >= 0.0) return

        ExcitLevel = FindBitExcitLevel(iLut,iLutHF, nel)
        if (ExcitLevel == NEl) then
            call BinSearchParts2(iLut(:), HistMinInd2(ExcitLevel),Det,PartIndex,tSuc)
            HistMinInd2(ExcitLevel) = PartIndex
        else if (ExcitLevel == 0) then
            PartIndex = 1
            tSuc = .true.
        else
            call BinSearchParts2(iLut(:), HistMinInd2(ExcitLevel), FCIDetIndex(ExcitLevel+1) - 1, PartIndex, tSuc)
            HistMinInd2(ExcitLevel) = PartIndex
        end if
        if (tSuc) then
            AvAnnihil(part_type,PartIndex) = AvAnnihil(part_type,PartIndex) + &
                2*(min(abs(Sign1(part_type)), abs(Sign2(part_type))))
            InstAnnihil(part_type,PartIndex) = InstAnnihil(part_type,PartIndex) + &
                2*(min(abs(Sign1(part_type)), abs(Sign2(part_type))))
        else
            call stop_all("CompressSpawnedList", "Cannot find corresponding FCI determinant when histogramming")
        end if

    end subroutine HistAnnihilEvent

    subroutine FindResidualParticle (cum_det, new_det, part_type, iter_data, Spawned_No, Parent_Array_Ind)

        ! This routine is called whilst compressing the spawned list during
        ! annihilation. It considers the sign and flags from two particles
        ! on the same determinant, and calculates the correct sign/flags
        ! for the compressed particle.
        !
        ! --> The information is stored within the first particle in a block
        ! --> Should be called for real/imaginary particles seperately

        integer(n_int), intent(inout) :: cum_det(0:nIfTot)
        integer(n_int), intent(in) :: new_det(0:niftot+nifdbo+2)
        integer, intent(in) :: part_type, Spawned_No 
        integer, intent(inout) :: Parent_Array_Ind
        type(fcimc_iter_data), intent(inout) :: iter_data
        real(dp) :: new_sgn, cum_sgn, updated_sign, sgn_prod
        integer :: run

        ! Obtain the signs and sign product. Ignore new particle if zero.
        new_sgn = extract_part_sign (new_det, part_type)

        if (new_sgn == 0.0_dp) then
            ! New sign is just an entry from SpawnedParts - this should only ever be zero
            ! in the complex case. 
            ! If it is 0 and we're not filling the RDM (and therefore filling up the 
            ! Spawned_Parents array), can just ignore the zero entry.
            if (.not. tFillingStochRDMonFly) return
        end if

        cum_sgn = extract_part_sign (cum_det, part_type)

        ! If the cumulative and new signs for this replica are both non-zero
        ! then there have been at least two spawning events to this site, so
        ! set the initiator flag.
        ! Also set the initiator flag if the new walker has its initiator flag
        ! set.
        if ((abs(cum_sgn) > 1.e-12_dp .and. abs(new_sgn) > 1.e-12_dp) .or. &
             test_flag(new_det, flag_initiator(part_type))) &
            call set_flag(cum_det, flag_initiator(part_type))

        sgn_prod = cum_sgn * new_sgn

        ! Update annihilation statistics.
        if (sgn_prod < 0.0_dp) then
            run = part_type_to_run(part_type)
            Annihilated = Annihilated(run) + 2*min(abs(cum_sgn), abs(new_sgn))
            iter_data%nannihil(part_type) = iter_data%nannihil(part_type)&
                + 2 * min(abs(cum_sgn), abs(new_sgn))
        end if

        ! Update the cumulative sign count
        updated_sign = cum_sgn + new_sgn
        call encode_part_sign (cum_det, updated_sign, part_type)

        ! Obviously only add the parent determinant into the parent array if it is 
        ! actually being stored - and is therefore not zero.
        if (((tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib)) .and. &
            (.not. DetBitZero(new_det(NIfTot+1:NIfTot+NIfDBO+1), NIfDBO)))) then
            if (new_sgn /= 0) then
                ! No matter what the final sign is, always want to add any Di stored in 
                ! SpawnedParts to the parent array.
                Spawned_Parents(0:NIfDBO+1,Parent_Array_Ind) = new_det(NIfTot+1:NIfTot+NIfDBO+2)
                Spawned_Parents(NIfDBO+2,Parent_Array_Ind) = part_type
                Parent_Array_Ind = Parent_Array_Ind + 1
                Spawned_Parents_Index(2,Spawned_No) = Spawned_Parents_Index(2,Spawned_No) + 1
            end if
        end if

    end subroutine FindResidualParticle

    subroutine deterministic_annihilation(iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: i, j, MinInd, MaxInd, PartInd
        integer :: nI(nel)
        real(dp), dimension(lenof_sign) :: SpawnedSign, CurrentSign, SignProd
        logical :: tSuccess

        ! Copy across the weights from partial_determ_vector (the result of the deterministic projection)
        ! to CurrentDets:
        do i = 1, determ_proc_sizes(iProcIndex)
            call extract_sign(CurrentDets(:, indices_of_determ_states(i)), CurrentSign)
            SpawnedSign = partial_determ_vector(:,i)
            call encode_sign(CurrentDets(:, indices_of_determ_states(i)), SpawnedSign + CurrentSign)

            ! Update stats:
            ! Number born:
            iter_data%nborn = iter_data%nborn + abs(SpawnedSign)
            ! Number annihilated:
            SignProd = CurrentSign*SpawnedSign
            do j = 1, lenof_sign
                if (SignProd(j) < 0.0_dp) iter_data%nannihil(j) = iter_data%nannihil(j) + &
                    2*(min(abs(CurrentSign(j)), abs(SpawnedSign(j))))
            end do
        end do

    end subroutine deterministic_annihilation
    
    subroutine AnnihilateSpawnedParts(ValidSpawned,TotWalkersNew, iter_data)

        ! In this routine we want to search through the list of spawned
        ! particles. For each spawned particle, we search the list of particles
        ! in the main walker list (CurrentDets) to see if annihilation events
        ! can occur. The annihilated particles are then removed from the
        ! spawned list to the whole list of spawned particles at the end of the
        ! routine. In the main list, we change the 'sign' element of the array
        ! to zero.  These will be deleted at the end of the total annihilation
        ! step.

        use nElRDMMod, only: check_fillRDM_DiDj

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, intent(inout) :: TotWalkersNew
        integer, intent(inout) :: ValidSpawned 
        integer :: MinInd,PartInd,i,j,ToRemove,DetsMerged,PartIndex
        real(dp), dimension(lenof_sign) :: CurrentSign, SpawnedSign, SignTemp
        real(dp), dimension(lenof_sign) :: TempCurrentSign
        real(dp), dimension(lenof_sign) :: SignProd, NewSignTemp
        real(dp) :: pRemove, r
        integer :: ExcitLevel, nJ(NEl),DetHash,FinalVal,clash,walkExcitLevel, dettemp(NEl)
        integer(kind=n_int), pointer :: PointTemp(:,:)
        logical :: tSuccess,tSuc,tPrevOcc, tDetermState
        character(len=*), parameter :: t_r = "AnnihilateSpawnedParts"
        integer :: comp, run
        type(ll_node), pointer :: TempNode

        ! Only node roots to do this.
        if (.not. bNodeRoot) return

        call set_timer(AnnMain_time, 30)

        ! The number of particles to annihilate.
        ToRemove = 0
        
        ! MinInd indicates the minimum bound of the main array in which the
        ! particle can be found (and is only relevant when not using the
        ! linear scaling algorithm). Since the SpawnedParts array is ordered in
        ! the same fashion as the main array, we can find the particle position
        ! in the main array by only searching a subset.
        MinInd = 1
        PartInd = 1

        if (tHistSpawn) HistMinInd2(1:nEl) = FCIDetIndex(1:nEl)

        call set_timer(BinSearch_time,45)

        do i = 1, ValidSpawned

            ! tSuccess will determine whether the particle has been found or not.
            if (tHashWalkerList) then
                call decode_bit_det(nJ, SpawnedParts(:,i)) 
                ! Search the hash table HashIndex for the determinant defined by
                ! nJ and SpawnedParts(:,i). If it is found, tSuccess will be
                ! returned .true. and PartInd will hold the position of the
                ! determinant in CurrentDets. Else, tSuccess will be returned
                ! .false. (and PartInd shouldn't be accessed).
                ! Also, the hash value, DetHash, is returned by this routine.
                call hash_table_lookup(nJ, SpawnedParts(:,i), NIfDBO, HashIndex, CurrentDets, PartInd, DetHash, tSuccess)
            else
                ! This will binary search the CurrentDets array to find the
                ! desired particle. It will also return the index of the
                ! position one below where the particle would be found if
                ! was in the list.
                call BinSearchParts(SpawnedParts(:,i), MinInd, TotWalkersNew, PartInd, tSuccess)
            end if

            tDetermState = .false.

            if (tSuccess) then

                ! Our SpawnedParts determinant is found in CurrentDets.

                call extract_sign(CurrentDets(:,PartInd),CurrentSign)
                call extract_sign(SpawnedParts(:,i),SpawnedSign)

                SignProd = CurrentSign*SpawnedSign

                tDetermState = test_flag(CurrentDets(:,PartInd), flag_deterministic)

                if (sum(abs(CurrentSign)) /= 0.0_dp .or. tDetermState) then
                    ! Transfer new sign across.
                    call encode_sign(CurrentDets(:,PartInd), SpawnedSign+CurrentSign)
                    call encode_sign(SpawnedParts(:,i), null_part)

                    ! The only way SpawnedSign can be zero is if we are
                    ! calculating the RDM. If this is the case, we would have
                    ! already added the SpawnedDet to Spawned_Parts_Zero when it
                    ! was compressed and all walkers were annihilated. This only
                    ! counts the walkers where the SpawnedSign has newly become
                    ! zero, by merging with CurrentDets.
                    if (sum(abs(SpawnedSign)) /= 0.0_dp) ToRemove = ToRemove + 1

                    do j = 1, lenof_sign
                        run = part_type_to_run(j)
#ifndef __CMPLX
                        if (CurrentSign(j) == 0.0_dp) then
                            ! This determinant is actually *unoccupied* for the
                            ! walker type/set we're considering. We need to
                            ! decide whether to abort it or not.
                            if (tTruncInitiator) then
                                if (.not. test_flag (SpawnedParts(:,i), flag_initiator(j)) .and. &
                                     .not. tDetermState) then
                                    ! Walkers came from outside initiator space.
                                    NoAborted(j) = NoAborted(j) + abs(SpawnedSign(j))
                                    iter_data%naborted(j) = iter_data%naborted(j) + abs(SpawnedSign(j))
                                    call encode_part_sign (CurrentDets(:,PartInd), 0.0_dp, j)
                                end if
                            end if
                        else if (SignProd(j) < 0) then
#else
                        if (SignProd(j) < 0) then
#endif
                            ! This indicates that the particle has found the
                            ! same particle of opposite sign to annihilate with.
                            ! In this case we just need to update some statistics:
                            Annihilated(run) = Annihilated(run) + 2*(min(abs(CurrentSign(j)),abs(SpawnedSign(j))))
                            iter_data%nannihil(j) = iter_data%nannihil(j) + &
                                2*(min(abs(CurrentSign(j)), abs(SpawnedSign(j))))

                            if (tHistSpawn) then
                                ! We want to histogram where the particle
                                ! annihilations are taking place.
                                ExcitLevel = FindBitExcitLevel(SpawnedParts(:,i), iLutHF, nel)
                                if (ExcitLevel == NEl) then
                                    call BinSearchParts2(SpawnedParts(:,i), HistMinInd2(ExcitLevel), Det, PartIndex, tSuc)
                                else if (ExcitLevel == 0) then
                                    PartIndex = 1
                                    tSuc = .true.
                                else
                                    call BinSearchParts2(SpawnedParts(:,i), HistMinInd2(ExcitLevel), &
                                            FCIDetIndex(ExcitLevel+1)-1, PartIndex, tSuc)
                                end if
                                HistMinInd2(ExcitLevel) = PartIndex
                                if (tSuc) then
                                    AvAnnihil(j,PartIndex) = AvAnnihil(j,PartIndex)+ &
                                    real(2*(min(abs(CurrentSign(j)), abs(SpawnedSign(j)))), dp)
                                    InstAnnihil(j,PartIndex) = InstAnnihil(j,PartIndex)+ &
                                    real(2*(min(abs(CurrentSign(j)), abs(SpawnedSign(j)))), dp)
                                else
                                    write(6,*) "***",SpawnedParts(0:NIftot,i)
                                    Call WriteBitDet(6,SpawnedParts(0:NIfTot,i), .true.)
                                    call stop_all("AnnihilateSpawnedParts","Cannot find corresponding FCI "&
                                        & //"determinant when histogramming")
                                end if
                            end if
                        end if

                    end do ! Over all components of the sign.

                    if (tHashWalkerList .and. (.not. tDetermState)) then
                        call extract_sign (CurrentDets(:,PartInd), SignTemp)
                        if (IsUnoccDet(SignTemp)) then
                            ! All walkers in this main list have been annihilated
                            ! away. Remove it from the hash index array so that
                            ! no others find it (it is impossible to have another
                            ! spawned walker yet to find this determinant).
                            call remove_hash_table_entry(HashIndex, nJ, PartInd)
                            ! Add to "freeslot" list so it can be filled in.
                            iEndFreeSlot=iEndFreeSlot+1
                            FreeSlot(iEndFreeSlot)=PartInd
                        end if
                    end if

                    if (tFillingStochRDMonFly .and. (.not.tNoNewRDMContrib)) then
                        call extract_sign(CurrentDets(:,PartInd), TempCurrentSign)
                        ! We must use the instantaneous value for the off-diagonal
                        ! contribution. However, we can't just use CurrentSign from
                        ! the previous iteration, as this has been subject to death
                        ! but not the new walkers. We must add on SpawnedSign, so
                        ! we're effectively taking the instantaneous value from the
                        ! next iter. This is fine as it's from the other population,
                        ! and the Di and Dj signs are already strictly uncorrelated.
                        call check_fillRDM_DiDj(i,CurrentDets(:,PartInd),TempCurrentSign)
                    end if 

                end if

            end if
                
            if ( (.not.tSuccess) .or. (tSuccess .and. sum(abs(CurrentSign)) == 0.0_dp .and. (.not. tDetermState)) ) then

                ! Determinant in newly spawned list is not found in CurrentDets.
                ! Usually this would mean the walkers just stay in this list and
                ! get merged later - but in this case we want to check where the
                ! walkers came from - because if the newly spawned walkers are
                ! from a parent outside the active space they should be killed,
                ! as they have been spawned on an unoccupied determinant.
                if (tTruncInitiator) then

                    call extract_sign (SpawnedParts(:,i), SignTemp)

                    ! If we abort these particles, we'll still need to add them
                    ! to ToRemove.
                    tPrevOcc=.false.
                    if (.not. IsUnoccDet(SignTemp)) tPrevOcc=.true.   
                        
                    do j = 1, lenof_sign
                        run = part_type_to_run(j)
                        if ( .not. test_flag (SpawnedParts(:,i), flag_initiator(j)) ) then
                            ! If this option is on, include the walker to be
                            ! cancelled in the trial energy estimate.
                            if (tIncCancelledInitEnergy) call add_trial_energy_contrib(SpawnedParts(:,i), SignTemp(j))

                            ! Walkers came from outside initiator space.
                            NoAborted(j) = NoAborted(j) + abs(SignTemp(j))
                            iter_data%naborted(j) = iter_data%naborted(j) + abs(SignTemp(j))
                            ! We've already counted the walkers where SpawnedSign
                            ! become zero in the compress, and in the merge, all
                            ! that's left is those which get aborted which are
                            ! counted here only if the sign was not already zero
                            ! (when it already would have been counted).
                            SignTemp(j) = 0.0_dp
                            call encode_part_sign (SpawnedParts(:,i), SignTemp(j), j)

                        end if
                        
                        if (tInitOccThresh .and. test_flag(CurrentDets(:,j), flag_has_been_initiator(1)))then
                            if ((abs(SignTemp(j)) > 0.0_dp).and.(abs(SignTemp(j)) < InitiatorOccupiedThresh)) then
                                pRemove=(InitiatorOccupiedThresh-abs(SignTemp(j)))/InitiatorOccupiedThresh
                                r = genrand_real2_dSFMT ()
                                if (pRemove > r) then
                                    ! Remove the determinant.
                                    NoRemoved(run) = NoRemoved(run)+abs(SignTemp(j))
                                    iter_data%nremoved(j) = iter_data%nremoved(j) &
                                    + abs(SignTemp(j))
                                    SignTemp(j) = 0.0_dp
                                    call nullify_ilut_part(SpawnedParts(:,i),j)
                                    ! Also cancel the has_been_initiator_flag.
                                    call clear_has_been_initiator(CurrentDets(:,j),flag_has_been_initiator(1))
                                else if (tEnhanceRemainder) then
                                    NoBorn(run) = NoBorn(run) + InitiatorOccupiedThresh - abs(SignTemp(j))
                                    iter_data%nborn(j) = iter_data%nborn(j) &
                                    + InitiatorOccupiedThresh - abs(SignTemp(j))
                                    SignTemp(j) = sign(InitiatorOccupiedThresh, SignTemp(j))
                                    call encode_part_sign (SpawnedParts(:,i), SignTemp(j), j)
                                end if
                            end if
                        else
                            ! Either the determinant has never been an initiator,
                            ! or we want to treat them all the same, as before.
                            if ((abs(SignTemp(j)) > 0.0_dp) .and. (abs(SignTemp(j)) < OccupiedThresh)) then
                                ! We remove this walker with probability 1-RealSignTemp
                                pRemove=(OccupiedThresh-abs(SignTemp(j)))/OccupiedThresh
                                r = genrand_real2_dSFMT ()
                                if (pRemove > r) then
                                    ! Remove this walker.
                                    NoRemoved(run) = NoRemoved(run) + abs(SignTemp(j))
                                    !Annihilated = Annihilated + abs(SignTemp(j))
                                    !iter_data%nannihil = iter_data%nannihil + abs(SignTemp(j))
                                    iter_data%nremoved(j) = iter_data%nremoved(j) &
                                                          + abs(SignTemp(j))
                                    SignTemp(j) = 0.0_dp
                                    call nullify_ilut_part (SpawnedParts(:,i), j)
                                else if (tEnhanceRemainder) then
                                    NoBorn(run) = NoBorn(run) + OccupiedThresh - abs(SignTemp(j))
                                    iter_data%nborn(j) = iter_data%nborn(j) &
                                              + OccupiedThresh - abs(SignTemp(j))
                                    SignTemp(j) = sign(OccupiedThresh, SignTemp(j))
                                    call encode_part_sign (SpawnedParts(:,i), SignTemp(j), j)
                                end if
                            end if
                        end if
                    end do

                    if (IsUnoccDet(SignTemp) .and. tPrevOcc) then
                        ! All particle 'types' have been aborted. The zero sign
                        ! has already been taken into account in Spawned_Parts_Zero,
                        ! if it was zero directly after the compress. Only add in
                        ! here if not already taken care of there.
                         ToRemove = ToRemove + 1
                    else if (tHashWalkerList.and.(.not.(IsUnoccDet(SignTemp)))) then
                        ! Walkers have not been aborted and so we should copy the
                        ! determinant straight over to the main list. We do not
                        ! need to recompute the hash, since this should be the
                        ! same one as was generated at the beginning of the loop.
                        call AddNewHashDet(TotWalkersNew, SpawnedParts(:,i), DetHash, nJ)
                    end if

                else
                    ! Running the full, non-initiator scheme.
                    ! Determinant in newly spawned list is not found in
                    ! CurrentDets. If coeff <1, apply removal criterion.
                    call extract_sign (SpawnedParts(:,i), SignTemp)
                    
                    ! If we abort these particles, we'll still need to add
                    ! them to ToRemove.
                    tPrevOcc = .false.
                    if (.not. IsUnoccDet(SignTemp)) tPrevOcc = .true. 
                    
                    do j = 1, lenof_sign
                        run = part_type_to_run(j)
                        if ((abs(SignTemp(j)) > 0.0_dp) .and. (abs(SignTemp(j)) < OccupiedThresh)) then
                            ! We remove this walker with probability 1-RealSignTemp.
                            pRemove = (OccupiedThresh-abs(SignTemp(j)))/OccupiedThresh
                            r = genrand_real2_dSFMT ()
                            if (pRemove  >  r) then
                                ! Remove this walker.
                                NoRemoved(run) = NoRemoved(run) + abs(SignTemp(j))
                                !Annihilated = Annihilated + abs(SignTemp(j))
                                !iter_data%nannihil = iter_data%nannihil + abs(SignTemp(j))
                                iter_data%nremoved(j) = iter_data%nremoved(j) &
                                                      + abs(SignTemp(j))
                                SignTemp(j) = 0
                                call nullify_ilut_part (SpawnedParts(:,i), j)
                            else if (tEnhanceRemainder) then
                                NoBorn(run) = NoBorn(run) + OccupiedThresh - abs(SignTemp(j))
                                iter_data%nborn(j) = iter_data%nborn(j) &
                                            + OccupiedThresh - abs(SignTemp(j))
                                SignTemp(j) = sign(OccupiedThresh, SignTemp(j))
                                call encode_part_sign (SpawnedParts(:,i), SignTemp(j), j)
                            end if
                        end if
                    end do
                    
                    if (IsUnoccDet(SignTemp) .and. tPrevOcc) then
                        ! All particle 'types' have been aborted. The zero sign
                        ! has already been taken into account in Spawned_Parts_Zero,
                        ! if it was zero directly after the compress. Only add in
                        ! here if not already taken care of there.
                         ToRemove = ToRemove + 1
                    else if (tHashWalkerList .and. (.not. (IsUnoccDet(SignTemp)))) then
                        ! Walkers have not been aborted and so we should copy the
                        ! determinant straight over to the main list. We do not
                        ! need to recompute the hash, since this should be the
                        ! same one as was generated at the beginning of the loop.
                        call AddNewHashDet(TotWalkersNew, SpawnedParts(:,i), DetHash, nJ)
                    end if
                end if

                if (tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib)) then
                    ! We must use the instantaneous value for the off-diagonal contribution.
                    call check_fillRDM_DiDj(i,SpawnedParts(0:NifTot,i),SignTemp)
                end if 
            end if

            ! Even if a corresponding particle wasn't found, we can still
            ! search a smaller list next time... so not all bad news then...
            MinInd = PartInd

        end do

        call halt_timer(BinSearch_time)

        if (tHashWalkerList) then
            ! Update remaining number of holes in list for walkers stats.
            if (iStartFreeSlot > iEndFreeSlot) then
                ! All slots filled.
                HolesInList = 0
            else
                HolesInList = iEndFreeSlot - (iStartFreeSlot-1)
            end if
        else

            ! Now we have to remove the annihilated particles from the spawned
            ! list. They will be removed from the main list at the end of the
            ! annihilation process.
            if ((ToRemove + Spawned_Parts_Zero) > 0) then
            ! Since reading and writing from the same array is slow, copy the
            ! information across to the other spawned array, and just swap
            ! the pointers around after.
                DetsMerged=0
                do i = 1, ValidSpawned
                    ! We want to move all the elements above this point down to
                    ! 'fill in' the annihilated determinant.
                    call extract_sign(SpawnedParts(:,i), SignTemp)
                    if (IsUnoccDet(SignTemp)) then
                        DetsMerged = DetsMerged + 1
                    else
                        SpawnedParts2(0:NIfTot,i-DetsMerged) = SpawnedParts(0:NIfTot,i)
                    end if
                end do
                ValidSpawned = ValidSpawned - DetsMerged
                if (DetsMerged /= (ToRemove+Spawned_Parts_Zero)) then
                    write(6,*) "***", Iter, DetsMerged, ToRemove, Spawned_Parts_Zero
                    call stop_all("AnnihilateSpawnedParts", "Incorrect number of particles &
                                                             &removed from spawned list")
                end if
                ! We always want to annihilate from the SpawedParts and
                ! SpawnedSign arrays, so swap them around.
                PointTemp => SpawnedParts2
                SpawnedParts2 => SpawnedParts
                SpawnedParts => PointTemp
            end if

        end if

        call halt_timer(AnnMain_time)

    end subroutine AnnihilateSpawnedParts

    subroutine AddNewHashDet(TotWalkersNew, iLutCurr, DetHash, nJ)

        ! Add a new determinant to the main list when tHashWalkerList is true.
        ! This involves updating the list length, copying it across, updating
        ! its flag, adding its diagonal helement(if neccessary). We also need
        ! to update the hash table to point at it correctly

        integer, intent(inout) :: TotWalkersNew 
        integer(n_int), intent(inout) :: iLutCurr(0:NIfTot)
        integer, intent(in) :: DetHash, nJ(nel)
        integer :: DetPosition
        HElement_t :: HDiag
        real(dp) :: trial_amp
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

        ! If using a trial wavefunction, search to see if this state is in either the trial or
        ! connected space. If so, bin_search_trial sets the correct flag and returns the corresponding
        ! amplitude, which is stored.
        if (tTrialWavefunction) then
            if (tTrialHash) then
                call hash_search_trial(CurrentDets(:,DetPosition), nJ, trial_amp)
            else
                call bin_search_trial(CurrentDets(:,DetPosition), trial_amp)
            end if
            current_trial_amps(DetPosition) = trial_amp
        end if

        ! Add the new determinant to the hash table.
        call add_hash_table_entry(HashIndex, DetPosition, DetHash)

    end subroutine AddNewHashDet

    subroutine CalcHashTableStats(TotWalkersNew, iter_data)

        use nElRDMMod, only: det_removed_fill_diag_rdm 
        use util_mod, only: abs_sign
        use CalcData, only: tCheckHighestPop

        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: i, j, AnnihilatedDet
        real(dp) :: CurrentSign(lenof_sign), SpawnedSign(lenof_sign)
        real(dp) :: pRemove, r
        integer :: nI(nel), run
        logical :: tIsStateDeterm
        character(*), parameter :: t_r = 'CalcHashTableStats'

        if (.not. bNodeRoot) return

        TotParts = 0.0_dp
        norm_psi_squared = 0.0_dp
        norm_semistoch_squared = 0.0_dp
        iHighestPop = 0
        AnnihilatedDet = 0
        tIsStateDeterm = .false.
        InstNoAtHf = 0.0_dp

        if (TotWalkersNew > 0) then
            do i=1,TotWalkersNew
                call extract_sign(CurrentDets(:,i),CurrentSign)
                if (DetBitEQ(CurrentDets(:,i), iLutHF_True, NIfDBO)) InstNoAtHF=CurrentSign
                if (tSemiStochastic) tIsStateDeterm = test_flag(CurrentDets(:,i), flag_deterministic)

                if (IsUnoccDet(CurrentSign) .and. (.not. tIsStateDeterm)) then
                    AnnihilatedDet = AnnihilatedDet + 1 
                else
                    do j=1, lenof_sign
                        run = part_type_to_run(j)
                        if (.not. tIsStateDeterm) then
                            if (tInitOccThresh.and.test_flag(CurrentDets(:,j), flag_has_been_initiator(1)))then
                                if ((abs(CurrentSign(j)) > 0.0) .and. (abs(CurrentSign(j)) < InitiatorOccupiedThresh)) then
                                    ! We remove this walker with probability 1-RealSignTemp.
                                    pRemove = (InitiatorOccupiedThresh-abs(CurrentSign(j)))/InitiatorOccupiedThresh
                                    r = genrand_real2_dSFMT ()
                                    if (pRemove > r) then
                                        ! Remove this walker.
                                        NoRemoved(run) = NoRemoved(run) + abs(CurrentSign(j))
                                        iter_data%nremoved(j) = iter_data%nremoved(j) &
                                                              + abs(CurrentSign(j))
                                        CurrentSign(j) = 0.0_dp
                                        call nullify_ilut_part(CurrentDets(:,i), j)
                                        call decode_bit_det(nI, CurrentDets(:,i))
                                        call clear_has_been_initiator(CurrentDets(:,j),flag_has_been_initiator(1))
                                        if (IsUnoccDet(CurrentSign)) then
                                            call remove_hash_table_entry(HashIndex, nI, i)
                                            iEndFreeSlot=iEndFreeSlot+1
                                            FreeSlot(iEndFreeSlot)=i
                                        end if
                                    else if (tEnhanceRemainder) then
                                        NoBorn(run) = NoBorn(run) + InitiatorOccupiedThresh - abs(CurrentSign(j))
                                        iter_data%nborn(j) = iter_data%nborn(j) &
                                             + InitiatorOccupiedThresh - abs(CurrentSign(j))
                                        CurrentSign(j) = sign(InitiatorOccupiedThresh, CurrentSign(j))
                                        call encode_part_sign (CurrentDets(:,i), CurrentSign(j), j)
                                    end if
                                end if
                            else
                                if ((abs(CurrentSign(j)) > 0.0) .and. (abs(CurrentSign(j)) < OccupiedThresh)) then
                                !We remove this walker with probability 1-RealSignTemp
                                pRemove=(OccupiedThresh-abs(CurrentSign(j)))/OccupiedThresh
                                r = genrand_real2_dSFMT ()
                                if (pRemove  >  r) then
                                    !Remove this walker
                                    NoRemoved(run) = NoRemoved(run) + abs(CurrentSign(j))
                                    iter_data%nremoved(j) = iter_data%nremoved(j) &
                                                          + abs(CurrentSign(j))
                                    CurrentSign(j) = 0.0_dp
                                    call nullify_ilut_part(CurrentDets(:,i), j)
                                    call decode_bit_det(nI, CurrentDets(:,i))
                                    if (IsUnoccDet(CurrentSign)) then
                                        call remove_hash_table_entry(HashIndex, nI, i)
                                        iEndFreeSlot=iEndFreeSlot+1
                                        FreeSlot(iEndFreeSlot)=i
                                    end if
                                else if (tEnhanceRemainder) then
                                    NoBorn(run) = NoBorn(run) + OccupiedThresh - abs(CurrentSign(j))
                                    iter_data%nborn(j) = iter_data%nborn(j) &
                                         + OccupiedThresh - abs(CurrentSign(j))
                                    CurrentSign(j) = sign(OccupiedThresh, CurrentSign(j))
                                    call encode_part_sign (CurrentDets(:,i), CurrentSign(j), j)
                                end if
                            end if
                            end if
                            !!!!
                        end if
                    end do

                    TotParts = TotParts + abs(CurrentSign)
#if defined(__CMPLX) || defined(__DOUBLERUN)
                    norm_psi_squared = norm_psi_squared + sum(CurrentSign**2)
                    if (tIsStateDeterm) norm_semistoch_squared = norm_semistoch_squared + sum(CurrentSign**2)
#else
                    norm_psi_squared = norm_psi_squared + CurrentSign**2
                    if (tIsStateDeterm) norm_semistoch_squared = norm_semistoch_squared + CurrentSign**2
#endif
                    
                    if (tCheckHighestPop) then
                        ! If this option is on, then we want to compare the 
                        ! weight on each determinant to the weight at the HF 
                        ! determinant.
                        !
                        ! Record the highest weighted determinant on each 
                        ! processor. If double run, only consider set 1 to keep things simple.
                        if (abs_sign(ceiling(CurrentSign)) > iHighestPop) then
                            iHighestPop = int(abs_sign(ceiling(CurrentSign)))
                            HighestPopDet(:)=CurrentDets(:,i)
                        end if
                    end if
                end if

                if (IsUnoccDet(CurrentSign) .and. (.not. tIsStateDeterm)) then
                    if (DetBitEQ(CurrentDets(:,i), iLutHF_True, NIfDBO)) then
                        ! We have to do this such that AvNoAtHF matches up with AvSign.
                        ! AvSign is extracted from CurrentH, and if the HFDet is unoccupied
                        ! at this moment during annihilation, it's CurrentH entry is removed
                        ! and the averaging information in it is lost.
                        ! In some cases (a successful spawning event) a CurrentH entry will
                        ! be recreated, but with AvSign 0, so we must match this here.
                        AvNoAtHF = 0.0_dp 
                        IterRDM_HF = Iter + 1 
                    end if
                end if

                if (tFillingStochRDMonFly .and. (.not. tIsStateDeterm)) then
                    if (inum_runs == 2) then

                        if ((CurrentSign(1) == 0 .and. get_iter_occ(i, 1) /= 0) .or. &
                            (CurrentSign(inum_runs) == 0 .and. get_iter_occ(i, 2) /= 0) .or. &
                            (CurrentSign(1) /= 0 .and. get_iter_occ(i, 1) == 0) .or. &
                            (CurrentSign(inum_runs) /= 0 .and. get_iter_occ(i, 2) == 0)) then
                               
                            ! At least one of the signs has just gone to zero or just become reoccupied
                            ! so we need to consider adding in diagonal elements and connections to HF
                            ! The block that's just ended was occupied in at least one population.
                            call det_removed_fill_diag_rdm(CurrentDets(:,i), i)
                        end if
                    else
                        if (IsUnoccDet(CurrentSign)) then
                            call det_removed_fill_diag_rdm(CurrentDets(:,i), i)
                        end if
                    end if
                end if

                if (IsUnoccDet(CurrentSign) .and. (.not. tIsStateDeterm) .and. tTruncInitiator) then
                    do j=1,lenof_sign
                        if (test_flag(CurrentDets(:,i),flag_initiator(j))) then
                            !determinant was an initiator...it obviously isn't any more...
                            NoAddedInitiators(j)=NoAddedInitiators(j)-1
                        end if
                    end do
                end if

            end do
        end if

        if (AnnihilatedDet /= HolesInList) then
            write(6,*) "TotWalkersNew: ",TotWalkersNew
            write(6,*) "AnnihilatedDet: ",AnnihilatedDet
            write(6,*) "HolesInList: ",HolesInList
            call stop_all(t_r, "Error in determining annihilated determinants")
        end if

    end subroutine CalcHashTableStats
    
    subroutine InsertRemoveParts(ValidSpawned, TotWalkersNew, iter_data)

        ! This routine will run through the total list of particles
        ! (TotWalkersNew in CurrentDets with sign CurrentSign) and the list of
        ! newly-spawned but non-annihilated particles (ValidSpawned in
        ! SpawnedParts and SpawnedSign) and move the  new particles into the
        ! correct place in the new list, while removing the particles with
        ! sign = 0 from CurrentDets.

        use util_mod, only: abs_sign
        use SystemData, only: tHPHF, tRef_Not_HF
        use bit_reps, only: NIfD
        use LoggingData, only: tRDMonFly, tExplicitAllRDM
        use nElRDMMod, only: det_removed_fill_diag_rdm 
        use CalcData, only: tCheckHighestPop, NMCyc

        integer, intent(in) :: ValidSpawned
        integer, intent(inout) :: TotWalkersNew
        real(dp) :: CurrentSign(lenof_sign), SpawnedSign(lenof_sign)
        real(dp) :: HDiag, pRemove, r
        integer :: i,DetsMerged,nJ(NEl),part_type, ExcitLevelCurr, j, run
        integer :: trial_merged, con_merged, i_trial, i_conn
        logical :: TestClosedShellDet
        logical :: tIsStateDeterm, tTrialState, tConState
        type(fcimc_iter_data), intent(inout) :: iter_data
        character(*), parameter :: t_r = 'InsertRemoveParts'
        HElement_t :: HDiagTemp

        ! It appears that the rest of this routine isn't thread-safe if ValidSpawned is zero.

        if (.not. bNodeRoot) return
            
        ! This logical is only used for the semi-stochastic code. If true then we don't
        ! try to remove the state.
        tIsStateDeterm = .false.

        ! Annihilated determinants first are removed from the main array (zero sign). 
        ! Surely we only need to perform this loop if the number of annihilated particles > 0?
        TotParts = 0.0
        norm_psi_squared = 0.0_dp
        norm_semistoch_squared = 0.0_dp
        DetsMerged = 0
        trial_merged = 0
        i_trial = 0
        con_merged = 0
        i_conn = 0
        iHighestPop = 0
        InstNoatHF = 0.0_dp

        if (TotWalkersNew > 0) then
            do i=1,TotWalkersNew
                call extract_sign(CurrentDets(:,i),CurrentSign)
                if (tSemiStochastic) tIsStateDeterm = test_flag(CurrentDets(:,i), flag_deterministic)
                do j=1, lenof_sign
                    run = part_type_to_run(j)
                    if (.not. tIsStateDeterm) then
                        if (tInitOccThresh .and. test_flag(CurrentDets(:,j),flag_has_been_initiator(1)))then
                            if ((abs(CurrentSign(j)) > 0.0) .and. (abs(CurrentSign(j)) < InitiatorOccupiedThresh)) then
                                ! We remove this walker with probability 1-RealSignTemp.
                                pRemove=(InitiatorOccupiedThresh-abs(CurrentSign(j)))/InitiatorOccupiedThresh
                                r = genrand_real2_dSFMT ()
                                if (pRemove  >  r) then
                                    ! Remove this walker.
                                    NoRemoved(run) = NoRemoved(run) + abs(CurrentSign(j))
                                    !Annihilated = Annihilated + abs(CurrentSign(j))
                                    !iter_data%nannihil = iter_data%nannihil + abs(CurrentSign(j))
                                    iter_data%nremoved(j) = iter_data%nremoved(j) &
                                                          + abs(CurrentSign(j))
                                    CurrentSign(j) = 0
                                    call nullify_ilut_part (CurrentDets(:,i), j)
                                    call clear_has_been_initiator(CurrentDets(:,j),flag_has_been_initiator(1))
                                else if (tEnhanceRemainder) then
                                    ! SDS: TODO: Account for the TotParts Changes
                                    NoBorn(run) = NoBorn(run) + InitiatorOccupiedThresh - abs(CurrentSign(j))
                                    iter_data%nborn(j) = iter_data%nborn(j) &
                                             + InitiatorOccupiedThresh - abs(CurrentSign(j))
                                    CurrentSign(j) = sign(InitiatorOccupiedThresh, CurrentSign(j))
                                    call encode_part_sign (CurrentDets(:,i), CurrentSign(j), j)
                                end if
                            end if
                        else
                            if ((abs(CurrentSign(j)) > 0.0) .and. (abs(CurrentSign(j)) < OccupiedThresh)) then
                                ! We remove this walker with probability 1-RealSignTemp.
                                pRemove=(OccupiedThresh-abs(CurrentSign(j)))/OccupiedThresh
                                r = genrand_real2_dSFMT ()
                                if (pRemove  >  r) then
                                    ! Remove this walker.
                                    NoRemoved(run) = NoRemoved(run) + abs(CurrentSign(j))
                                    !Annihilated = Annihilated + abs(CurrentSign(j))
                                    !iter_data%nannihil = iter_data%nannihil + abs(CurrentSign(j))
                                    iter_data%nremoved(j) = iter_data%nremoved(j) &
                                                          + abs(CurrentSign(j))
                                    CurrentSign(j) = 0
                                    call nullify_ilut_part (CurrentDets(:,i), j)
                                else if (tEnhanceRemainder) then
                                    ! SDS: TODO: Account for the TotParts Changes
                                    NoBorn(run) = NoBorn(run) + OccupiedThresh - abs(CurrentSign(j))
                                    iter_data%nborn(j) = iter_data%nborn(j) &
                                             + OccupiedThresh - abs(CurrentSign(j))
                                    CurrentSign(j) = sign(OccupiedThresh, CurrentSign(j))
                                    call encode_part_sign (CurrentDets(:,i), CurrentSign(j), j)
                                end if
                            end if
                        end if
                    end if
                end do
                
                if (DetBitEQ(CurrentDets(:,i), iLutHF_True, NIfDBO)) &
                    InstNoAtHF(1:lenof_sign) = CurrentSign

                ! Is this state a trial or connected state, or neither?
                if (tTrialWavefunction) then
                    if (test_flag(CurrentDets(:,i), flag_trial)) then
                        i_trial = i_trial + 1
                        tTrialState = .true.
                        tConState = .false.
                    else if (test_flag(CurrentDets(:,i), flag_connected)) then
                        i_conn = i_conn + 1
                        tConState = .true.
                        tTrialState = .false.
                    else
                        tTrialState = .false.
                        tConState = .false.
                    end if
                end if
                
                if (tFillingStochRDMonFly .and. (.not. tIsStateDeterm)) then
                    if (inum_runs == 2) then

                        if ((CurrentSign(1) == 0 .and. get_iter_occ(i, 1) /= 0) .or. &
                            (CurrentSign(inum_runs) == 0 .and. get_iter_occ(i, 2) /= 0) .or. &
                            (CurrentSign(1) /= 0 .and. get_iter_occ(i, 1) == 0) .or. &
                            (CurrentSign(inum_runs) /= 0 .and. get_iter_occ(i, 2) == 0)) then
                        
                            ! At least one of the signs has just gone to zero or just become reoccupied
                            ! so we need to consider adding in diagonal elements and connections to HF
                            ! The block that's just ended was occupied in at least one population.
                            call det_removed_fill_diag_rdm(CurrentDets(:,i), i)
                        end if
                    else
                        if (IsUnoccDet(CurrentSign)) then
                            call det_removed_fill_diag_rdm(CurrentDets(:,i), i)
                        end if
                    end if
                end if

                if (IsUnoccDet(CurrentSign) .and. (.not. tIsStateDeterm)) then

                    DetsMerged=DetsMerged+1
                    if (tTrialWavefunction) then
                        if (tTrialState) trial_merged = trial_merged + 1
                        if (tConState) con_merged = con_merged + 1
                    end if

                    if (i == HFInd) then
                        ! We have to do this such that AvNoAtHF matches up with
                        ! AvSign. AvSign is extracted from global_determinant_data,
                        ! and if the HFDet is unoccupied at this moment during
                        ! annihilation, it's global_determinant_data entry is removed
                        ! and the averaging information in it is lost. In some
                        ! cases (a successful spawning event) a global_determinant_data
                        ! entry will be recreated, but with AvSign 0, so we must match
                        ! this here.
                        AvNoAtHF=0.0_dp 
                        IterRDM_HF = Iter + 1 
                    end if
                    if (tTruncInitiator) then
                        do part_type=1,lenof_sign
                            if (test_flag(CurrentDets(:,i),flag_initiator(part_type))) then
                                ! Determinant was an initiator...
                                ! It obviously isn't any more...
                                NoAddedInitiators(part_type)=NoAddedInitiators(part_type)-1
                            end if
                        end do
                    end if
                else
                    ! We want to move all the elements above this point down
                    ! to 'fill in' the annihilated determinant.
                    if (DetsMerged /= 0) then
                        CurrentDets(0:NIfTot,i-DetsMerged)=CurrentDets(0:NIfTot,i)
                        ! It isn't normally acceptable to access this array
                        ! directly (it has accessor functions). But we are
                        ! merely shuffling data in blocks, so its OK.
                        global_determinant_data(:, i - DetsMerged) = &
                            global_determinant_data(:, i)

                        ! Move the elements in the occupied trial and connected
                        ! vectors to fill in the values of the annihilate
                        ! determinants.
                        if (tTrialWavefunction) then
                            if (tTrialState) occ_trial_amps(i_trial-trial_merged) = occ_trial_amps(i_trial)
                            if (tConState) occ_con_amps(i_conn-con_merged) = occ_con_amps(i_conn)
                        end if

                    end if

                    TotParts(1:lenof_sign) = TotParts(1:lenof_sign) &
                                           + abs(CurrentSign)
#ifdef __CMPLX
                    norm_psi_squared = norm_psi_squared + sum(CurrentSign**2)
                    if (tIsStateDeterm) norm_semistoch_squared = norm_semistoch_squared + sum(CurrentSign**2)
#else
                    norm_psi_squared(1:inum_runs) = &
                            norm_psi_squared(1:inum_runs) + CurrentSign**2
                    if (tIsStateDeterm) norm_semistoch_squared(1:inum_runs) = &
                        norm_semistoch_squared(1:inum_runs) + CurrentSign**2
#endif
                    
                    ! If this option is on, then we want to compare the weight
                    ! on each determinant to the weight at the HF determinant.
                    ! Record the highest weighted determinant on each processor.
                    if (tCheckHighestPop) then
                        if (abs_sign(ceiling(CurrentSign)) > iHighestPop) then
                            iHighestPop = int(abs_sign(ceiling(CurrentSign)))
                            HighestPopDet(:) = CurrentDets(:,i)
                        end if
                    end if
                end if
            end do
            TotWalkersNew=TotWalkersNew-DetsMerged
            ntrial_occ = ntrial_occ-trial_merged
            ncon_occ = ncon_occ-con_merged
        end if

        ! We now calculate the contribution to the total number of particles
        ! from the spawned lists. The list has previously been compressed.
        if (ValidSpawned > 0) then
            call extract_sign(SpawnedParts(:,1),SpawnedSign)
            TotParts(1:lenof_sign) = TotParts(1:lenof_sign) + abs(SpawnedSign)
#ifdef __CMPLX
            norm_psi_squared = norm_psi_squared + sum(SpawnedSign**2)
#else
            norm_psi_squared(1:inum_runs) = norm_psi_squared(1:inum_runs) &
                                          + SpawnedSign**2
#endif

        end if
        do i=2,ValidSpawned
            call extract_sign(SpawnedParts(:,i),SpawnedSign)
            TotParts(1:lenof_sign) = TotParts(1:lenof_sign) + abs(SpawnedSign)
#ifdef __CMPLX
            norm_psi_squared = norm_psi_squared + sum(SpawnedSign**2)
#else
            norm_psi_squared(1:inum_runs) = &
                         norm_psi_squared(1:inum_runs) + SpawnedSign**2
#endif
        end do

        ! TotWalkersNew is now the number of determinants in the main list left.
        ! We now want to merge the main list with the spawned list of
        ! non-annihilated spawned particles. The final list will be of length
        ! TotWalkersNew+ValidSpawned. This will be returned in the first element
        ! of MergeLists updated.
        if (TotWalkersNew + ValidSpawned > MaxWalkersPart) then
            write(6,*) "Non-annihilated old walkers:", TotWalkersNew
            write(6,*) "Non-annihilated spawned:", ValidSpawned
            write(6,*) "Total walkers to remain:", TotWalkersNew + ValidSpawned
            write(6,*) "Size of Particle List:", MaxWalkersPart
            call stop_all(t_r, "Not enough space in particle list for merge.") 
        end if
        if (TotWalkersNew == 0) then
            ! Merging algorithm will not work with no determinants in the main list.
            TotWalkersNew = ValidSpawned

            ! If using a trial wavefunction then call a routine to set the trial and connected
            ! flags for the determinants in SpawnedParts, where necessary, and to also count and
            ! return the number of trial and connected determinants. The corresponding trial
            ! and connected vector amplitudes are stored in trial_temp and con_temp on output.
            ! These are then copied across to occ_trial_amps and con_trial_amps.
            if (tTrialWavefunction) then
                if (tTrialHash) then
                    call find_trial_and_con_states_hash(int(ValidSpawned,8), &
                                SpawnedParts(0:NIfTot,1:ValidSpawned), ntrial_occ, ncon_occ)
                else
                    call find_trial_and_con_states_bin(int(ValidSpawned,8), &
                                SpawnedParts(0:NIfTot,1:ValidSpawned), ntrial_occ, ncon_occ)
                end if
                occ_trial_amps(1:ntrial_occ) = trial_temp(1:ntrial_occ)
                occ_con_amps(1:ncon_occ) = con_temp(1:ncon_occ)
            end if

            do i = 1, ValidSpawned
                CurrentDets(:,i) = SpawnedParts(:,i)
                ! We want to calculate the diagonal hamiltonian matrix element for the new particle to be merged.
                if (DetBitEQ(CurrentDets(:,i), iLutHF_True, NIfDBO)) call extract_sign(CurrentDets(:,i),InstNoatHF)
                if (DetBitEQ(CurrentDets(:,i), iLutRef, NIfDBO)) then
                    ! We know we are at HF - HDiag=0
                    HDiag = 0.0_dp
                else
                    call decode_bit_det (nJ, CurrentDets(:,i))
                    if (tHPHF) then
                        HDiagTemp = hphf_diag_helement (nJ, &
                                                        CurrentDets(:,i))
                    else
                        HDiagTemp = get_helement (nJ, nJ, 0)
                    end if
                    HDiag = (real(HDiagTemp,dp)) - Hii
                end if
                call set_det_diagH(i, HDiag)

                ! Store the iteration this particle is being created on
                call set_part_init_time(i, TotImagTime)

            end do
        else
            call MergeListswH(TotWalkersNew,ValidSpawned,SpawnedParts(0:NIfTot,1:ValidSpawned))
        end if

    end subroutine InsertRemoveParts

end module AnnihilationMod
