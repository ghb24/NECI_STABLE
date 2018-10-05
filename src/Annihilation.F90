#include "macros.h"

module AnnihilationMod

    use SystemData, only: NEl, tHPHF
    use CalcData, only:   tTruncInitiator, OccupiedThresh, tSemiStochastic, &
                          tTrialWavefunction, tKP_FCIQMC, tContTimeFCIMC, &
                          tContTimeFull, InitiatorWalkNo, tau, tEN2, tEN2Init, &
                          tEN2Started, tEN2Truncated, tInitCoherentRule, t_truncate_spawns, &
                          n_truncate_spawns, t_prone_walkers, t_truncate_unocc, &
                          tSpawnSeniorityBased, numMaxExLvlsSet, maxKeepExLvl
    use DetCalcData, only: Det, FCIDetIndex
    use Parallel_neci
    use dSFMT_interface, only: genrand_real2_dSFMT
    use FciMCData
    use DetBitOps, only: DetBitEQ, FindBitExcitLevel, ilut_lt, &
                         ilut_gt, DetBitZero, count_open_orbs
    use sort_mod
    use constants, only: n_int, lenof_sign, null_part, sizeof_int
    use bit_rep_data
    use bit_reps, only: decode_bit_det, &
                        encode_sign, test_flag, set_flag, &
                        flag_initiator, encode_part_sign, &
                        extract_part_sign, extract_bit_rep, &
                        nullify_ilut_part, clr_flag, get_num_spawns,&
                        encode_flags, bit_parent_zero, get_initiator_flag
    use hist_data, only: tHistSpawn, HistMinInd2
    use LoggingData, only: tNoNewRDMContrib
    use load_balance, only: DetermineDetNode, AddNewHashDet, &
                            CalcHashTableStats, get_diagonal_matel
    use searching
    use hash
    use global_det_data, only: det_diagH
    use procedure_pointers, only: scaleFunction
    use Determinants, only: get_helement
    use hphf_integrals, only: hphf_diag_helement
    use rdm_data, only: rdm_estimates, en_pert_main
    use rdm_data_utils, only: add_to_en_pert_t
    use fcimc_helper, only: CheckAllowedTruncSpawn

    implicit none

    contains

    subroutine DirectAnnihilation(TotWalkersNew, iter_data, tSingleProc)

        integer, intent(inout) :: TotWalkersNew
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: MaxIndex
        integer(kind=n_int), pointer :: PointTemp(:,:)
        logical, intent(in) :: tSingleProc
        type(timer), save :: Compress_time

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

        CALL set_timer(Sort_Time,30)
        call CalcHashTableStats(TotWalkersNew, iter_data) 
        CALL halt_timer(Sort_Time)

    end subroutine DirectAnnihilation

    subroutine SendProcNewParts(MaxIndex,tSingleProc)

        ! This routine is used for sending the determinants to the correct
        ! processors. 

        integer, intent(out) :: MaxIndex
        logical, intent(in) :: tSingleProc

        integer :: i, error
        integer(MPIArg), dimension(nProcessors) :: sendcounts, disps, &
                                                   recvcounts, recvdisps
        integer :: MaxSendIndex
        integer(MPIArg) :: SpawnedPartsWidth

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
              ! n.b. work around PGI bug.
              do i = 2, nNodes
                  disps(i) = int(ValidSpawnedList(1), MPIArg)
              end do
              !disps(2:nNodes)=int(ValidSpawnedList(1),MPIArg)
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

        ! Set this global data - the total number of spawned determants.
        nspawned = sum(recvcounts)

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
#ifdef __DEBUG
            write(6,*) MaxIndex,MaxSpawned
#else
            write(*,*) 'On task ',iProcIndex,': ',MaxIndex,MaxSpawned
#endif
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
        integer :: EndBlockDet, part_type, Parent_Array_Ind
        integer :: No_Spawned_Parents
        integer(kind=n_int), pointer :: PointTemp(:,:)
        integer(n_int) :: cum_det(0:nifbcast), temp_det(0:nifbcast)
        character(len=*), parameter :: t_r = 'CompressSpawnedList'
        type(timer), save :: Sort_time

        ! We want to sort the list of newly spawned particles, in order for
        ! quicker binary searching later on. They should remain sorted after
        ! annihilation between spawned.

        if (.not. bNodeRoot) return

        Sort_time%timer_name='Compress Sort interface'
        call set_timer(Sort_time, 20)

        call sort(SpawnedParts(0:NIfBCast,1:ValidSpawned), ilut_lt, ilut_gt)

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
                if (.not. (DetBitEQ(SpawnedParts(:, BeginningBlockDet), &
                                    SpawnedParts(:, CurrentBlockDet)))) exit
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

                ! If this one entry has no amplitude then don't add it to the
                ! compressed list, but just cycle.
                call extract_sign (SpawnedParts(:, BeginningBlockDet), temp_sign)
                if ( (sum(abs(temp_sign)) < 1.e-12_dp) .and. (.not. (tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib))) ) then
                    DetsMerged = DetsMerged + 1
                    BeginningBlockDet = CurrentBlockDet 
                    cycle
                end if

                ! Transfer all info to the other array.
                SpawnedParts2(:, VecInd) = SpawnedParts(:, BeginningBlockDet)

                if (tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib)) then
                    ! SpawnedParts contains the determinants spawned on (Dj),
                    ! and it's parent (Di) plus it's sign (Cj).
                    ! As in | Dj | Di | Ci |
                    ! We then compress multiple occurances of Dj, but these may
                    ! have come from different parents, and we want to keep
                    ! track of all Di's. As we compress SpawnedParts, we
                    ! therefore move all the parents (Di's) into Spawned_Parents.
                    ! If the compressed Dj is at position VecInd in SpawnedParts,
                    ! then Spawned_Parents_Index(1,VecInd) is the starting point
                    ! of it's parents (Di) in Spawned_Parents, and there are 
                    ! Spawned_Parents_Index(2,VecInd) entries corresponding to
                    ! this Dj.

                    if (.not. bit_parent_zero(SpawnedParts(:, BeginningBlockDet))) then
                    
                        ! If the parent determinant is null, the contribution to
                        ! the RDM is zero. No point in doing anything more with it.

                        ! Why is this length nifdbo+2? What is the extra bit?

                        Spawned_Parents(0:NIfDBO+1,Parent_Array_Ind) = &
                            SpawnedParts(nOffParent:nOffParent+nIfDBO+1, BeginningBlockDet)

                        call extract_sign (SpawnedParts(:,BeginningBlockDet), temp_sign)
                        
                        ! Search to see which sign is non-zero, and therefore
                        ! find which simulation the spawning occured from and to.
                        ! NOTE: it is safe to compare against zero exactly here,
                        ! because all other components will have been set to zero
                        ! exactly and can't have changed at all.
                        Spawned_Parents(NIfDBO+2,Parent_Array_Ind) = 0
                        do part_type = 1, lenof_sign
                            if (abs(temp_sign(part_type)) > 1.0e-12_dp) then
                                Spawned_Parents(NIfDBO+2,Parent_Array_Ind) = part_type
                                exit
                            end if
                        end do
                        
                        ! The first NIfDBO of the Spawned_Parents entry is the
                        ! parent determinant, the NIfDBO + 1 entry is the Ci.
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
                BeginningBlockDet = CurrentBlockDet 
                cycle ! Skip the rest of this block.
            end if

            ! Reset the cumulative determinant
            cum_det = 0_n_int
            cum_det (0:nifdbo) = SpawnedParts(0:nifdbo, BeginningBlockDet)
        
            if (tFillingStochRDMonFly .and. (.not.tNoNewRDMContrib)) then
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


            do i = BeginningBlockDet, EndBlockDet
               ! if logged, accumulate the number of spawn events
               if(tLogNumSpawns) then
                  cum_det(nSpawnOffset) = cum_det(nSpawnOffset) + &
                       SpawnedParts(nSpawnOffset,i)
               end if
               ! Annihilate in this block seperately for walkers of different types.
               do part_type = 1, lenof_sign
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
                ! we are concerned with the sign of Dj in the CurrentDets array,
                ! not the newly spawned sign.  We still want to check if Dj has
                ! a non-zero Cj in CurrentDets, so we need to carry this Dj
                ! through to the stage of checking CurrentDets regardless of
                ! the sign here.  Also getting rid of them here would make the
                ! biased sign of Ci slightly wrong.

                SpawnedParts2(0:NIfTot,VecInd) = cum_det(0:NIfTot)
                VecInd = VecInd + 1
                DetsMerged = DetsMerged + EndBlockDet - BeginningBlockDet

                ! Spawned_Parts_Zero is the number of spawned parts that are
                ! zero after compression of the spawned_parts list - and should
                ! have been removed from SpawnedParts if we weren't calculating
                ! the RDM - need this for a check later.
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

    end subroutine CompressSpawnedList

    subroutine HistAnnihilEvent(iLut, Sign1, Sign2, part_type)

        ! Histogram a possible annihilation event.

        integer(kind=n_int), intent(in) :: iLut(0:NIfTot)
        real(dp), dimension(lenof_sign), intent(in) :: Sign1, Sign2
        integer, intent(in) :: part_type
        integer :: ExcitLevel, PartIndex
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

        new_sgn = extract_part_sign (new_det, part_type)
        cum_sgn = extract_part_sign (cum_det, part_type)

        ! If the cumulative and new signs for this replica are both non-zero
        ! then there have been at least two spawning events to this site, so
        ! set the initiator flag.
        ! (There is now an option (tInitCoherentRule = .false.) to turn this
        ! coherent spawning rule off, mainly for testing purposes).
        ! Also set the initiator flag if the new walker has its initiator flag
        ! set.
        if (tTruncInitiator) then
            if (tInitCoherentRule) then
                if ((abs(cum_sgn) > 1.e-12_dp .and. abs(new_sgn) > 1.e-12_dp) .or. &
                     test_flag(new_det, get_initiator_flag(part_type))) &
                    call set_flag(cum_det, get_initiator_flag(part_type))
            else
                if (test_flag(new_det, get_initiator_flag(part_type))) &
                    call set_flag(cum_det, get_initiator_flag(part_type))
            end if
        end if

        ! set the multi-spawn flag if there has been more than one spawn
        if(abs(cum_sgn) > eps .and. abs(new_sgn) > eps) call set_flag(cum_det, flag_multi_spawn)

        sgn_prod = cum_sgn * new_sgn

        ! Update annihilation statistics.
        if (sgn_prod < 0.0_dp) then
            run = part_type_to_run(part_type)
            Annihilated(run) = Annihilated(run) + 2*min(abs(cum_sgn), abs(new_sgn))
            iter_data%nannihil(part_type) = iter_data%nannihil(part_type)&
                + 2 * min(abs(cum_sgn), abs(new_sgn))
        end if

        ! Update the cumulative sign count.
        updated_sign = cum_sgn + new_sgn
        call encode_part_sign (cum_det, updated_sign, part_type)

        ! Obviously only add the parent determinant into the parent array if it is 
        ! actually being stored - and is therefore not zero.
        if (((tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib)) .and. &
            (.not. DetBitZero(new_det(NIfTot+1:NIfTot+NIfDBO+1), NIfDBO)))) then
            if (abs(new_sgn) > 1.e-12_dp) then
                ! Add parent (Di) stored in SpawnedParts to the parent array.
                Spawned_Parents(0:NIfDBO+1,Parent_Array_Ind) = new_det(NIfTot+1:NIfTot+NIfDBO+2)
                Spawned_Parents(NIfDBO+2,Parent_Array_Ind) = part_type
                Parent_Array_Ind = Parent_Array_Ind + 1
                Spawned_Parents_Index(2,Spawned_No) = Spawned_Parents_Index(2,Spawned_No) + 1
            end if
        end if

    end subroutine FindResidualParticle

    subroutine deterministic_annihilation(iter_data)

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer :: i, j
        real(dp), dimension(lenof_sign) :: SpawnedSign, CurrentSign, SignProd

        ! Copy across the weights from partial_determ_vecs (the result of the deterministic projection)
        ! to CurrentDets:
        do i = 1, determ_sizes(iProcIndex)
            call extract_sign(CurrentDets(:, indices_of_determ_states(i)), CurrentSign)
            SpawnedSign = partial_determ_vecs(:,i)
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
    
    subroutine AnnihilateSpawnedParts(ValidSpawned, TotWalkersNew, iter_data)

        ! In this routine we want to search through the list of spawned
        ! particles. For each spawned particle, we search the list of particles
        ! in the main walker list (CurrentDets) to see if annihilation events
        ! can occur. The annihilated particles are then removed from the
        ! spawned list to the whole list of spawned particles at the end of the
        ! routine. In the main list, we change the 'sign' element of the array
        ! to zero.  These will be deleted at the end of the total annihilation
        ! step.

        use LoggingData, only: tOldRDMs
        use rdm_data_old, only: rdms, one_rdms_old
        use rdm_data, only: rdm_definitions, two_rdm_spawn, one_rdms
        use rdm_filling_old, only: check_fillRDM_DiDj_old
        use rdm_filling, only: check_fillRDM_DiDj

        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, intent(inout) :: TotWalkersNew
        integer, intent(inout) :: ValidSpawned 
        integer :: PartInd, i, j, PartIndex, run
        real(dp), dimension(lenof_sign) :: CurrentSign, SpawnedSign, SignTemp
        real(dp), dimension(lenof_sign) :: TempCurrentSign, SignProd
        real(dp) :: ScaledOccupiedThresh, scFVal, diagH
        integer :: ExcitLevel, DetHash, nJ(nel)
        logical :: tSuccess, tSuc, tDetermState
        logical :: abort(lenof_sign)
        logical :: tTruncSpawn, t_truncate_this_det

        ! Only node roots to do this.
        if (.not. bNodeRoot) return

        call set_timer(AnnMain_time, 30)

        if (tHistSpawn) HistMinInd2(1:nEl) = FCIDetIndex(1:nEl)

        call set_timer(BinSearch_time,45)

        do i = 1, ValidSpawned

            call decode_bit_det(nJ, SpawnedParts(:,i)) 
            ! Just to be sure
            CurrentSign = 0.0_dp
            ! Search the hash table HashIndex for the determinant defined by
            ! nJ and SpawnedParts(:,i). If it is found, tSuccess will be
            ! returned .true. and PartInd will hold the position of the
            ! determinant in CurrentDets. Else, tSuccess will be returned
            ! .false. (and PartInd shouldn't be accessed).
            ! Also, the hash value, DetHash, is returned by this routine.
            ! tSuccess will determine whether the particle has been found or not.
            call hash_table_lookup(nJ, SpawnedParts(:,i), NIfDBO, HashIndex, &
                                   CurrentDets, PartInd, DetHash, tSuccess)

            tDetermState = .false.

            ! for scaled walkers, truncation is done here
            t_truncate_this_det = t_truncate_spawns .and. .not. &
                 test_flag(SpawnedParts(:,i), flag_multi_spawn) .and. &
                 tEScaleWalkers


!            WRITE(6,*) 'i,SpawnedParts(:,i)',i,SpawnedParts(:,i)

            if (tSuccess) then

                ! Our SpawnedParts determinant is found in CurrentDets.

                call extract_sign(CurrentDets(:,PartInd),CurrentSign)
                call extract_sign(SpawnedParts(:,i),SpawnedSign)

                SignProd = CurrentSign*SpawnedSign

                ! truncate if requested
                if(t_truncate_this_det .and. .not. t_truncate_unocc) then
                   scFVal = scaleFunction(det_diagH(PartInd))
                   do j = 1, lenof_sign
                      call truncateSpawn(iter_data, SpawnedSign, i, j, scFVal, SignProd(j))
                   enddo
                endif

                tDetermState = test_flag(CurrentDets(:,PartInd), flag_deterministic)

                if (sum(abs(CurrentSign)) >= 1.e-12_dp .or. tDetermState) then
                    ! Transfer new sign across.
                    call encode_sign(CurrentDets(:,PartInd), SpawnedSign+CurrentSign)
                    call encode_sign(SpawnedParts(:,i), null_part)                    

                    ! If the sign changed, the adi check has to be redone
                    if(any(real(SignProd,dp) < 0.0_dp)) &
                         call clr_flag(CurrentDets(:,PartInd), flag_adi_checked)

                    ! this det is not prone anymore
                    if(t_prone_walkers) call clr_flag(CurrentDets(:,PartInd), flag_prone)

                    do j = 1, lenof_sign
                        run = part_type_to_run(j)

                        if (is_run_unnocc(CurrentSign,run)) then
                            ! This determinant is actually *unoccupied* for the
                            ! run we're considering. We need to
                            ! decide whether to abort it or not.
                            if (tTruncInitiator) then
                                if (.not. test_flag (SpawnedParts(:,i), get_initiator_flag(j)) .and. &
                                     .not. tDetermState) then
                                    ! Walkers came from outside initiator space.
                                    NoAborted(j) = NoAborted(j) + abs(SpawnedSign(j))
                                    iter_data%naborted(j) = iter_data%naborted(j) + abs(SpawnedSign(j))
                                    call encode_part_sign (CurrentDets(:,PartInd), 0.0_dp, j)
                                end if
                            end if
                        end if
                            

                        if (SignProd(j) < 0) then
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

                    if (.not. tDetermState) then
                        call extract_sign (CurrentDets(:,PartInd), SignTemp)
                        if (IsUnoccDet(SignTemp)) then
                            ! All walkers in this main list have been annihilated
                            ! away. Remove it from the hash index array so that
                            ! no others find it (it is impossible to have another
                            ! spawned walker yet to find this determinant).
                            call remove_hash_table_entry(HashIndex, nJ, PartInd)
                            ! Add to "freeslot" list so it can be filled in.
                            iEndFreeSlot = iEndFreeSlot + 1
                            FreeSlot(iEndFreeSlot) = PartInd
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
                        if (tOldRDMs) call check_fillRDM_DiDj_old(rdms, one_rdms_old, i, CurrentDets(:,PartInd), TempCurrentSign)
                        call check_fillRDM_DiDj(rdm_definitions, two_rdm_spawn, one_rdms, i, &
                                                CurrentDets(:,PartInd), TempCurrentSign)
                    end if 

                end if

            end if
                
            if ( (.not.tSuccess) .or. (tSuccess .and. IsUnoccDet(CurrentSign) .and. (.not. tDetermState)) ) then


                ! Determinant in newly spawned list is not found in CurrentDets.
                ! Usually this would mean the walkers just stay in this list and
                ! get merged later - but in this case we want to check where the
                ! walkers came from - because if the newly spawned walkers are
                ! from a parent outside the active space they should be killed,
                ! as they have been spawned on an unoccupied determinant.
                if (tTruncInitiator) then

                    call extract_sign (SpawnedParts(:,i), SignTemp)

                    ! Are we about to abort this spawn (on any replica) due to
                    ! initiator criterion?
                    do j = 1, lenof_sign
                        abort(j) = test_abort_spawn(SpawnedParts(:, i), j)
                    end do

                    ! If calculating an EN2 correction to initiator error,
                    ! check now if we should add anything.
                    if (tEN2Init) then
                        call add_en2_pert_for_init_calc(i, abort, nJ, SignTemp)
                    end if

                    do j = 1, lenof_sign

                        if (abort(j)) then

                            ! If this option is on, include the walker to be
                            ! cancelled in the trial energy estimate.
                            if (tIncCancelledInitEnergy) call add_trial_energy_contrib(SpawnedParts(:,i), SignTemp(j), j)

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
                        ! truncate to a minimum population given by the scale factor
                    end do

                    if(t_prone_walkers) then
                       if(.not. test_flag(SpawnedParts(:,i), flag_multi_spawn)) &
                            call set_flag(SpawnedParts(:,i), flag_prone)
                    endif

                    if (.not. IsUnoccDet(SignTemp)) then
                       
                       ! if we did not kill the walkers, get the scaling factor
                       call getEScale(nJ, i, diagH, scFVal, ScaledOccupiedThresh)
                       do j = 1, lenof_sign
                          call stochRoundSpawn(iter_data, SignTemp, i, j, scFVal, &
                               ScaledOccupiedThresh, t_truncate_this_det)
                       enddo

                       if(.not. IsUnoccDet(SignTemp)) then
                        ! Walkers have not been aborted and so we should copy the
                        ! determinant straight over to the main list. We do not
                        ! need to recompute the hash, since this should be the
                        ! same one as was generated at the beginning of the loop.
                          if(.not. tEScaleWalkers) diagH = get_diagonal_matel(nJ, SpawnedParts(:,i))
                          call AddNewHashDet(TotWalkersNew, SpawnedParts(0:NIfTot,i), DetHash, nJ, &
                               diagH)
                       end if
                    end if

                else
                    ! Running the full, non-initiator scheme.
                    ! Determinant in newly spawned list is not found in
                    ! CurrentDets. If coeff <1, apply removal criterion.
                    call extract_sign (SpawnedParts(:,i), SignTemp)

                    ! no chance to kill the spawn by initiator criterium
                    ! so get the diagH immediately
                    call getEScale(nJ, i, diagH, scFVal, ScaledOccupiedThresh)

                    ! If using an EN2 perturbation to correct a truncated
                    ! calculation, then this spawn may need to be truncated
                    ! away now. Check this here:
                    tTruncSpawn = .false.
                    if (tEN2Truncated) then
                        tTruncSpawn = .not. CheckAllowedTruncSpawn(0, nJ, SpawnedParts(:,i), 0)
                    end if

                    if (tTruncSpawn) then
                        ! Needs to be truncated away, and a contribution
                        ! added to the EN2 correction.
                        call add_en2_pert_for_trunc_calc(i, nJ, SignTemp)
                    else
                        do j = 1, lenof_sign
                            ! truncate the spawn if required
                            call stochRoundSpawn(iter_data, SignTemp, i, j, scFVal, &
                                 ScaledOccupiedThresh, t_truncate_this_det)
                        end do

                        if (.not. IsUnoccDet(SignTemp)) then
                            ! Walkers have not been aborted and so we should copy the
                            ! determinant straight over to the main list. We do not
                            ! need to recompute the hash, since this should be the
                            ! same one as was generated at the beginning of the loop.
                           if(.not. tEScaleWalkers) diagH = get_diagonal_matel(nJ, SpawnedParts(:,i))
                            call AddNewHashDet(TotWalkersNew, SpawnedParts(:,i), DetHash, nJ, &
                                 diagH)
                        end if
                    end if
                end if

                if (tFillingStochRDMonFly .and. (.not. tNoNewRDMContrib)) then
                    ! We must use the instantaneous value for the off-diagonal contribution.
                    if (tOldRDMs) call check_fillRDM_DiDj_old(rdms, one_rdms_old, i, SpawnedParts(0:NifTot,i), SignTemp)
                    call check_fillRDM_DiDj(rdm_definitions, two_rdm_spawn, one_rdms, i, SpawnedParts(0:NifTot,i), SignTemp)
                end if 
            end if

        end do

        call halt_timer(BinSearch_time)

        ! Update remaining number of holes in list for walkers stats.
        if (iStartFreeSlot > iEndFreeSlot) then
            ! All slots filled
            HolesInList = 0
        else
            HolesInList = iEndFreeSlot - (iStartFreeSlot-1)
        endif

        call halt_timer(AnnMain_time)

    end subroutine AnnihilateSpawnedParts

    subroutine stochRoundSpawn(iter_data, SignTemp, i, j, scFVal, ScaledOccupiedThresh, &
         tTruncate)
      implicit none
      type(fcimc_iter_data), intent(inout) :: iter_data
      real(dp), intent(inout) :: SignTemp(lenof_sign)
      integer, intent(in) :: i, j
      real(dp), intent(in) :: scFVal, ScaledOccupiedThresh
      logical, intent(in) :: tTruncate

      real(dp) :: pRemove, r
      integer :: run

      run = part_type_to_run(j)

      if ((abs(SignTemp(j)) > 1.e-12_dp) .and. (abs(SignTemp(j)) < ScaledOccupiedThresh)) then
         ! We remove this walker with probability OccupiedThresh - Sign/ScaleFactor
         pRemove=1.0_dp-abs(SignTemp(j))/(ScaledOccupiedThresh)
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
         else
            !Round up
            NoBorn(run) = NoBorn(run) + OccupiedThresh*scFVal - abs(SignTemp(j))
            iter_data%nborn(j) = iter_data%nborn(j) &
                 + scaledOccupiedThresh - abs(SignTemp(j))
            SignTemp(j) = sign(scaledOccupiedThresh, SignTemp(j))
            call encode_part_sign (SpawnedParts(:,i), SignTemp(j), j)
         end if
      else if(abs(SignTemp(j)) > eps) then
         ! truncate down to a minimum number of spawns to
         ! prevent blooms if requested
         if(tTruncate) then
            call truncateSpawn(iter_data,SignTemp,i,j,scFVal,1.0_dp)
            call encode_part_sign(SpawnedParts(:,i), SignTemp(j), j)
         endif
      end if

    end subroutine stochRoundSpawn

    subroutine truncateSpawn(iter_data, SignTemp, i, j, scFVal, SignProd)
      implicit none
      type(fcimc_iter_data), intent(inout) :: iter_data
      real(dp), intent(inout) :: SignTemp(lenof_sign)
      integer, intent(in) :: i, j ! i: index of ilut in SpawnedParts, j: part index
      real(dp), intent(in) :: scFVal, SignProd

      real(dp) :: maxSpawns

      ! we allow n_truncate_spawns unit walkers to be created per spawn event
      maxSpawns = n_truncate_spawns * scFVal * get_num_spawns(SpawnedParts(:,i))

      ! truncate the new walkers to a maximum value
      if(abs(SignTemp(j)) > maxSpawns) then
         iter_data%nremoved(j) = iter_data%nremoved(j) + &
              abs(SignTemp(j)) - sign(maxSpawns,SignProd)
         ! log the truncated weight
         truncatedWeight = truncatedWeight + abs(SignTemp(j)) - maxSpawns
         ! reduce the sign to maxSpawns
         SignTemp(j) = sign(maxSpawns, SignTemp(j))
      endif

    end subroutine truncateSpawn

    subroutine getEScale(nJ, i, diagH, scFVal, ScaledOccupiedThresh)
      implicit none
      integer, intent(in) :: nJ(nel), i
      real(dp), intent(out) :: diagH, scFVal, ScaledOccupiedThresh

      if(tEScaleWalkers) then
         ! the diagonal element of H is needed anyway ONLY for scaled walkers
         ! if we dont scale walkers, we might not need it, if the round kills the 
         ! walkers
         diagH = get_diagonal_matel(nJ, SpawnedParts(:,i))
         ! evaluate the scaling function
         scFVal = scaleFunction(diagH - Hii)
      else
         scFVal = 1.0_dp
      endif
      ScaledOccupiedThresh = scFVal * OccupiedThresh
    end subroutine getEScale

    pure function test_abort_spawn(ilut_spwn, part_type) result(abort)

        ! Should this spawn be aborted (according to the initiator
        ! criterion).
        !
        ! This function accepts an initiator spawn with probability
        ! (1.0 - alpha(ilut_spwn, part_type)), which can be artibrarily
        ! complicated.

        integer(n_int), intent(in) :: ilut_spwn(0:nIfBCast)
        integer, intent(in) :: part_type
        logical :: abort
        integer :: maxExLvl, nopen

        ! If a particle comes from a site marked as an initiator, then it can
        ! live

        abort = .not. test_flag(ilut_spwn, get_initiator_flag(part_type))
        
        ! optionally keep spawns up to a given seniority level + excitaion level
        if(abort .and. tSpawnSeniorityBased) then
           ! get the seniority level
           nopen = count_open_orbs(ilut_spwn) 
           if(nopen < numMaxExLvlsSet) then
              maxExLvl = 0
              ! get the corresponding max excitation level
              do while(maxExLvl == 0)
                 maxExLvl = maxKeepExLvl(nopen+1)
                 nopen = nopen + 1
                 if(nopen >= numMaxExLvlsSet) exit
              end do
              ! if we are below this level, keep the spawn anyway
              if(FindBitExcitLevel(ilutHF, ilut_spwn) <= maxExLvl) abort = .false.
           end if
        end if

    end function test_abort_spawn

    subroutine add_en2_pert_for_init_calc(ispawn, abort, nJ, SpawnedSign)

        ! Add a contribution to the second-order Epstein-Nesbet correction to
        ! initiator error.

        ! This adds a correction for the determinant at position ispawn in
        ! the spawning array, and nJ is an array holding the occupied
        ! electrons on this determinant. SpawnedSign is the array of
        ! amplitudes spawned onto this determinant, and abort is array
        ! specifying whether each of the spawnings are about to be aborted
        ! due to the initiator criterion.

        integer, intent(in) :: ispawn
        logical, intent(in) :: abort(lenof_sign)
        integer, intent(in) :: nJ(nel)
        real(dp), intent(in):: SpawnedSign(lenof_sign)

        integer :: istate
        real(dp) :: contrib_sign(en_pert_main%sign_length)
        logical :: pert_contrib(en_pert_main%sign_length)
        real(dp) :: trial_contrib(lenof_sign)
        HElement_t(dp) :: h_diag

        ! RDM-energy-based estimate:
        ! Only add a contribution if we've started accumulating this estimate.
        if (tEN2Started) then
            pert_contrib = .false.

            do istate = 1, en_pert_main%sign_length
                ! Was a non-zero contribution aborted on *both* replicas for
                ! a given state?
                if (abort(2*istate-1) .and. abort(2*istate) .and. &
                      abs(SpawnedSign(2*istate-1)) > 1.e-12_dp .and. &
                      abs(SpawnedSign(2*istate)) > 1.e-12_dp) &
                    pert_contrib(istate) = .true.
            end do

            if (any(pert_contrib)) then
                contrib_sign = 0.0_dp
                do istate = 1, en_pert_main%sign_length
                    if (pert_contrib(istate)) then
                        contrib_sign(istate) = SpawnedSign(2*istate-1)*SpawnedSign(2*istate) / (tau**2)
                    end if
                end do

                call add_to_en_pert_t(en_pert_main, nJ, SpawnedParts(:,ispawn), contrib_sign)
            end if
        end if

    end subroutine add_en2_pert_for_init_calc

    subroutine add_en2_pert_for_trunc_calc(ispawn, nJ, SpawnedSign)

        ! Add a contribution to the second-order Epstein-Nesbet correction to
        ! error due to truncation of the Hilbert space being sampled.

        ! This adds a correction for the determinant at position ispawn in
        ! the spawning array, and nJ is an array holding the occupied
        ! electrons on this determinant. SpawnedSign is the array of
        ! amplitudes spawned onto this determinant.

        integer, intent(in) :: ispawn
        integer, intent(in) :: nJ(nel)
        real(dp), intent(inout):: SpawnedSign(lenof_sign)

        integer :: j, istate
        real(dp) :: contrib_sign(en_pert_main%sign_length)
        logical :: pert_contrib(en_pert_main%sign_length)

        ! Only add a contribution if we've started accumulating this estimate.
        if (tEN2Started) then
            pert_contrib = .false.

            do istate = 1, en_pert_main%sign_length
                if (abs(SpawnedSign(2*istate-1)) > 1.e-12_dp .and. &
                      abs(SpawnedSign(2*istate)) > 1.e-12_dp) then
                    pert_contrib(istate) = .true.
                end if
            end do

            if (any(pert_contrib)) then
                contrib_sign = 0.0_dp
                do istate = 1, en_pert_main%sign_length
                    if (pert_contrib(istate)) then
                        contrib_sign(istate) = SpawnedSign(2*istate-1)*SpawnedSign(2*istate) / (tau**2)
                    end if
                end do

                call add_to_en_pert_t(en_pert_main, nJ, SpawnedParts(:,ispawn), contrib_sign)
            end if
        end if

        ! Remove the spawning
        do j = 1, lenof_sign
            SpawnedSign(j) = 0.0_dp
            call encode_part_sign (SpawnedParts(:,ispawn), SpawnedSign(j), j)
        end do

    end subroutine add_en2_pert_for_trunc_calc

end module AnnihilationMod
