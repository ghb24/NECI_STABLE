#include "macros.h"

! Some general procedures, created for the semi-stochastic (and trial wavefunction) code.
! Some are just used for initialisation and others in the main FCIQMC loop.

module semi_stoch_procs

    use bit_rep_data, only: flag_deterministic, nIfDBO, NIfD, NIfTot, test_flag
    use bit_reps, only: decode_bit_det, get_initiator_flag_by_run
    use CalcData
    use constants
    use FciMCData, only: determ_sizes, determ_displs, determ_space_size, &
                         SpawnedParts, TotWalkers, CurrentDets, core_space, &
                         MaxSpawned,indices_of_determ_states, ilutRef, determ_last
    use Parallel_neci, only: iProcIndex, nProcessors, MPIArg
    use sparse_arrays, only: sparse_core_ham, approx_ham
    use procedure_pointers, only: shiftFactorFunction
    use SystemData, only: nel
    use timing_neci
    use adi_data, only: tSignedRepAv
    use global_det_data, only: set_tot_acc_spawns, set_apvals, tAccumEmptyDet
    use LoggingData, only: tAccumPopsActive

    implicit none

contains

    subroutine determ_projection()

        ! This subroutine gathers together partial_determ_vecs from each processor so
        ! that the full vector for the whole deterministic space is stored on each processor.
        ! It then performs the deterministic multiplication of the projector on this full vector.

        use FciMCData, only: partial_determ_vecs, full_determ_vecs, SemiStoch_Comms_Time
        use FciMCData, only: SemiStoch_Multiply_Time
        use Parallel_neci, only: MPIBarrier, MPIAllGatherV
        use DetBitOps, only: DetBitEQ

        integer :: i, j, ierr, run, part_type
        real(dp) :: scaledDiagSft(inum_runs)
        call MPIBarrier(ierr)

        call set_timer(SemiStoch_Comms_Time)

        call MPIAllGatherV(partial_determ_vecs, full_determ_vecs, &
                            determ_sizes, determ_displs)

        call halt_timer(SemiStoch_Comms_Time)

        call set_timer(SemiStoch_Multiply_Time)

        if (determ_sizes(iProcIndex) >= 1) then

            ! For the moment, we're only adding in these contributions when we need the energy
            ! This will need refinement if we want to continue with the option of inst vs true full RDMs
            ! (as in another CMO branch).

            ! Perform the multiplication.

            partial_determ_vecs = 0.0_dp

#ifdef __CMPLX
            do i = 1, determ_sizes(iProcIndex)
                do j = 1, sparse_core_ham(i)%num_elements
                    do run = 1, inum_runs
                        partial_determ_vecs(min_part_type(run),i) = partial_determ_vecs(min_part_type(run),i) - &
                            Real(sparse_core_ham(i)%elements(j))*full_determ_vecs(min_part_type(run),sparse_core_ham(i)%positions(j)) +&
                            Aimag(sparse_core_ham(i)%elements(j))*full_determ_vecs(max_part_type(run),sparse_core_ham(i)%positions(j))
                        partial_determ_vecs(max_part_type(run),i) = partial_determ_vecs(max_part_type(run),i) - &
                            Aimag(sparse_core_ham(i)%elements(j))*full_determ_vecs(min_part_type(run),sparse_core_ham(i)%positions(j)) -&
                            Real(sparse_core_ham(i)%elements(j))*full_determ_vecs(max_part_type(run),sparse_core_ham(i)%positions(j))
                    end do
                end do
            end do
#else
            do i = 1, determ_sizes(iProcIndex)
                do j = 1, sparse_core_ham(i)%num_elements
                    partial_determ_vecs(:,i) = partial_determ_vecs(:,i) - &
                        sparse_core_ham(i)%elements(j)*full_determ_vecs(:,sparse_core_ham(i)%positions(j))
                end do
            end do
#endif

            ! Now add shift*full_determ_vecs to account for the shift, not stored in
            ! sparse_core_ham.
#ifdef __CMPLX
            do i = 1, determ_sizes(iProcIndex)
               ! Only scale the shift for the corespace when the option is set
               if(tCoreAdaptiveShift .and. tAdaptiveShift) then
                  do part_type = 1, inum_runs
                     ! scale the shift using the abs of this run's complex coefficient
                     scaledDiagSft(part_type) = &
                     shiftFactorFunction(&
                          indices_of_determ_states(i), part_type,&
                          sqrt(full_determ_vecs(min_part_type(part_type),&
                          i+determ_displs(iProcIndex))**2 + &
                          full_determ_vecs(max_part_type(part_type),&
                          i+determ_displs(iProcIndex))**2)) *  DiagSft(part_type)
                  end do
               else
                  scaledDiagSft = DiagSft
               endif

               do part_type  = 1, lenof_sign
                  partial_determ_vecs(part_type,i) = partial_determ_vecs(part_type,i) + &
                       scaledDiagSft(part_type_to_run(part_type)) * full_determ_vecs(part_type,i+determ_displs(iProcIndex))
               enddo
            end do
#else
            do i = 1, determ_sizes(iProcIndex)
               ! Only scale the shift for the corespace when the option is set
               if(tCoreAdaptiveShift .and. tAdaptiveShift) then
                  ! get the re-scaled shift accounting for undersampling error
                  do  part_type = 1, inum_runs
                     scaledDiagSft(part_type) = DiagSft(part_type) * shiftFactorFunction(&
                          indices_of_determ_states(i), part_type,&
                          abs(full_determ_vecs(part_type,i+determ_displs(iProcIndex))))
                  end do
               else
                  scaledDiagSft = DiagSft
               endif
               partial_determ_vecs(:,i) = partial_determ_vecs(:,i) + &
                    scaledDiagSft * full_determ_vecs(:,i+determ_displs(iProcIndex))
            end do
#endif

            ! Now multiply the vector by tau to get the final projected vector.
            partial_determ_vecs = partial_determ_vecs * tau

            do i = 1, determ_sizes(iProcIndex)
                do part_type  = 1, lenof_sign
                    run = part_type_to_run(part_type)
                    if(tSkipRef(run) .and. DetBitEQ(CurrentDets(:,indices_of_determ_states(i)),iLutRef(:,run),nIfD)) then
                        partial_determ_vecs(part_type, i) = 0.0_dp
                    end if
                end do
            end do
        end if

        call halt_timer(SemiStoch_Multiply_Time)

    end subroutine determ_projection

    subroutine determ_projection_kp_hamil(partial_vecs, full_vecs, determ_sizes, determ_disps)

        use FciMCData, only: Hii, SemiStoch_Comms_Time, SemiStoch_Multiply_Time
        use Parallel_neci, only: MPIBarrier, MPIAllGatherV

        real(dp), allocatable, intent(inout) :: partial_vecs(:,:)
        real(dp), allocatable, intent(inout) :: full_vecs(:,:)
        integer(MPIArg), allocatable, intent(in) :: determ_sizes(:), determ_disps(:)

        integer :: i, j, ierr

        call MPIBarrier(ierr)

        call set_timer(SemiStoch_Comms_Time)

        call MPIAllGatherV(partial_vecs, full_vecs, determ_sizes, determ_disps)

        call halt_timer(SemiStoch_Comms_Time)

        call MPIBarrier(ierr)

        call set_timer(SemiStoch_Multiply_Time)

        if (determ_sizes(iProcIndex) >= 1) then
            ! Start with this because sparse_core_hamil has Hii taken off, but actually we
            ! don't want the projected Hamiltonian to be relative to the HF determinant.
            partial_vecs = Hii*full_vecs(:, determ_disps(iProcIndex)+1:&
                                            determ_disps(iProcIndex)+determ_sizes(iProcIndex))

            do i = 1, determ_sizes(iProcIndex)
                do j = 1, sparse_core_ham(i)%num_elements
                    partial_vecs(:,i) = partial_vecs(:,i) + &
                        sparse_core_ham(i)%elements(j)*full_vecs(:,sparse_core_ham(i)%positions(j))
                end do
            end do
        end if

        call halt_timer(SemiStoch_Multiply_Time)

    end subroutine determ_projection_kp_hamil

    subroutine determ_projection_no_death()

        ! This subroutine gathers together partial_determ_vecs from each processor so
        ! that the full vector for the whole deterministic space is stored on each processor.
        ! It then performs the deterministic multiplication of the projector on this full vector.

        use DetBitOps, only: DetBitEQ
        use FciMCData, only: partial_determ_vecs, full_determ_vecs, SemiStoch_Comms_Time
        use FciMCData, only: SemiStoch_Multiply_Time, core_ham_diag
        use Parallel_neci, only: MPIBarrier, MPIAllGatherV

        integer :: i, j, ierr, run, part_type

        call MPIBarrier(ierr)

        call set_timer(SemiStoch_Comms_Time)

        call MPIAllGatherV(partial_determ_vecs, full_determ_vecs, &
                            determ_sizes, determ_displs)

        call halt_timer(SemiStoch_Comms_Time)

        call set_timer(SemiStoch_Multiply_Time)

        if (determ_sizes(iProcIndex) >= 1) then

            ! Perform the multiplication.

            partial_determ_vecs = 0.0_dp

#ifdef __CMPLX
            do i = 1, determ_sizes(iProcIndex)
                do j = 1, sparse_core_ham(i)%num_elements
                    do run = 1, inum_runs
                        partial_determ_vecs(min_part_type(run),i) = partial_determ_vecs(min_part_type(run),i) - &
                            Real(sparse_core_ham(i)%elements(j))*full_determ_vecs(min_part_type(run),sparse_core_ham(i)%positions(j)) +&
                            Aimag(sparse_core_ham(i)%elements(j))*full_determ_vecs(max_part_type(run),sparse_core_ham(i)%positions(j))
                        partial_determ_vecs(max_part_type(run),i) = partial_determ_vecs(max_part_type(run),i) - &
                            Aimag(sparse_core_ham(i)%elements(j))*full_determ_vecs(min_part_type(run),sparse_core_ham(i)%positions(j)) -&
                            Real(sparse_core_ham(i)%elements(j))*full_determ_vecs(max_part_type(run),sparse_core_ham(i)%positions(j))
                    end do
                end do
            end do
#else
            do i = 1, determ_sizes(iProcIndex)
                do j = 1, sparse_core_ham(i)%num_elements
                    partial_determ_vecs(:,i) = partial_determ_vecs(:,i) - &
                        sparse_core_ham(i)%elements(j)*full_determ_vecs(:,sparse_core_ham(i)%positions(j))
                end do
            end do
#endif

            ! Remove contribution from the diagonal elements, since the
            ! propagator has no diagonal.
            do i = 1, determ_sizes(iProcIndex)
                partial_determ_vecs(:,i) = partial_determ_vecs(:,i) + &
                   core_ham_diag(i) * full_determ_vecs(:,i+determ_displs(iProcIndex))
            end do

            ! Now multiply the vector by tau to get the final projected vector.
            partial_determ_vecs = partial_determ_vecs * tau

            do i = 1, determ_sizes(iProcIndex)
                do part_type  = 1, lenof_sign
                    run = part_type_to_run(part_type)
                    if(tSkipRef(run) .and. DetBitEQ(CurrentDets(:,indices_of_determ_states(i)),iLutRef(:,run),nIfD)) then
                        partial_determ_vecs(part_type, i) = 0.0_dp
                    end if
                end do
            end do

        end if

        call halt_timer(SemiStoch_Multiply_Time)

    end subroutine determ_projection_no_death

    subroutine determ_proj_approx()

        ! This subroutine gathers together partial_determ_vecs from each processor so
        ! that the full vector for the whole deterministic space is stored on each processor.
        ! It then performs the deterministic multiplication of the projector on this full vector.

        use FciMCData, only: partial_determ_vecs, full_determ_vecs, SemiStoch_Comms_Time
        use FciMCData, only: SemiStoch_Multiply_Time
        use Parallel_neci, only: MPIBarrier, MPIAllGatherV
        use DetBitOps, only: DetBitEQ

        integer :: i, j, ierr, run, part_type

        call MPIBarrier(ierr)

        call set_timer(SemiStoch_Comms_Time)

        call MPIAllGatherV(partial_determ_vecs, full_determ_vecs, &
                            determ_sizes, determ_displs)

        call halt_timer(SemiStoch_Comms_Time)

        call set_timer(SemiStoch_Multiply_Time)

        if (determ_sizes(iProcIndex) >= 1) then

            ! For the moment, we're only adding in these contributions when we need the energy
            ! This will need refinement if we want to continue with the option of inst vs true full RDMs
            ! (as in another CMO branch).

            ! Perform the multiplication.

            partial_determ_vecs = 0.0_dp

#ifdef __CMPLX
            do i = 1, determ_sizes(iProcIndex)
                do j = 1, approx_ham(i)%num_elements
                    do run = 1, inum_runs
                        partial_determ_vecs(min_part_type(run),i) = partial_determ_vecs(min_part_type(run),i) - &
                            Real(approx_ham(i)%elements(j))*full_determ_vecs(min_part_type(run),approx_ham(i)%positions(j)) +&
                            Aimag(approx_ham(i)%elements(j))*full_determ_vecs(max_part_type(run),approx_ham(i)%positions(j))
                        partial_determ_vecs(max_part_type(run),i) = partial_determ_vecs(max_part_type(run),i) - &
                            Aimag(approx_ham(i)%elements(j))*full_determ_vecs(min_part_type(run),approx_ham(i)%positions(j)) -&
                            Real(approx_ham(i)%elements(j))*full_determ_vecs(max_part_type(run),approx_ham(i)%positions(j))
                    end do
                end do
            end do
#else
            do i = 1, determ_sizes(iProcIndex)
                do j = 1, approx_ham(i)%num_elements
                    partial_determ_vecs(:,i) = partial_determ_vecs(:,i) - &
                        approx_ham(i)%elements(j)*full_determ_vecs(:,approx_ham(i)%positions(j))
                end do
            end do
#endif

            ! Now add shift*full_determ_vecs to account for the shift, not stored in
            ! approx_ham.
#ifdef __CMPLX
            do i = 1, determ_sizes(iProcIndex)
                do part_type  = 1, lenof_sign
                    partial_determ_vecs(part_type,i) = partial_determ_vecs(part_type,i) + &
                       DiagSft(part_type_to_run(part_type)) * full_determ_vecs(part_type,i+determ_displs(iProcIndex))
                enddo
            end do
#else
            do i = 1, determ_sizes(iProcIndex)
                partial_determ_vecs(:,i) = partial_determ_vecs(:,i) + &
                   DiagSft * full_determ_vecs(:,i+determ_displs(iProcIndex))
            end do
#endif

            ! Now multiply the vector by tau to get the final projected vector.
            partial_determ_vecs = partial_determ_vecs * tau

            do i = 1, determ_sizes(iProcIndex)
                do part_type  = 1, lenof_sign
                    run = part_type_to_run(part_type)
                    if(tSkipRef(run) .and. DetBitEQ(CurrentDets(:,indices_of_determ_states(i)),iLutRef(:,run),nIfD)) then
                        partial_determ_vecs(part_type, i) = 0.0_dp
                    end if
                end do
            end do
        end if

        call halt_timer(SemiStoch_Multiply_Time)

    end subroutine determ_proj_approx

    subroutine average_determ_vector()

        use FciMCData, only: Iter, IterRDMStart
        use FciMCData, only: full_determ_vecs, full_determ_vecs_av, PreviousCycles
        use LoggingData, only: RDMEnergyIter

        real(dp) :: iter_curr, iter_start_av

        ! If this condition is met then RDM energies were added in on the
        ! previous iteration. We now want to start a new averaging block so
        ! that the same contributions aren't added in again later.
        if (mod(Iter+PreviousCycles-IterRDMStart, RDMEnergyIter) == 0) then
            full_determ_vecs_av = 0.0_dp
            write(6,*) "Reset fdv av at iteration ", iter
        end if

        ! The current iteration, converted to a double precision real.
        iter_curr = real(Iter+PreviousCycles, dp)
        ! The iteration that this averaging block started on.
        iter_start_av = real(RDMEnergyIter*((Iter+PreviousCycles - IterRDMStart)/RDMEnergyIter) + IterRDMStart, dp)

        ! Add in the current deterministic vector to the running average.
        full_determ_vecs_av = (((iter_curr - iter_start_av)*full_determ_vecs_av) + full_determ_vecs)/&
                                 (iter_curr - iter_start_av + 1.0_dp)

    end subroutine average_determ_vector

    function is_core_state(ilut, nI) result (core_state)

        use FciMCData, only: determ_space_size_int
        use hash, only: FindWalkerHash
        use sparse_arrays, only: core_ht

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: nI(:)
        integer :: i, hash_val
        logical :: core_state

        core_state = .false.

        hash_val = FindWalkerHash(nI, determ_space_size_int)

        do i = 1, core_ht(hash_val)%nclash
            if (all(ilut(0:NIfDBO) == core_space(0:NIfDBO,core_ht(hash_val)%ind(i)) )) then
                core_state = .true.
                return
            end if
        end do

    end function is_core_state

    function core_space_pos(ilut, nI) result (pos)

        use FciMCData, only: ll_node, determ_space_size_int
        use hash, only: FindWalkerHash
        use sparse_arrays, only: core_ht

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: nI(:)
        integer :: i, hash_val
        character(len=*), parameter :: t_r = "core_space_pos"

        integer :: pos

        pos = 0

        hash_val = FindWalkerHash(nI, determ_space_size_int)

        do i = 1, core_ht(hash_val)%nclash
            if (all(ilut(0:NIfDBO) == core_space(0:NIfDBO,core_ht(hash_val)%ind(i)) )) then
                pos = core_ht(hash_val)%ind(i)
                return
            end if
        end do

        if (pos == 0) then
            call stop_all(t_r, "State not found in core hash table.")
        end if

    end function core_space_pos

    function check_determ_flag(ilut) result (core_state)

        ! The reason for using this instead of just using test_flag is that test_flag
        ! crashes if flags are not being used. Calling this function therefore makes
        ! things neater!

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        logical :: core_state

        if (tSemiStochastic) then
            core_state = test_flag(ilut, flag_deterministic)
        else
            core_state = .false.
        end if

    end function check_determ_flag

    subroutine recalc_core_hamil_diag(old_Hii, new_Hii)

        use FciMCData, only: core_ham_diag

        real(dp) :: old_Hii, new_Hii
        real(dp) :: Hii_shift
        integer :: i, j

        ! Only attempt this if we have already performed the semi-stochastic
        ! initialisation, in which case determ_sizes will have been allocated.
        if (allocated(determ_sizes)) then
            write(6,'(a56)') "Recalculating diagonal elements of the core Hamiltonian."

            Hii_shift = old_Hii - new_Hii

            do i = 1, determ_sizes(iProcIndex)
                do j = 1, sparse_core_ham(i)%num_elements
                    if (sparse_core_ham(i)%positions(j) == i + determ_displs(iProcIndex)) then
                        sparse_core_ham(i)%elements(j) = sparse_core_ham(i)%elements(j) + Hii_shift
                    end if
                end do
            end do

            core_ham_diag = core_ham_diag + Hii_shift
        end if

    end subroutine recalc_core_hamil_diag

    subroutine generate_core_connections()

        use DetBitOps, only: FindBitExcitLevel
        use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc
        use Parallel_neci, only: MPIAllGatherV
        use sparse_arrays, only: core_connections

        integer :: i, j, ic, counter, ierr
        integer :: Ex(2,nel)
        logical :: tSign
        integer(n_int), allocatable, dimension(:,:) :: temp_store
        integer(TagIntType) :: TempStoreTag
        character(len=*), parameter :: t_r = "calculate_det_hamiltonian_sparse"

        allocate(core_connections(determ_sizes(iProcIndex)))

        allocate(temp_store(0:NIfTot, determ_space_size), stat=ierr)
        call LogMemAlloc('temp_store', maxval(determ_sizes)*(NIfTot+1), 8, t_r, &
                         TempStoreTag, ierr)

        ! Stick together the deterministic states from all processors, on all processors.
        call MPIAllGatherV(SpawnedParts(0:NIfTot,1:determ_sizes(iProcIndex)), temp_store, &
                       determ_sizes, determ_displs)

        ! Over all core states on this processor.
        do i = 1, determ_sizes(iProcIndex)

            ! The number of non-zero elements in this array will be almost the same as in
            ! the core Hamiltonian array, except the diagonal element is not considered,
            ! so there will actually be one less.
            allocate(core_connections(i)%elements(sparse_core_ham(i)%num_elements-1))
            allocate(core_connections(i)%positions(sparse_core_ham(i)%num_elements-1))

            ! The total number of non-zero elements in row i.
            core_connections(i)%num_elements = sparse_core_ham(i)%num_elements-1

            counter = 0
            do j = 1, sparse_core_ham(i)%num_elements
                ! If not the diagonal element.
                if (sparse_core_ham(i)%positions(j) /= i + determ_displs(iProcIndex)) then
                    Ex = 0
                    Ex(1,1) = nel
                    counter = counter + 1
                    ! The positions of the non-zero and non-diagonal elements in this row i.
                    core_connections(i)%positions(counter) = sparse_core_ham(i)%positions(j)

                    ic = FindBitExcitLevel(SpawnedParts(:,i), temp_store(:, sparse_core_ham(i)%positions(j)))
                    call GetBitExcitation(SpawnedParts(0:NIfD,i), temp_store(0:NIfD, &
                                          sparse_core_ham(i)%positions(j)),Ex,tSign)
                    if (tSign) then
                        ! Odd number of permutations. Minus the excitation level.
                        core_connections(i)%elements(counter) = -ic
                    else
                        ! Even number of permutations. The excitation level.
                        core_connections(i)%elements(counter) = ic
                    end if
                end if
            end do

        end do

        deallocate(temp_store, stat=ierr)
        call LogMemDealloc(t_r, TempStoreTag, ierr)

    end subroutine generate_core_connections

    subroutine store_whole_core_space()

        use FciMCData, only: CoreSpaceTag
        use MemoryManager, only: LogMemAlloc
        use Parallel_neci, only: MPIAllGatherV

        integer :: ierr
        character(len=*), parameter :: t_r = "store_whole_core_space"

        allocate(core_space(0:NIfTot, determ_space_size), stat=ierr)
        call LogMemAlloc('core_space', maxval(determ_sizes)*(NIfTot+1), 8, t_r, &
                         CoreSpaceTag, ierr)
        core_space = 0_n_int

        ! Give explicit limits for SpawnedParts slice, as NIfTot is not nesc.
        ! equal to NIfBCast. (It may be longer)
        call MPIAllGatherV(SpawnedParts(0:NIfTot, 1:determ_sizes(iProcIndex)),&
                           core_space, determ_sizes, determ_displs)

    end subroutine store_whole_core_space

    subroutine initialise_core_hash_table(ilut_list, space_size, hash_table)

        use bit_reps, only: decode_bit_det
        use hash, only: FindWalkerHash
        use sparse_arrays, only: core_hashtable
        use SystemData, only: nel

        integer(n_int), intent(in) :: ilut_list(0:,:)
        integer, intent(in) :: space_size
        type(core_hashtable), allocatable, intent(out) :: hash_table(:)

        integer :: nI(nel)
        integer :: i, ierr, hash_val

        allocate(hash_table(space_size), stat=ierr)

        do i = 1, space_size
            hash_table(i)%nclash = 0
        end do

        ! Count the number of states with each hash value.
        do i = 1, space_size
            call decode_bit_det(nI, ilut_list(:,i))
            hash_val = FindWalkerHash(nI, int(space_size,sizeof_int))
            hash_table(hash_val)%nclash = hash_table(hash_val)%nclash + 1
        end do

        do i = 1, space_size
            allocate(hash_table(i)%ind(hash_table(i)%nclash), stat=ierr)
            hash_table(i)%ind = 0
            ! Reset this for now.
            hash_table(i)%nclash = 0
        end do

        ! Now fill in the indices of the states in the space.
        do i = 1, space_size
            call decode_bit_det(nI, ilut_list(:,i))
            hash_val = FindWalkerHash(nI, int(space_size, sizeof_int))
            hash_table(hash_val)%nclash = hash_table(hash_val)%nclash + 1
            hash_table(hash_val)%ind(hash_table(hash_val)%nclash) = i
        end do

    end subroutine initialise_core_hash_table

    subroutine remove_high_energy_orbs(ilut_list, num_states, target_num_states, tParallel)

        use Parallel_neci, only: MPISumAll
        use sort_mod, only: sort
        use SystemData, only: nBasis, BRR, Arr

        integer, intent(inout) :: num_states
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot, 1:num_states)
        integer, intent(in) :: target_num_states
        logical, intent(in) :: tParallel
        integer, allocatable, dimension(:) :: orbitals_rmvd
        integer :: i, j, orb, num_orbs_rmvd
        integer :: bit, elem
        integer :: states_rmvd_this_proc, counter
        integer(MPIArg) :: tot_num_states, states_rmvd_all_procs
        logical :: occupied

        states_rmvd_this_proc = 0
        num_orbs_rmvd = 0

        if (tParallel) then
            call MPISumAll(int(num_states, MPIArg), tot_num_states)
        else
            tot_num_states = int(num_states,MPIArg)
        end if

        if (tot_num_states <= target_num_states) return

        write(6,'(a32)') "Removing high energy orbitals..."

        ! Loop through all orbitals, from highest energy to lowest.
        do i = nBasis, 1, -1

            num_orbs_rmvd = nBasis-i+1

            orb = BRR(i)

            bit = mod(orb-1, bits_n_int)
            elem = (orb-1-bit)/bits_n_int

            ! Loop through all states and remove those states with orbital orb
            ! occupied.
            do j = 1, num_states
                occupied = btest(ilut_list(elem, j), bit)
                if (occupied) then
                    ilut_list(:, j) = 0_n_int
                    states_rmvd_this_proc = states_rmvd_this_proc + 1
                end if
            end do

            if (tParallel) then
                call MPISumAll(int(states_rmvd_this_proc, MPIArg), states_rmvd_all_procs)
            else
                states_rmvd_all_procs = int(states_rmvd_this_proc,MPIArg)
            end if

            ! If there are degenerate orbitals, then cycle to remove the
            ! degenerate orbitals too, before giving the program chance to quit.
            if (i > 1) then
                ! If the orbitals energies are the same:
                if (abs(Arr(i, 1) - Arr((i-1), 1)) < 1.0e-12_dp ) cycle
            end if

            if (tot_num_states-states_rmvd_all_procs <= target_num_states) exit

        end do

        ! Loop through the list and shuffle states down to fill in the gaps
        ! created above.
        counter = 0
        do i = 1, num_states
            ! If the state wasn't set to 0:
            if (.not. all(ilut_list(:,i) == 0_n_int)) then
                counter = counter + 1
                ilut_list(:, counter) = ilut_list(:, i)
            end if
        end do
        if (counter /= num_states-states_rmvd_this_proc) &
            call stop_all("remove_high_energy_orbs", &
                          "Wrong number of states found.")
        num_states = counter

        ! Finally, output information:
        write(6,'(a36)',advance='no') "The following orbitals were removed:"
        allocate(orbitals_rmvd(num_orbs_rmvd))
        orbitals_rmvd = BRR(nBasis-num_orbs_rmvd+1:nBasis)
        call sort(orbitals_rmvd)
        do i = 1, num_orbs_rmvd
            write(6,'(i5)',advance='no') orbitals_rmvd(i)
            if (i == num_orbs_rmvd) write(6,'()',advance='yes')
        end do
        deallocate(orbitals_rmvd)
        write(6,'(i6,1x,a29)') states_rmvd_all_procs, "states were removed in total."
        write(6,'(i6,1x,a17)') tot_num_states-states_rmvd_all_procs, "states were kept."
        call neci_flush(6)

    end subroutine remove_high_energy_orbs

    subroutine sort_states_by_energy(ilut_list, num_states, tSortDoubles)

        ! Note: If requested, keep all doubles at the top, then sort by energy.

        use bit_reps, only: decode_bit_det
        use DetBitOps, only: FindBitExcitLevel
        use Determinants, only: get_helement
        use FciMCData, only: ilutHF
        use hphf_integrals, only: hphf_diag_helement
        use sort_mod, only: sort
        use SystemData, only: tHPHF

        integer, intent(in) :: num_states
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot, 1:num_states)
        logical, intent(in) :: tSortDoubles
        integer(n_int) :: temp_ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: i, excit_level, num_sing_doub, block_size
        real(dp) :: energies(num_states)

        do i = 1, num_states
            call decode_bit_det(nI, ilut_list(:,i))
            if (tHPHF) then
                energies(i) = hphf_diag_helement(nI, ilut_list(:,i))
            else
                energies(i) = get_helement(nI, nI, 0)
            end if
        end do

        ! Sort the states in order of energy, from smallest to largest.
        call sort(energies, ilut_list(0:NIfTot, 1:num_states))

        ! If requested, sort singles and doubles to the top, keeping the rest of
        ! the ordering the same.
        if (tSortDoubles) then
            num_sing_doub = 0
            block_size = 0
            do i = 1, num_states
                excit_level = FindBitExcitLevel(ilut_list(:,i), ilutHF)
                if (excit_level <= 2) then
                    num_sing_doub = num_sing_doub + 1
                    temp_ilut = ilut_list(:,i)

                    ! Move block of non-singles or doubles down one in ilut_list.
                    ilut_list(:, num_sing_doub+1:num_sing_doub+block_size) = &
                        ilut_list(:, num_sing_doub:num_sing_doub+block_size-1)

                    ! Then move the single or double just found into the space
                    ! created above this block.
                    ilut_list(:, num_sing_doub) = temp_ilut
                else
                    block_size = block_size + 1
                end if
            end do
        end if

    end subroutine sort_states_by_energy

    subroutine sort_space_by_proc(ilut_list, ilut_list_size, num_states_procs)

        ! And also output the number of states on each processor in the space.

        use load_balance_calcnodes, only: DetermineDetNode
        use MemoryManager, only: TagIntType, LogMemAlloc, LogMemDealloc

        integer, intent(in) :: ilut_list_size
        integer(n_int), intent(inout) :: ilut_list(0:,:)
        integer(MPIArg), intent(out) :: num_states_procs(0:nProcessors-1)
        integer(n_int), allocatable, dimension(:,:) :: temp_list
        integer, allocatable, dimension(:) :: proc_list
        integer :: nI(nel)
        integer :: i, ierr
        integer :: width, max_ind
        integer :: counter(0:nProcessors-1)
        integer(TagIntType) :: TempConTag, ProcListTag
        character(len=*), parameter :: t_r = "sort_space_by_proc"

        width = int(size(ilut_list, 1), MPIArg)
        max_ind = width - 1

        allocate(proc_list(ilut_list_size), stat=ierr)
        call LogMemAlloc('proc_list', ilut_list_size, sizeof_int, t_r, ProcListTag, ierr)

        allocate(temp_list(0:max_ind, ilut_list_size), stat=ierr)
        call LogMemAlloc('temp_list', ilut_list_size*width, size_n_int, t_r, &
                         TempConTag, ierr)

        num_states_procs = 0

        ! Create a list, proc_list, with the processor numbers of the corresponding iluts.
        do i = 1, ilut_list_size
            call decode_bit_det(nI, ilut_list(0:NIfTot,i))
            proc_list(i) = DetermineDetNode(nel,nI,0)
            num_states_procs(proc_list(i)) = int(num_states_procs(proc_list(i)) + 1, MPIArg)
        end do

        counter(0) = 0
        do i = 1, nProcessors-1
            counter(i) = counter(i-1) + num_states_procs(i-1)
        end do

        do i = 1, ilut_list_size
            counter(proc_list(i)) = counter(proc_list(i)) + 1
            temp_list(0:NIfTot, counter(proc_list(i))) = ilut_list(0:NIfTot,i)
        end do

        ilut_list(:,1:ilut_list_size) = temp_list(:,1:ilut_list_size)

        deallocate(temp_list, stat=ierr)
        deallocate(proc_list, stat=ierr)
        call LogMemDealloc(t_r, TempConTag, ierr)
        call LogMemDealloc(t_r, ProcListTag, ierr)

    end subroutine sort_space_by_proc

    subroutine fill_in_diag_helements()

        use Determinants, only: get_helement
        use FciMCData, only: Hii
        use global_det_data, only: set_det_diagH, tAccumEmptyDet
        use hphf_integrals, only: hphf_diag_helement
        use SystemData, only: tHPHF
        use bit_reps, only: extract_sign

        integer :: i
        integer :: nI(nel)
        real(dp) :: tmpH
        real(dp) :: sgn(lenof_sign)

        do i = 1, int(TotWalkers)
            call extract_sign(CurrentDets(:,i), sgn)
            if (IsUnoccDet(sgn) .and. (.not. tAccumEmptyDet(i))) cycle
            call decode_bit_det(nI, CurrentDets(:,i))

            if (tHPHF) then
                tmpH = hphf_diag_helement(nI, CurrentDets(:,i)) - Hii
            else
                tmpH = get_helement(nI, nI, 0) - Hii
            end if
            call set_det_diagh(i, tmpH)

        end do

    end subroutine fill_in_diag_helements

    subroutine write_core_space()

        use Parallel_neci, only: MPIBarrier
        use ParallelHelper, only: root
        use util_mod, only: get_free_unit

        integer :: i, k, iunit, ierr

        write(6,'(a35)') "Writing the core space to a file..."

        iunit = get_free_unit()

        ! Only let the root process write the states.
        if (iProcIndex == root) then
            open(iunit, file='CORESPACE', status='replace')

            do i = 1, determ_space_size
                do k = 0, NIfDBO
                    write(iunit, '(i24)', advance='no') core_space(k,i)
                end do
                write(iunit, '()')
            end do

            call neci_flush(iunit)
            close(iunit)
        end if

        call MPIBarrier(ierr)

    end subroutine write_core_space

    subroutine add_core_states_currentdets()

        ! And if the state is already present, simply set its flag.
        ! Also sort the states afterwards.

        use bit_rep_data, only: flag_initiator
        use bit_reps, only: set_flag
        use DetBitOps, only: ilut_lt, ilut_gt, DetBitLT
        use searching, only: BinSearchParts
        use sort_mod, only: sort

        integer :: i, run, comp, MinInd, PartInd, nwalkers
        logical :: tSuccess

        MinInd = 1
        nwalkers = int(TotWalkers,sizeof_int)

        do i = 1, determ_sizes(iProcIndex)

            if (nwalkers > 0) then
                ! If there is only one state in CurrentDets to check then BinSearchParts doesn't
                ! return the desired value for PartInd, so do this separately...
                if (MinInd == nwalkers) then
                    comp = DetBitLT(CurrentDets(:,MinInd), SpawnedParts(0:NIfTot,i), NIfDBO)
                    if (comp == 0) then
                        tSuccess = .true.
                        PartInd = MinInd
                    else if (comp == 1) then
                        tSuccess = .false.
                        PartInd = MinInd
                    else if (comp == -1) then
                        tSuccess = .false.
                        PartInd = MinInd - 1
                    end if
                else
                    call BinSearchParts(SpawnedParts(:,i), MinInd, nwalkers, PartInd, tSuccess)
                end if
            else
                tSuccess = .false.
                PartInd = 0
            end if

            if (tSuccess) then
                call set_flag(CurrentDets(:,PartInd), flag_deterministic)
                if (tTruncInitiator) then
                    do run = 1, inum_runs
                        call set_flag(CurrentDets(:,PartInd), get_initiator_flag_by_run(run))
                    end do
                end if
                MinInd = PartInd
            else
                ! Move all states below PartInd down one and insert the new state in the slot.
                CurrentDets(:, PartInd+2:nwalkers+1) = CurrentDets(:, PartInd+1:nwalkers)
                CurrentDets(:, PartInd+1) = SpawnedParts(0:NIfTot,i)
                nwalkers = nwalkers + 1
                MinInd = PartInd + 1
            end if

        end do

        call sort(CurrentDets(:,1:nwalkers), ilut_lt, ilut_gt)

        TotWalkers = int(nwalkers, int64)

    end subroutine add_core_states_currentdets

    subroutine add_core_states_currentdet_hash()

        ! This routine adds the core states in SpawnedParts into CurrentDets. For all
        ! such states already in CurrentDets, we want to keep the amplitude (which
        ! may have come from a popsfile).

        ! This routine is for when the hashed walker main list. In this case,
        ! as all core states are always kept in the list, it is beneficial to keep
        ! them at the top always. So, in this routine, we move the non-core states
        ! in CurrentDets to the end and add the new core states in the gaps.

        ! WARNING: If there are any determinants in CurrentDets on input which are
        ! unoccupied then, for this function to work correctly, the determiant
        ! *must* have an entry in the hash table. Otherwise, these determinants
        ! will end up being repeated in CurrentDets is they are core determinants.
        ! This isn't ideal because when the FCIQMC calculation starts, such
        ! unoccupied determinants should *not* be in the hash table. During this
        ! routine, such determinants will be removed from the hash table and so
        ! on output, everything will be fine and ready for the FCIQMC calculation
        ! to start.

        use bit_reps, only: set_flag, extract_sign, encode_sign
        use FciMCData, only: ll_node, indices_of_determ_states, HashIndex, nWalkerHashes
        use hash, only: clear_hash_table, FindWalkerHash

        integer :: i, hash_val, PartInd, nwalkers, i_non_core
        integer :: nI(nel)
        real(dp) :: walker_sign(lenof_sign)
        type(ll_node), pointer :: temp_node
        logical :: tSuccess
        integer :: ierr
        character(*), parameter :: this_routine = 'add_core_states_currentdet'
        real(dp), allocatable :: fvals(:,:)
        real(dp), allocatable :: apvals(:,:)

        nwalkers = int(TotWalkers,sizeof_int)

        ! Test that SpawnedParts is going to be big enough
        if (determ_sizes(iProcIndex) > MaxSpawned) then
#ifdef __DEBUG
            write(6,*) 'Spawned parts array will not be big enough for &
                       &Semi-Stochastic initialisation'
            write(6,*) 'Please increase MEMORYFACSPAWN'
#else
            write(iout,*) 'Spawned parts array will not be big enough for &
                       &Semi-Stochastic initialisation on task ', iProcIndex
            write(iout,*) 'Please increase MEMORYFACSPAWN'
#endif
            call stop_all(this_routine, "Insufficient memory assigned")
        end if

        ! we need to reorder the adaptive shift data, too
        if(tAutoAdaptiveShift) then
           ! the maximally required buffer size is the current size of the
           ! determinant list plus the size of the semi-stochastic space (in case
           ! all core-dets are new)
           allocate(fvals(2*inum_runs,(nwalkers+determ_sizes(iProcIndex))), stat = ierr)
           if(ierr.ne.0) call stop_all(this_routine, &
                "Failed to allocate buffer for adaptive shift data")
           allocate(apvals(lenof_sign+1,(nwalkers+determ_sizes(iProcIndex))), stat = ierr)
           if(ierr.ne.0) call stop_all(this_routine, &
                "Failed to allocate buffer for accumulated population data")
        endif

        ! First find which CurrentDet states are in the core space.
        ! The warning above refers to this bit of code: If a core determinant is not in the
        ! hash table then they won't be found here and the deterministic flag won't be set!
        do i = 1, determ_sizes(iProcIndex)

            tSuccess = .false.
            call decode_bit_det (nI, SpawnedParts(:,i))
            hash_val = FindWalkerHash(nI, nWalkerHashes)
            temp_node => HashIndex(hash_val)
            if (temp_node%ind /= 0) then
                do while (associated(temp_node))
                    if ( all(SpawnedParts(0:NIfDBO, i) == CurrentDets(0:NIfDBO, temp_node%ind)) ) then
                        tSuccess = .true.
                        PartInd = temp_node%ind
                        exit
                    end if
                    temp_node => temp_node%next
                end do
            end if
            nullify(temp_node)

            ! Core state i is in CurrentDets.
            if (tSuccess) then
                call set_flag(CurrentDets(:,PartInd), flag_deterministic)
                ! Copy the amplitude of the state across to SpawnedParts.
                call extract_sign(CurrentDets(:,PartInd), walker_sign)
                call encode_sign(SpawnedParts(:,i), walker_sign)
                if(tAutoAdaptiveShift) call cache_fvals(i,PartInd)
                if(tAccumPopsActive) call cache_apvals(i,PartInd)
            else
                ! This will be a new state added to CurrentDets.
                nwalkers = nwalkers + 1
                ! no auto-adaptive shift data available
                if(tAutoAdaptiveShift) fvals(:,i) = 0.0_dp
                if(tAccumPopsActive) apvals(:,i) = 0.0_dp
            end if

        end do
        ! Next loop through CurrentDets and move all non-core states to after the last
        ! core state slot in SpawnedParts.
        i_non_core = determ_sizes(iProcIndex)
        do i = 1, int(TotWalkers,sizeof_int)
            if (.not. test_flag(CurrentDets(:,i), flag_deterministic)) then
                i_non_core = i_non_core + 1

                ! Add a quick test in, to ensure that we don't overflow the
                ! spawned parts array...
                if (i_non_core > MaxSpawned) then
#ifdef __DEBUG
                    write(6,*) 'Spawned parts array too small for &
                               &semi-stochastic initialisation'
                    write(6,*) 'Please increase MEMORYFACSPAWN'
#else
                    write(iout,*) 'Spawned parts array too small for &
                               &semi-stochastic initialisation on task ', iProcIndex
                    write(iout,*) 'Please increase MEMORYFACSPAWN'
#endif
                    call stop_all(this_routine, 'Insufficient memory assigned')
                end if

                SpawnedParts(0:NIfTot,i_non_core) = CurrentDets(:,i)
                if(tAutoAdaptiveShift) call cache_fvals(i_non_core,i)
                if(tAccumPopsActive) call cache_apvals(i_non_core,i)
            end if
        end do
        ! Now copy all the core states in SpawnedParts into CurrentDets.
        ! Note that the amplitude in CurrentDets was copied across, so this is fine.
        do i = 1, nwalkers
            CurrentDets(:,i) = SpawnedParts(0:NIfTot,i)
            ! also re-order the adaptive shift data if auto-adapive shift is used
        end do
        if(tAutoAdaptiveShift) call set_tot_acc_spawns(fvals, nwalkers)
        if(tAccumPopsActive) call set_apvals(apvals, nwalkers)

        call clear_hash_table(HashIndex)

        ! Finally, add the indices back into the hash index array.
        do i = 1, nwalkers
            call extract_sign(CurrentDets(:,i), walker_sign)
            ! Don't add the determinant to the hash table if its unoccupied and not
            ! in the core space and not accumulated.
            if (IsUnoccDet(walker_sign) .and. (.not. test_flag(CurrentDets(:,i), flag_deterministic)) .and. .not. tAccumEmptyDet(i)) cycle
            call decode_bit_det(nI, CurrentDets(:,i))
            hash_val = FindWalkerHash(nI,nWalkerHashes)
            temp_node => HashIndex(hash_val)
            ! If the first element in the list has not been used.
            if (temp_node%ind == 0) then
                temp_node%ind = i
            else
                do while (associated(temp_node%next))
                    temp_node => temp_node%next
                end do
                allocate(temp_node%next)
                nullify(temp_node%next%next)
                temp_node%next%ind = i
            end if
            nullify(temp_node)

            ! These core states will always stay in the same position.
            if (i <= determ_sizes(iProcIndex)) indices_of_determ_states(i) = i
        end do

        TotWalkers = int(nwalkers, int64)

        contains

          subroutine cache_fvals(i,j)
            use global_det_data, only: get_acc_spawns, get_tot_spawns
            implicit none
            integer, intent(in) :: i,j
            integer :: run

            do run = 1, inum_runs
               fvals(run,i) = get_acc_spawns(j,run)
               fvals(run+inum_runs,i) = get_tot_spawns(j,run)
            end do
          end subroutine cache_fvals

          subroutine cache_apvals(i,j)
            use global_det_data, only: get_pops_sum, get_pops_iter
            implicit none
            integer, intent(in) :: i,j
            integer :: k

            do k = 1, lenof_sign
               apvals(k,i) = get_pops_sum(j,k)
            end do
            apvals(lenof_sign+1,i) = get_pops_iter(j)

          end subroutine cache_apvals
    end subroutine add_core_states_currentdet_hash

    subroutine return_most_populated_states(n_keep, largest_walkers, norm)

        ! Return the most populated states in CurrentDets on *this* processor only.
        ! Also return the norm of these states, if requested.

        use bit_reps, only: extract_sign
        use DetBitOps, only: sign_lt, sign_gt
        use sort_mod, only: sort

        integer, intent(in) :: n_keep
        integer(n_int), intent(out) :: largest_walkers(0:NIfTot, n_keep)
        real(dp), intent(out), optional :: norm
        integer :: i, j, smallest_pos, part_type
        real(dp) :: smallest_sign, sign_curr_real
        real(dp), dimension(lenof_sign) :: sign_curr, low_sign

        largest_walkers = 0_n_int
        smallest_sign = 0.0_dp
        smallest_pos = 1
        if (present(norm)) norm = 0.0_dp

        ! Run through all walkers on process.
        do i = 1, int(TotWalkers,sizeof_int)
            call extract_sign(CurrentDets(:,i), sign_curr)

#ifdef __CMPLX
            sign_curr_real = sqrt(sum(abs(sign_curr(1::2)))**2 + sum(abs(sign_curr(2::2)))**2)
#else
            if(tSignedRepAv) then
               sign_curr_real = real(abs(sum(sign_curr)),dp)
            else
               sign_curr_real = sum(real(abs(sign_curr),dp))
            endif
#endif
            if (present(norm)) norm = norm + (sign_curr_real**2.0)

            ! Is this determinant more populated than the smallest? First in
            ! the list is always the smallest.
            if (sign_curr_real > smallest_sign) then
                largest_walkers(:,smallest_pos) = CurrentDets(:,i)

                ! Instead of resorting, just find new smallest sign and position.
                call extract_sign(largest_walkers(:,1),low_sign)

#ifdef __CMPLX
                smallest_sign = sqrt(real(low_sign(1),dp)**2+real(low_sign(lenof_sign),dp)**2)
#else
                smallest_sign = sum(real(abs(low_sign),dp))
#endif

                smallest_pos = 1
                do j = 2, n_keep
                    call extract_sign(largest_walkers(:,j), low_sign)
#ifdef __CMPLX
                    sign_curr_real = sqrt(sum(real(low_sign**2,dp)))
#else
                    sign_curr_real = sum(real(abs(low_sign),dp))
#endif
                    if (sign_curr_real < smallest_sign .or. all(largest_walkers(:,j) == 0_n_int)) then
                        smallest_pos = j
                        smallest_sign = sign_curr_real
                    end if
                end do

            endif

        end do

        call sort(largest_walkers(:,1:n_keep), sign_lt, sign_gt)

    end subroutine return_most_populated_states

    subroutine return_largest_indices(n_keep, list_size, list, largest_indices)

        ! Return the indices of the largest elements in list.

        integer, intent(in) :: n_keep, list_size
        real(dp), intent(in) :: list(list_size)
        integer, intent(out) :: largest_indices(n_keep)
        integer :: i, j, ind, smallest_pos
        real(dp) :: smallest_sign, sign_curr, sign_curr_abs, low_sign

        largest_indices = 0
        smallest_sign = 0.0_dp
        smallest_pos = 1

        do i = 1, list_size
            sign_curr = list(i)
            sign_curr_abs = abs(sign_curr)

            if (sign_curr_abs > smallest_sign) then
                largest_indices(smallest_pos) = i

                low_sign = list(largest_indices(1))

                smallest_sign = abs(low_sign)

                smallest_pos = 1
                do j = 2, n_keep
                    ind = largest_indices(j)
                    if (ind == 0) then
                        low_sign = 0.0_dp
                    else
                        low_sign = list(ind)
                    end if

                    sign_curr_abs = abs(low_sign)

                    if (sign_curr_abs < smallest_sign .or. largest_indices(j) == 0) then
                        smallest_pos = j
                        smallest_sign = sign_curr_abs
                    end if
                end do
            end if
        end do

    end subroutine return_largest_indices

    subroutine start_walkers_from_core_ground(tPrintInfo)

        use bit_reps, only: encode_sign
        use davidson_semistoch, only: davidson_ss, perform_davidson_ss, destroy_davidson_ss
        use Parallel_neci, only: MPISumAll

        logical, intent(in) :: tPrintInfo
        integer :: i, counter, ierr
        real(dp) :: eigenvec_pop, eigenvec_pop_tot, pop_sign(lenof_sign)
        character(len=*), parameter :: t_r = "start_walkers_from_core_ground"
        type(davidson_ss) :: dc

        if (tPrintInfo) then
            write(6,'(a69)') "Using the deterministic ground state as initial walker configuration."
            write(6,'(a34)') "Performing Davidson calculation..."
            call neci_flush(6)
        end if

        ! Call the Davidson routine to find the ground state of the core space.
        call perform_davidson_ss(dc, .true.)

        if (tPrintInfo) then
            write(6,'(a30)') "Davidson calculation complete."
            write(6,'("Deterministic subspace correlation energy:",1X,f15.10)') dc%davidson_eigenvalue
            call neci_flush(6)
        end if

        ! We need to normalise this vector to have the correct 'number of walkers'.
        eigenvec_pop = 0.0_dp
        do i = 1, determ_sizes(iProcIndex)
            eigenvec_pop = eigenvec_pop + abs(dc%davidson_eigenvector(i))
        end do
        call MPISumAll(eigenvec_pop, eigenvec_pop_tot)

        if (tStartSinglePart) then
            dc%davidson_eigenvector = dc%davidson_eigenvector*InitialPart/eigenvec_pop_tot
        else
            dc%davidson_eigenvector = dc%davidson_eigenvector*InitWalkers/eigenvec_pop_tot
        end if

        ! Then copy these amplitudes across to the corresponding states in CurrentDets.
        counter = 0
        do i = 1, int(TotWalkers, sizeof_int)
            if (test_flag(CurrentDets(:,i), flag_deterministic)) then
                counter = counter + 1
                pop_sign = dc%davidson_eigenvector(counter)
                call encode_sign(CurrentDets(:,i), pop_sign)
            end if
        end do

        call destroy_davidson_ss(dc)

    end subroutine start_walkers_from_core_ground

    subroutine copy_core_dets_to_spawnedparts()

        ! This routine will copy all the core determinants *ON THIS PROCESS
        ! ONLY* to the SpawnedParts array.

        use load_balance_calcnodes, only: DetermineDetNode

        integer :: i, ncore, proc
        integer :: nI(nel)
        character (len=*), parameter :: t_r = "copy_core_dets_to_spawnedparts"

        ncore = 0
        SpawnedParts = 0_n_int

        do i = 1, determ_space_size
            call decode_bit_det(nI, core_space(:,i))
            proc = DetermineDetNode(nel,nI,0)
            if (proc == iProcIndex) then
                ncore = ncore + 1
                SpawnedParts(0:NIfTot,ncore) = core_space(:,i)
            end if
        end do

        if (ncore /= determ_sizes(iProcIndex)) call stop_all(t_r, "The number of &
            &core determinants counted is less than was previously counted.")

    end subroutine copy_core_dets_to_spawnedparts

    subroutine return_mp1_amp_and_mp2_energy(nI, ilut, ex, tParity, amp, energy_contrib)

        ! For a given determinant (input as nI), find the amplitude of it in the MP1 wavefunction.
        ! Also return the contribution from this determinant in the MP2 energy.

        ! To use this routine, generate an excitation from the Hartree-Fock determinant using the
        ! GenExcitations3 routine. This will return nI, ex and tParity which can be input into this
        ! routine.

        use Determinants, only: get_helement, GetH0Element3, GetH0Element4
        use FciMCData, only: ilutHF, HFDet, Fii
        use hphf_integrals, only: hphf_off_diag_helement
        use SystemData, only: tHPHF, tUEG

        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: ex(2,2)
        logical, intent(in) :: tParity
        real(dp), intent(out) :: amp, energy_contrib
        integer :: ic
        real(dp) :: hel, H0tmp, denom

        amp = 0.0_dp
        energy_contrib = 0.0_dp

        if (ex(1,2) == 0) then
            ic = 1
        else
            ic = 2
        end if

        if (tHPHF) then
            ! Assume since we are using HPHF that the alpha and
            ! beta orbitals of the same spatial orbital have the same
            ! fock energies, so can consider either.
            hel = hphf_off_diag_helement(HFDet, nI, iLutHF, ilut)
        else
            hel = get_helement(HFDet, nI, ic, ex, tParity)
        end if

        if (tUEG) then
            ! This will calculate the MP2 energies without having to use the fock eigenvalues.
            ! This is done via the diagonal determinant hamiltonian energies.
            H0tmp = getH0Element4(nI, HFDet)
        else
            H0tmp = getH0Element3(nI)
        end if

        ! If the relevant excitation from the Hartree-Fock takes electrons from orbitals
        ! (i,j) to (a,b), then denom will be equal to
        ! \epsilon_a + \epsilon_b - \epsilon_i - \epsilon_j
        ! as required in the denominator of the MP1 amplitude and MP2 energy.
        denom = Fii - H0tmp

        if (.not. (abs(denom) > 0.0_dp)) then
            call warning_neci("return_mp1_amp_and_mp2_energy", &
            "One of the determinants under consideration for the MP1 wave function is degenerate &
            &with the Hartree-Fock determinant. Degenerate perturbation theory has not been &
            &considered, but the amplitude of this determinant will be returned as huge(0.0_dp) &
            &so that it should be included in the space.")
            amp = huge(0.0_dp)
        else
            amp = hel/denom
            energy_contrib = (hel**2)/denom
        end if

    end subroutine return_mp1_amp_and_mp2_energy

    subroutine reinit_current_trial_amps()

        ! Recreate current trial amps, without using arrays such as trial_space
        ! and trial_wfs, which are deallocated after the first init_trial_wf
        ! call.

        use bit_rep_data, only: flag_trial, flag_connected
        use bit_reps, only: decode_bit_det, set_flag
        use FciMCData, only: CurrentDets, TotWalkers, tTrialHash, current_trial_amps, ntrial_excits
        use searching, only: hash_search_trial, bin_search_trial
        use SystemData, only: nel

        integer(int64) :: i
        integer :: nI(nel)
        HElement_t(dp) :: trial_amps(ntrial_excits)
        logical :: tTrial, tCon

        ! Don't do anything if this is called before the trial wave function
        ! initialisation.
        if (.not. allocated(current_trial_amps)) return

        current_trial_amps = 0.0_dp

        do i = 1, TotWalkers
            if (tTrialHash) then
                call decode_bit_det(nI, CurrentDets(:,i))
                call hash_search_trial(CurrentDets(:,i), nI, trial_amps, tTrial, tCon)
            else
                call bin_search_trial(CurrentDets(:,i), trial_amps, tTrial, tCon)
            end if

            ! Set the appropraite flag (if any). Unset flags which aren't
            ! appropriate, just in case.
            if (tTrial) then
                call set_flag(CurrentDets(:,i), flag_trial, .true.)
                call set_flag(CurrentDets(:,i), flag_connected, .false.)
            else if (tCon) then
                call set_flag(CurrentDets(:,i), flag_trial, .false.)
                call set_flag(CurrentDets(:,i), flag_connected, .true.)
            else
                call set_flag(CurrentDets(:,i), flag_trial, .false.)
                call set_flag(CurrentDets(:,i), flag_connected, .false.)
            end if

            ! Set the amplitude (which may be zero).
            current_trial_amps(:,i) = trial_amps
        end do

    end subroutine reinit_current_trial_amps

    subroutine end_semistoch()

        use FciMCData, only: partial_determ_vecs, full_determ_vecs, full_determ_vecs_av
        use FciMCData, only: PDetermTag, FDetermTag, FDetermAvTag, IDetermTag
        use FciMCData, only: indices_of_determ_states, core_ham_diag, hamiltonian
        use FciMCData, only: core_space, determ_sizes, determ_displs, HamTag
        use FciMCData, only: CoreSpaceTag
        use MemoryManager, only: LogMemDealloc
        use sparse_arrays, only: SparseCoreHamilTags, deallocate_sparse_ham, core_ht
        use sparse_arrays, only: core_connections, deallocate_core_hashtable
        use sparse_arrays, only: deallocate_sparse_matrix_int

        character(len=*), parameter :: t_r = "end_semistoch"
        integer :: ierr

        call deallocate_sparse_ham(sparse_core_ham, 'sparse_core_ham', SparseCoreHamilTags)

        call deallocate_core_hashtable(core_ht)

        call deallocate_sparse_matrix_int(core_connections)

        if (allocated(partial_determ_vecs)) then
            deallocate(partial_determ_vecs, stat=ierr)
            call LogMemDealloc(t_r, PDetermTag, ierr)
        end if
        if (allocated(full_determ_vecs)) then
            deallocate(full_determ_vecs, stat=ierr)
            call LogMemDealloc(t_r, FDetermTag, ierr)
        end if
        if (allocated(full_determ_vecs_av)) then
            deallocate(full_determ_vecs_av, stat=ierr)
            call LogMemDealloc(t_r, FDetermAvTag, ierr)
        end if
        if (allocated(indices_of_determ_states)) then
            deallocate(indices_of_determ_states, stat=ierr)
            call LogMemDealloc(t_r, IDetermTag, ierr)
        end if
        if (allocated(core_ham_diag)) then
            deallocate(core_ham_diag, stat=ierr)
!            call LogMemDealloc(t_r, IDetermTag, ierr)
        end if
        if (allocated(core_space)) then
            deallocate(core_space, stat=ierr)
            call LogMemDealloc(t_r, CoreSpaceTag, ierr)
        end if
        if (allocated(hamiltonian)) then
            deallocate(hamiltonian, stat=ierr)
            call LogMemDealloc(t_r, HamTag, ierr)
        end if
        if (allocated(determ_sizes)) then
            deallocate(determ_sizes, stat=ierr)
            if (ierr /= 0) write(6,'("Error when deallocating determ_sizes:",1X,i8)') ierr
        end if
        if (allocated(determ_displs)) then
            deallocate(determ_displs, stat=ierr)
            if (ierr /= 0) write(6,'("Error when deallocating determ_displs:",1X,i8)') ierr
        end if
        if (allocated(determ_last)) then
            deallocate(determ_last, stat=ierr)
            if (ierr /= 0) write(6,'("Error when deallocating determ_last:",1X,i8)') ierr
        end if

    end subroutine end_semistoch

end module semi_stoch_procs
