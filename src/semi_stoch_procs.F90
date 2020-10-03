#include "macros.h"

! Some general procedures, created for the semi-stochastic (and trial wavefunction) code.
! Some are just used for initialisation and others in the main FCIQMC loop.

module semi_stoch_procs

    use bit_rep_data, only: flag_deterministic, NIfD, NIfTot, test_flag, &
                            test_flag_multi, IlutBits

    use bit_reps, only: decode_bit_det, get_initiator_flag_by_run

    use CalcData

    use constants

    use orb_idx_mod, only: SpinOrbIdx_t

    use FciMCData, only: SpawnedParts, TotWalkers, CurrentDets, &
                         MaxSpawned, ilutRef, &
                         t_global_core_space, core_run

    use core_space_util, only: core_space_t, cs_replicas, min_pt, max_pt, &
                               deallocate_sparse_ham

    use Parallel_neci, only: iProcIndex, nProcessors, MPIArg

    use SystemData, only: nel, tHPHF, tGUGA, t_non_hermitian

    use procedure_pointers, only: shiftScaleFunction

    use hphf_integrals, only: hphf_diag_helement, hphf_off_diag_helement

    use Determinants, only: get_helement

    use sparse_arrays, only: approx_ham

    use procedure_pointers, only: shiftFactorFunction

    use timing_neci

    use bit_reps, only: encode_sign

    use hamiltonian_linalg, only: parallel_sparse_hamil_type

    use davidson_neci, only: DavidsonCalcType, perform_davidson, DestroyDavidsonCalc

    use FciMCData, only: DavidsonTag

    use MemoryManager, only: LogMemAlloc, LogMemDealloc, TagIntType

    use Parallel_neci, only: MPIScatterV

    use ParallelHelper, only: root

    use sparse_arrays, only: sparse_ham, hamil_diag, HDiagTag

    use shared_rhash, only: shared_rhash_t, shared_rht_lookup

    use sparse_arrays, only: SparseHamilTags, allocate_sparse_ham_row

    use matrix_util, only: print_matrix, eig

    use adi_data, only: tSignedRepAv

    use global_det_data, only: readFVals, readAPVals

    use LoggingData, only: tAccumPopsActive

    use gdata_io, only: gdata_io_t

    use LoggingData, only: t_print_core_info

    use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi

    use guga_bitrepops, only: init_csf_information

    implicit none

    ! Distinguishing value for 'use all runs'
    integer, parameter :: GLOBAL_RUN = -45

contains

    subroutine determ_projection()

        ! This subroutine gathers together partial_determ_vecs from each processor so
        ! that the full vector for the whole deterministic space is stored on each processor.
        ! It then performs the deterministic multiplication of the projector on this full vector.

        use FciMCData, only: SemiStoch_Comms_Time
        use FciMCData, only: SemiStoch_Multiply_Time
        use Parallel_neci, only: MPIBarrier, MPIAllGatherV
        use DetBitOps, only: DetBitEQ
        integer :: run
        integer :: i, j, part_type, c_run
        integer :: ierr
        integer(MPIArg) :: MPIerr
        real(dp) :: scaledDiagSft(inum_runs)

        do run = 1, size(cs_replicas)
            associate(rep => cs_replicas(run))
                call MPIBarrier(ierr)

                call set_timer(SemiStoch_Comms_Time)

                call MPIAllGatherV(rep%partial_determ_vecs, rep%full_determ_vecs, &
                                   rep%determ_sizes, rep%determ_displs)

                call halt_timer(SemiStoch_Comms_Time)

                call set_timer(SemiStoch_Multiply_Time)

                if (rep%determ_sizes(iProcIndex) >= 1) then

                    ! For the moment, we're only adding in these contributions when we need the energy
                    ! This will need refinement if we want to continue with the option of inst vs true full RDMs
                    ! (as in another CMO branch).

                    ! Perform the multiplication.

                    rep%partial_determ_vecs = 0.0_dp

#ifdef CMPLX_
                    block
                        integer :: r_pt, i_pt
                        do i = 1, rep%determ_sizes(iProcIndex)
                            do j = 1, rep%sparse_core_ham(i)%num_elements
                                do r_pt = rep%min_part(), rep%max_part(), 2
                                    i_pt = r_pt + 1
                                    rep%partial_determ_vecs(r_pt, i) = rep%partial_determ_vecs(r_pt, i) - &
                              Real(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs(r_pt, rep%sparse_core_ham(i)%positions(j)) + &
                                 Aimag(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs(i_pt, rep%sparse_core_ham(i)%positions(j))
                                    rep%partial_determ_vecs(i_pt, i) = rep%partial_determ_vecs(i_pt, i) - &
                             Aimag(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs(r_pt, rep%sparse_core_ham(i)%positions(j)) - &
                                  Real(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs(i_pt, rep%sparse_core_ham(i)%positions(j))
                                end do
                            end do
                        end do
                    end block
#else

                    do i = 1, rep%determ_sizes(iProcIndex)
                        do j = 1, rep%sparse_core_ham(i)%num_elements
                            rep%partial_determ_vecs(:, i) = rep%partial_determ_vecs(:, i) - &
                                           rep%sparse_core_ham(i)%elements(j) * rep%full_determ_vecs(:, rep%sparse_core_ham(i)%positions(j))
                        end do
                    end do
#endif

                    ! Now add shift*full_determ_vecs to account for the shift, not stored in
                    ! sparse_core_ham.
#ifdef CMPLX_
                    do i = 1, rep%determ_sizes(iProcIndex)
                        ! Only scale the shift for the corespace when the option is set
                        if (tCoreAdaptiveShift .and. tAdaptiveShift) then
                            ! scale the shift using the abs of this run's complex coefficient
                            scaledDiagSft(run) = &
                                shiftFactorFunction( &
                                rep%indices_of_determ_states(i), run, &
                                sqrt(rep%full_determ_vecs(min_pt, &
                                                          i + rep%determ_displs(iProcIndex))**2 + &
                                     rep%full_determ_vecs(max_pt, &
                                                          i + rep%determ_displs(iProcIndex))**2)) * DiagSft(run)
                        else
                            scaledDiagSft = DiagSft
                        end if

                        do part_type = 1, size(rep%partial_determ_vecs, dim=1)
                            ! Convert the index along partial_determ_vecs into a part_type
                            rep%partial_determ_vecs(part_type, i) = rep%partial_determ_vecs(part_type, i) + &
                                                                    scaledDiagSft(part_type_to_run(rep%min_part() + part_type - 1)) &
                                                                    * rep%full_determ_vecs(part_type, i + rep%determ_displs(iProcIndex))
                        end do
                    end do
#else
                    do i = 1, rep%determ_sizes(iProcIndex)
                        ! Only scale the shift for the corespace when the option is set
                        if (tCoreAdaptiveShift .and. tAdaptiveShift) then
                            ! Here, translate between positins in the full_determ_vecs
                            ! and replicas
                            if (rep%t_global) then
                                c_run = run
                            else
                                c_run = 1
                            end if
                            ! get the re-scaled shift accounting for undersampling error
                            scaledDiagSft(run) = DiagSft(run) * shiftFactorFunction( &
                                                 rep%indices_of_determ_states(i), run, &
                                                 abs(rep%full_determ_vecs(c_run, i + rep%determ_displs(iProcIndex))))
                        else
                            scaledDiagSft = DiagSft
                        end if
                        rep%partial_determ_vecs(:, i) = rep%partial_determ_vecs(:, i) + &
                                                        scaledDiagSft(run) * rep%full_determ_vecs(:, i + rep%determ_displs(iProcIndex))
                    end do
#endif

                    ! Now multiply the vector by tau to get the final projected vector.
                    rep%partial_determ_vecs = rep%partial_determ_vecs * tau

                    do i = 1, rep%determ_sizes(iProcIndex)
                        if (tSkipRef(run) .and. DetBitEQ(CurrentDets(:, rep%indices_of_determ_states(i)), iLutRef(:, run), nIfD)) then
                            rep%partial_determ_vecs(:, i) = 0.0_dp
                        end if
                    end do
                end if
            end associate

            call halt_timer(SemiStoch_Multiply_Time)
        end do
    end subroutine determ_projection

    subroutine determ_projection_kp_hamil(partial_vecs, full_vecs, rep)

        use FciMCData, only: Hii, SemiStoch_Comms_Time, SemiStoch_Multiply_Time
        use Parallel_neci, only: MPIBarrier, MPIAllGatherV

        real(dp), allocatable, intent(inout) :: partial_vecs(:, :)
        real(dp), allocatable, intent(inout) :: full_vecs(:, :)
        type(core_space_t), intent(inout) :: rep

        integer :: i, j, ierr, run

        call MPIBarrier(ierr)

        call set_timer(SemiStoch_Comms_Time)

        call MPIAllGatherV(partial_vecs, full_vecs, rep%determ_sizes, rep%determ_displs)

        call halt_timer(SemiStoch_Comms_Time)

        call MPIBarrier(ierr)

        call set_timer(SemiStoch_Multiply_Time)

        if (rep%determ_sizes(iProcIndex) >= 1) then
            ! Start with this because sparse_core_hamil has Hii taken off, but actually we
            ! don't want the projected Hamiltonian to be relative to the HF determinant.
            partial_vecs = Hii * full_vecs(:, rep%determ_displs(iProcIndex) + 1: &
                                           rep%determ_displs(iProcIndex) + rep%determ_sizes(iProcIndex))

            do i = 1, rep%determ_sizes(iProcIndex)
                do j = 1, rep%sparse_core_ham(i)%num_elements
                    partial_vecs(:, i) = partial_vecs(:, i) + &
                                         rep%sparse_core_ham(i)%elements(j) * full_vecs(:, rep%sparse_core_ham(i)%positions(j))
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
        use FciMCData, only: SemiStoch_Comms_Time
        use FciMCData, only: SemiStoch_Multiply_Time
        use Parallel_neci, only: MPIBarrier, MPIAllGatherV

        integer :: i, j, ierr, run, part_type

        do run = 1, size(cs_replicas)
            associate(rep => cs_replicas(run))
                call MPIBarrier(ierr)

                call set_timer(SemiStoch_Comms_Time)

                call MPIAllGatherV(rep%partial_determ_vecs, rep%full_determ_vecs, &
                                   rep%determ_sizes, rep%determ_displs)

                call halt_timer(SemiStoch_Comms_Time)

                call set_timer(SemiStoch_Multiply_Time)

                if (rep%determ_sizes(iProcIndex) >= 1) then

                    ! Perform the multiplication.

                    rep%partial_determ_vecs = 0.0_dp

#ifdef CMPLX_
                    block
                        integer :: r_pt, i_pt
                        do i = 1, rep%determ_sizes(iProcIndex)
                            do j = 1, rep%sparse_core_ham(i)%num_elements
                                do r_pt = rep%min_part(), rep%max_part(), 2
                                    i_pt = r_pt + 1
                                    rep%partial_determ_vecs(r_pt, i) = rep%partial_determ_vecs(r_pt, i) - &
                              Real(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs(r_pt, rep%sparse_core_ham(i)%positions(j)) + &
                                 Aimag(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs(i_pt, rep%sparse_core_ham(i)%positions(j))
                                    rep%partial_determ_vecs(i_pt, i) = rep%partial_determ_vecs(i_pt, i) - &
                             Aimag(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs(r_pt, rep%sparse_core_ham(i)%positions(j)) - &
                                  Real(rep%sparse_core_ham(i)%elements(j)) * rep%full_determ_vecs(i_pt, rep%sparse_core_ham(i)%positions(j))
                                end do
                            end do
                        end do
                    end block
#else
                    do i = 1, rep%determ_sizes(iProcIndex)
                        do j = 1, rep%sparse_core_ham(i)%num_elements
                            rep%partial_determ_vecs(:, i) = rep%partial_determ_vecs(:, i) - &
                                           rep%sparse_core_ham(i)%elements(j) * rep%full_determ_vecs(:, rep%sparse_core_ham(i)%positions(j))
                        end do
                    end do
#endif

                    ! Remove contribution from the diagonal elements, since the
                    ! propagator has no diagonal.
                    do i = 1, rep%determ_sizes(iProcIndex)
                        rep%partial_determ_vecs(:, i) = rep%partial_determ_vecs(:, i) + &
                                                        rep%core_ham_diag(i) * rep%full_determ_vecs(:, i + rep%determ_displs(iProcIndex))
                    end do

                    ! Now multiply the vector by tau to get the final projected vector.
                    rep%partial_determ_vecs = &
                        rep%partial_determ_vecs * tau

                    do i = 1, rep%determ_sizes(iProcIndex)
                        if (tSkipRef(run) .and. DetBitEQ(CurrentDets(:, rep%indices_of_determ_states(i)), iLutRef(:, run), nIfD)) then
                            rep%partial_determ_vecs(:, i) = 0.0_dp
                        end if
                    end do

                end if

                call halt_timer(SemiStoch_Multiply_Time)

            end associate
        end do
    end subroutine determ_projection_no_death

    subroutine determ_proj_approx()

        ! This subroutine gathers together partial_determ_vecs from each processor so
        ! that the full vector for the whole deterministic space is stored on each processor.
        ! It then performs the deterministic multiplication of the projector on this full vector.

        use FciMCData, only: SemiStoch_Comms_Time
        use FciMCData, only: SemiStoch_Multiply_Time
        use Parallel_neci, only: MPIBarrier, MPIAllGatherV
        use DetBitOps, only: DetBitEQ

        integer :: i, j, ierr, run, part_type, c_run
        character(*), parameter :: t_r = "determ_proj_approx"

        if (.not. t_global_core_space) then
            call stop_all(t_r, "Cannot do approximate projection with core-space replicas")
        end if
        associate(rep => cs_replicas(core_run))
            call MPIBarrier(ierr)

            call set_timer(SemiStoch_Comms_Time)

            call MPIAllGatherV(rep%partial_determ_vecs, rep%full_determ_vecs, &
                               rep%determ_sizes, rep%determ_displs)

            call halt_timer(SemiStoch_Comms_Time)

            call set_timer(SemiStoch_Multiply_Time)

            if (rep%determ_sizes(iProcIndex) >= 1) then

                ! For the moment, we're only adding in these contributions when we need the energy
                ! This will need refinement if we want to continue with the option of inst vs true full RDMs
                ! (as in another CMO branch).

                ! Perform the multiplication.

                rep%partial_determ_vecs = 0.0_dp

#ifdef CMPLX_
                do i = 1, rep%determ_sizes(iProcIndex)
                    do j = 1, approx_ham(i)%num_elements
                        rep%partial_determ_vecs(min_pt, i) = rep%partial_determ_vecs(min_pt, i) - &
                                              Real(approx_ham(i)%elements(j)) * rep%full_determ_vecs(min_pt, approx_ham(i)%positions(j)) + &
                                                 Aimag(approx_ham(i)%elements(j)) * rep%full_determ_vecs(max_pt, approx_ham(i)%positions(j))
                        rep%partial_determ_vecs(max_pt, i) = rep%partial_determ_vecs(max_pt, i) - &
                                             Aimag(approx_ham(i)%elements(j)) * rep%full_determ_vecs(min_pt, approx_ham(i)%positions(j)) - &
                                                  Real(approx_ham(i)%elements(j)) * rep%full_determ_vecs(max_pt, approx_ham(i)%positions(j))
                    end do
                end do
#else
                do i = 1, rep%determ_sizes(iProcIndex)
                    do j = 1, approx_ham(i)%num_elements
                        rep%partial_determ_vecs(:, i) = rep%partial_determ_vecs(:, i) - &
                                                        approx_ham(i)%elements(j) * rep%full_determ_vecs(:, approx_ham(i)%positions(j))
                    end do
                end do
#endif

                ! Now add shift*full_determ_vecs to account for the shift, not stored in
                ! approx_ham.
#ifdef CMPLX_
                do i = 1, rep%determ_sizes(iProcIndex)
                    do part_type = 1, lenof_sign
                        rep%partial_determ_vecs(part_type, i) = rep%partial_determ_vecs(part_type, i) + &
                                                           DiagSft(run) * rep%full_determ_vecs(part_type, i + rep%determ_displs(iProcIndex))
                    end do
                end do
#else
                do i = 1, rep%determ_sizes(iProcIndex)
                    rep%partial_determ_vecs(:, i) = rep%partial_determ_vecs(:, i) + &
                                                    DiagSft * rep%full_determ_vecs(:, i + rep%determ_displs(iProcIndex))
                end do
#endif

                ! Now multiply the vector by tau to get the final projected vector.
                rep%partial_determ_vecs = rep%partial_determ_vecs * tau

                do i = 1, rep%determ_sizes(iProcIndex)
                    do part_type = 1, rep_size
                        if (tSkipRef(run) .and. DetBitEQ(CurrentDets(:, rep%indices_of_determ_states(i)), iLutRef(:, run), nIfD)) then
                            rep%partial_determ_vecs(part_type, i) = 0.0_dp
                        end if
                    end do
                end do
            end if

            call halt_timer(SemiStoch_Multiply_Time)
        end associate

    end subroutine determ_proj_approx

    subroutine average_determ_vector()

        use FciMCData, only: Iter, IterRDMStart
        use FciMCData, only: PreviousCycles
        use LoggingData, only: RDMEnergyIter

        real(dp) :: iter_curr, iter_start_av
        integer :: run
        ! If this condition is met then RDM energies were added in on the
        ! previous iteration. We now want to start a new averaging block so
        ! that the same contributions aren't added in again later.
        if (mod(Iter + PreviousCycles - IterRDMStart, RDMEnergyIter) == 0) then
            do run = 1, size(cs_replicas)
                cs_replicas(run)%full_determ_vecs_av = 0.0_dp
            end do
            write(6, *) "Reset fdv av at iteration ", iter
        end if

        ! The current iteration, converted to a double precision real.
        iter_curr = real(Iter + PreviousCycles, dp)
        ! The iteration that this averaging block started on.
        iter_start_av = real(RDMEnergyIter * ((Iter + PreviousCycles - IterRDMStart) / RDMEnergyIter) + IterRDMStart, dp)

        ! Add in the current deterministic vector to the running average.
        do run = 1, size(cs_replicas)
            associate(rep => cs_replicas(run))
                rep%full_determ_vecs_av = (((iter_curr - iter_start_av) * rep%full_determ_vecs_av) + rep%full_determ_vecs) / &
                                          (iter_curr - iter_start_av + 1.0_dp)
            end associate
        end do

    end subroutine average_determ_vector

    !> Check whether an ilut belongs to the core space
    !> @param[in] ilut  ilut we want to check
    !> @param[in] nI  determinant corresponding to this ilut. Redundant, but is passed
    !!                for performance reasons (decoding is expensive and we likely
    !!                already know nI at this point)
    !> @return t_core  true if and only if ilut is in the core space
    function is_core_state(ilut, nI, run_) result(t_core)
        use hash, only: FindWalkerHash

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: nI(:)
        integer, intent(in), optional :: run_
        integer :: i, run
        logical :: t_core

        def_default(run, run_, core_run)
        if (t_global_core_space) run = core_run

        associate(rep => cs_replicas(run))
            call shared_rht_lookup(rep%core_ht, ilut, nI, rep%core_space, i, t_core)
        end associate
    end function is_core_state

    !> Check where an ilut is in the core space
    !> @param[in] ilut  ilut we want to check
    !> @param[in] nI  determinant corresponding to this ilut. Redundant, but is passed
    !!                for performance reasons (decoding is expensive and we likely
    !!                already know nI at this point)
    !> @return pos  position of ilut in the core space, 0 if ilut is not in the core space
    function core_space_pos(ilut, nI, run_) result(pos)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: nI(:)
        integer, intent(in), optional :: run_
        character(len=*), parameter :: t_r = "core_space_pos"

        integer :: pos, run
        logical :: t_core

        def_default(run, run_, core_run)
        if (t_global_core_space) run = core_run
        associate(rep => cs_replicas(run))
            call shared_rht_lookup(rep%core_ht, ilut, nI, rep%core_space, pos, t_core)
        end associate
        if (pos == 0) then
            call stop_all(t_r, "State not found in core hash table.")
        end if
    end function core_space_pos

    function check_determ_flag(ilut, run_) result(core_state)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in), optional :: run_
        logical :: core_state
        integer :: run

        def_default(run, run_, GLOBAL_RUN)
        if (t_global_core_space) run = core_run
        if (tSemiStochastic) then
            if (run == GLOBAL_RUN) then
                core_state = test_flag_multi(ilut, flag_deterministic)
            else
                core_state = test_flag(ilut, flag_deterministic(run))
            end if
        else
            core_state = .false.
        end if

    end function check_determ_flag

    subroutine recalc_core_hamil_diag(old_Hii, new_Hii)

        real(dp) :: old_Hii, new_Hii
        real(dp) :: Hii_shift
        integer :: i, j, run

        ! Only attempt this if we have already performed the semi-stochastic
        ! initialisation, in which case determ_sizes will have been allocated.
        do run = 1, size(cs_replicas)
            associate(rep => cs_replicas(run))
                if (allocated(rep%determ_sizes)) then
                    write(6, '(a56)') "Recalculating diagonal elements of the core Hamiltonian."

                    Hii_shift = old_Hii - new_Hii

                    do i = 1, rep%determ_sizes(iProcIndex)
                        do j = 1, rep%sparse_core_ham(i)%num_elements
                            if (rep%sparse_core_ham(i)%positions(j) == i + rep%determ_displs(iProcIndex)) then
                                rep%sparse_core_ham(i)%elements(j) = rep%sparse_core_ham(i)%elements(j) + Hii_shift
                            end if
                        end do
                    end do

                    rep%core_ham_diag = rep%core_ham_diag + Hii_shift
                end if
            end associate
        end do

    end subroutine recalc_core_hamil_diag

    subroutine generate_core_connections(rep)

        use DetBitOps, only: FindBitExcitLevel
        use Parallel_neci, only: MPIAllGatherV
        type(core_space_t), intent(inout) :: rep
        integer :: i, j, ic, counter, ierr
        integer :: Ex(2, nel)
        logical :: tSign
        integer(n_int), allocatable, dimension(:, :) :: temp_store
        integer(TagIntType) :: TempStoreTag
        character(len=*), parameter :: t_r = "calculate_det_hamiltonian_sparse"

        allocate(rep%core_connections(rep%determ_sizes(iProcIndex)))

        allocate(temp_store(0:NIfTot, rep%determ_space_size), stat=ierr)
        call LogMemAlloc('temp_store', maxval(rep%determ_sizes) * (NIfTot + 1), 8, t_r, &
                         TempStoreTag, ierr)

        ! Stick together the deterministic states from all processors, on all processors.
        call MPIAllGatherV(SpawnedParts(0:NIfTot, 1:rep%determ_sizes(iProcIndex)), temp_store, &
                           rep%determ_sizes, rep%determ_displs)

        ! Over all core states on this processor.
        do i = 1, rep%determ_sizes(iProcIndex)

            ! The number of non-zero elements in this array will be almost the same as in
            ! the core Hamiltonian array, except the diagonal element is not considered,
            ! so there will actually be one less.
            allocate(rep%core_connections(i)%elements(rep%sparse_core_ham(i)%num_elements - 1))
            allocate(rep%core_connections(i)%positions(rep%sparse_core_ham(i)%num_elements - 1))

            ! The total number of non-zero elements in row i.
            rep%core_connections(i)%num_elements = rep%sparse_core_ham(i)%num_elements - 1

            counter = 0
            do j = 1, rep%sparse_core_ham(i)%num_elements
                ! If not the diagonal element.
                if (rep%sparse_core_ham(i)%positions(j) /= i + rep%determ_displs(iProcIndex)) then
                    Ex = 0
                    Ex(1, 1) = nel
                    counter = counter + 1
                    ! The positions of the non-zero and non-diagonal elements in this row i.
                    rep%core_connections(i)%positions(counter) = rep%sparse_core_ham(i)%positions(j)

                    ! for the GUGA implementation this has to be changed in the
                    ! future. but since this routine is only called if we calc.
                    ! RDMs on the fly, i can postpone that until then.. todo
                    ic = FindBitExcitLevel(SpawnedParts(:, i), temp_store(:, rep%sparse_core_ham(i)%positions(j)))
                    call GetBitExcitation(SpawnedParts(0:NIfD, i), temp_store(0:NIfD, &
                                                                              rep%sparse_core_ham(i)%positions(j)), Ex, tSign)
                    if (tSign) then
                        ! Odd number of permutations. Minus the excitation level.
                        rep%core_connections(i)%elements(counter) = -ic
                    else
                        ! Even number of permutations. The excitation level.
                        rep%core_connections(i)%elements(counter) = ic
                    end if
                end if
            end do

        end do

        deallocate(temp_store, stat=ierr)
        call LogMemDealloc(t_r, TempStoreTag, ierr)

    end subroutine generate_core_connections

    subroutine store_whole_core_space(rep)

        use FciMCData, only: CoreSpaceTag
        use MemoryManager, only: LogMemAlloc
        use Parallel_neci, only: MPIAllGatherV, iProcIndex_inter, iProcIndex_intra, &
                                 mpi_comm_inter, mpi_comm_intra, nNodes, NodeLengths, iNodeIndex, MPIAllGather
        type(core_space_t), intent(inout) :: rep
        integer(MPIArg) :: MPIerr
        integer :: ierr
        character(len=*), parameter :: t_r = "store_whole_core_space"
        integer(MPIArg), allocatable :: sizes_this_node(:)
        integer(MPIArg) :: total_size_this_node
        integer(MPIArg), allocatable :: node_offsets(:), sizes_per_node(:)
        integer(MPIArg) :: proc_offset
        integer(MPIArg) :: core_width
        integer :: i
        integer(MPIArg) :: node_size, num_nodes
        integer(MPIArg) :: global_offset

        call mpi_comm_size(mpi_comm_intra, node_size, MPIerr)
        call mpi_comm_size(mpi_comm_inter, num_nodes, MPIerr)

        allocate(sizes_this_node(0:node_size - 1))
        allocate(node_offsets(0:num_nodes - 1))
        allocate(sizes_per_node(0:num_nodes - 1))

        call shared_allocate_mpi(rep%core_space_win, rep%core_space_direct, &
                                 (/int(1 + NIfTot, int64), int(rep%determ_space_size, int64)/))
        ! Convert from 1-based first dimension to 0-based first dimension as used in iluts
        rep%core_space(0:, 1:) => rep%core_space_direct(1:, 1:)
        call LogMemAlloc('core_space', maxval(rep%determ_sizes) * (NIfTot + 1), 8, t_r, &
                         rep%CoreSpaceTag, ierr)
        if (iProcIndex_intra == 0) rep%core_space = 0_n_int

        ! Write the core-space on this node into core_space
        call MPI_AllGather(rep%determ_sizes(iProcIndex), 1, MPI_INTEGER4, &
                           sizes_this_node, 1, MPI_INTEGER4, mpi_comm_intra, MPIerr)

        ! Get the intra-node offset
        proc_offset = 0
        do i = 1, iProcIndex_intra
            proc_offset = proc_offset + sizes_this_node(i - 1)
        end do

        ! Sum the size on this node
        total_size_this_node = sum(sizes_this_node)

        call MPI_AllGather(total_size_this_node, 1, MPI_INTEGER4, sizes_per_node, 1, MPI_INTEGER4, &
                           mpi_comm_inter, MPIerr)

        ! Get the inter-node offset
        node_offsets(0) = 0
        do i = 1, num_nodes - 1
            node_offsets(i) = node_offsets(i - 1) + sizes_per_node(i - 1)
        end do

        global_offset = node_offsets(iProcIndex_inter) + proc_offset

        call MPI_Win_Sync(rep%core_space_win, MPIerr)
        call MPI_Barrier(mpi_comm_intra, MPIerr)
        rep%core_space(0:NIfTot, (global_offset + 1): &
                       (global_offset + rep%determ_sizes(iProcIndex))) = &
            SpawnedParts(0:NIfTot, 1:rep%determ_sizes(iProcIndex))
        call MPI_Win_Sync(rep%core_space_win, MPIerr)
        call MPI_Barrier(mpi_comm_intra, MPIerr)
        call MPI_Win_Sync(rep%core_space_win, MPIerr)

        ! Multiply with message width (1+NIfTot)
        core_width = int(size(rep%core_space, dim=1), MPIArg)
        node_offsets = node_offsets * core_width
        total_size_this_node = total_size_this_node * core_width
        sizes_per_node = sizes_per_node * core_width

        ! Give explicit limits for SpawnedParts slice, as NIfTot is not nesc.
        ! equal to NIfBCast. (It may be longer)
        call MPI_AllGatherV(MPI_IN_PLACE, total_size_this_node, MPI_INTEGER8, rep%core_space, &
                            sizes_per_node, node_offsets, MPI_INTEGER8, mpi_comm_inter, MPIerr)

        ! And sync the shared window
        call MPI_Win_Sync(rep%core_space_win, MPIerr)
        call MPI_Barrier(mpi_comm_intra, MPIerr)
        ! Communicate the indices in the full vector at which the various processors take over, relative
        ! to the first index position in the vector (i.e. the array displs in MPI routines).
        call MPIAllGather(global_offset, rep%determ_displs, ierr)

        deallocate(sizes_per_node)
        deallocate(node_offsets)
        deallocate(sizes_this_node)

    end subroutine store_whole_core_space

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
            tot_num_states = int(num_states, MPIArg)
        end if

        if (tot_num_states <= target_num_states) return

        write(6, '(a32)') "Removing high energy orbitals..."

        ! Loop through all orbitals, from highest energy to lowest.
        do i = nBasis, 1, -1

            num_orbs_rmvd = nBasis - i + 1

            orb = BRR(i)

            bit = mod(orb - 1, bits_n_int)
            elem = (orb - 1 - bit) / bits_n_int

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
                states_rmvd_all_procs = int(states_rmvd_this_proc, MPIArg)
            end if

            ! If there are degenerate orbitals, then cycle to remove the
            ! degenerate orbitals too, before giving the program chance to quit.
            if (i > 1) then
                ! If the orbitals energies are the same:
                if (abs(Arr(i, 1) - Arr((i - 1), 1)) < 1.0e-12_dp) cycle
            end if

            if (tot_num_states - states_rmvd_all_procs <= target_num_states) exit

        end do

        ! Loop through the list and shuffle states down to fill in the gaps
        ! created above.
        counter = 0
        do i = 1, num_states
            ! If the state wasn't set to 0:
            if (.not. all(ilut_list(:, i) == 0_n_int)) then
                counter = counter + 1
                ilut_list(:, counter) = ilut_list(:, i)
            end if
        end do
        if (counter /= num_states - states_rmvd_this_proc) &
            call stop_all("remove_high_energy_orbs", &
                          "Wrong number of states found.")
        num_states = counter

        ! Finally, output information:
        write(6, '(a36)', advance='no') "The following orbitals were removed:"
        allocate(orbitals_rmvd(num_orbs_rmvd))
        orbitals_rmvd = BRR(nBasis - num_orbs_rmvd + 1:nBasis)
        call sort(orbitals_rmvd)
        do i = 1, num_orbs_rmvd
            write(6, '(i5)', advance='no') orbitals_rmvd(i)
            if (i == num_orbs_rmvd) write(6, '()', advance='yes')
        end do
        deallocate(orbitals_rmvd)
        write(6, '(i6,1x,a29)') states_rmvd_all_procs, "states were removed in total."
        write(6, '(i6,1x,a17)') tot_num_states - states_rmvd_all_procs, "states were kept."
        call neci_flush(6)

    end subroutine remove_high_energy_orbs

    subroutine sort_states_by_energy(ilut_list, num_states, tSortDoubles)

        ! Note: If requested, keep all doubles at the top, then sort by energy.

        use bit_reps, only: decode_bit_det
        use DetBitOps, only: FindBitExcitLevel
        use FciMCData, only: ilutHF
        use sort_mod, only: sort

        integer, intent(in) :: num_states
        integer(n_int), intent(inout) :: ilut_list(0:NIfTot, 1:num_states)
        logical, intent(in) :: tSortDoubles
        integer(n_int) :: temp_ilut(0:NIfTot)
        integer :: nI(nel)
        integer :: i, excit_level, num_sing_doub, block_size
        real(dp) :: energies(num_states)

        do i = 1, num_states
            call decode_bit_det(nI, ilut_list(:, i))
            if (tHPHF) then
                energies(i) = hphf_diag_helement(nI, ilut_list(:, i))
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
                ! for GUGA this would have to be changed, but apparently this
                ! function is never called in the rest of the code so
                ! ignore it for now..
                excit_level = FindBitExcitLevel(ilut_list(:, i), ilutHF)
                if (excit_level <= 2) then
                    num_sing_doub = num_sing_doub + 1
                    temp_ilut = ilut_list(:, i)

                    ! Move block of non-singles or doubles down one in ilut_list.
                    ilut_list(:, num_sing_doub + 1:num_sing_doub + block_size) = &
                        ilut_list(:, num_sing_doub:num_sing_doub + block_size - 1)

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

        integer, intent(in) :: ilut_list_size
        integer(n_int), intent(inout) :: ilut_list(0:, :)
        integer(MPIArg), intent(out) :: num_states_procs(0:nProcessors - 1)
        integer(n_int), allocatable, dimension(:, :) :: temp_list
        integer, allocatable, dimension(:) :: proc_list
        integer :: nI(nel)
        integer :: i, ierr
        integer :: width, max_ind
        integer :: counter(0:nProcessors - 1)
        integer(TagIntType) :: TempConTag, ProcListTag
        character(len=*), parameter :: t_r = "sort_space_by_proc"

        width = int(size(ilut_list, 1), MPIArg)
        max_ind = width - 1

        allocate(proc_list(ilut_list_size), stat=ierr)
        call LogMemAlloc('proc_list', ilut_list_size, sizeof_int, t_r, ProcListTag, ierr)

        allocate(temp_list(0:max_ind, ilut_list_size), stat=ierr)
        call LogMemAlloc('temp_list', ilut_list_size * width, size_n_int, t_r, &
                         TempConTag, ierr)

        num_states_procs = 0

        ! Create a list, proc_list, with the processor numbers of the corresponding iluts.
        do i = 1, ilut_list_size
            call decode_bit_det(nI, ilut_list(0:NIfTot, i))
            proc_list(i) = DetermineDetNode(nel, nI, 0)
            num_states_procs(proc_list(i)) = int(num_states_procs(proc_list(i)) + 1, MPIArg)
        end do

        counter(0) = 0
        do i = 1, nProcessors - 1
            counter(i) = counter(i - 1) + num_states_procs(i - 1)
        end do

        do i = 1, ilut_list_size
            counter(proc_list(i)) = counter(proc_list(i)) + 1
            temp_list(0:NIfTot, counter(proc_list(i))) = ilut_list(0:NIfTot, i)
        end do

        ilut_list(:, 1:ilut_list_size) = temp_list(:, 1:ilut_list_size)

        deallocate(temp_list, stat=ierr)
        deallocate(proc_list, stat=ierr)
        call LogMemDealloc(t_r, TempConTag, ierr)
        call LogMemDealloc(t_r, ProcListTag, ierr)

    end subroutine sort_space_by_proc

    subroutine fill_in_diag_helements()

        use FciMCData, only: Hii
        use global_det_data, only: set_det_diagH

        integer :: i
        integer :: nI(nel)
        real(dp) :: tmpH
        real(dp) :: sgn(lenof_sign)

        do i = 1, int(TotWalkers)
            call decode_bit_det(nI, CurrentDets(:, i))

            if (tHPHF) then
                tmpH = hphf_diag_helement(nI, CurrentDets(:, i)) - Hii
            else
                tmpH = get_helement(nI, nI, 0) - Hii
            end if
            call set_det_diagh(i, tmpH)

        end do

    end subroutine fill_in_diag_helements

    subroutine write_core_space(rep)

        use Parallel_neci, only: MPIBarrier
        use util_mod, only: get_free_unit
        type(core_space_t), intent(in) :: rep
        integer :: i, k, iunit, ierr

        write(6, '(a35)') "Writing the core space to a file..."

        iunit = get_free_unit()

        ! Only let the root process write the states.
        if (iProcIndex == root) then
            open(iunit, file='CORESPACE', status='replace')

            do i = 1, rep%determ_space_size
                do k = 0, nifd
                    write(iunit, '(i24)', advance='no') rep%core_space(k, i)
                end do
                write(iunit, '()')
            end do

            call neci_flush(iunit)
            close(iunit)
        end if

        call MPIBarrier(ierr)

    end subroutine write_core_space

    subroutine add_core_states_currentdets(run)

        ! And if the state is already present, simply set its flag.
        ! Also sort the states afterwards.

        use bit_reps, only: set_flag
        use DetBitOps, only: ilut_lt, ilut_gt, DetBitLT
        use searching, only: BinSearchParts
        use sort_mod, only: sort
        integer, intent(in) :: run
        integer :: i, comp, MinInd, PartInd, nwalkers
        logical :: tSuccess

        associate(rep => cs_replicas(run))
            MinInd = 1
            nwalkers = int(TotWalkers, sizeof_int)

            do i = 1, rep%determ_sizes(iProcIndex)

                if (nwalkers > 0) then
                    ! If there is only one state in CurrentDets to check then BinSearchParts doesn't
                    ! return the desired value for PartInd, so do this separately...
                    if (MinInd == nwalkers) then
                        comp = DetBitLT(CurrentDets(:, MinInd), SpawnedParts(0:NIfTot, i), nifd)
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
                        call BinSearchParts(SpawnedParts(:, i), MinInd, nwalkers, PartInd, tSuccess)
                    end if
                else
                    tSuccess = .false.
                    PartInd = 0
                end if

                if (tSuccess) then
                    call set_flag(CurrentDets(:, PartInd), flag_deterministic(run))
                    if (tTruncInitiator .and. t_core_inits) then
                        call set_flag(CurrentDets(:, PartInd), get_initiator_flag_by_run(run))
                    end if
                    MinInd = PartInd
                else
                    ! Move all states below PartInd down one and insert the new state in the slot.
                    CurrentDets(:, PartInd + 2:nwalkers + 1) = CurrentDets(:, PartInd + 1:nwalkers)
                    CurrentDets(:, PartInd + 1) = SpawnedParts(0:NIfTot, i)
                    nwalkers = nwalkers + 1
                    MinInd = PartInd + 1
                end if

            end do
        end associate
        call sort(CurrentDets(:, 1:nwalkers), ilut_lt, ilut_gt)

        TotWalkers = int(nwalkers, int64)

    end subroutine add_core_states_currentdets

    subroutine add_core_states_currentdet_hash(run)

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

        use bit_reps, only: set_flag, extract_sign
        use FciMCData, only: ll_node, HashIndex, nWalkerHashes
        use hash, only: clear_hash_table, FindWalkerHash
        use DetBitOps, only: tAccumEmptyDet

        integer, intent(in) :: run
        integer :: i, hash_val, PartInd, nwalkers, i_non_core
        integer :: nI(nel)
        real(dp) :: walker_sign(lenof_sign)
        type(ll_node), pointer :: temp_node
        logical :: tSuccess
        integer :: ierr
        character(*), parameter :: this_routine = 'add_core_states_currentdet'
        real(dp), allocatable :: gdata_buf(:, :)
        type(gdata_io_t) :: reorder_handler

        nwalkers = int(TotWalkers, sizeof_int)

        associate(rep => cs_replicas(run))
            ! Test that SpawnedParts is going to be big enough
            if (rep%determ_sizes(iProcIndex) > MaxSpawned) then
#ifdef DEBUG_
                write(6, *) 'Spawned parts array will not be big enough for &
                    &Semi-Stochastic initialisation'
                write(6, *) 'Please increase MEMORYFACSPAWN'
#else
                write(iout, *) 'Spawned parts array will not be big enough for &
                    &Semi-Stochastic initialisation on task ', iProcIndex
                write(iout, *) 'Please increase MEMORYFACSPAWN'
#endif
                call stop_all(this_routine, "Insufficient memory assigned")
            end if

            call reorder_handler%init_gdata_io(tAutoAdaptiveShift, tScaleBlooms, tAccumPopsActive, &
                                               2 * inum_runs, 1, lenof_sign + 1)
            ! we need to reorder the adaptive shift data, too
            ! the maximally required buffer size is the current size of the
            ! determinant list plus the size of the semi-stochastic space (in case
            ! all core-dets are new)
            allocate(gdata_buf(reorder_handler%entry_size(), (nwalkers + rep%determ_sizes(iProcIndex))), &
                      stat=ierr)
            if (ierr /= 0) call stop_all(this_routine, &
                                         "Failed to allocate buffer for global det data")

            ! First find which CurrentDet states are in the core space.
            ! The warning above refers to this bit of code: If a core determinant is not in the
            ! hash table then they won't be found here and the deterministic flag won't be set!
            do i = 1, rep%determ_sizes(iProcIndex)

                tSuccess = .false.
                call decode_bit_det(nI, SpawnedParts(:, i))
                hash_val = FindWalkerHash(nI, nWalkerHashes)
                temp_node => HashIndex(hash_val)
                if (temp_node%ind /= 0) then
                    do while (associated(temp_node))
                        if (all(SpawnedParts(0:nifd, i) == CurrentDets(0:nifd, temp_node%ind))) then
                            tSuccess = .true.
                            PartInd = temp_node%ind
                            exit
                        end if
                        temp_node => temp_node%next
                    end do
                end if
                nullify (temp_node)

                ! Core state i is in CurrentDets.
                if (tSuccess) then
                    call set_flag(CurrentDets(:, PartInd), flag_deterministic(run))
                    ! Copy the amplitude of the state across to SpawnedParts.
                    call extract_sign(CurrentDets(:, PartInd), walker_sign)
                    call encode_sign(SpawnedParts(:, i), walker_sign)
                    ! Add up the already set flags to those to be set
                    SpawnedParts(IlutBits%ind_flag, i) = ior(CurrentDets(IlutBits%ind_flag, PartInd), &
                                                             SpawnedParts(IlutBits%ind_flag, i))
                    ! Cache the accumulated global det data
                    call reorder_handler%write_gdata(gdata_buf, 1, PartInd, i)
                else
                    ! This will be a new state added to CurrentDets.
                    nwalkers = nwalkers + 1
                    ! no auto-adaptive shift data available
                    gdata_buf(:, i) = 0.0_dp
                end if

            end do
            ! Next loop through CurrentDets and move all non-core states to after the last
            ! core state slot in SpawnedParts.
            i_non_core = rep%determ_sizes(iProcIndex)
            do i = 1, int(TotWalkers, sizeof_int)
                if (.not. check_determ_flag(CurrentDets(:, i), run)) then
                    i_non_core = i_non_core + 1

                    ! Add a quick test in, to ensure that we don't overflow the
                    ! spawned parts array...
                    if (i_non_core > MaxSpawned) then
#ifdef DEBUG_
                        write(6, *) 'Spawned parts array too small for &
                            &semi-stochastic initialisation'
                        write(6, *) 'Please increase MEMORYFACSPAWN'
#else
                        write(iout, *) 'Spawned parts array too small for &
                            &semi-stochastic initialisation on task ', iProcIndex
                        write(iout, *) 'Please increase MEMORYFACSPAWN'
#endif
                        call stop_all(this_routine, 'Insufficient memory assigned')
                    end if

                    SpawnedParts(0:NIfTot, i_non_core) = CurrentDets(:, i)
                    call reorder_handler%write_gdata(gdata_buf, 1, i, i_non_core)
                end if
            end do
            ! Now copy all the core states in SpawnedParts into CurrentDets.
            ! Note that the amplitude in CurrentDets was copied across, so this is fine.
            do i = 1, nwalkers
                CurrentDets(:, i) = SpawnedParts(0:NIfTot, i)
            end do
            ! Re-assign the reordered global det data cached in gdata_buf
            call reorder_handler%read_gdata(gdata_buf, nwalkers)

            call clear_hash_table(HashIndex)

            ! Finally, add the indices back into the hash index array.
            do i = 1, nwalkers
                call extract_sign(CurrentDets(:, i), walker_sign)
                ! Don't add the determinant to the hash table if its unoccupied and not
                ! in the core space and not accumulated.
         if (IsUnoccDet(walker_sign) .and. (.not. check_determ_flag(CurrentDets(:, i))) .and. .not. tAccumEmptyDet(CurrentDets(:, i))) cycle
                call decode_bit_det(nI, CurrentDets(:, i))
                hash_val = FindWalkerHash(nI, nWalkerHashes)
                temp_node => HashIndex(hash_val)
                ! If the first element in the list has not been used.
                if (temp_node%ind == 0) then
                    temp_node%ind = i
                else
                    do while (associated(temp_node%next))
                        temp_node => temp_node%next
                    end do
                    allocate(temp_node%next)
                    nullify (temp_node%next%next)
                    temp_node%next%ind = i
                end if
                nullify (temp_node)

                ! These core states will always stay in the same position.
                if (i <= rep%determ_sizes(iProcIndex)) rep%indices_of_determ_states(i) = i
            end do
        end associate
        TotWalkers = int(nwalkers, int64)
    end subroutine add_core_states_currentdet_hash

    subroutine proc_most_populated_states(n_keep, run, &
                                          largest_walkers, opt_source, opt_source_size, norm)
        ! Return the most populated states in CurrentDets on *this* processor only.
        ! Also return the norm of these states, if requested.

        use bit_reps, only: extract_sign
        use DetBitOps, only: sign_lt, sign_gt
        use sort_mod, only: sort

        integer(int64), intent(in), optional :: opt_source_size
        integer(n_int), intent(in), optional :: opt_source(0:, 1:)
        integer, intent(in) :: n_keep, run
        integer(n_int), intent(out) :: largest_walkers(0:NIfTot, n_keep)
        real(dp), intent(out), optional :: norm
        integer :: i, j, smallest_pos, part_type
        real(dp) :: smallest_sign, sign_curr_real
        real(dp), dimension(lenof_sign) :: sign_curr, low_sign
        character(*), parameter :: this_routine = "return_most_populated_states"

        integer(n_int), allocatable :: loc_source(:, :)
        integer(int64) :: source_size

        if (present(opt_source)) then
            ASSERT(present(opt_source_size))

            source_size = opt_source_size
            ! ask Kai if I have to allocate
            allocate(loc_source(0:niftot, 1:source_size), &
                      source=opt_source(0:NIfTot, 1:source_size))
!             loc_source => opt_source

        else
            source_size = TotWalkers
            allocate(loc_source(0:niftot, 1:source_size), &
                      source=CurrentDets(0:NIfTot, 1:source_size))
!             loc_source => CurrentDets
        end if

        largest_walkers = 0_n_int
        smallest_sign = 0.0_dp
        smallest_pos = 1
        if (present(norm)) norm = 0.0_dp

        ! Run through all walkers on process.
        do i = 1, int(source_size, sizeof_int)
            call extract_sign(loc_source(:, i), sign_curr)

            sign_curr_real = core_space_weight(sign_curr, run)
            if (present(norm)) norm = norm + (sign_curr_real**2.0)

            ! Is this determinant more populated than the smallest? First in
            ! the list is always the smallest.
            if (sign_curr_real > smallest_sign) then
                largest_walkers(:, smallest_pos) = loc_source(:, i)

                ! Instead of resorting, just find new smallest sign and position.
                call extract_sign(largest_walkers(:, 1), low_sign)

                smallest_sign = core_space_weight(low_sign, run)

                smallest_pos = 1
                do j = 2, n_keep
                    call extract_sign(largest_walkers(:, j), low_sign)
                    sign_curr_real = core_space_weight(low_sign, run)
                    if (sign_curr_real < smallest_sign .or. all(largest_walkers(:, j) == 0_n_int)) then
                        smallest_pos = j
                        smallest_sign = sign_curr_real
                    end if
                end do

            end if

        end do

        call sort(largest_walkers(:, 1:n_keep), sign_lt_run, sign_gt_run)

    contains

        pure function sign_lt_run(ilutI, ilutJ) result(bLt)

            ! This is a comparison function between two bit strings of length
            ! 0:NIfTot, and will return true if absolute value of the sign of
            ! ilutI is less than ilutJ

            integer(n_int), intent(in) :: iLutI(0:), iLutJ(0:)
            logical :: bLt
            real(dp) :: SignI(lenof_sign), SignJ(lenof_sign)

            call extract_sign(ilutI, SignI)
            call extract_sign(ilutJ, SignJ)

            bLt = core_space_weight(signI, run) < core_space_weight(signJ, run)

        end function sign_lt_run

        pure function sign_gt_run(ilutI, ilutJ) result(bGt)

            ! This is a comparison function between two bit strings of length
            ! 0:NIfTot, and will return true if the abs sign of ilutI is greater
            ! than ilutJ

            integer(n_int), intent(in) :: iLutI(0:), iLutJ(0:)
            logical :: bGt
            real(dp) :: SignI(lenof_sign), SignJ(lenof_sign)

            call extract_sign(ilutI, SignI)
            call extract_sign(ilutJ, SignJ)

            bGt = core_space_weight(signI, run) > core_space_weight(signJ, run)

        end function sign_gt_run

    end subroutine proc_most_populated_states

!>  @brief
!>      Return the most populated states over all processors.
!>
!>  @author Oskar Weser
!>
!>  @details
!>  Reducing version of `proc_most_populated_states`, which works per process.
!>  Return as many states as the size of largest_walkers.
!>  Returns the norm as well, if requested.
!>  @param[out] largest_walkers, Array of most `n_keep` most populated states.
    subroutine global_most_populated_states(n_keep, run, largest_walkers, norm, rank_of_largest)
        use Parallel_neci, only: MPISumAll, MPIAllReduceDatatype, MPIBCast
        use bit_reps, only: extract_sign

        integer, intent(in) :: n_keep, run
        integer(n_int), intent(out) :: largest_walkers(0:NIfTot, n_keep)
        real(dp), intent(out), optional :: norm
        integer, intent(out), optional :: rank_of_largest(n_keep)
        character(*), parameter :: this_routine = 'global_most_populated_states'

        integer(n_int), allocatable :: proc_largest_walkers(:, :)
        integer, allocatable :: rank_of_largest_(:)

        largest_walkers = 0_n_int
        allocate(proc_largest_walkers(0:NIfTot, n_keep), source=0_n_int)
        block
            real(dp) :: proc_norm, all_norm
            call proc_most_populated_states( &
                n_keep, run, proc_largest_walkers, CurrentDets, TotWalkers, proc_norm)
            if (present(norm)) then
                call MpiSumAll(proc_norm, all_norm)
                norm = sqrt(all_norm)
            end if
        end block

        allocate(rank_of_largest_(n_keep))
        block
            real(dp) :: high_sign, curr_sign(lenof_sign)
            integer :: high_pos
            integer :: i, j

            fill_largest_walkers: do i = 1, n_keep
                high_sign = 0.0_dp
                high_pos = 1
                find_largest_sign_per_proc: do j = n_keep, 1, -1
                    if (any(proc_largest_walkers(:, j) /= 0)) then
                        call extract_sign(proc_largest_walkers(:, j), curr_sign)
                        high_pos = j
#ifdef CMPLX_
                        high_sign = sqrt(sum(abs(curr_sign(1::2)))**2 &
                                         + sum(abs(curr_sign(2::2)))**2)
#else
                        high_sign = sum(real(abs(curr_sign), dp))
#endif
                        exit find_largest_sign_per_proc
                    end if
                end do find_largest_sign_per_proc

                block
                    real(dp) :: reduce_in(2, 1), reduce_out(2, 1)
                    reduce_in = reshape([high_sign, real(iProcIndex, dp)], shape(reduce_in))
                    call MPIAllReduceDatatype( &
                        reduce_in, size(reduce_in, 2), MPI_MAXLOC, MPI_2DOUBLE_PRECISION, reduce_out)
                    ! Now, reduce_out(2, :) has the rank of the largest weighted determinant
                    rank_of_largest_(i) = nint(reduce_out(2, 1))
                end block

                block
                    integer(n_int) :: HighestDet(0:NIfTot)
                    if (iProcIndex == rank_of_largest_(i)) then
                        HighestDet(0:NIfTot) = proc_largest_walkers(:, high_pos)
                    end if
                    call MPIBCast(HighestDet(0:NIfTot), size(HighestDet), rank_of_largest_(i))
                    largest_walkers(0:NIfTot, i) = HighestDet(:)
                end block

                ! Zeroing essentially deletes the element because we search
                ! for first nonzero from the end. Also no resorting is required.
                if (iProcIndex == rank_of_largest_(i)) then
                    proc_largest_walkers(:, high_pos) = 0_n_int
                end if
            end do fill_largest_walkers
        end block

        if (present(rank_of_largest)) rank_of_largest = rank_of_largest_

    end subroutine

    !> Weight function for picking the most populated states. Trivial in
    !! single run mode, but multiple options exist in mneci
    pure function core_space_weight(sign_curr, run) result(sign_curr_real)
        real(dp), intent(in) :: sign_curr(lenof_sign)
        integer, intent(in) :: run
        real(dp) :: sign_curr_real
#ifdef CMPLX_
        sign_curr_real = sqrt(sum(abs(sign_curr(1::2)))**2 + sum(abs(sign_curr(2::2)))**2)
#else
        if (tSignedRepAv) then
            sign_curr_real = real(abs(sum(sign_curr)), dp)
        else
            if (t_global_core_space .or. run == GLOBAL_RUN) then
                sign_curr_real = real(sum(abs(sign_curr)), dp)
            else
                sign_curr_real = mag_of_run(sign_curr, run)
            end if
        end if
#endif

    end function core_space_weight

    !> Specialized routine that returns the number of determinants that are going into
    !! the core-space on this processor. Once the leading determinants have been obtained on each
    !! processor, this function requires the minimal and maximal populations among these for all
    !! processors (requires previous MPI_Gather), then the procedure to determine the core-space
    !! size on this processor is as follows:
    !! 1) Get the minimum of the maximal populations, then count the number of determinants
    !!    above this population. This count is then broadcasted to the other procs, and the
    !!    total number of determinants above the smallest maximum is determined. If it is smaller
    !!    than the core-space size, these determinants are put into the core-space, else
    !!    we repeat with the second smallest of the maximal populations, and so on.
    !! 2) Get the maximum of the minimal populations, then count the number of determinants below
    !!    this population. This count is then broadcasted to the other procs, and the total number
    !!    of determinants below the largest minimum is determined. If the number of determinants
    !!    that are remaining (i.e. larger than the largest minimum and smaller than the smallest
    !!    maximum) is sufficient to fill up the core-space (in particular, the smallest maximum
    !!    has to be bigger than the largest minimum), the small determinants are discarded. Else,
    !!    we repeat this with the second largest minimum, and so on.
    !! 3) From the remaining determinants, each processor contributes a share that equals to the
    !!    ratio of the remaining determinants on this proc to the total remaining determinants
    !> @param[in] n_keep  core-space size
    !> @param[in] min_vals  minimal population of the canditates per processor
    !> @param[in] max_vals  maximal population of the canditates per processor
    !> @param[in] lengths  number of candidates per processor
    !> @param[in] list  candidates on this processor
    !> @param[out] n_dets_this_proc  number of core-space determinants on this processor
    subroutine return_proc_share(n_keep, min_vals, max_vals, lengths, list, n_dets_this_proc)
        use util_mod, only: binary_search_first_ge
        real(dp), intent(inout) :: max_vals(0:nProcessors - 1), min_vals(0:nProcessors - 1)
        real(dp), intent(in) :: list(:)
        integer, intent(in) :: n_keep, lengths(0:nProcessors - 1)
        integer, intent(out) :: n_dets_this_proc

        integer :: sum_max, sum_min
        real(dp) :: min_max, max_min
        integer :: dets_left, pool_left
        real(dp) :: total_pool, ip_ratio
        integer :: n_max(0:nProcessors - 1), n_min(0:nProcessors - 1)
        integer :: missing, n_full, i

        dets_left = -1
        ! Increase the cutoff until our selection is small enough
        do while (dets_left < 0)
            call get_pp_ex(min_max, n_max, sum_max, max_vals)
            ! Number of determinants left when keeping the maximal ones
            dets_left = n_keep - sum_max
        end do

        ! Definitely take these determinants
        n_dets_this_proc = n_max(iProcIndex)

        max_min = min_max + 1
        total_pool = -1
        ! Reduce the cutoff until we are below the min_max
        do while (max_min > min_max .or. total_pool < dets_left)
            call get_pp_ex(max_min, n_min, sum_min, min_vals, t_max=.true.)
            ! Size of the pool left when not keeping the minimal ones ( has to be at least
            ! big enough to fill the core-space)
            total_pool = sum(lengths) - sum_max - sum_min
        end do

        ! If the corespace consists of all chosen determinants, the remaining pool might be 0
        ! -> no further action, take all determinants
        if (total_pool > 0) then
            ! Number of available dets on this proc after removing min/max
            pool_left = lengths(iProcIndex) - n_max(iProcIndex) - n_min(iProcIndex)
            ! Ratio of available dets on this proc vs. in totap
            ip_ratio = pool_left / real(total_pool, dp)
            ! If any further dets have to be picked, get them from all procs weighted with the pool sizes
            n_dets_this_proc = n_dets_this_proc + int(ip_ratio * dets_left)
        end if

    contains

        subroutine get_pp_ex(ex, n_ex, sum_ex, vals, t_max)
            use Parallel_neci, only: MPIAllGather
            real(dp), intent(out) :: ex
            integer, intent(out) :: n_ex(0:nProcessors - 1), sum_ex
            real(dp), intent(inout) :: vals(0:nProcessors - 1)
            logical, intent(in), optional :: t_max
            integer :: ex_ind, n_ex_loc
            integer :: ierr
            real(dp) :: pre
            logical :: t_max_

            def_default(t_max_, t_max, .false.)

            if (t_max_) then
                pre = -1.0
            else
                pre = 1.0
            end if
            ! Get the smallest value of the per-proc max
            ex_ind = minloc(pre * vals, dim=1) - 1
            ex = vals(ex_ind)
            ! Invalidate this value, such that the next call finds the second smallest value and so on
            vals(ex_ind) = pre * sum(abs(vals))
            ! Now, get the location of the first element above the extremum
            if (size(list) > 0) then
                n_ex_loc = binary_search_first_ge(list, ex)
            else
                ! it might be possible that a proc is empty (has no candidates)
                ! in this case, never find anything above/below global extremal values
                n_ex_loc = -1
            end if
            ! If no such element exists, return 0 on this proc
            if (n_ex_loc < 0) then
                n_ex_loc = 0
            else if (t_max_) then
                ! From the position, get the number of elements below the extremum (for max_min)
                n_ex_loc = n_ex_loc - 1
            else
                ! Or above the extremum (for min_max)
                n_ex_loc = lengths(iProcIndex) - n_ex_loc + 1
            end if

            call MPIAllGather(n_ex_loc, n_ex, ierr)
            ! Check if the maximum pop is already sufficient
            sum_ex = sum(n_ex)

        end subroutine get_pp_ex

    end subroutine return_proc_share

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

    subroutine start_walkers_from_core_ground(tPrintInfo, run)
        use davidson_semistoch, only: davidson_ss, perform_davidson_ss, destroy_davidson_ss
        use Parallel_neci, only: MPISumAll
        implicit none

        logical, intent(in) :: tPrintInfo
        integer, intent(in) :: run
        integer :: nI(nel)
        integer :: i, counter, ierr
        real(dp) :: eigenvec_pop, eigenvec_pop_tot, pop_sign(lenof_sign)
        integer(n_int) :: tmp
        character(len=*), parameter :: t_r = "start_walkers_from_core_ground"

        type(davidson_ss) :: dc

        if (tPrintInfo) then
            write(6, '(a69)') "Using the deterministic ground state as initial walker configuration."
            write(6, '(a34)') "Performing Davidson calculation..."
            call neci_flush(6)
        end if

        ! Call the Davidson routine to find the ground state of the core space.
        call perform_davidson_ss(dc, .true., run)

        if (tPrintInfo) then
            write(6, '(a30)') "Davidson calculation complete."
            write(6, '("Deterministic subspace correlation energy:",1X,f15.10)') dc%davidson_eigenvalue
            call neci_flush(6)
        end if

        associate(rep => cs_replicas(run))
            ! We need to normalise this vector to have the correct 'number of walkers'.
            eigenvec_pop = 0.0_dp
            do i = 1, rep%determ_sizes(iProcIndex)
                eigenvec_pop = eigenvec_pop + abs(dc%davidson_eigenvector(i))
            end do
            call MPISumAll(eigenvec_pop, eigenvec_pop_tot)

            if (tStartSinglePart) then
                dc%davidson_eigenvector = dc%davidson_eigenvector * InitialPart / eigenvec_pop_tot
            else
                dc%davidson_eigenvector = dc%davidson_eigenvector * InitWalkers / eigenvec_pop_tot
            end if

            ! Then copy these amplitudes across to the corresponding states in CurrentDets.
            counter = 0
            do i = 1, int(TotWalkers, sizeof_int)
                if (check_determ_flag(CurrentDets(:, i), run)) then
                    counter = counter + 1
                    pop_sign = dc%davidson_eigenvector(counter)
                    call decode_bit_det(nI, CurrentDets(:, i))
                    tmp = transfer(pop_sign(rep%min_part():rep%max_part()), tmp)
                    CurrentDets(IlutBits%ind_pop + rep%min_part() - 1:IlutBits%ind_pop + rep%max_part() - 1, i) = tmp
                end if
            end do
        end associate

        call destroy_davidson_ss(dc)

    end subroutine start_walkers_from_core_ground

    subroutine start_walkers_from_core_ground_nonhermit(tPrintInfo, run)
        use bit_reps, only: encode_sign
        implicit none

        logical, intent(in) :: tPrintInfo
        integer, intent(in) :: run
        integer :: i, counter, ierr
        integer :: nI(nel)
        real(dp), allocatable :: e_values(:)
        HElement_t(dp), allocatable :: e_vectors(:, :)
        integer(n_int) :: tmp
        real(dp) :: eigenvec_pop, pop_sign(lenof_sign)
        character(len=*), parameter :: t_r = "start_walkers_from_core_ground_nonhermit"

        if (tPrintInfo) then
            write(6, '(a69)') "Using the deterministic ground state as initial walker configuration."
            write(6, '(a53)') "Performing diagonalization of non-Hermitian matrix..."
            call neci_flush(6)
        end if

        associate(rep => cs_replicas(run))
            ! Call the non-Hermitian diagonalizer to find the ground state of the core space.
            call diagonalize_core_non_hermitian(e_values, e_vectors, rep)

            if (tPrintInfo) then
                write(6, '("Energies of the deterministic subspace:")')
                write(6, *) e_values(1:rep%determ_space_size)
                call neci_flush(6)
            end if

            ! We need to normalise this vector to have the correct 'number of walkers'.
            eigenvec_pop = 0.0_dp
            do i = 1, rep%determ_space_size
                eigenvec_pop = eigenvec_pop + abs(e_vectors(i, 1))
            end do

            if (tStartSinglePart) then
                e_vectors(:, 1) = e_vectors(:, 1) * InitialPart / eigenvec_pop
            else
                e_vectors(:, 1) = e_vectors(:, 1) * InitWalkers / eigenvec_pop
            end if

            write(6, *) 'The ground state vector:'
            write(6, *) e_vectors(:, 1)
            ! Then copy these amplitudes across to the corresponding states in CurrentDets.
            counter = 0
            do i = 1, iProcIndex
                counter = counter + rep%determ_sizes(i - 1)
            end do
            do i = 1, rep%determ_space_size !int(TotWalkers, sizeof_int)
                if (check_determ_flag(CurrentDets(:, i), run)) then
                    counter = counter + 1
                    pop_sign = e_vectors(counter, 1)
                    tmp = transfer(pop_sign(min_pt:max_pt), tmp)
                    CurrentDets(IlutBits%ind_pop + min_part_type(run) - 1:IlutBits%ind_pop + max_part_type(run) - 1, i) = tmp
                end if
            end do
        end associate

        deallocate(e_values, e_vectors)

    end subroutine start_walkers_from_core_ground_nonhermit

    subroutine diagonalize_core(e_value, e_vector, rep)
        real(dp), intent(out)  :: e_value
        HElement_t(dp), intent(out), allocatable :: e_vector(:)
        type(core_space_t), intent(in) :: rep
        type(DavidsonCalcType) :: davidsonCalc
        integer :: ierr
        character(*), parameter :: t_r = "diagonalize_core"

        call create_sparse_ham_from_core(rep)

        ! Call the Davidson routine to find the ground state of the core space.
        call perform_davidson(davidsonCalc, parallel_sparse_hamil_type, .true.)

        e_value = davidsonCalc%davidson_eigenvalue
        allocate(e_vector(rep%determ_space_size))
        e_vector = davidsonCalc%davidson_eigenvector

        write(6, '(a30)') "Davidson calculation complete."
        write(6, '("Deterministic subspace correlation energy:",1X,f15.10)') &
            e_value

        call neci_flush(6)

        call DestroyDavidsonCalc(davidsonCalc)
        ! call LogMemDealloc(t_r, DavidsonTag, ierr)
        deallocate(hamil_diag, stat=ierr)
        call LogMemDealloc(t_r, HDiagTag, ierr)
        call deallocate_sparse_ham(sparse_ham, SparseHamilTags)

    end subroutine diagonalize_core

    subroutine create_sparse_ham_from_core(rep)
        type(core_space_t), intent(in) :: rep
        character(*), parameter :: t_r = "create_sparse_ham_from_core"
        integer :: ierr, i

        ! Create the arrays used by the Davidson routine.
        ! First, the whole Hamiltonian in sparse form.
        if (allocated(sparse_ham)) deallocate(sparse_ham)
        if (allocated(SparseHamilTags)) deallocate(SparseHamilTags)
        allocate(sparse_ham(rep%determ_sizes(iProcIndex)))
        allocate(SparseHamilTags(2, rep%determ_sizes(iProcIndex)))
        do i = 1, rep%determ_sizes(iProcIndex)
            call allocate_sparse_ham_row(sparse_ham, i, rep%sparse_core_ham(i)%num_elements, "sparse_ham", SparseHamilTags(:, i))
            sparse_ham(i)%elements = rep%sparse_core_ham(i)%elements
            sparse_ham(i)%positions = rep%sparse_core_ham(i)%positions
            sparse_ham(i)%num_elements = rep%sparse_core_ham(i)%num_elements
        end do

        ! Next create the diagonal used by Davidson by copying the core one.
        if (allocated(hamil_diag)) deallocate(hamil_diag)
        allocate(hamil_diag(rep%determ_sizes(iProcIndex)), stat=ierr)
        call LogMemAlloc('hamil_diag', int(rep%determ_sizes(iProcIndex), sizeof_int), 8, t_r, HDiagTag, ierr)
        hamil_diag = rep%core_ham_diag

    end subroutine create_sparse_ham_from_core

    subroutine diagonalize_core_non_hermitian(e_values, e_vectors, rep)
        type(core_space_t), intent(in) :: rep
        real(dp), allocatable, intent(out) :: e_values(:)
        HElement_t(dp), allocatable :: e_vectors(:, :)
        HElement_t(dp), allocatable :: full_H(:, :)
        integer i, nI(nel), space_size

        ! if the Hamiltonian is non-hermitian we cannot use the
        ! standard Lanzcos or Davidson routines. so:
        ! build the full Hamiltonian
        call calc_determin_hamil_full(full_H, rep)

        if (t_print_core_info) then
            root_print "The determinants are"
            root_print "semistochastic basis:"
            if_root
            do i = 1, rep%determ_space_size
                call decode_bit_det(nI, rep%core_space(:, i))
                print *, nI
            end do
            end_if_root

            root_print "deterministic hamiltonian:"
            if_root
            call print_matrix(full_H)
            end_if_root
        end if

        allocate(e_values(size(full_H, 1)))
        allocate(e_vectors(size(full_H, 1), size(full_H, 1)))
        e_values = 0.0_dp
        e_vectors = 0.0_dp

        call eig(full_H, e_values, e_vectors)

        ! maybe we also want to start from a different eigenvector in
        ! this case? this would be practial for the hubbard problem case..
        root_print "Full diagonalisation for non-hermitian Hamiltonian completed!"

    end subroutine diagonalize_core_non_hermitian

    subroutine calc_determin_hamil_full(hamil, rep)
        use guga_data, only: ExcitationInformation_t
        use guga_excitations, only: calc_guga_matrix_element
        type(core_space_t) :: rep
        type(ExcitationInformation_t) :: excitInfo

        HElement_t(dp), allocatable, intent(out) :: hamil(:, :)
        integer :: i, j, nI(nel), nJ(nel)

        allocate(hamil(rep%determ_space_size, rep%determ_space_size))

        hamil = h_cast(0.0_dp)

        do i = 1, rep%determ_space_size
            call decode_bit_det(nI, rep%core_space(:, i))

            call init_csf_information(rep%core_space(0:nifd, i))

            if (tHPHF) then
                hamil(i, i) = hphf_diag_helement(nI, rep%core_space(:, i))
            else
                hamil(i, i) = get_helement(nI, nI, 0)
            end if

            do j = 1, rep%determ_space_size

                if (i == j) cycle

                call decode_bit_det(nJ, rep%core_space(:, j))

                if (tHPHF) then
                    hamil(i, j) = hphf_off_diag_helement(nI, nJ, &
                                                         rep%core_space(:, i), rep%core_space(:, j))
                else if (tGUGA) then
                    call calc_guga_matrix_element(rep%core_space(:, i), rep%core_space(:, j), &
                                                  excitInfo, hamil(i, j), .true., 1)
                else
                    hamil(i, j) = get_helement(nI, nJ, rep%core_space(:, i), rep%core_space(:, j))
                end if

            end do
        end do

    end subroutine calc_determin_hamil_full

    subroutine copy_core_dets_to_spawnedparts(rep)

        ! This routine will copy all the core determinants *ON THIS PROCESS
        ! ONLY* to the SpawnedParts array.

        use load_balance_calcnodes, only: DetermineDetNode
        type(core_space_t), intent(in) :: rep
        integer :: i, ncore, proc
        integer :: nI(nel)
        character(len=*), parameter :: t_r = "copy_core_dets_to_spawnedparts"

        ncore = 0
        SpawnedParts = 0_n_int

        do i = 1, rep%determ_space_size
            call decode_bit_det(nI, rep%core_space(:, i))
            proc = DetermineDetNode(nel, nI, 0)
            if (proc == iProcIndex) then
                ncore = ncore + 1
                SpawnedParts(0:NIfTot, ncore) = rep%core_space(:, i)
            end if
        end do

        if (ncore /= rep%determ_sizes(iProcIndex)) call stop_all(t_r, "The number of &
            &core determinants counted is less than was previously counted.")

    end subroutine copy_core_dets_to_spawnedparts

    subroutine return_mp1_amp_and_mp2_energy(nI, ilut, ex, tParity, amp, energy_contrib)

        ! For a given determinant (input as nI), find the amplitude of it in the MP1 wavefunction.
        ! Also return the contribution from this determinant in the MP2 energy.

        ! To use this routine, generate an excitation from the Hartree-Fock determinant using the
        ! GenExcitations3 routine. This will return nI, ex and tParity which can be input into this
        ! routine.

        use Determinants, only: GetH0Element3, GetH0Element4
        use FciMCData, only: ilutHF, HFDet, Fii
        use SystemData, only: tUEG
        use SystemData, only: tGUGA
        use guga_matrixElements, only: calcDiagMatEleGUGA_nI
        use guga_excitations, only: calc_guga_matrix_element
        use guga_data, only: ExcitationInformation_t
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: ex(2, maxExcit)
        logical, intent(in) :: tParity
        real(dp), intent(out) :: amp, energy_contrib
        integer :: ic
        HElement_t(dp) :: hel, H0tmp, denom
        type(ExcitationInformation_t) :: excitInfo

        amp = 0.0_dp
        energy_contrib = 0.0_dp

        if (ex(1, 2) == 0) then
            ic = 1
        else
            ic = 2
        end if

        if (tHPHF) then
            ! Assume since we are using HPHF that the alpha and
            ! beta orbitals of the same spatial orbital have the same
            ! fock energies, so can consider either.
            hel = hphf_off_diag_helement(HFDet, nI, iLutHF, ilut)
        else if (tGUGA) then
            ! i am not sure if the ref_stepvector thingies are set up for
            ! the ilutHF in this case..
            call calc_guga_matrix_element(ilut, ilutHF, excitInfo, hel, &
                                          .true., 2)
        else
            hel = get_helement(HFDet, nI, ic, ex, tParity)
        end if

        if (tUEG) then
            ! This will calculate the MP2 energies without having to use the fock eigenvalues.
            ! This is done via the diagonal determinant hamiltonian energies.
            H0tmp = getH0Element4(nI, HFDet)
        else if (tGUGA) then
            ! do i have a routine to calculate the diagonal and double
            ! contributions for GUGA csfs? yes!
            H0tmp = calcDiagMatEleGUGA_nI(nI)
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
            amp = hel / denom
            energy_contrib = (hel**2) / denom
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
                call decode_bit_det(nI, CurrentDets(:, i))
                call hash_search_trial(CurrentDets(:, i), nI, trial_amps, tTrial, tCon)
            else
                call bin_search_trial(CurrentDets(:, i), trial_amps, tTrial, tCon)
            end if

            ! Set the appropraite flag (if any). Unset flags which aren't
            ! appropriate, just in case.
            if (tTrial) then
                call set_flag(CurrentDets(:, i), flag_trial, .true.)
                call set_flag(CurrentDets(:, i), flag_connected, .false.)
            else if (tCon) then
                call set_flag(CurrentDets(:, i), flag_trial, .false.)
                call set_flag(CurrentDets(:, i), flag_connected, .true.)
            else
                call set_flag(CurrentDets(:, i), flag_trial, .false.)
                call set_flag(CurrentDets(:, i), flag_connected, .false.)
            end if

            ! Set the amplitude (which may be zero).
            current_trial_amps(:, i) = trial_amps
        end do

    end subroutine reinit_current_trial_amps

    subroutine end_semistoch()

        use FciMCData, only: HamTag
        use MemoryManager, only: LogMemDealloc

        character(len=*), parameter :: t_r = "end_semistoch"
        integer :: ierr
        integer :: run

        do run = 1, size(cs_replicas)
            call cs_replicas(run)%dealloc()
        end do
        deallocate(cs_replicas)

    end subroutine end_semistoch

end module semi_stoch_procs
