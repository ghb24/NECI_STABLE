#include "macros.h"

module initiator_space_procs

    use bit_rep_data, only: nifd, NIfTot
    use bit_reps, only: decode_bit_det
    use CalcData
    use constants
    use FciMCData, only: ilutHF
    use gndts_mod, only: gndts, gndts_all_sym_this_proc
    use Parallel_neci, only: iProcIndex, nProcessors, MPIArg, &
                             MPIAllGather, MPIBarrier
    use semi_stoch_gen
    use semi_stoch_procs
    use sparse_arrays
    use timing_neci
    use shared_rhash, only: initialise_shared_rht

    implicit none

    integer(n_int), allocatable :: initiator_space(:, :)

    integer(MPIArg), allocatable :: initiator_sizes(:)
    integer(MPIArg), allocatable :: initiator_displs(:)

    integer(MPIArg) :: initiator_space_size
    integer :: initiator_space_size_int

    type(shared_rhash_t) :: initiator_ht

contains

    subroutine init_initiator_space(space_in)

        use DetBitOps, only: ilut_lt, ilut_gt
        use DeterminantData, only: write_det
        use FciMCData, only: SpawnedParts, InitSpace_Init_Time
        use sort_mod, only: sort
        use SystemData, only: nel

        type(subspace_in) :: space_in

        integer :: i, ierr
        integer :: nI(nel)
        integer(MPIArg) :: mpi_temp
        character(len=*), parameter :: t_r = "init_initiator_space"

        call MPIBarrier(ierr, tTimeIn=.false.)

        call set_timer(InitSpace_Init_Time)

        write(stdout, '(/,12("="),1x,a30,1x,12("="))') "Initiator space initialisation"; call neci_flush(6)

        allocate(initiator_sizes(0:nProcessors - 1))
        allocate(initiator_displs(0:nProcessors - 1))
        initiator_sizes = 0_MPIArg
        initiator_displs = 0_MPIArg

        if (.not. (tStartCAS .or. space_in%tPops .or. space_in%tDoubles .or. space_in%tCAS .or. space_in%tRAS .or. &
                   space_in%tOptimised .or. space_in%tLowE .or. space_in%tRead .or. space_in%tMP1 .or. &
                   space_in%tFCI .or. space_in%tHeisenbergFCI .or. space_in%tHF)) then
            call stop_all("init_initiator_space", "You have not selected an initiator space to use.")
        end if

        ! Call the enumerating subroutines to create all excitations and add these states to
        ! SpawnedParts on the correct processor. As they do this, they count the size of the
        ! deterministic space (on their own processor only).
        write(stdout, '("Generating the initiator space...")'); call neci_flush(6)
        call generate_initiator_space(space_in)

        ! So that all procs store the size of the deterministic spaces on all procs.
        mpi_temp = initiator_sizes(iProcIndex)
        call MPIAllGather(mpi_temp, initiator_sizes, ierr)

        initiator_space_size = sum(initiator_sizes)
        initiator_space_size_int = int(initiator_space_size, sizeof_int)

        write(stdout, '("Total size of initiator space:",1X,i8)') initiator_space_size
        write(stdout, '("Size of initiator space on this processor:",1X,i8)') initiator_sizes(iProcIndex)
        call neci_flush(6)

        ! Calculate the indices in the full vector at which the various processors take over, relative
        ! to the first index position in the vector (i.e. the array disps in MPI routines).
        initiator_displs(0) = 0
        do i = 1, nProcessors - 1
            initiator_displs(i) = initiator_displs(i - 1) + initiator_sizes(i - 1)
        end do

        call sort(SpawnedParts(0:NIfTot, 1:initiator_sizes(iProcIndex)), ilut_lt, ilut_gt)

        ! Do a check that no states are in the initiator space twice. The list is sorted
        ! already so simply check states next to each other in the list.
        do i = 2, initiator_sizes(iProcIndex)
            if (all(SpawnedParts(0:nifd, i - 1) == SpawnedParts(0:nifd, i))) then
                call decode_bit_det(nI, SpawnedParts(:, i))
                write(stdout, '("State found twice:")')
                write(stdout, *) SpawnedParts(:, i)
                call write_det(6, nI, .true.)
                call stop_all(t_r, "The same state has been found twice in the initiator space.")
            end if
        end do

        ! Store every determinant from all processors on all processors, in initiator_space.
        call store_whole_initiator_space()
        ! Create the hash table to address the initiator determinants.
        call initialise_shared_rht(initiator_space, initiator_space_size_int, initiator_ht)

        call set_initiator_space_flags()

        SpawnedParts = 0_n_int

        ! Call MPIBarrier here so that InitSpace_Init_Time will give the
        ! initialisation time for all processors to finish.
        call MPIBarrier(ierr, tTimeIn=.false.)

        call halt_timer(InitSpace_Init_Time)

        write(stdout, '("Initialisation of initiator space complete.")')
        write(stdout, '("Total time (seconds) taken for initiator space initialisation:", f9.3, /)') &
            get_total_time(InitSpace_Init_Time)
        call neci_flush(6)

    end subroutine init_initiator_space

    subroutine generate_initiator_space(space_in)

        ! A wrapper to call the correct generating routine.

        use bit_rep_data, only: flag_initiator
        use bit_reps, only: set_flag, encode_sign
        use FciMCData, only: SpawnedParts
        use searching, only: remove_repeated_states
        use SystemData, only: tAllSymSectors

        type(subspace_in) :: space_in

        integer :: ndets_this_proc, run
        !real(dp) :: zero_sign(lenof_sign)
        character(len=*), parameter :: t_r = "generate_initiator_space"

        ndets_this_proc = 0

        ! Call the requested generating routines.
        if (space_in%tHF) call add_state_to_space(ilutHF, SpawnedParts, ndets_this_proc)
        if (space_in%tPops) call generate_space_most_populated(space_in%npops, &
                           space_in%tApproxSpace, space_in%nApproxSpace, SpawnedParts, ndets_this_proc, GLOBAL_RUN, CurrentDets, TotWalkers)
        if (space_in%tRead) call generate_space_from_file(space_in%read_filename, SpawnedParts, ndets_this_proc)
        if (space_in%tDoubles) call generate_sing_doub_determinants(SpawnedParts, ndets_this_proc, space_in%tHFConn)
        if (space_in%tCAS) call generate_cas(space_in%occ_cas, space_in%virt_cas, SpawnedParts, ndets_this_proc)
        if (space_in%tRAS) call generate_ras(space_in%ras, SpawnedParts, ndets_this_proc)
        if (space_in%tOptimised) call generate_optimised_space(space_in%opt_data, space_in%tLimitSpace, &
                                                               SpawnedParts, ndets_this_proc, space_in%max_size)
        if (space_in%tMP1) call generate_using_mp1_criterion(space_in%mp1_ndets, SpawnedParts, ndets_this_proc)
        if (space_in%tFCI) then
            if (tAllSymSectors) then
                call gndts_all_sym_this_proc(SpawnedParts, .false., ndets_this_proc)
            else
                call generate_fci_core(SpawnedParts, ndets_this_proc)
            end if
            !else if (space_in%tHeisenbergFCI) then
            !call generate_heisenberg_fci(SpawnedParts, ndets_this_proc)
        end if

        ! If two different spaces have been called then there may be some
        ! repeated states. We don't want repeats, so remove them and update
        ! ndets_this_proc accordingly.
        call remove_repeated_states(SpawnedParts, ndets_this_proc)

        !zero_sign = 0.0_dp
        !do i = 1, ndets_this_proc
        !    call encode_sign(SpawnedParts(:,i), zero_sign)

        !    if (tTruncInitiator) then
        !        do run = 1, inum_runs
        !            call set_flag(SpawnedParts(:,i), get_initiator_flag_by_run(run))
        !        end do
        !    end if
        !end do

        ! Set the initiator space size for this process.
        initiator_sizes(iProcIndex) = int(ndets_this_proc, MPIArg)

        ! If requested, remove high energy orbitals so that the space size is below some max.
        if (space_in%tLimitSpace) then
            call remove_high_energy_orbs(SpawnedParts(0:NIfTot, 1:ndets_this_proc), ndets_this_proc, &
                                         space_in%max_size, .true.)
            initiator_sizes(iProcIndex) = int(ndets_this_proc, MPIArg)
        end if

    end subroutine generate_initiator_space

    subroutine store_whole_initiator_space()

        use Parallel_neci, only: MPIAllGatherV

        integer :: ierr
        character(len=*), parameter :: t_r = "store_whole_initiator_space"

        allocate(initiator_space(0:NIfTot, initiator_space_size), stat=ierr)
        initiator_space = 0_n_int

        call MPIAllGatherV(SpawnedParts(0:NIfTot, 1:initiator_sizes(iProcIndex)), &
                           initiator_space, initiator_sizes, initiator_displs)

    end subroutine store_whole_initiator_space

    function is_in_initiator_space(ilut, nI) result(initiator_state)

        use hash, only: FindWalkerHash

        integer(n_int), intent(in) :: ilut(0:NIfTot)

        integer, intent(in) :: nI(:)
        integer(int64) :: hash_val, pos

        logical :: initiator_state

        hash_val = FindWalkerHash(nI, initiator_space_size_int)

        call initiator_ht%callback_lookup(hash_val, pos, initiator_state, loc_verify)

    contains

        function loc_verify(ind) result(match)
            integer(int64), intent(in) :: ind
            logical :: match
            match = (all(ilut(0:nifd) == initiator_space(0:nifd, ind)))
        end function loc_verify

    end function is_in_initiator_space

    subroutine set_initiator_space_flags()

        use bit_rep_data, only: flag_static_init
        use bit_reps, only: set_flag, get_initiator_flag_by_run
        use FciMCData, only: CurrentDets, TotWalkers
        integer(int64) :: i
        integer :: j
        integer :: nI(nel)
        logical :: tInitiatorDet

        tInitiatorDet = .false.

        do i = 1_int64, TotWalkers

            call decode_bit_det(nI, CurrentDets(:, i))
            tInitiatorDet = is_in_initiator_space(CurrentDets(:, i), nI)

            if (tInitiatorDet) then
                do j = 1, inum_runs
                    call set_flag(CurrentDets(:, i), get_initiator_flag_by_run(j))
                    call set_flag(CurrentDets(:, i), flag_static_init(j))
                end do
            end if

        end do

    end subroutine set_initiator_space_flags

    subroutine set_conn_init_space_flags_slow(ilut_list, list_size)

        use bit_reps, only: set_flag, get_initiator_flag
        use DetBitOps, only: CountBits, TestClosedShellDet

        ! Set initiator flags for any determinants in ilut_list which are
        ! connected to states within the initiator space

        integer(n_int), intent(inout) :: ilut_list(0:, :)
        integer, intent(in) :: list_size

        integer :: i, j, part_type, IC
        integer(n_int) :: tmp(0:NIfD)

        do i = 1, list_size
            ! If the initiator flag is already set, then we don't need to check
            ! whether or not to set it
            if (test_flag(ilut_list(:, i), get_initiator_flag(1))) cycle

            ! First the naive way: just check every single state in the initiator space...
            do j = 1, initiator_space_size
                tmp = ieor(ilut_list(0:NIfD, i), initiator_space(0:NIfD, j))
                tmp = iand(ilut_list(0:NIfD, i), tmp)
                IC = CountBits(tmp, NIfD)
                if (IC <= 2) then
                    do part_type = 1, lenof_sign
                        call set_flag(ilut_list(:, i), get_initiator_flag(part_type))
                    end do
                    exit
                end if
            end do
        end do

    end subroutine set_conn_init_space_flags_slow

end module initiator_space_procs
