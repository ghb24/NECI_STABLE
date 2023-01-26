module test_shared_array_mod
    use constants, only: int64, MPIArg, stdout, stderr
    use shared_array, only: shared_array_int64_t
    use util_mod, only: stop_all
    use Parallel_neci, only: mpi_comm_intra, iProcIndex_intra, MPI_LOGICAL, MPI_LAND, &
        mpi_comm_size
#ifndef IFORT_
    use MPI_Wrapper, only: MPI_Allreduce
#endif
    use fruit, only: assert_true
    implicit none
    private
    public :: test_large_array

    type :: Interval_t
        integer(int64) :: start, end, length
    end type

contains
    !> This test should test if MPI can allocate shared arrays larger than 4Gb.
    subroutine test_large_array()
        ! The maximum number of int64 that fit in an array,
        ! that is indexed with uint32_t (32bit system). (4Gb limit)
        integer(int64) :: i
        integer(MPIArg) :: ierr, node_size
        type(shared_array_int64_t) :: shared_arr
        type(Interval_t) :: my_workarea

        write(stdout, *) 'modify the source of test_shared_array.F90 '
        write(stdout, *) 'if you want to test if MPI supports 64bit indexed shared arrays.'
        write(stdout, *) 'We had to switch the test off, otherwise the unit tests take too long.'
        call shared_arr%shared_alloc(9_int64)

        call mpi_comm_size(mpi_comm_intra, node_size, ierr)

        my_workarea = calculate_my_workarea(size(shared_arr%ptr, kind=int64), iProcIndex_intra, node_size)

        do i = my_workarea%start, my_workarea%end
            shared_arr%ptr(i) = i
        end do


        call shared_arr%sync()

        block
            logical :: correct(1), all_correct(1)
            correct(1) = .true.
            do i = 1_int64, size(shared_arr%ptr, kind=int64)
                if (shared_arr%ptr(i) /= i) then
                    correct(1) = .false.
                    exit
                end if
            end do
#ifdef USE_MPI
            call MPI_Allreduce(correct, all_correct, 1_MPIArg, MPI_LOGICAL, MPI_LAND, mpi_comm_intra, ierr)
#endif
            call assert_true(all_correct(1))
        end block
        call shared_arr%shared_dealloc()

        contains

        elemental function calculate_my_workarea(problem_size, rank, comm_size) result(res)
            integer(int64), intent(in) :: problem_size
            integer(MPIArg), intent(in) :: rank, comm_size
            type(Interval_t) :: res

            integer(int64) :: L, overhead
            L = problem_size / int(comm_size, kind=int64)
            overhead = mod(problem_size, int(comm_size, kind=int64))

            res = Interval_t(start=rank * L + min(rank, overhead) + 1_int64, &
                             end=(rank + 1_int64) * L + min(rank + 1_int64, overhead), &
                             length= merge(L + 1_int64, L, (rank + 1_int64) <= overhead))
        end function
    end subroutine test_large_array

end module


program test_shared_array

    use Parallel_neci, only: MPIInit, MPIEnd
    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count, run_test_case
    use util_mod, only: stop_all

    use test_shared_array_mod, only: test_large_array

    implicit none

    integer :: failed_count

    block

        call MPIInit(.false.)

        call init_fruit()

        call shared_array_test_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_shared_array', 'Some tests failed')

        call MPIEnd(.false.)
    end block

contains

    subroutine shared_array_test_driver()

        call run_test_case(test_large_array, "first_test")

    end subroutine shared_array_test_driver


end program test_shared_array
