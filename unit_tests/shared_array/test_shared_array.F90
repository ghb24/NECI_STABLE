module test_shared_array_mod
    use mpi
    use constants, only: int64, MPIArg
    use shared_array, only: shared_array_int64_t
    use Parallel_neci, only: mpi_comm_intra, iProcIndex_intra
    use fruit
    implicit none
    private
    public :: first_test

contains
    subroutine first_test()
        integer(int64), parameter :: val = 1_int64
        ! The maximum number of real64 floats that fit in an array, that is
        ! indexed with int32
        integer(int64), parameter :: max_size = 2_int64**32 / sizeof(val)
        integer(int64) :: i
        integer(MPIArg) :: ierr
        type(shared_array_int64_t) :: shared_arr

        call shared_arr%shared_alloc(max_size + 100)

        if (iProcIndex_intra == 0) then
            do i = 1_int64, size(shared_arr%ptr, kind=int64)
                shared_arr%ptr(i) = i
            end do
        end if
        call shared_arr%sync()

        block
            logical :: correct, all_correct
            correct = .true.
            do i = 1_int64, size(shared_arr%ptr, kind=int64)
                if (shared_arr%ptr(i) /= 2) then
                    correct = .false.
                    exit
                end if
            end do
            call MPI_Allreduce(correct, all_correct, 1, MPI_LOGICAL, MPI_LAND, mpi_comm_intra, ierr)
            call assert_true(all_correct)
        end block
        call shared_arr%shared_dealloc()
    end subroutine first_test

end module


program test_shared_array

    use Parallel_neci, only: MPIInit, MPIEnd
    use util_mod, only: stop_all
    use fruit

    use test_shared_array_mod
    ! use shared_array

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

        call run_test_case(first_test, "first_test")

    end subroutine shared_array_test_driver


end program test_shared_array

