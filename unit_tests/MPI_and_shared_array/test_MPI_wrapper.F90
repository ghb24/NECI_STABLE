#include "macros.h"
module test_MPI_wrapper_mod
    use constants, only: int64, MPIArg, stdout, stderr
    use Parallel_neci, only: mpi_comm_intra, iProcIndex_intra, &
        mpi_comm_size, Node, MPIBcast
    use fruit, only: assert_true, run_test_case
    better_implicit_none
    private
    public :: test_driver

contains

    subroutine test_driver()
        call run_test_case(test_MPIBcast, "test_MPIBcast")
    end subroutine

    subroutine test_MPIBcast()
        block
            logical :: communicated_even
            integer :: rank, ierr, node_size
            call mpi_comm_size(mpi_comm_intra, node_size, ierr)
            associate(i_am_even => mod(iProcIndex_intra, 2) == 0)

                do rank = 0, node_size - 1
                    communicated_even = i_am_even
                    call MPIBCast(communicated_even, iProcIndex_intra == rank, Node)
                    call assert_true(communicated_even .eqv. mod(rank, 2) == 0)
                end do
            end associate
        end block
    end subroutine



end module


program test_MPI_wrapper

    use Parallel_neci, only: MPIInit, MPIEnd
    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count
    use util_mod, only: stop_all

    use test_MPI_wrapper_mod, only: test_driver

    implicit none

    integer :: failed_count

    block

        call MPIInit(.false.)

        call init_fruit()

        call test_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_shared_array', 'Some tests failed')

        call MPIEnd(.false.)
    end block


end program test_MPI_wrapper
