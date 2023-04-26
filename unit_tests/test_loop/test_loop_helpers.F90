module test_loop_helpers
    ! use mpi
    use util_mod, only: get_free_unit
    use unit_test_helper_excitgen, only: FciDumpWriter_t, InputWriter_t, Writer_t, delete_file
    use Parallel_neci, only: MPI_Comm_Rank, mpi_comm_world
    implicit none
#include "macros.h"
#include "NECICore.h"
    private
    public :: test_loop_factory

contains

    subroutine test_loop_factory(input_writer, fcidump_writer, additional_writers)
        class(InputWriter_t), intent(in) :: input_writer
        class(FciDumpWriter_t), intent(in) :: fcidump_writer
        class(Writer_t), intent(in), optional :: additional_writers(:)
        integer :: i, myrank, err

        call mpi_comm_rank(MPI_COMM_WORLD, myrank, err)

        if (myrank == 0) then
            call input_writer%write()
            call fcidump_writer%write()
            if (present(additional_writers)) then
                do i = lbound(additional_writers, 1), ubound(additional_writers, 1)
                    call additional_writers(i)%write()
                end do
            end if
        end if

        do i = 1, 3
            call NECICore(filename_in=input_writer%filepath, &
                          & int_name=fcidump_writer%filepath, &
                          & call_as_lib=.true.)
        end do

        if (myrank == 0) then
            call delete_file(input_writer%filepath)
            call delete_file(fcidump_writer%filepath)
            if (present(additional_writers)) then
                do i = lbound(additional_writers, 1), ubound(additional_writers, 1)
                    call delete_file(additional_writers(i)%filepath)
                end do
            end if
        end if
    end subroutine test_loop_factory

end module test_loop_helpers

