#include "macros.h"

module test_parser_mod
    use fruit
    use constants, only: dp, n_int
    use input_parser_mod, only: FileReader_t, TokenIterator_t, tokenize
    use fortran_strings, only: Token_t
    better_implicit_none
    private
    public :: test_parser_driver



contains

    subroutine test_parser_driver()
        call run_test_case(test_open_close, "test_open_close")
    end subroutine

    subroutine test_open_close()
        block
            type(FileReader_t) :: file_reader
            ! character(:), allocatable :: line
            type(TokenIterator_t) :: tokens

            file_reader = FileReader_t('H_30.FciInp')

            do while (file_reader%nextline(tokens))
                write(*, *)
                do while (tokens%remaining_items() >= 1)
                    write(*, *) tokens%get_char()
                end do
            end do

        end block
    end subroutine

end module test_parser_mod

program test_parser_program

    use mpi
    use fruit
    use Parallel_neci, only: MPIInit, MPIEnd
    use test_parser_mod, only: test_parser_driver
    use util_mod, only: stop_all

    better_implicit_none
    integer :: failed_count
    block

        call MPIInit(.false.)

        call init_fruit()

        call test_parser_driver()

        call fruit_summary()
        call fruit_finalize()
        call get_failed_count(failed_count)

        if (failed_count /= 0) call stop_all('test_parser_program', 'failed_tests')

        call MPIEnd(.false.)
    end block

end program test_parser_program
