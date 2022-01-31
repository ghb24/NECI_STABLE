#include "macros.h"

module test_parser_mod
    use fruit
    use constants, only: dp, n_int, int64, stdout
    use input_parser_mod, only: FileReader_t, TokenIterator_t, tokenize
    use fortran_strings, only: Token_t
    better_implicit_none
    private
    public :: test_parser_driver



contains

    subroutine test_parser_driver()
        call run_test_case(test_tokenize, "test_tokenize")
        call run_test_case(test_open_close, "test_open_close")
    end subroutine

    subroutine test_open_close()
        block
            type(FileReader_t) :: file_reader
            ! character(:), allocatable :: line
            type(TokenIterator_t) :: tokens

            file_reader = FileReader_t('H_30.FciInp')

            do while (file_reader%nextline(tokens))
                write(stdout, *)
                do while (tokens%remaining_items() >= 1)
                    write(stdout, *) tokens%get_char()
                end do
            end do

        end block
    end subroutine

    subroutine test_tokenize()
        type(Token_t), allocatable :: calculated(:), expected(:)
        type(TokenIterator_t) :: tokens

        calculated = tokenize('A B C')
        expected = [Token_t('A'), Token_t('B'), Token_t('C')]
        call assert_equals(size(calculated), 3)
        call assert_true(all(calculated == expected))


        calculated = tokenize('')
        expected = [Token_t ::]
        call assert_true(size(calculated) == 0)
        call assert_true(all(calculated == expected))

        tokens = TokenIterator_t(tokenize('A 123 B # sadfsfd'))
        call assert_equals(3, tokens%remaining_items())
        call assert_equals(tokens%get_char(), 'A')
        call assert_equals(tokens%get_int64(), 123_int64)
        call assert_equals(tokens%get_char(), 'B')
        call assert_equals(0, tokens%remaining_items())
        !
        ! expected = [Token_t('A'), Token_t('B'), Token_t('C')]
        ! call assert_true(size(calculated) == 3)
        ! call assert_true(all(calculated == expected))

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
