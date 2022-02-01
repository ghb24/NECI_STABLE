#include "macros.h"

module test_parser_mod
    use fruit
    use constants, only: dp, n_int, int64, stdout
    use input_parser_mod, only: ManagingFileReader_t, TokenIterator_t, tokenize
    use fortran_strings, only: Token_t
    ! use util_mod, only: remove
    better_implicit_none
    private
    public :: test_parser_driver



contains

    subroutine test_parser_driver()
        call run_test_case(test_tokenize, "test_tokenize")
        call run_test_case(test_open_close, "test_open_close")
        call run_test_case(test_range, "test_range")
    end subroutine

    subroutine test_open_close()
        integer :: file_id
        character(*), parameter :: file_name = 'test.FciInp'

        open(file=file_name, newunit=file_id, action='write')
            write(file_id, '(A)') 'System read'
            write(file_id, '(A)') ''
            write(file_id, '(A)') '    # electrons           30'
            write(file_id, '(A)') '    GAS-CI GENERAL-PCHB # This is a comment'
            write(file_id, '(A)') '    GAS-SPEC LOCAL 30 +++'
            write(file_id, '(A)') '       1 4 9 +++'
            write(file_id, '(A)') '       16 25 36 # This is a comment'
        close(file_id)


        block
            type(ManagingFileReader_t) :: file_reader
            type(TokenIterator_t) :: tokens
            file_reader = ManagingFileReader_t(file_name)

            call assert_true(file_reader%nextline(tokens))
            call assert_equals(2, tokens%remaining_items())
            call assert_equals('System', tokens%get_char())
            call assert_equals('read', tokens%get_lower())
            call assert_equals(0, tokens%remaining_items())

            call assert_true(file_reader%nextline(tokens))
            call assert_equals(0, tokens%remaining_items())

            call assert_true(file_reader%nextline(tokens))
            call assert_equals(0, tokens%remaining_items())

            call assert_true(file_reader%nextline(tokens))
            call assert_equals(2, tokens%remaining_items())
            call assert_equals('GAS-CI', tokens%get_upper())
            call assert_equals('GENERAL-PCHB', tokens%get_upper())
            call assert_equals(0, tokens%remaining_items())

            call assert_true(file_reader%nextline(tokens))
            call assert_equals(9, tokens%remaining_items())
            call assert_equals('GAS-SPEC', tokens%get_upper())
            call assert_equals('LOCAL', tokens%get_upper())
            call assert_equals(30_int64, tokens%get_int64())
            block
                integer :: i
                do i = 1, 6
                    call assert_equals(i**2, tokens%get_int())
                end do
            end block
            call assert_equals(0, tokens%remaining_items())

            call assert_false(file_reader%nextline(tokens))
            ! Note that the FileReader_t closes the file automatically,
            ! when leaving scope.
        end block

        ! Delete the file
        open(file=file_name, newunit=file_id, action='read')
        close(file_id, status='delete')

    end subroutine

    subroutine test_range()
        integer, allocatable :: expected(:), calculated(:)

        expected = [1]
        calculated = get_range('1')
        call assert_equals(expected, calculated)
    end subroutine test_range

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
