#include "macros.h"

module test_parser_mod
    use fruit
    use constants, only: dp, n_int, int64, stdout
    use input_parser_mod, only: ManagingFileReader_t, TokenIterator_t, tokenize, get_range
    use fortran_strings, only: Token_t, to_lower, to_upper, to_int, to_int64
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
            write(file_id, '(A)') '    GAS-SPEC LOCAL 30 \'
            write(file_id, '(A)') '       1 4 9 \'
            write(file_id, '(A)') '       16 25 36 # This is a comment'
        close(file_id)


        block
            type(ManagingFileReader_t) :: file_reader
            type(TokenIterator_t) :: tokens
            file_reader = ManagingFileReader_t(file_name)

            call assert_equals(file_name, file_reader%get_file_name())

            call assert_true(file_reader%nextline(tokens))
            call assert_equals(2, tokens%remaining_items())
            call assert_equals('System', tokens%next())
            call assert_equals('read', to_lower(tokens%next()))
            call assert_equals(1, file_reader%get_current_line())
            call assert_equals(0, tokens%remaining_items())

            call assert_true(file_reader%nextline(tokens))
            call assert_equals(0, tokens%remaining_items())

            call assert_true(file_reader%nextline(tokens))
            call assert_equals(0, tokens%remaining_items())

            call assert_true(file_reader%nextline(tokens))
            call assert_equals(2, tokens%remaining_items())
            call assert_equals('GAS-CI', to_upper(tokens%next()))
            call assert_equals('GENERAL-PCHB', to_upper(tokens%next()))
            call assert_equals(0, tokens%remaining_items())

            call assert_true(file_reader%nextline(tokens))
            call assert_equals(7, file_reader%get_current_line())
            call assert_equals(9, tokens%remaining_items())
            call assert_equals('GAS-SPEC', to_upper(tokens%next()))
            call assert_equals('LOCAL', to_upper(tokens%next()))
            call assert_equals(30_int64, to_int64(tokens%next()))
            block
                integer :: i
                do i = 1, 6
                    call assert_equals(i**2, to_int(tokens%next()))
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
        call assert_equals(size(expected), size(calculated))
        call assert_equals(expected, calculated, size(expected))

        expected = [1]
        calculated = get_range('1-1')
        call assert_equals(size(expected), size(calculated))
        call assert_equals(expected, calculated, size(calculated))

        expected = [1, 2, 3, 4]
        calculated = get_range('1-4')
        call assert_equals(size(expected), size(calculated))
        call assert_equals(expected, calculated, size(expected))

        expected = [integer :: ]
        calculated = get_range('4-1')
        call assert_equals(size(expected), size(calculated))
        call assert_equals(expected, calculated, size(calculated))
    end subroutine test_range

    subroutine test_tokenize()
        type(Token_t), allocatable :: calculated(:), expected(:)
        type(TokenIterator_t) :: tokens

        calculated = tokenize('A B C')
        expected = [Token_t('A'), Token_t('B'), Token_t('C')]
        call assert_equals(size(calculated), 3)
        call assert_true(all(calculated == expected))


        calculated = tokenize('')
        deallocate(expected); allocate(expected(0))
        call assert_true(size(calculated) == 0)
        call assert_true(all(calculated == expected))

        tokens = TokenIterator_t(tokenize('A 123 B # sadfsfd'))
        call assert_equals(3, tokens%remaining_items())
        call assert_equals(tokens%next(), 'A')
        call assert_equals(to_int64(tokens%next()), 123_int64)
        call assert_equals(tokens%next(), 'B')
        call assert_equals(0, tokens%remaining_items())
        call assert_equals(tokens%next(if_exhausted='asdf'), 'asdf')

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
