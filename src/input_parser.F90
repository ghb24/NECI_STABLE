#include "macros.h"

module input_parser_mod
    use constants, only: sp, dp, int32, int64, stderr, stdout
    use util_mod, only: stop_all, operator(.div.)
    use fortran_strings, only: Token_t, to_int, to_int32, to_int64, &
        to_realsp, to_realdp, split, to_upper, to_lower, operator(.in.)
    use growing_buffers, only: buffer_token_1D_t
    better_implicit_none
    private
    public :: FileReader_t, ManagingFileReader_t, AttachedFileReader_t, TokenIterator_t, tokenize, get_range

    integer, parameter :: max_line_length = 1028

    character(*), parameter :: delimiter = ' ', comment = '#', alt_comment = '(', concat = '+++'

    type, abstract :: FileReader_t
        private
        integer :: file_id
        integer, allocatable :: echo_lines
    contains
        private
        procedure :: raw_nextline
        procedure, public :: nextline
        procedure, public :: rewind => my_rewind
        procedure, public :: set_echo_lines
    end type

    type, extends(FileReader_t) :: ManagingFileReader_t
        private
        character(:), allocatable :: file_name
    contains
        private
        procedure, public :: close => my_close
        final :: automatic_finalize
        procedure, public :: is_open
    end type

    type, extends(FileReader_t) :: AttachedFileReader_t
        private
    contains
        private
    end type

    type :: TokenIterator_t
        private
        type(Token_t), allocatable, public :: tokens(:)
        integer :: i_curr_token = 1
    contains
        private
        procedure, public :: size => size_TokenIterator_t
        procedure, public :: remaining_items
        procedure, public :: next
        procedure, public :: reset
    end type

    interface ManagingFileReader_t
        module procedure construct_ManagingFileReader_t
    end interface

    interface AttachedFileReader_t
        module procedure construct_AttachedFileReader_t
    end interface

contains

    function construct_ManagingFileReader_t(file_name, echo_lines, err) result(res)
        character(*), intent(in) :: file_name
        integer, intent(in), optional :: echo_lines
        integer, intent(out), optional :: err
        type(ManagingFileReader_t) :: res
        res%file_name = file_name
        if (present(echo_lines)) res%echo_lines = echo_lines
        if (present(err)) then
            open(file=res%file_name, newunit=res%file_id, action='read', status='old', form='formatted', iostat=err)
        else
            open(file=res%file_name, newunit=res%file_id, action='read', status='old', form='formatted')
        end if
    end function

    function construct_AttachedFileReader_t(file_id, echo_lines) result(res)
        integer, intent(in) :: file_id
        integer, intent(in), optional :: echo_lines
        type(AttachedFileReader_t) :: res
        res%file_id = file_id
        if (present(echo_lines)) res%echo_lines = echo_lines
    end function

    impure elemental subroutine my_close(this, delete)
        class(ManagingFileReader_t), intent(inout) :: this
        logical, intent(in), optional :: delete
        deallocate(this%file_name)
        if (present(delete)) then
            if (delete) then
                close(this%file_id, status='delete')
            end if
        else
            close(this%file_id)
        end if
    end subroutine

    impure elemental subroutine automatic_finalize(this)
        type(ManagingFileReader_t), intent(inout) :: this
        if (this%is_open()) call this%close()
    end subroutine

    logical elemental function is_open(this)
        class(ManagingFileReader_t), intent(in) :: this
        is_open = allocated(this%file_name)
    end function


    logical function raw_nextline(this, line)
        !! Return if the next line can be read and return it
        !!
        !! Note that it just reads the next line and does not
        !! know about line-continuation etc.
        class(FileReader_t), intent(inout) :: this
        character(:), allocatable, intent(out) :: line
        character(*), parameter :: this_routine = 'nextline'
        character(max_line_length) :: buffer
        integer :: iread
        read(this%file_id, '(A)', iostat=iread) buffer
        raw_nextline = .false.
        if (iread > 0) then
            call stop_all(this_routine, 'Error in nextline')
        else if (is_iostat_end(iread)) then
            line = ''
        else
            raw_nextline = .true.
            line = trim(buffer)
            if (allocated(this%echo_lines)) write(this%echo_lines, '(A)') line
        end if
    end function

    logical function nextline(this, tokenized_line)
        !! Return if the next line can be read and get it tokenized.
        !!
        !! Note that it reads the next **logical** line,
        !! so if there are two lines connected by a line-continuation
        !! symbol, the two lines will be read.
        class(FileReader_t), intent(inout) :: this
        type(TokenIterator_t), intent(out) :: tokenized_line
        character(*), parameter :: this_routine = 'nextline'

        type(Token_t), allocatable :: tokens(:)
        character(:), allocatable :: line
        logical :: await_new_line


        nextline = this%raw_nextline(line)
        if (line /= '') then
            tokens = tokenize(line)
            if (size(tokens) /= 0) then
                await_new_line = tokens(size(tokens))%str == concat
                do while (await_new_line)
                    if (this%raw_nextline(line)) then
                        tokens = [tokens(: size(tokens) - 1), tokenize(line)]
                        await_new_line = tokens(size(tokens))%str == concat
                    else
                        call stop_all(this_routine, 'Open line continuation, but EOF reached.')
                    end if
                end do
            end if
            tokenized_line = TokenIterator_t(tokens)
        else
            allocate(tokenized_line%tokens(0))
        end if
    end function

    subroutine my_rewind(this)
        class(FileReader_t), intent(inout) :: this
        rewind(this%file_id)
    end subroutine

    subroutine set_echo_lines(this, echo_lines)
        class(FileReader_t), intent(inout) :: this
        integer, intent(in), optional :: echo_lines
        if (present(echo_lines)) then
            this%echo_lines = echo_lines
        else
            if (allocated(this%echo_lines)) deallocate(this%echo_lines)
        end if
    end subroutine

    function tokenize(line) result(res)
        !! Tokenize a line.
        !!
        !! If a token starts with the `comment` symbol,
        !!  this token and every following token is removed from the output.
        character(*), intent(in) :: line
        type(Token_t), allocatable :: res(:)

        integer :: i

        res = split(line, delimiter)
        do i = 1, size(res)
            if (res(i)%str(1 : len(comment)) == comment) then
                exit
            end if
            if (res(i)%str(1 : len(alt_comment)) == alt_comment) then
                write(stderr, '(A)') 'The usage of "' // alt_comment // '" as comment is deprecated.'
                write(stdout, '(A)') 'The usage of "' // alt_comment // '" as comment is deprecated.'
                exit
            end if
        end do
        res = res(: i - 1)
    end function

    integer elemental function remaining_items(this)
        class(TokenIterator_t), intent(in) :: this
        character(*), parameter :: this_routine = 'remaining_items'
        remaining_items = this%size() - this%i_curr_token + 1
        ASSERT(remaining_items >= 0)
    end function

    integer elemental function size_TokenIterator_t(this)
        class(TokenIterator_t), intent(in) :: this
        size_TokenIterator_t = size(this%tokens)
    end function

    function next(this, if_exhausted) result(res)
        class(TokenIterator_t), intent(inout) :: this
        character(*), intent(in), optional :: if_exhausted
        character(:), allocatable :: res
        character(*), parameter :: this_routine = 'next'
        integer :: i
        if (this%remaining_items() == 0) then
            if (present(if_exhausted)) then
                res = if_exhausted
            else
                write(stderr, *) 'There are no tokens remaining and the next item was requested.'
                write(stderr, *) 'The tokens are:'
                call this%reset()
                do i = 1, this%size()
                    write(stderr, *) this%tokens(i)%str
                end do
                res = ''
                call stop_all(this_routine, 'No tokens for next remaining.')
            end if
        else
            res = this%tokens(this%i_curr_token)%str
            this%i_curr_token = this%i_curr_token + 1
        end if
    end function

    elemental subroutine reset(this, k)
        class(TokenIterator_t), intent(inout) :: this
        integer, intent(in), optional :: k
        character(*), parameter :: this_routine = 'reset'
        if (present(k)) then
            if (k < 0 .and. (this%i_curr_token + k) >= 1) then
                this%i_curr_token = this%i_curr_token + k
            else
                call stop_all(this_routine, 'k has to be smaller 0 and one cannot reset past the beginning.')
            end if
        else
            this%i_curr_token = 1
        end if
    end subroutine

    pure function get_range(str_range) result(res)
        character(*), intent(in) :: str_range
        integer, allocatable :: res(:)

        character(*), parameter :: this_routine = 'get_range'
        type(Token_t), allocatable :: tokens(:)
        integer :: i

        tokens = split(str_range, '-')

        if (size(tokens) == 1) then
            res = [to_int(tokens(1)%str)]
        else if (size(tokens) == 2) then
            res = [(i, i = to_int(tokens(1)%str), to_int(tokens(2)%str))]
        else
            call stop_all(this_routine, 'Invalid input: '//str_range)
        end if
    end function

end module
