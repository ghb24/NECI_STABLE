#include "macros.h"

module input_parser_mod
    use constants, only: sp, dp, int32, int64, stderr, stdout
    use util_mod, only: stop_all, operator(.div.)
    use fortran_strings, only: Token_t, to_int, to_int32, to_int64, &
        to_realsp, to_realdp, split, to_upper, to_lower
    use growing_buffers, only: buffer_token_1D_t
    better_implicit_none
    private
    public :: FileReader_t, TokenIterator_t, tokenize

    integer, parameter :: max_line_length = 1028

    character(*), parameter :: delimiter = ' ', comment = '#', alt_comment = '(', concat = '+++'

    type :: FileReader_t
        private
        character(:), allocatable :: file_name
        integer :: file_id
        integer :: current_line = 0

    contains
        private
        procedure, public :: close => my_close
        final :: automatic_finalize

        procedure, public :: is_open
        procedure, public :: raw_nextline
        procedure, public :: nextline
    end type

    type :: TokenIterator_t
        private
        type(Token_t), allocatable, public :: tokens(:)
        integer :: i_curr_token = 1
    contains
        private
        procedure, public :: remaining_items
        procedure, public :: get_char
        procedure, public :: get_upper
        procedure, public :: get_lower
        procedure, public :: get_realsp
        procedure, public :: get_realdp
        procedure, public :: get_int
        procedure, public :: get_int32
        procedure, public :: get_int64
        procedure, public :: reset
    end type

    interface FileReader_t
        module procedure construct_File_t
    end interface

contains

    function construct_File_t(file_name, err) result(res)
        character(*), intent(in) :: file_name
        integer, intent(out), optional :: err
        type(FileReader_t) :: res
        res%file_name = file_name
        if (present(err)) then
            open(file=res%file_name, newunit=res%file_id, action='read', form='formatted', iostat=err)
        else
            open(file=res%file_name, newunit=res%file_id, action='read', form='formatted')
        end if
    end function

    impure elemental subroutine my_close(this)
        class(FileReader_t), intent(inout) :: this
        deallocate(this%file_name)
        close(this%file_id)
    end subroutine

    impure elemental subroutine automatic_finalize(this)
        type(FileReader_t), intent(inout) :: this
        call this%close()
    end subroutine

    logical elemental function is_open(this)
        class(FileReader_t), intent(in) :: this
        is_open = allocated(this%file_name)
    end function


    integer elemental function get_current_line(this)
        class(FileReader_t), intent(in) :: this
        get_current_line = this%current_line
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
        ASSERT(this%is_open())
        read(this%file_id, '(A)', iostat=iread) buffer
        raw_nextline = .false.
        if (iread > 0) then
            call stop_all(this_routine, 'Error in nextline')
        else if (is_iostat_end(iread)) then
            line = ''
        else
            raw_nextline = .true.
            line = trim(buffer)
            this%current_line = this%current_line + 1
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
        remaining_items = size(this%tokens) - this%i_curr_token + 1
        ASSERT(remaining_items >= 0)
    end function

    function get_char(this) result(res)
        class(TokenIterator_t), intent(inout) :: this
        character(:), allocatable :: res
        character(*), parameter :: this_routine = 'get_char'
        ASSERT(this%remaining_items() >= 1)
        res = this%tokens(this%i_curr_token)%str
        this%i_curr_token = this%i_curr_token + 1
    end function

    function get_lower(this) result(res)
        class(TokenIterator_t), intent(inout) :: this
        character(:), allocatable :: res
        res = to_lower(this%get_char())
    end function

    function get_upper(this) result(res)
        class(TokenIterator_t), intent(inout) :: this
        character(:), allocatable :: res
        res = to_upper(this%get_char())
    end function

    integer impure elemental function get_int(this)
        class(TokenIterator_t), intent(inout) :: this
        get_int = to_int(this%get_char())
    end function

    integer(int32) impure elemental function get_int32(this)
        class(TokenIterator_t), intent(inout) :: this
        get_int32 = to_int32(this%get_char())
    end function

    integer(int64) impure elemental function get_int64(this)
        class(TokenIterator_t), intent(inout) :: this
        get_int64 = to_int64(this%get_char())
    end function

    real(sp) impure elemental function get_realsp(this)
        class(TokenIterator_t), intent(inout) :: this
        get_realsp = to_realsp(this%get_char())
    end function

    real(dp) impure elemental function get_realdp(this)
        class(TokenIterator_t), intent(inout) :: this
        get_realdp = to_realdp(this%get_char())
    end function

    elemental subroutine reset(this)
        class(TokenIterator_t), intent(inout) :: this
        this%i_curr_token = 1
    end subroutine

end module
