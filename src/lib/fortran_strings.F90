#include "macros.h"

module fortran_strings
    use constants, only: int32, int64, sp, dp
    implicit none
    save
    private
    public :: str, to_lower, to_upper, operator(.in.), split, Token_t, &
        count_char, join, to_int32, to_int64, to_realsp, to_realdp

!>  @brief
!>    Convert to Fortran string
!>
!>  @author Oskar Weser
!>
!>  @details
!>  It is a generic procedure that accepts int32 or int64.
!>
!>  @param[in] An int32 or int64.
    interface str
        module procedure int32_to_str, int64_to_str
    end interface

    interface operator(.in.)
        module procedure contains
    end interface

character(*), parameter ::  &
    UPPERCASE_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ',&
    lowercase_chars = 'abcdefghijklmnopqrstuvwxyz'


    type :: Token_t
        character(len=:), allocatable :: str
    end type

contains

    function int32_to_str(i) result(str)
        character(:), allocatable :: str
        integer(int32), intent(in) :: i
        character(range(i) + 2) :: tmp
        write(tmp, '(I0)') I
        str = trim(tmp)
    end function


    function int64_to_str(i) result(str)
        character(:), allocatable :: str
        integer(int64), intent(in) :: i
        character(range(i) + 2) :: tmp
        write(tmp, '(I0)') I
        str = trim(tmp)
    end function


    !> Changes a string to upper case
    pure function to_upper (in_str) result (string)
        character(*), intent(in) :: in_str
        character(len(in_str)) :: string
        integer :: ic, i, L

        L = len_trim(in_str)
        do i = 1, L
            ic = index(lowercase_chars, in_str(i:i))
            if (ic > 0) then
                string(i:i) = UPPERCASE_chars(ic:ic)
            else
                string(i:i) = in_str(i:i)
            end if
        end do
        string(L + 1: ) = ' '
    end function to_upper

    !> Changes a string to lower case
    pure function to_lower (in_str) result (string)
        character(*), intent(in) :: in_str
        character(len(in_str)) :: string
        integer :: ic, i, L

        L = len_trim(in_str)
        do i = 1, L
            ic = index(UPPERCASE_chars, in_str(i:i))
            if (ic > 0) then
                string(i:i) = lowercase_chars(ic:ic)
            else
                string(i:i) = in_str(i:i)
            end if
        end do
        string(L + 1: ) = ' '
    end function to_lower

    logical pure function contains(substring, string)
        character(*), intent(in) :: string, substring

        contains = index(string, substring) /= 0
    end function

    !> @brief
    !> Split string by delimiter (defaults to space).
    pure function split(expr, delimiter) result(res)
        character(*), intent(in) :: expr
        character(1), intent(in), optional :: delimiter
        type(Token_t), allocatable :: res(:)
        type(Token_t), allocatable :: tmp(:)
        character(len=1) :: delimiter_

        integer :: n, low, high

        def_default(delimiter_, delimiter, ' ')

        allocate(tmp(len(expr) / 2 + 1))
        low = 1; n = 0
        do while (low <= len(expr))
            do while (expr(low : low) == delimiter_)
                low = low + 1
                if (low > len(expr)) exit
            end do
            if (low > len(expr)) exit

            high = low
            if (high < len(expr)) then
                do while (expr(high + 1 : high + 1) /= delimiter_)
                    high = high + 1
                    if (high == len(expr)) exit
                end do
            end if
            n = n + 1
            tmp(n)%str = expr(low : high)
            low = high + 2
        end do
        res = tmp(: n)
    end function

    !> Join an array of tokens into one string
    pure function join(tokens, delimiter) result(res)
        type(Token_t), intent(in) :: tokens(:)
        character(*), intent(in) :: delimiter
        character(:), allocatable :: res
        integer :: i
        res = ''
        do i = 1, size(tokens) - 1
            res = res // tokens(i)%str // delimiter
        end do
        res = res // tokens(size(tokens))%str
    end function

    !> @brief
    !> Count the occurence of a character in a string.
    pure function count_char(str, char) result(c)
        character(len=*), intent(in) :: str
        character(len=1), intent(in) :: char
        integer :: c
        integer :: i

        c = 0
        do i = 1, len(str)
            if (str(i : i) == char) c = c + 1
        end do
    end function

    integer(int32) elemental function to_int32(str)
        character(*), intent(in) :: str
        read(unit=str, fmt=*) to_int32
    end function

    integer(int64) elemental function to_int64(str)
        character(*), intent(in) :: str
        read(unit=str, fmt=*) to_int64
    end function

    real(sp) elemental function to_realsp(str)
        character(*), intent(in) :: str
        read(unit=str, fmt=*) to_realsp
    end function

    real(dp) elemental function to_realdp(str)
        character(*), intent(in) :: str
        read(unit=str, fmt=*) to_realdp
    end function
end module fortran_strings
