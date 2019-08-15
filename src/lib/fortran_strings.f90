module fortran_strings
    use constants, only: int32, int64
    implicit none
    save
    private
    public :: str, to_lower, to_upper, operator(.in.)

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
end module fortran_strings
