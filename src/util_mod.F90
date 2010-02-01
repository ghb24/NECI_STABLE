module util_mod
    implicit none

    ! An elemental routine to swap specified data.
    interface swap
        module procedure swap_int
        module procedure swap_doub
    end interface
    
contains

    logical function int_arr_eq (a, b, len)

        ! If two specified integer arrays are equal, return true. Otherwise
        ! return false.
        !
        ! In:  a, b - The arrays to consider
        !      len  - The maximum length to consider. Otherwise will use whole
        !             length of array

        integer, intent(in), dimension(:) :: a, b
        integer, intent(in), optional :: len
        
        integer llen, i

        ! Obtain the lengths of the arrays if a bound is not specified.
        ! Return false if mismatched sizes and not specified.
        if (present(len)) then
            llen = len
        else
            if (size(a) /= size(b)) then
                int_arr_eq = .false.
                return
            endif
            llen = size(a)
        endif

        ! Run through the arrays. Return if they differ at any point.
        do i=1,llen
            if (a(i) /= b(i)) then
                int_arr_eq = .false.
                return
            endif
        enddo

        ! If we get this far, they are equal
        int_arr_eq = .true.
    end function

    elemental subroutine swap_doub (a, b)
        
        ! Swap the doubles A, B

        real*8, intent(inout) :: a, b
        real*8 :: tmp

        tmp = a
        a = b
        b = tmp
    end subroutine

    elemental subroutine swap_int (a, b)

        ! Swap the integers A, B

        integer, intent(inout) :: a, b
        integer :: tmp

        tmp = a
        a = b
        b = tmp
    end subroutine

end module
