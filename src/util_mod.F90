!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generic routine macro definitions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Swap the element a with the element b via a temporary variable.
#define swap_def(type,name) elemental subroutine name (a, b); \
    type, intent(inout) :: a, b; \
    type :: tmp; \
    tmp = a; \
    a = b; \
    b = tmp; \
end subroutine

! Make a comparison we can sort integer arrays by. Return true if the first
! differing items of a, b is such that a(i) > b(i).
!
! In:  a, b - The arrays to compare
! Ret:      - a > b
#define arr_gt_def(type,name) pure function name (a, b) result (bGt); \
    type, intent(in), dimension(:) :: a, b; \
    logical :: bGt; \
    integer :: i, length; \
\
    length = min(size(a), size(b)); \
\
    /* Sort by the first item first ... */ \
    do i = 1, length; \
        if (a(i) /= b(i)) exit; \
    enddo; \
\
    /* Make the comparison */ \
    if (i > length) then; \
        bGt = .false.; \
    else; \
        bGt = a(i) > b(i); \
    endif; \
end function

! Make a comparison we can sort integer arrays by. Return true if the first
! differing items of a, b is such that a(i) < b(i).
!
! In:  a, b - The arrays to compare
! Ret:      - a < b
#define arr_lt_def(type,name) pure function name (a, b) result (bLt); \
    type, intent(in), dimension(:) :: a, b; \
    logical :: bLt; \
    integer :: i, length; \
\
    length = min(size(a), size(b)); \
\
    /* Sort by the first item first ... */ \
    do i = 1, length; \
        if (a(i) /= b(i)) exit; \
    enddo; \
\
    /* Make the comparison */ \
    if (i > length) then; \
        bLt = .false.; \
    else; \
        bLt = a(i) < b(i); \
    endif; \
end function



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module starts here
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module util_mod
    implicit none
    private

    public :: swap, arr_lt, arr_gt, operator(.arrlt.), operator(.arrgt.)
    public factrl, choose, int_fmt
    

    ! Note that in these interfaces we do NOT include a version for reals, as
    ! we are using -r8 in the compile scripts --> real == real*8.
    ! If we choose to change this, then simply add a new version.

    ! An elemental routine to swap specified data.
    interface swap
        module procedure swap_int
        module procedure swap_doub
        module procedure swap_cplx
    end interface

    ! Provide operators for array comparison. Would like to use arr_gt, arr_lt
    ! but fortran requires operators to only consist of letters
    ! Not defined for complex data types, as ordering not defined.
    interface operator (.arrgt.)
        module procedure arr_gt_int
        module procedure arr_gt_doub
    end interface

    interface operator (.arrlt.)
        module procedure arr_lt_int
        module procedure arr_lt_doub
    end interface

    ! And a function interface to the same
    interface arr_gt
        module procedure arr_gt_int
        module procedure arr_gt_doub
    end interface

    interface arr_lt
        module procedure arr_lt_int
        module procedure arr_lt_doub
    end interface

contains

    ! Swap routines for different variable types.
    swap_def(integer, swap_int)
    swap_def(real*8, swap_doub)
    swap_def(complex*16, swap_cplx)

    ! Comparison routines for arrays to allow sorting
    arr_gt_def(integer, arr_gt_int)
    arr_gt_def(real*8, arr_gt_doub)

    arr_lt_def(integer, arr_lt_int)
    arr_lt_def(real*8, arr_lt_doub)

    elemental real*8 function factrl (n)

        ! Return the factorial on n, i.e. n!
        ! This is not done in the most efficient way possible (i.e. use with
        ! care if N is large, or if called many times!).
        ! If a more efficient procedure is required, refer to:
        ! http://www.luschny.de/math/factorial/FastFactorialFunctions.htm.

        integer, intent(in) :: n
        integer :: i

        factrl = 1
        do i = 2, n
            factrl = factrl * i
        enddo
    end function factrl

    elemental real*8 function choose (n, r)
        
        ! Return the binomail coefficient nCr

        integer, intent(in) :: n, r
        integer :: i, k

        if (r > n) then
            choose = 0
        else
            ! Always use the smaller possibility
            if (r > (n / 2)) then
                k = n - r
            else
                k = r
            endif

            choose = 1
            do i = 0, k-1
                choose = (choose * (n - i)) / (i + 1)
            enddo
        endif
    end function choose


    elemental function int_fmt(i, padding) result(fmt1)

        ! In:
        !    i: an integer
        !    padding (optional): amount of padding to add to format statement.
        !        The default amount is 2.  The padding is used to include the
        !        sign if i is negative.
        ! Returns:
        !    fmt1: a format statement for an integer field which will hold
        !        i perfectly plus an amount of padding.

        ! This does take i/o formatting to a slightly OCD level addmittedly...

        character(2) :: fmt1
        integer, intent(in) :: i
        integer, intent(in), optional :: padding
        integer :: p
        real :: r

        if (present(padding)) then
            p = padding
        else
            p  = 2
        end if

        if (i == 0 .or. i==1) then
            r = 1.0
        else
            r = log10(real(abs(i)+1))
        end if

        if (r < 10) then
            write (fmt1,'("i",i1)') ceiling(r+p)
        else if (r < 100) then
            write (fmt1,'("i",i2)') ceiling(r+p)
        else
            ! By this point we'll have hit integer overflow anyway...
            write (fmt1,'("i",i3)') ceiling(r+p)
        end if

    end function int_fmt

    pure function binary_search (arr, val, data_size, num_items) result(pos)

        integer, intent(in) :: data_size, num_items
        integer, intent(in) :: arr(data_size, num_items)
        integer, intent(in) :: val(data_size)
        integer :: pos

        integer :: hi, lo

        ! The search range
        lo = 1
        hi = num_items

        ! Narrow the search range down in steps.
        do while (hi /= lo)
            pos = int(real(hi + lo) / 2)

            if (all(arr(:,pos) == val)) then
                exit
            else if (arr_gt(val, arr(:,pos))) then
                ! val is "greater" than arr(:,pos).
                ! The lowest position val can take is hence pos + 1 (i.e. if
                ! val is greater than pos by smaller than pos + 1).
                lo = pos + 1
            else
                ! arr(:,pos) is "greater" than val.
                ! The highest position val can take is hence pos (i.e. if val is
                ! smaller than pos but greater than pos - 1).  This is why
                ! we differ slightly from a standard binary search (where lo
                ! is set to be pos+1 and hi to be pos-1 accordingly), as
                ! a standard binary search assumes that the element you are
                ! searching for actually appears in the array being
                ! searched...
                hi = pos
            endif
        enddo

        ! If we have narrowed down to one position, and it is not the item, 
        ! then return -pos to indicate that the item is not present, but that
        ! this is the location it should be in.
        if (hi == lo) then
            if (all(arr(:,hi) == val)) then
                pos = hi
            else if (arr_gt(val, arr(:,hi))) then
                pos = -hi - 1
            else
                pos = -hi
            endif
        endif

    end function

end module
