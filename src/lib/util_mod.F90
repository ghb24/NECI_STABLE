module util_mod
    use util_mod_comparisons
    use util_mod_cpts
    use constants, only: dp
    implicit none

    ! sds: It would be nice to use a proper private/public interface here,
    !      BUT PGI throws a wobbly on using the public definition on
    !      a new declared operator. --> "Empty Operator" errors!
    !      to fix when compilers work!
!    private

!    public :: swap, arr_lt, arr_gt, operator(.arrlt.), operator(.arrgt.)
!    public :: factrl, choose, int_fmt, binary_search
!    public :: append_ext, get_unique_filename, get_nan, isnan
    
contains

!--- Array utilities ---

      SUBROUTINE NECI_ICOPY(N,A,IA,B,IB)
         ! Copy elements from integer array A to B.
         ! Simple version of BLAS routine ICOPY, which isn't always implemented
         ! in BLAS.
         ! Fortran 90 array features allow this to be done in one line of
         ! standard fortran, so this is just for legacy purposes.
         ! In:
         !    N: number of elements in A.
         !    A: vector to be copied.
         !    IA: increment between elements to be copied in A.  
         !        IA=1 for continuous data blocks.
         !    IB: increment between elements to be copied to in B.  
         !        IB=1 for continuous data blocks.
         ! Out:
         !    B: result vector.
         IMPLICIT NONE
!        Arguments
         INTEGER, INTENT(IN) :: N,IA,IB
         INTEGER, INTENT(IN) :: A(IA*N)
         INTEGER, INTENT(OUT) :: B(IB*N)
!        Variables
         INTEGER I,IAX,IBX
     
         DO I=1,N
           IAX=(I-1)*IA + 1
           IBX=(I-1)*IB + 1
           B(IBX) = A(IAX)
         ENDDO
     
         RETURN
      END SUBROUTINE NECI_ICOPY

!--- Numerical utilities ---

    ! If all of the compilers supported ieee_arithmetic
    ! --> could use ieee_value(1.0_dp, ieee_quiet_nan)
    real(dp) function get_nan ()
        real(dp) :: a, b
        a = 1
        b = 1
        get_nan = log (a-2*b)
    end function

    ! If all of the compilers supported ieee_arithmetic
    ! --> could use ieee_is_nan (r)
    elemental logical function isnan (r)
        real(dp), intent(in) :: r

        if ( (r == 0) .and. (r * 1 == 1) ) then
            isnan = .true.
        else
            isnan = .false.
        endif
    end function

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

!--- Comparison of subarrays ---
    
    logical pure function det_int_arr_gt (a, b, len)
        use constants, only: n_int

        ! Make a comparison we can sort determinant integer arrays by. Return true if the
        ! first differing integer of a, b is such that a(i) > b(i).
        !
        ! In:  a, b - The arrays to compare
        !      len  - An optional argument to specify the size to consider.
        !             If not provided, then min(size(a), size(b)) is used.
        ! Ret:      - If a > b
        !NOTE: These will sort by the bit-string integer length, n_int.
        !Therefore, these may be 32 or 64 bit integers and should
        !only be used as such.
    
        integer(kind=n_int), intent(in), dimension(:) :: a, b
        integer, intent(in), optional :: len

        integer llen, i

        if (present(len)) then
            llen = len
        else
            llen = min(size(a), size(b))
        endif

        ! Sort by the first integer first ...
        i = 1
        do i = 1, llen
            if (a(i) /= b(i)) exit
        enddo

        ! Make the comparison
        if (i > llen) then
            det_int_arr_gt = .false.
        else
            if (a(i) > b(i)) then
                det_int_arr_gt = .true.
            else
                det_int_arr_gt = .false.
            endif
        endif
    end function det_int_arr_gt
        

    logical pure function det_int_arr_eq (a, b, len)
        use constants, only: n_int

        ! If two specified integer arrays are equal, return true. Otherwise
        ! return false.
        !
        ! In:  a, b - The arrays to consider
        !      len  - The maximum length to consider. Otherwise will use whole
        !             length of array
        !NOTE: These will sort by the bit-string integer length, n_int.
        !Therefore, these may be 32 or 64 bit integers and should
        !only be used as such.

        integer(kind=n_int), intent(in), dimension(:) :: a, b
        integer, intent(in), optional :: len
        
        integer llen, i

        ! Obtain the lengths of the arrays if a bound is not specified.
        ! Return false if mismatched sizes and not specified.
        if (present(len)) then
            llen = len
        else
            if (size(a) /= size(b)) then
                det_int_arr_eq = .false.
                return
            endif
            llen = size(a)
        endif

        ! Run through the arrays. Return if they differ at any point.
        do i=1,llen
            if (a(i) /= b(i)) then
                det_int_arr_eq = .false.
                return
            endif
        enddo

        ! If we get this far, they are equal
        det_int_arr_eq = .true.
    end function det_int_arr_eq

!--- Output utilties ---

    elemental function int_fmt(i, padding) result(fmt1)

        ! In:
        !    i: an integer
        !    padding (optional): amount of padding to add to format statement.
        !        The default amount is 2.  The padding is used to include the
        !        sign if i is negative.
        ! Returns:
        !    fmt1: a format statement for an integer field which will hold
        !        i perfectly plus an amount of padding.

        ! This does take i/o formatting to a slightly OCD level, admittedly...

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

!--- Searching ---

    ! NOTE: This can only be used for binary searching determinant bit 
    !       strings now. We can template it if it wants to be more general 
    !       in the future if needed.
    pure function binary_search (arr, val, data_size, num_items) result(pos)
        use constants, only: n_int

        integer, intent(in) :: data_size, num_items
        integer(kind=n_int), intent(in) :: arr(data_size, num_items)
        integer(kind=n_int), intent(in) :: val(data_size)
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

!--- File utilities ---

    integer function record_length(bytes)
       ! Some compilers use record lengths in units of bytes.
       ! Some compilers use record lengths in units of words.
       ! This is an utter *pain* for reading unformatted files,
       ! where you must specify the record length.
       !
       ! In:
       !    bytes: number of bytes in record type of interest (should
       !    be a multiple of 4).
       !
       ! Returns:
       !    record_length: size of record in units of the compiler's
       !    choice.
       integer, intent(in) :: bytes
       inquire(iolength=record_length) bytes
       record_length = (bytes/4)*record_length 
    end function record_length

    subroutine append_ext(stem, n, s)

        ! Returns stem.n in s.

        character(*), intent(in) :: stem
        integer, intent(in) :: n
        character(*), intent(out) :: s
        character(10) :: ext

        write (ext,'('//int_fmt(n,0)//')') n
        s = stem//'.'//ext

    end subroutine append_ext

    subroutine get_unique_filename(stem, tincrement, tnext, istart, filename)

        ! Find a filename which is either the "newest" or the next to be used.
        ! The filename is assumed to be stem.x, where x is an integer.

        ! In:
        !    stem: stem of the filename.
        !    tincrement: the filename is given as stem.x if true, otherwise the
        !        filename is simply set to be equal to stem.
        !    tnext: the next unused filename is found if true, else the
        !        filename is set to be stem.x where stem.x exists and stem.x+1
        !        doesn't and x is greater than istart or unless the file
        !        stem exists, then the filename is set to be stem (with no
        !        extension).
        !    istart: the integer of the first x value to check.
        !        If istart is negative, then the filename is set to be stem.x,
        !        where x = |istart+1|.  This overrides everything else.
        ! Out:
        !    filename.

        character(*), intent(in) :: stem
        logical, intent(in) :: tincrement, tnext
        integer, intent(in) :: istart
        character(*), intent(out) :: filename

        integer :: i
        logical :: exists

        if (tincrement) then
            i = istart
            exists = .true.
            do while (exists)
                call append_ext(stem, i, filename)
                inquire(file=filename,exist=exists)
                i = i + 1
            end do
            if (.not.tnext) then
                ! actually want the last file which existed.
                ! this will return stem.istart if stem.istart doesn't exist.
                i = max(istart,i - 2)
                call append_ext(stem, i, filename)
            end if
        else
            filename = stem
        end if

        if (.not.tnext) then
            inquire(file=filename,exist=exists)
            if (.not.exists) then
                inquire(file=stem,exist=exists)
                if (exists) filename = stem
            end if
        end if

        if (istart < 0) then
            call append_ext(stem, abs(i+1), filename)
        end if

    end subroutine get_unique_filename

    function get_free_unit() result(free_unit)

        ! Returns:
        !    The first free file unit above 10 and less than or equal to
        !    the paramater max_unit (currently set to 200).

        integer, parameter :: max_unit = 100
        integer :: free_unit
        integer :: i
        logical :: t_open, t_exist

        do i = 10, max_unit
            inquire(unit=i, opened=t_open, exist=t_exist)
            if (.not.t_open .and. t_exist) then
                free_unit = i
                exit
            end if
        end do
        if (i == max_unit+1) call stop_all('get_free_unit','Cannot find a free unit below max_unit.')

    end function get_free_unit

end module
