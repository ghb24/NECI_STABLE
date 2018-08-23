module util_mod
    use util_mod_comparisons
    use util_mod_numerical
    use util_mod_byte_size
    use util_mod_cpts
    use fmt_utils
    use dSFMT_interface, only: genrand_real2_dSFMT
    use constants
    use iso_c_hack

    ! We want to use the builtin etime intrinsic with ifort to 
    ! work around some broken behaviour.
#ifdef __IFORT
    use ifport, only: etime
#endif
    implicit none

    interface
        pure function strlen_wrap (str) result(len) bind(c)
            use iso_c_hack
            implicit none
            character(c_char), intent(in) :: str(*)
            integer(c_int) :: len
        end function
        pure function erf_local (x) result(e) bind(c, name='erf')
            use iso_c_hack
            implicit none
            real(c_double), intent(in) :: x
            real(c_double) :: e
        end function
        pure function erfc_local (x) result(ec) bind(c, name='erfc')
            use iso_c_hack
            implicit none
            real(c_double), intent(in) :: x
            real(c_double) :: ec
        end function
    end interface

    interface abs_sign
        module procedure abs_int4_sign
        module procedure abs_int8_sign
        module procedure abs_real_sign
    end interface

    interface abs_l1
        module procedure abs_l1_dp
        module procedure abs_l1_sp
        module procedure abs_l1_cdp
        module procedure abs_l1_csp
    end interface


    ! sds: It would be nice to use a proper private/public interface here,
    !      BUT PGI throws a wobbly on using the public definition on
    !      a new declared operator. --> "Empty Operator" errors!
    !      to fix when compilers work!
!    private

!    public :: swap, arr_lt, arr_gt, operator(.arrlt.), operator(.arrgt.)
!    public :: factrl, choose, int_fmt, binary_search
!    public :: append_ext, get_unique_filename, get_nan, isnan_neci

contains

    function stochastic_round (r) result(i)

        ! Stochastically round the supplied real value to an integer. This is
        ! the primary method of introducing the monte-carlo nature of spawning
        ! or death into the algorithm.
        ! --> Probably nicer to use a centralised implementation than a bunch
        !     of hacked-in ones all over the place...
        !
        ! Unfortunately, we cannot make this pure, as we would need to have
        ! a mutable variable in genrand_real2_dSFMT...

        real(dp), intent(in) :: r
        integer :: i
        real(dp) :: res

        i = int(r)
        res = r - real(i, dp)
        
        if (abs(res) >= 1.0e-12_dp) then
            if (abs(res) > genrand_real2_dSFMT()) &
                i = i + nint(sign(1.0_dp, r))
        end if

    end function

    function stochastic_round_r (num, r) result(i)

        ! Perform the stochastic rounding of the above function where the 
        ! random number is already specified.

        real(dp), intent(in) :: num, r
        integer :: i
        real(dp) :: res

        i = int(num)
        res = num - real(i, dp)

        if (abs(res) >= 1.0e-12_dp) then
            if (abs(res) > r) &
                i = i + nint(sign(1.0_dp, num))
        end if

    end function

    subroutine print_cstr (str) bind(c, name='print_cstr')

        ! Write a string outputted by calling fort_printf in C.
        ! --> Ensure that it is redirected to the same place as the normal
        !     STDOUT within fortran.

        character(c_char), intent(in) :: str(*)
        integer :: l

        l = strlen_wrap(str)
        call print_cstr_local (str, l)

    end subroutine

    subroutine print_cstr_local (str, l)

        character(c_char), intent(in) :: str(*)
        integer, intent(in) :: l
        character(len=l) :: tmp_s

        tmp_s = transfer(str(1:l), tmp_s)
        write(6, '(a)', advance='no') tmp_s

    end subroutine

    ! routine to calculation the absolute magnitude of a complex integer 
    ! variable (to nearest integer)
    pure real(dp) function abs_int4_sign(sgn)
        integer(int32), intent(in) :: sgn(lenof_sign/inum_runs)

#ifdef __CMPLX
            abs_int4_sign=real(int(sqrt(real(sgn(1),dp)**2+real(sgn(2),dp)**2)),dp)
            ! The integerisation here is an approximation, but one that is 
            ! used in the integer algorithm, so is retained in this real 
            ! version of the algorithm
#else
            abs_int4_sign=abs(sgn(1))
#endif
    end function abs_int4_sign

!routine to calculation the absolute magnitude of a complex integer(int64) variable (to nearest integer)
    pure integer(kind=int64) function abs_int8_sign(wsign)
        integer(kind=int64), dimension(lenof_sign/inum_runs), intent(in) :: wsign

#ifdef __CMPLX
            abs_int8_sign=nint(sqrt(real(wsign(1),dp)**2+real(wsign(2),dp)**2),int64)
#else
            abs_int8_sign=abs(wsign(1))
#endif
    end function abs_int8_sign

    pure real(dp) function abs_real_sign (sgn)
        real(dp), intent(in) :: sgn(lenof_sign/inum_runs)
#ifdef __CMPLX
            abs_real_sign = real(nint(sqrt(sum(sgn ** 2))), dp)
#else
            abs_real_sign = abs(sgn(1))
#endif
    end function

! --------- L1 norm function
! These return the absolute L1 norm of the specified value
!
! --> for complex numbers this is not sqrt(r**2 + i**2), but is the sum
!     of the absolute values of the terms
!----------------------------

    pure function abs_l1_sp (val) result (ret)

        real(sp), intent(in) :: val
        real(sp) :: ret

        ret = abs(val)

    end function

    pure function abs_l1_dp (val) result (ret)

        real(dp), intent(in) :: val
        real(dp) :: ret

        ret = abs(val)

    end function

    pure function abs_l1_csp (val) result (ret)

        complex(sp), intent(in) :: val
        real(sp) :: ret

        ret = abs(real(val, sp)) + abs(aimag(val))

    end function

    pure function abs_l1_cdp (val) result (ret)

        complex(dp), intent(in) :: val
        real(dp) :: ret

        ret = abs(real(val, dp)) + abs(aimag(val))

    end function

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
    elemental logical function isnan_neci (r)
        real(dp), intent(in) :: r

#ifdef __GFORTRAN__
        isnan_neci = isnan(r)
#else
        if ( (r == 0) .and. (r * 1 == 1) ) then
            isnan_neci = .true.
        else
            isnan_neci = .false.
        endif
#endif
    end function

    elemental real(dp) function factrl (n)

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

    elemental real(dp) function choose (n, r)

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

!--- Searching ---

    ! NOTE: This can only be used for binary searching determinant bit
    !       strings now. We can template it if it wants to be more general
    !       in the future if needed.
    function binary_search (arr, val, cf_len) &
                                 result(pos)
        use constants, only: n_int

        integer(kind=n_int), intent(in) :: arr(:,:)
        integer(kind=n_int), intent(in) :: val(:)
        integer, intent(in), optional :: cf_len
        integer :: data_lo, data_hi, val_lo, val_hi
        integer :: pos

        integer :: hi, lo

        ! The search range
        lo = lbound(arr,2)
        hi = ubound(arr,2)

        ! Account for poor usage (i.e. array len == 0)
        if (hi < lo) then
            pos = -lo
            return
        endif

        ! Have we specified how much to look at?
        data_lo = lbound(arr, 1)
        val_lo = lbound(val, 1)
        if (present(cf_len)) then
            data_hi = data_lo + cf_len - 1
            val_hi = val_lo + cf_len - 1
        else
            data_hi = ubound(arr, 1)
            val_hi = ubound(val, 1)
        endif

        ! Narrow the search range down in steps.
        do while (hi /= lo)
            pos = int(real(hi + lo,sp) / 2_sp)

            if (all(arr(data_lo:data_hi,pos) == val(val_lo:val_hi))) then
                exit
            else if (arr_gt(val(val_lo:val_hi), arr(data_lo:data_hi,pos))) then
                ! val is "greater" than arr(:len,pos).
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
            if (all(arr(data_lo:data_hi,hi) == val(val_lo:val_hi))) then
                pos = hi
            else if (arr_gt(val(val_lo:val_hi), arr(data_lo:data_hi,hi))) then
                pos = -hi - 1
            else
                pos = -hi
            endif
        endif

    end function binary_search

    function binary_search_int (arr, val) result(pos) 
        ! W.D.: also write a binary search for "normal" lists of ints
        integer, intent(in) :: arr(:)
        integer, intent(in) :: val
        integer :: pos 

        integer :: hi, lo

        lo = lbound(arr,1)
        hi = ubound(arr,1)

        if (hi < lo) then 
            pos  = -lo
            return
        end if

        do while (hi /= lo)
            pos = int(real(hi + lo, sp) / 2.0_dp)

            if (arr(pos) == val) then 
                exit 
            else if (val > arr(pos)) then 
                lo = pos + 1
            else 
                hi = pos
            end if
        end do

        if (hi == lo) then 
            if (arr(hi) == val) then 
                pos = hi 
            else if (val > arr(hi)) then 
                pos = -hi - 1 

            else 
                pos = -hi
            end if
        end if

    end function binary_search_int


    function binary_search_real (arr, val, thresh) &
                                 result(pos)

        real(dp), intent(in) :: arr(:)
        real(dp), intent(in) :: val
        real(dp), intent(in) :: thresh
        integer :: pos

        integer :: hi, lo

        ! The search range
        lo = lbound(arr,1)
        hi = ubound(arr,1)

        ! Account for poor usage (i.e. array len == 0)
        if (hi < lo) then
            pos = -lo
            return
        endif

        ! Narrow the search range down in steps.
        do while (hi /= lo)
            pos = int(real(hi + lo,sp) / 2_sp)

            if (abs(arr(pos) - val) < thresh) then
                exit
            else if (val > arr(pos)) then
                ! val is "greater" than arr(:len,pos).
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
            if (abs(arr(hi) - val) < thresh) then
                pos = hi
            else if (val > arr(hi)) then
                pos = -hi - 1
            else
                pos = -hi
            endif
        endif

    end function binary_search_real

    function binary_search_custom (arr, val, cf_len, custom_gt) &
                                                     result(pos)
        use constants, only: n_int
        use DetBitOps, only: DetBitLt

        interface
            pure function custom_gt(a, b) result (ret)
                use constants, only: n_int
                logical :: ret
                integer(kind=n_int), intent(in) :: a(:), b(:)
            end function
        end interface

        integer(kind=n_int), intent(in) :: arr(:,:)
        integer(kind=n_int), intent(in) :: val(:)
        integer, intent(in), optional :: cf_len
        integer :: data_lo, data_hi, val_lo, val_hi
        integer :: pos

        integer :: hi, lo

        ! The search range
        lo = lbound(arr,2)
        hi = ubound(arr,2)

        ! Account for poor usage (i.e. array len == 0)
        if (hi < lo) then
            pos = -lo
            return
        endif

        ! Have we specified how much to look at?
        data_lo = lbound(arr, 1)
        val_lo = lbound(val, 1)
        if (present(cf_len)) then
            data_hi = data_lo + cf_len - 1
            val_hi = val_lo + cf_len - 1
        else
            data_hi = ubound(arr, 1)
            val_hi = ubound(val, 1)
        endif

        ! Narrow the search range down in steps.
        do while (hi /= lo)
            pos = int(real(hi + lo,sp) / 2_sp)

            if (DetBitLT(arr(data_lo:data_hi,pos), val(val_lo:val_hi)) == 0) then
                exit
            else if (custom_gt(val(val_lo:val_hi), arr(data_lo:data_hi,pos))) then
                ! val is "greater" than arr(:len,pos).
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
            if (DetBitLT(arr(data_lo:data_hi,hi), val(val_lo:val_hi)) == 0) then
                pos = hi
            else if (custom_gt(val(val_lo:val_hi), arr(data_lo:data_hi,hi))) then
                pos = -hi - 1
            else
                pos = -hi
            endif
        endif

    end function binary_search_custom

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
       integer :: record_length_loc
       inquire(iolength=record_length_loc) bytes
!       record_length = (bytes/4)*record_length
       record_length = (bytes/sizeof_int)*int(record_length_loc,sizeof_int)
! 8 indicates 8-byte words I think
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

    subroutine get_unique_filename(stem, tincrement, tnext, istart, filename, &
                                   ext)

        ! Find a filename which is either the "newest" or the next to be used.
        ! The filename is assumed to be stem.x, where x is an integer.

        ! In:
        !    stem: stem of the filename.
        !    tincrement: the filename is given as stem.x if true, otherwise the
        !        filename is simply set to be equal to stem.
        !    tnext: the next unused filename is found if true, else the
        !        filename is set to be stem.x where stem.x exists and stem.x+1
        !        doesn't and x is greater than istart
        !    istart: the integer of the first x value to check.
        !        If istart is negative, then the filename is set to be stem.x,
        !        where x = |istart+1|.  This overrides everything else.
        !    ext: The file extension. Appended after the numbers.
        ! Out:
        !    filename.

        character(*), intent(in) :: stem
        logical, intent(in) :: tincrement, tnext
        integer, intent(in) :: istart
        character(*), intent(out) :: filename
        character(*), intent(in), optional :: ext

        integer :: i
        logical :: exists

        if (tincrement) then
            i = istart
            exists = .true.
            do while (exists)
                call append_ext(stem, i, filename)
                if (present(ext)) filename = trim(filename) // ext
                inquire(file=filename,exist=exists)
                i = i + 1
            end do
            if (.not.tnext) then
                ! actually want the last file which existed.
                ! this will return stem.istart if stem.istart doesn't exist.
                i = max(istart,i - 2)
                call append_ext(stem, i, filename)
                if (present(ext)) filename = trim(filename) // ext
            end if
        else
            filename = stem
            if (present(ext)) filename = trim(filename) // ext
        end if

        if (.not.tnext) then
            inquire(file=filename,exist=exists)
            if (.not.exists) then
                inquire(file=stem,exist=exists)
                if (exists) then
                    filename = stem
                    if (present(ext)) filename = trim(filename) // ext
                endif
            end if
        end if

        if (istart < 0) then
            call append_ext(stem, abs(i+1), filename)
            if (present(ext)) filename = trim(filename) // ext
        end if

    end subroutine get_unique_filename

    function get_free_unit() result(free_unit)

        ! Returns:
        !    The first free file unit above 10 and less than or equal to
        !    the paramater max_unit (currently set to 200).
        !
        !    If max_unit is exceeded, the function returns -1

        integer, parameter :: max_unit = 100
        integer :: free_unit
        integer :: i
        logical :: t_open, t_exist

        free_unit = -1
        do i = 10, max_unit
            inquire(unit=i, opened=t_open, exist=t_exist)
            if (.not.t_open .and. t_exist) then
                free_unit = i
                exit
            end if
        end do
        if (i == max_unit+1) call stop_all('get_free_unit','Cannot find a free unit below max_unit.')

    end function get_free_unit

    function error_function_c(argument) result (res)

        use constants, only: dp
        use iso_c_hack
        implicit none 
      
        real(dp), intent(in) :: argument
        real(dp) :: res

        !interface
        !    pure function erfc_lm(x) bind(c, name='erfc') result (ret)
        !        use iso_c_hack
        !        implicit none
        !        real(c_double) :: ret
        !            real(c_double), intent(in), value :: x
        !    end function erfc_lm
        !!end interface

        !res = erfc_lm(real(argument, c_double))
        res = erfc_local (real(argument, c_double))
    end function error_function_c



    function error_function(argument) result(res)
        
        use constants, only: dp
        use iso_c_hack
        implicit none 
      
        real(dp), intent(in) :: argument
        real(dp) :: res

!        interface
!                pure function erf_lm(x) bind(c, name='erf') result(ret)
!                use iso_c_hack
!                implicit none
!                real(c_double) :: ret
!                real(c_double), intent(in), value :: x
!            end function erf_lm
!        end interface
!        res = erf_lm(real(argument, c_double))
        res = erf_local(real(argument, c_double))

    end function error_function

    pure subroutine find_next_comb(comb, k, n, finish)

        integer, intent(in) :: k, n
        integer, intent(inout) :: comb(k)
        logical, intent(out) :: finish
        integer :: i

        if (k == 0 .or. n == 0) then
            finish = .true.
            return
        else if (comb(1) > n-k) then
            finish = .true.
            return
        else
            finish = .false.
        end if

        i = k
        comb(i) = comb(i) + 1

        do
            if (i < 1 .or. comb(i) < n-k+i+1) exit
            i = i-1
            comb(i) = comb(i) + 1
        end do

        do i = i+1, k
            comb(i) = comb(i-1) + 1
        end do

    end subroutine find_next_comb

    function neci_etime(time) result(ret)

        ! Return elapsed time for timing and calculation ending purposes.

        real(dp), intent(out) :: time(2)
        real(dp) :: ret

#ifdef __IFORT
        ! intels etime takes a real(4)
        real(4) :: ioTime(2)
        ! Ifort defines etime directly in its compatibility modules. 
        ! Avoid timing inaccuracies from using cpu_time on cerebro.
        ret = real(etime(ioTime),dp)
        time = real(ioTime,dp)
#else
#ifdef BLUEGENE_HACKS
        time = 0.0_dp
        ret = 0.0_dp
#else
        ! use MPI_WTIME - etime returns wall-clock time on multi-processor
        ! environments, so keep it consistent
        ret = MPI_WTIME()
        time(1) = ret
        time(2) = real(0.0,dp)
#endif
#endif

    end function neci_etime

end module

!Hacks for compiler specific system calls.

    integer function neci_iargc()
        implicit none
#if defined(CBINDMPI) && !defined(MOLPRO)
        interface
            function c_argc () result(ret) bind(c)
                use iso_c_hack
                integer(c_int) :: ret
            end function
        end interface
        neci_iargc = c_argc()
#elif defined(MOLPRO_f2003)
        integer command_argument_count
        neci_iargc=command_argument_count()
#elif defined(MOLPRO)
        integer iargc
        neci_iargc=iargc()
#else
        integer :: command_argument_count
        neci_iargc = command_argument_count ()
#endif
    end function



    subroutine neci_getarg (i, str)

#ifdef NAGF95
        use f90_unix_env, only: getarg
#endif
        use constants
        use iso_c_hack
        use util_mod
        implicit none
        integer, intent(in) :: i
        character(len=*), intent(out) :: str

#if defined(__OPEN64__) || defined(__PATHSCALE__)
        integer(int32) :: j
#else
        integer :: j
#endif

        ! Eliminate compiler warnings
        j = i

#if defined(CBINDMPI) && !defined(MOLPRO)
        ! Define interfaces that we need
        interface
            pure function c_getarg_len (i) result(ret) bind(c)
                use iso_c_hack
                integer(c_int), intent(in), value :: i
                integer(c_int) :: ret
            end function
            pure subroutine c_getarg (i, str) bind(c)
                use iso_c_hack
                integer(c_int), intent(in), value :: i
                character(c_char), intent(out) :: str
            end subroutine
        end interface
        character(len=c_getarg_len(int(i, c_int))) :: str2

        call c_getarg (int(i, c_int), str2)

        str = str2

#elif defined NAGF95
        call getarg(i, str)
#elif defined(MOLPRO) && !defined(MOLPRO_f2003)
        call getarg(i, str)
#elif defined(BLUEGENE_HACKS)
        call getarg(int(i, 4), str)
#elif defined(__OPEN64__) || defined(__PATHSCALE__)
        j = i
        call get_command_argument (j, str)
#else
        call get_command_argument (i, str)
#endif

    end subroutine neci_getarg


    subroutine neci_flush(un)
#ifdef MOLPRO
    implicit none
    integer, intent(in) :: un
    flush(un)
#else
#ifdef NAGF95
    USe f90_unix, only: flush
    use constants, only: int32
#endif
    implicit none
    integer, intent(in) :: un
#ifdef NAGF95
    integer(kind=int32) :: dummy
#endif
#ifdef BLUEGENE_HACKS
        call flush_(un)
#else
#ifdef NAGF95
        dummy=un
        call flush(dummy)
#else
        call flush(un)
#endif
#endif
#endif
    end subroutine neci_flush


    integer function neci_system(str)
#ifdef NAGF95
    Use f90_unix_proc, only: system
#endif
        character(*), intent(in) :: str
#ifndef NAGF95
        integer :: system
        neci_system=system(str)
#else
        call system(str)
        neci_system=0
#endif
    end function neci_system

! Hacks for the IBM compilers on BlueGenes.
! --> The compiler intrinsics are provided as flush_, etime_, sleep_ etc.
! --> We need to either change the names used in the code, or provide wrappers
#ifdef BLUEGENE_HACKS
! I presume that the function cpu_time will work here?
! If not, simply add BLUEGENE_HACKS to the neci_etime above.
!    function etime (t) result(ret)
!        implicit none
!        real(4) :: t(2), etime_, ret
!        ret = etime_(t)
!    end function
    function hostnm (nm) result(ret)
        implicit none
        integer :: ret, hostnm_
        character(255) :: nm
        ret = hostnm_(nm)
    end function

#endif

#ifdef CRAY_ETIME

    function etime (tarr) result (tret)
        implicit none
        real(4) :: tarr(2), tret, second

        tret = second()
        tarr = tret
    end function

#endif

#ifdef __GFORTRAN__
    function g_loc (var) result(addr)

        use iso_c_binding

        integer, target :: var
        type(c_ptr) :: addr
        
        addr = c_loc(var)

    end function
#endif

