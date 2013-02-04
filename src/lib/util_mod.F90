module util_mod
    use util_mod_comparisons
    use util_mod_cpts
    use constants, only: dp, lenof_sign,sizeof_int
    use dSFMT_interface, only: genrand_real2_dSFMT
    use iso_c_hack
    implicit none

    interface
        pure function strlen_wrap (str) result(len) bind(c)
            use iso_c_hack
            implicit none
            character(c_char), intent(in) :: str(*)
            integer(c_int) :: len
        end function
        pure function erf (x) result(e) bind(c)
            use iso_c_hack
            implicit none
            real(c_double), intent(in) :: x
            real(c_double) :: e
        end function
        pure function erfc (x) result(ec) bind(c)
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

        if (res /= 0) then
            if (abs(res) > genrand_real2_dSFMT()) &
                i = i + nint(sign(1.0_dp, r))
        end if

    end function

    subroutine print_cstr (str) bind(c, name='print_cstr')

        ! Write a string outputted by calling fort_printf in C.
        ! --> Ensure that it is redirected to the same place as the normal
        !     STDOUT within fortran.

        character(c_char), intent(in) :: str(*)
        integer :: l, i

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
        integer(int32), intent(in) :: sgn(lenof_sign)

        if(lenof_sign.eq.1) then
            abs_int4_sign = abs(sgn(1))
        else
            abs_int4_sign = real(int(sqrt(real(sgn(1),dp)**2 + real(sgn(lenof_sign),dp)**2)), dp)
            ! The integerisation here is an approximation, but one that is 
            ! used in the integer algorithm, so is retained in this real 
            ! version of the algorithm
        endif
    end function abs_int4_sign

!routine to calculation the absolute magnitude of a complex integer(int64) variable (to nearest integer)
    pure integer(kind=int64) function abs_int8_sign(wsign)
        integer(kind=int64), dimension(lenof_sign), intent(in) :: wsign

        if(lenof_sign.eq.1) then
            abs_int8_sign=abs(wsign(1))
        else
            abs_int8_sign=nint(sqrt(real(wsign(1),dp)**2+real(wsign(lenof_sign),dp)**2),int64)
        endif
    end function abs_int8_sign

    pure real(dp) function abs_real_sign (sgn)
        real(dp), intent(in) :: sgn(lenof_sign)

        if (lenof_sign == 1) then
            abs_real_sign = abs(sgn(1))
        else
            abs_real_sign = real(nint(sqrt(sum(sgn ** 2))), dp)
        end if
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

        if ( (r == 0) .and. (r * 1 == 1) ) then
            isnan_neci = .true.
        else
            isnan_neci = .false.
        endif
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
        p = ceiling(r+p)

        if (p < 10) then
            write (fmt1,'("i",i1)') p
        else if (r < 100) then
            write (fmt1,'("i",i2)') p
        else
            ! By this point we'll have hit integer overflow anyway...
            write (fmt1,'("i",i3)') p
        end if

    end function int_fmt

!--- Searching ---

    function binary_search_first_ge (arr, val) result(pos)

        ! Find the first element in an array which is >= val, in an array
        ! which has been sorted.
        !
        ! If there is no such element, the function returns -1.
        ! TODO: Should this only happen in debug mode?

        integer, intent(in) :: val
        integer, intent(in) :: arr(:)
        integer :: pos

        integer :: hi, lo

        ! The search range
        lo = lbound(arr, 1)
        hi = ubound(arr, 1)

        ! Test if such an element exists
        if (arr(hi) < val) then
            pos = -1
            return
        endif

        do while (hi /= lo)
            pos = int(real(hi + lo) / 2)

            if (arr(pos) >= val) then
                hi = pos
            else
                lo = pos + 1
            endif
        enddo

        ! Return the converged value.
        pos = hi

    end function

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
            pos = int(real(hi + lo) / 2)

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

    function binary_search_custom (arr, val, cf_len, custom_gt) &
                                                     result(pos)
        !use bit_reps, only: NIfD
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
            pos = int(real(hi + lo) / 2)

            if (DetBitLT(arr(data_lo:data_hi,pos), val(val_lo:val_hi), &
                    use_flags_opt = .false.) == 0) then
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
            if (DetBitLT(arr(data_lo:data_hi,hi), val(val_lo:val_hi), &
                    use_flags_opt = .false.) == 0) then
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
        res = erfc (real(argument, c_double))
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
        res = erf(real(argument, c_double))

    end function error_function


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
#elif defined(__OPEN64__) || defined(__PATHSCALE__)
        j = i
        call get_command_argument (j, str)
#else
        call get_command_argument (i, str)
#endif

    end subroutine neci_getarg


    subroutine neci_flush(un)
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
    end subroutine neci_flush


    function neci_etime(time) result(ret)
      use constants, only: sp
      real(sp) :: ret
      real(sp) :: time(2)
      call cpu_time(ret)
      time(1) = ret
      time(2) = real(0.0,sp)
    end function neci_etime

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
        character(8) :: nm
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

