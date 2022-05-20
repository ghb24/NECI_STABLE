#include "macros.h"
#:include "../algorithms.fpph"
#:include "../macros.fpph"

#:set countable_types = {'integer': {'int32', 'int64'}}
#:set primitive_types = {'integer': {'int32', 'int64'}, 'real': {'sp', 'dp'}, 'complex': {'sp', 'dp'}}
#:set log_entry = {'logical':{''}}
#:set extended_types = dict(primitive_types, **log_entry)
#:set ops = {'integer': '==', 'real': '.isclose.', 'complex': '.isclose.', 'logical': '.eqv.'}

module util_mod
    use util_mod_comparisons, only: operator(.arrgt.), operator(.arrlt.), arr_gt, arr_lt
    use util_mod_numerical
    use util_mod_cpts
    use util_mod_epsilon_close
    use binomial_lookup, only: factrl => factorial, binomial_lookup_table_i64
#ifdef GFORTRAN_
    use binomial_lookup, only: binomial_lookup_table_i128
#endif
    use fmt_utils
    use dSFMT_interface, only: genrand_real2_dSFMT
    use constants
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_double

    ! We want to use the builtin etime intrinsic with ifort to
    ! work around some broken behaviour.
#if defined(IFORT_) || defined(INTELLLVM_)
    use ifport, only: etime
#endif
    implicit none

#if defined(IFORT_) || defined(INTELLLVM_)
    public
#else
    private
#endif

    public :: get_nan, isnan_neci, factrl, choose_i64, NECI_icopy, operator(.implies.), &
        abs_l1, abs_sign, near_zero, operator(.isclose.), operator(.div.), &
        stochastic_round, stochastic_round_r
#ifdef GFORTRAN_
    public :: choose_i128
#endif

    public :: error_function, error_function_c, stop_all, toggle_lprof
    public :: arr_2d_ptr, ptr_abuse_1d, ptr_abuse_2d, ptr_abuse_scalar
    public :: neci_etime, get_free_unit, int_fmt,&
        strlen_wrap, record_length, open_new_file, &
        append_ext, get_unique_filename, neci_flush, print_cstr_local, &
        stats_out
    public :: arr_lt, arr_gt, operator(.arrlt.), operator(.arrgt.), &
        find_next_comb, binary_search_ilut, binary_search_custom, binary_search_first_ge, &
        cumsum, pairswap, swap, lex_leq, lex_geq, &
        get_permutations, custom_findloc, addToIntArray, fuseIndex, linearIndex, &
        getSpinIndex, binary_search_int, binary_search_real

    public :: EnumBase_t


    interface
        ! NOTE: A stop all is of course state-changing, but even
        !   the Fortran standard allows an `error stop`.
        pure subroutine stop_all(sub_name, error_msg)
            character(*), intent(in) :: sub_name, error_msg
        end subroutine

        subroutine neci_flush(n)
            integer, intent(in) :: n
        end subroutine
    end interface

    interface binary_search_int
    #:for kind in primitive_types['integer']
        module procedure binary_search_int_${kind}$
    #:endfor
    end interface

    interface operator(.implies.)
        module procedure implies
    end interface

    interface choose_i64
    #:for kind in primitive_types['integer']
        module procedure choose_i64_${kind}$
    #:endfor
    end interface

#ifdef GFORTRAN_
    interface choose_i128
    #:for kind in primitive_types['integer']
        module procedure choose_i128_${kind}$
    #:endfor
    end interface
#endif


    interface
        pure function strlen_wrap(str) result(len) bind(c)
            import :: c_char, c_int
            implicit none
            character(c_char), intent(in) :: str(*)
            integer(c_int) :: len
        end function
        pure function erf_local(x) result(e) bind(c, name='erf')
            import :: c_double
            implicit none
            real(c_double), intent(in) :: x
            real(c_double) :: e
        end function
        pure function erfc_local(x) result(ec) bind(c, name='erfc')
            import :: c_double
            implicit none
            real(c_double), intent(in) :: x
            real(c_double) :: ec
        end function
        subroutine toggle_lprof() bind(C, name='toggle_lprof')
            !! This call toggles the profiling using MAQAO
            !!
            !! It can be for example used to switch on profiling
            !! right before the main loop
            !! and to switch it off directly afterwards.
        end subroutine
    end interface

    interface operator(.div.)
        module procedure div_int32, div_int64
#ifdef GFORTRAN_
        module procedure div_int128
#endif
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
    end interface abs_l1

    interface fuseIndex
        module procedure fuseIndex_int32
        module procedure fuseIndex_int64
    end interface fuseIndex

    interface swap
        module procedure swap_int64
        module procedure swap_int32
    end interface swap

    interface custom_findloc
        #:for type, kinds in extended_types.items()
        #:for kind in kinds
        module procedure custom_findloc_${type}$_${kind}$
        #:endfor
        #:endfor
    end interface custom_findloc

    interface cumsum
        #:for type, kinds in primitive_types.items()
        #:for kind in kinds
        module procedure cumsum_${type}$_${kind}$
        #:endfor
        #:endfor
    end interface

    type, abstract :: EnumBase_t
        integer :: val
    contains
        private
        procedure :: eq_EnumBase_t
        procedure :: neq_EnumBase_t
        generic, public :: operator(==) => eq_EnumBase_t
        generic, public :: operator(/=) => neq_EnumBase_t
    end type


contains

    function stochastic_round(r) result(i)

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

        if(abs(res) >= 1.0e-12_dp) then
            if(abs(res) > genrand_real2_dSFMT()) &
                i = i + nint(sign(1.0_dp, r))
        end if

    end function

    elemental function stochastic_round_r(num, r) result(i)

        ! Perform the stochastic rounding of the above function where the
        ! random number is already specified.

        real(dp), intent(in) :: num, r
        integer :: i
        real(dp) :: res

        i = int(num)
        res = num - real(i, dp)

        if(abs(res) >= 1.0e-12_dp) then
            if(abs(res) > r) &
                i = i + nint(sign(1.0_dp, num))
        end if

    end function stochastic_round_r

    subroutine print_cstr_local(str, l)

        character(c_char), intent(in) :: str(*)
        integer, intent(in) :: l
        character(len=l) :: tmp_s

        tmp_s = transfer(str(1:l), tmp_s)
        write(6, '(a)', advance='no') tmp_s

    end subroutine

    ! routine to calculation the absolute magnitude of a complex integer
    ! variable (to nearest integer)
    pure real(dp) function abs_int4_sign(sgn)
        integer(int32), intent(in) :: sgn(lenof_sign / inum_runs)

#ifdef CMPLX_
        abs_int4_sign = real(int(sqrt(real(sgn(1), dp)**2 + real(sgn(2), dp)**2)), dp)
        ! The integerisation here is an approximation, but one that is
        ! used in the integer algorithm, so is retained in this real
        ! version of the algorithm
#else
        abs_int4_sign = abs(sgn(1))
#endif
    end function abs_int4_sign

!routine to calculation the absolute magnitude of a complex integer(int64) variable (to nearest integer)
    pure integer(kind=int64) function abs_int8_sign(wsign)
        integer(kind=int64), dimension(lenof_sign/inum_runs), intent(in) :: wsign

#ifdef CMPLX_
        abs_int8_sign = nint(sqrt(real(wsign(1), dp)**2 + real(wsign(2), dp)**2), int64)
#else
        abs_int8_sign = abs(wsign(1))
#endif
    end function abs_int8_sign

    pure real(dp) function abs_real_sign(sgn)
        real(dp), intent(in) :: sgn(lenof_sign / inum_runs)
#ifdef CMPLX_
        abs_real_sign = real(nint(sqrt(sum(sgn**2))), dp)
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

    pure function abs_l1_sp(val) result(ret)

        real(sp), intent(in) :: val
        real(sp) :: ret

        ret = abs(val)

    end function

    pure function abs_l1_dp(val) result(ret)

        real(dp), intent(in) :: val
        real(dp) :: ret

        ret = abs(val)

    end function

    pure function abs_l1_csp(val) result(ret)

        complex(sp), intent(in) :: val
        real(sp) :: ret

        ret = abs(real(val, sp)) + abs(aimag(val))

    end function

    pure function abs_l1_cdp(val) result(ret)

        complex(dp), intent(in) :: val
        real(dp) :: ret

        ret = abs(real(val, dp)) + abs(aimag(val))

    end function

!--- Array utilities ---

    SUBROUTINE NECI_ICOPY(N, A, IA, B, IB)
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
        INTEGER, INTENT(IN) :: N, IA, IB
        INTEGER, INTENT(IN) :: A(IA * N)
        INTEGER, INTENT(OUT) :: B(IB * N)
!        Variables
        INTEGER I, IAX, IBX

        DO I = 1, N
            IAX = (I - 1) * IA + 1
            IBX = (I - 1) * IB + 1
            B(IBX) = A(IAX)
        ENDDO

        RETURN
    END SUBROUTINE NECI_ICOPY

    subroutine addToIntArray(arr, ind, elem)
        integer, intent(inout), allocatable :: arr(:)
        integer, intent(in) :: ind, elem

        integer, allocatable :: tmp(:)
        integer :: nelems

        if(allocated(arr)) then
            nelems = size(arr)

            if(ind > nelems) then
                ! resize the array
                allocate(tmp(nelems))
                tmp = arr
                deallocate(arr)
                allocate(arr(ind), source=0)
                arr(1:nelems) = tmp(1:nelems)
            endif
        else
            allocate(arr(ind), source=0)
        endif

        arr(ind) = elem

    end subroutine addToIntArray

    !------------------------------------------------------------------------------------------!

    !> Custom implementation of the findloc intrinsic (with somewhat reduced functionality)
    !! as it requires fortran2008 support and is thus not available for some relevant compilers
    #:for type, kinds in extended_types.items()
    #:for kind in kinds
    #:set this_op = ops[type]
    pure function custom_findloc_${type}$_${kind}$(arr, val, back) result(loc)
        @{get_decl(${type}$, ${kind}$,0)}@, intent(in) :: arr(:)
        @{get_decl(${type}$, ${kind}$,0)}@, intent(in) :: val
        logical, intent(in), optional :: back
        integer :: loc

        integer :: i, first, last, step
        logical :: back_

        def_default(back_, back, .false.)

        if(back_) then
            first = size(arr)
            last = 1
            step = -1
        else
            first = 1
            last = size(arr)
            step = 1
        end if

        loc = 0
        do i = first, last, step
            if(arr(i) ${this_op}$ val) then
                loc = i
                return
            endif
        end do
    end function custom_findloc_${type}$_${kind}$
    #:endfor
    #:endfor

!--- Indexing utilities

    pure function fuseIndex_int32(q, p) result(ind)
        ! fuse p,q into one symmetric index
        ! the resulting index is not contigious in p or q
        ! Input: p,q - 2d-array indices
        ! Output: ind - 1d-array index assuming the array is symmetric w.r. p<->q
        integer(int32), intent(in) :: p, q
        integer(int32) :: ind

        ! qp and pq are considered to be the same index
        ! -> permutational symmetry
        ! implemented in terms of fuseIndex_int64
        ind = int(fuseIndex_int64(int(q, int64), int(p, int64)))
    end function fuseIndex_int32

!------------------------------------------------------------------------------------------!

    pure function fuseIndex_int64(x, y) result(xy)
        ! create a composite index out of two indices, assuming they are unordered
        ! i.e. their ordering does not matter
        ! Input: p,q - 2d-array indices
        ! Output: ind - 1d-array index assuming the array is symmetric w.r. p<->q
        integer(int64), intent(in) :: x, y
        integer(int64) :: xy

        if(x < y) then
            xy = x + y * (y - 1) / 2
        else
            xy = y + x * (x - 1) / 2
        endif
    end function fuseIndex_int64

!------------------------------------------------------------------------------------------!

    elemental subroutine swap_int32(a, b)
        ! exchange the value of two integers a,b
        ! Input: a,b - integers to swapp (on return, a has the value of b on call and vice versa)
        integer(int32), intent(inout) :: a, b
        integer(int32) :: tmp

        tmp = a
        a = b
        b = tmp
    end subroutine swap_int32

!------------------------------------------------------------------------------------------!

    elemental subroutine swap_int64(a, b)
        ! exchange the value of two integers a,b
        ! Input: a,b - integers to swapp (on return, a has the value of b on call and vice versa)
        integer(int64), intent(inout) :: a, b
        integer(int64) :: tmp

        tmp = a
        a = b
        b = tmp
    end subroutine swap_int64

!------------------------------------------------------------------------------------------!

    pure subroutine pairSwap(a, i, b, j)
        ! exchange a pair of integers
        integer(int64), intent(inout) :: a, i, b, j

        call swap(a, b)
        call swap(i, j)
    end subroutine pairSwap

!------------------------------------------------------------------------------------------!

    function linearIndex(p, q, dim) result(ind)
        ! fuse p,q into one contiguous index
        ! the resulting index is contiguous in q
        ! Input: p,q - 2d-array indices
        !        dim - dimension of the underlying array in q-direction
        ! Output: ind - contiguous 1d-array index
        integer, intent(in) :: p, q, dim
        integer :: ind

        ind = q + (p - 1) * dim
    end function linearIndex

    pure elemental function getSpinIndex(orb) result(ms)
        ! return a spin index of the orbital orb which can be used to address arrays
        ! Input: orb - spin orbital
        ! Output: ms - spin index of orb with the following values:
        !              0 - alpha
        !              1 - beta
        integer, intent(in) :: orb
        integer :: ms

        ms = mod(orb, 2)
    end function getSpinIndex

!--- Numerical utilities ---

    ! If all of the compilers supported ieee_arithmetic
    ! --> could use ieee_value(1.0_dp, ieee_quiet_nan)
    real(dp) function get_nan()
        real(dp) :: a, b
        a = 1
        b = 1
        get_nan = log(a - 2 * b)
    end function

    ! If all of the compilers supported ieee_arithmetic
    ! --> could use ieee_is_nan (r)
    elemental logical function isnan_neci(r)
        real(dp), intent(in) :: r

#ifdef GFORTRAN_
        isnan_neci = isnan(r)
#else
        if((r == 0) .and. (r * 1 == 1)) then
            isnan_neci = .true.
        else
            isnan_neci = .false.
        endif
#endif
    end function

        !> @brief
        !> Calculate 1 + ... + n
        integer elemental function gauss_sum(n)
            integer, intent(in) :: n
            gauss_sum = (n * (n + 1)) .div. 2
        end function

        !> @brief
        !> Get the index in the binomial_lookup_table
        integer elemental function get_index(n, k)
            integer, intent(in) :: n, k
            get_index = gauss_sum((n - 3) .div. 2) + gauss_sum((n - 4) .div. 2) + k - 1
        end function

#:for kind in primitive_types['integer']
    ! Unfortunately there are no recursive elemental functions in Fortran.
    recursive pure function choose_i64_${kind}$(n, r, signal_overflow) result(res)
        !! Return the binomail coefficient nCr(n, r)
        integer(${kind}$), intent(in) :: n, r
        logical, intent(in), optional :: signal_overflow
            !! If true then the function returns -1 instead of aborting
            !! when overflow is encountered.
        integer(int64) :: res
        integer(int64) :: k

        character(*), parameter :: this_routine = "choose_i64"

        ! NOTE: This is highly optimized. If you change something, please time it!

        @:pure_ASSERT(n >= 0_${kind}$, import_stop_all=False)
        @:pure_ASSERT(r >= 0_${kind}$, import_stop_all=False)

        if(r > n) then
            res = 0_int64
            return
        end if

        k = int(merge(r, n - r, r <= n - r), kind=int64)

        if (k == 0) then
            res = 1_int64
        else if (k == 1) then
            res = int(n, int64)
        else if (n <= 66) then
            ! use lookup table
            res = binomial_lookup_table_i64(get_index(int(n), int(k)))
        else
            block
                integer(int64) :: prev
                prev = choose_i64_${kind}$(n - 1, r - 1, signal_overflow)
            ! Note that the recursion stops at n = 66
                res = (prev * n) .div. k
                check_for_overflow: if (prev < 0 .or. res < 0) then
                    if (present(signal_overflow)) then
                        if (signal_overflow) then
                            res = -1
                        else
#if defined(IFORT_) || defined(INTELLLVM_)
                            error stop 'Binomial coefficient exceeds range of int64.'
#else
                            call stop_all(this_routine, 'Binomial coefficient exceeds range of int64.')
#endif
                        end if
                    else
#if defined(IFORT_) || defined(INTELLLVM_)
                            error stop 'Binomial coefficient exceeds range of int64.'
#else
                            call stop_all(this_routine, 'Binomial coefficient exceeds range of int64.')
#endif
                    end if
                end if check_for_overflow
            end block
        end if
    end function

#ifdef GFORTRAN_
    recursive pure function choose_i128_${kind}$(n, r, signal_overflow) result(res)
        !! Return the binomail coefficient nCr(n, r)
        integer(${kind}$), intent(in) :: n, r
        logical, intent(in), optional :: signal_overflow
            !! If true then the function returns -1 instead of aborting
            !! when overflow is encountered.
        integer(int128) :: res
        integer(int128) :: k
        character(*), parameter :: this_routine = "choose_i128"

        ! NOTE: This is highly optimized. If you change something, please time it!

        @:pure_ASSERT(n >= 0_${kind}$, import_stop_all=False)
        @:pure_ASSERT(r >= 0_${kind}$, import_stop_all=False)

        if(r > n) then
            res = 0_int128
            return
        end if

        k = int(merge(r, n - r, r <= n - r), kind=int128)

        if (k == 0) then
            res = 1_int128
        else if (k == 1) then
            res = int(n, int128)
        else if (n <= 130) then
            ! use lookup table
            res = binomial_lookup_table_i128(get_index(int(n), int(k)))
        else
            ! Note that the recursion stops at n = 130
            block
                integer(int128) :: prev
                prev = choose_i128_${kind}$(n - 1, r - 1, signal_overflow)
            ! Note that the recursion stops at n = 66
                res = (prev * n) .div. k
                check_for_overflow: if (prev < 0 .or. res < 0) then
                    if (present(signal_overflow)) then
                        if (signal_overflow) then
                            res = -1
                        else
#if defined(IFORT_) || defined(INTELLLVM_)
                            error stop 'Binomial coefficient exceeds range of int128.'
#else
                            call stop_all(this_routine, 'Binomial coefficient exceeds range of int128.')
#endif
                        end if
                    else
#if defined(IFORT_) || defined(INTELLLVM_)
                            error stop 'Binomial coefficient exceeds range of int128.'
#else
                            call stop_all(this_routine, 'Binomial coefficient exceeds range of int128.')
#endif
                    end if
                end if check_for_overflow
            end block
        end if
    end function
#endif
#:endfor

    elemental integer(int32) function div_int32(a, b)
        integer(int32), intent(in) :: a, b
#ifdef WARNING_WORKAROUND_
        div_int32 = int(real(a, kind=sp) / real(b, kind=sp), kind=int32)
#else
        div_int32 = a / b
#endif
    end function

    elemental integer(int64) function div_int64(a, b)
        integer(int64), intent(in) :: a, b
#ifdef WARNING_WORKAROUND_
        div_int64 = int(real(a, kind=dp) / real(b, kind=dp), kind=int64)
#else
        div_int64 = a / b
#endif
    end function

#ifdef GFORTRAN_
    elemental integer(int128) function div_int128(a, b)
        integer(int128), intent(in) :: a, b
#ifdef WARNING_WORKAROUND_
        div_int128 = int(real(a, kind=dp) / real(b, kind=dp), kind=int128)
#else
        div_int128 = a / b
#endif
    end function
#endif

!--- Comparison of subarrays ---

    logical pure function det_int_arr_gt(a, b, len)
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

        if(present(len)) then
            llen = len
        else
            llen = min(size(a), size(b))
        endif

        ! Sort by the first integer first ...
        i = 1
        do i = 1, llen
            if(a(i) /= b(i)) exit
        enddo

        ! Make the comparison
        if(i > llen) then
            det_int_arr_gt = .false.
        else
            if(a(i) > b(i)) then
                det_int_arr_gt = .true.
            else
                det_int_arr_gt = .false.
            endif
        endif
    end function det_int_arr_gt

    logical pure function det_int_arr_eq(a, b, len)
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
        if(present(len)) then
            llen = len
        else
            if(size(a) /= size(b)) then
                det_int_arr_eq = .false.
                return
            endif
            llen = size(a)
        endif

        ! Run through the arrays. Return if they differ at any point.
        do i = 1, llen
            if(a(i) /= b(i)) then
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
    function binary_search_ilut(arr, val, cf_len) result(pos)
        use constants, only: n_int

        integer(kind=n_int), intent(in) :: arr(:, :)
        integer(kind=n_int), intent(in) :: val(:)
        integer, intent(in), optional :: cf_len
        integer :: data_lo, data_hi, val_lo, val_hi
        integer :: pos

        integer :: hi, lo

        ! The search range
        lo = lbound(arr, 2)
        hi = ubound(arr, 2)

        ! Account for poor usage (i.e. array len == 0)
        if(hi < lo) then
            pos = -lo
            return
        endif

        ! Have we specified how much to look at?
        data_lo = lbound(arr, 1)
        val_lo = lbound(val, 1)
        if(present(cf_len)) then
            data_hi = data_lo + cf_len - 1
            val_hi = val_lo + cf_len - 1
        else
            data_hi = ubound(arr, 1)
            val_hi = ubound(val, 1)
        endif

        ! Narrow the search range down in steps.
        do while(hi /= lo)
            pos = int(real(hi + lo, sp) / 2_sp)

            if(all(arr(data_lo:data_hi, pos) == val(val_lo:val_hi))) then
                exit
            else if(arr_gt(val(val_lo:val_hi), arr(data_lo:data_hi, pos))) then
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
        if(hi == lo) then
            if(all(arr(data_lo:data_hi, hi) == val(val_lo:val_hi))) then
                pos = hi
            else if(arr_gt(val(val_lo:val_hi), arr(data_lo:data_hi, hi))) then
                pos = -hi - 1
            else
                pos = -hi
            endif
        endif

    end function binary_search_ilut

    #:for kind in primitive_types['integer']
        pure function binary_search_int_${kind}$(arr, val) result(pos)
            integer(${kind}$), intent(in) :: arr(:)
            integer(${kind}$), intent(in) :: val
            integer(int64) :: pos

            integer(int64) :: hi, lo

            lo = 1
            hi = size(arr)

            if(hi < lo) then
                pos = -lo
                return
            end if

            do while(hi /= lo)
                pos = int((hi + lo) / 2.0_dp, kind=int64)

                if(arr(pos) == val) then
                    exit
                else if(val > arr(pos)) then
                    lo = pos + 1
                else
                    hi = pos
                end if
            end do

            if(hi == lo) then
                if(arr(hi) == val) then
                    pos = hi
                else
                    pos = -1
                end if
            end if

        end function binary_search_int_${kind}$
    #:endfor

    function binary_search_real(arr, val, thresh) &
        result(pos)

        real(dp), intent(in) :: arr(:)
        real(dp), intent(in) :: val
        real(dp), intent(in) :: thresh
        integer :: pos

        integer :: hi, lo

        ! The search range
        lo = lbound(arr, 1)
        hi = ubound(arr, 1)

        ! Account for poor usage (i.e. array len == 0)
        if(hi < lo) then
            pos = -lo
            return
        endif

        ! Narrow the search range down in steps.
        do while(hi /= lo)
            pos = int(real(hi + lo, sp) / 2_sp)

            if(abs(arr(pos) - val) < thresh) then
                exit
            else if(val > arr(pos)) then
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
        if(hi == lo) then
            if(abs(arr(hi) - val) < thresh) then
                pos = hi
            else if(val > arr(hi)) then
                pos = -hi - 1
            else
                pos = -hi
            endif
        endif

    end function binary_search_real

    function binary_search_custom(arr, val, cf_len, custom_gt) &
        result(pos)
        use constants, only: n_int
        use DetBitOps, only: DetBitLt

        interface
            pure function custom_gt(a, b) result(ret)
                use constants, only: n_int
                logical :: ret
                integer(kind=n_int), intent(in) :: a(:), b(:)
            end function
        end interface

        integer(kind=n_int), intent(in) :: arr(:, :)
        integer(kind=n_int), intent(in) :: val(:)
        integer, intent(in), optional :: cf_len
        integer :: data_lo, data_hi, val_lo, val_hi
        integer :: pos

        integer :: hi, lo

        ! The search range
        lo = lbound(arr, 2)
        hi = ubound(arr, 2)

        ! Account for poor usage (i.e. array len == 0)
        if(hi < lo) then
            pos = -lo
            return
        endif

        ! Have we specified how much to look at?
        data_lo = lbound(arr, 1)
        val_lo = lbound(val, 1)
        if(present(cf_len)) then
            data_hi = data_lo + cf_len - 1
            val_hi = val_lo + cf_len - 1
        else
            data_hi = ubound(arr, 1)
            val_hi = ubound(val, 1)
        endif

        ! Narrow the search range down in steps.
        do while(hi /= lo)
            pos = int(real(hi + lo, sp) / 2_sp)

            if(DetBitLT(arr(data_lo:data_hi, pos), val(val_lo:val_hi)) == 0) then
                exit
            else if(custom_gt(val(val_lo:val_hi), arr(data_lo:data_hi, pos))) then
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
        if(hi == lo) then
            if(DetBitLT(arr(data_lo:data_hi, hi), val(val_lo:val_hi)) == 0) then
                pos = hi
            else if(custom_gt(val(val_lo:val_hi), arr(data_lo:data_hi, hi))) then
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
        record_length = (bytes / sizeof_int) * int(record_length_loc)
! 8 indicates 8-byte words I think
    end function record_length

    subroutine append_ext(stem, n, s)

        ! Returns stem.n in s.

        character(*), intent(in) :: stem
        integer, intent(in) :: n
        character(*), intent(out) :: s
        character(10) :: ext

        write(ext, '('//int_fmt(n, 0)//')') n
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

        if(tincrement) then
            i = istart
            exists = .true.
            do while(exists)
                call append_ext(stem, i, filename)
                if(present(ext)) filename = trim(filename)//ext
                inquire(file=filename, exist=exists)
                i = i + 1
            end do
            if(.not. tnext) then
                ! actually want the last file which existed.
                ! this will return stem.istart if stem.istart doesn't exist.
                i = max(istart, i - 2)
                call append_ext(stem, i, filename)
                if(present(ext)) filename = trim(filename)//ext
            end if
        else
            filename = stem
            if(present(ext)) filename = trim(filename)//ext
        end if

        if(.not. tnext) then
            inquire(file=filename, exist=exists)
            if(.not. exists) then
                inquire(file=stem, exist=exists)
                if(exists) then
                    filename = stem
                    if(present(ext)) filename = trim(filename)//ext
                endif
            end if
        end if

        if(istart < 0) then
            call append_ext(stem, abs(i + 1), filename)
            if(present(ext)) filename = trim(filename)//ext
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
            if(.not. t_open .and. t_exist) then
                free_unit = i
                exit
            end if
        end do
        if(i == max_unit + 1) call stop_all('get_free_unit', 'Cannot find a free unit below max_unit.')

    end function get_free_unit

    function error_function_c(argument) result(res)

        use constants, only: dp

        real(dp), intent(in) :: argument
        real(dp) :: res

        res = erfc_local(real(argument, c_double))
    end function error_function_c

    function error_function(argument) result(res)

        use constants, only: dp

        real(dp), intent(in) :: argument
        real(dp) :: res

        res = erf_local(real(argument, c_double))

    end function error_function

    pure subroutine find_next_comb(comb, k, n, finish)

        integer, intent(in) :: k, n
        integer, intent(inout) :: comb(k)
        logical, intent(out) :: finish
        integer :: i

        if(k == 0 .or. n == 0) then
            finish = .true.
            return
        else if(comb(1) > n - k) then
            finish = .true.
            return
        else
            finish = .false.
        end if

        i = k
        comb(i) = comb(i) + 1

        do
            if(i < 1 .or. comb(i) < n - k + i + 1) exit
            i = i - 1
            comb(i) = comb(i) + 1
        end do

        do i = i + 1, k
            comb(i) = comb(i - 1) + 1
        end do

    end subroutine find_next_comb

    function neci_etime(time) result(ret)
#if !defined(IFORT_) && !defined(INTELLLVM_)
        use mpi
#endif
        ! Return elapsed time for timing and calculation ending purposes.

        real(dp), intent(out) :: time(2)
        real(dp) :: ret

#if defined(IFORT_) || defined(INTELLLVM_)
        ! intels etime takes a real(4)
        real(sp) :: ioTime(2)
        ! Ifort defines etime directly in its compatibility modules.
        ! Avoid timing inaccuracies from using cpu_time on cerebro.
        ret = real(etime(ioTime), dp)
        time = real(ioTime, dp)
#else
#ifdef BLUEGENE_HACKS
        time = 0.0_dp
        ret = 0.0_dp
#else
        ! use MPI_WTIME - etime returns wall-clock time on multi-processor
        ! environments, so keep it consistent
        ret = MPI_WTIME()
        time(1) = ret
        time(2) = real(0.0, dp)
#endif
#endif

    end function neci_etime

    subroutine open_new_file(funit, filename)
        integer, intent(in) :: funit
        character(*), intent(in) :: filename
        logical :: exists
        integer :: ierr, i
        character(43) :: filename2
        character(12) :: num
        character(*), parameter :: t_r = 'open_new_file'

        ! If we are doing a normal calculation, move existing fciqmc_stats
        ! files so that they are not overwritten, and then create a new one
        inquire(file=filename, exist=exists)

        if(exists) then

            ! Loop until we find an available spot to move the existing
            ! file to.
            i = 1
            do while(exists)
                write(num, '(i12)') i
                filename2 = trim(adjustl(filename))//"."// &
                            trim(adjustl(num))
                inquire(file=filename2, exist=exists)
                if(i > 10000) &
                    call stop_all(t_r, 'Error finding free fciqmc_stats.*')
                i = i + 1
            end do

            ! Move the file
            call rename(filename, filename2)

        end if

        ! And finally open the file
        open(funit, file=filename, status='unknown', iostat=ierr)

    end subroutine open_new_file

    #:for type, kinds in primitive_types.items()
    #:for kind in kinds
    pure function cumsum_${type}$_${kind}$(X) result(Y)
        ${type}$(${kind}$), intent(in) :: X(:)
        ${type}$(${kind}$) :: Y(lbound(X, 1):ubound(X, 1))

        integer :: i

        if(size(X) /= 0) then
            Y(lbound(X, 1)) = X(lbound(X, 1))
            do i = lbound(X, 1) + 1, ubound(X, 1)
                Y(i) = Y(i - 1) + X(i)
            end do
        end if
    end function
    #:endfor
    #:endfor

    pure function lex_leq(lhs, rhs) result(res)
        integer, intent(in) :: lhs(:), rhs(size(lhs))
        logical :: res
        integer :: i

        res = .true.
        do i = 1, size(lhs)
            if (lhs(i) == rhs(i)) then
                cycle
            else if (lhs(i) < rhs(i)) then
                return
            else if (lhs(i) > rhs(i)) then
                res = .false.
                return
            end if
        end do
    end function

    pure function lex_geq(lhs, rhs) result(res)
        integer, intent(in) :: lhs(:), rhs(size(lhs))
        logical :: res
        integer :: i

        res = .true.
        do i = 1, size(lhs)
            if (lhs(i) == rhs(i)) then
                cycle
            else if (lhs(i) > rhs(i)) then
                return
            else if (lhs(i) < rhs(i)) then
                res = .false.
                return
            end if
        end do
    end function

    !> @brief
    !> Create all possible permutations of [1, ..., n]
    pure function get_permutations(n) result(res)
        integer, intent(in) :: n
        integer :: res(n, factrl(n))

        integer :: tmp(n), i, j, f

        tmp = [(i, i = 1, n)]

        res(:, 1) = tmp
        do f = 2, size(res, 2)
            i = 2
            do while (tmp(i - 1) > tmp(i))
                i = i + 1
            end do
            j = 1
            do while (tmp(j) > tmp(i))
                j = j + 1
            end do
            call swap(tmp(i), tmp(j))

            i = i - 1
            j = 1
            do while (j < i)
                call swap(tmp(i), tmp(j))
                i = i - 1
                j = j + 1
            end do
            res(:, f) = tmp
        end do
    end function

    !> @brief
    !> The logical operator P => Q
    !>
    !> @details
    !>    P  |  Q   |  P => Q   | ¬ P ∨ Q
    !>    -------------------------------
    !>    T  |  T   |     T     |     T
    !>    T  |  F   |     F     |     F
    !>    F  |  T   |     T     |     T
    !>    F  |  F   |     T     |     T
    logical elemental function implies(P, Q)
        logical, intent(in) :: P, Q
        implies = .not. P .or. Q
    end function

    logical elemental function eq_EnumBase_t(this, other)
        class(EnumBase_t), intent(in) :: this, other
        if (.not. SAME_TYPE_AS(this, other)) error stop 'Can only compare objects of same type'
        eq_EnumBase_t = this%val == other%val
    end function

    logical elemental function neq_EnumBase_t(this, other)
        class(EnumBase_t), intent(in) :: this, other
        if (.not. SAME_TYPE_AS(this, other)) error stop 'Can only compare objects of same type'
        neq_EnumBase_t = this%val /= other%val
    end function

end module

!Hacks for compiler specific system calls.

integer function neci_iargc()
    implicit none
    integer :: command_argument_count
    neci_iargc = command_argument_count()
end function

subroutine neci_getarg(i, str)

#ifdef NAGF95
    use f90_unix_env, only: getarg
#endif
    use constants
    use util_mod
    implicit none
    integer, intent(in) :: i
    character(len=*), intent(out) :: str

#if defined(__OPEN64__) || defined(__PATHSCALE__)
    integer(int32) :: j
#else
    integer :: j
#endif

#ifdef WARNING_WORKAROUND_
    j = i
#endif

#if defined NAGF95
    call getarg(i, str)
#elif defined(BLUEGENE_HACKS)
    call getarg(int(i, 4), str)
#elif defined(__OPEN64__) || defined(__PATHSCALE__)
    j = i
    call get_command_argument(j, str)
#else
    call get_command_argument(i, str)
#endif

end subroutine neci_getarg


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
function hostnm(nm) result(ret)
    implicit none
    integer :: ret, hostnm_
    character(255) :: nm
    ret = hostnm_(nm)
end function

#endif

#ifdef CRAY_ETIME

function etime(tarr) result(tret)
    implicit none
    real(4) :: tarr(2), tret, second

    tret = second()
    tarr = tret
end function

#endif

subroutine neci_flush(un)
#ifdef NAGF95
    use f90_unix, only: flush
    use constants, only: int32
#endif
    integer, intent(in) :: un
#ifdef NAGF95
    integer(kind=int32) :: dummy
#endif
#ifdef BLUEGENE_HACKS
    call flush_(un)
#else
#ifdef NAGF95
    dummy = un
    call flush(dummy)
#else
    call flush(un)
#endif
#endif
end subroutine neci_flush

subroutine warning_neci(sub_name,error_msg)
    != Print a warning message in a (helpfully consistent) format.
    !=
    != In:
    !=    sub_name:  calling subroutine name.
    !=    error_msg: error message.
    use, intrinsic :: iso_fortran_env, only: stderr => error_unit
    character(*), intent(in) :: sub_name, error_msg

    write (stderr,'(/a)') 'WARNING.  Error in '//adjustl(sub_name)
    write (stderr,'(a/)') adjustl(error_msg)
end subroutine warning_neci
