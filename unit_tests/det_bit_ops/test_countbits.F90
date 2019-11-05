#include "test_macros.h"

program countbits_tests

    use DetBitOps, only: CountBits_sparse, CountBits_nifty, Countbits, &
                         CountBits_elemental
    use constants
    use fruit
    implicit none

    abstract interface
        function count_bits_t (ilut, nlast, nbitsmax) result(bits)
            use constants
            implicit none
            integer, intent(in) :: nlast
            integer(int64), intent(in) :: ilut(0:nlast)
            integer, intent(in), optional :: nbitsmax
            integer :: bits
        end function
    end interface

    integer :: failed_count

    call init_fruit
    call countbits_drive_tests
    call fruit_finalize

    call get_failed_count(failed_count)
    if (failed_count /= 0) then
        stop -1
    end if

contains

    subroutine countbits_drive_tests
        ! We have multiple variants of the countbits function. Test that they
        ! all pass the tests!
        procedure(count_bits_t), pointer :: proc
        proc => CountBits_sparse
        call test_countbits_fn(proc, "sparse")
        proc => CountBits_nifty
        call test_countbits_fn(proc, "nifty")
        proc => CountBits_elemental
        call test_countbits_fn(proc, "elemental")
    end subroutine

    subroutine test_countbits_fn(proc, subname)
        ! This calls each of the relevant tests on the supplied function
        procedure(count_bits_t), pointer :: proc
        character(*), intent(in) :: subname

        TEST1(test_zero_bits, proc)
        TEST1(test_max_bits, proc)
        TEST1(test_nlast_restriction, proc)
        TEST1(test_max_bits2, proc)

    end subroutine

    subroutine test_zero_bits(proc)
        ! Check that bits are counted correctly
        procedure(count_bits_t), pointer :: proc
        integer :: bits
        bits = proc(int([1234, 5678], int64), 1)
        call assert_equals(bits, 12)
    end subroutine

    subroutine test_max_bits(proc)
        ! If we set a maximum number of bits, it should specify the maximum
        ! number that the routine _must_ count.
        ! --> i.e. if it is more efficient to just count them, then it does.
        procedure(count_bits_t), pointer :: proc
        integer :: bits
        bits = proc(int([1234, 5678], int64), 1, 4)
        call assert_true(bits >= 4)
    end subroutine

    subroutine test_max_bits2(proc)
        ! If we set a maximum number of bits, then it should have no impact
        ! if the actual number of bits is lower.
        procedure(count_bits_t), pointer :: proc
        integer :: bits
        bits = proc(int([1234, 5678], int64), 1, 12345)
        call assert_equals(bits, 12)
    end subroutine

    subroutine test_nlast_restriction(proc)
        ! If we specify a final integer to look at in the array, check that
        ! further integers are ignored.
        ! n.b. defined such that array is from 0:nLast
        procedure(count_bits_t), pointer :: proc
        integer :: bits
        bits = proc(int([1234, 5678, 9101], int64), 1)
        call assert_equals(bits, 12)
    end subroutine
end
