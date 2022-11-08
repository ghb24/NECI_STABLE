#include "macros.h"
module CDF_sampler_test_mod
    use constants, only: dp
    use util_mod, only: operator(.isclose.), stop_all
    use fruit, only: assert_true
    use dSFMT_interface, only: genrand_real2_dSFMT, dSFMT_init
    use CDF_sampling_mod, only: CDF_Sampler_t
    better_implicit_none
    private
    public :: CDF_sampler_test_driver

contains

    !------------------------------------------------------------------------------------------!

    subroutine CDF_sampler_test_driver()
        call test_CDF_Sampler()
    end subroutine

    !------------------------------------------------------------------------------------------!


    subroutine test_CDF_Sampler()
        type(CDF_Sampler_t) :: sampler
        integer :: i

        check_empty_set: block
            real(dp) :: p
            integer :: val
            sampler = CDF_Sampler_t([real(dp)::])

            do i = 1, 10
                call sampler%sample(val, p)
                call assert_true(val == 0)
                call assert_true(p .isclose. 1.0_dp)
            end do
        end block check_empty_set

        check_actual_sampling: block
            integer, parameter :: sample_iterations(4) = [100, 1000, 10000, 1000000], number_probs = 10
            real(dp) :: diffs(size(sample_iterations))
            real(dp), parameter :: diff_tol = 3e-3_dp
            real(dp), allocatable :: probs(:)

            probs = generate_probs(number_probs)
            sampler = CDF_Sampler_t(probs)
            call assert_true(all(probs .isclose. sampler%get_prob([(i, i = 1, sampler%size())])))


            do i = 1, size(sample_iterations)
                diffs(i) = get_diff(get_hist(sampler, sample_iterations(i)), probs)
            end do

            ! ensure that the error goes monotonically down, when increasing sample size
            call assert_true(all(diffs - eoshift(diffs, +1) > 0))
            ! ensure that the last difference is below our error threshold
            call assert_true(diffs(size(sample_iterations)) < diff_tol)
        end block check_actual_sampling

    end subroutine


    function generate_probs(N) result(p)
        integer, intent(in) :: N
        real(dp), allocatable :: p(:)
        integer :: i
        allocate(p(N))
        do i = 1, size(p)
            p(i) = genrand_real2_dSFMT()
        end do
        p = p / sum(p)
    end function

    function get_hist(sampler, N) result(hist)
        type(CDF_Sampler_t), intent(in) :: sampler
        integer, intent(in) :: N
        integer :: hist(sampler%size())
        integer :: i, val
        real(dp) :: p
        hist(:) = 0
        do i = 1, N
            call sampler%sample(val, p)
            hist(val) = hist(val) + 1
        end do
    end function


    pure function get_diff(hist, w) result(diff)
        integer, intent(in) :: hist(:)
        real(dp), intent(in) :: w(:)
        integer :: nSamples
        character(*), parameter :: this_routine = 'get_diff'
        real(dp) :: diff
        ! compare the drawn numbers with the distribtuion
        if (size(w) /= size(hist)) then
            call stop_all(this_routine, "Error: size mismatch in getDiff")
        end if
        nSamples = sum(hist)
        diff = sum(abs(hist / real(nSamples, dp) - w))
    end function get_diff


end module

program test_CDF_sampler
    use fruit, only: fruit_summary, fruit_finalize, get_failed_count, init_fruit
    use Parallel_neci, only: MPIInit, MPIEnd
    use util_mod, only: stop_all, near_zero, operator(.isclose.)
    use dSFMT_interface, only: dSFMT_init

    use CDF_sampler_test_mod, only: CDF_sampler_test_driver
    better_implicit_none

    call MPIInit(.false.)
    call dSFMT_init(8)
    call init_fruit()
    call CDF_sampler_test_driver()
    call fruit_summary()
    block
        integer :: failed_count
        call get_failed_count(failed_count)
        if (failed_count /= 0) call stop_all('test_CDF_sampler', 'failed_tests')
    end block
    call fruit_finalize()
    call MPIEnd(.false.)

end program
