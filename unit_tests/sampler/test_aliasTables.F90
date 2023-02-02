! Do unit tests for the aliasSampling module
#include "macros.h"
program test_aliasTables
    use constants, only: dp, stdout
    use fruit, only: assert_true, init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count
    use aliasSampling, only: aliasTable_t, aliasSampler_t
    use Parallel_neci, only: MPIInit, MPIEnd
    use dSFMT_interface, only: genrand_real2_dSFMT, dSFMT_init
    use util_mod, only: stop_all, near_zero, operator(.isclose.)
    better_implicit_none

    call init_fruit()
    call aliasSampling_test_driver()
    call fruit_summary()
    call fruit_finalize()
    block
        integer :: failed_count
        call get_failed_count(failed_count)
        if (failed_count /= 0) call stop_all('test_aliasTable', 'failed_tests')
    end block

contains

    !------------------------------------------------------------------------------------------!

    subroutine aliasSampling_test_driver()
        call MPIInit(.false.)
        ! initialize the RNG
        call dSFMT_init(8)
        ! do two test series: one for the aliasTable class
        call test_aliasTable()
        ! one for the aliasSampler class - this is a wrapper for the aliasTable class that
        ! also knows the original probability distribtuion
        call test_aliasSampler()

        call test_CDF_Sampler()
        call MPIEnd(.false.)
    end subroutine aliasSampling_test_driver

    !------------------------------------------------------------------------------------------!

    function genrand_aliasTable(huge_number) result(diff)
        integer, intent(in) :: huge_number
        real(dp) :: diff
        type(aliasTable_t) :: test_table
        integer, parameter :: tSize = 1000
        real(dp) :: w(tSize)
        integer :: hist(tSize)
        integer :: i, r

        ! use random probabilities to seed the aliasTable
        call create_rand_probs(w)
        ! init the table
        call test_table%setupTable(w)
        ! draw a huge number of random numbers from the table
        hist = 0
        do i = 1, huge_number
            r = test_table%getRand()
            hist(r) = hist(r) + 1
        end do

        diff = get_diff(hist, w)
    end function genrand_aliasTable

    !------------------------------------------------------------------------------------------!

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

    !------------------------------------------------------------------------------------------!

    subroutine test_aliasTable()
        integer :: nSamples
        real(dp), parameter :: diff_tol = 1e-2
        integer, parameter :: nRuns = 4
        real(dp) :: diff(nRuns)
        integer :: i

        ! check the total sampling error for different numbers of samples
        nSamples = 1000000
        do i = 1, nRuns
            diff(i) = genrand_aliasTable(nSamples)
            nSamples = 2 * nSamples
        end do
        ! it should only ever go down if we increase the number of samples
        do i = 1, nRuns - 1
            call assert_true(diff(i) > diff(i + 1))
        end do
        ! and we expect some threshold to be beaten
        call assert_true(diff(nRuns) < diff_tol)

!    print *, "Diff with 1e9 samples", genrand_aliasTable(1000000000) ( this is 7.2e-4 )
    end subroutine test_aliasTable

    !------------------------------------------------------------------------------------------!

    subroutine test_aliasSampler()
        ! we already know that the sampling is correct, just check that the sampler has the
        ! right probabilities
        use timing_neci, only: timer, set_timer, halt_timer, get_total_time

        type(aliasSampler_t) :: sampler
        integer, parameter :: huge_number = 100000000
        integer, parameter :: tSize = 400
        real(dp), parameter :: diff_tol = 3e-3_dp
        real(dp) :: w(tSize), p, probs(tSize)
        integer :: hist(tSize)
        type(timer) :: const_sample_timer
        integer :: r
        integer :: i
        real(dp) :: diff

        call create_rand_probs(w)
        call sampler%setupSampler(w)

        call assert_true(all(w.isclose.sampler%getProb([(i, i=1, tSize)])))

        hist = 0
        call set_timer(const_sample_timer)
        do i = 1, huge_number
            call sampler%sample(r, p)
            hist(r) = hist(r) + 1
            probs(r) = p
        end do
        call halt_timer(const_sample_timer)
        write(stdout, *) 'Full Alias sample', get_total_time(const_sample_timer)

        diff = get_diff(hist, w)
        ! is the sampling reasonable?

        call assert_true(diff < diff_tol)
        ! are the probabilities claimed by the sampler the ones we put in?
        call assert_true(near_zero(sum(abs(probs - w))))

    end subroutine test_aliasSampler

    subroutine test_CDF_Sampler()
        use timing_neci, only: timer, set_timer, halt_timer, get_total_time
        use CDF_sampling_mod, only: CDF_Sampler_t

        type(aliasSampler_t) :: alias_sampler

        integer, parameter :: huge_number = 10000000
        integer, parameter :: tSize = 10
        type(timer) :: full_sampler, const_sample_timer
        real(dp), parameter :: diff_tol = 3e-3_dp
        real(dp) :: renorm
        integer, allocatable :: contain(:), exclude(:)
        real(dp) :: w(tSize), p, probs(tSize)
        integer :: hist(tSize)
        integer :: r
        integer :: i
        real(dp) :: diff

        call create_rand_probs(w)
        contain = [(i, i=1, tSize / 2)]
        exclude = [integer::]

        call alias_sampler%setupSampler(w)
        call assert_true(all(w.isclose.alias_sampler%getProb([(i, i=1, tSize)])))

        hist = 0
        call set_timer(full_sampler)
        do i = 1, huge_number
            call alias_sampler%sample(r, p)
            hist(r) = hist(r) + 1
            probs(r) = p
        end do
        call halt_timer(full_sampler)
        write(stdout, *) 'Full Alias sample', get_total_time(full_sampler)

        diff = get_diff(hist(contain), w(contain) / sum(w(contain)))

        ! is the sampling reasonable?
        call assert_true(diff < diff_tol)
        ! are the probabilities claimed by the sampler the ones we put in?

        hist = 0
        call set_timer(const_sample_timer)
        do i = 1, huge_number
            renorm = sum(w(contain))
            block
                integer :: pos
                call alias_sampler%constrained_sample(contain, renorm, pos, r, p)
                call assert_true(contain(pos) == r)
            end block
            hist(r) = hist(r) + 1
            probs(r) = p
        end do
        call halt_timer(const_sample_timer)
        write(stdout, *) 'constrained Alias sample', get_total_time(const_sample_timer)

        diff = get_diff(hist(contain), w(contain) / sum(w(contain)))

        ! is the sampling reasonable?
        call assert_true(diff < diff_tol)
        ! are the probabilities claimed by the sampler the ones we put in?

        call alias_sampler%setupSampler(w)
        call assert_true(all(w.isclose.alias_sampler%getProb([(i, i=1, tSize)])))

    end subroutine

    !------------------------------------------------------------------------------------------!

    subroutine create_rand_probs(w)
        ! fill an array with a normalized distribution of uniform random numbers
        ! Input: w - array of reals, on return filled with random numbers
        real(dp), intent(out) :: w(:)
        integer :: i
        ! create some random probabilities
        do i = 1, size(w)
            w(i) = genrand_real2_dSFMT()
        end do
        ! normalize the probs
        w = w / sum(w)
    end subroutine create_rand_probs

end program test_aliasTables
