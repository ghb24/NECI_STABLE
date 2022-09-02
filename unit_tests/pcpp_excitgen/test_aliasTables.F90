! Do unit tests for the aliasSampling module
program test_aliasTables
    use constants, only: dp, EPS
    use fruit, only: init_fruit, fruit_summary, fruit_finalize, &
                     get_failed_count, run_test_case, assert_true
    use util_mod, only: stop_all
    use aliasSampling, only: aliasSampler_t, aliasTable_t
    use Parallel_neci, only: MPIInit, MPIEnd
    use dSFMT_interface, only: genrand_real2_dSFMT, dSFMT_init
    implicit none

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

        diff = get_diff(hist, huge_number, w)
    end function genrand_aliasTable

    !------------------------------------------------------------------------------------------!

    function get_diff(hist, nSamples, w) result(diff)
        integer, intent(in) :: hist(:), nSamples
        real(dp), intent(in) :: w(:)
        real(dp) :: diff
        integer :: i
        ! compare the drawn numbers with the distribtuion
        diff = 0
        if (size(w) /= size(hist)) then
            print *, "Error: size mismatch in getDiff"
            return
        end if
        do i = 1, size(hist)
            diff = diff + abs(hist(i) / real(nSamples) - w(i))
        end do

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
        type(aliasSampler_t) :: sampler
        integer, parameter :: huge_number = 1000000
        integer, parameter :: tSize = 1000
        real(dp), parameter :: diff_tol = 3e-2
        real(dp) :: w(tSize), p, probs(tSize)
        integer :: hist(tSize)
        integer :: r
        integer :: i
        real(dp) :: diff

        call create_rand_probs(w)
        call sampler%setupSampler(w)

        hist = 0
        do i = 1, huge_number
            call sampler%sample(r, p)
            hist(r) = hist(r) + 1
            probs(r) = p
        end do

        diff = get_diff(hist, huge_number, w)
        ! is the sampling reasonable?
        call assert_true(diff < diff_tol)
        ! are the probabilities claimed by the sampler the ones we put in?
        call assert_true(sum(abs(probs - w)) < eps)
    end subroutine test_aliasSampler

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
