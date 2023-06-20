! Do unit tests for the aliasSampling module
#include "macros.h"
program test_aliasTables
    use constants, only: dp
    use fruit, only: assert_true, init_fruit, fruit_summary, fruit_finalize, &
        get_failed_count
    use aliasSampling, only: AliasTable_t, AliasSampler_t
    use Parallel_neci, only: MPIInit, MPIEnd, MPIBarrier, iProcIndex_intra, root, MPIBcast, Node
    use dSFMT_interface, only: genrand_real2_dSFMT, dSFMT_init
    use util_mod, only: stop_all, near_zero, operator(.isclose.)
    better_implicit_none
    integer :: ierr

    call MPIInit(.false.)
    call init_fruit()
    call aliasSampling_test_driver()
    call MPIBarrier(ierr)
    call fruit_summary()
    call fruit_finalize()
    block
        integer :: failed_count
        call get_failed_count(failed_count)
        if (failed_count /= 0) call stop_all('test_aliasTable', 'failed_tests')
    end block
    call MPIEnd(.false.)

contains

    !------------------------------------------------------------------------------------------!

    subroutine aliasSampling_test_driver()
        ! initialize the RNG
        call dSFMT_init(8 + iProcIndex_intra)
        ! do two test series: one for the aliasTable class
        call test_aliasTable()
        !
        ! ! one for the aliasSampler class - this is a wrapper for the aliasTable class that
        ! ! also knows the original probability distribtuion
        call test_aliasSampler()
        !
        call test_if_works_from_one_proc()

        call test_CDF_Sampler()
    end subroutine aliasSampling_test_driver

    !------------------------------------------------------------------------------------------!

    function genrand_aliasTable(n_probs, n_sample) result(diff)
        integer, intent(in) :: n_probs, n_sample
        real(dp) :: diff
        type(AliasTable_t) :: test_table
        real(dp) :: w(n_probs)
        integer :: hist(n_probs)
        integer :: i, r

        ! use random probabilities to seed the aliasTable
        if (iProcIndex_intra == root) w = create_rand_probs(n_probs)
        call MPIBcast(w, iProcIndex_intra == root, Node)
        ! init the table
        call test_table%setup(root, w)
        ! draw a huge number of random numbers from the table
        hist = 0
        do i = 1, n_sample
            r = test_table%sample()
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
        integer :: n_samples
        real(dp), parameter :: diff_tol = 1e-2
        integer, parameter :: nRuns = 4, n_probs = 1000
        real(dp) :: diff(nRuns)
        integer :: i

        ! check the total sampling error for different numbers of samples
        n_samples = 1000000
        do i = 1, nRuns
            diff(i) = genrand_aliasTable(n_probs, n_samples)
            n_samples = 2 * n_samples
        end do
        ! it should only ever go down if we increase the number of samples
        do i = 1, nRuns - 1
            call assert_true(diff(i) > diff(i + 1))
        end do
        ! and we expect some threshold to be beaten
        call assert_true(diff(nRuns) < diff_tol)

!    print *, "Diff with 1e9 samples", genrand_aliasTable(1000000000) ( this is 7.2e-4 )
    end subroutine test_aliasTable


    subroutine test_if_works_from_one_proc()
        !! Here we test, if it is possible to
        !! initialize the sampler with just one weights array per node.
        !! The non-root procs can have an array of zero-length
        integer, parameter :: n_probs = 4, huge_number = 10**1
        type(AliasSampler_t) :: sampler
        real(dp), allocatable :: w_in(:)
        real(dp) :: probs(n_probs), diff, p, w(n_probs)
        integer :: i, r, hist(n_probs)

        ! use random probabilities to seed the aliasTable
        if (iProcIndex_intra == root) then
            w_in = create_rand_probs(n_probs)
        else
            w_in = [real(dp) :: ]
        end if

        ! init the table
        call sampler%setup(root, w_in)

        if (iProcIndex_intra == root) w = w_in
        call MPIBcast(w, iProcIndex_intra == root, Node)
        call assert_true(all(w .isclose. sampler%get_prob([(i, i=1, n_probs)])))

        hist = 0
        do i = 1, huge_number
            call sampler%sample(r, p)
            hist(r) = hist(r) + 1
            probs(r) = p
        end do

        diff = get_diff(hist, w_in)
    end subroutine

    !------------------------------------------------------------------------------------------!

    subroutine test_aliasSampler()
        ! we already know that the sampling is correct, just check that the sampler has the
        ! right probabilities
        use timing_neci, only: timer, set_timer, halt_timer, get_total_time

        type(AliasSampler_t) :: sampler
        integer, parameter :: huge_number = 100000000
        integer, parameter :: n_probs = 400
        real(dp), parameter :: diff_tol = 3e-3_dp
        real(dp) :: w(n_probs), p, probs(n_probs)
        integer :: hist(n_probs)
        type(timer) :: const_sample_timer
        integer :: r
        integer :: i
        real(dp) :: diff

        if (iProcIndex_intra == root) w = create_rand_probs(n_probs)
        call MPIBcast(w, iProcIndex_intra == root, Node)

        call sampler%setup(root, w)

        call assert_true(all(w .isclose. sampler%get_prob([(i, i=1, n_probs)])))

        hist = 0
        call set_timer(const_sample_timer)
        do i = 1, huge_number
            call sampler%sample(r, p)
            hist(r) = hist(r) + 1
            probs(r) = p
        end do
        call halt_timer(const_sample_timer)

        diff = get_diff(hist, w)
        ! is the sampling reasonable?

        call assert_true(diff < diff_tol)
        ! are the probabilities claimed by the sampler the ones we put in?
        call assert_true(near_zero(sum(abs(probs - w))))

    end subroutine test_aliasSampler

    subroutine test_CDF_Sampler()
        use timing_neci, only: timer, set_timer, halt_timer, get_total_time
        use CDF_sampling_mod, only: CDF_Sampler_t

        type(AliasSampler_t) :: alias_sampler

        integer, parameter :: n_samples = 10000000
        integer, parameter :: n_probs = 10
        type(timer) :: full_sampler, const_sample_timer
        real(dp), parameter :: diff_tol = 3e-3_dp
        real(dp) :: renorm
        integer, allocatable :: contain(:)
        real(dp) :: w(n_probs), p, probs(n_probs)
        integer :: hist(n_probs)
        integer :: r
        integer :: i
        real(dp) :: diff

        if (iProcIndex_intra == root) w = create_rand_probs(n_probs)
        call MPIBcast(w, iProcIndex_intra == root, Node)

        contain = [(i, i=1, n_probs / 2)]

        call alias_sampler%setup(root, w)
        call assert_true(all(w .isclose. alias_sampler%get_prob([(i, i=1, n_probs)])))

        hist = 0
        call set_timer(full_sampler)
        do i = 1, n_samples
            call alias_sampler%sample(r, p)
            hist(r) = hist(r) + 1
            probs(r) = p
        end do
        call halt_timer(full_sampler)

        diff = get_diff(hist(contain), w(contain) / sum(w(contain)))

        ! is the sampling reasonable?
        call assert_true(diff < diff_tol)

        ! are the probabilities claimed by the sampler the ones we put in?
        hist = 0
        call set_timer(const_sample_timer)
        do i = 1, n_samples
            renorm = sum(w(contain))
            block
                integer :: pos
                call alias_sampler%constrained_sample(contain, renorm, pos, r, p)
            end block
            hist(r) = hist(r) + 1
            probs(r) = p
        end do
        call halt_timer(const_sample_timer)

        diff = get_diff(hist(contain), w(contain) / sum(w(contain)))

        ! is the sampling reasonable?
        call assert_true(diff < diff_tol)
        ! are the probabilities claimed by the sampler the ones we put in?

        call alias_sampler%setup(root, w)
        call assert_true(all(w .isclose. alias_sampler%get_prob([(i, i=1, n_probs)])))

    end subroutine

    !------------------------------------------------------------------------------------------!

    function create_rand_probs(n_probs) result(w)
        ! fill an array with a normalized distribution of uniform random numbers
        ! Input: w - array of reals, on return filled with random numbers
        integer, intent(in) :: n_probs
        real(dp) :: w(n_probs)
        integer :: i
        ! create some random probabilities
        do i = 1, size(w)
            w(i) = genrand_real2_dSFMT()
        end do
        ! normalize the probs
        w = w / sum(w)
    end function

end program test_aliasTables
