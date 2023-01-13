#include "macros.h"
#:include "macros.fpph"
module CDF_sampling_mod
    !! This module implements the sampling of non-uniform probability
    !! distributions that are constructed on the fly
    !! using the cumulated distribution function method.
    !!
    !! This sampling is faster to construct than the alias-sampling,
    !! but it is slower to sample.
    !! For distributions that are known to be constant during
    !! a calculation use the alias sampling.
    !! This is done for example in PCHB.
    use constants, only: dp
    use util_mod, only: binary_search_first_ge, stop_all, cumsum, near_zero, stop_all, operator(.isclose.)
    use dSFMT_interface, only: genrand_real2_dSFMT, dSFMT_init
    better_implicit_none
    private
    public :: CDF_Sampler_t

    type :: CDF_Sampler_t
        private
        real(dp), allocatable :: p(:)
            !! The probabilities.
        real(dp), allocatable :: cum_p(:)
            !! The cumulated probabilities.
        integer :: my_size
            !! We store the size in an additional integer,
            !! because `p` and `cum_p` are allocated empty,
            !! if all probabilities are zero.
    contains
        procedure, public :: sample
        procedure, public :: get_prob
        procedure, public :: all_zero
        procedure, public :: size => get_size
    end type

    interface CDF_Sampler_t
        module procedure construct_CDF_sampler_t, construct_CDF_sampler_with_total_t
    end interface

contains

    function construct_CDF_sampler_t(w) result(res)
        !! Construct a CDF sampler from given weights.
        !!
        !! The weights do not have to be normalized.
        !! If all weights are zero, or the array is empty,
        !! the sampler will return the (non-existent) index 0 with probability 1.
        real(dp), intent(in) :: w(:)
        type(CDF_Sampler_t) :: res
        res = CDF_Sampler_t(w, sum(w))
    end function


    function construct_CDF_sampler_with_total_t(w, total) result(res)
        !! Construct a CDF sampler from given weights.
        !!
        !! The weights do not have to be normalized.
        !! If all weights are zero, or the array is empty,
        !! the sampler will return the (non-existent) index 0 with probability 1.
        real(dp), intent(in) :: w(:), total
        type(CDF_Sampler_t) :: res
        debug_function_name("construct_CDF_sampler_with_total_t")
        res%my_size = size(w)

        @:ASSERT(all(w >= 0._dp))
        @:ASSERT(sum(w) .isclose. total)
        if (near_zero(total)) then
            ! allocate as empty sets
            allocate(res%p(0), res%cum_p(0))
        else
            @:ASSERT(sum(w) .isclose. total, sum(w), total, w)
            res%p = w(:) / total
            res%cum_p = cumsum(res%p)
            @:ASSERT(res%cum_p(size(res%cum_p)) .isclose. 1._dp)
        end if
    end function

    logical elemental function all_zero(this)
        !! Return if all probabilities are zero, or the set of probabilities is empty.
        class(CDF_Sampler_t), intent(in) :: this
        all_zero = size(this%p) == 0
    end function

    real(dp) elemental function get_prob(this, val)
        !! Get the probability of a given value `val`.
        class(CDF_Sampler_t), intent(in) :: this
        integer, intent(in) :: val
        debug_function_name("get_p")
        @:pure_ASSERT(1 <= val .and. val <= size(this%p))
        if (this%all_zero()) then
            get_prob = 0.0
        else
            get_prob = this%p(val)
        end if
    end function

    subroutine sample(this, val, p)
        !! Return randomly a value `val` and its probability `p`.
        class(CDF_Sampler_t), intent(in) :: this
        integer, intent(out) :: val
        real(dp), intent(out) :: p
        real(dp) :: r
        if (this%all_zero()) then
            val = 0
            p = 1.0_dp
        else
            r = genrand_real2_dSFMT()
            val = binary_search_first_ge(this%cum_p, r)
            p = this%p(val)
        end if
    end subroutine

    elemental integer function get_size(this)
        !! Return the number of probabilites.
        class(CDF_Sampler_t), intent(in) :: this
        get_size = this%my_size
    end function

end module
