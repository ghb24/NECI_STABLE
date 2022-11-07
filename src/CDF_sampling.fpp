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
    use constants, only: int32, int64, sp, dp
    use util_mod, only: binary_search_first_ge, stop_all, cumsum, near_zero
    use dSFMT_interface, only: genrand_real2_dSFMT, dSFMT_init
    better_implicit_none
    private
    public :: CDF_Sampler_t

    type :: CDF_Sampler_t
        private
        real(dp), allocatable :: p(:)
            !! The probabilities
        real(dp), allocatable :: cum_p(:)
            !! The cumulated probabilities.
    contains
        procedure, public :: sample
        procedure, public :: get_p
    end type

    interface CDF_Sampler_t
        module procedure construct_CDF_sampler_t
    end interface

contains

    pure function construct_CDF_sampler_t(w) result(res)
        !! Construct a CDF sampler from given weights.
        !!
        !! The weights do not have to be normalized.
        real(dp), intent(in) :: w(:)
        type(CDF_Sampler_t) :: res
        debug_function_name("construct_CDF_sampler_t")
        real(dp), allocatable :: p(:)
        ASSERT(all(w >= 0._dp) .and. .not. near_zero(sum(w)))
        p = w / sum(w)
        res = CDF_Sampler_t(p, cumsum(p))
    end function

    real(dp) elemental function get_p(this, val)
        !! Get the probability of a given value `val`.
        class(CDF_Sampler_t), intent(in) :: this
        integer, intent(in) :: val
        get_p = this%p(val)
    end function

    subroutine sample(this, val, p)
        !! Return randomly a value `val` and its probability `p`.
        class(CDF_Sampler_t), intent(in) :: this
        integer, intent(out) :: val
        real(dp), intent(out) :: p
        real(dp) :: r
        r = genrand_real2_dSFMT()
        val = binary_search_first_ge(this%cum_p, r)
        p = this%p(val)
    end subroutine

end module
