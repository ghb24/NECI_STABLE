#include "macros.h"
#:include "macros.fpph"

module CDF_sampling_mod
    !! This module implements the sampling from nonuniform probality
    !!  distributions, via the Cumulative Distribution functions.
    !!
    !! The algorithm works for drawing from distributions
    !! that are constructed on the fly.
    !! For probality distributions that are fully known
    !! at initialization time of the calculation it
    !! is better to use alias sampling as it is done for e.g. PCHB.
    use constants, only: int32, int64, sp, dp
    use util_mod, only: binary_search_first_ge, stop_all, cumsum
    use dSFMT_interface, only: genrand_real2_dSFMT, dSFMT_init
    better_implicit_none
    private
    public :: CDF_Sampler_t

    type :: CDF_Sampler_t
        private
        real(dp), allocatable :: p(:)
        real(dp), allocatable :: cum_p(:)
    contains
        procedure, public :: sample
        procedure, public :: get_p
    end type

    interface CDF_Sampler_t
        module procedure construct_CDF_sampler_t
    end interface

contains

    pure function construct_CDF_sampler_t(w) result(res)
        !! Initialize the CDF sampler with given weights.
        !!
        !! The weights `w` do not have to be normalized.
        real(dp), intent(in) :: w(:)
        type(CDF_Sampler_t) :: res
        real(dp), allocatable :: p(:)
        p = w / sum(w)
        res = CDF_Sampler_t(p, cumsum(p))
    end function

    real(dp) elemental function get_p(this, val)
        !! Return the probability of drawing `val`
        class(CDF_Sampler_t), intent(in) :: this
        integer, intent(in) :: val
        get_p = this%p(val)
    end function

    subroutine sample(this, val, p)
        !! Return a random value `val` and its probality `p`.
        class(CDF_Sampler_t), intent(in) :: this
        integer, intent(out) :: val
        real(dp), intent(out) :: p
        real(dp) :: r
        r = genrand_real2_dSFMT()
        val = binary_search_first_ge(this%cum_p, r)
        p = this%p(val)
    end subroutine

end module
