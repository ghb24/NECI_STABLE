#:set OrbIdxTypes = ['SpinOrbIdx_t', 'SpatOrbIdx_t']

module orb_idx_mod
    implicit none
    private
    public :: OrbIdx_t, SpinOrbIdx_t, SpatOrbIdx_t

    type, abstract :: OrbIdx_t
        integer, allocatable :: idx(:)
    end type

    !> We assume order [beta_1, alpha_1, beta_2, alpha_2, ...]
    type, extends(OrbIdx_t) :: SpinOrbIdx_t
    end type

    type, extends(OrbIdx_t) :: SpatOrbIdx_t
    end type

#:for orb_idx_type in OrbIdxTypes
    interface ${orb_idx_type}$
        module procedure construction_from_array_${orb_idx_type}$
    end interface
#:endfor

    contains


#:for orb_idx_type in OrbIdxTypes
    pure function construction_from_array_${orb_idx_type}$(idx) result(res)
        integer, intent(in) :: idx(:)
        type(${orb_idx_type}$) :: res
        res%idx = idx
    end function
#:endfor
end module
