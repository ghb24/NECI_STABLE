#include "macros.h"
module rdm_integral_fns

    use constants
    use excitation_types, only: Excite_2_t
    use SystemData, only: t_heisenberg_model, exchange_j 
    use lattice_mod, only: lat

    implicit none

contains

    pure function one_elec_int(i, j) result(integral)

        use OneEInts, only: GetTMatEl

        integer, intent(in) :: i, j
        real(dp) :: integral

        if(t_heisenberg_model)then
            integral = 0.0_dp
        else            
            integral = GetTMatEl(i, j)
        endif

    end function one_elec_int

    function two_elec_int_heisenberg(i, j, k, l) result(integral)

        integer, intent(in) :: i, j, k, l
        real(dp) :: integral

        integer :: x, y, z, w !site index
        logical :: a, b, c, d !spin
        integer, allocatable :: neighbors(:)
        logical :: are_neighbors
        integer :: idx

        !Spin-Orb Excitation: i->k, j-> l
        !Site Excitation: x->z, y-> w

        x = get_spatial(i)
        y = get_spatial(j)
        z = get_spatial(k)
        w = get_spatial(l)

        if((x/=z) .or. (y/=w))then
            ! exciting to  a different site
            integral = 0.0_dp
        else
            neighbors = lat%get_neighbors(x)
            are_neighbors = .false.
            do idx=1, size(neighbors)
                if(y == neighbors(idx))then
                    are_neighbors = .true.
                    exit
                endif
            enddo

            if(.not. are_neighbors)then
                integral = 0.0_dp
            else
                !Spin Excitation: a->c, b-> d
                a = is_beta(i)
                b = is_beta(j)
                c = is_beta(k)
                d = is_beta(l)

                if((a.eqv.c) .and. (b.eqv.d))then
                    !both spins stayed the same

                    if (a.eqv.b)then
                        integral = exchange_j/4.0_dp 
                    else 
                        integral = -exchange_j/4.0_dp
                    endif

                elseif((a.neqv.c) .and. (b.neqv.d))then

                    !both spins flipped
                    integral = exchange_j/2.0_dp 
                else
                    !one spin flipped but not the other
                    integral = 0.0_dp
                endif

            endif
            deallocate(neighbors)
        endif

    end function two_elec_int_heisenberg

    function two_elec_int(i, j, k, l) result(integral)

        use sltcnd_mod, only: sltcnd_2_kernel
        use UMatCache, only: gtID

        integer, intent(in) :: i, j, k, l
        real(dp) :: integral

        if(t_heisenberg_model)then
            integral = two_elec_int_heisenberg(i, j, k, l)
        else
            integral = real(sltcnd_2_kernel(Excite_2_t(j, l, i, k)), dp)
        endif

    end function two_elec_int

    function GetPropInts(i, j, iprop) result(integral)

        use OneEInts, only: OneEPropInts
        integer, intent(in) :: i, j, iprop
        real(dp) :: integral

        integral = OneEPropInts(i, j, iprop)

    end function GetPropInts

end module rdm_integral_fns
