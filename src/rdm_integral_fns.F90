! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
module rdm_integral_fns

    use constants

    implicit none

contains

    pure function one_elec_int(i, j) result(integral)

        use OneEInts, only: TMat2D

        integer, intent(in) :: i, j
        real(dp) :: integral

        integral = TMat2D(i,j)

    end function one_elec_int

    function two_elec_int(i, j, k, l) result(integral)

        use sltcnd_mod, only: sltcnd_2

        integer, intent(in) :: i, j, k, l
        integer :: ex(2,2)
        real(dp) :: integral

        ex(1,1) = j
        ex(1,2) = i
        ex(2,1) = l
        ex(2,2) = k

        integral = real(sltcnd_2(ex, .false.), dp)

    end function two_elec_int

end module rdm_integral_fns
