module vasp_neci_interface

use constants, only: dp

contains
   subroutine construct_ijab_one(i,j,a,b,integral)
   implicit none
   integer :: i,j,a,b
   complex(dp) :: integral(1,1)
   end subroutine construct_ijab_one
end module vasp_neci_interface
