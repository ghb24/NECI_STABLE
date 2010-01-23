module vasp_neci_interface

integer, parameter :: q=kind(0.d0)

contains
   subroutine construct_ijab_one(i,j,a,b,integral)
   implicit none
   integer :: i,j,a,b
   complex(q) :: integral(1,1)
   end subroutine construct_ijab_one
end module vasp_neci_interface
