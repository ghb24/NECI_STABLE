module en2molcas
   implicit NONE
   save
   real(8) :: NECI_E
end module en2molcas

subroutine NECImain(NECIen)
  use en2molcas, only : NECI_E
  character(64) :: dummy1, dummy2
  real(8), intent (out) :: NECIen
  write(6,*) "STARTING NECI from Molcas"
  dummy1 = ' '
  dummy2 = ' '
 ! Indicate not called by CPMD, VASP, Molpro
  call NECICore(0, .false., .false., .false., .true., dummy1, dummy2)
 ! Once we got excited states energies we will add them here to ENER array.
!        write(6,*) 'NECI_E in necimain', NECI_E
  NECIen = NECI_E
end subroutine NECImain
