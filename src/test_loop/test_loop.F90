program neci_loop_test
  implicit none
  integer :: i
  character(64) :: dummy1, dummy2
  integer, parameter :: nLoop = 3

  dummy1 = ' '
  dummy2 = ' '

  do i = 1, nLoop
     call NECICore(0, .false., .false., .false., .true., dummy1, dummy2)
  end do
  
end program neci_loop_test
