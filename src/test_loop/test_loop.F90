program neci_loop_test
  use mpi
  use constants
  implicit none
  integer :: i
  integer, parameter :: nLoop = 3
  integer(MPIArg) :: err

  ! run neci nLoop times - this reveals initialization/finalization issues
  do i = 1, nLoop
     call NECICore(0, .false., .false., .false., .true.)
  end do
  ! we run in external mode to mimick the behaviour of embedded neci
  ! this requires manual finalization of MPI
  call MPI_Finalize(err)

end program neci_loop_test
