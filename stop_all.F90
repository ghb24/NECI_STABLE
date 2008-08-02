subroutine stop_all(sub_name,error_msg)
! Stop calculation due to an error.
! Exit with code 999.
!
! In:
!    sub_name:  calling subroutine name.
!    error_msg: error message.

#ifdef PARALLEL
use Parallel, only: iProcIndex,MPIStopAll
#endif
implicit none
character(*) :: sub_name,error_msg

! It seems that giving STOP a string is far more portable.
! MPI_Abort requires an integer though.
integer, parameter :: error_code=999
character(3), parameter :: error_str='999'

write (6,'(a7)') 'ERROR.'
write (6,'(a27,a)') 'NECI stops in subroutine: ',trim(sub_name)
write (6,'(a9,18X,a)') 'Reason: ',trim(error_msg)
write (6,'(a11)') 'EXITING...'

#ifdef PARALLEL
write (6,'(a12,15X,i3)') 'Processor: ',iProcIndex
call MPIStopAll(error_str)
#else
stop error_str
#endif

return
end subroutine stop_all
