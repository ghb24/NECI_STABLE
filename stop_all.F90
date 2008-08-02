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

integer, parameter :: error_code=999

write (6,'(a7)') 'ERROR.'
write (6,'(a27,a)') 'NECI stops in subroutine: ',trim(sub_name)
write (6,'(a9,18X,a)') 'Reason: ',trim(error_msg)
write (6,'(a11)') 'EXITING...'

#ifdef PARALLEL
write (6,'(a12,15X,i3)') 'Processor: ',iProcIndex
call MPIStopAll(error_code)
#else
stop 999
#endif

return
end subroutine stop_all
