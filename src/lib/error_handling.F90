subroutine stop_all(sub_name,error_msg)
!= Stop calculation due to an error.
!= Exit with code 999.
!=
!= In:
!=    sub_name:  calling subroutine name.
!=    error_msg: error message.

#ifdef PARALLEL
use Parallel, only: iProcIndex,MPIStopAll
#endif
implicit none
character(*), intent(in) :: sub_name,error_msg

! It seems that giving STOP a string is far more portable.
! MPI_Abort requires an integer though.
integer, parameter :: error_code=999
character(3), parameter :: error_str='999'

write (6,'(/a7)') 'ERROR.'
write (6,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
write (6,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
write (6,'(a11)') 'EXITING...'
CALL FLUSH(6)

#ifdef PARALLEL
write (6,'(a12,15X,i3)') 'Processor: ',iProcIndex
call MPIStopAll(error_str)
#else
stop error_str
#endif

return
end subroutine stop_all



subroutine warning(sub_name,error_msg)
!= Print a warning message in a (helpfully consistent) format.
!= I was bored of typing the same formatting in different places. ;-)
!=
!= In:
!=    sub_name:  calling subroutine name.
!=    error_msg: error message.

implicit none
character(*), intent(in) :: sub_name,error_msg

write (6,'(/a)') 'WARNING.  Error in '//adjustl(sub_name)
write (6,'(a/)') adjustl(error_msg)

return
end subroutine warning



subroutine quiet_stop(msg)
!= Exit without making any noise.  Useful for when there's no error, but you
!= still want to exit midway through a calculation (e.g. for testing purposes,
!= or for use with the SOFTEXIT functionality).
!= In:
!=    msg (optional) : Print msg before exiting if msg is present.
#ifdef PARALLEL
use Parallel, only: MPIStopAll
#endif

implicit none
character(*), intent(in), optional :: msg

if (present(msg)) then
    write (6,'(X,a)') adjustl(msg)
    CALL FLUSH(6)
end if

#ifdef PARALLEL
call MPIStopAll(msg)
#else
stop
#endif

end subroutine quiet_stop
