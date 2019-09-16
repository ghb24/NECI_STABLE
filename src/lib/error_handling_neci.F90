subroutine stop_all_c (sub_name, error_msg) bind(c)
    use iso_c_hack
    use util_mod, only: strlen_wrap
    implicit none

    character(c_char), target, intent(in) :: sub_name(*), error_msg(*)
    character(len=strlen_wrap(sub_name)), target :: sub_name_tmp
    character(len=strlen_wrap(error_msg)), target :: error_msg_tmp

    ! Convert from C character to standard fortran character string.
    ! Note that strlen does not include the null character at the end of the
    ! C string.  This is the behaviour we want.
    sub_name_tmp = transfer(sub_name(:strlen_wrap(sub_name)), sub_name_tmp)
    error_msg_tmp = transfer(error_msg(:strlen_wrap(error_msg)), error_msg_tmp)

    call stop_all (sub_name_tmp, error_msg_tmp)

end subroutine

subroutine stop_all (sub_name, error_msg)

    ! Stop calculation due to an error. Exit with code 999?
    !
    ! In: sub_name    - Calling routine
    !     error_msg   - Error message

#ifdef PARALLEL
    use Parallel_neci, only: iProcIndex, MPIStopAll
#endif
    implicit none

    interface
        subroutine print_backtrace_neci () bind(c)
        end subroutine
    end interface

    character(*), intent(in) :: sub_name, error_msg

    ! It seems that giving STOP a string is far more portable.
    ! MPI_Abort requires an integer though.
    character(3), parameter :: error_str='999'

#ifdef __DEBUG
    write (6,'(/a7)') 'ERROR.'
    write (6,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write (6,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef PARALLEL
    write (6,'(a12,15X,i5)') 'Processor: ',iProcIndex
#endif
    write (6,'(a11)') 'EXITING...'
    call neci_flush (6)
#else
    write (iout,'(/a7)') 'ERROR.'
    write (iout,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write (iout,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef PARALLEL
    write (iout,'(a12,15X,i5)') 'Processor: ',iProcIndex
#endif
    write (iout,'(a11)') 'EXITING...'
#endif

    ! Also push this to the stderr unit, so it hopefully ends up somewhere
    ! more useful.
    write (7,'(/a7)') 'ERROR.'
    write (7,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write (7,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef PARALLEL
    write (7,'(a12,15X,i5)') 'Processor: ',iProcIndex
#endif
    write (7,'(a11)') 'EXITING...'

    call print_backtrace_neci()

#ifdef PARALLEL
    call MPIStopAll(error_str)
#else
!    stop error_str
    stop
#endif

end subroutine stop_all



subroutine warning_neci(sub_name,error_msg)
!= Print a warning message in a (helpfully consistent) format.
!= I was bored of typing the same formatting in different places. ;-)
!=
!= In:
!=    sub_name:  calling subroutine name.
!=    error_msg: error message.

implicit none
character(*), intent(in) :: sub_name,error_msg

#ifdef __DEBUG
write (6,'(/a)') 'WARNING.  Error in '//adjustl(sub_name)
write (6,'(a/)') adjustl(error_msg)
#else
write (iout,'(/a)') 'WARNING.  Error in '//adjustl(sub_name)
write (iout,'(a/)') adjustl(error_msg)
#endif

return
end subroutine warning_neci



subroutine quiet_stop(msg)
!= Exit without making any noise.  Useful for when there's no error, but you
!= still want to exit midway through a calculation (e.g. for testing purposes,
!= or for use with the SOFTEXIT functionality).
!= In:
!=    msg (optional) : Print msg before exiting if msg is present.
#ifdef PARALLEL
use Parallel_neci, only: MPIStopAll
#endif

implicit none
character(*), intent(in), optional :: msg

if (present(msg)) then
    write (6,'(1X,a)') adjustl(msg)
    CALL neci_flush(6)
end if

#ifdef PARALLEL
call MPIStopAll(msg)
#else
stop
#endif

end subroutine quiet_stop
