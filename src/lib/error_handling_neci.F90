subroutine pp_stop_all(sub_name, error_msg, file_name, line_number)
    ! additional pre-processor info (pp) of file and line-number in
    ! stop_all call
    implicit none

    character(*), intent(in) :: sub_name, error_msg, file_name
    integer, intent(in) :: line_number


    ! print to stdout
    write (6,'(/a10,17X,a)')  ' In file: ', adjustl(file_name)
    write (6,'(a10,17X,i6)') ' At line: ', line_number

    ! print to stderr
    write (7,'(/a10,17X,a)')  ' In file: ', adjustl(file_name)
    write (7,'(a10,17X,i6)') ' At line: ', line_number

    call stop_all(sub_name, error_msg)

end subroutine pp_stop_all


subroutine stop_all(sub_name, error_msg)

    ! Stop calculation due to an error. Exit with code 222?
    !
    ! In: sub_name    - Calling routine
    !     error_msg   - Error message

#ifdef USE_MPI
    use Parallel_neci, only: iProcIndex, MPIStopAll
#endif
    implicit none

    interface
        subroutine print_backtrace_neci () bind(c)
        end subroutine
    end interface

    character(*), intent(in) :: sub_name, error_msg

    ! I found problems when the error code is larger 2^8 - 1 == 255
    integer, parameter :: error_code = 222

    write (6,'(/a7)') 'ERROR.'
    write (6,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write (6,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef USE_MPI
    write (6,'(a12,15X,i5)') 'Processor: ',iProcIndex
#endif
    write (6,'(a11)') 'EXITING...'
#ifdef DEBUG_
    call neci_flush (6)
#else
    write (6,'(/a7)') 'ERROR.'
    write (6,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write (6,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef USE_MPI
    write (6,'(a12,15X,i5)') 'Processor: ',iProcIndex
#endif
    write (6,'(a11)') 'EXITING...'
#endif

    ! Also push this to the stderr unit, so it hopefully ends up somewhere
    ! more useful.
    write (7,'(/a7)') 'ERROR.'
    write (7,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write (7,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef USE_MPI
    write (7,'(a12,15X,i5)') 'Processor: ',iProcIndex
#endif
    write (7,'(a11)') 'EXITING...'

    call print_backtrace_neci()

#ifdef USE_MPI
    call MPIStopAll(error_code)
#else
    stop error_code
#endif

end subroutine stop_all



subroutine quiet_stop(msg)
!= Exit without making any noise.  Useful for when there's no error, but you
!= still want to exit midway through a calculation (e.g. for testing purposes,
!= or for use with the SOFTEXIT functionality).
!= In:
!=    msg (optional) : Print msg before exiting if msg is present.
#ifdef USE_MPI
use Parallel_neci, only: MPIStopAll
#endif

    implicit none
    character(*), intent(in), optional :: msg

    if (present(msg)) then
        write (6,'(1X,a)') adjustl(msg)
        CALL neci_flush(6)
    end if

#ifdef USE_MPI
call MPIStopAll(0)
#else
    stop
#endif

end subroutine quiet_stop


! subroutine stop_all_c (sub_name, error_msg) bind(c)
!     use iso_c_hack
!     use util_mod, only: strlen_wrap
!     implicit none
!
!     character(c_char), target, intent(in) :: sub_name(*), error_msg(*)
!     character(len=strlen_wrap(sub_name)), target :: sub_name_tmp
!     character(len=strlen_wrap(error_msg)), target :: error_msg_tmp
!
!     ! Convert from C character to standard fortran character string.
!     ! Note that strlen does not include the null character at the end of the
!     ! C string.  This is the behaviour we want.
!     sub_name_tmp = transfer(sub_name(:strlen_wrap(sub_name)), sub_name_tmp)
!     error_msg_tmp = transfer(error_msg(:strlen_wrap(error_msg)), error_msg_tmp)
!
!     call stop_all (sub_name_tmp, error_msg_tmp)
!
! end subroutine


