subroutine stop_all(sub_name, error_msg)

    ! Stop calculation due to an error. Exit with code 222?
    !
    ! In: sub_name    - Calling routine
    !     error_msg   - Error message
    use constants, only: stdout, stderr

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

    write (stdout,'(/a7)') 'ERROR.'
    write (stdout,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write (stdout,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef USE_MPI
    write (stdout,'(a12,15X,i5)') 'Processor: ',iProcIndex
#endif
    write (stdout,'(a11)') 'EXITING...'
#ifdef DEBUG_
    call neci_flush (stdout)
#else
    write (stdout,'(/a7)') 'ERROR.'
    write (stdout,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write (stdout,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef USE_MPI
    write (stdout,'(a12,15X,i5)') 'Processor: ',iProcIndex
#endif
    write (stdout,'(a11)') 'EXITING...'
#endif

    ! Also push this to the stderr unit, so it hopefully ends up somewhere
    ! more useful.
    write (stderr,'(/a7)') 'ERROR.'
    write (stderr,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write (stderr,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef USE_MPI
    write (stderr,'(a12,15X,i5)') 'Processor: ',iProcIndex
#endif
    write (stderr,'(a11)') 'EXITING...'

    call print_backtrace_neci()

#ifdef USE_MPI
    call MPIStopAll(error_code)
#else
    stop error_code
#endif

end subroutine stop_all





