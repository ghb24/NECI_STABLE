#include "macros.h"
submodule (error_handling_neci) error_handling_neci_impls
    use constants, only: stdout, stderr

#ifdef USE_MPI
    use Parallel_neci, only: iProcIndex, MPIStopAll
#endif
#ifdef NAGF95
    use f90_unix, only: flush
    use constants, only: int32
#endif
    better_implicit_none

contains


    subroutine neci_flush(un)
        integer, intent(in) :: un
#ifdef NAGF95
        integer(kind=int32) :: dummy
#endif
#ifdef BLUEGENE_HACKS
        call flush_(un)
#else
#ifdef NAGF95
        dummy = un
        call flush(dummy)
#else
        call flush(un)
#endif
#endif
    end subroutine neci_flush

    subroutine warning_neci(sub_name,error_msg)
        != Print a warning message in a (helpfully consistent) format.
        !=
        != In:
        !=    sub_name:  calling subroutine name.
        !=    error_msg: error message.
        character(*), intent(in) :: sub_name, error_msg

        write (stderr,'(/a)') 'WARNING.  Error in '//adjustl(sub_name)
        write (stderr,'(a/)') adjustl(error_msg)
    end subroutine warning_neci


    module subroutine neci_backtrace()
#ifdef GFORTRAN_
    call backtrace()
#elif IFORT_
    call tracebackqq()
#else
    write(stdout, *) 'Traceback not implemented for this compiler'
#endif
    end subroutine
end submodule


subroutine hidden_stop_all(sub_name, error_msg)
    ! Do not call this directly
    use constants, only: stdout, stderr
    use error_handling_neci, only: neci_flush, neci_backtrace
    use Parallel_neci, only: iProcIndex, MPIStopAll

    ! Stop calculation due to an error. Exit with code 222?
    !
    ! In: sub_name    - Calling routine
    !     error_msg   - Error message
    character(*), intent(in) :: sub_name, error_msg

    ! I found problems when the error code is larger 2^8 - 1 == 255
    integer, parameter :: error_code = 222

    call neci_flush(stdout)
    write(stdout,'(/a7)') 'ERROR.'
    write(stdout,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write(stdout,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef USE_MPI
    write(stdout,'(a12,15X,i5)') 'Processor: ', iProcIndex
#endif
    write(stdout,'(a11)') 'EXITING...'

    ! Also push this to the stderr unit, so it hopefully ends up somewhere
    ! more useful.
    call neci_flush(stderr)
    write(stderr,'(/a7)') 'ERROR.'
    write(stderr,'(a27,a)') 'NECI stops in subroutine: ',adjustl(sub_name)
    write(stderr,'(a9,18X,a)') 'Reason: ',adjustl(error_msg)
#ifdef USE_MPI
    write(stderr,'(a12,15X,i5)') 'Processor: ', iProcIndex
#endif
    write(stderr,'(a11)') 'EXITING...'

    call neci_backtrace()

#ifdef USE_MPI
    call MPIStopAll(error_code)
#else
    stop error_code
#endif

end subroutine hidden_stop_all


