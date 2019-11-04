#include "macros.h"
module neci_signals

    ! Here we manage the NECI signal handlers.
    !
    ! The actual legwork is done in C/C++ as the relevant constants are only
    ! accessible in signal.h.

    use iso_c_hack
    implicit none
    private

    integer, volatile :: sigint_count

    interface
        subroutine init_signals_helper() bind(c)
        end subroutine
        subroutine clear_signals() bind(c)
        end subroutine
    end interface

    public :: init_signals, clear_signals

contains

    subroutine init_signals()

        ! Set the signal handlers, and set the call count to zero

        sigint_count = 0
        call init_signals_helper()

    end subroutine

    subroutine neci_sigint(signo) bind(c, name='neci_sigint')

        use soft_exit, only: tSoftExitFound

        ! Parameter unused. Required by POSIX
        integer(c_int), intent(in), value :: signo
        character(*), parameter :: t_r = 'neci_sigint'

        unused_var(signo)

        ! Flush existing output in the stdout buffer
        ! --> Try and avoid issues if we happen to Ctrl-C during a write.
        call neci_flush(6)
        write(6,*)
        write(6,*) '----------------------------------------'
        write(6,*) 'NECI SIGINT (Ctrl-C) handler'
        write(6,*)
        sigint_count = sigint_count + 1

        if (sigint_count == 1) then
            write(6,*) 'Calculation will cleanly exit at the next update cycle.'
            write(6,*) 'To kill the application immediately, resend signal again &
                       & (Ctrl-C)'
            write(6,*) '----------------------------------------'
            write(6,*)

            ! Trigger soft exit
            tSoftExitFound = .true.
        else
            write(6,*) 'Killing calculation'
            write(6,*) '----------------------------------------'
            write(6,*)
            call neci_flush(6)
            call stop_all(t_r, "User requested")
        end if

    end subroutine

end module
