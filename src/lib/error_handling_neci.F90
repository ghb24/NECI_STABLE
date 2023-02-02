#include "macros.h"
module error_handling_neci
    better_implicit_none
    private
    public :: stop_all, neci_backtrace, neci_flush, warning_neci

    interface
        pure subroutine hidden_stop_all(sub_name, error_msg)
            character(*), intent(in) :: sub_name, error_msg
        end subroutine

        module subroutine neci_backtrace()
        end subroutine

        module subroutine neci_flush(un)
            integer, intent(in) :: un
        end subroutine neci_flush

        module subroutine warning_neci(sub_name, error_msg)
            != Print a warning message in a (helpfully consistent) format.
            !=
            != In:
            !=    sub_name:  calling subroutine name.
            !=    error_msg: error message.
            character(*), intent(in) :: sub_name, error_msg
        end subroutine warning_neci
    end interface

contains

    pure subroutine stop_all(sub_name, error_msg)
        character(*), intent(in) :: sub_name, error_msg
        call hidden_stop_all(sub_name, error_msg)
    end subroutine stop_all
end module


