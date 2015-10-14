module fmt_utils

    use constants
    implicit none
    private

    interface int_fmt
        module procedure int64_fmt
        module procedure int32_fmt
    end interface

    public :: int_fmt

contains

    function int32_fmt(i, padding) result(fmt_str)

        ! Return the format string that exactly contains the specified
        ! integer

        integer(int32), intent(in) :: i
        integer, intent(in), optional :: padding
        character(4) :: fmt_str
        real(dp) :: logi

        if (i == 0 .or. i == 1) then
            logi = 1.0
        else
            logi = log10(real(abs(i)+1.0, dp))
        end if
        if (i < 0) logi = logi + 1

        fmt_str = int_fmt_local(logi, padding)

    end function

    function int64_fmt(i, padding) result(fmt_str)

        ! Return the format string that exactly contains the specified
        ! integer

        integer(int64), intent(in) :: i
        integer, intent(in), optional :: padding
        character(4) :: fmt_str
        real(dp) :: logi

        if (i == 0 .or. i == 1) then
            logi = 1.0
        else
            logi = log10(real(abs(i)+1.0, dp))
        end if
        if (i < 0) logi = logi + 1

        fmt_str = int_fmt_local(logi, padding)

    end function

    function int_fmt_local(logi, padding) result(fmt_str)

        ! Return the format string that exactly contains the specified
        ! integer (local worker function)

        real(dp), intent(in) :: logi
        integer, intent(in), optional :: padding
        character(4) :: fmt_str
        integer :: ndigit

        ndigit = int(ceiling(logi))

        if (present(padding)) ndigit = ndigit + padding

        if (ndigit < 10) then
            write(fmt_str, '("i",i1)') ndigit
        else if (ndigit < 100) then
            write(fmt_str, '("i",i2)') ndigit
        else
            ! n.b. should never hit this. Will hit integer overflow first
            write(fmt_str, '("i",i3)') ndigit
        end if

    end function


end module
