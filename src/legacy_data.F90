module legacy_data
    use constants, only: int64

! Module to store legacy data that was once upon a time in evil include files.

! Data originally in irat.inc.
! irat is the 64/(number of bits per integer).
    integer, parameter, private :: irat_test = 0
    integer, parameter :: irat = 64 / bit_size(irat_test)

end module legacy_data
