module legacy_data
    use constants, only: dp, int64

! Module to store legacy data that was once upon a time in evil include files.

! Data originally in irat.inc.
! irat is the 64/(number of bits per integer).
    integer, parameter, private :: irat_test = 0
    integer, parameter :: irat = 64 / bit_size(irat_test)

! Data originally in calcp.inc.
! Probably only Alex knows what there are.
    INTEGER CALCP_N(0:1023), CALCP_NT(0:1023, 15, 2)
    real(dp) CALCP_SUM(0:1023, 3), CALCP_HSUM(0:1023)

! Data originally in csf.inc.
! The offset of the CSF open shell dets
    integer, parameter :: CSF_NBSTART = 100000

end module legacy_data

subroutine legacy_data_dummy
    implicit none
    return
end subroutine legacy_data_dummy
