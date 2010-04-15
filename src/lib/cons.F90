module constants

! Constant data.

integer, parameter :: sp = selected_real_kind(6,37)
integer, parameter :: dp = selected_real_kind(15,307)
integer, parameter :: qp = selected_real_kind(33,4931)
integer, parameter :: int32 = selected_int_kind(8)
integer, parameter :: int64 = selected_int_kind(18)

real(dp), parameter ::  PI    = 3.1415926535897932384626433832795028841971693993751_dp
real(dp), parameter ::  PI2   = 9.8696044010893586188344909998761511353136994072408_dp
real(dp), parameter ::  THIRD = 0.3333333333333333333333333333333333333333333333333_dp

integer, parameter :: sizeof_int = 4
integer, parameter :: sizeof_dp = 8
integer, parameter :: sizeof_sp = 4
#ifdef __CMPLX
integer, parameter :: sizeof_helement = 16
#else
integer, parameter :: sizeof_helement = 8
#endif

end module constants
