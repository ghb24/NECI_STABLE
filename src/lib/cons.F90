module constants

#ifdef PARALLEL
uSE mpi
#endif

! Constant data.

integer, parameter :: dp = selected_real_kind(15,307)
integer, parameter :: int64 = selected_int_kind(18)
integer, parameter :: int32 = selected_int_kind(8)

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

#ifdef __INT64

integer, parameter :: n_int=selected_int_kind(18)
integer, parameter :: size_n_int = 8
integer, parameter :: bits_n_int = 64
#ifdef PARALLEL
integer, parameter :: MpiDetInt=MPI_INTEGER8
#endif

#else

integer, parameter :: n_int=selected_int_kind(8)
integer, parameter :: size_n_int = 4
integer, parameter :: bits_n_int = 32
#ifdef PARALLEL
integer, parameter :: MpiDetInt=MPI_INTEGER
#endif

#endif

end module constants
