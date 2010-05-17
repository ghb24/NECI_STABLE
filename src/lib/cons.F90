module constants

#ifdef PARALLEL
uSE mpi
#endif

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
integer, parameter :: sizeof_complexdp = 16
integer, parameter :: sizeof_sp = 4
#ifdef __CMPLX
integer, parameter :: sizeof_helement = 16
integer, parameter :: lenof_sign = 2
#else
integer, parameter :: sizeof_helement = 8
integer, parameter :: lenof_sign = 1
#endif

#ifdef __INT64

! Kind parameter for 64-bit integers.
integer, parameter :: n_int=selected_int_kind(18)

! MPI integer kind associated with n_int.
#ifdef PARALLEL
integer, parameter :: MpiDetInt=MPI_INTEGER8
#endif

#else

! Kind parameter for 32-bit integers.
integer, parameter :: n_int=selected_int_kind(8)

! MPI integer kind associated with n_int.
#ifdef PARALLEL
integer, parameter :: MpiDetInt=MPI_INTEGER
#endif

#endif

! Number of bits in an n_int integer.
! Note that PGI (at least in 10.3) has a bug which causes
! bit_size(int(0,n_int)) to return an incorrect value.
integer, parameter :: bits_n_int = bit_size(0_n_int)
! Number of bytes in an n_int integer.
integer, parameter :: size_n_int = bits_n_int/8
! Index of last bit in an n_int integer (bits are indexed 0,1,...,bits_n_int-1).
integer, parameter :: end_n_int = bits_n_int - 1

#ifndef PARALLEL
! This should not be used in serial.  Set to a nonsense value.
integer, parameter :: MpiDetInt=-1
#endif

end module constants
