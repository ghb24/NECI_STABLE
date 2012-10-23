module constants

!All use of mpi routines come from this module
#ifdef PARALLEL
#ifndef CBINDMPI 
uSE mpi
#endif
#endif
implicit none

! Constant data.

integer, parameter :: sp = selected_real_kind(6,37)
integer, parameter :: dp = selected_real_kind(15,307)
integer, parameter :: qp = selected_real_kind(33,4931)
integer, parameter :: int32 = selected_int_kind(8)
integer, parameter :: int64 = selected_int_kind(15)

real(dp), parameter ::  PI    = 3.1415926535897932384626433832795028841971693993751_dp
real(dp), parameter ::  PI2   = 9.8696044010893586188344909998761511353136994072408_dp
real(dp), parameter ::  THIRD = 0.3333333333333333333333333333333333333333333333333_dp
real(dp), parameter ::  Root2 = 1.4142135623730950488016887242096980785696718753769_dp
!real(dp), parameter ::  Root2 = sqrt(2.0_dp)   !Removed since sun comiler didn't like this: bug 3853 

integer :: temp
integer, parameter :: sizeof_int = kind(temp)   !Default integer size (not necessarily = no. bytes)
!potential hack for molpro, which seems to support a compiler which doesn't like the kind() intrinsic..?
!integer, parameter :: sizeof_int = selected_int_kind(digits(huge(temp)))   !Default integer size (not necessarily = no. bytes)

logical :: temp2=.true.
integer, parameter :: sizeof_log = kind(temp2) 

integer, parameter :: sizeof_int32 = 4
integer, parameter :: sizeof_int64 = 8
integer, parameter :: sizeof_dp = 8
integer, parameter :: sizeof_complexdp = 16
integer, parameter :: sizeof_sp = 4
#ifdef __CMPLX
integer, parameter :: sizeof_helement = 16
integer, parameter :: lenof_sign = 2
integer, dimension(2), parameter :: null_part = 0
#else
integer, parameter :: sizeof_helement = 8
integer, parameter :: lenof_sign = 1
integer, dimension(1), parameter :: null_part = 0
#endif

!This is the integer type which is used in MPI call arguments
!This should normally be integer(4)'s.
integer, parameter :: MPIArg=int32

#ifdef __INT64

! Kind parameter for 64-bit integers.
integer, parameter :: n_int=int64

#else

! Kind parameter for 32-bit integers.
integer, parameter :: n_int=int32

#endif

! Number of bits in an n_int integer.
! Note that PGI (at least in 10.3) has a bug which causes
! bit_size(int(0,n_int)) to return an incorrect value.
integer(n_int) :: temp3=0
integer, parameter :: bits_n_int = bit_size(temp3)
! Number of bytes in an n_int integer.
integer, parameter :: size_n_int = bits_n_int/8
! Index of last bit in an n_int integer (bits are indexed 0,1,...,bits_n_int-1).
integer, parameter :: end_n_int = bits_n_int - 1


end module constants
