module constants

implicit none

! Constant data.

integer, parameter :: sp = selected_real_kind(6,37)
integer, parameter :: dp = selected_real_kind(15,307)
integer, parameter :: qp = selected_real_kind(33,4931)
integer, parameter :: int32 = selected_int_kind(8)
integer, parameter :: int64 = selected_int_kind(15)
integer, parameter :: cLineSize = 8

real(dp), parameter ::  PI    = 3.1415926535897932384626433832795028841971693993751_dp
real(dp), parameter ::  PI2   = 9.8696044010893586188344909998761511353136994072408_dp
real(dp), parameter ::  THIRD = 0.3333333333333333333333333333333333333333333333333_dp
real(dp), parameter ::  Root2 = 1.4142135623730950488016887242096980785696718753769_dp
real(dp), parameter ::  OverR2 = 1.0_dp/Root2
real(dp), parameter :: EPS = 0.0000000000001_dp
real(dp), parameter :: INFINITY = huge(1.0_dp)

integer :: temp
integer, parameter :: sizeof_int = kind(temp)   !Default integer size (not necessarily = no. bytes)
integer, parameter :: bits_int = bit_size(temp)
integer, parameter :: bytes_int = bits_int/8
!potential hack for molpro, which seems to support a compiler which doesn't like the kind() intrinsic..?
!integer, parameter :: sizeof_int = selected_int_kind(digits(huge(temp)))   !Default integer size (not necessarily = no. bytes)

logical :: temp2=.true.
integer, parameter :: sizeof_log = kind(temp2)

integer, parameter :: sizeof_int32 = 4
integer, parameter :: sizeof_int64 = 8
integer, parameter :: sizeof_dp = 8
integer, parameter :: sizeof_complexdp = 16
integer, parameter :: sizeof_sp = 4

! number of possible excitations per step
integer, parameter :: maxExcit = 3

#if defined(CMPLX_)
! The ratio of lenof_sign / inum_runs
integer, parameter :: rep_size = 2
#else
integer, parameter :: rep_size = 1
#endif
! Give ourselves the option of lenof_sign/inum_runs being a runtime
! variable, rather than a compile-time constant
#if defined(PROG_NUMRUNS_)
#if defined(CMPLX_)
!Complex integrals, (arbitrary, run-time) multiple replicas
    integer :: nreplicas = 1    !1 or 2   (for replica sampling, not multiple states)
    integer :: lenof_sign       !2 x inum_runs (2 for complex x number of seperate wavefuncs sampled)
    integer :: inum_runs        !nreplicas x nstates
    integer :: lenof_sign_kp
    integer, parameter :: lenof_sign_max = 16
    integer, parameter :: inum_runs_max = 16
    integer, parameter :: sizeof_helement = 16
    HElement_t(dp), parameter :: HEl_zero = cmplx(0.0_dp, 0.0_dp, dp)
#else
!Real integrals, (arbitrary, run-time) multiple replicas
!Also, define MULTI_RUN below, to mean real (not complex) double or multiple run code
#define MULTI_RUN
    integer :: nreplicas = 1
    integer :: lenof_sign
    integer :: inum_runs
    integer :: lenof_sign_kp
    integer, parameter :: lenof_sign_max = 16
    integer, parameter :: inum_runs_max = 16
    integer, parameter :: sizeof_helement = 16
    real(dp), parameter :: HEl_zero = 0.0_dp
#endif

#elif defined(DOUBLERUN_)
#if defined(CMPLX_)
!Complex integrals, double replica
    integer, parameter :: nreplicas = 2
    integer, parameter :: lenof_sign = 4
    integer, parameter :: inum_runs = 2
    integer, parameter :: lenof_sign_kp = 4
    integer, parameter :: lenof_sign_max = lenof_sign
    integer, parameter :: inum_runs_max = inum_runs
    integer, parameter :: sizeof_helement = 16
    HElement_t(dp), parameter :: HEl_zero = cmplx(0.0_dp, 0.0_dp, dp)
#else
!Real integrals, double replica
#define MULTI_RUN
    integer, parameter :: nreplicas = 2
    integer, parameter :: lenof_sign = 2
    integer, parameter :: inum_runs = lenof_sign
    integer, parameter :: lenof_sign_kp = 2
    integer, parameter :: lenof_sign_max = lenof_sign
    integer, parameter :: inum_runs_max = inum_runs
    integer, parameter :: sizeof_helement = 8
    real(dp), parameter :: HEl_zero = 0.0_dp
#endif

#else
#if defined(CMPLX_)
!Complex integrals, single replica
    integer, parameter :: nreplicas = 1
    integer, parameter :: lenof_sign = 2
    integer, parameter :: inum_runs = 1
    integer, parameter :: lenof_sign_kp = 2
    integer, parameter :: lenof_sign_max = lenof_sign
    integer, parameter :: inum_runs_max = inum_runs
    integer, parameter :: sizeof_helement = 16
    complex(dp), parameter :: HEl_zero = cmplx(0.0_dp, 0.0_dp, dp)
#else
!Real, single replica
    integer, parameter :: nreplicas = 1
    integer, parameter :: lenof_sign = 1
    integer, parameter :: inum_runs = 1
    integer, parameter :: lenof_sign_kp = 1
    integer, parameter :: lenof_sign_max = lenof_sign
    integer, parameter :: inum_runs_max = inum_runs
    integer, parameter :: sizeof_helement = 8
    real(dp), parameter :: HEl_zero = 0.0_dp
#endif
#endif

real(dp), dimension(lenof_sign_max), parameter :: null_part = 0.0_dp

!This is the integer type which is used in MPI call arguments
!This should normally be integer(4)'s.
integer, parameter :: MPIArg=int32

#ifdef INT64_

! Kind parameter for 64-bit integers.
integer, parameter :: n_int=int64
logical :: build_64bit = .true.

integer, parameter :: int_rdm=int64

#else

! Kind parameter for 32-bit integers.
integer, parameter :: n_int=int32
logical :: build_64bit = .false.

integer, parameter :: int_rdm=int32

#endif

! Number of bits in an n_int integer.
! Note that PGI (at least in 10.3) has a bug which causes
! bit_size(int(0,n_int)) to return an incorrect value.
integer(n_int) :: temp3=0
integer, parameter :: bits_n_int = bit_size(temp3)
! to avoid too long lines in GUGA macros.h introduce a shorter named bits_n_int
integer, parameter :: bni_ = bits_n_int
! also define a bits_n_int/2 to further reduce the character count
integer, parameter :: bn2_ = bits_n_int/2
! Number of bytes in an n_int integer.
integer, parameter :: size_n_int = bits_n_int/8
! Index of last bit in an n_int integer (bits are indexed 0,1,...,bits_n_int-1).
integer, parameter :: end_n_int = bits_n_int - 1

! Number of bits in an int_rdm integer.
integer(int_rdm) :: temp4=0
integer, parameter :: bits_int_rdm = bit_size(temp4)
! Number of bytes in an int_rdm integer.
integer, parameter :: size_int_rdm = bits_int_rdm/8
! Index of last bit in an int_rdm integer (bits are indexed 0,1,...,bits_n_int-1).
integer, parameter :: end_int_rdm = bits_int_rdm - 1

    integer, parameter :: iout = 6

    ! Internal state storage for the stats_out integration
    ! n.b. This shouldn't be here, but there is nowhere els eto put it
    type write_state_t
        integer :: funit, cols, cols_mc
        logical :: init, mc_out, prepend
    end type

! Typedef for HDF5 variables
integer, parameter :: hdf_err = int32
integer, parameter :: hdf_log = int32

integer, parameter :: NEL_UNINITIALIZED = -1

end module constants
