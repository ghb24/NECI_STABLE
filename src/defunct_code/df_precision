module precision
!define single and double precision REALs
!USAGE: 
!   USE precision
!   real(dp) :: X
!   X = 1.0e-12_dp
integer, parameter :: sp = selected_real_kind(6,37)
integer, parameter :: dp = selected_real_kind(15,307)
!Same for integers:
integer, parameter :: si = selected_int_kind(9)  !ints between +-10^9
integer, parameter :: li = selected_int_kind(18) !ints between +-10^18
!
!Number of bytes in reals and integers
integer, parameter :: int_bytes = si
integer, parameter :: real_bytes = dp
!integer(li), parameter :: int_bytes_li = int(si,li)
!integer(li), parameter :: real_bytes_li = int(dp,li)
end module precision


