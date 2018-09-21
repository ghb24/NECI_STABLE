module tc_three_body_data
  use constants
  implicit none

  ! biases for 3-body excitation generation
  real(dp) :: pTriples, p0B, p1B, p2B
  ! number of empty alpha/beta electrons
  integer :: nUnoccA, nUnoccB
  
  HElement_t(dp), allocatable :: LMat(:)

end module tc_three_body_data
