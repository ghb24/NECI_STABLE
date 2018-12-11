module tc_three_body_data
  use constants
  implicit none

  ! biases for 3-body excitation generation
  real(dp) :: pTriples, p0B, p0A, p2B, p1B
  ! precomputed probabilities
  real(dp) :: pgen0B, pgen1B, pgen2B, pgen3B
  ! number of empty alpha/beta electrons
  integer :: nUnoccAlpha, nUnoccBeta
  
  HElement_t(dp), allocatable :: LMat(:)

end module tc_three_body_data
