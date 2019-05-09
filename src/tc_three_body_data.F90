module tc_three_body_data
  use constants
  implicit none

  ! biases for 3-body excitation generation
  real(dp) :: pTriples = 0.0_dp
  real(dp) :: p0B, p0A, p2B, p1B
  logical :: tReadPTriples = .false.
  ! precomputed probabilities
  real(dp) :: pgen0B, pgen1B, pgen2B, pgen3B
  ! number of empty alpha/beta electrons
  integer :: nUnoccAlpha, nUnoccBeta

  ! option to reduce the k-matrix element for same-spin excitations
  logical :: tDampKMat
  ! the two-index tc integrals K
  real(dp), pointer :: kMat(:), kMatLin(:), kMatQuad(:)
  integer(MPIArg) :: kMat_win, kMatLin_win, kMatQuad_win
  
  ! if kMatrix is used
  logical :: tUseKMat

end module tc_three_body_data
