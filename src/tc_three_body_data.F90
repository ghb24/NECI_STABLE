module tc_three_body_data
  use constants
  use FciMCData, only: ll_node
  implicit none

  ! biases for 3-body excitation generation
  real(dp) :: pTriples = 0.0_dp
  real(dp) :: p0B, p0A, p2B, p1B
  real(dp) :: lMatEps = 0.0_dp
  logical :: tReadPTriples = .false.
  ! precomputed probabilities
  real(dp) :: pgen0B, pgen1B, pgen2B, pgen3B
  ! number of empty alpha/beta electrons
  integer :: nUnoccAlpha, nUnoccBeta

  ! option to reduce the k-matrix element for same-spin excitations
  logical :: tDampKMat = .false.
  logical :: tDampLMat = .false.
  
  ! if kMatrix is used
  logical :: tUseKMat
  
  ! if lMat does not have permutational symmetry (required for spin-projection)
  logical :: tSymBrokenLMat = .false.
  logical :: tSpinCorrelator = .false.

  ! using a sparse format to store the 6-index integrals
  logical :: tSparseLMat

  type lMat_t
     HElement_t(dp), pointer :: LMatPtr(:)
     integer(int64), pointer :: indexPtr(:)
     integer :: tag, indexTag
     integer(MPIArg) :: shm_win, index_win
     integer :: nInts
     type(ll_node), pointer :: hTable(:)
     integer :: htSize
  end type lMat_t


end module tc_three_body_data
