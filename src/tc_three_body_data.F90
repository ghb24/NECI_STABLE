module tc_three_body_data
  use constants
  use FciMCData, only: ll_node
  use procedure_pointers, only: lMatInd_t
  use lMat_indexing, only: lMatIndSym  
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
  ! storage options for hdf5
  logical :: tHDF5LMat = .false.

  ! Number of entries in a UMat obj.
  integer(int64) :: twoIndexSize, fourIndexSize

  type lMat_t
     HElement_t(dp), pointer :: LMatPtr(:)
     integer(int64), pointer :: indexPtr(:)
     integer :: tag, indexTag
     integer(MPIArg) :: shm_win, index_win
     integer(int64) :: nInts
     type(ll_node), pointer :: hTable(:)
     integer(int64) :: htSize
     procedure(lMatInd_t), nopass, pointer :: indexFunc => lMatIndSym
  end type lMat_t  
  
  

  logical :: tLMatCalc !Calculate LMat on the fly. tcfactors.h5 is required.
  real :: lMatCalcHFactor !Size of hash table used to cache LMat values as a fraction of total LMat size.
  integer(int64) :: lMatCalcHSize, lMatCalcHUsed !Size of hash table and number of entries used
  integer :: lMatCalcTot, lMatCalcHit !Total calls of lMatCalc and the number of times the value found in cache.
  integer :: lMatCalcStatsIters = 100 !How often to print lMatCalc statistics.

  real :: lMatABCalcHFactor !Size of hash table used to cache LMatAB values as a fraction of total LMatAB size.
  integer(int64) :: lMatABCalcHSize, lMatABCalcHUsed !Size of hash table and number of entries used
  integer :: lMatABCalcTot, lMatABCalcHit !Total calls of lMatABCalc and the number of times the value found in cache.
end module tc_three_body_data
