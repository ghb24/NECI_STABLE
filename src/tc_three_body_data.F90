module tc_three_body_data
  use constants
  use FciMCData, only: ll_node
  use lMat_indexing, only: lMatIndSym
  use shared_array
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

  ! Is the permutational symmetry of the 6-index integrals broken?
  logical :: tSymBrokenLMat = .false.
  
  ! using a sparse format to store the 6-index integrals
  logical :: tSparseLMat
  ! storage options for hdf5
  logical :: tHDF5LMat = .false.

  ! Number of entries in a UMat obj.
  integer(int64) :: twoIndexSize, fourIndexSize  

  logical :: tLMatCalc !Calculate LMat on the fly. tcfactors.h5 is required.
  real :: lMatCalcHFactor !Size of hash table used to cache LMat values as a fraction of total LMat size.
  integer(int64) :: lMatCalcHSize, lMatCalcHUsed !Size of hash table and number of entries used
  integer :: lMatCalcTot, lMatCalcHit !Total calls of lMatCalc and the number of times the value found in cache.
  integer :: lMatCalcStatsIters = 100 !How often to print lMatCalc statistics.

end module tc_three_body_data
