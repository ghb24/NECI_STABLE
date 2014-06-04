module hist_data

    use constants
    use MemoryManager, only: TagIntType
    implicit none

    ! Should we histogram spawning
    logical :: tHistSpawn

    ! How many energy bins are there for energy histogramming
    integer :: iNoBins
    real(dp) :: BinRange

    ! Should we histogram the distribution of spin dets within a given
    ! spatial structure --> Analyse spin development
    integer(n_int), allocatable, target :: hist_spin_dist(:,:)
    real(dp), allocatable :: hist_csf_coeffs(:,:)
    integer(TagIntType) :: tag_spindist=0, tag_histcsfs=0

    real(dp), allocatable :: Histogram(:,:), AllHistogram(:,:)
    real(dp), allocatable :: BeforeNormHist(:)
    real(dp), allocatable :: InstHist(:,:), AllInstHist(:,:)
    real(dp), allocatable :: HistogramEnergy(:), AllHistogramEnergy(:)

    integer, allocatable :: HistMinInd(:), HistMinInd2(:)

    ! Store data about where spawns are going to/from
    real(dp), allocatable :: hist_excit_tofrom(:,:)
    integer(TagIntType) :: tag_hist_excit
    integer :: excit_tofrom_unit

end module
