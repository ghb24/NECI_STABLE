module GreensFuncData
    use constants, only: dp, int64, n_int
    use global_utilities
    implicit none
    save

    real(dp) , allocatable :: G_ret(:,:,:), G_adv(:,:,:)
    real(dp) , allocatable :: G_ret_all(:,:,:), G_adv_all(:,:,:)
    real(dp) , allocatable :: G_ret_sq(:,:,:), G_adv_sq(:,:,:)
    integer(n_int) , allocatable :: CurrentDetsSaved(:,:)
    integer , allocatable :: HashIndexSaved(:,:)
    integer :: TotWalkersSaved

    integer :: nTimePnts_ret,nTimePnts_adv  !Input parameters (default to 10?)
    integer :: GFSamples    !Input parameter (Independant samples of GF)
    integer :: nIterEquilGF !Input parameter (equilibration iterations)
    logical :: tJustAdvGF   !Input parameter (def: F)
    logical :: tJustRetGF   !Input parameter (def: F)
    real(dp), dimension(2) :: GreensImTime  !Input parameter
    integer, dimension(2) :: GreensIters    !Conversion of the above

    logical :: tGreensFuncs !Whether to calculate greensfunctions
    logical :: tNoDiag      !Whether to calculate diagonal elements of the GF
    logical :: tSpecificMatEls  !Manually specify the specific GF matrix elements to calculate
    integer :: iSpecificMatEls  !The number of matrix elements to calculate
    integer, allocatable :: SpecificMatEls(:)   !The specific matrix element pairs (triangular indexed)
    integer, allocatable :: SpecificMatEls_temp(:,:)  !The specific matrix element pairs (expanded for reading in)
    integer, allocatable :: Allowed_j(:)    !The j's that have allowed i's to give matrix elements
    integer :: Contributing_Js              !The number of J's that we need to run over

    integer :: nel_orig !the original number of electrons
    integer :: nOccAlpha_orig   !...alpha electrons
    integer :: nOccBeta_orig    !... and beta electrons

    type(timer) :: Greens_Time

    !!  MEM OPTIONS  !!
    logical :: tJustSpectrum    !Skip calculation of GF, and just go straight for the spectrum?
    logical :: tCalcSpectrum    !Whether to calculate the spectrum or not at the end of the calc.

end module GreensFuncData
