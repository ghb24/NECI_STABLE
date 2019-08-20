module IntegralsData

use iso_c_hack
use constants, only: dp, MPIArg
use MemoryManager, only: TagIntType

IMPLICIT NONE
save

LOGICAL :: TQUADRHO,TEXPRHO,THFBASIS,THFCALC,TCALCREALPROD
LOGICAL :: TRHF,TReadTUMat,TReadHF,TQuadValMax,TQuadVecMax
LOGICAL :: TSUMPROD,TCALCRHOPROD,TDISCONODES,TCalcExcitStar
LOGICAL :: TJustQuads,TNoDoubs,TDiagStarStars,TExcitStarsRootChange
LOGICAL :: TRmRootExcitStarsRootChange,TLinRootChange,tPartFreezeCore,tPartFreezeVirt

INTEGER :: NTAY(2),nHFit,NFROZEN,NTFROZEN,NFROZENIN,NTFROZENIN
INTEGER :: NRSTEPSMAX,IHFMETHOD
INTEGER :: NPartFrozen,NHolesFrozen,NVirtPartFrozen,NElVirtFrozen


real(dp) :: NRCONV,RFCONV,OrbOrder2(8)
real(dp) :: HFMix,HFEDelta,HFCDelta
real(dp) :: HFRand
real(dp) :: DMatEpsilon !  The cutoff for density matrix elements
Logical :: tPostFreezeHF ! Do we do HF after freezing
logical :: tDumpFCIDUMP !Do we write out an FCIDUMP file (after freezing)

!  From NECI.F
! UMAT stores 4-index, 2-electron integrals.  Lookup is via the UMATIND function
! (in the UMatCache module).
! Warning! It is best (ie only safe) to access 4-index integrals via the
! get_umat_el function.  This enables the appropriate transformations, complex
! conjugations and storage implementations (eg with UMAT2D/UMatCache/etc) to be
! used as appropriate.
HElement_t(dp), pointer :: UMAT(:)
INTEGER(MPIArg) :: umat_win
INTEGER(TagIntType) :: tagUMat=0
COMPLEX(dp) ,pointer :: FCK(:) => null()
INTEGER(TagIntType) :: tagFCK=0
INTEGER :: NMAX
real(dp) :: CST

! from Calc      
real(dp) :: ChemPot
! Set if we need virtual orbitals  (usually set).  Will be unset (by Calc readinput) if I_VMAX=1 and TENERGY is false
Logical :: tNeedsVirts  

! Details to permit reversing freezing for convenience.
integer :: nel_pre_freezing
integer, allocatable :: frozen_orb_list(:), frozen_orb_reverse_map(:)
integer(TagIntType) :: tagFrozen=0, tagFrozenMap=0

end module IntegralsData
