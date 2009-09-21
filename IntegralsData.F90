module IntegralsData

use HElem, only: HElement

IMPLICIT NONE
save

LOGICAL :: TQUADRHO,TEXPRHO,THFBASIS,THFCALC,TCALCREALPROD
LOGICAL :: TRHF,TReadTUMat,TReadHF,TQuadValMax,TQuadVecMax
LOGICAL :: TSUMPROD,TCALCRHOPROD,TDISCONODES,TCalcExcitStar
LOGICAL :: TJustQuads,TNoDoubs,TDiagStarStars,TExcitStarsRootChange
LOGICAL :: TRmRootExcitStarsRootChange,TLinRootChange,tPartFreezeCore

INTEGER :: NTAY(2),nHFit,NFROZEN,NTFROZEN,NFROZENIN,NTFROZENIN
INTEGER :: NRSTEPSMAX,IHFMETHOD
INTEGER :: NPartFrozen,NHolesFrozen


REAL*8 :: NRCONV,RFCONV,OrbOrder2(8)
REAL*8 :: HFMix,HFEDelta,HFCDelta
REAL*8 :: HFRand
REAL*8 :: DMatEpsilon !  The cutoff for density matrix elements
Logical :: tPostFreezeHF ! Do we do HF after freezing

!  From NECI.F
TYPE(HElement), pointer :: UMAT(:)      
INTEGER :: tagUMat=0
COMPLEX*16,pointer :: FCK(:)
INTEGER :: tagFCK=0
INTEGER :: NMAX
REAL*8 :: CST

! from Calc      
real*8 :: ChemPot
Logical :: tNeedsVirts  ! Set if we need virtual orbitals  (usually set).  Will be unset (by Calc readinput) if I_VMAX=1 and TENERGY is false

end module IntegralsData
