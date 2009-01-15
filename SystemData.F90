module SystemData

implicit none

save

LOGICAL :: TSTARBIN,TREADINT,THFORDER,TDFREAD,TPBC,TUEG,TCPMD,THUB
LOGICAL :: TSPN,TCSF,TPARITY,TUSEBRILLOUIN,TEXCH,TREAL,TTILT
LOGICAL :: TALPHA,TSTOREASEXCITATIONS,TBIN,tStarStore,tVASP
LOGICAL :: TNoRenormRandExcits,tAssumeSizeExcitgen,tCycleOrbs
LOGICAL :: tNonUniRandExcits,tNoSymGenRandExcits
INTEGER :: LMS,STOT,IPARITY(5),NMAXX,NMAXY,NMAXZ,NMSH,COULDAMPORB,ElecPairs
INTEGER :: iPeriodicDampingType,ISTATE,NEL,ITILTX,ITILTY,nOccAlpha,nOccBeta
REAL*8 :: BOX,BOA,COA,FUEGRS,fRc,FCOUL,OrbECutoff,UHUB,BHUB
REAL*8 :: ALPHA,FCOULDAMPBETA,FCOULDAMPMU

! Used to be stored in Integrals
INTEGER :: ORBORDER(8,2)

! From NECICB
integer :: lmsBasis

TYPE Symmetry
   SEQUENCE
   INTEGER*8 :: S
END TYPE

integer, PARAMETER :: SymmetrySize=2
integer, PARAMETER :: SymmetrySizeB=SymmetrySize*8
TYPE BasisFN
   SEQUENCE
   INTEGER :: k(3)
   INTEGER :: Ms
   TYPE(Symmetry) :: sym
END TYPE

integer, PARAMETER :: BasisFNSize=SymmetrySize+4
integer, PARAMETER :: BasisFNSizeB=BasisFNSize*8


TYPE(BASISFN) :: SymRestrict
INTEGER :: nBasisMax(5,7)
REAL*8 :: ALAT(5)
REAL*8 :: ECore
INTEGER :: nBasis
integer :: nMax
integer :: nnr
integer :: nocc
REAL*8 :: OMEGA
logical :: tSpinPolar
INTEGER :: iSpinSkip ! Often referred to as ISS.

!From Calc  
REAL*8 :: Beta
        
!Renewed for compile

! List of orbital energies.
! Arr(:,1)
!   ordered by energy.  Moreover when we have a set of degenerate orbitals with
!   different symmetry, they are ordered within that set such that orbitals with
!   the same symmetry are next to each other, and then by Ms=-1,1.
!   Arr(10,1) is the energy of the 10th lowest energy
!   spin-orbital.
! Arr(:,2)
!     ordered by spin-orbital index.  Arr(10,1) is the energy of the 10th
!     spin-orbital (given the index scheme in use).
! Reallocated with the correct (new) size during freezing.
REAL*8, pointer :: Arr(:,:) 
INTEGER :: tagArr

! Lists orbitals in energy order. i.e. Brr(1) is the lowest energy orbital
INTEGER, pointer :: BRR(:) 
INTEGER :: tagBrr

Type(BasisFN), pointer :: G1(:)  ! Info about the basis functions.
INTEGER :: tagG1

INTEGER :: LMS2

!  Set if we turn symmetry off
LOGICAL :: lNoSymmetry

end module SystemData
