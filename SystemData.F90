module SystemData

implicit none

save

LOGICAL :: TSTARBIN,TREADINT,THFORDER,TDFREAD,TPBC,TUEG,TCPMD,THUB,tHPHF,tHPHFInts,tUHF
LOGICAL tRIIntegrals  !Read in RI 2-e integrals from RIDUMP file
! Why is so little of this commented.  'tis horrific.  AJWT.
LOGICAL :: TSPN,TCSF,TPARITY,TUSEBRILLOUIN,TEXCH,TREAL,TTILT,tUMatEps,tOneElIntMax,tOnePartOrbEnMax,tROHF,tNoBrillouin,tCSFOLD
LOGICAL :: tStoreSpinOrbs   !This is set when the orbitals are stored in spin-orbital notation
LOGICAL :: tVirtCoulombMax,tVirtExchangeMin,tHijSqrdMin,tDiagonalizehij,tHFSingDoubExcMax,tSpinOrbs,tReadInCoeff,tUseMP2VarDenMat
LOGICAL :: TALPHA,TSTOREASEXCITATIONS,TBIN,tStarStore,tVASP,tOffDiagSqrdMin,tOffDiagSqrdMax,tOffDiagMax,tShakeDelay
LOGICAL :: tSeparateOccVirt,tMerTwist,tExactSizeSpace,tRotatedOrbs,tImportanceSample,tERLocalization,tOffDiagMin,tFindCINatOrbs
LOGICAL :: TNoRenormRandExcits,tAssumeSizeExcitgen,tCycleOrbs,tROIteration,tShakeIter,tRotateOccOnly,tDoubExcMin,tUseHFOrbs
LOGICAL :: tNonUniRandExcits,tNoSymGenRandExcits,tRotateOrbs,tLagrange,tShake,tShakeApprox,tRotateVirtOnly,tMaxHLGap,tCacheFCIDUMPInts
LOGICAL :: tNoFailAb, tLatticeGens ! These are temporary inputs for UEG
INTEGER :: LMS,STOT,IPARITY(5),NMAXX,NMAXY,NMAXZ,NMSH,COULDAMPORB,ElecPairs,ROIterMax,iRanLuxLev,DiagMaxMinFac,OneElMaxMinFac
INTEGER :: iPeriodicDampingType,ISTATE,NEL,ITILTX,ITILTY,nOccAlpha,nOccBeta,ShakeIterMax,ShakeStart,MaxMinFac,MaxABPairs
REAL*8 :: BOX,BOA,COA,FUEGRS,fRc,FCOUL,OrbECutoff,UHUB,BHUB,DiagWeight,OffDiagWeight,OrbEnMaxAlpha
REAL*8 :: ALPHA,FCOULDAMPBETA,FCOULDAMPMU,TimeStep,ConvergedForce,ShakeConverged,UMatEps,OneElWeight

LOGICAL :: tListDets    !Means that a list of allowed determinants in FciMC will be read in an particles are only allowed here.

LOGICAL :: tMCSizeSpace
INTEGER*8 :: CalcDetPrint,CalcDetCycles   !parameters for the MC determination of the FCI determinant space size.

! Used to be stored in Integrals
INTEGER :: ORBORDER(8,2)

LOGICAL :: tFixLz   !This indicates that in FCIMC, the Lz of the determinants is fixed at LzTot
INTEGER :: LzTot,iMaxLz     !LzTot is the total Ml quantum number of the state to converge upon. iMaxLz is the abs(maximum Ml basis function).

! NIfTot indicates upper bound of determinants in bit form
! NIfD is final byte used to represent determinants (NIfD+1 bytes total)
! NIfY is the number of bytes used to represent a Yamanouchi symbol
! NIfP is an integer that is added on to the end of NIfD in CASSTAR calculations - keeps track of the parent determinant of spawned walkers.
! NB. bit representations are zero indexed.
INTEGER :: NIfD, NIfY,NIfP,NIfDBO
integer :: NIfTot

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
   INTEGER :: Ml            !This is the Ml symmetry of the orbital
   INTEGER :: spacer    ! The spacer is there to make sure we have a structure which is a multiple of 8-bytes for 64-bit machines.
   TYPE(Symmetry) :: sym
END TYPE

! Empty basis function is used in many places.
! This is useful so if BasisFn changes, we don't have to go
! through the code and change the explicit null statements.
type(BasisFn) :: NullBasisFn=BasisFn((/0,0,0/),0,0,0,Symmetry(0))

integer, PARAMETER :: BasisFNSize=SymmetrySize+6
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

!Set if we don't reorder the HF orbitals on entry
LOGICAL :: tHFNoOrder
! When set, ignore differences in orbital energies between pairs of orbitals (which should be beta/alpha)
!  and group them under the same symrep
LOGICAL :: tSymIgnoreEnergies

end module SystemData
