module SystemData

    use constants, only: n_int

implicit none

save

! Why is so little of this commented.  'tis horrific.  AJWT.
! Agree.
! But...
! I blame the original authors for setting the precedent... JSS
! ;-)

logical :: tStarBin, tReadInt, tHFOrder, tDFRead, tPBC, tUEG, tCPMD, tHUB
logical :: tHPHF, tHPHFInts, tUHF, tSPN, tParity, tUseBrillouin, tExch, tReal
logical :: tTilt, tUmatEps, tOneElIntMax, tOnePartOrbEnMax, tROHF, tBrillouinsDefault
logical :: tNoBrillouin, tVirtCoulombMax, tVirtExchangeMin, tHijSqrdMin
logical :: tDiagonalizehij, tHFSingDoubExcMax, tSpinOrbs, tReadInCoeff
logical :: tUseMP2VarDenMat, tAlpha, tStoreAsExcitations, tBin, tStarStore
logical :: tVASP, tOffDiagSqrdMin, tOffDiagSqrdMax, tOffDiagmax, tShakeDelay
logical :: tSeparateOccVirt, tExactSizeSpace, tRotatedOrbs
logical :: tImportanceSample, tERLocalization, tOffDiagMin, tFindCINatOrbs
logical :: tNoRenormRandExcits, tAssumeSizeExcitgen, tCycleOrbs, tROIteration
logical :: tShakeIter, tRotateOccOnly, tDoubExcMin, tUseHFOrbs, tRotateOrbs
logical :: tNonUniRandExcits, tNoSymGenRandExcits, tLagrange, tShakeApprox
logical :: tShake, tRotateVirtOnly, tMaxHLGap, tCacheFCIDUMPInts, tNoRODump
logical :: tKPntSym        !Are we using KPoint symmetry?
logical :: tRotatedOrbsReal     !This means we are reading in a complex FCIDUMP, but all 
                                !orbitals have been rotated to be real. This requires all
                                !kpoints to be at the gamma point or BZ boundary.
                                !At the reading in, the integrals will be converted to reals,
                                !but kpoint symmetry can still be used.

logical :: tRIIntegrals   ! Read in RI 2-e integrals from RIDUMP file
logical :: tStoreSpinOrbs ! This is set when the orbitals are stored in 
                          ! spin-orbital notation

logical :: tPickVirtUniform ! Use the 3rd generation, random excitation
                            ! generators, which pick pairs of virtual orbitals
                            ! at random (but uniformly)

logical :: tISKFuncs      ! Only for use in systems where the kpoint mesh has inversion symmetry,this ensures all
                          ! integrals are real.
logical :: tOddS_HPHF     !If this is true, and you are using HPHF, then it will converge onto an Odd S HPHF state.
integer :: iParity(5), nMaxX, nMaxY, nMaxZ, nMSH, coulDampOrb, elecPairs
integer :: roIterMax, iRanLuxLev, DiagMaxMinFac, OneElmaxMinFac, iState
integer :: iTiltX, iTiltY, nOccAlpha, nOccBeta, ShakeIterMax, ShakeStart
integer :: MaxMinFac, MaxABPairs
real*8 :: BOX, BOA, COA, fUEGRs, fRc, fCoul, OrbECutoff, UHUB, BHUB
real*8 :: Diagweight, OffDiagWeight, OrbEnMaxAlpha, Alpha, fCoulDampBeta
real*8 :: fCoulDampMu, TimeStep, ConvergedForce, ShakeConverged, UMATEps
real*8 :: OneElWeight


integer :: nEl             ! Number of (non-frozen) electrons in the system
integer :: Stot            ! Restrict S to Stot when using CSFs
integer :: LMS             ! Restrict determinants/CSFs to Ms == LMS
integer :: csf_trunc_level ! Max nopen for CSFs if tTruncateCSF enabled. Above
                           ! this, switch to using determinants.

! Inputs for CSFs
logical :: tCSFOld        ! Use (Alex's) old CSF code
logical :: tCSF           ! Use CSFs
logical :: tTruncateCSF   ! Use determinants not CSFs for nopen > 
                          ! csf_trunc_level

! Calculate size of FCI determinant space using MC
logical :: tMCSizeSpace 
integer*8 :: CalcDetPrint, CalcDetCycles   ! parameters

! Inputs for the UEG
logical :: tUEGTrueEnergies ! This is the logical for use of unscaled energies in the UEG calculation; will normally break spawning
logical :: tLatticeGens   ! Use new UEG excitation generators
logical :: tNoFailAb
logical :: tUEGOffset     ! Use twisted boundary conditions
real*8 :: k_offset(3)      ! UEG parameter for twist-averaging
logical :: tUEGSpecifyMomentum ! UEG parameter to allow specification of total momentum
integer :: k_momentum(3) ! UEG parameter for total momentum
logical :: tOrbECutoff ! Whether we're using a spherical cutoff in momentum space or not

! For the UEG, we damp the exchange interactions.
!    0 means none
!    1 means screened (using an erfc)
!    2 means hard spherical cut-off (at a distance Rc=ALAT(4))
!      [see JSS, ASA PRB 77, 193110 (2008)]
integer :: iPeriodicDampingType

! Used to be stored in Integrals
INTEGER :: ORBORDER(8,2)

logical :: tFixLz   ! Fix Lz of determinants in FCIMC. Lz == LzTot
integer :: LzTot    ! Total Ml quantum number of state to converge on.
integer :: iMaxLz   ! abs(maximum Ml basis function))


integer :: nIfP     ! Size appended to nIfD in CASSTAR calculations. Keeps
                    ! track of the parent determinant of spawned walkers.

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
!     ordered by spin-orbital index.  Arr(10,2) is the energy of the 10th
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

    ! These should really be in hist.F90, but we get circular dependencies
    ! These are bad.
    logical :: tHistSpinDist
    integer(n_int), allocatable :: ilut_spindist(:)
    integer :: hist_spin_dist_iter
    integer, allocatable :: nI_spindist(:)

! Operators for type(symmetry)
interface assignment (=)
    module procedure SymAssign
end interface
interface operator (.eq.)
    module procedure SymEq
end interface
interface operator (.ne.)
    module procedure SymNEq
end interface
interface operator (.gt.)
    module procedure SymGt
end interface
interface operator (.lt.)
    module procedure SymLt
end interface

contains
    ! Operations on type(symmetry)
    elemental subroutine SymAssign (lhs, rhs)
        type(Symmetry), intent(out) :: lhs
        type(Symmetry), intent(in) :: rhs
        lhs%S = rhs%S
    end subroutine
    elemental logical function SymEq (a, b)
        type(Symmetry), intent(in) :: a, b
        SymEq = a%S .eq. b%S
    end function
    elemental logical function SymNEq (a, b)
        type(Symmetry), intent(in) :: a, b
        SymNEq = a%S .ne. b%S
    end function
    elemental logical function SymGt (a, b)
        type(Symmetry), intent(in) :: a, b
        SymGt = a%S .gt. b%S
    end function
    elemental logical function SymLt (a, b)
        type(Symmetry), intent(in) :: a, b
        SymLt = a%S .lt. b%S
    end function

end module SystemData
