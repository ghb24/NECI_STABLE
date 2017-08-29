module SystemData

    use constants, only: n_int,int64,dp
    use MemoryManager, only: TagIntType

implicit none

save

integer :: symmax ! Max number of irreps to deal with. Value computed in GETFCIBASIS readint.f90
logical :: tMolpro,tMolcas  !True if the code has been called from Molpro or Molcas
logical :: tMolproMimic !True if the code is being run from standalone neci, but designed to mimic the runtime 
                        !behaviour of molpro
character(12) :: MolproID

logical :: tNoSingExcits    !True if there are no single excitations in the system

logical :: tStarBin, tReadInt, tHFOrder, tDFRead, tPBC, tUEG, tUEG2, tCPMD, tHUB, tHeisenberg
logical :: tHPHF, tHPHFInts, tUHF, tSPN, tParity, tUseBrillouin, tExch, tReal
logical :: tTilt, tOneElIntMax, tOnePartOrbEnMax, tROHF, tBrillouinsDefault
logical :: tNoBrillouin, tVirtCoulombMax, tVirtExchangeMin, tHijSqrdMin
logical :: tDiagonalizehij, tHFSingDoubExcMax, tSpinOrbs, tReadInCoeff
logical :: tUseMP2VarDenMat, tAlpha, tStoreAsExcitations, tBin
logical :: tVASP, tOffDiagSqrdMin, tOffDiagSqrdMax, tOffDiagmax, tShakeDelay
logical :: tSeparateOccVirt, tExactSizeSpace, tRotatedOrbs
logical :: tImportanceSample, tERLocalization, tOffDiagMin, tFindCINatOrbs
logical :: tNoRenormRandExcits, tAssumeSizeExcitgen, tCycleOrbs, tROIteration
logical :: tShakeIter, tRotateOccOnly, tDoubExcMin, tUseHFOrbs, tRotateOrbs
logical :: tNonUniRandExcits, tNoSymGenRandExcits, tLagrange, tShakeApprox
logical :: tShake, tRotateVirtOnly, tMaxHLGap, tCacheFCIDUMPInts
logical :: tKPntSym        !Are we using KPoint symmetry?
logical :: tRotatedOrbsReal     !This means we are reading in a complex FCIDUMP, but all 
                                !orbitals have been rotated to be real. This requires all
                                !kpoints to be at the gamma point or BZ boundary.
                                !At the reading in, the integrals will be converted to reals,
                                !but kpoint symmetry can still be used.
logical :: tReadFreeFormat,tReltvy

logical :: tRIIntegrals   ! Read in RI 2-e integrals from RIDUMP file
logical :: tStoreSpinOrbs ! This is set when the orbitals are stored in 
                          ! spin-orbital notation

logical :: tPickVirtUniform ! Use the 3rd generation, random excitation
                            ! generators, which pick pairs of virtual orbitals
                            ! at random (but uniformly)

logical :: tISKFuncs      ! Only for use in systems where the kpoint mesh has inversion symmetry,this ensures all
                          ! integrals are real.
logical :: tOddS_HPHF     !If this is true, and you are using HPHF, then it will converge onto an Odd S HPHF state.
logical :: tAntisym_MI    !Antisymmetric MI functions.
logical :: tComplexOrbs_RealInts    !We are using complex orbitals, but real integrals. 
                                    !Therefore, check the mom sym before looking up integral, 
                                    !since we only have 4x perm sym.
logical :: tComplexWalkers_RealInts !We are using complex walkers, but real integrals and real
                                    !basis functions, such that there is still 8x perm sym.
integer :: iParity(5), nMaxX, nMaxY, nMaxZ, nMSH, coulDampOrb, elecPairs
integer :: roIterMax, iRanLuxLev, DiagMaxMinFac, OneElmaxMinFac, iState
integer :: iTiltX, iTiltY, nOccAlpha, nOccBeta, ShakeIterMax, ShakeStart
integer :: MaxMinFac, MaxABPairs
real(dp) :: BOX, BOA, COA, fUEGRs, fRc, OrbECutoff, UHUB, BHUB, btHub
real(dp) :: Diagweight, OffDiagWeight, OrbEnMaxAlpha, Alpha, fCoulDampBeta
real(dp) :: fCoulDampMu, TimeStep, ConvergedForce, ShakeConverged, UMATEps
real(dp) :: OneElWeight

integer :: AA_elec_pairs, BB_elec_pairs, par_elec_pairs, AB_elec_pairs
integer :: AA_hole_pairs, BB_hole_pairs, par_hole_pairs, AB_hole_pairs
integer :: hole_pairs, nholes_a, nholes_b, nholes


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
logical :: tMCSizeSpace,tMCSizeTruncSpace
integer :: iMCCalcTruncLev
integer(int64) :: CalcDetPrint, CalcDetCycles   ! parameters

! Inputs for the UEG
logical :: tUEGTrueEnergies ! This is the logical for use of unscaled energies in the UEG calculation; will normally break spawning
logical :: tLatticeGens   ! Use new UEG excitation generators
logical :: tNoFailAb
logical :: tUEGOffset     ! Use twisted boundary conditions
real(dp) :: k_offset(3)      ! UEG parameter for twist-averaging
logical :: tUEGSpecifyMomentum ! UEG parameter to allow specification of total momentum
integer :: k_momentum(3) ! UEG parameter for total momentum
logical :: tOrbECutoff ! Whether we're using a spherical cutoff in momentum space or not
logical :: tgCutoff ! Whether we're using a spherical cutoff for the momentum transfer vector
real(dp) :: gCutoff ! Spherical cutoff for the momentum transfer vector
logical :: tMP2UEGRestrict ! Restricts the MP2 sum over a single electron pair, specified by: 
integer :: kiRestrict(3), kjRestrict(3) ! ki/kj pair
integer :: kiMsRestrict, kjMsRestrict ! and their spins
logical :: tMadelung ! turning on self-interaction term
real(dp) :: Madelung ! variable storage for self-interaction term
logical :: tUEGFreeze ! Freeze core electrons for the UEG, a crude hack for this to work-around freezing not working for UEG
real(dp) :: FreezeCutoff
logical :: tRef_Not_HF

! Inputs for the UEG2
character(len=3) :: real_lattice_type ! type of reciprocal lattice (eg. fcc, sc, bcc, hcp)
integer :: k_lattice_vectors(3,3)
real(dp) :: k_lattice_constant
integer, allocatable :: kvec(:,:)
integer ::  Highest_orb_index
integer :: dimen

! Matrix elements for the modified hubbard model (including breathing effect)
real(dp), allocatable :: breathingCont(:)
integer, allocatable :: momIndexTable(:,:,:,:)
logical :: tmodHub

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
   integer(int64) :: S
END TYPE

integer, PARAMETER :: SymmetrySize=2
integer, PARAMETER :: SymmetrySizeB=SymmetrySize*8
TYPE BasisFN
   TYPE(Symmetry) :: sym
   INTEGER :: k(3)
   INTEGER :: Ms
   INTEGER :: Ml            !This is the Ml symmetry of the orbital
   INTEGER :: Dummy         !Rather than use SEQUENCE which has caused endless bother...
END TYPE

! Empty basis function is used in many places.
! This is useful so if BasisFn changes, we don't have to go
! through the code and change the explicit null statements.
type(BasisFn) :: NullBasisFn=BasisFn(Symmetry(0),(/0,0,0/),0,0,0)

integer, PARAMETER :: BasisFNSize=SymmetrySize+5
integer, PARAMETER :: BasisFNSizeB=BasisFNSize*8


TYPE(BASISFN) :: SymRestrict
INTEGER :: nBasisMax(5,7)
real(dp) :: ALAT(5)
real(dp) :: ECore
INTEGER :: nBasis
integer :: nMax
integer :: nnr
integer :: nocc
real(dp) :: OMEGA
logical :: tSpinPolar
INTEGER :: iSpinSkip ! Often referred to as ISS.

!From Calc  
real(dp) :: Beta
        
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
real(dp), pointer :: Arr(:,:) 
INTEGER(TagIntType) :: tagArr

! Lists orbitals in energy order. i.e. Brr(1) is the lowest energy orbital
INTEGER, pointer :: BRR(:) 
INTEGER(TagIntType) :: tagBrr

Type(BasisFN), pointer :: G1(:)  ! Info about the basis functions.
INTEGER(TagIntType) :: tagG1

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

! Should we use |K| for FCIQMC?
logical :: modk_offdiag

! True if we are performing a calculation in all symmetry sectors at once.
! This is used in finite-temperature KP-FCIQMC calculations.
logical :: tAllSymSectors

    logical :: tGenHelWeighted, tGen_4ind_weighted, tGen_4ind_reverse
    logical :: tUEGNewGenerator, tGen_4ind_part_exact, tGen_4ind_lin_exact
    logical :: tGen_4ind_2, tGen_4ind_2_symmetric

! Are we using multiple replicas?
logical :: tMultiReplicas

! Has the user set the symmetry using the 'SYM' option?
logical :: tSymSet = .false.

logical :: tGiovannisBrokenInit

! twisted boundary implementation for the hubbard model: 
! use keyword twisted-bc [real, real] in System Block of input 
! twist value is in values of periodicity (2*pi/L) for cubic and (pi/L) for 
! tilted lattice
logical :: t_twisted_bc = .false. 
real(dp) :: twisted_bc(2) = 0.0_dp

! flags for the use of open boundary conditions in the real-space 
! hubbard model. 
! for the cubic lattice the can be set separately, for the tilted only 
! full open BC are implemented
logical :: t_open_bc_x = .false.
logical :: t_open_bc_y = .false.

! use an intermediate flag for a new implementation of the newest excitation
! generator
logical :: tGen_4ind_unbound = .false.

! also implement a next-nearest neighbor Hubbard model implementation: 
! for k-space hubbard, this only affects the diagonal part! 
real(dp) :: nn_bhub = 0.0_dp

! do a quick test with different weightings of picking orbitals (a) 
logical :: t_iiaa = .false., t_ratio = .false. 
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
