module SystemData

    use constants, only: n_int,int64,dp
    use MemoryManager, only: TagIntType

implicit none

save

logical :: tMolpro  !True if the code has been called from Molpro
logical :: tMolproMimic !True if the code is being run from standalone neci, but designed to mimic the runtime 
                        !behaviour of molpro

logical :: tNoSingExcits    !True if there are no single excitations in the system

logical :: tStarBin, tReadInt, tHFOrder, tDFRead, tPBC, tUEG, tUEG2, tCPMD, tHUB
logical :: tHPHF, tHPHFInts, tUHF, tSPN, tParity, tUseBrillouin, tExch, tReal
logical :: tTilt, tUmatEps, tOneElIntMax, tOnePartOrbEnMax, tROHF, tBrillouinsDefault
logical :: tNoBrillouin, tVirtCoulombMax, tVirtExchangeMin, tHijSqrdMin, tMomInv
logical :: tDiagonalizehij, tHFSingDoubExcMax, tSpinOrbs, tReadInCoeff
logical :: tUseMP2VarDenMat, tAlpha, tStoreAsExcitations, tBin, tStarStore
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
logical :: tReadFreeFormat

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
integer :: iParity(5), nMaxX, nMaxY, nMaxZ, nMSH, coulDampOrb, elecPairs
integer :: roIterMax, iRanLuxLev, DiagMaxMinFac, OneElmaxMinFac, iState
integer :: iTiltX, iTiltY, nOccAlpha, nOccBeta, ShakeIterMax, ShakeStart
integer :: MaxMinFac, MaxABPairs
real(dp) :: BOX, BOA, COA, fUEGRs, fRc, OrbECutoff, UHUB, BHUB
real(dp) :: Diagweight, OffDiagWeight, OrbEnMaxAlpha, Alpha, fCoulDampBeta
real(dp) :: fCoulDampMu, TimeStep, ConvergedForce, ShakeConverged, UMATEps
real(dp) :: OneElWeight


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

! Options relating to the semi-stochastic code.
logical :: tSemiStochastic ! Performing a semi-stochastic simulation if true.
! Options regarding splitting the space into core and non-core elements. Needed, for example when performing a
! semi-stochastic simulation, to specify the deterministic space.
logical :: tCSFCore ! Use CSFs for the core states.
logical :: tDeterminantCore ! Use determinants for the core states.
logical :: tOptimisedCore ! Generate an optimised deterministic space by diagonalising part of the space.
logical :: tDoublesCore ! Use single and double excitations for the core states.
logical :: tCASCore ! Use Determinants where orbitals within an active space can differ from the Hartree-Fock for core states.
logical :: tRASCore ! Use a RAS space for the core space (see ras.F90 for definition).
logical :: tLowECore ! Like the optimised core space, but instead of diagonalising the space each iteration to find which states to keep, we keep the states with the lowest energies.
! cas_determ_bitmask has all bits that refer to the active space set, and all other bits unset.
! cas_not_determ_bitmask is simply the result after the not operation is applied to cas_determ_bitmask.
integer(n_int), allocatable, dimension(:) :: cas_determ_bitmask
integer(n_int), allocatable, dimension(:) :: cas_determ_not_bitmask
! Bitmasks with all bits corresponding to orbitals in RAS1 and RAS3, repectively, set.
integer(n_int), allocatable, dimension(:) :: core_ras1_bitmask
integer(n_int), allocatable, dimension(:) :: core_ras3_bitmask
! When using a CAS deterministic space, these integers store the number of orbitals above and below the Fermi energy to
! include in the CAS active space (the occupied and virtual orbitals).
integer :: OccDetermCASOrbs
integer :: VirtDetermCASOrbs
! When using tOptimisedCore, this specifies the minimum amplitude that a basis state should have in the ground state
! (of the deterministic space generated) in order to be kept in the next deterministic space. The first component
! refers to the first iteration, the second component to the second iteration...
real(dp), allocatable, dimension(:) :: determ_space_cutoff_amp
! When using tOptimisedCore, this specifies how many basis states (of the deterministic space generated) should be kept
! in the next deterministic space. The first component refers to the first iteration, the second component to the second
! iteration... i.e. if determ_space_cutoff_num(1) = 5 then in the first iteration, the 5 basis states with the largest
! amplitudes in the ground state are kept.
integer, allocatable, dimension(:) :: determ_space_cutoff_num
! If this logical is true, then the cutoff criterion is done using the amplitude. If false, it is done using a fixed number.
logical :: tDetermAmplitudeCutoff
! When using the optimised core option for semi-stochastic simulations, this option specifies how many iterations of
! the generation procedure should be performed.
integer :: num_det_generation_loops
! If this is true then set a limit on the maximum deterministic space size.
logical :: tLimitDetermSpace
! This is maximum number of elements in the deterministic space, if tLimitDetermSpace is true.
integer :: max_determ_size
! This option gives the maximum excitation level to go up to when generating the low energy deterministic space.
integer :: low_e_core_excit
! This integer specifies the number of states to keep for each iteration of the low energy core generation.
integer :: low_e_core_num_keep
! When using tLowECore, if this option is true then all doubles will be kept.
logical :: tLowECoreAllDoubles

! Options relating to the trial wavefunction.
logical :: tTrialWavefunction ! Use a trial wavefunction-based energy estimator.
logical :: tDoublesTrial ! Use single and double exciations for the trial space.
logical :: tCASTrial ! Use a CAS space for the trial space.
logical :: tOptimisedTrial ! Generate an optimised trial space by diagonalisaing part of the space.
! As for determ_space_cutoff_amp and determ_space_cutoff_num above, but the following two quantities refer to the trial space
! generation rather than the deterministic space generation.
logical :: tLowETrial ! Like the optimised trial space, but instead of diagonalising the space each iteration to find which states to keep, we keep the states with the lowest energies.
real(dp), allocatable, dimension(:) :: trial_space_cutoff_amp
integer, allocatable, dimension(:) :: trial_space_cutoff_num
! When using a CAS trial space, these integers store the number of orbitals above and below the Fermi energy to
! include in the CAS active space (the occupied and virtual orbitals).
integer :: OccTrialCASOrbs
integer :: VirtTrialCASOrbs
! If this logical is true, then the cutoff criterion is done using the amplitude. If false, it is done using a fixed number.
logical :: tTrialAmplitudeCutoff
! As for num_determ_generation_loops above, but here for the trial wavefunction generation.
integer :: num_trial_generation_loops
! If this is true then set a limit on the maximum trial space size.
logical :: tLimitTrialSpace
! This is maximum number of elements in the trial space, if tLimitDetermSpace is true.
integer :: max_trial_size
! This option gives the maximum excitation level to go up to when generating the low energy trial space.
integer :: low_e_trial_excit
! This integer specifies the number of states to keep for each iteration of the low energy trial generation.
integer :: low_e_trial_num_keep
! When using tLowETrial, if this option is true then all doubles will be kept.
logical :: tLowETrialAllDoubles

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
   SEQUENCE
   INTEGER :: k(3)
   INTEGER :: Ms
   INTEGER :: Ml            !This is the Ml symmetry of the orbital
   TYPE(Symmetry) :: sym
END TYPE

! Empty basis function is used in many places.
! This is useful so if BasisFn changes, we don't have to go
! through the code and change the explicit null statements.
type(BasisFn) :: NullBasisFn=BasisFn((/0,0,0/),0,0,Symmetry(0))

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
