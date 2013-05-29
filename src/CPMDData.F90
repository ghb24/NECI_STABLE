! Copyright (c) 2013, Ali Alavi
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
module CPMDData

! Module for data required for the CPMD interface.
! Data passed from CPMD.
      use MemoryManager ,only: TagIntType
      use constants, only: dp,int64

      implicit none

      save

!     NSTATES:   number of states.
!     NCOMPLETE: number of states that are part of complete degenerate
!                sets.  Used only in CPMD routines. The highest energy
!                incomplete degenerate sets are automatically detected in the
!                gamma-point symmetry routines.  nStates is set to be nComplete.
!     NROTOP:    number of rotational symmetry operations.
!     NTRANSLAT: number of translational symmetry operations.
!     NSpGpOp:   number of space group operations of the cell.
!     NKPS: Number of k-points (excluding symmetry equivalent ones.)
!     NKPTRANS:  number of translational operations of supercell
!                equivalent to the k-point mesh.
!     CELLSIZES: dimensions of supercell in terms of the primitive cell.  Used
!                only in gamma-point calculations when treating translational
!                symmetry separately from rotational.
      INTEGER :: NSTATES,NROTOP,NTRANSLAT,NCOMPLETE,NKPS,NKPTRANS,NSpGpOp,CellSizes(3)

!     BIGCELL:   true if the cell is non-primitive.
!     PROJECT:   true if states are projected according to
!     translational symmetry.
!     TSPACEGP:  true if using all space group operations, rather
!                than considering rotational and translaitonal
!                symmetries separately.
!     TKP:       true if doing a k-point calculation rather than gamma-point.
      LOGICAL :: BIGCELL,PROJECT,TSPACEGP,TKP


!     XI:        finite size correction.  See paper by Fraser et al.  Currently
!                not used---superseded by periodic image integrals.
!     EIONION:   nuclear-nuclear interactions.
!     DEGENTOL:  Tolerance (in eV) below which eigenvalues are taken to be
!                degenerate.
      real(dp) :: XI,EIONION,DEGENTOL

!     ROT_CHAR:  characters of (the degenerate sets of) the
!                wavefunctions.
!     EIGENVALUES:eigenvalues of the states, in rydberg.
!     PIInt: Periodic image integrals, <ii|ii>_PI.
!            Interaction of a charge density, i*i with it's periodic images
!            but not itself. <ii|ii>_PI=<ii|ii>_space - <ii|ii>_cell
!                                      =~<ii|ii>_space - <ii|ii>_attenuated
!     TRANS_CHAR:characters of (the degenerate sets of) the
!                wavefunctions.  If the translational irreducible
!                representations are real, then these are the
!                representations of the (1D) symmetry states.
!     GROUP_CHAR:characters due to the operations of the full
!                space group (or at least a large subgroup thereof.
!     Trans_Char_Table: character table of translational symmetry group.
!                Used when considering rotational symmetry separately from
!                rotational symmetry.
!     KPnt_Char: Translational representation spanned by each k-point. Used when
!                considering a k-point mesh to be equivalent to an effective
!                supercell.  Legacy.  Superseded by directly working with k-point
!                vectors (i.e. directly conserving crystal momenta) with vector
!                propogators and Abelian symmetry.
!     IMPROPER_ROT: true if operation is improper.
!     K_VECTORS: Stores the momentum of each translational label
!     ROT/TRANS_LABEL : label characters of states by unique
!                representations.
!     KPNT_CHAR: Representation of a state with a given k-point
!                the translational symmetry operations of the
!                equivalent supercell.
!     KPNTIND:   K-point index of each state (using a combined index).
      real(dp), ALLOCATABLE :: ROT_CHAR(:,:)     ! shape: (NSTATES,NROTOP)
      real(dp), ALLOCATABLE :: EIGENVALUES(:)    ! shape: (NSTATES)
      real(dp), ALLOCATABLE :: PIInt(:)          ! shape: (NStates)
      complex(dp), ALLOCATABLE :: TRANS_CHAR(:,:)       ! shape: (NSTATES,NTRANSLAT)
      complex(dp), ALLOCATABLE :: GROUP_CHAR(:,:)       ! shape: (nSpGpOp,NSTATES)
      complex(dp), ALLOCATABLE :: TRANS_CHAR_TABLE(:,:) ! shape: (NTRANSLAT,NTRANSLAT)
      complex(dp), ALLOCATABLE :: KPNT_CHAR(:,:)        ! shape: (NKPTRANS,NKPS)
      LOGICAL, ALLOCATABLE :: IMPROPER_OP(:)    ! shape: (NROTOP*NTRANSLAT)
      INTEGER, ALLOCATABLE :: K_VECTORS(:,:)    ! shape: (3,NTRANSLAT)
      INTEGER, ALLOCATABLE :: ROT_LABEL(:)      ! shape: (NSTATES)
      INTEGER, ALLOCATABLE :: TRANS_LABEL(:)    ! shape: (NSTATES)
      INTEGER, ALLOCATABLE :: KPNTIND(:)        ! shape: (2*nkps*nStates)

      ! Memory logging tags.
      INTEGER(TagIntType) :: tagROT_CHAR=0
      INTEGER(TagIntType) :: tagEIGENVALUES=0
      INTEGER(TagIntType) :: tagPIInt=0
      INTEGER(TagIntType) :: tagTRANS_CHAR=0
      INTEGER(TagIntType) :: tagGROUP_CHAR=0
      INTEGER(TagIntType) :: tagTRANS_CHAR_TABLE=0
      INTEGER(TagIntType) :: tagKPNT_CHAR=0
      INTEGER(TagIntType) :: tagIMPROPER_OP=0
      INTEGER(TagIntType) :: tagK_VECTORS=0
      INTEGER(TagIntType) :: tagROT_LABEL=0
      INTEGER(TagIntType) :: tagTRANS_LABEL=0
      INTEGER(TagIntType) :: tagKPNTIND=0

end module CPMDData
