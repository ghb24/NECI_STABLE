! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
MODULE SymExcitDataMod
    use constants, only: dp
    IMPLICIT NONE
    SAVE

    real(dp) :: pDoubNew
    INTEGER , ALLOCATABLE :: SymLabelList2(:),SymLabelCounts2(:,:)

    ! What are the upper bounds for the scratch arrays neede for excitation
    ! generation.
    integer :: ScratchSize

    ! Often ScratchSize1-4 will == ScratchSize, but not necessarily.
    integer :: ScratchSize1 = 0, ScratchSize2 = 0, ScratchSize3 = 0

    ! Do we want to extract lists of occupied / unoccupied orbitals when
    ! decoding bit determinants in the main FCIMC loop?
    logical :: tBuildOccVirtList = .false.
!This is set up in SpinOrbSymSetup, and is a default ClassCount excitation generator, 
!from which it is then easier to set up the determinant specific ones.
    INTEGER , ALLOCATABLE :: OrbClassCount(:)  
    !This is set up in SpinOrbSymSetup for the hubbard model, indicies are kx, ky, kz and a spin index value
    INTEGER , ALLOCATABLE :: kPointToBasisFn(:,:,:,:) 
    INTEGER :: kTotal(3) !This is the total momentum of the reference configuration

    INTEGER , ALLOCATABLE :: SpinOrbSymLabel(:)        !Find symmetry label (for symexcit routines: 0 -> nSymLabels-1) from BasisFn
    INTEGER , ALLOCATABLE :: SymInvLabel(:)  !Find inverse symmetry label (0 -> nSymLabels-1)
    INTEGER , ALLOCATABLE :: SymTableLabels(:,:)    !Symmetry table for symexcit labels (not syms themselves)
    INTEGER , ALLOCATABLE :: KPntInvSymOrb(:)       !Lookup table to find the inverse-k point equivalent spin orbital

    integer, allocatable :: MomInvSymOrb(:) !This is for when using MI functions, and gives the momentum mapping between orbs.

      ! Excitation generator stored information
      ! --> Due to the allocatables, we can add as many things to here as
      !     we like without causing any problems :-).
    type excit_gen_store_type
        ! These next two need to be pointers, rather than allocatable, to
        ! work with the GetNextSpawner stuff in CCMC.F90. It isn't really
        ! desirable
        ! --> Must ensure that init/clean_excit_gen_store is NEVER called on
        !     a store object being used in that way.
        integer, pointer :: ClassCountOcc(:) => null()
        integer, pointer :: ClassCountUnocc(:) => null()
        integer, pointer :: scratch3(:) => null()
        integer, pointer :: occ_list(:,:) => null()
        integer, pointer :: virt_list(:,:) => null()
        logical :: tFilled
        integer, pointer :: dorder_i (:) => null()
        integer, pointer :: dorder_j (:) => null()
        integer :: nopen
    end type

END MODULE SymExcitDataMod
