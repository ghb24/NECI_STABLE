MODULE SymExcitDataMod
    IMPLICIT NONE

    REAL*8 :: pDoubNew
    INTEGER , ALLOCATABLE :: SymLabelList2(:),SymLabelCounts2(:,:)
    INTEGER , SAVE :: ScratchSize          !This indicates the upper bound of the arrays needed for the excitation generation. The array bounds are ScratchSize.
    INTEGER , ALLOCATABLE :: OrbClassCount(:)  !This is set up in SpinOrbSymSetup, and is a default ClassCount excitation generator, from which it is then easier to set up the determinant specific ones.
    INTEGER , ALLOCATABLE :: kPointToBasisFn(:,:,:,:) !This is set up in SpinOrbSymSetup for the hubbard model, indicies are kx, ky, kz and a spin index value
    INTEGER :: kTotal(3) !This is the total momentum of the reference configuration

    INTEGER , ALLOCATABLE :: SpinOrbSymLabel(:)        !Find symmetry label (for symexcit routines: 0 -> nSymLabels-1) from BasisFn
    INTEGER , ALLOCATABLE :: SymInvLabel(:)  !Find inverse symmetry label (0 -> nSymLabels-1)
    INTEGER , ALLOCATABLE :: SymTableLabels(:,:)    !Symmetry table for symexcit labels (not syms themselves)
    INTEGER , ALLOCATABLE :: KPntInvSymOrb(:)       !Lookup table to find the inverse-k point equivalent spin orbital

END MODULE SymExcitDataMod
