MODULE SymExcitDataMod
    IMPLICIT NONE
    SAVE

    REAL*8 :: pDoubNew
    INTEGER , ALLOCATABLE :: SymLabelList2(:),SymLabelCounts2(:,:)
    INTEGER :: MaxABPairs
    LOGICAL :: tNoSingsPossible
    INTEGER :: ScratchSize          !This indicates the upper bound of the arrays needed for the excitation generation. The array bounds are ScratchSize.
    INTEGER , ALLOCATABLE :: OrbClassCount(:)  !This is set up in SpinOrbSymSetup, and is a default ClassCount excitation generator, from which it is then easier to set up the determinant specific ones.

END MODULE SymExcitDataMod
