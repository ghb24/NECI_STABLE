MODULE RotateOrbsData
! This just contains those variables which need to be used in other modules - in particular NatOrbsMod.
    USE global_utilities
    IMPLICIT NONE
    SAVE


    REAL*8 , ALLOCATABLE :: CoeffT1(:,:)  ! This contains the transformation matrix which rotates the HF orbitals into their new basis.
    INTEGER :: CoeffT1Tag,MemAllocRot

! These are the labelling arrays which allow us to separate the occupied and virtual orbitals, mix all together, use spin or spatial 
! orbitals and maitain symmetry etc.
    INTEGER , ALLOCATABLE :: SymLabelList2(:),SymLabelCounts2(:,:),SymLabelListInv(:),SymLabelList3(:),SymOrbs(:)
    INTEGER :: SymLabelList2Tag,SymLabelCounts2Tag,SymLabelListInvTag,SymLabelList3Tag,SymOrbsTag

! NoOrbs is either nBasis or SpatOrbs depending on whether we are using spin or spatial orbitals.    
    INTEGER :: NoOrbs,SpatOrbs,NoFrozenVirt

    TYPE(timer) , save :: FillOneRDM_Time,FillMP2VDM_Time,DiagNatOrbMat_Time,OrderCoeff_Time,FillCoeff_Time

END MODULE RotateOrbsData
