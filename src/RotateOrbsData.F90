MODULE RotateOrbsData
! This just contains those variables which need to be used in other modules - in particular NatOrbsMod.
    USE global_utilities
    use MemoryManager, only: TagIntType
    use constants, only: dp
    IMPLICIT NONE
    SAVE


    real(dp) , ALLOCATABLE :: CoeffT1(:,:)  ! This contains the transformation matrix which 
                                            !rotates the HF orbitals into their new basis.
    INTEGER(TagIntType) :: CoeffT1Tag
    INTEGER :: MemAllocRot

! These are the labelling arrays which allow us to separate the occupied and virtual 
!orbitals, mix all together, use spin or spatial 
! orbitals and maitain symmetry etc.
    INTEGER , ALLOCATABLE :: SymLabelList2_rot(:),SymLabelCounts2_rot(:,:), &
                             SymLabelListInv_rot(:),SymLabelList3_rot(:), &
                             SymOrbs_rot(:)
    INTEGER(TagIntType)  :: SymLabelList2_rotTag,SymLabelCounts2_rotTag, &
                            SymLabelListInv_rotTag,SymLabelList3_rotTag, &
                            SymOrbs_rotTag,EvaluesTruncTag
    real(dp) , ALLOCATABLE :: EvaluesTrunc(:)

! NoOrbs is either nBasis or SpatOrbs depending on whether we are using spin or spatial orbitals.    
    INTEGER :: NoOrbs,SpatOrbs,NoFrozenVirt,NoRotOrbs
    real(dp) :: TruncEval

    ! unit to write the TRANSFORM file to.
    integer :: transform_unit

    TYPE(timer) :: FillOneRDM_Time,FillMP2VDM_Time,DiagNatOrbMat_Time, &
                   OrderCoeff_Time,FillCoeff_Time

    LOGICAL :: tTurnStoreSpinOff

END MODULE RotateOrbsData
