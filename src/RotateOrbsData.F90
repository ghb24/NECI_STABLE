module RotateOrbsData

    ! This contains variables which need to be used in particular in NatOrbsMod.

    use global_utilities
    use MemoryManager, only: TagIntType
    use constants, only: dp

    implicit none
    save

    ! This contains the transformation matrix which rotates the HF orbitals
    ! into their new basis.
    real(dp), allocatable :: CoeffT1(:,:)  
    integer(TagIntType) :: CoeffT1Tag

    integer :: MemAllocRot

    ! These are the labelling arrays which allow us to separate the occupied and virtual 
    ! orbitals, mix all together, use spin or spatial 
    ! orbitals and maintain symmetry etc.
    integer :: NoSymLabelCounts
    integer, allocatable :: SymLabelList2_rot(:), SymLabelCounts2_rot(:,:)
    integer, allocatable :: SymLabelListInv_rot(:), SymLabelList3_rot(:)
    integer, allocatable :: SymOrbs_rot(:)

    integer(TagIntType) :: SymLabelList2_rotTag, SymLabelCounts2_rotTag
    integer(TagIntType) :: SymLabelListInv_rotTag, SymLabelList3_rotTag
    integer(TagIntType) :: SymOrbs_rotTag, EvaluesTruncTag

    real(dp), allocatable :: EvaluesTrunc(:)

    ! NoOrbs is either nBasis or SpatOrbs depending on whether we are using
    ! spin or spatial orbitals.
    integer :: NoOrbs, SpatOrbs, NoFrozenVirt, NoRotOrbs
    real(dp) :: TruncEval

    ! Unit to write the TRANSFORM file to.
    integer :: transform_unit

    type(timer) :: FillOneRDM_Time, FillMP2VDM_Time, DiagNatOrbMat_Time
    type(timer) :: OrderCoeff_Time, FillCoeff_Time

    logical :: tTurnStoreSpinOff

end module RotateOrbsData
