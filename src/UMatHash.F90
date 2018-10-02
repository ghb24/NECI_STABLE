module UMatHash
  use constants
  use SystemData, only: nel, nBasis
  use dSFMT_interface, only: genrand_real2_dSFMT
  use procedure_pointers, only: get_umat_el
  use sltcnd_mod, only: sltcnd_2
  use IntegralsData, only: UMat
  use util_mod_numerical
  implicit none

  ! the cumulative sum only contains the absolute values, it 
  ! is always real
  real(dp), allocatable :: CumSparseUMat(:), TotalCumWeights(:)
  ! copies of the same stuff for parallel spin excitations
  real(dp), allocatable :: CumSparseUMatPar(:), TotalCumWeightsPar(:)
  ! number of nonzero elements per pq (using pq as index) / number of large elements
  integer, allocatable :: nPQ(:), nLargePQ(:)
  integer, allocatable :: nPQPar(:), nLargePQPar(:)
  ! rs indices for each pq
  integer, allocatable :: rsPQ(:,:) 
  integer, allocatable :: rsPQPar(:,:) 
  ! position of first element of pq in CumSparseUMat and rspq
  integer, allocatable :: posPQ(:)
  integer, allocatable :: posPQPar(:)
  ! number of basis functions used for storage (number of (spin)-orbitals)
  integer :: nStoreBas
  integer :: numPQPairs
  integer :: numPQPairsPar

  contains 

    subroutine initializeSparseUMat()
      implicit none

      call setupCumulativeSparseUMat(CumSparseUMat,TotalCumWeights,nPQ,nLargePQ,rsPQ,&
           posPQ, numPQPairs, .false.)
      call setupCumulativeSparseUMat(CumSparseUMatPar, TotalCumWeightsPar, nPQPar, nLargePQPar, &
           rsPQPar, posPQPar, numPQPairsPar, .true.)
    end subroutine initializeSparseUMat

    !------------------------------------------------------------------------------------------!

    subroutine setupCumulativeSparseUMat(CumSparseUMatLoc,TotalCumWeightsLoc, &
         nPQLoc,nLargePQLoc,rsPQLoc,posPQLoc,numPQPairsLoc,tPar)
      implicit none
      real(dp), intent(out), allocatable :: CumSparseUMatLoc(:), TotalCumWeightsLoc(:)
      integer, intent(out), allocatable :: nPQLoc(:), nLargePQLoc(:), rsPQLoc(:)
      integer, intent(out), allocatable :: posPQLoc(:), numPQPairsLoc
      logical, intent(in) :: tPar      
      integer :: j, umatsize, nnz
      integer :: p,q,r,s,ex(2,2),pq
      ! current matrix element
      real(dp) :: matel

      call GetUMatSize(nBasis,nel,umatsize)

      ! number of pq-values
      numPQPairsLoc = fuseIndex(nStoreBasis,nStoreBasis-1)

      allocate(nPQLoc(numPQPairsLoc))
      allocate(posPQLoc(numPQPairsLoc))
      allocate(TotalCumWeightsLoc(numPQPairsLoc))

      ! count the number of nonzero elements in umat
      nnz = nnz + 1
      do j = 1, umatsize
         if(abs(UMat(j)) > eps) nnz = nnz + 1
      end do
      allocate(rsPQLoc(2,nnz))
      allocate(CumSparseUMatLoc(nnz))

      nnz = 0
      do pq = 1, numPQPairsLoc
         ! recover the indices p,q from pq
         call splitIndex(pq,p,q)
         ! store the starting position of the values for this pq
         posPQLoc(pq) = nnz + 1
         do r = 1, nStoreBasis
            do s = r+1, nStoreBasis
               if(tPar) then
                  matel = abs(UMat(UMatInd(p,q,r,s)) - UMat(UMatInd(p,q,s,r)))
               else
                  matel = abs(UMat(UMatInd(p,q,r,s)))
               endif
               if(matel > eps) then
                  nnz = nnz + 1
                  CumSparseUMatLoc(nnz) = matel
                  rsPQLoc(1,nnz) = r
                  rsPQLoc(2,nnz) = s
               endif
            end do
         end do

         ! number of matrix elements for this pq
         if(r == 1 .and. p==2) then
            nPQLoc(pq) = nnz
         else
            nPQLoc(pq) = nnz - posPQLoc(pq-1) - 1
         endif
         ! so far, we got the unordered matrix elements for pq
         ! with their respective r,s
         ! now, sort the entries + rsPQ entries
         ! remember to only act on the slice belonging to the current PQ
         ! last index belonging to this pq is nnz
         call sortParallel(CumSparseUMatLoc(posPQLoc(pq):nnz),&
              rsPQLoc(posPQLoc(pq):endThisPQ),nPQLoc(pq))

         ! accumulate the values of CumSparseUMat
         call accumulateValues(CumSparseUMatLoc(posPQLoc(pq):nnz))

         nLargePQLoc(pq) = binary_search_first_ge(CumSparseUMatLoc(posPQLoc(pq):nnz),&
              0.5*CumSparseUMatLoc(nnz))

         TotalCumWeightsLoc(pq) = CumSparseUMatLoc(nnz)
      end do

    end subroutine setupCumulativeSparseUMat

    !------------------------------------------------------------------------------------------!

    subroutine sortParallel(valArray, indexArray, nelems)
      ! sorts the arrays valArray and indexArray simultaneously using the comparison
      ! of valArray
      implicit none
       ! we require valArray and indexArray to have the same number of elements
      integer, intent(in) :: nelems
      real(dp), intent(inout) :: valArray(nelems)
      integer, intent(inout) :: indexArray(2,nelems)

      integer :: i
      integer, allocatable :: auxIndices(:)

      allocate(auxIndices(nelems))

      do i = 1, nelems
         auxIndices(i) = i
      end do

      ! sort the auxiliary indices using the comparison operator of valArray
      ! (sort in ascending order)
      call sort(auxIndices,matel_gt,matel_lt)
      
      ! reorder valArray and indexArray
      valArray(1:nelems) = valArray(auxIndices)
      indexArray(:,1:nelems) = indexArray(:,auxIndices)

      deallocate(auxIndices)

      contains 
        
        pure function matel_gt(mI, mJ) result(bGt)
          integer, intent(in) :: mI,mJ
          logical :: bGt
          
          bGt = CumSparseUMat(mI) > CumSparseUMat(mJ)
        end function matel_gt

        pure function matel_lt(mI, mJ) result(bGt)
          integer, intent(in) :: mI,mJ
          logical :: bGt
          
          bGt = CumSparseUMat(mI) < CumSparseUMat(mJ)
        end function matel_gt

    end subroutine sortParallel

    !------------------------------------------------------------------------------------------!

    ! convert an array of values into an array of accumulated sums of those values
    subroutine accumulateValues(valArray)
      real(dp), intent(inout) :: valArray(:)
      integer, intent(in) :: nelems
      
      integer :: i

      ! add up the elements of valArray
      do i = 2, size(valArray)
         valArray(i) = valArray(i) + valArray(i-1)
      end do
    end subroutine accumulateValues

    !-----------------------------------------------------------------------------------------!-

    subroutine freeCumulativeSparseUMat()
      implicit none

      if(allocated(CumSparseUMat)) deallocate(CumSparseUMat)
      if(allocated(TotalCumWeights)) deallocate(TotalCumWeights)
      if(allocated(posPQ)) deallocate(posPQ)
      if(allocated(rsPQ)) deallocate(rsPQ)
      if(allocated(nPQ)) deallocate(nPQ)

      if(allocated(CumSparseUMatPar)) deallocate(CumSparseUMatPar)
      if(allocated(TotalCumWeightsPar)) deallocate(TotalCumWeightsPar)
      if(allocated(posPQPar)) deallocate(posPQPar)
      if(allocated(rsPQPar)) deallocate(rsPQPar)
      if(allocated(nPQPar)) deallocate(nPQPar)
           
    end subroutine freeCumulativeSparseUMat

    !------------------------------------------------------------------------------------------!

    subroutine gen_excit_hel_cached(nI, ilutI, nJ, ilutJ, exFlag, ic, ex, tParity, &
         pgen, helgen, store, part_type)
      implicit none
      
      integer, intent(in) :: nI(nel), exFlag
      integer(n_int), intent(in) :: ilutI(0:NIfTot)
      integer, intent(out) :: nJ(nel), ic, ex(2,2)
      logical, intent(out) :: tParity
      real(dp), intent(out) :: pGen
      HElement_t(dp), intent(out) :: HElGen
      type(excit_gen_store_type), intent(inout), target :: store
      integer, intent(in), optional :: part_type

      

    end subroutine gen_excit_hel_cached

    !------------------------------------------------------------------------------------------!

    ! may not be needed, or neads tweaking since we need to guarantee that
    ! p,q are occupied
    subroutine selectPQFromCSUM(p,q,tLarge)
      implicit none
      integer, intent(out) :: p,q
      logical, intent(in) :: tLarge

      real(dp) :: randThresh
      integer :: pos
      
      randThresh = genrand_real2_dSFMT()*TotalCumWeights(numPQPairs)
      pos = binary_search_first_ge(TotalCumWeights,randThresh)

      call splitIndex(pos,p,q)

    end subroutine selectPQFromCSUM

    !------------------------------------------------------------------------------------------!

    ! r,s, might be occupied and thus no valid excitation is generated
    ! we will answer this by generating multiple excitations, such that 1 valid one
    ! is generated on average
    subroutine selectRSFromCSUM(p,q,r,s,pgen,CumSparseUMatLoc,TotalCumWeightsLoc, &
         nPQLoc,rsPQLoc,posPQLoc)
      implicit none
      ! input p,q
      integer, intent(in) :: p,q
      ! output r,s
      integer, intent(out) :: r,s
      real(dp), intent(inout) :: pgen
      real(dp), intent(in) :: CumSparseUMatLoc(:), TotalCumWeightsLoc(:)
      integer, intent(in) :: nPQLoc(:), rsPQLoc(:), posPQLoc(:)

      real(dp) :: randThresh
      ! positions and aux indices
      integer :: startPQ, endPQ, pos, pq

      ! get the fused index
      pq = fuseIndex(p,q)
      ! and the related positions in the cache
      startPQ = posPQLoc(pq)
      endPQ = startPQ + nPQLoc(pq) - 1
      
      ! choose the HElement
      ! random factor between 0 and 1 times total cumulated value
      randThresh = genrand_real2_dSFMT()*CumSparseUMat(endPQ)
      
      pos = binary_search_first_ge(CumSparseUMat(startPQ:endPQ), randThresh) + startPQ
      
      if(pos > 1) then
         pgen = pgen * (CumSparseUMat(pos) - CumSparseUMat(pos-1))/CumSparseUMat(endPQ)
      else
         pgen = pgen * CumSparseUMat(pos)/CumSparseUMat(endPQ)
      end if

      r = rsPQ(1,pos)
      s = rsPQ(2,pos)
      
    end subroutine selectRSFromCSUM

    !------------------------------------------------------------------------------------------!

    function fuseIndex(q,p) result(ind)
      implicit none
      integer, intent(in) :: p,q
      integer :: ind

      ! qp and pq are considered to be the same index
      ! -> permutational symmetry

      if(p > q) then
         ind = q + p*(p-1)/2
      else
         ind = p + q*(q-1)/2
      end if
    end function fuseIndex

    !------------------------------------------------------------------------------------------!

    subroutine splitIndex(ind,p,q)
      ! inversion of fuseIndex()
      implicit none
      integer, intent(in) :: ind
      integer, intent(out) :: p,q
      integer :: i, tmp
      ! additionally, p > q, so reorder if necessary
      
      tmp = 0
      do i = 1, nBasis
         if(ind > tmp) then
            exit
         end if
         tmp = tmp + i
      end do

      q = ind - tmp
      ! i has the value of the larger index
      p = i
    end subroutine splitIndex

end module UMatHash
