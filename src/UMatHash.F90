module UMatHash
  use constants
  use SystemData, only: nel, nBasis
  use dSFMT_interface, only: genrand_real2_dSFMT
  use procedure_pointers, only: get_umat_el
  use sltcnd_mod, only: sltcnd_2
  use IntegralsData, only: UMat
  use util_mod_numerical
  use shared_memory_mpi
  use UMatCache, only: UMatInd, getUMatSize, numBasisIndices
  use sort_mod
  implicit none

  ! the cumulative sum only contains the absolute values, it 
  ! is always real
  real(dp), pointer :: CumSparseUMat(:), TotalCumWeights(:)
  ! copies of the same stuff for parallel spin excitations
  real(dp), pointer :: CumSparseUMatPar(:), TotalCumWeightsPar(:)
  ! number of nonzero elements per pq (using pq as index) / number of large elements
  integer, pointer :: nPQ(:), nLargePQ(:)
  integer, pointer :: nPQPar(:), nLargePQPar(:)
  ! rs indices for each pq
  integer, pointer :: rsPQ(:,:) 
  integer, pointer :: rsPQPar(:,:) 
  ! position of first element of pq in CumSparseUMat and rspq
  integer, pointer :: posPQ(:)
  integer, pointer :: posPQPar(:)
  ! number of basis functions used for storage (number of (spin)-orbitals)
  integer :: nStoreBasis
  integer :: numPQPairs
  integer :: numPQPairsPar
  ! mpi windows for shared memory allocation
  integer, parameter :: nShmWins = 6 ! number of shared memory windows
  integer(MPIArg) :: PQshmWins(nShmWins), PQshmWinsPar(nShmWins)

  contains 

    subroutine initializeSparseUMat()
      implicit none

      nStoreBasis = numBasisIndices(nBasis)

      call setupCumulativeSparseUMat(CumSparseUMat,TotalCumWeights,nPQ,nLargePQ,rsPQ,&
           posPQ, numPQPairs, PQShmWins, .false.)
      call setupCumulativeSparseUMat(CumSparseUMatPar, TotalCumWeightsPar, nPQPar, nLargePQPar, &
           rsPQPar, posPQPar, numPQPairsPar, PQshmWinsPar, .true.)

    end subroutine initializeSparseUMat

    !------------------------------------------------------------------------------------------!

    subroutine setupCumulativeSparseUMat(CumSparseUMatLoc,TotalCumWeightsLoc, &
         nPQLoc,nLargePQLoc,rsPQLoc,posPQLoc,numPQPairsLoc, shmWins, tPar)
      implicit none
      real(dp), intent(out), pointer :: CumSparseUMatLoc(:), TotalCumWeightsLoc(:)
      integer, intent(out), pointer :: nPQLoc(:), nLargePQLoc(:), rsPQLoc(:,:)
      integer, intent(out), pointer :: posPQLoc(:)
      integer, intent(out) :: numPQPairsLoc
      integer(MPIArg), intent(out) :: shmWins(nShmWins)
      logical, intent(in) :: tPar      
      integer(int64) :: nnz
      integer :: p,q,r,s,pq
      integer(int64) :: nPQP_64
      ! current matrix element
      real(dp) :: matel

      ! temporary counters for aligning the sub-arrays of CumSparseUMatLoc
      integer :: nnzPQ
     
      ! number of pq-values
      numPQPairsLoc = fuseIndex(nStoreBasis,nStoreBasis)
      ! 64-bit version
      nPQP_64 = int(numPQPairsLoc,int64)

      call shared_allocate_mpi(shmWins(1),nPQLoc,(/nPQP_64/))
      call shared_allocate_mpi(shmWins(2),nLargePQLoc,(/nPQP_64/))
      call shared_allocate_mpi(shmWins(3),posPQLoc,(/nPQP_64/))
      call shared_allocate_mpi(shmWins(4),TotalCumWeightsLoc,(/nPQP_64/))

      ! count the number of nonzero elements in umat
      nnz = 0
      do pq = 1, numPQPairsLoc
         ! add empty entries to align to 64 byte
         call roundTo64(nnz)
         ! get the indices p,q
         call splitIndex(pq,p,q)
         do r = 1, nStoreBasis
            do s = minInd(r), nStoreBasis
               if(umatel(p,q,r,s) > eps) nnz = nnz + 1
            end do
         end do
      end do
      call shared_allocate_mpi(shmWins(5),rsPQLoc,(/2_int64,nnz/))
      call shared_allocate_mpi(shmWins(6),CumSparseUMatLoc,(/nnz/))
      
      allocate(rsPQLoc(2,nnz))
      allocate(CumSparseUMatLoc(nnz))

      nnz = 0
      do pq = 1, numPQPairsLoc
         ! align to 64 byte
         call roundTo64(nnz)
         ! recover the indices p,q from pq
         call splitIndex(pq,p,q)         
         ! store the starting position of the values for this pq
         posPQLoc(pq) = nnz + 1

         do r = 1, nStoreBasis
            ! as we consider spatial orbs in general, we allow for s==r
            do s = minInd(r), nStoreBasis
               matel = umatel(p,q,r,s)
               if(matel > eps) then
                  nnz = nnz + 1
                  CumSparseUMatLoc(nnz) = matel
                  rsPQLoc(1,nnz) = r
                  rsPQLoc(2,nnz) = s
               endif
            end do
         end do
        
         ! number of matrix elements for this pq
         nPQLoc(pq) = nnz - posPQLoc(pq) + 1
         ! so far, we got the unordered matrix elements for pq
         ! with their respective r,s
         ! now, sort the entries + rsPQ entries
         ! remember to only act on the slice belonging to the current PQ
         ! last index belonging to this pq is nnz
         if(nPQLoc(pq) > 0) then
            call sort(CumSparseUMatLoc(posPQLoc(pq):nnz),rsPQLoc(:,posPQLoc(pq):nnz))

            ! ordering is increasing now, but we want it to be decreasing
            CumSparseUMatLoc(posPQLoc(pq):nnz) = CumSparseUMatLoc(nnz:posPQLoc(pq):-1)
            rsPQLoc(:,posPQLoc(pq):nnz) = rsPQLoc(:,nnz:posPQLoc(pq):-1)

            ! accumulate the values of CumSparseUMat
            call accumulateValues(CumSparseUMatLoc(posPQLoc(pq):nnz))

            nLargePQLoc(pq) = binary_search_first_ge(CumSparseUMatLoc(posPQLoc(pq):nnz),&
                 0.5*CumSparseUMatLoc(nnz))

            TotalCumWeightsLoc(pq) = CumSparseUMatLoc(nnz)
         else
            ! if there are no excitations with this pq, posPQ(pq) is set to the position
            ! of the next pq with valid excitations and the number nPQ/nPQLarge is set to 0
            nPQLoc(pq) = 0
            nLargePQLoc(pq) = 0
            TotalCumWeightsLoc(pq) = 0.0_dp
         end if
      end do

    contains

      function umatel(p,q,r,s) result(matel)
        integer, intent(in) :: p,q,r,s
        real(dp) :: matel

        if(tPar) then
           matel = abs(UMat(UMatInd(p,q,r,s)) - UMat(UMatInd(p,q,s,r)))
        else
           matel = abs(UMat(UMatInd(p,q,r,s)))
        endif
      end function umatel

      function minInd(ind) result(lowerBound)
        ! for parallel excitations, rs are ordered, for antiparallel, they are
        ! not, since r and p and s and q have to have the same spin
        implicit none
        integer, intent(in) :: ind
        integer :: lowerBound
        
        if(tPar) then
           lowerBound = ind
        else
           lowerBound = 1
        end if
      end function minInd

      subroutine roundTo64(pos)
        implicit none
        integer(int64), intent(inout) :: pos

        if(mod(pos,cLineSize) .ne. 0) pos = pos + cLineSize - mod(pos,cLineSize)
      end subroutine roundTo64

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
      !call sort(auxIndices,matel_gt,matel_lt)
      
      ! reorder valArray and indexArray
      valArray(1:nelems) = valArray(auxIndices)
      indexArray(:,1:nelems) = indexArray(:,auxIndices)

      deallocate(auxIndices)

      contains 
        
        pure function matel_gt(mI, mJ) result(bGt)
          use constants
          integer, intent(in) :: mI,mJ
          logical :: bGt
          
          bGt = CumSparseUMat(mI) > CumSparseUMat(mJ)
        end function matel_gt

        pure function matel_lt(mI, mJ) result(bGt)
          use constants
          integer, intent(in) :: mI,mJ
          logical :: bGt
          
          bGt = CumSparseUMat(mI) < CumSparseUMat(mJ)
        end function matel_lt

    end subroutine sortParallel

    !------------------------------------------------------------------------------------------!

    ! convert an array of values into an array of accumulated sums of those values
    subroutine accumulateValues(valArray)
      real(dp), intent(inout) :: valArray(:)
      
      integer :: i

      ! add up the elements of valArray
      do i = 2, size(valArray)
         valArray(i) = valArray(i) + valArray(i-1)
      end do
    end subroutine accumulateValues

    !-----------------------------------------------------------------------------------------!-

    subroutine freeCumulativeSparseUMat()
      implicit none

      if(associated(CumSparseUMat)) deallocate(CumSparseUMat)
      if(associated(TotalCumWeights)) deallocate(TotalCumWeights)
      if(associated(posPQ)) deallocate(posPQ)
      if(associated(rsPQ)) deallocate(rsPQ)
      if(associated(nPQ)) deallocate(nPQ)

      if(associated(CumSparseUMatPar)) deallocate(CumSparseUMatPar)
      if(associated(TotalCumWeightsPar)) deallocate(TotalCumWeightsPar)
      if(associated(posPQPar)) deallocate(posPQPar)
      if(associated(rsPQPar)) deallocate(rsPQPar)
      if(associated(nPQPar)) deallocate(nPQPar)
           
    end subroutine freeCumulativeSparseUMat

    !------------------------------------------------------------------------------------------!
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
    subroutine selectRSFromCSUM(p,q,rs,pgen,CumSparseUMatLoc, &
         nPQLoc,rsPQLoc,posPQLoc)
      implicit none
      ! input p,q
      integer, intent(in) :: p,q
      ! output r,s
      integer, intent(out) :: rs(2)
      real(dp), intent(inout) :: pgen
      real(dp), intent(in) :: CumSparseUMatLoc(:)
      integer, intent(in) :: nPQLoc(:), rsPQLoc(:,:), posPQLoc(:)

      real(dp) :: randThresh
      ! positions and aux indices
      integer :: startPQ, endPQ, pos, pq

      ! get the fused index
      pq = fuseIndex(p,q)
      ! and the related positions in the cache
      startPQ = posPQLoc(pq)
      endPQ = startPQ + nPQLoc(pq) - 1
      if(endPQ > startPQ) then
         ! choose the HElement
         ! random factor between 0 and 1 times total cumulated value
         randThresh = genrand_real2_dSFMT()*CumSparseUMatLoc(endPQ)

         pos = linearSearch_old(CumSparseUMatLoc(startPQ:endPQ), randThresh) + startPQ - 1

         if(pos > startPQ) then
            pgen = pgen * (CumSparseUMatLoc(pos) - CumSparseUMatLoc(pos-1))/CumSparseUMatLoc(endPQ)
         else
            pgen = pgen * CumSparseUMatLoc(pos)/CumSparseUMatLoc(endPQ)
         end if

         rs = rsPQLoc(:,pos)
      else
         ! no non-zero matrix elements for this pq
         rs = 0
      end if
      
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
      ! additionally, q > p, so reorder if necessary
      
      tmp = 0
      do i = 1, nStoreBasis
         tmp = tmp + i
         if(ind <= tmp) exit
      end do

      p = ind - tmp + i
      ! i has the value of the larger index
      q = i
    end subroutine splitIndex

    !------------------------------------------------------------------------------------------!

    function linearSearch(arr, thresh) result(ind)
      implicit none
      real(dp), intent(in) :: arr(:)
      real(dp), intent(in) :: thresh
      integer :: ind

      ind = 1
      do 
         if(arr(ind) > thresh) return
         ind = ind + 1
      end do

    end function linearSearch

    !------------------------------------------------------------------------------------------!

    function linearSearch_old(arr, thresh) result(ind)
      implicit none
      real(dp), intent(in) :: arr(:)
      real(dp), intent(in) :: thresh
      integer :: ind
      
      ! left and right point of extrapolation
      integer :: a, b
      real(dp) :: fa, fb, fx
      real(dp), parameter :: scale = 0.5
      
      a = 0
      b = size(arr)
      ! for 1-sized arrays, the search does not work
      if(size(arr) < 2) then
         ind = 1
         return
      endif
      fa = -thresh
      fb = arr(b) - thresh
      
      do 
         if(abs(a-b)==1) then
            ind = max(a,b)
            return
         endif
         ind = ceiling(a - (fa)/(fb-fa)*(b-a))
         fx = arr(ind) - thresh
         if(fx*fb < 0) then
            a = b
            fa = fb
         else
            fa = fa * scale
         endif
         b = ind
         fb = fx
      end do
    end function linearSearch_old

    !------------------------------------------------------------------------------------------!

    function cacheLineSearch(arr, thresh) result(ind)
      implicit none
      real(dp), intent(in) :: arr(:)
      real(dp), intent(in) :: thresh
      integer :: ind

      integer :: i, iStartLinSearch, modulus, nReads
      real(dp) :: fb, fa

      iStartLinSearch = 1
      nReads = cLineSize - 1
      do 
         ! do a secant search using cache-line length blocks of data
         ind = iStartLinSearch
         fa = arr(ind)
         if(fa > thresh) return
         do i = 1, nReads
            ind = i + iStartLinSearch
            fb = arr(ind)
            if(fb > thresh) return
         end do
         iStartLinSearch = (thresh - fa)/(fb-fa)*(i-1) + iStartLinSearch
         
         ! align to 64 byte
         modulus = mod(iStartLinSearch-1,cLineSize)
         if(modulus .ne. 0) then
            nReads = 7 - modulus
            if(nReads .eq. 0) nReads = cLineSize
         endif
      end do
   

    end function cacheLineSearch

end module UMatHash
