module fiveIndexLMat
  ! this module takes care of those entries of the L-matrix that have
  ! at least one repeated index
  ! special treatment is only for performance reasons
  use SystemData, only: G1
  use UMatCache, only: gtID
  use tc_three_body_data, only: lMat_t, nBI, twoIndexSize, fourIndexSize, tSpinCorrelator,&
       tDampLMat, tDampKMat, lMat
  use procedure_pointers, only: lMatInd
  use LMat_aux, only: fuseIndex, intswap, pairSwap, dampLMatel
  use constants
  implicit none
  
contains

!------------------------------------------------------------------------------------------!
  
  function get_lmat_el_five_ind_sparse(a,b,i,j,n) result(matel)
    ! Input: a,b - indices of orbitals to excite to
    !        i,j - indices of orbitals to excite from
    !        n   - index of unchanged orbital
    ! Output: matel - matrix element of this excitation
    implicit none
    integer, intent(in) :: a,b,i,j,n
    HElement_t(dp) :: matel
    integer(int64) :: ida, idb, idi, idj, idn
    logical :: tSameSpin
    
    tSameSpin = .not. tSpinCorrelator .or. (G1(a)%MS==G1(b)%ms .and. G1(a)%MS==G1(n)%MS)
    ! convert to spatial orbs if required
    ida = gtID(a)
    idb = gtID(b)
    idi = gtID(i)
    idj = gtID(j)
    idn = gtID(n)

    matel = 0

    call addExchangeTerm(a,b,n,ida,idb,idn,3,1)
    call addExchangeTerm(b,n,a,idb,idn,ida,2,1)
    call addExchangeTerm(n,a,b,idn,ida,idb,1,1)
    call addExchangeTerm(b,a,n,idb,ida,idn,3,-1)
    call addExchangeTerm(a,n,b,ida,idn,idb,2,-1)
    call addExchangeTerm(n,b,a,idn,idb,ida,1,-1)    
    if(tDampKMat .and. .not. tDampLMat) call dampLMatel(i,j,n,matel)

  contains

    subroutine addExchangeTerm(p,q,r,idp,idq,idr,pos,sgn)
      implicit none
      integer, intent(in) :: p,q,r,sgn,pos
      integer(int64), intent(in) :: idp,idq,idr
      integer(int64) :: index, exch
      integer(int64), parameter :: largerEx(2) = (/4,5/)
      logical :: tSpinAllowed
      
      tSpinAllowed = G1(i)%MS .eq. G1(p)%MS .and. G1(j)%MS .eq. G1(q)%MS .and. &
           G1(r)%MS .eq. G1(n)%MS
      
      if(tSpinAllowed) then
         exch = pos
         if(pos.eq.1) then
            if(idp > idi) exch = largerEx(pos)
         else if(pos.eq.2) then
            if(idq > idj) exch = largerEx(pos)
         endif
         
         index = LMatFiveInd(idp,idq,idi,idj,idr,exch)

         if(tSameSpin) then
            matel = matel + sgn * real(lMatAccessFiveIndex(lMat,index),dp)
         endif            
      endif
    end subroutine addExchangeTerm    
    
  end function get_lmat_el_five_ind_sparse

!------------------------------------------------------------------------------------------!  

  pure function LMatFiveInd(a,b,i,j,n,exch) result(index)
    ! the indexing function for the five-index integrals (ie contracted 6-index ints)
    ! this is the index of the matrix element L_{ijk}^{abn} where n is one of ijk
    ! Input: a,b - orbitals to excite to
    !        i,j,k - orbitals to excite from
    !        n   - orbital that remains unchanged
    !        exch - which exchange term do we ask for - this is the position of the repeated index
    ! Here, we assume that the orbitals to excite to are sorted with the unchanged one being
    ! the last one
    implicit none
    integer(int64), intent(in) :: a,b,n,i,j,exch
    integer(int64) :: index
    integer(int64) :: ai, bj
    integer, parameter :: maxEx = 5

    ai = fuseIndex(a,i)
    bj = fuseIndex(b,j)

    ! the ordering of ai/bj is important
    index = exch + maxEx*(ai-1+twoIndexSize*(bj-1)+fourIndexSize*(n-1))
      
  end function LMatFiveInd

!------------------------------------------------------------------------------------------!

    pure function lMatFiveIndUnsymmetric(a,b,n,i,j,k) result(index)
      implicit none
      integer(int64), intent(in) :: a,b,n,i,j,k
      integer(int64) :: index

      integer(int64) :: pos
      integer(int64), parameter :: maxPos = 5

      if(a.eq.n) then
         pos = 1
      else if(b.eq.n) then
         pos = 2
      else if(i.eq.n) then
         pos = 3
      else if(j.eq.n) then
         pos = 4
      else if(k.eq.n) then
         pos =5
      endif

      ! linear, contiguous, unsymmetrix 5-index function
      index = pos + maxPos*( (a-1) + nBI*( (b-1) + nBI*( (i-1) + nBI*( (j-1) + nBI*(k-1) ) ) ) )
    end function lMatFiveIndUnsymmetric

!------------------------------------------------------------------------------------------!
    
!------------------------------------------------------------------------------------------!

    ! only assign if there is one repeated index
    subroutine assignFiveIndexElem(arr,val,indicesInp)
      implicit none
      HElement_t(dp), intent(inout) :: arr(:)
      HElement_t(dp), intent(in) :: val
      integer(int64), intent(in) :: indicesInp(6)
      integer(int64) :: index, indices(6)
      integer(int64) :: nRep(6), pRep, nFound, iRep, marked, iEx, nEx
      ! we pair the indices into pairs, 1-4, 2-5 and 3-6 are a pair each, as they
      ! can be swapped without changing the matrix element
      integer(int64), parameter :: pairInd(6) = (/4,5,6,1,2,3/)
      ! there is a sixth, auxiliary index, which is given by the position of the repeated index
      integer(int64) :: exch(5)
      ! the action of swapping two electrons on the exchange index
      integer(int64), parameter :: mapExch(5) = (/2,1,3,5,4/)
      integer(int64), parameter :: exchTable(5) = (/1,2,3,1,2/)
      ! we need to clarify which index of a pair is the duplicate one, therefore,
      ! the auxiliary exchange index contains that information:
      ! if its the larger one, the exchange index is 4/5, else its 1/2 (for 3, they are equal)
      integer(int64), parameter :: largerEx(3) = (/4,5,3/)

      ! the indices have to be brought into a determined ordering now
      ! first, sort the indices into the canonical ordering      
      nRep = 0
      nFound = 0
      ! search for all repeated index pairs
      do pRep = 1, size(indicesInp)-1
         if(any(indicesInp(pRep).eq.indicesInp(pRep+1:))) then
            ! distinguish between an index appearing twice and one appearing more often
            if(any(nRep(:nFound).eq.pRep)) then
               ! this index is already logged, only add one entry
               nFound = nFound + 1
            else
               ! this is a new repeated index
               nFound = nFound + 2
               nRep(nFound-1) = pRep
            endif
            ! store the position of the second index in this pair of repeated indices
            ! this sum just casts an array of size 1 to a scalar
            ! findloc returns with an offset of pRep, correct for this
            nRep(nFound) = pRep + sum(findloc(indicesInp(pRep+1:),indicesInp(pRep)))            
         end if
      end do

      ! we must go through all repeated indices and store them      
      do pRep = 1, nFound
         ! operate on a copy of the input indices
         indices = indicesInp
         ! consider the pRep-th repeated index
         marked = nRep(pRep)
         ! marked is now the index of the chosen repeated index
         ! make sure the marked index is an upper one (i.e. > 4)
         if(marked < 4) then
            call intswap(indices(marked),indices(pairInd(marked)))
            ! now, the marked index is its previous pair index
            marked = pairInd(marked)
         end if

         ! now, we have to move the marked index to position 6
         if(marked.ne.6) then
            call pairSwap(indices(marked),indices(pairInd(marked)),indices(6),indices(pairInd(6)))
            marked = 6
         endif

         ! now, select the auxiliary exchange index - its given by where the other repeated index is
         ! if there are multiple candidates, each one has to be assigned
         iEx = 0
         exch = 0
         do iRep = 1, size(indices)-1
            if(indices(iRep).eq.indices(marked)) then
               iEx = iEx + 1
               exch(iEx) = exchTable(iRep)
               ! if the paired index is smaller, set the corresponding exchange index
               if(indices(iRep) > indices(pairInd(iRep))) exch(iEx) = largerEx(exch(iEx))
            endif
         end do

         ! assign all valid exchange terms
         nEx = iEx
         do iEx = 1, nEx
            ! the indices have a determined ordering now

            ! assign all possible permutations with these indices and the given exchange index
            index = lMatFiveInd(indices(1),indices(2),indices(4),indices(5),indices(3),exch(iEx))
            arr(index) = val

            ! permuting two indices might change which exchange index we are looking at
            index = lMatFiveInd(indices(2),indices(1),indices(5),indices(4),&
                 indices(3),mapExch(exch(iEx)))
            arr(index) = val
         end do

      end do
    end subroutine assignFiveIndexElem
      
!------------------------------------------------------------------------------------------!

  subroutine initFiveIndexAccess()
    implicit none

    ! set the stride for the last dimension of the five-index index function
    twoIndexSize = fuseIndex(nBI,nBI)
    fourIndexSize = twoIndexSize**2
  end subroutine initFiveIndexAccess

!------------------------------------------------------------------------------------------!

  pure function lMatAccessFiveIndex(lMatObj, index) result(matel)
    ! return the entry of the five-index part of the full six-index LMat
    ! that has the given index
    ! Input: lMatObj - container of the six-index integrals
    !        index   - linear index of the desired entry
    ! Output: matel  - matrix element of lMatObj
    use tc_three_body_data, only: lMat_t
    implicit none
    type(lMat_t), intent(in) :: lMatObj
    integer(int64), intent(in) :: index
    HElement_t(dp) :: matel

    ! the five-index part is stored separately
    matel = lMatObj%fiveIndexPtr(index)
  end function lMatAccessFiveIndex

!------------------------------------------------------------------------------------------!
  
  pure function isFiveIndex(indices) result(fI)
    ! given five orbital indices, check if there is at least one repeated index
    implicit none
    integer, intent(in) :: indices(:)
    logical :: fI
    integer :: i, nPairs

    nPairs = 0
    do i = 1, size(indices)-1
       if(any(indices(i).eq.indices(i+1:))) then
          nPairs = nPairs + 1
       endif
    end do
    ! its five index if there is at least one repeated index (three repeated indices are boring)
    fI = nPairs > 0

  end function isFiveIndex

!------------------------------------------------------------------------------------------!
  
  pure subroutine orderRepeated(q,r)
    ! Given two repeated indices and their paired indices, re-order them such that the
    ! one paired with the bigger one is behind
    ! Input: p,q - Orbital indices of the first pair, with p being a repeated index
    !        r,s - Orbital indices of the second pair, with s being a repeated index
    ! On return, p,q is the pair with the larger non-repeated index and p is the
    ! repeated index
    implicit none
    integer(int64), intent(inout) :: q,r

    if(q < r) then
       ! exchange the pairs if they have the wrong order
       call intswap(q,r)
    endif
  end subroutine orderRepeated
  
end module fiveIndexLMat
