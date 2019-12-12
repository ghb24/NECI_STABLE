module lMat_indexing
  use constants, only: int64
  use util_mod, only: intSwap, fuseIndex
  use SystemData, only: nBI
  implicit none

  ! for the symmetry broken index function
  integer :: strideInner, strideOuter

  contains

!------------------------------------------------------------------------------------------!
!  Index functions for the six-index addressing
!------------------------------------------------------------------------------------------!

    pure function lMatIndSym(a,b,c,i,j,k) result(index)
      ! Indexing function implementing 48-fold symmetry:
      ! Symmetric with respect to permutation of pairs (a,i), (b,j) and (c,k) as well as
      ! with respect to exchange of a<->i, b<->j and c<->k
      ! Input: a,b,c - orbital indices of electrons
      !        i,j,k - orbital indices of holes
      ! Output: index - contiguous index I(a,b,c,i,j,k) with the aforementioned symmetry
      implicit none
      integer(int64), value, intent(in) :: a,b,c ! occupied orb indices
      integer(int64), value, intent(in) :: i,j,k ! unoccupied orb
      integer(int64) :: index

      integer(int64) :: ai,bj,ck

      ai = fuseIndex(a,i)
      bj = fuseIndex(b,j)
      ck = fuseIndex(c,k)

      ! sort the indices
      if(ai > bj) call intswap(ai,bj)
      if(bj > ck) call intswap(bj,ck)
      if(ai > bj) call intswap(ai,bj)

      index = ai + bj*(bj-1)/2 + ck*(ck-1)*(ck+1)/6
    end function lMatIndSym

!------------------------------------------------------------------------------------------!

    pure function oldLMatInd(aI,bI,cI,iI,jI,kI) result(index)
      ! Indexing function with a 12-fold symmetry: symmetric with respect to
      ! permuting (a,i), (b,j) and (c,k) and with exchange (a,b,c)<->(i,j,k)
      ! Input: a,b,c - orbital indices of electrons
      !        i,j,k - orbital indices of holes
      ! Output: index - contiguous index I(a,b,c,i,j,k) with the aforementioned symmetry
      implicit none
      integer(int64), value, intent(in) :: aI,bI,cI ! occupied orb indices
      integer(int64), value, intent(in) :: iI,jI,kI ! unoccupied orb
      integer(int64) :: index
      integer(int64) :: a,b,c,i,j,k

      ! Somehow, I cannot directly use aI-kI, even though they are passed by value
      a = aI
      b = bI
      c = cI
      i = iI
      j = jI
      k = kI

      ! we store the permutation where a < b < c (regardless of i,j,k)
      ! or i < j < k, depending on (permuted) a < i
      ! sort such that the ordered indices start with the smallest index
      if(minval((/a,b,c/)) > minval((/i,j,k/))) then
         call intswap(a,i)
         call intswap(b,j)
         call intswap(c,k)
      endif
      ! -> create the ordered permutation on ap,bp,cp
      call sort2Els(a,b,i,j)
      call sort2Els(b,c,j,k)
      call sort2Els(a,b,i,j)

      ! indexing function: there are three ordered indices (ap,bp,cp)
      ! and three larger indices (ip,jp,kp)
      ! the last larger index kp, it is the contigous index, then follow (jp,cp)
      ! then (ip,bp) and then the smallest index ap
      index = k + nBI*(j-1) + nBI**2*(i-1) + nBI**3*(a-1) + nbI**3*(b-1)*b/2+&
           nBI**3*(c+1)*(c-1)*c/6

      contains

        ! sorts the indices a,b and i,j with respect to the
        ! ordering selected in iPermute
        pure subroutine sort2Els(r,s,p,q)
          implicit none
          integer(int64), intent(inout) :: r,s,p,q

          if(r > s) then
             call intswap(r,s)
             call intswap(p,q)
          end if
        end subroutine sort2Els

    end function oldLMatInd

!------------------------------------------------------------------------------------------!

    pure function lMatIndSymBroken(a,b,c,i,j,k) result(index)
      ! broken-symmetry index function that operates on LMat without permutational
      ! symmetry between ai, bj, ck
      ! A 6-fold symmetry remains: swapping of a<->i, b<->j and c<->k
      ! Input: a,b,c - orbital indices of electrons
      !        i,j,k - orbital indices of holes
      ! Output: index - contiguous index I(a,b,c,i,j,k) with the aforementioned symmetry
      implicit none
      integer(int64), value, intent(in) :: a,b,c ! occupied orb indices
      integer(int64), value, intent(in) :: i,j,k ! unoccupied orb
      integer(int64) :: index

      integer(int64) :: ai,bj,ck

      ai = fuseIndex(a,i)
      bj = fuseIndex(b,j)
      ck = fuseIndex(c,k)

      index = ai + strideInner*bj + strideOuter * ck

    end function lMatIndSymBroken

!------------------------------------------------------------------------------------------!

    pure function lMatIndSpin(i,j,k,a,b,c) result(index)
      ! Index functions that is symmetric with respect to swapping a<->i, b<->j and c<->k,
      ! as well as swapping (b,j)<->(c,k), but does not have the full ai,bj,ck, permutational
      ! symmetry
      ! Input: a,b,c - orbital indices of electrons
      !        i,j,k - orbital indices of holes
      ! Output: index - contiguous index I(a,b,c,i,j,k) with the aforementioned symmetry
      integer(int64), value, intent(in) :: i,j,k
      integer(int64), value, intent(in) :: a,b,c
      integer(int64) :: index
      integer(int64) :: ai, bj, ck

      ai = fuseIndex(a,i)
      bj = fuseIndex(b,j)
      ck = fuseIndex(c,k)

      index = ai + strideInner * fuseIndex(bj,ck)
    end function lMatIndSpin

end module lMat_indexing
