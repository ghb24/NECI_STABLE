module lMat_indexing
  use constants, only: int64
  use util_mod, only: intSwap, fuseIndex
  use tc_three_body_data, only: nBI
  implicit none

  ! for the symmetry broken index function
  integer :: strideInner, strideOuter

  contains

!------------------------------------------------------------------------------------------!    
!  Index functions for the six-index addressing      
!------------------------------------------------------------------------------------------!

    pure function lMatIndSym(a,b,c,i,j,k) result(index)
      implicit none
      integer(int64), value :: a,b,c ! occupied orb indices
      integer(int64), value :: i,j,k ! unoccupied orb
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
      implicit none
      integer(int64), value :: aI,bI,cI ! occupied orb indices
      integer(int64), value :: iI,jI,kI ! unoccupied orb
      integer(int64) :: index
      integer(int64) :: a,b,c,i,j,k

      ! guarantee pass-by-value without changing the signature to value
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
      implicit none
      integer(int64), value :: a,b,c ! occupied orb indices
      integer(int64), value :: i,j,k ! unoccupied orb
      integer(int64) :: index

      integer(int64) :: ai,bj,ck

      ai = fuseIndex(a,i)
      bj = fuseIndex(b,j)
      ck = fuseIndex(c,k)

      index = ai + strideInner*bj + strideOuter * ck
      
    end function lMatIndSymBroken
end module lMat_indexing
