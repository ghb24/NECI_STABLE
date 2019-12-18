module LMat_aux
  ! This module contains common routines used in handling multi-index arrays, like
  ! functions used to create indexing functions. It is used in the LMat modules
  use constants
  use SystemData, only: G1
  
  implicit none

  contains

!------------------------------------------------------------------------------------------!
    
    pure subroutine dampLMatel(a,b,c,matel)
      ! spin-projector for the tc terms - apply heuristically for three-body terms
      ! Input: a,b,c - Target orbital indices
      ! Input/Output: matel - on input, the bare matrix element, on return, the
      !                       matrix element including spin-correction
      integer, intent(in) :: a,b,c
      HElement_t(dp), intent(inout) :: matel


      ! same-spin contributions are divided by 4 (this is exact)
      if(G1(a)%MS .eq. G1(b)%MS .and. G1(b)%MS .eq. G1(c)%MS) then
         matel = matel/4.0_dp
      else
         ! opposite-spin contributions are divided by 2 (this is a guess, the
         ! exact form has an admixture of exchange terms)
         matel = matel/2.0_dp
      endif
    end subroutine dampLMatel

!------------------------------------------------------------------------------------------!

    function diffSpinPos(i,j,k,a,b,c) result(pos)
      ! given three excitations (i,a), (j,b), (k,c), find the position of the one with
      ! different spin in the ordering a'<b'<c' where a'=min(a,i) etc.
      ! Input: i,j,k - occupied orbital indices
      !        a,b,c - unoccupied orbital indices
      ! Outpu: pos - 0 if all electrons have the same spin, the position of the one
      !              with different spin else
      implicit none
      integer, intent(in) :: i,j,k,a,b,c
      integer :: pos
      integer :: ap, bp, cp
      integer :: tmp(3)

      ! the minimum of each pair to be sorted
      ap = min(a,i)
      bp = min(b,j)
      cp = min(c,k)
      tmp = (/ap,bp,cp/)

      ! at this point, both indices of a pair have the same spin, so we just use the one
      ! of the primed indices
      if(G1(minval(tmp))%MS.ne.G1(maxval(tmp))%MS) then
         ! either the first or the last is different
         ! do a check if the first entry has different spin: subtract its spin from the
         ! total spin
         if(abs(G1(ap)%MS+G1(bp)%MS+G1(cp)%MS-G1(minval(tmp))%MS)==2) then
            ! if this is 2, minval(tmp) has different spin than the other two
            pos = 1
         else
            ! else, maxval(tmp) has the different spin
            pos = 3
         endif
      else
         pos = 2
      endif

    end function diffSpinPos

end module LMat_aux
