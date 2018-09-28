module LMat_mod
  use constants
  use HElem, only: HElement_t_SizeB
  use SystemData, only: tStoreSpinOrbs, nBasis
  use MemoryManager, only: LogMemAlloc, LogMemDealloc
  use util_mod, only: get_free_unit
  implicit none

  ! this is likely to be stored in a hashtable long term
  HElement_t(dp), allocatable :: LMat(:)
  integer :: nBI, LMatTag

  contains

    function get_lmat_el(a,b,c,i,j,k) result(matel)
      ! Gets an entry of the 3-body tensor L:
      ! L_{abc}^{ijk} - triple excitation from abc to ijk
      implicit none
      integer, intent(in) :: a,b,c
      integer, intent(in) :: i,j,k
      HElement_t(dp) :: matel

      matel = LMat(LMatInd(a,b,c,i,j,k)) &
           + LMat(LMatInd(a,b,c,j,k,i)) &
           + LMat(LMatInd(a,b,c,k,i,j)) &
           - LMat(LMatInd(a,b,c,j,i,k)) &
           - LMat(LMatInd(a,b,c,i,k,j)) &
           - LMat(LMatInd(a,b,c,k,j,i)) 
    end function get_lmat_el

!------------------------------------------------------------------------------------------!

    function LMatInd(a,b,c,i,j,k) result(index)
      implicit none
      integer, intent(in) :: a,b,c ! occupied orb indices
      integer, intent(in) :: i,j,k ! unoccupied orb
      integer :: index, ap, bp, cp, ip, jp, kp
      logical :: iPermute

      ! we store the permutation where a < b < c (regardless of i,j,k)
      ! or i < j < k, depending on (permuted) a < i
      iPermute = minVal([a,b,c]) > minVal([i,j,k])
      if(iPermute) then
         ! permute (a,b,c) <-> (i,j,k) if required
         ap = i
         bp = j
         cp = k
         ip = a
         jp = b
         kp = c
      else
         ! or not
         ap = a
         bp = b
         cp = c
         ip = i
         jp = j
         kp = k
      endif
      
      ! -> create that permutation
      call sort2Els(ap,bp,ip,jp)
      call sort2Els(bp,cp,jp,kp)
      call sort2Els(ap,bp,ip,jp)

      ! indexing function: there are three ordered indices (ap,bp,cp)
      ! and three unrestricted indices (ip,jp,kp)
      ! the last unrestricted index is the contiguous one, then follow the other unrestricted
      ! ones
      ! then, the restricted ones follow, limited by the previous index (the last ordered
      ! index has no restriction
      index = kp + nBI * (jp-1) + nBI**2 * (ip-1) + nBI**3 * (ap-1) + &
           nBI**3 * (bp-1)*bp/2 + nBI**3 * (cp+1)*(cp-1)*cp/6

      contains

        ! sorts the indices a,b and i,j with respect to the 
        ! ordering selected in iPermute
        pure subroutine sort2Els(r,s,p,q)
          implicit none
          integer, intent(inout) :: r,s,p,q

          if(r > s) then
             call intswap(r,s)
             call intswap(p,q)
          end if
        end subroutine sort2Els
      
    end function LMatInd

!------------------------------------------------------------------------------------------!

    pure subroutine intswap(a,b)
      integer, intent(inout) :: a,b
      integer :: tmp
      
      tmp = a
      a = b 
      b = tmp
    end subroutine intswap

!------------------------------------------------------------------------------------------!

    subroutine readLMat()
      implicit none
      character(*), parameter :: LMatFileName = "TCDUMP"
      integer :: iunit, ierr
      integer :: LMatSize
      integer :: a,b,c,i,j,k
      HElement_t(dp) :: matel
      character(*), parameter :: t_r = "readLMat"

      if(tStoreSpinOrbs) then
         nBI = nBasis
      else
         nBI = nBasis / 2
      endif

      ! The size is given by the largest index (LMatInd is monotonous in all arguments)
      LMatSize = LMatInd(nBI,nBI,nBI,nBI,nBI,nBI)

      ! allocate LMat
      allocate(LMat(LMatSize), stat = ierr)
      LMat = 0.0_dp
      call LogMemAlloc("LMat", LMatSize, HElement_t_SizeB, t_r, LMatTag, ierr)

      iunit = get_free_unit()
      open(iunit,file = LMatFileName,status = 'old')

      do
         read(iunit,*,iostat = ierr) matel, a,b,c,i,j,k
         ! end of file reached?
         if(ierr < 0) then
            exit
         else if(ierr > 0) then
            ! error while reading?
            call stop_all(t_r,"Error reading TCDUMP file")
         else
            ! else assign the matrix element
            LMat(LMatInd(a,b,c,i,j,k)) = matel
         endif
         
      end do
    end subroutine readLMat

!------------------------------------------------------------------------------------------!

    subroutine freeLMat()
      implicit none
      character(*), parameter :: t_r = "freeLMat"
      
      if(allocated(LMat)) then
         deallocate(LMat)
         call LogMemDealloc(t_r, LMatTag)
      end if
      
    end subroutine freeLMat

end module LMat_mod
