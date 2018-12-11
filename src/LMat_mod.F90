module LMat_mod
  use constants
  use HElem, only: HElement_t_SizeB
  use SystemData, only: tStoreSpinOrbs, nBasis
  use MemoryManager, only: LogMemAlloc, LogMemDealloc
  use util_mod, only: get_free_unit
  use shared_memory_mpi
  use ParallelHelper, only: iProcIndex_intra
  implicit none

  ! this is likely to be stored in a hashtable long term
  HElement_t(dp), pointer :: LMat(:)
  integer :: LMatTag
  integer, allocatable :: n2Ind(:), n3Ind(:)
  integer(MPIarg) :: LMatWin
  integer(int64) :: nBI
  logical :: tDebugLMat

  contains

    function get_lmat_el(a,b,c,i,j,k) result(matel)
      use SystemData, only: G1
      use UMatCache, only: gtID
      ! Gets an entry of the 3-body tensor L:
      ! L_{abc}^{ijk} - triple excitation from abc to ijk
      implicit none
      integer, value :: a,b,c
      integer, intent(in) :: i,j,k
      HElement_t(dp) :: matel
      integer(int64) :: ida, idb, idc, idi, idj, idk

      ! convert to spatial orbs if required
      ida = gtID(a)
      idb = gtID(b)
      idc = gtID(c)
      idi = gtID(i)
      idj = gtID(j)
      idk = gtID(k)

      matel = 0
      ! only add the contribution if the spins match
      if(tDebugLMat) then
         call addAllMatelConts(a,b,c,1)
         call addAllMatelConts(b,c,a,1)
         call addAllMatelConts(c,a,b,1)
         call addAllMatelConts(b,a,c,-1)
         call addAllMatelConts(a,c,b,-1)
         call addAllMatelConts(c,b,a,-1)
         matel = matel/2.0_dp
      else
         call addMatelContribution(i,j,k,idi,idj,idk,1)
         call addMatelContribution(j,k,i,idj,idk,idi,1)
         call addMatelContribution(k,i,j,idk,idi,idj,1)
         call addMatelContribution(j,i,k,idj,idi,idk,-1)
         call addMatelContribution(i,k,j,idi,idk,idj,-1)
         call addMatelContribution(k,j,i,idk,idj,idi,-1)
      endif
      contains 

        subroutine addAllMatelConts(aI,bI,cI,sgn)
          implicit none
          integer, value :: aI,bI,cI
          integer, intent(in) :: sgn
          integer :: ab, bb, cb
          
          ab = a
          bb = b
          cb = c
          
          a = aI
          b = bI
          c = cI
          
          ida = gtID(a)
          idb = gtID(b)
          idc = gtID(c)
          call addMatelContribution(i,j,k,idi,idj,idk,sgn)
          call addMatelContribution(j,k,i,idj,idk,idi,sgn)
          call addMatelContribution(k,i,j,idk,idi,idj,sgn)
          call addMatelContribution(j,i,k,idj,idi,idk,-sgn)
          call addMatelContribution(i,k,j,idi,idk,idj,-sgn)
          call addMatelContribution(k,j,i,idk,idj,idi,-sgn)    
          c = cb
          a = ab
          b = bb
        end subroutine addAllMatelConts

        subroutine addMatelContribution(p,q,r,idp,idq,idr,sgn)
          implicit none
          integer, value :: idp,idq,idr,p,q,r
          integer, intent(in) :: sgn
          integer(int64) :: ai,bj,ck
          
          if(G1(p)%ms == G1(a)%ms .and. G1(q)%ms == G1(b)%ms .and. G1(r)%ms == G1(c)%ms) then
             matel = matel + sgn * LMat(LMatInd(int(ida,int64),int(idb,int64),&
                  int(idc,int64),int(idp,int64),int(idq,int64),int(idr,int64)))
          endif
        end subroutine addMatelContribution
        
    end function get_lmat_el

!------------------------------------------------------------------------------------------!

    function LMatInd(a,b,c,i,j,k) result(index)
      implicit none
      integer(int64), intent(in) :: a,b,c ! occupied orb indices
      integer(int64), intent(in) :: i,j,k ! unoccupied orb
      integer(int64) :: index

      integer(int64) :: ai,bj,ck

      if(tDebugLMat) then
         index = k + nBI*j + nBI**2*i + nBI**3*c + nBI**4*b + nBI**5*a
      else
      ai = fuseIndex(a,i)
      bj = fuseIndex(b,j)
      ck = fuseIndex(c,k)

      index = LMatIndFused(ai,bj,ck)
   endif
    end function LMatInd

    function LMatIndFused(ai,bj,ck) result(index)
      integer(int64) :: ai,bj,ck
      integer(int64) :: index

      ! sort the indices
      if(ai > bj) call intswap(ai,bj)
      if(bj > ck) call intswap(bj,ck)
      if(ai > bj) call intswap(ai,bj)

      index = ai + bj*(bj-1)/2 + ck*(ck-1)*(ck+1)/6
      
    end function LMatIndFused

    function fuseIndex(x,y) result(xy)
      implicit none
      integer(int64), intent(in) :: x,y
      integer(int64) :: xy

      if(x < y) then
         xy = x + y*(y-1)/2
      else
         xy = y + x*(x-1)/2
      endif
    end function fuseIndex

    function oldLMatInd(a,b,c,i,j,k) result(index)
      implicit none
      integer(int64), intent(in) :: a,b,c ! occupied orb indices
      integer(int64), intent(in) :: i,j,k ! unoccupied orb
      integer(int64) :: index
      integer(int64) :: ap, bp, cp, ip, jp, kp

      ! we store the permutation where a < b < c (regardless of i,j,k)
      ! or i < j < k, depending on (permuted) a < i, b < j, c < k

      ap = min(a,i)
      ip = max(a,i)
      bp = min(j,b)
      jp = max(j,b)
      cp = min(c,k)
      kp = max(c,k)

      ! -> create that permutation
      call sort2Els(ap,bp,ip,jp)
      call sort2Els(bp,cp,jp,kp)
      call sort2Els(ap,bp,ip,jp)

      ! indexing function: there are three ordered indices (ap,bp,cp)
      ! and three larger indices (ip,jp,kp)
      ! the last larger index kp, it is the contigous index, then follow (jp,cp) 
      ! then (ip,bp) and then the smallest index ap
      index = kp + nBI*(jp-1) + nBI**2*(ip-1) + nBI**3*(ap-1) + nbI**3*(bp-1)*bp/2+&
           nBI**3*(cp+1)*(cp-1)*cp/6

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

    pure subroutine intswap(a,b)
      integer(int64), intent(inout) :: a,b
      integer(int64) :: tmp
      
      tmp = a
      a = b 
      b = tmp
    end subroutine intswap

!------------------------------------------------------------------------------------------!

    subroutine readLMat()
      implicit none
      character(*), parameter :: LMatFileName = "TCDUMP"
      integer :: iunit, ierr
      integer(int64) :: LMatSize
      integer(int64) :: a,b,c,i,j,k
      HElement_t(dp) :: matel
      character(*), parameter :: t_r = "readLMat"
      integer :: counter
      real(dp) :: fac

      if(tStoreSpinOrbs) then
         nBI = nBasis
      else
         nBI = nBasis / 2
      endif

      call initializeNInd()

      ! The size is given by the largest index (LMatInd is monotonous in all arguments)
      LMatSize = LMatInd(nBI,nBI,nBI,nBI,nBI,nBI)

      ! allocate LMat (shared memory)
      call shared_allocate_mpi(LMatWin, LMat, (/LMatSize/))
      LMat = 0.0_dp
      call LogMemAlloc("LMat", int(LMatSize), HElement_t_SizeB, t_r, LMatTag)

      if(iProcIndex_intra .eq. 0) then
         iunit = get_free_unit()
         open(iunit,file = LMatFileName,status = 'old')
         counter = 0
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
               if(LMatInd(a,b,c,i,j,k) > LMatSize) then
                  counter = LMatInd(a,b,c,i,j,k)
                  write(iout,*) "Warning, exceeding size" 
               endif
               LMat(LMatInd(a,b,c,i,j,k)) = 3.0_dp * matel
               if(abs(matel)> 0.0_dp) counter = counter + 1
            endif

         end do

         counter = counter / 12

         write(iout, *), "Sparsity of LMat", real(counter)/real(LMatSize)
         write(iout, *), "Nonzero elements in LMat", counter
         write(iout, *), "Allocated size of LMat", LMatSize
      endif
      
    end subroutine readLMat

!------------------------------------------------------------------------------------------!

    subroutine initializeNInd()
      implicit none
      integer :: i,j      
      ! prepare an array containing the offset of certain blocks in the
      ! LMat array
      
      allocate(n2Ind(nBI))
      n2Ind(1) = 0
      do i = 2, nBI
         n2Ind(i) = n2Ind(i-1) + (nBI + 2 - i)**2
      end do

      allocate(n3Ind(nBI))
      n3Ind(1) = 1
      do i = 1, nBI-1
         n3Ind(nBI-i) = n3Ind(nBI-i+1) + (nBI+1-i)**2
      end do
    end subroutine initializeNInd

!------------------------------------------------------------------------------------------!

    subroutine freeLMat()
      implicit none
      character(*), parameter :: t_r = "freeLMat"
      
      if(associated(LMat)) then
         call shared_deallocate_mpi(LMatWin,LMat)
         call LogMemDealloc(t_r, LMatTag)
         LMAt => null()
      end if
      deallocate(n3Ind)
      deallocate(n2Ind)
      
    end subroutine freeLMat

end module LMat_mod
