module LMat_mod
  use constants
  use FciMCData, only: ll_node
  use HElem, only: HElement_t_SizeB
  use SystemData, only: tStoreSpinOrbs, nBasis, tHDF5LMat, t12FoldSym
  use MemoryManager, only: LogMemAlloc, LogMemDealloc
  use util_mod, only: get_free_unit
  use shared_memory_mpi
  use sort_mod
  use hash, only: add_hash_table_entry, clear_hash_table
  use ParallelHelper, only: iProcIndex_intra
  use tc_three_body_data, only: tDampKMat, tDampLMat, tSymBrokenLMat, tSpinCorrelator, LMatEps, &
       lMat_t, tSparseLMat
  use procedure_pointers, only: LMatInd, get_lmat_el, get_lmat_el_symInternal
  use LoggingData, only: tHistLMat
#ifdef __USE_HDF5
  use hdf5
#endif
  implicit none


  type(lMat_t) :: LMat, LMatAB
  integer(int64) :: nBI
  logical :: tDebugLMat

  ! for the symmetry broken index function
  integer :: strideInner, strideOuter

  ! access function for lMat
  abstract interface
     function lmat_access_t(lMatObj, index) result(matel)
       use constants
       use tc_three_body_data, only: lMat_t
       implicit none
       type(lMat_t), intent(in) :: lMatObj
       integer(int64), intent(in) :: index
       HElement_t(dp) :: matel
     end function lmat_access_t
  end interface

  procedure(lmat_access_t), pointer :: lMatAccess


  contains

    function get_lmat_el_base(a,b,c,i,j,k) result(matel)
      use SystemData, only: G1
      use UMatCache, only: gtID
      ! Gets an entry of the 3-body tensor L:
      ! L_{abc}^{ijk} - triple excitation from abc to ijk
      implicit none
      integer, value :: a,b,c
      integer, intent(in) :: i,j,k
      HElement_t(dp) :: matel
      integer :: ida, idb, idc, idi, idj, idk
      logical :: tSameSpin

      ! initialize spin-correlator check: if all spins are the same, use LMat, else LMatAB
      ! without spin-dependent correlator, always use LMat
      tSameSpin = .not. tSpinCorrelator .or. (G1(i)%MS==G1(j)%ms .and. G1(i)%MS==G1(k)%MS)
      
      ! convert to spatial orbs if required
      ida = gtID(a)
      idb = gtID(b)
      idc = gtID(c)
      idi = gtID(i)
      idj = gtID(j)
      idk = gtID(k)

      matel = 0
      ! only add the contribution if the spins match
     
      call addMatelContribution(i,j,k,idi,idj,idk,1)         
      call addMatelContribution(j,k,i,idj,idk,idi,1)
      call addMatelContribution(k,i,j,idk,idi,idj,1)
      call addMatelContribution(j,i,k,idj,idi,idk,-1)
      call addMatelContribution(i,k,j,idi,idk,idj,-1)
      call addMatelContribution(k,j,i,idk,idj,idi,-1)

      ! spin-projector for the tc terms - apply heuristically for three-body terms
      if(tDampKMat) then
         ! same-spin contributions are divided by 4 (this is exact)
         if(G1(a)%MS .eq. G1(b)%MS .and. G1(b)%MS .eq. G1(c)%MS) then
            matel = matel/4.0_dp
         else
            ! opposite-spin contributions are divided by 2 (this is a guess, the
            ! exact form has an admixture of exchange terms)
            matel = matel/2.0_dp
         endif
      endif
      contains 

        subroutine addMatelContribution(p,q,r,idp,idq,idr,sgn)
          implicit none
          integer, value :: idp,idq,idr,p,q,r
          integer, intent(in) :: sgn
          
          if(G1(p)%ms == G1(a)%ms .and. G1(q)%ms == G1(b)%ms .and. G1(r)%ms == G1(c)%ms) then
             if(tSameSpin) then
                matel = matel + sgn * real(lMatAccess(lMat,LMatInd(int(ida,int64),int(idb,int64),&
                     int(idc,int64),int(idp,int64),int(idq,int64),int(idr,int64))),dp)
             else
                matel = matel + sgn * real(lMatAccess(lMatAB,LMatInd(int(ida,int64),int(idb,int64),&
                     int(idc,int64),int(idp,int64),int(idq,int64),int(idr,int64))),dp)
             endif
          endif
        end subroutine addMatelContribution
        
      end function get_lmat_el_base

!------------------------------------------------------------------------------------------!

      function get_lmat_el_symmetrized(a,b,c,i,j,k) result(matel)
        ! post-symmetrized access when storing non-symmetrized lMat
        implicit none
        integer, value :: a,b,c
        integer, intent(in) :: i,j,k
        HElement_t(dp) :: matel

        matel = 0.0_dp
        
        ! get_lmat_el_base is not symmetric in this case => symmetrize w.r. 
        ! to exchange of electrons
        matel = matel + get_lmat_el_symInternal(a,b,c,i,j,k)
        matel = matel + get_lmat_el_symInternal(b,c,a,j,k,i)
        matel = matel + get_lmat_el_symInternal(c,a,b,k,i,j)
        ! + spin correction (get_lmat_el_base(ap,bp,cp,ip,jp,kp) with ap etc being the 
        ! spin-swapped indices)
        matel = matel / 3.0_dp
        
      end function get_lmat_el_symmetrized

!------------------------------------------------------------------------------------------!
      
      function get_lmat_el_spinProj(a,b,c,i,j,k) result(matel)
        ! get the spin-projected matel
        implicit none
        integer, value :: a,b,c
        integer, intent(in) :: i,j,k
        HElement_t(dp) :: matel

        ! auxiliary, spin-swapped indices
        integer :: ap, bp, cp, ip, jp, kp
        ! prefactors for the different parts
        real(dp), parameter :: directFac = 0.5625_dp
        real(dp), parameter :: swapFac = -0.1875_dp
        real(dp), parameter :: permFac = 0.0625_dp
        
        matel = 0.0_dp
        matel = matel + directFac * get_lmat_el_base(a,b,c,i,j,k)
        call resetAux()
        call swapSpins(ap,bp,ip,jp)
        matel = matel + swapFac * get_lmat_el_base(ap,bp,cp,ip,jp,kp)
        call resetAux()
        call swapSpins(ap,cp,ip,kp)
        matel = matel + swapFac * get_lmat_el_base(ap,bp,cp,ip,jp,kp)
        call resetAux()
        call swapSpins(ap,cp,ip,kp)
        call swapSpins(ap,bp,ip,jp)
        matel = matel + swapFac * get_lmat_el_base(ap,bp,cp,ip,jp,kp)
      
      contains

        subroutine resetAux()
          implicit none
          ! reset the auxiliary variables
          ap = a
          bp = b
          cp = c
          ip = i
          jp = j
          kp = k
        end subroutine resetAux

        subroutine swapSpins(src1, src2, tgt1, tgt2)
          implicit none
          integer, intent(inout) :: src1, src2, tgt1, tgt2
          integer :: mst1, mst2, mss1, mss2

          ! get the spin values
          mss1 = mod(src1,2)
          mss2 = mod(src2,2)
          mst1 = mod(tgt1,2)
          mst2 = mod(tgt2,2)

          ! swap two spins
          src1 = src1 + mss1 - mss2
          src2 = src2 + mss2 - mss1
          tgt1 = tgt1 + mst1 - mst2
          tgt2 = tgt2 + mst2 - mst1
        end subroutine swapSpins
      end function get_lmat_el_spinProj

!------------------------------------------------------------------------------------------!

    function LMatIndSym(a,b,c,i,j,k) result(index)
      implicit none
      integer(int64), intent(in) :: a,b,c ! occupied orb indices
      integer(int64), intent(in) :: i,j,k ! unoccupied orb
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
    end function LMatIndSym

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
      integer(int64), value :: a,b,c ! occupied orb indices
      integer(int64), value :: i,j,k ! unoccupied orb
      integer(int64) :: index

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

    function lMatIndSymBroken(a,b,c,i,j,k) result(index)
      ! broken-symmetry index function that operates on LMat without permutational
      ! symmetry between ai, bj, ck
      implicit none
      integer(int64), intent(in) :: a,b,c ! occupied orb indices
      integer(int64), intent(in) :: i,j,k ! unoccupied orb
      integer(int64) :: index

      integer(int64) :: ai,bj,ck

      ai = fuseIndex(a,i)
      bj = fuseIndex(b,j)
      ck = fuseIndex(c,k)

      index = ai + strideInner*bj + strideOuter * ck
      
    end function lMatIndSymBroken

!------------------------------------------------------------------------------------------!

    subroutine initializeLMatInd()
      implicit none
     
      ! set the LMatInd function pointer
      if(t12FoldSym) then
         LMatInd => oldLMatInd
      else if(tSymBrokenLMat) then
         ! also need to set the size of the blocks
         strideInner = fuseIndex(nBI,nBI)
         strideOuter = strideInner**2
         LMatInd => lMatIndSymBroken
      else
         LMatInd => LMatIndSym
      endif

      ! set the get_lmat_el function pointer
      if(tSymBrokenLMat) then
         get_lmat_el => get_lmat_el_symmetrized
      else
         get_lmat_el => get_lmat_el_base
      endif

      ! the internal pointer of get_lmat_el_symmetrized
      if(tDampLMat) then
         get_lmat_el_symInternal => get_lmat_el_spinProj
      else
         get_lmat_el_symInternal => get_lmat_el_base
      endif

      if(tSparseLMat) then
         lMatAccess => lMatHashedAccess
      else
         lMatAccess => lMatDirectAccess
      endif
    end subroutine initializeLMatInd

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
      use SystemData, only: nel
      implicit none

      ! we need at least three electrons to make use of the six-index integrals
      ! => for less electrons, this can be skipped
      if(nel<=2) return

      if(tStoreSpinOrbs) then
         nBI = nBasis
      else
         nBI = nBasis / 2
      endif

      call initializeLMatInd()

      ! now, read lmat from file
      call readLMatArray(LMat,"TCDUMP","tcdump.h5")
      ! for spin-dependent LMat, also read the opp. spin matrix
      if(tSpinCorrelator) call readLMatArray(LMatAB,"TCDUMPAB","tcdumpab.h5")
    end subroutine readLMat

!------------------------------------------------------------------------------------------!

    subroutine readLMatArray(LMatLoc, LMatFileName, h5Filename)
      use ParallelHelper, only: mpi_comm_intra
      implicit none
      type(lMat_t), intent(inout) :: LMatLoc
      character(*), intent(in) :: LMatFileName
      character(*), intent(in) :: h5Filename
      integer :: iunit, ierr
      integer(int64) :: LMatSize
      integer(int64) :: a,b,c,i,j,k
      HElement_t(dp) :: matel
      character(*), parameter :: t_r = "readLMat"
      integer :: counter
      integer(int64) :: iChunk, chunkSize
      integer :: iNodeSize, iProcInNode
      integer :: dummy(6)
      integer(MPIArg) :: err

      if(.not.tSparseLMat) then
         ! for sparse storage, we first need to get the number of integrals
         ! this works differently for hdf5 and non-hdf5 dumpfiles, so it is
         ! done directly in the respective code

         ! The size is given by the largest index (LMatInd is monotonous in all arguments)
         LMatSize = LMatInd(nBI,nBI,nBI,nBI,nBI,nBI)

         call allocLMat(LMatLoc,LMatSize)
      endif

      if(tHDF5LMat) then
#ifdef __USE_HDF
         call readHDF5LMat(LMatLoc,h5filename)
#else
         call stop_all(t_r, "HDF5 integral files disabled at compile time")
#endif
      else

         if(tSparseLMat) call stop_all(t_r,"Sparse storage requires hdf5-integrals")
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
                  if(abs(3.0_dp*matel) > LMatEps) &
                       LMatLoc%LMatPtr(LMatInd(a,b,c,i,j,k)) = 3.0_dp * matel
                  if(abs(matel)> 0.0_dp) counter = counter + 1
               endif

            end do

            counter = counter / 12

            write(iout, *), "Sparsity of LMat", real(counter)/real(LMatSize)
            write(iout, *), "Nonzero elements in LMat", counter
            write(iout, *), "Allocated size of LMat", LMatSize
         endif
         LMatLoc%nInts = counter
      end if

      if(tHistLMat) call histogramLMat(LMatLoc)
      
    end subroutine readLMatArray

!------------------------------------------------------------------------------------------!

    subroutine allocLMat(LMatLoc,LMatSize)
      use HElem, only: Helement_t_sizeB
      implicit none
      type(lMat_t) :: LMatLoc
      integer(int64), intent(in) :: LMatSize
      integer :: k
      character(*), parameter :: t_r = "allocLMat"
      
      write(iout,*) "Allocating LMat, memory required:", LMatSize*HElement_t_sizeB/(2.0_dp**30), "GB"
      LMatLoc%nInts = LMatSize
      ! allocate LMat (shared memory)
      call shared_allocate_mpi(LMatLoc%shm_win, LMatLoc%LMatPtr, (/LMatSize/))
      call LogMemAlloc("LMat", int(LMatSize), HElement_t_SizeB, t_r, LMatLoc%tag)
      if(iProcIndex_intra.eq.0) then
         do k = 1, LMatSize
            LMatLoc%LMatPtr(k) = 0.0_dp
         end do
      endif
      write(iout,*) "Successfully allocated LMat"
    end subroutine allocLMat

!------------------------------------------------------------------------------------------!

    subroutine freeLMat()
      implicit none
      character(*), parameter :: t_r = "freeLMat"

      call deallocLMatArray(LMat)
      call deallocLMatArray(LMatAB)

      contains 

        subroutine deallocLMatArray(LMatLoc)
          implicit none
          type(LMat_t) :: LMatLoc

          if(associated(LMatLoc%LMatPtr)) then
             call shared_deallocate_mpi(LMatLoc%shm_win,LMatLoc%LMatPtr)
             call LogMemDealloc(t_r, LMatLoc%tag)
             LMatLoc%LMatPtr => null()
          end if

          if(associated(LMatLoc%indexPtr)) then
             call shared_deallocate_mpi(LMatLoc%index_win,LMatLoc%indexPtr)
             call LogMemDealloc(t_r, LMatLoc%indexTag)
             LMatLoc%indexPtr => null()
          endif

          if(associated(LMatLoc%htable)) then
             call clear_hash_table(LMatLoc%htable)
             deallocate(LMatLoc%htable)
          endif

        end subroutine deallocLMatArray

      end subroutine freeLMat

!------------------------------------------------------------------------------------------!
! For huge LMats, a sparse storage scheme is required. Since we still need to access
! the matrix elements, we use a hash table and re-use the hash.F90 module
!------------------------------------------------------------------------------------------!

      subroutine initLMatHash(lMatCtr)
        implicit none
        type(LMat_t), intent(in) :: lMatCtr
        integer :: i
        integer(int64) :: hashVal

        allocate(lMatCtr%htable(lMatCtr%htSize))

        ! for each entry, store the index at the position in the hashtable given
        ! by its hash value
        do i = 1, lMatCtr%nInts
           hashVal = LMatHash(lMatCtr%indexPtr(i), lMatCtr%htSize)
           call add_hash_table_entry(lMatCtr%htable, i, hashVal)
        end do
      end subroutine initLMatHash


!------------------------------------------------------------------------------------------!

      function LMatHash(index, htSize) result(hashVal)
        implicit none
        integer(int64), intent(in) :: index
        integer, intent(in) :: htSize
        integer(int64) :: hashVal

! TODO: Implement an actual hash function
        hashVal = mod(index-1,htSize)+1
      end function LMatHash

!------------------------------------------------------------------------------------------!

      function lMatHashedAccess(lMatCtr,index) result(matel)
        implicit none
        type(LMat_t), intent(in) :: lMatCtr
        integer(int64), intent(in) :: index
        HElement_t(dp) :: matel

        integer(int64) :: hashVal
        type(ll_node), pointer :: tmp_node
        logical :: found
        
        hashVal = LMatHash(index, lMatCtr%htSize)
        ! this is the kernel of hash_table_lookup from the hash module
        found = .false.
        tmp_node => lMatCtr%htable(hashVal)
        if(tmp_node%ind /= 0) then
           do while(associated(tmp_node))
              if(index==lMatCtr%indexPtr(tmp_node%ind)) then
                 found = .true.
                 exit
              endif
              tmp_node => tmp_node%next
           end do
        end if
        
        ! now, tmp_node%ind is the index of the element we wanted to have
        if(found) then 
           matel = lMatCtr%LMatPtr(tmp_node%ind)
        else
           ! if not found, the matrix element is 0
           matel = 0.0_dp
        endif

        nullify(tmp_node)
      end function lMatHashedAccess

!------------------------------------------------------------------------------------------!

      ! direct (dense) access function
      function lMatDirectAccess(lMatCtr, index) result(matel)
        implicit none
        type(lMat_t), intent(in) :: lMatCtr
        integer(int64), intent(in) :: index
        HElement_t(dp) :: matel
        
        matel = lMatCtr%lMatPtr(index)
      end function lMatDirectAccess

!------------------------------------------------------------------------------------------!
! HDF5 Functionality
!------------------------------------------------------------------------------------------!


#ifdef __USE_HDF
    subroutine readHDF5LMat(LMatLoc, filename)
      use hdf5
      use hdf5_util
      use ParallelHelper, only: mpi_comm_inter, mpiInfoNull, iProcIndex, &
           mpi_comm_intra
      use Parallel_neci, only: nProcessors
      implicit none
      type(lMat_t) :: LMatLoc
      character(*), intent(in) :: filename
      integer :: proc, i
      integer(hid_t) :: err, file_id, plist_id, grp_id, ds_vals, ds_inds
      integer(hsize_t) :: nInts, rest, blocksize, blockstart, blockend, this_blocksize, countsEnd
      integer(hsize_t), allocatable :: counts(:), offsets(:)
      integer(hsize_t), allocatable :: indices(:,:), entries(:,:)
      character(*), parameter :: nm_grp = "tcdump", nm_nInts = "nInts", nm_vals = "values", nm_indices = "indices"
      integer(MPIArg) :: procs_per_node, ierr
      integer :: sparseBlockStart
      real(dp) :: rVal
      logical :: running, any_running
      integer :: this_nEntries, counter, sparseBlock, allCounter
      integer, allocatable :: all_nEntries(:)
      character(*), parameter :: t_r = "readHDF5LMat"
      rVal = 0.0_dp

      call h5open_f(err)
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, err)
      call h5pset_fapl_mpio_f(plist_id, mpi_comm_inter, mpiInfoNull, err)

      ! open the file
      call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, err, access_prp=plist_id)

      call h5gopen_f(file_id, nm_grp, grp_id, err)

      ! get the number of integrals
      call read_int64_attribute(grp_id, nm_nInts, nInts, required=.true.)
      write(iout,*) "Reading", nInts, "integrals"

      ! if the LMat is stored sparse, we can allocate now
      if(tSparseLMat) then
         call allocLMat(LMatLoc,nInts)
         ! also, allocate the index array
         call shared_allocate_mpi(LMatLoc%index_win, LMatLoc%indexPtr, (/int(nInts,int64)/))
         call LogMemAlloc("LMat Indices", int(nInts), sizeof_int64, t_r, LMatLoc%indexTag)
         ! come up with some reasonable size
         LMatLoc%htSize = nInts
         ! start filling LMat at the first entry
         sparseBlockStart = 0
      endif

      ! how many entries does each proc get?
      call MPI_Comm_Size(mpi_comm_intra, procs_per_node, ierr)
      allocate(counts(0:procs_per_node-1))
      allocate(offsets(0:procs_per_node-1))
      ! the number of entries read in per node - required for sparse storage
      if(tSparseLMat) allocate(all_nEntries(0:procs_per_node-1))
      counts = nInts / int(procs_per_node, hsize_t)
      rest = mod(nInts, procs_per_node)
      if(rest>0) counts(0:rest-1) = counts(0:rest-1) + 1

      offsets(0) = 0
      do proc = 1, procs_per_node - 1
         offsets(proc) = offsets(proc-1) + counts(proc-1)
      end do
      call h5dopen_f(grp_id, nm_vals, ds_vals, err)
      call h5dopen_f(grp_id, nm_indices, ds_inds, err)

      ! reserve max. 100MB buffer size for dumpfile I/O
      blocksize = 100000000/(7*sizeof(LMatLoc%LMatPtr(1)))
      blockstart = offsets(iProcIndex_intra)

      ! the last element to read on each proc
      if(iProcIndex_intra.eq.nProcessors - 1) then
         countsEnd = nInts - 1
      else
         countsEnd = offsets(iProcIndex_intra + 1) - 1
      end if
      blockend = min(blockstart + blocksize - 1, countsEnd)
      any_running = .true.
      running = .true.
      do while(any_running)
         if(running) then
            ! the number of elements to read in this block
            this_blocksize = blockend - blockstart + 1
         else
            this_blocksize = 0
         end if
         ! number of integrals read in this block
         counter = 0

         allocate(indices(6,this_blocksize), source = 0_int64)
         allocate(entries(1,this_blocksize), source = 0_int64)

         ! read in the data
         call read_2d_multi_chunk(&
              ds_vals, entries, H5T_NATIVE_REAL_8, &
              [1_hsize_t, this_blocksize],&
              [0_hsize_t, blockstart],&
              [0_hsize_t, 0_hsize_t])

         call read_2d_multi_chunk(&
              ds_inds, indices, H5T_NATIVE_INTEGER_8, &
              [6_hsize_t, this_blocksize], &
              [0_hsize_t, blockstart], &
              [0_hsize_t, 0_hsize_t])

         ! for a sparse format, we need to know NOW, how many nonzero entries each
         ! processor contributes
         if(tSparseLMat) then
            this_nEntries = count_entries(entries)
            ! communicate the number of nonzeros
            call MPI_AllGather(this_nEntries,1,MPI_INT,all_nEntries,1,MPI_INT,mpi_comm_intra,ierr)
            ! the window we write to with this proc starts with this offset
            sparseBlock = sparseBlockStart
            do i = 0, iProcIndex_intra-1
               sparseBlock = sparseBlock + all_nEntries(i)
            end do
         endif

         ! assign LMat entries
         do i = 1, this_blocksize
            ! truncate down to lMatEps
            rVal = 3.0_dp * transfer(entries(1,i),rVal)
            if(abs(rVal)>lMatEps) then
               if(tSparseLMat) then
                  ! increase the number of read-in integrals
                  counter = counter + 1
                  ! write to the local window within the shared memory
                  LMatLoc%LMatPtr(sparseBlock+counter) = rVal
                  LMatLoc%indexPtr(sparseBlock+counter) = LMatInd(int(indices(1,i),int64),int(indices(2,i),int64),&
                       int(indices(3,i),int64),&
                       int(indices(4,i),int64),int(indices(5,i),int64),int(indices(6,i),int64))
               else
                  LMatLoc%LMatPtr(LMatInd(int(indices(1,i),int64),int(indices(2,i),int64),&
                       int(indices(3,i),int64),&
                       int(indices(4,i),int64),int(indices(5,i),int64),int(indices(6,i),int64))) &
                       = rVal
               endif
            endif
         end do

         ! set the size/offset for the next block
         if(running) then
            blockstart = blockend + 1
            blockend = min(blockstart + blocksize - 1, countsEnd)
            if(blockstart > countsEnd) running = .false.
         end if

         deallocate(entries)
         deallocate(indices)

         call MPIAllLORLogical(running, any_running)

         ! communicate how many nonzero entries have been read and set the starting point 
         ! for the next write
         if(tSparseLMat) then
            call MPISumAll(counter,allCounter)
            sparseBlockStart = sparseBlockStart + allCounter
         endif
      end do
      deallocate(offsets)
      deallocate(counts)
      call h5dclose_f(ds_vals, err)
      call h5dclose_f(ds_inds, err)
      ! close the file, finalize hdf5
      call h5gclose_f(grp_id, err)
      call h5pclose_f(plist_id, err)
      call h5fclose_f(file_id, err)
      call h5close_f(err)

      if(tSparseLMat) then
         LMatLoc%nInts = sparseBlockStart
         call initLMatHash(LMatLoc)
      endif

      contains 

        function count_entries(entries) result(nEntries)
          integer(hsize_t), intent(in) :: entries(:,:)
          integer :: blocksize, j
          integer :: nEntries

          blocksize = size(entries,2)
          nEntries = 0
          ! loop over the matrix elements in entries and count the number of
          ! non-truncated matrix elements
          do j = 1, blocksize
             if(abs(3.0*transfer(entries(1,j),rVal)) > lMatEps) nEntries = nEntries + 1
          end do
             
        end function count_entries

    end subroutine readHDF5LMat
#endif

!------------------------------------------------------------------------------------------!

    subroutine histogramLMat(lMatObj)
      implicit none
      type(lMat_t), intent(in) :: lMatObj
      integer :: i,thresh
      integer, parameter :: minExp = 10
      integer :: histogram(0:minExp)
      real :: ratios(0:minExp)
      
      histogram = 0
      do i = 1, lMatObj%nInts
         do thresh = minExp,1,-1
            ! in each step, count all matrix elements that are below the threshold and
            ! have not been counted yet
            if(abs(lMatObj%LMatPtr(i))<0.1**(thresh)) then
               histogram(thresh) = histogram(thresh) + 1
               ! do not count this one again
               exit
            endif
            ! the last check has a different form: everything that is bigger than 0.1 counts here
            if(abs(lMatObj%LMatPtr(i)) > 0.1) histogram(0) = histogram(0) + 1
         end do
      end do

      ratios(:) = real(histogram(:))/real(lMatObj%nInts)
      ! print the ratios
      write(iout,*) "Matrix elements below", 0.1**(minExp), ":", ratios(minExp)
      do i = minExp-1,1,-1
         write(iout,*) "Matrix elements from", 0.1**(i+1),"to",0.1**(i),":",ratios(i)
      end do
      write(iout,*) "Matrix elements above", 0.1,":",ratios(0)
      write(iout,*), "Total number of logged matrix elements", lMatObj%nInts
    end subroutine histogramLMat

end module LMat_mod
