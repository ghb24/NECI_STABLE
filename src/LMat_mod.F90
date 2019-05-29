module LMat_mod
  use constants
  use HElem, only: HElement_t_SizeB
  use SystemData, only: tStoreSpinOrbs, nBasis, tHDF5LMat, t12FoldSym
  use MemoryManager, only: LogMemAlloc, LogMemDealloc
  use util_mod, only: get_free_unit
  use shared_memory_mpi
  use sort_mod
  use ParallelHelper, only: iProcIndex_intra
  use tc_three_body_data, only: tDampKMat, tDampLMat, tSymBrokenLMat, tSpinCorrelator
  use procedure_pointers, only: LMatInd, get_lmat_el, get_lmat_el_symInternal
#ifdef __USE_HDF5
  use hdf5
#endif
  implicit none

  HElement_t(dp), pointer :: LMat(:), LMatAB(:)
  integer :: LMatTag, LMatABTag
  integer(MPIarg) :: LMatWin, LMatABWin
  integer(int64) :: nBI
  logical :: tDebugLMat

  ! for the symmetry broken index function
  integer :: strideInner, strideOuter

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
                matel = matel + sgn * real(LMat(LMatInd(int(ida,int64),int(idb,int64),&
                     int(idc,int64),int(idp,int64),int(idq,int64),int(idr,int64))),dp)
             else
                matel = matel + sgn * real(LMatAB(LMatInd(int(ida,int64),int(idb,int64),&
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
      call readLMatArray(LMat,LMatTag,LMatWin,"TCDUMP","tcdump.h5")
      ! for spin-dependent LMat, also read the opp. spin matrix
      if(tSpinCorrelator) call readLMatArray(LMatAB,LMatABTag,LMatABWin,"TCDUMPAB","tcdumpAB.h5")
    end subroutine readLMat

!------------------------------------------------------------------------------------------!

    subroutine readLMatArray(LMatLoc, LMatLocTag, LMatLocWin, LMatFileName, h5Filename)
      use ParallelHelper, only: mpi_comm_intra
      use HElem, only: Helement_t_sizeB
      implicit none
      HElement_t(dp), pointer, intent(inout) :: LMatLoc(:)
      integer, intent(inout) :: LMatLocTag
      integer(MPIArg), intent(inout) :: LMatLocWin
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


      ! The size is given by the largest index (LMatInd is monotonous in all arguments)
      LMatSize = LMatInd(nBI,nBI,nBI,nBI,nBI,nBI)

      write(iout,*) "Allocating LMat, memory required:", LMatSize*HElement_t_sizeB/(2.0_dp**30), "GB"

      ! allocate LMat (shared memory)
      call shared_allocate_mpi(LMatLocWin, LMatLoc, (/LMatSize/))
      call LogMemAlloc("LMat", int(LMatSize), HElement_t_SizeB, t_r, LMatLocTag)
!      call MPI_Comm_Size(mpi_comm_intra,iNodeSize,err)
!      call MPI_Comm_Rank(mpi_comm_intra,iProcInNode,err)
!      chunkSize = LMatSize / int(iNodeSize,int64)
!      iChunk = int(iProcInNode,int64)*chunkSize
!      print *, "On proc", iProcInNode, "on node", iNodeSize, "slice", iChunk+1, iChunk+chunkSize 

      if(iProcIndex_intra.eq.0) then
         do k = 1, LMatSize
            LMatLoc(k) = 0.0_dp
         end do
      endif

      write(iout,*) "Successfully allocated LMat"

      call MPI_Barrier(mpi_comm_intra, err)

      if(tHDF5LMat) then
#ifdef __USE_HDF
         call readHDF5LMat(LMatLoc,LMatLocWin,LMatLocTag,h5filename)
#else
         call stop_all(t_r, "HDF5 integral files disabled at compile time")
#endif
      else

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

                  LMatLoc(LMatInd(a,b,c,i,j,k)) = 3.0_dp * matel
                  if(abs(matel)> 0.0_dp) counter = counter + 1
               endif

            end do

            counter = counter / 12

            write(iout, *), "Sparsity of LMat", real(counter)/real(LMatSize)
            write(iout, *), "Nonzero elements in LMat", counter
            write(iout, *), "Allocated size of LMat", LMatSize
         endif
      end if

      call MPI_Barrier(mpi_comm_intra, err)
      
    end subroutine readLMatArray

!------------------------------------------------------------------------------------------!

    subroutine freeLMat()
      implicit none
      character(*), parameter :: t_r = "freeLMat"

      call deallocLMatArray(LMat,LMatWin,LMatTag)
      call deallocLMatArray(LMatAB,LMatABWin,LMatABTag)

      contains 

        subroutine deallocLMatArray(LMatPtr,shm_win,tag)
          implicit none
          HElement_t(dp), pointer :: LMatPtr(:)
          integer(MPIArg) :: shm_win
          integer :: tag

          if(associated(LMat)) then
             call shared_deallocate_mpi(LMatWin,LMat)
             call LogMemDealloc(t_r, LMatTag)
             LMAt => null()
          end if
        end subroutine deallocLMatArray
      end subroutine freeLMat

!------------------------------------------------------------------------------------------!

#ifdef __USE_HDF
    subroutine readHDF5LMat(LMatPtr, filename)
      use hdf5
      use hdf5_util
      use ParallelHelper, only: mpi_comm_inter, mpiInfoNull, iProcIndex, &
           mpi_comm_intra
      use Parallel_neci, only: nProcessors
      implicit none
      HElement_t(dp), pointer :: LMatPtr(:)
      character(*), intent(in) :: filename
      integer :: proc, i
      integer(hid_t) :: err, file_id, plist_id, grp_id, ds_vals, ds_inds
      integer(hsize_t) :: nInts, rest, blocksize, blockstart, blockend, this_blocksize, countsEnd
      integer(hsize_t), allocatable :: counts(:), offsets(:)
      integer(hsize_t), allocatable :: indices(:,:), entries(:,:)
      character(*), parameter :: nm_grp = "tcdump", nm_nInts = "nInts", nm_vals = "values", nm_indices = "indices"
      integer(MPIArg) :: procs_per_node, ierr
      real(dp) :: rVal
      logical :: running, any_running
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

      ! how many entries does each proc get?
      call MPI_Comm_Size(mpi_comm_intra, procs_per_node, ierr)
      allocate(counts(0:procs_per_node-1))
      allocate(offsets(0:procs_per_node-1))
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
      blocksize = 100000000/(7*sizeof(LMatPtr(1)))
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

         ! assign LMat entries
         do i = 1, this_blocksize
            LMatPtr(LMatInd(int(indices(1,i),int64),int(indices(2,i),int64),int(indices(3,i),int64),&
                 int(indices(4,i),int64),int(indices(5,i),int64),int(indices(6,i),int64))) &
                 = 3.0_dp * transfer(entries(1,i),rVal)
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
      
    end subroutine readHDF5LMat
#endif

end module LMat_mod
