module LMat_mod
  use constants
  use FciMCData, only: ll_node
  use HElem, only: HElement_t_SizeB
  use SystemData, only: nBasis, t12FoldSym, G1, t_mol_3_body, nel, nBI, tStoreSpinOrbs, tContact
  use MemoryManager, only: LogMemAlloc, LogMemDealloc
  use util_mod, only: get_free_unit, fuseIndex
  use gen_coul_ueg_mod, only: get_lmat_ueg, get_lmat_ua
  use shared_memory_mpi
  use sort_mod
  use hash, only: add_hash_table_entry, clear_hash_table
  use ParallelHelper, only: iProcIndex_intra
  use tc_three_body_data, only: tDampKMat, tDampLMat, tSpinCorrelator, lMatEps, &
       tHDF5LMat, tSymBrokenLMat, tSparseLMat, lMat_t, tLMatCalc
  use procedure_pointers, only: get_lmat_el, get_lmat_el_symInternal
  use LoggingData, only: tHistLMat
  use LMat_aux, only: diffSpinPos, dampLMatel
  use LMat_indexing, only: lMatIndSym, lMatIndSymBroken, oldLMatInd, strideInner, strideOuter, &
       lMatIndSpin
  use LMat_calc, only: readlMatFactors, freelMatFactors, lMatCalc, lMatABCalc
#ifdef __USE_HDF5
  use hdf5
#endif
  implicit none

  ! access function for lMat
  abstract interface
     function lMatAccess_t(lMatObj, index) result(matel)
       use constants
       use tc_three_body_data, only: lMat_t
       implicit none
       type(lMat_t), intent(in) :: lMatObj
       integer(int64), intent(in) :: index
       HElement_t(dp) :: matel
     end function lMatAccess_t
     
  end interface

  ! actual objects storing the 6-index integrals
  type(lMat_t), target :: LMat, LMatAB

  procedure(lMatAccess_t), pointer :: lMatAccess

  contains

!------------------------------------------------------------------------------------------!    
! Access function for six-index integrals: get a matrix element given the spin-orbitals
!------------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------------!  
    
    ! this is the common logic of all 6-index lmat-acceses
    function get_lmat_el_base(a,b,c,i,j,k) result(matel)
      ! Input: a,b,c - indices of orbitals to excite to
      !        i,j,k - indices of orbitals to excite from
      ! Output: matel - matrix element of this excitation, including all exchange terms
      use UMatCache, only: gtID
      ! Gets an entry of the 3-body tensor L:
      ! L_{abc}^{ijk} - triple excitation from abc to ijk
      implicit none
      integer, value :: a,b,c
      integer, value :: i,j,k
      HElement_t(dp) :: matel
      integer(int64) :: ida, idb, idc, idi, idj, idk
      logical :: tSameSpin

      ! initialize spin-correlator check: if all spins are the same, use LMat
      ! without spin-dependent correlator, always use LMat
      
      ! for matrix elements involving different spins, there are three cases:
      ! each of the orbitals i,j,k can have the different spin
      ! The position is important because the spinCorrelator breaks permutation symmetry
      ! w.r.t spin -> we fix the differing spin
      tSameSpin = .not. tSpinCorrelator .or. (G1(a)%MS==G1(b)%ms .and. G1(a)%MS==G1(c)%MS)
      
      ! convert to spatial orbs if required
      ida = gtID(a)
      idb = gtID(b)
      idc = gtID(c)
      idi = gtID(i)
      idj = gtID(j)
      idk = gtID(k)

      matel = 0
      ! only add the contribution if the spins match

      ! TODO: Use sameSpin to determine which contributions can appear - no need for
      ! any further checks (remember to tweak for sameSpin!= possible without tSpinCorrelator)
      ! here, we add up all the exchange terms
      call addMatelContribution(i,j,k,idi,idj,idk,1)
      call addMatelContribution(j,k,i,idj,idk,idi,1)
      call addMatelContribution(k,i,j,idk,idi,idj,1)
      call addMatelContribution(j,i,k,idj,idi,idk,-1)
      call addMatelContribution(i,k,j,idi,idk,idj,-1)
      call addMatelContribution(k,j,i,idk,idj,idi,-1)
      ! if a heuristic spin-projection is done, it happens here
      if(tDampKMat .and. .not. tDampLMat) call dampLMatel(a,b,c,matel)

      contains 

        subroutine addMatelContribution(p,q,r,idp,idq,idr,sgn)
          ! get a single entry of the LMat array and add it to the matrix element
          implicit none
          integer(int64), value :: idp,idq,idr
          integer, value :: p,q,r
          integer, intent(in) :: sgn
          integer(int64) :: index
          integer :: spinPos
          type(lMat_t), pointer :: lMatPtr
          real(dp) :: lMatVal
     !     integer(int64) :: ai,bj,ck
          
          if(G1(p)%ms == G1(a)%ms .and. G1(q)%ms == G1(b)%ms .and. G1(r)%ms == G1(c)%ms) then

             if(tContact) then
                    lMatVal = get_lmat_ueg(ida,idb,idc,idp,idq,idr)
             else if(tLMatCalc)then
                if(tSameSpin) then
                    lMatVal = lMatCalc(ida,idb,idc,idp,idq,idr)
                else
                    spinPos = diffSpinPos(p,q,r,a,b,c)
                    lMatVal = lMatABCalc(ida,idb,idc,idp,idq,idr, spinPos)
                end if
             else
                 ! pick the lMat object used here according to the spin-relation
                 if(tSameSpin) then                
                    lMatPtr => lMat                
                    ! the indexing function is contained in the lMat object                      
                    index = lMatPtr%indexFunc(ida,idb,idc,idp,idq,idr)                
                 else
                    ! for different spins, check which one is the different one and
                    ! call the index function accordingly
                    lMatPtr => lMatAB
                    spinPos = diffSpinPos(p,q,r,a,b,c)
                    ! the different-spin LMat assumes the first electron in the
                    ! index function has the different spin (the order of the other two does
                    ! not matter)
                    select case(spinPos)
                    case(1)
                       index = lMatPtr%indexFunc(ida,idb,idc,idp,idq,idr)
                    case(2)
                       index = lMatPtr%indexFunc(idb,ida,idc,idq,idp,idr)
                    case(3)
                       index = lMatPtr%indexFunc(idc,idb,ida,idr,idq,idp)
                    end select
                 endif
                 lMatVal = real(lMatAccess(lMatPtr,index),dp)
             endif
             matel = matel + sgn * lMatVal
          endif

        end subroutine addMatelContribution
        
      end function get_lmat_el_base

!------------------------------------------------------------------------------------------!

      function get_lmat_el_symmetrized(a,b,c,i,j,k) result(matel)
        ! post-symmetrized access when storing non-symmetrized lMat
        implicit none
        integer, value :: a,b,c
        integer, value :: i,j,k
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
        integer, value :: i,j,k
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
        matel = matel + permFac * get_lmat_el_base(ap,bp,cp,ip,jp,kp)
      
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
! Auxiliary functions for indexing and accessing the LMat    
!------------------------------------------------------------------------------------------!

    subroutine initializeLMatPtrs()
      implicit none

      ! some typical array dimensions useful in the indexing functions
      strideInner = fuseIndex(int(nBI),int(nBI))
      strideOuter = strideInner**2
      ! set the LMatInd function pointer
      if(t12FoldSym) then
         lMat%indexFunc => oldLMatInd
      else if(tSymBrokenLMat) then
         ! also need to set the size of the blocks
         lMat%indexFunc => lMatIndSymBroken
      else
         lMat%indexFunc => lMatIndSym
      endif
      ! set the spin-correlator index functions
      lMatAB%indexFunc => lMatIndSpin

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
    end subroutine initializeLMatPtrs

!------------------------------------------------------------------------------------------!
! For huge LMats, a sparse storage scheme is required. Since we still need to access
! the matrix elements, we use a hash table and re-use the hash.F90 module
!------------------------------------------------------------------------------------------!

      subroutine initLMatHash(lMatCtr)
        implicit none
        type(LMat_t), intent(inout) :: lMatCtr
        integer :: i
        integer :: hashVal

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
        integer :: hashVal

! TODO: Implement an actual hash function
        hashVal = mod(int(index)-1,htSize)+1
      end function LMatHash

!------------------------------------------------------------------------------------------!
! Functions that access single entries of the LMat array      
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
! Six-index integral I/O functions    
!------------------------------------------------------------------------------------------!

    subroutine readLMat()
      use SystemData, only: nel
      implicit none

      ! we need at least three electrons to make use of the six-index integrals
      ! => for less electrons, this can be skipped
      if(nel<=2) return

      call initializeLMatPtrs()

      if(tLMatCalc) then
           call readLMatFactors("tcfactors.h5")
      else
          ! now, read lmat from file
          call readLMatArray(LMat,"TCDUMP","tcdump.h5")
          ! for spin-dependent LMat, also read the opp. spin matrices 
          if(tSpinCorrelator) then
             ! they break permutational symmetry w.r.t spin, so the cheapest solution is
             ! to have three instances, one for each spin-permutation
             ! (the alternatives are: spin-orbitals with full symmetry or spatial orbitals
             ! without spin symmetry - both more expensive)
             call readLMatArray(LMatAB,"TCDUMPAB","tcdumpab.h5")
          end if
      end if
    end subroutine readLMat

!------------------------------------------------------------------------------------------!

    subroutine readLMatArray(LMatLoc, LMatFileName, h5Filename)
      implicit none
      type(lMat_t), intent(inout) :: LMatLoc
      character(*), intent(in) :: LMatFileName
      character(*), intent(in) :: h5Filename
      integer :: iunit, ierr
      integer(int64) :: a,b,c,i,j,k
      HElement_t(dp) :: matel
      integer(int64) :: LMatSize      
      character(*), parameter :: t_r = "readLMat"
      integer(int64) :: counter

      if(.not.tSparseLMat) then
         ! for sparse storage, we first need to get the number of integrals
         ! this works differently for hdf5 and non-hdf5 dumpfiles, so it is
         ! done directly in the respective code

         ! The size is given by the largest index (LMatInd is monotonous in all arguments)
         LMatSize = lMatLoc%indexFunc(nBI,nBI,nBI,nBI,nBI,nBI)

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
                  if(lMatLoc%indexFunc(a,b,c,i,j,k) > LMatSize) then
                     counter = lMatLoc%indexFunc(a,b,c,i,j,k)
                     write(iout,*) "Warning, exceeding size" 
                  endif
                  if(abs(3.0_dp*matel) > LMatEps) &
                       LMatLoc%LMatPtr(lMatLoc%indexFunc(a,b,c,i,j,k)) = 3.0_dp * matel
                  if(abs(matel)> 0.0_dp) counter = counter + 1
               endif

            end do

            counter = counter / 12

            write(iout, *), "Sparsity of LMat", real(counter)/real(LMatSize)
            write(iout, *), "Nonzero elements in LMat", counter
            write(iout, *), "Allocated size of LMat", LMatSize
         endif
         call MPIBcast(counter)
         LMatLoc%nInts = counter
      end if

      if(tHistLMat) call histogramLMat(LMatLoc)

    end subroutine readLMatArray

!------------------------------------------------------------------------------------------!
! Memory management for the LMat_t objects    
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

      if(tSparseLMat) then
         ! also, allocate the index array
         call shared_allocate_mpi(LMatLoc%index_win, LMatLoc%indexPtr, (/LMatSize/))
         call LogMemAlloc("LMat Indices", int(LMatSize), sizeof_int64, t_r, LMatLoc%indexTag)
         ! come up with some reasonable size
         LMatLoc%htSize = LMatSize
      endif

      write(iout,*) "Successfully allocated LMat"
    end subroutine allocLMat

!------------------------------------------------------------------------------------------!

    subroutine freeLMat()
      implicit none
      character(*), parameter :: t_r = "freeLMat"

      if(tLMatCalc) then
          call freeLMatFactors()
      else
          call deallocLMatArray(LMatAB)
          call deallocLMatArray(LMat)
      end if

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
      call h5pset_fapl_mpio_f(plist_id, mpi_comm_intra, mpiInfoNull, err)

      ! open the file
      call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, err, access_prp=plist_id)

      call h5gopen_f(file_id, nm_grp, grp_id, err)

      ! get the number of integrals
      call read_int64_attribute(grp_id, nm_nInts, nInts, required=.true.)
      write(iout,*) "Reading", nInts, "integrals"

      ! if the LMat is stored sparse, we can allocate now
      if(tSparseLMat) then
         call allocLMat(LMatLoc,nInts)
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
      if(iProcIndex_intra.eq.procs_per_node - 1) then
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
                  LMatLoc%indexPtr(sparseBlock+counter) = lMatLoc%indexFunc(int(indices(1,i),int64),int(indices(2,i),int64),&
                       int(indices(3,i),int64),&
                       int(indices(4,i),int64),int(indices(5,i),int64),int(indices(6,i),int64))
               else
                  LMatLoc%LMatPtr(lMatLoc%indexFunc(int(indices(1,i),int64),int(indices(2,i),int64),&
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

         ! once all procs on this node are done reading, we can exit
         call MPI_ALLREDUCE(running, any_running, 1, MPI_LOGICAL, MPI_LOR, mpi_comm_intra, ierr)

         ! communicate how many nonzero entries have been read and set the starting point 
         ! for the next write
         if(tSparseLMat) then
            call MPI_ALLREDUCE(counter,allCounter,1,MPI_INT,MPI_SUM,mpi_comm_intra,ierr)
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
!------------------------------------------------------------------------------------------------
!functions for contact interaction

    function get_lmat_el_ua(a,b,c,i,j,k) result(matel)
      use SystemData, only: G1
      use UMatCache, only: gtID
      ! Gets an entry of the 3-body tensor L:
      ! L_{abc}^{ijk} - triple excitation from abc to ijk
      implicit none
      integer, value :: a,b,c
      integer :: a2,b2,c2
      integer, intent(in) :: i,j,k
      HElement_t(dp) :: matel

      ! convert to spatial orbs if required

      matel = 0
      
      if(G1(a)%ms == G1(b)%ms .and. G1(a)%ms.ne.G1(c)%ms ) then
         a2=a
         b2=b
         c2=c
      elseif(G1(a)%ms == G1(c)%ms .and. G1(a)%ms.ne.G1(b)%ms ) then
         a2=c
         b2=a
         c2=b
      elseif(G1(b)%ms == G1(c)%ms .and. G1(a)%ms.ne.G1(b)%ms ) then
         a2=b
         b2=c
         c2=a
      else
        return
      endif

      ! only add the contribution if the spins match
         call addMatelContribution_ua(i,j,k,1)
         call addMatelContribution_ua(j,k,i,1)
         call addMatelContribution_ua(k,i,j,1)
         call addMatelContribution_ua(j,i,k,-1)
         call addMatelContribution_ua(i,k,j,-1)
         call addMatelContribution_ua(k,j,i,-1)
      contains

        subroutine addMatelContribution_ua(p,q,r,sgn)
          implicit none
          integer, value :: p,q,r
          integer, intent(in) :: sgn
     !     integer(int64) :: ai,bj,ck

          if(G1(p)%ms == G1(a2)%ms .and. G1(q)%ms == G1(b2)%ms .and. G1(r)%ms ==G1(c2)%ms) then
             matel = matel + 2.d0 * sgn * get_lmat_ua(a2,b2,c2,p,q,r)
          endif

        end subroutine addMatelContribution_ua

    end function get_lmat_el_ua

!------------------------------------------------------------------------------------------!

    subroutine histogramLMat(lMatObj)
      implicit none
      type(lMat_t), intent(in) :: lMatObj
      integer(int64) :: i
      integer :: thresh
      integer, parameter :: minExp = 10
      integer :: histogram(0:minExp)
      real :: ratios(0:minExp)
      
      histogram = 0
      do i = 1, size(lMatObj%LMatPtr)
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

!------------------------------------------------------------------------------------------!
    
    function getPCWeightOffset() result(offset)
      ! get the offset used for generating precomputed weights
      ! Output: offset - constant to be added to each weight of the precomputed-excitation generators
      implicit none
      real(dp) :: offset
      integer :: i
      real(dp) :: maxInts(nel-2)

      if(t_mol_3_body.and.nel>2) then
         maxInts = 0.0_dp
         write(iout,*) "Scanning", lMat%nInts, "integrals"
         ! we estimate the 6-index contribution with the maximum possible contribution
         do i = 1, size(lMat%lMatPtr)
            if(abs(lMat%LMatPtr(i)) > minval(maxInts)) then
               maxInts(minloc(maxInts)) = abs(lMat%LMatPtr(i))
            endif
         end do
         offset = sum(maxInts)
         write(iout,*) "Excitgen is offset of", offset
      else
         offset = 0.0_dp
      endif
    end function getPCWeightOffset

!------------------------------------------------------------------------------------------!

    ! last resort in debugging: Write out the lmat from memory
    subroutine write_lmat_debug()
      implicit none
      integer :: i
      character(*), parameter :: lmatfilename = "LMAT_FILE"
      integer :: iunit
      
      iunit = get_free_unit()
      open(iunit,file=lmatfilename,status='UNKNOWN')
      do i = 1, size(LMat%LMatPtr)
         write(iunit,*) LMat%LMatPtr(i)
      end do
      close(iunit)
      
    end subroutine write_lmat_debug

end module LMat_mod

