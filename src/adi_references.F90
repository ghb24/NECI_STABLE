module adi_references
use Parallel_neci
use FciMCData, only: ilutRefAdi, ilutRef, nRefsCurrent
use bit_rep_data, only: niftot, nifdbo
use constants

contains

  subroutine read_in_refs(nRefs, filename)
    use util_mod, only: get_free_unit
    implicit none
    integer, intent(in) :: nRefs
    character(255), intent(in) :: filename
    integer :: iunit, i, run, stat
    logical :: exists
    character(*), parameter :: this_routine = "read_in_refs"

    ! All procs read in the file, as all need all references
    inquire(file=trim(filename), exist = exists)
    if(.not. exists) call stop_all(this_routine, "No "//trim(filename)//" file detected.")
    iunit = get_free_unit()
    open(iunit, file=trim(filename), status = 'old')
    ! Read in at most nRefs references
    do i = 1, nRefs
       ! Note that read-in references always have precedence over generated references
       read(iunit, *, iostat=stat) ilutRefAdi(:,1,i)
       ! If there are no more dets to be read, exit
       if(stat < 0) exit
    enddo

    ! Copy the references to all runs
    do i = 1, nRefs
       do run = 2, inum_runs
          ilutRefAdi(:,run,i) = ilutRefAdi(:,1,i)
       end do
    enddo

    close(iunit)
    nRefsCurrent = nRefs
  end subroutine read_in_refs

!------------------------------------------------------------------------------------------!

    subroutine generate_ref_space(nRefs)
      use semi_stoch_gen, only: generate_space_most_populated
      use LoggingData, only: ref_filename, tWriteRefs
      use CalcData, only: tProductReferences, nExProd
      use SystemData, only: tHPHF
      implicit none
      integer, intent(inout) :: nRefs
      integer(MPIArg) :: mpi_refs_found
      integer :: ierr, i, all_refs_found, refs_found
      integer(MPIArg) :: refs_found_per_proc(0:nProcessors-1), refs_displs(0:nProcessors-1)
      integer(n_int) :: ref_buf(0:NIfTot,nRefs), mpi_buf(0:NIfTot,nRefs)
      character(*), parameter :: this_routine = "generate_ref_space"
      
      write(6,*) "Getting reference determinants for all-doubs-initiators"

      ! Get the nRefs most populated determinants
      refs_found = 0
      call generate_space_most_populated(nRefs, .false., 1, ref_buf, refs_found)
      ! Communicate the refs_found info
      mpi_refs_found = int(refs_found,MPIArg)
      call MPIAllGather(mpi_refs_found, refs_found_per_proc, ierr)
      all_refs_found = sum(refs_found_per_proc)
      if(all_refs_found .ne. nRefs) then
         write(6,*) "all_refs_found = ", all_refs_found
         call stop_all(this_routine, &
           "Number of references found differs from requested")
      endif
      refs_displs(0) = 0
      do i = 1, nProcessors - 1
         refs_displs(i) = sum(refs_found_per_proc(0:i-1))
      enddo
      ! Store them on all processors
      call MPIAllGatherV(ref_buf(0:NIfTot, 1:refs_found), mpi_buf, refs_found_per_proc, refs_displs)

      ! And write the so-merged references to ilutRef
      call set_additional_references(mpi_buf, nRefs)

      if(tProductReferences) call add_product_references(nRefs, nExProd)
      if(tHPHF) call spin_symmetrize_references(nRefs)

      if(tWriteRefs) call output_reference_space(ref_filename)
      nRefsCurrent = nRefs
    end subroutine generate_ref_space

!------------------------------------------------------------------------------------------!

    subroutine set_additional_references(mpi_buf,nRefs)
      use DetBitOps, only: DetBitEQ, sign_lt, sign_gt
      use sort_mod, only: sort
      use bit_rep_data, only: extract_sign
      implicit none
      integer, intent(in) :: nRefs
      integer(n_int), intent(inout) :: mpi_buf(0:NIfTot,nRefs)
      integer :: i, k, run
      integer(n_int) :: tmp_ilut(0:NIfTot)
      real(dp) :: minOcc, tmp_sgn(lenof_sign)

      ! Check if the current reference is in mpi_buf
      ! k is the index of the current reference (-1 if not present)
      k = -1
      do i = 1, nRefs
         if(DetBitEQ(ilutRef(:,1),mpi_buf(:,i),NIfDBO)) k = i
      enddo
      
      ! If the reference is not at pos. 1, move it there
      if(k .ne. -1) then 
         ! If k==1, everything is already alright, nothing to do here
         if(k > 1) then
            tmp_ilut = mpi_buf(:,k)
            mpi_buf(:,k) = mpi_buf(:,1)
            mpi_buf(:,1) = tmp_ilut
         endif
      else
         ! if it is not even in mpi_buf, replace the lowest populated det
         ! in the buffer with it
         k = 1
         call extract_sign(mpi_buf(:,1), tmp_sgn)
         minOcc = sum(abs(tmp_sgn))
         do i = 2, nRefs
            call extract_sign(mpi_buf(:,i),tmp_sgn)
            if(sum(abs(tmp_sgn)) < minOcc) k = i
         enddo
         tmp_ilut = mpi_buf(:,1)
         mpi_buf(:,1) = ilutRef(:,1)
         if(k > 1) mpi_buf(:,k) = tmp_ilut
      endif
      ! Now we have the main reference at position 1

      ! Then sort the buffer with respect to sign (but keep the original reference at pos. 1
      call sort(mpi_buf(:,2:), sign_gt, sign_lt)

      ! Copy the buffer to ilutRef
      do i = 1, nRefs
         do run = 1, inum_runs
            ilutRefAdi(:,run,i) = mpi_buf(:,i)
         enddo
      enddo
    end subroutine set_additional_references

!------------------------------------------------------------------------------------------!

    subroutine output_reference_space(filename)
      use ParallelHelper, only: root
      use util_mod, only: get_free_unit
      use FciMCData, only: nRefs
      implicit none
      character(255), intent(in) :: filename
      integer :: iunit, i, ierr
      
      if(iProcIndex == root) then
         iunit = get_free_unit()
         open(iunit, file = filename, status = 'replace')
         ! write the references into the file
         do i = 1, nRefs
            write(iunit, *) ilutRefAdi(:,1,i)
         enddo
         call neci_flush(iunit)
         close(iunit)
      endif
      call MPIBarrier(ierr)
    end subroutine output_reference_space

!------------------------------------------------------------------------------------------!

    subroutine print_reference_notification(nRefs)
      use bit_rep_data, only: extract_sign
      implicit none
      integer, intent(in) :: nRefs
      integer :: i, j
      real(dp) :: temp_sgn(lenof_sign)
      
      if(iProcIndex==root) then
         print *, "References for all-doubs-initiators set as"
         do i=1, nRefs
            ! First output the reference determinant
            call WriteDetBit(iout,ilutRefAdi(:,1,i),.false.)
            ! And then the sign
            call extract_sign(ilutRefAdi(:,1,i),temp_sgn)
            do j = 1, lenof_sign
               write(iout,"(G16.7)",advance='no') temp_sgn(j)
            enddo
            write(iout,'()')
         enddo
      endif
    end subroutine print_reference_notification

!------------------------------------------------------------------------------------------!

    subroutine spin_symmetrize_references(nRefs)
      use DetBitOps, only: spin_flip, DetBitEQ
      ! If using HPHF, we want to have the spinflipped version of each reference, too
      ! This is a bit cumbersome to implement, as the spinflipped version might or
      ! might not be already in the references
      implicit none
      integer, intent(inout) :: nRefs
      integer :: iRef, jRef, nRefs_new, cRef, run
      integer(n_int) :: ilutRef_new(0:NIfTot,2*nRefs), tmp_ilut(0:NIfTot)
      logical :: missing
      
      cRef = 1
      do iRef = 1 ,nRefs
         ! First, copy the current ilut to a new buffer
         ilutRef_new(:,cRef) = ilutRefAdi(:,1,iRef)
         cRef = cRef + 1
         ! get the spinflipped determinant
         tmp_ilut = spin_flip(ilutRefAdi(:,1,iRef))
         missing = .true.
         ! check if it is already present
         do jRef = 1, nRefs
            if(DetBitEQ(tmp_ilut,ilutRefAdi(:,1,jRef),NIfDBO)) missing = .false.
         end do
         ! If not, also add it to the list
         if(missing) then
            ilutRef_new(:,cRef) = tmp_ilut
            cRef = cRef + 1
         endif
      enddo

      ! Resize the ilutRef array
      deallocate(ilutRefAdi)
      allocate(ilutRefAdi(0:NIfTot,inum_runs,cRef))
      nRefs = cRef
      ! Now, copy the newly constructed references back to the original array
      do iRef = 1, cRef
         do run = 1, inum_runs
            ilutRefAdi(:,run,iRef) = ilutRef_new(:,iRef)
         end do
      end do
    end subroutine spin_symmetrize_references

!------------------------------------------------------------------------------------------!

    subroutine add_product_references(nRefs, prodLvl)
      implicit none
      integer, intent(inout) :: nRefs
      integer, intent(in) :: prodLvl
      integer(n_int) :: ilut_tmp(0:NIfTot)
      integer(n_int), allocatable :: prod_buffer(:,:)
      integer :: nProdsMax, nProds, i, run, cLvl
      logical :: t_is_valid

      ! Get the maximum number of product excitations
      nProdsMax = maxProdEx(nRefs, prodLvl)

      ! Setup the buffer
      allocate(prod_buffer(0:NifTot,nProdsMax))
      prod_buffer(:,1:nRefs) = ilutRefAdi(:,1,1:nRefs)
      nProds = nRefs

      do cLvl = 2, prodLvl
         ! For each product level, get all possible product excitations
         do i = 1, nRefs**cLvl
            ! Get the excitation product corresponding to i
            call get_product_excitation(i, cLvl, nRefs, ilut_tmp, t_is_valid)
            
            ! Check, if it is already in the list (e.g. if both tau operators are those of
            ! the reference)
            if(t_is_valid) t_is_valid = ilut_not_in_list(ilut_tmp, prod_buffer, nProds)
         
            ! If all excitations were valid and the new ilut is not already in the buffer, 
            ! add it
            if(t_is_valid) then
               nProds = nProds + 1
               prod_buffer(:,nProds) = ilut_tmp
            endif
         end do
         ! Do this for each targeted product level
      end do
            
      ! Write the buffer to ilutRefAdi
      deallocate(ilutRefAdi)
      allocate(ilutRefAdi(0:NIfTot,inum_runs,nProds))
      ! Reassign the number of references
      nRefs = nProds
      do run = 1, inum_runs
         ilutRefAdi(:,run,:) = prod_buffer(:,1:nProds)
      end do
    end subroutine add_product_references

!------------------------------------------------------------------------------------------!

    subroutine get_product_excitation(i, cLvl, nRefs, ilut_tmp, isValid)
      use DetBitOps, only: CountBits
      use bit_rep_data, only: NIfD
      use SystemData, only: nEl
      implicit none
      integer, intent(in) :: i, cLvl, nRefs
      integer(n_int), intent(out) :: ilut_tmp(0:NIfTot)
      logical, intent(out) :: isValid
      integer(n_int) :: tau(0:NIfTot)
      integer :: cp, tmp

      ilut_tmp = ilutRef(:,1)
      isValid = .true.
      ! By getting the cLvl excitation operators that correspond to i
      do cp = 1, cLvl
         ! Store the excitation operators in tau
         tau = IEOR(ilutRefAdi(:,1,cpIndex(i,nRefs,cp)),ilutRef(:,1))
         tmp = cpIndex(i,nRefs,cp)
         ! And apply it on the temporary iLut
         ilut_tmp = IEOR(ilut_tmp,tau)

         ! Check if it's valid, if not, go to the next value of i
         if(CountBits(ilut_tmp,NIfD) .ne. nEl) then
            isValid = .false.
            exit
         endif

      end do
    end subroutine get_product_excitation

!------------------------------------------------------------------------------------------!

    pure function ilut_not_in_list(ilut, buffer, buffer_size) result(t_is_new)
      use DetBitOps, only: DetBitEQ
      use bit_rep_data, only: NIfD
      implicit none
      integer, intent(in) :: buffer_size
      integer(n_int), intent(in) :: ilut(0:NIfTot), buffer(0:NIfTot,buffer_size)
      logical :: t_is_new
      integer :: i_ilut

      t_is_new = .true.     
      do i_ilut = 1, buffer_size
         if(DetBitEQ(ilut, buffer(:,i_ilut), NIfD)) t_is_new = .false.
      enddo
    end function ilut_not_in_list

!------------------------------------------------------------------------------------------!

    function cpIndex(i,nRefs,cp) result(index)
      implicit none
      integer, intent(in) :: i, nRefs, cp
      integer :: index
      
      index = mod(i/(nRefs**(cp-1)),nRefs)+1
    end function cpIndex

!------------------------------------------------------------------------------------------!

    function maxProdEx(nRefs, prodLvl) result(nProdsMax)
      implicit none
      integer, intent(in) :: nRefs, prodLvl
      integer :: nProdsMax
      integer :: cLvl

      ! Gets the maximum number of possible product excitations of nRefs states up 
      ! to a product level of prodLvl
      nProdsMax = 0
      do cLvl = 2, prodLvl
         nProdsMax = nProdsMax + nRefs**cLvl
      enddo

    end function maxProdEx

!------------------------------------------------------------------------------------------!

    subroutine enable_adi
      use CalcData, only: tAllDoubsInitiators
      use FciMCData, only: nRefsCurrent
      implicit none
      integer :: i

      write(iout,'()') 
      write(iout,*) "Setting all double excitations to initiators"
      write(iout,*) "Using static references"
      write(iout,*) "Further notification on additional references will be given"
      write(iout,'()')
      tAllDoubsInitiators = .true.
      
      ! tSinglePartPhase = .false.
      
    end subroutine enable_adi

!------------------------------------------------------------------------------------------!

    function giovannis_check(ilut) result(init)
      use CalcData, only: g_markers, g_markers_num
      use DetBitOps, only: CountBits
      use bit_rep_data, only: NIfD
      ! Very much like DetBitOps::FindBitExcitLevel, but it does the check for the active-space
      ! initiator criterium
      implicit none
      integer(n_int), intent(in) :: ilut(0:NIfD)
      integer(n_int) :: ilut_tmp(0:NIfD)
      integer :: nOrbs
      logical :: init

      init = .false.
      ! Get the bitwise equivalence of the input with the reference
      ilut_tmp = NOT(IEOR(ilut,ilutRef(:,1)))
      ! Now compare ilut_tmp with the markers (these are 0 for the core orbitals and 1 else)
      ilut_tmp = IAND(ilut_tmp, g_markers)
      nOrbs = CountBits(ilut_tmp,NIfD)
      ! If nOrbs is the markersize, it is an initiator
      init = (nOrbs == g_markers_num)

    end function giovannis_check
    
end module adi_references
