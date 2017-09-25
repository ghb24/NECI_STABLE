module adi_references
use Parallel_neci
use FciMCData, only: ilutRef
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
       read(iunit, *, iostat=stat) ilutRef(:,1,i)
       ! If there are no more dets to be read, exit
       if(stat < 0) exit
    enddo

    ! Copy the references to all runs
    do i = 1, nRefs
       do run = 2, inum_runs
          ilutRef(:,run,i) = ilutRef(:,1,i)
       end do
    enddo

    close(iunit)
  end subroutine read_in_refs

!------------------------------------------------------------------------------------------!

    subroutine generate_ref_space(nRefs)
      use semi_stoch_gen, only: generate_space_most_populated
      use LoggingData, only: ref_filename, tWriteRefs
      implicit none
      integer, intent(in) :: nRefs
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

      if(tWriteRefs) call output_reference_space(ref_filename)
            
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
         if(DetBitEQ(ilutRef(:,1,1),mpi_buf(:,i),NIfDBO)) k = i
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
         mpi_buf(:,1) = ilutRef(:,1,1)
         if(k > 1) mpi_buf(:,k) = tmp_ilut
      endif
      ! Now we have the main reference at position 1

      ! Then sort the buffer with respect to sign (but keep the original reference at pos. 1
      call sort(mpi_buf(:,2:), sign_gt, sign_lt)

      ! Copy the buffer to ilutRef
      do i = 1, nRefs
         do run = 1, inum_runs
            ilutRef(:,run,i) = mpi_buf(:,i)
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
            write(iunit, *) ilutRef(:,1,i)
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
            call WriteDetBit(iout,ilutRef(:,1,i),.false.)
            ! And then the sign
            call extract_sign(ilutRef(:,1,i),temp_sgn)
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
      integer, intent(in) :: nRefs
      integer :: iRef, jRef, nRefs_new, cRef, run
      integer(n_int) :: ilutRef_new(0:NIfTot,2*nRefs), tmp_ilut(0:NIfTot)
      logical :: missing
      
      cRef = 1
      do iRef = 1 ,nRefs
         ! First, copy the current ilut to a new buffer
         ilutRef_new(:,cRef) = ilutRef(:,1,iRef)
         cRef = cRef + 1
         ! get the spinflipped determinant
         tmp_ilut = spin_flip(ilutRef(:,1,iRef))
         missing = .true.
         ! check if it is already present
         do jRef = 1, nRefs
            if(DetBitEQ(tmp_ilut,ilutRef(:,1,jRef),NIfDBO)) missing = .false.
         end do
         ! If not, also add it to the list
         if(missing) then
            ilutRef_new(:,cRef) = tmp_ilut
            cRef = cRef + 1
         endif
      enddo

      ! Resize the ilutRef array
      deallocate(ilutRef)
      allocate(ilutRef(0:NIfTot,inum_runs,cRef))
      ! Now, copy the newly constructed references back to the original array
      do iRef = 1, cRef
         do run = 1, inum_runs
            ilutRef(:,run,iRef) = ilutRef_new(:,iRef)
         end do
      end do
    end subroutine spin_symmetrize_references

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
      ilut_tmp = NOT(IEOR(ilut,ilutRef(:,1,1)))
      ! Now compare ilut_tmp with the markers (these are 0 for the core orbitals and 1 else)
      ilut_tmp = IAND(ilut_tmp, g_markers)
      nOrbs = CountBits(ilut_tmp,NIfD)
      ! If nOrbs is the markersize, it is an initiator
      init = (nOrbs == g_markers_num)

    end function giovannis_check
    
end module adi_references
