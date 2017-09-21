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
      use DetBitOps, only: DetBitEQ
      use bit_rep_data, only: extract_sign
      implicit none
      integer, intent(in) :: nRefs
      integer(n_int), intent(in) :: mpi_buf(0:NIfTot,nRefs)
      integer :: i, k, run
      real(dp) :: minOcc, tmp_sgn(lenof_sign)

      ! Check if the current reference is in mpi_buf
      ! k is the index of the current reference (-1 if not present)
      k = -1
      do i = 1, nRefs
         if(DetBitEQ(ilutRef(:,1,1),mpi_buf(:,i),NIfDBO)) k = i
      enddo
      ! TODO: Sort the references here
      if(k .ne. -1) then
         do run = 1, inum_runs
            ! We still do this to equalize the reference across replicas, 
            ! if necessary
            ilutRef(:,run,1) = mpi_buf(:,k)
            do i = 2, nRefs
               if(i==k) then
                  ilutRef(:,run,i) = mpi_buf(:,1)
               else
                  ilutRef(:,run,i) = mpi_buf(:,i)
               endif
            enddo
         enddo
      else
         ! If the old reference is not in the new references, discard the lowest 
         ! populated of the new references
         k = 1
         call extract_sign(mpi_buf(:,1), tmp_sgn)
         minOcc = sum(abs(tmp_sgn))
         do i = 2, nRefs
            call extract_sign(mpi_buf(:,i),tmp_sgn)
            if(sum(abs(tmp_sgn)) < minOcc) k = i
         enddo
         ! k is now the index of the lowest populated det in mpi_buf
         do run = 1, inum_runs
            ! For now, all replicas have the same reference
            if(run > 1) ilutRef(:,run,1) = ilutRef(:,1,1)
            do i = 2, nRefs
               if(i .ne. k) ilutRef(:,run,i) = mpi_buf(:,i)
               ! If k!=1, then we take the first entry, else we never hit this point
               if(i .eq. k) ilutRef(:,run,i) = mpi_buf(:,1)
            enddo
         enddo
      endif

    end subroutine set_additional_references

  !------------------------------------------------------------------------------------------!

    subroutine output_reference_space(filename)
      use ParallelHelper, only: root
      use util_mod, only: get_free_unit
      use FciMCData, only: nRefs
      implicit none
      character(255), intent(in) :: filename
      integer :: iunit, i, k, ierr
      
      if(iProcIndex == root) then
         iunit = get_free_unit()
         open(iunit, file = filename, status = 'replace')
         ! write the references into the file
         do i = 1, nRefs
            do k = 1, NIfTot
               write(iunit, '(i24)', advance = 'no') ilutRef(k,1,i)
            enddo
            ! newline
            write(iunit, '()')
         enddo
         call neci_flush(iunit)
         close(iunit)
      endif
      call MPIBarrier(ierr)
    end subroutine output_reference_space

end module adi_references
