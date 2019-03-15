module shared_memory_mpi
  use mpi
  use Parallel_neci
  implicit none

  public:: shared_allocate_mpi, shared_deallocate_mpi, shared_sync_mpi, MPIBCast_inter_byte
  private
contains

#ifdef __SHARED_MEM

      subroutine shared_allocate_mpi (win_shm, p_shm, dims)
      use HElem
      integer(MPIArg):: win_shm
      HElement_t(dp), pointer :: p_shm(:)
      integer(int64):: dims(1)

      integer(MPIArg):: disp_unit
      integer(MPIArg) :: ierr, jerr, errorclass
      integer(MPIArg) :: length
      integer(kind=mpi_address_kind):: wsize
!      character, allocatable :: string
      character(255) :: string

      TYPE(C_PTR):: cptr_shm

      if (iProcIndex_intra.eq.0) then
         wsize=dims(1)*HElement_t_sizeB
      else
         wsize=0
      end if

      call mpi_win_allocate_shared(wsize,int(HElement_t_sizeB,MPIArg),MPI_INFO_NULL,mpi_comm_intra,&
           cptr_shm,win_shm,ierr)
      if (ierr /= MPI_SUCCESS) then
          call mpi_error_class(ierr, errorclass, jerr)
          call mpi_error_string(errorclass, string, length, jerr)
          call stop_all('shared_allocate_mpi', string)
      end if

      call mpi_win_shared_query(win_shm,0_MPIArg,wsize,disp_unit,cptr_shm,ierr)

      !map to Fortran array pointer
      call c_f_pointer(cptr_shm,p_shm,dims)

      !start read/write epoch for this window
      call mpi_win_lock_all(MPI_MODE_NOCHECK,win_shm,ierr)

    end subroutine shared_allocate_mpi


    subroutine shared_deallocate_mpi(win_shm,p_shm)
      integer(MPIArg):: win_shm
      HElement_t(dp), pointer :: p_shm(:)
      integer(MPIArg):: ierr

      nullify(p_shm)
      call mpi_win_unlock_all(win_shm,ierr)
      call mpi_win_free(win_shm,ierr)
      if (ierr /= MPI_SUCCESS) then
          call stop_all('shared_deallocate_mpi', 'Could not free win.')
      end if

    end subroutine shared_deallocate_mpi

    subroutine shared_sync_mpi(win_shm)
      integer(MPIArg):: win_shm
      integer(MPIArg):: ierr

      call mpi_win_sync(win_shm,ierr)
      call mpi_barrier(mpi_comm_intra,ierr)
    end	subroutine shared_sync_mpi

#else
    subroutine shared_allocate_mpi (win_shm, p_shm, dims)
      use HElem
      integer(MPIArg):: win_shm
      HElement_t(dp), pointer :: p_shm(:)
      integer(int64):: dims(1)

      integer:: disp_unit
      integer(MPIArg) :: ierr
      integer(kind=mpi_address_kind):: wsize
      TYPE(C_PTR):: cptr_shm

      allocate(p_shm(dims(1)))

    end subroutine shared_allocate_mpi

    subroutine shared_deallocate_mpi(win_shm,p_shm)
      integer(MPIArg):: win_shm
      HElement_t(dp), pointer :: p_shm(:)

      deallocate(p_shm)
    end subroutine shared_deallocate_mpi

    subroutine shared_sync_mpi(win_shm)
      integer(MPIArg):: win_shm
    end subroutine shared_sync_mpi
#endif

    subroutine MPIBCast_inter_byte(p_shm,nbytes)
      HElement_t(dp):: p_shm
      integer:: nbytes
      integer(MPIArg):: ierr

      !only task 0 of each shared memory task range does the MPI communication
      if (iProcIndex_intra.eq.0) then
         call mpi_bcast(p_shm,int(nbytes,MPIArg), &
              MPI_BYTE,0_MPIArg,mpi_comm_inter,ierr)
      end if

    end subroutine MPIBCast_inter_byte

end module shared_memory_mpi
