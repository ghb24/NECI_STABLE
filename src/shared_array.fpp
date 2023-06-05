#include "macros.h"

! fpp types
#:set data_types = [['real(dp)', 'real'], ['integer(int64)', 'int64'], ['integer(int32)','int32'], ['complex(dp)','cmplx'], ['logical','bool']]
module shared_array
    use constants, only: int32, int64, dp, MPIArg
    use shared_memory_mpi, only: shared_allocate_mpi, shared_deallocate_mpi
    use mpi, only: MPI_Barrier, MPI_Win_Sync
    use MPI_wrapper, only: mpi_comm_intra
    use MemoryManager, only: LogMemALloc, LogMemDealloc, TagIntType
    better_implicit_none
    private
    ! Shared memory array types are defined here
#:for data_type, data_name in data_types
    public :: shared_array_${data_name}$_t
    type :: shared_array_${data_name}$_t
        ! They contain a ptr (access to the array)
        ! WARNING: DO NOT MANUALLY RE-ASSIGN ptr, THIS WILL MOST LIKELY BREAK STUFF
        ${data_type}$, pointer :: ptr(:) => null()
        ! and an MPI window
        integer(MPIArg) :: win
        ! Tag for the NECI memory manager
        integer(TagIntType) :: tag = 0
    contains
        ! allocation and deallocation routines
        procedure :: shared_alloc => safe_shared_memory_alloc_${data_name}$
        procedure :: shared_dealloc => safe_shared_memory_dealloc_${data_name}$
        procedure :: sync => sync_${data_name}$
    end type shared_array_${data_name}$_t
#:endfor
contains

    !------------------------------------------------------------------------------------------!
    ! Auxiliary functions to prevent code duplication
    !------------------------------------------------------------------------------------------!

#:for data_type, data_name in data_types
    !> Wrapper for shared_allocate_mpi that tests if the pointer is associated
    !> @param[out] win  MPI shared memory window for internal MPI usage
    !> @param[out] ptr  pointer to be allocated, on return points to a shared memory segment of given size
    !> @param[in] size  size of the memory segment to be allocated
    subroutine safe_shared_memory_alloc_${data_name}$ (this, size, name)

        class(shared_array_${data_name}$_t) :: this
        integer(int64), intent(in) :: size
        character(*), intent(in), optional :: name
        character(*), parameter :: t_r = "shared_alloc"

        ! if pointer was allocated prior, re-allocate the probabilities
        ! WARNING: DO NOT MANUALLY RE-ASSIGN ptr, THIS WILL MOST LIKELY BREAK STUFF
        call safe_shared_memory_dealloc_${data_name}$ (this)
        call shared_allocate_mpi(this%win, this%ptr, (/size/))

        ! If a name is given, log the allocation
        if (associated(this%ptr) .and. present(name)) &
            call LogMemAlloc(name, size, sizeof(this%ptr(1)), t_r, this%tag)
    end subroutine safe_shared_memory_alloc_${data_name}$

    !------------------------------------------------------------------------------------------!

    !> wrapper for shared_deallocate_mpi that tests if the pointer is associated
    !> @param[in,out] win  MPI shared memory window for internal MPI usage
    !> @param[in,out] ptr  pointer to be deallocated (if associated)
    ! WARNING: THIS ASSUMES THAT IF ptr IS ASSOCIATED, IT POINTS TO AN MPI SHARED MEMORY
    !          WINDOW win
    subroutine safe_shared_memory_dealloc_${data_name}$ (this)
        class(shared_array_${data_name}$_t) :: this
        character(*), parameter :: t_r = "shared_dealloc"

        ! assume that if ptr is associated, it points to mpi shared memory
        if (associated(this%ptr)) call shared_deallocate_mpi(this%win, this%ptr)
        ! First, check if we have to log the deallocation
        if (this%tag /= 0) call LogMemDealloc(t_r, this%tag)
    end subroutine safe_shared_memory_dealloc_${data_name}$

    !> callls MPI_Win_Sync on the array's shared memory window to sync rma
    !! This has to be called between read/write epochs to ensure all tasks of a node are
    !! looking at the same shared data
    subroutine sync_${data_name}$ (this)
        class(shared_array_${data_name}$_t) :: this
        integer(MPIArg) :: ierr

        call MPI_Win_Sync(this%win, ierr)
        call MPI_Barrier(mpi_comm_intra, ierr)

    end subroutine sync_${data_name}$
#:endfor

end module shared_array
