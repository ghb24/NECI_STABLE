#include "macros.h"

! fpp types
#:set data_types = [['real(dp)', 'real'], ['integer(int64)', 'int64'], ['integer','int32'], ['complex(dp)','cmplx']]
module shared_array
    use constants
    use shared_memory_mpi
    implicit none
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

    contains
        ! allocation and deallocation routines
        procedure :: shared_alloc => safe_shared_memory_alloc_${data_name}$
        procedure :: shared_dealloc => safe_shared_memory_dealloc_${data_name}$
    end type shared_array_${data_name}$_t
#:endfor
contains
    
  !------------------------------------------------------------------------------------------!
  ! Auxiliary functions to prevent code duplication
  !------------------------------------------------------------------------------------------!

#:for data_type, data_name in data_types    
  subroutine safe_shared_memory_alloc_${data_name}$(this, size)
    ! wrapper for shared_allocate_mpi that tests if the pointer is associated
    ! Input: win - MPI shared memory window for internal MPI usage
    !        ptr - pointer to be allocated, on return points to a shared memory segment of given size
    !        size - size of the memory segment to be allocated
      implicit none
      class(shared_array_${data_name}$_t) :: this
      integer(int64), intent(in) :: size

      ! if pointer was allocated prior, re-allocate the probabilities
      ! WARNING: DO NOT MANUALLY RE-ASSIGN ptr, THIS WILL MOST LIKELY BREAK STUFF
      call safe_shared_memory_dealloc_${data_name}$(this)
      call shared_allocate_mpi(this%win, this%ptr, (/size/))
  end subroutine safe_shared_memory_alloc_${data_name}$

  !------------------------------------------------------------------------------------------!

  subroutine safe_shared_memory_dealloc_${data_name}$(this)
    ! wrapper for shared_deallocate_mpi that tests if the pointer is associated
    ! Input: win - MPI shared memory window for internal MPI usage
    !        ptr - pointer to be deallocated (if associated)
    ! WARNING: THIS ASSUMES THAT IF ptr IS ASSOCIATED, IT POINTS TO AN MPI SHARED MEMORY
    !          WINDOW win
    implicit none
    class(shared_array_${data_name}$_t) :: this
    
    ! assume that if ptr is associated, it points to mpi shared memory
    if(associated(this%ptr)) call shared_deallocate_mpi(this%win, this%ptr)
    this%ptr => null()
  end subroutine safe_shared_memory_dealloc_${data_name}$
#:endfor

    
end module shared_array
