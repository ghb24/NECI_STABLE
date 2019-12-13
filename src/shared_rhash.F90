#include "macros.h"

module shared_rhash
    ! The linked list hashtable implementation of hash.F90 is not suited for intra-node shared memory
    ! A read-only hashtable can be much more efficient when implemented with a fixed size
    ! This makes hashed integrals usable in shared memory and increases their efficiency
    use constants
    use shared_array
    use ParallelHelper, only: iProcIndex_intra
    implicit none

    private
    public :: shared_rhash_t
    ! We always hash integer values (typically indices)
    ! -> the array we address can have any type
    type :: shared_rhash_t
        ! All arrays here are mpi-3 shared memory
        private
        ! The hashed indices are stored in one contiguous array
        type(shared_array_int64_t) :: indices

        ! Alongside, we store the offset of the hash values
        type(shared_array_int64_t) :: hval_offsets
        ! The range of the hash function
        integer(int64) :: hval_range

        ! auxiliary 
        integer, allocatable :: mult(:)

    contains
        procedure :: hash_func
        ! Allocate the memory
        procedure :: alloc
        procedure :: dealloc
        ! Fill up the indices - since this is memory critical, we allow direct write to them
        procedure :: count_index
        procedure :: add_index
        ! Set up the table. It is read-only, so this is the only way to set it up
        procedure :: finalize_setup
        ! After counting the indices, we have to get the offsets (doing so on the fly is
        ! horiffically slow)
        procedure :: setup_offsets

        ! Look up an index. Returns the position in the contiguous array
        procedure :: lookup
    end type shared_rhash_t

contains

    function hash_func(this, index) result(hval)
        class(shared_rhash_t), intent(in) :: this 
        integer(int64), intent(in) :: index
        integer(int64) :: hval

! TODO: Implement an actual hash function
        hval = mod(int(index)-1,this%hval_range)+1
    end function hash_func

    !------------------------------------------------------------------------------------------!
    ! Memory management
    !------------------------------------------------------------------------------------------!   

    subroutine alloc(this, n_elem, htsize)
        class(shared_rhash_t), intent(inout) :: this
        integer(int64), intent(in) :: n_elem
        integer(int64), intent(in) :: htsize
        
        this%hval_range = htsize
        ! Store all the indices with nonzero entries
        call this%indices%shared_alloc(n_elem)
        ! For each possible hash value, there will be on offset
        ! Add one additional offset at the end for easier initialization
        call this%hval_offsets%shared_alloc(htsize+1)        
        ! Only on node-root, the multiplicity of each hash value is counted during setup
        if(iProcIndex_intra == 0) then
            allocate(this%mult(this%hval_range))
            ! Initialize with 0
            this%hval_offsets%ptr = 0
            this%mult = 1
        end if
    end subroutine alloc

    !------------------------------------------------------------------------------------------!

    subroutine dealloc(this)
        class(shared_rhash_t), intent(inout) :: this

        call this%indices%shared_dealloc()
        call this%hval_offsets%shared_dealloc()
        if(allocated(this%mult)) deallocate(this%mult)
    end subroutine dealloc

    !------------------------------------------------------------------------------------------!
    ! Initialisation routines
    !------------------------------------------------------------------------------------------!  

    subroutine count_index(this, index)
        class(shared_rhash_t), intent(inout) :: this
        integer(int64), intent(in) :: index

        integer(int64) :: n_hval

        ! all following entries get their offset increased by one
        if(iProcIndex_intra == 0) then
            n_hval = this%hash_func(index) + 1
            this%hval_offsets%ptr(n_hval) = this%hval_offsets%ptr(n_hval) + 1
        endif
    end subroutine count_index

    !------------------------------------------------------------------------------------------!

    subroutine setup_offsets(this)
        class(shared_rhash_t), intent(inout) :: this

        integer(int64) :: i

        ! The first offset stays unchanged
        if(iProcIndex_intra == 0) then
            do i = 2, this%hval_range+1
                this%hval_offsets%ptr(i) = this%hval_offsets%ptr(i) + this%hval_offsets%ptr(i-1)
            end do
        endif
    end subroutine setup_offsets

    !------------------------------------------------------------------------------------------!

    subroutine add_index(this, index, pos)
        class(shared_rhash_t), intent(inout) :: this
        integer(int64), intent(in) :: index
        integer(int64), intent(out) :: pos

        integer(int64) :: hval

        if(iProcIndex_intra == 0) then
            hval = this%hash_func(index)
            pos = this%hval_offsets%ptr(hval)+this%mult(hval)
            this%indices%ptr(pos) = index
            this%mult(hval) = this%mult(hval) + 1
        end if
    end subroutine add_index

    !------------------------------------------------------------------------------------------!

    subroutine finalize_setup(this)
        class(shared_rhash_t), intent(inout) :: this

        if(iProcIndex_intra == 0) then            
            deallocate(this%mult)
        endif
    end subroutine finalize_setup    

    !------------------------------------------------------------------------------------------!
    ! Read access
    !------------------------------------------------------------------------------------------!
    
    subroutine lookup(this, index, pos, t_found)
        class(shared_rhash_t), intent(in) :: this
        integer(int64), intent(in) :: index
        integer(int64), intent(out) :: pos
        logical, intent(out) :: t_found

        integer(int64) :: hval, lower, upper, i
        
        hval = this%hash_func(index)
        lower = this%hval_offsets%ptr(hval) + 1
        upper = this%hval_offsets%ptr(hval+1)

        t_found = .false.
        pos = 0
        do i = lower, upper
            if(this%indices%ptr(i) == index) then
                pos = i
                t_found = .true.
                return
            end if
        end do
    end subroutine lookup
    
end module shared_rhash
