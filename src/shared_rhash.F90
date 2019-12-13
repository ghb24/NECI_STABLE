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

    !> The shared read-only hash table stores a given number of arbitrary input values
    !! in one contiguous array and addresses this contiguous array using a hashtable
    !! The input values are stored in order of ascending hash value, with conflicts stored
    !! adjacently. For each hash value, the position of the first value with that hash value
    !! is stored. The lookup then searches for a given value between the first and the last
    !! stored value with the same hash value. 
    type :: shared_rhash_t             
        ! All arrays here are mpi-3 shared memory
        private
        ! The hashed data is stored in one contiguous array
        ! Typically, we will be storing indices of some array
        type(shared_array_int64_t) :: indices

        ! Alongside, we store the offset of the hash values - the position of the first
        ! index for each hash value
        type(shared_array_int64_t) :: hval_offsets
        ! The range of the hash function
        integer(int64) :: hval_range

        ! auxiliary array for initialisation. This stores how many input values we already
        ! added for each hash value
        integer, allocatable :: mult(:)

        ! are the conflicts of the hashtable already counted (have the offsets been set?)
        logical :: t_conflicts_known = .false.

    contains
        ! The hash function used for storing indices
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

        ! Tell if the conflicts have been counted
        procedure :: known_conflicts
    end type shared_rhash_t

contains

    !> Get the hash value for an arbitrary input value
    !> @param index input value to get the hash value for
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

    !> Allocate the internal (shared) memory
    !> @param[in] n_elem number of distinct values to store
    !> @param[in] htsize range of the hash function
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

    !> Deallocate all arrays associated with this hash table object
    subroutine dealloc(this)
        class(shared_rhash_t), intent(inout) :: this

        call this%indices%shared_dealloc()
        call this%hval_offsets%shared_dealloc()
        if(allocated(this%mult)) deallocate(this%mult)
    end subroutine dealloc

    !------------------------------------------------------------------------------------------!
    ! Initialisation routines
    !------------------------------------------------------------------------------------------!  

    !> Log the occurence of this index in the set of indices to be stored
    !! Does not add it, only updates the offsets
    !> @param[in] index value to be logged
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

    !> For performance reasons, we cannot directly calculate the offsets, but instead
    !! first count the number of conflicts per hash value. Then, we sum these up cumulatively
    !! Directly counting the offsets is horrifically slow
    subroutine setup_offsets(this)
        class(shared_rhash_t), intent(inout) :: this

        integer(int64) :: i

        ! The first offset stays unchanged
        if(iProcIndex_intra == 0) then
            do i = 2, this%hval_range+1
                this%hval_offsets%ptr(i) = this%hval_offsets%ptr(i) + this%hval_offsets%ptr(i-1)
            end do
        endif
        this%t_conflicts_known = .true.
    end subroutine setup_offsets

    !------------------------------------------------------------------------------------------!

    !> Add an input value to the stored values, assuming we already know the offsets
    !> @param[in] index value to be stored
    !> @param[out] pos on return, the position where this value was stored
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

    !> Dealloates temporary arrays used for initialisation
    subroutine finalize_setup(this)
        class(shared_rhash_t), intent(inout) :: this

        if(iProcIndex_intra == 0) then            
            deallocate(this%mult)
        endif
    end subroutine finalize_setup    

    !------------------------------------------------------------------------------------------!
    ! Read access
    !------------------------------------------------------------------------------------------!

    !> Look up a value in this hash table. Returns whether the value is stored and if yes, where
    !> @param[in] index value to be looked up
    !> @param[out] pos on return, the position of index if found, else 0
    !> @param[out] t_found on return, true if and only if index was found
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

    !------------------------------------------------------------------------------------------!

    !> During initialisation, we can only start writing values once the offsets are known.
    !! This requires knowledge about the number of conflicts per hash value. This function
    !! tells us whether the conflicts have already been counted.
    !> @return t_kc true if and only if the conflicts have already been counted.
    function known_conflicts(this) result(t_kc)
        class(shared_rhash_t), intent(in) :: this
        logical :: t_kc

        t_kc = this%t_conflicts_known
    end function known_conflicts
    
end module shared_rhash
