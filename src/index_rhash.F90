module index_rhash
    ! A wrapper for the shared read-only hashtable that supplies its own hash function
    ! and for which the input of the hash function and the stored indices are identical
    use shared_rhash, only: shared_rhash_t
    use constants
    private
    public :: index_rhash_t

    type :: index_rhash_t
        private

        type(shared_rhash_t) :: shared_ht

    contains
        ! The hash function used for storing indices
        procedure :: hash_function
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

        ! Synchronize the table between tasks
        procedure :: sync
    end type index_rhash_t

contains

    !> Get the hash value for an arbitrary input value
    !> @param index input value to get the hash value for
    function hash_function(this, index) result(hval)
        class(index_rhash_t), intent(in) :: this
        integer(int64), intent(in) :: index
        integer(int64) :: hval

! TODO: Implement an actual hash function
        hval = mod(index - 1, this%shared_ht%val_range()) + 1_int64
    end function hash_function

    !------------------------------------------------------------------------------------------!
    ! Memory management
    !------------------------------------------------------------------------------------------!

    !> Allocate the internal (shared) memory
    !> @param[in] n_elem number of distinct values to store
    !> @param[in] htsize range of the hash function
    subroutine alloc(this, n_elem, htsize)
        class(index_rhash_t), intent(inout) :: this
        integer(int64), intent(in) :: n_elem
        integer(int64), intent(in) :: htsize

        call this%shared_ht%alloc(n_elem, htsize)
    end subroutine alloc

    !------------------------------------------------------------------------------------------!

    !> Deallocate all arrays associated with this hash table object
    subroutine dealloc(this)
        class(index_rhash_t), intent(inout) :: this

        call this%shared_ht%dealloc()
    end subroutine dealloc

    !------------------------------------------------------------------------------------------!
    ! Initialisation routines
    !------------------------------------------------------------------------------------------!

    !> Log the occurence of this index in the set of indices to be stored
    !! Does not add it, only updates the offsets
    !> @param[in] index index value to be logged
    subroutine count_index(this, index)
        class(index_rhash_t), intent(inout) :: this
        integer(int64), intent(in) :: index

        integer(int64) :: hval

        hval = this%hash_function(index)
        call this%shared_ht%count_value(hval)
    end subroutine count_index

    !------------------------------------------------------------------------------------------!

    !> For performance reasons, we cannot directly calculate the offsets, but instead
    !! first count the number of conflicts per hash value. Then, we sum these up cumulatively
    !! Directly counting the offsets is horrifically slow
    subroutine setup_offsets(this)
        class(index_rhash_t), intent(inout) :: this

        call this%shared_ht%setup_offsets()
    end subroutine setup_offsets

    !------------------------------------------------------------------------------------------!

    !> Add an input value to the stored values, assuming we already know the offsets
    !> @param[in] index value to be stored
    !> @param[out] pos on return, the position where this value was stored
    subroutine add_index(this, index, pos)
        class(index_rhash_t), intent(inout) :: this
        integer(int64), intent(in) :: index
        integer(int64), intent(out) :: pos

        integer(int64) :: hval

        hval = this%hash_function(index)
        call this%shared_ht%add_value(hval, index, pos)
    end subroutine add_index

    !------------------------------------------------------------------------------------------!

    !> Dealloates temporary arrays used for initialisation
    subroutine finalize_setup(this)
        class(index_rhash_t) :: this

        call this%shared_ht%finalize_setup()
    end subroutine finalize_setup

    !------------------------------------------------------------------------------------------!
    ! Read access
    !------------------------------------------------------------------------------------------!

    !> Look up a value in this hash table. Returns whether the value is stored and if yes, where
    !> @param[in] index value to be looked up
    !> @param[out] pos on return, the position of index if found, else 0
    !> @param[out] t_found on return, true if and only if index was found
    subroutine lookup(this, index, pos, t_found)
        class(index_rhash_t), intent(in) :: this
        integer(int64), intent(in) :: index
        integer(int64), intent(out) :: pos
        logical, intent(out) :: t_found

        integer(int64) :: hval

        hval = this%hash_function(index)
        call this%shared_ht%direct_lookup(hval, index, pos, t_found)
    end subroutine lookup

    !------------------------------------------------------------------------------------------!

    !> During initialisation, we can only start writing values once the offsets are known.
    !! This requires knowledge about the number of conflicts per hash value. This function
    !! tells us whether the conflicts have already been counted.
    !> @return t_kc true if and only if the conflicts have already been counted.
    function known_conflicts(this) result(t_kc)
        class(index_rhash_t), intent(in) :: this
        logical :: t_kc

        t_kc = this%shared_ht%known_conflicts()
    end function known_conflicts

    !------------------------------------------------------------------------------------------!

    !> For a MPI-3 shared memory array, synchronization is required after/before each read/write epoch
    subroutine sync(this)
        class(index_rhash_t), intent(in) :: this

        call this%shared_ht%sync()
    end subroutine sync

end module index_rhash
