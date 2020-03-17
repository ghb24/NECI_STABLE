#include "macros.h"

module shared_rhash
    ! The linked list hashtable implementation of hash.F90 is not suited for intra-node shared memory
    ! A read-only hashtable can be much more efficient when implemented with a fixed size
    ! This makes hashed integrals usable in shared memory and increases their efficiency
    use constants
    use shared_array
    use ParallelHelper, only: iProcIndex_intra, mpi_comm_intra
    implicit none

    private
    public :: shared_rhash_t, initialise_shared_rht, shared_rht_lookup

    interface initialise_shared_rht
        module procedure initialise_shared_rht_impl
        module procedure initialise_shared_rht_expl
    end interface initialise_shared_rht

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

        ! Allocate the memory
        procedure :: alloc
        procedure :: dealloc
        ! Fill up the indices - since this is memory critical, we allow direct write to them
        procedure :: count_value
        procedure :: add_value
        ! Set up the table. It is read-only, so this is the only way to set it up
        procedure :: finalize_setup
        ! After counting the indices, we have to get the offsets (doing so on the fly is
        ! horiffically slow)
        procedure :: setup_offsets

        ! Look up an index. Returns the position in the contiguous array
        procedure :: direct_lookup
        procedure :: callback_lookup

        ! Tell if the conflicts have been counted
        procedure :: known_conflicts
        ! Check how large the hash table shall be
        procedure :: val_range
    end type shared_rhash_t

contains

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
        call this%hval_offsets%shared_alloc(this%hval_range+1)
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

    !> Log the occurence of this hash value in the set of values to be stored
    !! Does not add it, only updates the offsets
    !> @param[in] hval hash value to be logged
    subroutine count_value(this, hval)
        class(shared_rhash_t), intent(inout) :: this
        integer(int64), intent(in) :: hval

        integer(int64) :: n_hval

        ! all following entries get their offset increased by one
        if(iProcIndex_intra == 0) then
            n_hval = hval + 1
            this%hval_offsets%ptr(n_hval) = this%hval_offsets%ptr(n_hval) + 1
        endif
    end subroutine count_value

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
    !> @param[in] hval value to be stored
    !> @param[in] index index belonging to this value
    !> @param[out] pos on return, the position where this value was stored
    subroutine add_value(this, hval, index, pos)
        class(shared_rhash_t), intent(inout) :: this
        integer(int64), intent(in) :: hval, index
        integer(int64), intent(out) :: pos

        if(iProcIndex_intra == 0) then
            pos = this%hval_offsets%ptr(hval)+this%mult(hval)
            this%indices%ptr(pos) = index
            this%mult(hval) = this%mult(hval) + 1
        end if
    end subroutine add_value

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
    !> @param[in] hval hash value of the index to look up
    !> @param[in] index value to be looked up
    !> @param[out] pos on return, the position of index if found, else 0
    !> @param[out] t_found on return, true if and only if index was found
    subroutine direct_lookup(this, hval, index, pos, t_found)
        class(shared_rhash_t), intent(in) :: this
        integer(int64), intent(in) :: index, hval
        integer(int64), intent(out) :: pos
        logical, intent(out) :: t_found

        integer(int64) :: lower, upper, i
        
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
    end subroutine direct_lookup

    !------------------------------------------------------------------------------------------!

    !> Generic lookup routine, using an external routine for verification
    !! DOES NOT TO THE SAME AS direct_lookup
    !> @param[in] hval hash value of the index to look up
    !> @param[out] pos on return, the matching entry
    !> @param[out] t_found on return, true if and only if index was found
    !> @param[in] verify  function to check if an entry matches    
    subroutine callback_lookup(this, hval, pos, t_found, loc_verify)
        class(shared_rhash_t), intent(in) :: this
        integer(int64), intent(in) :: hval
        integer(int64), intent(out) :: pos
        logical, intent(out) :: t_found

        integer(int64) :: lower, upper, i

        interface
            function loc_verify(i) result(match)
                use constants
                integer(int64), intent(in) :: i
                logical :: match
            end function loc_verify
        end interface

        lower = this%hval_offsets%ptr(hval) + 1
        upper = this%hval_offsets%ptr(hval+1)

        t_found = .false.
        pos = 0
        do i = lower, upper
            if(loc_verify(this%indices%ptr(i))) then
                pos = this%indices%ptr(i)
                t_found = .true.
                return
            end if
        end do
    end subroutine callback_lookup

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

    !------------------------------------------------------------------------------------------!

    !> Get the range of hash table values of this ht
    !> @return h_range  maximum possible hash value of this ht
    function val_range(this) result(h_range)
        class(shared_rhash_t), intent(in) :: this
        integer(int64) :: h_range

        h_range = this%hval_range
    end function val_range

    !------------------------------------------------------------------------------------------!
    ! Non-member function for global utility
    !------------------------------------------------------------------------------------------!

    ! Default the determinant size
    subroutine initialise_shared_rht_impl(ilut_list, space_size, hash_table, ht_size)    
        use SystemData, only: nel
        integer(n_int), intent(in) :: ilut_list(0:,:)
        integer, intent(in) :: space_size
        type(shared_rhash_t), intent(out) :: hash_table
        integer, intent(in), optional :: ht_size
        integer :: ht_size_

        def_default(ht_size_, ht_size, space_size)        

        call initialise_shared_rht_expl(ilut_list, space_size, hash_table, nel, ht_size_)
    end subroutine initialise_shared_rht_impl

    !------------------------------------------------------------------------------------------!    

    subroutine initialise_shared_rht_expl(ilut_list, space_size, hash_table, det_size, ht_size)
        use bit_reps, only: decode_bit_det
        use hash, only: FindWalkerHash

        integer(n_int), intent(in) :: ilut_list(0:,:)
        integer, intent(in) :: space_size
        type(shared_rhash_t), intent(out) :: hash_table
        integer, intent(in) :: det_size
        ! ht_size cannot be defaulted anymore as this would be ambigious        
        integer, intent(in) :: ht_size

        integer :: nI(det_size)
        integer :: i, ierr
        integer(int64) :: hash_val, pos       
        
        call hash_table%alloc(int(space_size,int64), int(ht_size,int64))

        ! Count the number of states with each hash value.
        do i = 1, space_size
            call decode_bit_det(nI, ilut_list(:,i))
            hash_val = FindWalkerHash(nI, ht_size)
            call hash_table%count_value(hash_val)
        end do

        call hash_table%setup_offsets()
        ! Now fill in the indices of the states in the space.
        do i = 1, space_size
            call decode_bit_det(nI, ilut_list(:,i))
            hash_val = FindWalkerHash(nI, ht_size)
            call hash_table%add_value(hash_val, int(i,int64), pos)
        end do

        ! Synchronize the node afterwards to keep tasks from using the un-initialized ht
        call MPI_Barrier(mpi_comm_intra, ierr)
    end subroutine initialise_shared_rht_expl

    !------------------------------------------------------------------------------------------!    

    subroutine shared_rht_lookup(core_ht, ilut, nI, tgt_space, i, core_state)
        use hash, only: FindWalkerHash
        use FciMCData, only: determ_space_size_int
        use bit_rep_data, only: NIfTot, NIfDBO
        type(shared_rhash_t), intent(in) :: core_ht
        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer, intent(in) :: nI(:)
        integer(int64), intent(in) :: tgt_space(0:,1:)

        integer, intent(out) :: i
        logical, intent(out) :: core_state
        integer(int64) :: hash_val, i_tmp

        hash_val = FindWalkerHash(nI, int(core_ht%val_range()))
        
        call core_ht%callback_lookup(hash_val, i_tmp, core_state, loc_verify)
        ! cast down to int32
        i = int(i_tmp)
        
    contains

        function loc_verify(ind) result(match)
            integer(int64), intent(in) :: ind
            logical :: match

            match = all(ilut(0:NIfDBO) == tgt_space(0:NIfDBO,ind))
        end function loc_verify

    end subroutine shared_rht_lookup   
    
end module shared_rhash
