
! Two-dimensional contiguous array class with non-uniform second dimension stored in shared memory. It is implemented this way to
! yield maximal memory efficiency

! Templated for these types
#:set data_types = [['real(dp)', 'real'], ['integer(int64)', 'int64'], ['integer(int32)','int32'], ['complex(dp)','cmplx'], ['logical','bool']]
! And these integer kinds can be used for indexing and allocation
#:set index_types = ['int32','int64']
module shared_ragged_array
    use constants
    use shared_array
    implicit none

    private
    
#:for data_type, data_name in data_types
    public :: shared_ragged_array_${data_name}$_t

    ! This type only serves the purpose to create an addressing array for quasi-2D-access -> allows to directly get pointers to
    ! sub-arrays
    type :: auxiliary_${data_name}$_t
        ${data_type}$, pointer :: res(:)
    end type auxiliary_${data_name}$_t

    !> Shared memory 2-D array template with non-uniform 2nd dimension ("ragged") of type ${data_name}$
    type :: shared_ragged_array_${data_name}$_t
        private

        ! Actual shared resource
        type(shared_array_${data_name}$_t) :: data_array
        ! Indexing array to hold pointers to the sub-arrays
        type(auxiliary_${data_name}$_t), allocatable :: ptr(:)
        ! Sizes of the sub-arrays
        integer, allocatable :: store_sizes(:)
    contains

        ! These functions work with different integer kinds as index values (index_types)
#:for index_type in index_types

        ! Generic interfaces
        ! 
        generic :: shared_alloc => shared_alloc_${data_name}$_${index_type}$
        ! sub returns a pointer to a 1d subarray (sub1d) or a specific entry (sub2d)
        generic :: sub => pos_1d_${data_name}$_${index_type}$, pos_2d_${data_name}$_${index_type}$
        ! Set val sets a specific value of the array
        generic :: set_val => set_val_${data_name}$_${index_type}$

        ! Index kind-specific implementations
        procedure :: shared_alloc_${data_name}$_${index_type}$
        procedure :: pos_1d_${data_name}$_${index_type}$
        procedure :: pos_2d_${data_name}$_${index_type}$
        procedure :: set_val_${data_name}$_${index_type}$
#:endfor
        procedure :: shared_dealloc => shared_dealloc_${data_name}$
        ! Sync synchronizes the content of the array within the inter-node communicator
        procedure :: sync => sync_${data_name}$
        ! For setting up the subarray pointers
        procedure :: reassign_pointers => reassign_pointers_${data_name}$
    end type shared_ragged_array_${data_name}$_t
#:endfor

contains

#:for data_type, data_name in data_types
    #:for index_type in index_types
    subroutine shared_alloc_${data_name}$_${index_type}$(this, sizes)
        class(shared_ragged_array_${data_name}$_t), intent(inout) :: this
        integer(${index_type}$), intent(in) :: sizes(:)

        integer :: n_entries

        ! Allocate the shared resource
        call this%data_array%shared_alloc(int(sum(sizes), int64))

        ! Assign the pointers
        n_entries = size(sizes)
        allocate(this%ptr(n_entries))

        ! Keep a local copy of sizes (fortran 2003 automatic allocation)
        allocate(this%store_sizes(n_entries))
        this%store_sizes(1:n_entries) = sizes(1:n_entries)

        ! Set the internal pointers
        call this%reassign_pointers()
    end subroutine shared_alloc_${data_name}$_${index_type}$
    #:endfor

    !------------------------------------------------------------------------------------------!

    subroutine shared_dealloc_${data_name}$(this)
        class(shared_ragged_array_${data_name}$_t), intent(inout) :: this

        call this%data_array%shared_dealloc()
        if (allocated(this%ptr)) deallocate(this%ptr)
        if (allocated(this%store_sizes)) deallocate(this%store_sizes)
    end subroutine shared_dealloc_${data_name}$

    !------------------------------------------------------------------------------------------!

    subroutine reassign_pointers_${data_name}$(this)
        class(shared_ragged_array_${data_name}$_t), intent(inout) :: this
        integer :: n_entries
        integer :: i, win_start, win_end

        n_entries = size(this%store_sizes)
        win_start = 1
        do i = 1, n_entries
            win_end = win_start - 1 + this%store_sizes(i)
            this%ptr(i)%res => this%data_array%ptr(win_start:win_end)
            win_start = win_end + 1
        end do
    end subroutine reassign_pointers_${data_name}$

    !------------------------------------------------------------------------------------------!
#:for index_type in index_types
    subroutine set_val_${data_name}$_${index_type}$(this, i, j, val)
        class(shared_ragged_array_${data_name}$_t), intent(inout) :: this
        integer(${index_type}$), intent(in) :: i, j
        ${data_type}$, intent(in) :: val

        this%ptr(i)%res(j) = val
    end subroutine set_val_${data_name}$_${index_type}$
#:endfor
    !------------------------------------------------------------------------------------------!

    subroutine sync_${data_name}$(this)
        class(shared_ragged_array_${data_name}$_t), intent(inout) :: this

        call this%data_array%sync()
    end subroutine sync_${data_name}$

    !------------------------------------------------------------------------------------------!

#:for index_type in index_types
    function pos_2d_${data_name}$_${index_type}$(this, i, j) result(val)
        class(shared_ragged_array_${data_name}$_t), intent(inout) :: this
        integer(${index_type}$), intent(in) :: i, j
        ${data_type}$ :: val

        val = this%ptr(i)%res(j)
    end function pos_2d_${data_name}$_${index_type}$

    !------------------------------------------------------------------------------------------!

    function pos_1d_${data_name}$_${index_type}$(this, i) result(pt)
        class(shared_ragged_array_${data_name}$_t), intent(inout) :: this
        integer(${index_type}$), intent(in) :: i
        ${data_type}$, pointer :: pt(:)

        pt => this%ptr(i)%res
    end function pos_1d_${data_name}$_${index_type}$
#:endfor
#:endfor

end module shared_ragged_array
