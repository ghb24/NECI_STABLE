module shared_ragged_array
    use constants
    use shared_array
    implicit none

    private
    public :: shared_ragged_array_t

    type :: auxiliary_t
        integer(int32), pointer :: res(:)
    end type auxiliary_t

    type :: shared_ragged_array_t
        private

        type(shared_array_int32_t) :: data_array
        type(auxiliary_t), allocatable :: ptr(:)
        integer, allocatable :: store_sizes(:)
    contains
        procedure :: shared_alloc
        procedure :: shared_dealloc
        procedure :: reassign_pointers
        procedure :: pos_1d
        procedure :: pos_2d
        generic :: sub => pos_1d, pos_2d
        procedure :: set_val
    end type shared_ragged_array_t

contains

    subroutine shared_alloc(this, sizes)
        class(shared_ragged_array_t), intent(inout) :: this
        integer, intent(in) :: sizes(:)

        integer :: n_entries, i

        ! Allocate the shared resource
        call this%data_array%shared_alloc(int(sum(sizes),int64))

        ! Assign the pointers
        n_entries = size(sizes)
        allocate(this%ptr(n_entries))

        ! Keep a local copy of sizes (fortran 2003 automatic allocation)
        this%store_sizes = sizes

        ! Set the internal pointers
        call this%reassign_pointers()
    end subroutine shared_alloc

    !------------------------------------------------------------------------------------------!    
    
    subroutine shared_dealloc(this)
        class(shared_ragged_array_t), intent(inout) :: this

        call this%data_array%shared_dealloc()
        if(allocated(this%ptr)) deallocate(this%ptr)
        if(allocated(this%store_sizes)) deallocate(this%store_sizes)        
    end subroutine shared_dealloc

    !------------------------------------------------------------------------------------------!    

    subroutine reassign_pointers(this)
        class(shared_ragged_array_t), intent(inout) :: this
        integer :: n_entries
        integer :: i, win_start, win_end
        
        n_entries = size(this%store_sizes)
        win_start = 1
        do i = 1, n_entries
            win_end = win_start - 1 + this%store_sizes(i)
            this%ptr(i)%res => this%data_array%ptr(win_start:win_end)
            win_start = win_end + 1
        end do
    end subroutine reassign_pointers

    !------------------------------------------------------------------------------------------!

    subroutine set_val(this, i, j, val)
        class(shared_ragged_array_t), intent(inout) :: this        
        integer, intent(in) :: i, j
        integer, intent(in) :: val

        this%ptr(i)%res(j) = val
    end subroutine set_val

    !------------------------------------------------------------------------------------------!        

    function pos_2d(this, i, j) result(val)
        class(shared_ragged_array_t), intent(inout) :: this        
        integer, intent(in) :: i, j
        integer :: val

        val = this%ptr(i)%res(j)
    end function pos_2d

    !------------------------------------------------------------------------------------------!    

    function pos_1d(this, i) result(pt)
        class(shared_ragged_array_t), intent(inout) :: this        
        integer, intent(in) :: i
        integer, pointer :: pt(:)

        pt => this%ptr(i)%res
    end function pos_1d

end module shared_ragged_array
