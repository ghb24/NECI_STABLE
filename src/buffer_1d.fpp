! fpp types
#:set data_types = [['real(dp)', 'real'], ['integer(int64)', 'int64'], ['integer','int32'], ['HElement_t(dp)', 'hel']]

module buffer_1d
    use constants
    implicit none

    private
#:for data_type, data_name in data_types    
    public :: buffer_${data_name}$_t

    type buffer_${data_name}$_t
        private

        ${data_type}$, allocatable :: buf(:)
        ! Internal position of the buffer
        integer(int64) :: pos
        ! Blocksize
        integer(int64) :: block_size
        ! Buffer size
        integer(int64) :: buf_size
    contains

        procedure :: init => init_${data_name}$
        procedure :: finalize => finalize_${data_name}$

        procedure :: add_val => add_val_${data_name}$
        procedure :: dump => dump_${data_name}$
        procedure :: num_elements => num_elements_${data_name}$
    end type buffer_${data_name}$_t
#:endfor
contains

#:for data_type, data_name in data_types    
    subroutine init_${data_name}$(this, block_size_)
        class(buffer_${data_name}$_t), intent(inout) :: this
        integer(int64) :: block_size_

        ! Set the initial buffersize to be one block
        this%block_size = block_size_
        this%buf_size = this%block_size
        ! allocate one block
        allocate(this%buf(this%buf_size))
        this%pos = 0
        
    end subroutine init_${data_name}$

    !------------------------------------------------------------------------------------------!

    subroutine finalize_${data_name}$(this)
        class(buffer_${data_name}$_t), intent(inout) :: this

        if(allocated(this%buf)) deallocate(this%buf)
    end subroutine finalize_${data_name}$

    !------------------------------------------------------------------------------------------!    

    subroutine add_val_${data_name}$(this, val)
        class(buffer_${data_name}$_t), intent(inout) :: this
        ${data_type}$, intent(in) :: val

        ! If the buffer still has room, add the entry
        if(this%pos < this%buf_size) then
            this%pos = this%pos + 1            
            this%buf(this%pos) = val
        else
            ! else, expand the buffer by another block
            block
              ${data_type}$, allocatable :: tmp(:)
              integer(int64) :: new_buf_size
              ! Fortran 2003 automatic allocation/assignment
              tmp = this%buf
              
              deallocate(this%buf)
              new_buf_size = this%buf_size + this%block_size              
              allocate(this%buf(new_buf_size))
              
              this%buf(1:this%buf_size) = tmp(1:this%buf_size)
              
              this%buf_size = new_buf_size
              
              this%pos = this%pos + 1
              this%buf(this%pos) = val
            end block
        endif
    end subroutine add_val_${data_name}$
    
    !------------------------------------------------------------------------------------------!

    function num_elements_${data_name}$(this) result(n_els)
        class(buffer_${data_name}$_t), intent(in) :: this
        integer(int64) :: n_els

        n_els = this%pos
    end function num_elements_${data_name}$

    !------------------------------------------------------------------------------------------!

    !> Dump the buffer to an allocatable array, deleting it afterwards
    subroutine dump_${data_name}$(this, tgt)
        class(buffer_${data_name}$_t), intent(inout) :: this
        ${data_type}$, intent(out), allocatable :: tgt(:)

        integer(int64) :: n_els

        ! Transfer the buffer's content to an external allocatable
        n_els = this%num_elements()
        allocate(tgt(n_els))
        tgt(1:n_els) = this%buf(1:n_els)

        ! Reset the buffer
        call this%finalize()
        call this%init(this%block_size)
    end subroutine dump_${data_name}$
#:endfor
    
end module buffer_1d

