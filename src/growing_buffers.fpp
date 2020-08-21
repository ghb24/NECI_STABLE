#:include "algorithms.fpph"

! fpp types
#:set DATA_TYPES = {('real', 'dp'): 'real', ('HElement_t', 'dp'): 'hel', &
                    &('integer', ''): 'int', ('integer', 'int32'): 'int32', ('integer', 'int64'): 'int64', }

#:set RANKS = [1, 2]

#:def _get_name(data_name, rank)
$:  '{}_{}D'.format(data_name, rank)
#:enddef

#:def shape_like_except_along(rank, along, var, size_along)
$:  ', '.join(size_along if i == int(along) else 'size({var}, {i})'.format(var=var, i =i) for i in range(1, int(rank) + 1))
#:enddef


module growing_buffers
    use constants
    implicit none

    private
#:for rank in RANKS
#:for type_spec, data_name in DATA_TYPES.items()
#:set type, kind = type_spec
    public :: buffer_${_get_name(data_name, rank)}$_t

    !> Re-sizeable array type that can be filled elementwise to build up a contiguous data chunk
    !! that can then be dumped to an allocatable
    type :: buffer_${_get_name(data_name, rank)}$_t
        private

        @{get_decl(${type}$, ${kind}$, ${rank}$, allocatable=True)}@ :: buf
        ! Internal position of the buffer
        integer(int64) :: pos
        ! Blocksize
        integer(int64) :: block_size
        ! Buffer size
        integer(int64) :: buf_size
    contains

        procedure :: init => init_${_get_name(data_name, rank)}$
        procedure :: finalize => finalize_${_get_name(data_name, rank)}$

        procedure :: add_val => add_val_${_get_name(data_name, rank)}$
        procedure :: dump => dump_${_get_name(data_name, rank)}$
        procedure :: num_elements => num_elements_${_get_name(data_name, rank)}$
    end type buffer_${_get_name(data_name, rank)}$_t
#:endfor
#:endfor
contains

#:for type_spec, data_name in DATA_TYPES.items()
#:set type, kind = type_spec

    !> Set up the re-sizeable array (buffer) with a fixed blocksize
    !> @param[in] block_size_  size of the blocks in which the buffer is allocated. Also the
    !!                         initial size
#:set rank = 1
    subroutine init_${_get_name(data_name, rank)}$ (this, block_size_)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this
        integer(int64) :: block_size_

        ! Set the initial buffersize to be one block
        this%block_size = block_size_
        this%buf_size = this%block_size
        ! allocate one block
        allocate(this%buf(this%buf_size))
        this%pos = 0

    end subroutine init_${_get_name(data_name, rank)}$

#:set rank = 2
    !> Set up the re-sizeable array (buffer) with a fixed blocksize
    !> @param[in] block_size_  size of the blocks in which the buffer is allocated. Also the
    !!                         initial size
    subroutine init_${_get_name(data_name, rank)}$ (this, rows, block_size_)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this
        integer, intent(in) :: rows
        integer(int64), intent(in) :: block_size_

        ! Set the initial buffersize to be one block
        this%block_size = block_size_
        this%buf_size = this%block_size
        ! allocate one block
        allocate(this%buf(rows, this%buf_size))
        this%pos = 0

    end subroutine init_${_get_name(data_name, rank)}$

    !------------------------------------------------------------------------------------------!

#:for rank in RANKS
    #:set select = functools.partial(_select, rank, rank)
    !> Deallocate the resource. This is automatically called when dumping the buffer, with a
    !! following re-initialization
    subroutine finalize_${_get_name(data_name, rank)}$ (this)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this

        if (allocated(this%buf)) deallocate(this%buf)
    end subroutine finalize_${_get_name(data_name, rank)}$

    !------------------------------------------------------------------------------------------!

    !> Append a value to the buffer, expanding it by a block when necessary
    !> @param[in] val  Value to be added
    subroutine add_val_${_get_name(data_name, rank)}$ (this, val)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this
        @{get_decl(${type}$, ${kind}$, ${int(rank) - 1}$)}@, intent(in) :: val

        ! If the buffer still has room, add the entry
        if (this%pos < this%buf_size) then
            this%pos = this%pos + 1
            @{select(this%buf, this%pos)}@ = val
        else
            ! else, expand the buffer by another block
            block
                @{get_decl(${type}$, ${kind}$, ${rank}$, allocatable=True)}@ :: tmp
                integer(int64) :: new_buf_size
                ! Fortran 2003 automatic allocation/assignment
                tmp = this%buf

                deallocate(this%buf)
                new_buf_size = this%buf_size + this%block_size
                allocate(this%buf(@{shape_like_except_along(${rank}$, ${rank}$, tmp, new_buf_size)}@))

                @{select(this%buf, 1 : this%buf_size)}@ = @{select(tmp, 1 : this%buf_size)}@

                this%buf_size = new_buf_size

                this%pos = this%pos + 1
                @{select(this%buf, this%pos)}@ = val
            end block
        end if
    end subroutine add_val_${_get_name(data_name, rank)}$

    !------------------------------------------------------------------------------------------!

    !> Returns the number of already stored elements in the buffer
    !> @return n_els  number of elements already added to the buffer
    function num_elements_${_get_name(data_name, rank)}$ (this) result(n_els)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(in) :: this
        integer(int64) :: n_els

        n_els = this%pos
    end function num_elements_${_get_name(data_name, rank)}$

    !------------------------------------------------------------------------------------------!

    !> Dump the buffer to an allocatable array, deleting it afterwards
    !> @param[out] tgt  Allocatable array (reset upon entry), contains the stored elements of
    !!                  the buffer on return. The buffer has to be reinitialized if used again.
    subroutine dump_${_get_name(data_name, rank)}$ (this, tgt)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this
        @{get_decl(${type}$, ${kind}$, ${rank}$, allocatable=True)}@, intent(out) :: tgt

        integer(int64) :: n_els

        ! Transfer the buffer's content to an external allocatable
        n_els = this%num_elements()
        allocate(tgt(@{shape_like_except_along(${rank}$, ${rank}$, this%buf, n_els)}@))

        @{select(tgt, : n_els)}@ = @{select(this%buf, : n_els)}@

        call this%finalize()
    end subroutine dump_${_get_name(data_name, rank)}$
#:endfor
#:endfor

end module growing_buffers
