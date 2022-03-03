#:include "macros.fpph"
#:include "algorithms.fpph"

! fpp types
#:set DATA_TYPES = {('real', 'dp'): 'real', ('HElement_t', 'dp'): 'hel', &
                    &('integer', ''): 'int', ('integer', 'int32'): 'int32', ('integer', 'int64'): 'int64', &
                    &('integer', 'int128'): 'int128', ('Token_t', ''): 'token'}

#:set RANKS = [1, 2]

#:def _get_name(data_name, rank)
$:  '{}_{}D'.format(data_name, rank)
#:enddef

#:def shape_like_except_along(rank, along, var, size_along)
$:  ', '.join(size_along if i == int(along) else 'size({var}, {i})'.format(var=var, i =i) for i in range(1, int(rank) + 1))
#:enddef


module growing_buffers
    use constants
    use fortran_strings, only: Token_t
    implicit none

    private
#:for rank in RANKS
#:for type_spec, data_name in DATA_TYPES.items()
#:set type, kind = type_spec
    public :: buffer_${_get_name(data_name, rank)}$_t

    !>  @brief
    !>  Re-sizeable array type that can be filled elementwise to build up a contiguous data chunk
    !>  that can then be dumped to an allocatable
    !>
    !>  @details
    !>  For multidimensional buffers only the last dimension can grow.
    !>  (i.e. it is only possible to add columns.)
    !>
    !>  The buffer has to be initiliazed before first use.
    !>  After dumping (`dump_reset`) it is automatically resetted.
    type :: buffer_${_get_name(data_name, rank)}$_t
        private

        @{get_decl(${type}$, ${kind}$, ${rank}$, allocatable=True)}@ :: buf
        ! Internal position of the buffer
        integer(int64) :: pos
        real(dp) :: grow_factor = 1.5_dp
        integer(int64) :: start_size = 100_int64

    contains

        procedure :: init => init_${_get_name(data_name, rank)}$
        procedure :: finalize => finalize_${_get_name(data_name, rank)}$
        procedure, private :: reset => reset_${_get_name(data_name, rank)}$

        procedure :: size => get_size_${_get_name(data_name, rank)}$
        procedure :: capacity => get_capacity_${_get_name(data_name, rank)}$

        procedure :: push_back => add_val_${_get_name(data_name, rank)}$

        procedure, private :: dump => dump_${_get_name(data_name, rank)}$
        procedure :: dump_reset => dump_reset_${_get_name(data_name, rank)}$
    end type buffer_${_get_name(data_name, rank)}$_t
#:endfor
#:endfor
contains

#:for type_spec, data_name in DATA_TYPES.items()
#:set type, kind = type_spec

#:set rank = 1
    !>  @brief
    !>  Set up the re-sizeable array (buffer) with a given start size and grow_factor.
    !>
    !>  @details
    !>  Has to be called before first use and can be called any time.
    !>
    !> @param[in] start_size Initial size of the buffer.
    !> @param[in] grow_factor Factor about which to grow the buffer, if the capacity is not sufficient.
    pure subroutine init_${_get_name(data_name, rank)}$ (this, grow_factor, start_size)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this
        real(dp), optional, intent(in) :: grow_factor
        integer(int64), optional, intent(in) :: start_size
        character(*), parameter :: this_routine = 'buffer::init'

        if (present(grow_factor)) this%grow_factor = grow_factor
        if (present(start_size)) this%start_size = start_size
        @:pure_ASSERT(this%grow_factor > 1.0_dp)
        @:pure_ASSERT(this%start_size >= 0_int64)

        if (.not. allocated(this%buf)) allocate(this%buf(this%start_size))
        this%pos = 0_int64
    end subroutine

    !>  @brief
    !>  Reset an already initiliazed buffer.
    pure subroutine reset_${_get_name(data_name, rank)}$ (this)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this

        deallocate(this%buf)
        allocate(this%buf(this%start_size))
        this%pos = 0_int64
    end subroutine

#:set rank = 2
    !>  @brief
    !>  Set up the re-sizeable array (buffer) with a given start size and grow_factor.
    !>
    !>  @details
    !>  Has to be called before first use and can be called any time.
    !>
    !> @param[in] rows Number of rows in the first dimension.
    !> @param[in] start_size Initial size of the buffer along the last dimension.
    !> @param[in] grow_factor Factor about which to grow the buffer along the last dimension,
    !>              if the capacity is not sufficient.
    pure subroutine init_${_get_name(data_name, rank)}$ (this, rows, grow_factor, start_size)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this
        integer, intent(in) :: rows
        real, optional, intent(in) :: grow_factor
        integer, optional, intent(in) :: start_size
        character(*), parameter :: this_routine = 'buffer::init'

        if (present(grow_factor)) this%grow_factor = grow_factor
        if (present(start_size)) this%start_size = start_size
        @:pure_ASSERT(this%grow_factor > 1.0_dp)
        @:pure_ASSERT(this%start_size >= 0_int64)

        allocate(this%buf(rows, this%start_size))
        this%pos = 0_int64
    end subroutine init_${_get_name(data_name, rank)}$

    !>  @brief
    !>  Reset an already initiliazed buffer.
    pure subroutine reset_${_get_name(data_name, rank)}$ (this)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this

        integer(int64) :: rows

        rows = size(this%buf, 1)
        deallocate(this%buf)
        allocate(this%buf(rows, this%start_size))
        this%pos = 0_int64
    end subroutine

    !------------------------------------------------------------------------------------------!

#:for rank in RANKS
    #:set select = functools.partial(_select, rank, rank)
    !>  @brief
    !>  Deallocate the resource.
    pure subroutine finalize_${_get_name(data_name, rank)}$ (this)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this

        if (allocated(this%buf)) deallocate(this%buf)
    end subroutine finalize_${_get_name(data_name, rank)}$

    !------------------------------------------------------------------------------------------!

    !>  @brief
    !>  Returns the number of already stored elements in the buffer along the last dimension.
    !>
    !>  @return n_els Number of elements already added to the buffer.
    pure function get_size_${_get_name(data_name, rank)}$ (this) result(n_els)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(in) :: this
        integer(int64) :: n_els

        n_els = this%pos
    end function

    !------------------------------------------------------------------------------------------!

    !>  @brief
    !>  Returns the capacity of the buffer along the last dimension.
    !>
    !>  @return n_els Number of elements already added to the buffer.
    pure function get_capacity_${_get_name(data_name, rank)}$ (this) result(capacity)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(in) :: this
        integer(int64) :: capacity

        capacity = size(this%buf, ${rank}$)
    end function

    !------------------------------------------------------------------------------------------!

    !>  @brief
    !>  Append a value to the buffer, expanding the capacity if necessary.
    !>
    !>  @param[in] val Value to be added
    pure subroutine add_val_${_get_name(data_name, rank)}$ (this, val)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this
        @{get_decl(${type}$, ${kind}$, ${int(rank) - 1}$)}@, intent(in) :: val


        ! If the buffer still has room, add the entry
        if (this%pos < size(this%buf, ${rank}$)) then
            this%pos = this%pos + 1_int64
            @{select(this%buf, this%pos)}@ = val
        else
            ! else, expand the buffer by another block
            block

                @{get_decl(${type}$, ${kind}$, ${rank}$, allocatable=True)}@ :: tmp
                integer(int64) :: new_buf_size
                ! Fortran 2003 automatic allocation/assignment
                tmp = this%buf

                deallocate(this%buf)
                ! We add a constant offset to allow growth if start_size == 0.
                ! The grow_factor then takes over for larger numbers and prevents the O(n^2) scaling.
                new_buf_size = ceiling(real(size(this%buf, ${rank}$), kind=dp) * this%grow_factor, kind=int64) + 10_int64
                allocate(this%buf(@{shape_like_except_along(${rank}$, ${rank}$, tmp, new_buf_size)}@))

                @{select(this%buf, : size(tmp, ${rank}$))}@ = tmp

                this%pos = this%pos + 1_int64
                @{select(this%buf, this%pos)}@ = val
            end block
        end if
    end subroutine add_val_${_get_name(data_name, rank)}$


    !------------------------------------------------------------------------------------------!

    !>  @brief
    !>  Dump the buffer to an allocatable array.
    !>
    !>  @param[out] tgt Allocatable array (reset upon entry), contains the stored elements of
    !>                   the buffer on return. The buffer has to be reinitialized if used again.
    pure subroutine dump_${_get_name(data_name, rank)}$ (this, tgt)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this
        @{get_decl(${type}$, ${kind}$, ${rank}$, allocatable=True)}@, intent(out) :: tgt

        tgt = @{select(this%buf, : this%size())}@
    end subroutine

    !>  Dump the buffer to an allocatable array and reset the buffer.
    !>
    !>  @param[out] tgt Allocatable array (reset upon entry), contains the stored elements of
    !>                   the buffer on return. The buffer is writable afterwards.
    pure subroutine dump_reset_${_get_name(data_name, rank)}$ (this, tgt)
        class(buffer_${_get_name(data_name, rank)}$_t), intent(inout) :: this
        @{get_decl(${type}$, ${kind}$, ${rank}$, allocatable=True)}@, intent(out) :: tgt

        call this%dump(tgt)
        call this%reset()
    end subroutine
#:endfor
#:endfor

end module growing_buffers
