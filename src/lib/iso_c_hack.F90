
module iso_c_hack

#ifndef __ISO_C_HACK
    use, intrinsic :: iso_c_binding
#else

    use constants
    implicit none

    integer, parameter :: c_int32_t = int32
    integer, parameter :: c_int = c_int32_t
    integer, parameter :: c_bool = c_int
    integer, parameter :: c_char = 1
    integer, parameter :: c_double = dp
    character(1), parameter :: c_null_char = '\0'

#ifdef POINTER8
    integer, parameter :: c_size_t = int64
    type c_ptr
        integer(int64) :: p
    end type
    integer(int64), parameter :: C_NULL_PTR = 0
#else
    integer, parameter :: c_size_t = int32
    type c_ptr
        integer(int32) :: p
    end type
    integer(int32), parameter :: C_NULL_PTR = 0
#endif

#endif

end module

#ifdef __ISO_C_HACK
subroutine c_f_pointer (a, b, dims)

    integer, intent(in) :: dims(:)
    integer, dimension(:), pointer :: a
    integer, dimension(:), pointer :: b

    b => a

end subroutine
#endif

