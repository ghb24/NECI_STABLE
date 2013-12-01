! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
! This routine intentionally does not have a declared interface. This does
! leave it a little bit open to abuse, but that is what we want in this case.
!
! --> This is only a worker routine.
subroutine assign_proc (ptr, proc)

    ! In:   proc - A function pointer (effectively). n.b. passed by value.
    ! Out:  A type(c_ptr) which
    use iso_c_hack
    use constants, only: int32, int64
    implicit none

    ! sds: This is to work around a bug in Pathscale which causes the 
    !      compiler to crash on valid code.
#if defined __PATHSCALE__ || defined __PGI
#ifdef POINTER8
    integer(int64), intent(in), value :: proc
    integer(int64), intent(inout) :: ptr
#else
    integer(int32), intent(in), value :: proc
    integer(int32), intent(inout) :: ptr
#endif
#else
    type(c_ptr), intent(in), value :: proc
    type(c_ptr), intent(inout) :: ptr
#endif

    ptr = proc
end subroutine

! This routine intentionally does not have a declared interface. This does
! leave it a little bit open to abuse, but that is what we want in this case.
!
! --> This is only a worker routine.
! --> Duplicate of assign_proc to avoid over zealous type checking by the
!     NAG compiler.
subroutine assign_data (ptr, proc)

    ! In:   proc - A function pointer (effectively). n.b. passed by value.
    ! Out:  A type(c_ptr) which
    use iso_c_hack
    use constants, only: int32, int64
    implicit none

    ! sds: This is to work around a bug in Pathscale which causes the 
    !      compiler to crash on valid code.
#if defined __PATHSCALE__ || defined __PGI
#ifdef POINTER8
    integer(int64), intent(in), value :: proc
    integer(int64), intent(inout) :: ptr
#else
    integer(int32), intent(in), value :: proc
    integer(int32), intent(inout) :: ptr
#endif
#else
    type(c_ptr), intent(in), value :: proc
    type(c_ptr), intent(inout) :: ptr
#endif

    ptr = proc
end subroutine

subroutine sub_dispatcher_1 (fn, arg1)
    use iso_c_hack
    implicit none

    interface
        subroutine fn (arg)
            use iso_c_hack
            implicit none
            type(c_ptr), intent(in), value :: arg
        end subroutine
    end interface

    type(c_ptr), intent(in) :: arg1

    call fn(arg1)
end subroutine

subroutine sub_dispatcher_10 (fn, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
    use iso_c_hack
    implicit none

    interface
        subroutine fn (arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
            use iso_c_hack
            implicit none
            type(c_ptr), intent(in), value :: arg1, arg2, arg3, arg4, arg5
            type(c_ptr), intent(in), value :: arg6, arg7, arg8, arg9
        end subroutine
    end interface

    type(c_ptr), intent(in) :: arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9

    call fn(arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
end subroutine

subroutine print_pointer (ptr)
    use constants, only: int64, int32
    implicit none

#ifdef POINTER8
    integer(int64), intent(inout) :: ptr
#else
    integer(int32), intent(inout) :: ptr
#endif
    write (6, '("Pointer: ", i16)') ptr

end subroutine





