! This routine intentionally does not have a declared interface. This does
! leave it a little bit open to abuse, but that is what we want in this case.
!
! --> This is only a worker routine.
subroutine assign_proc (ptr, proc)

    ! In:   proc - A function pointer (effectively). n.b. passed by value.
    ! Out:  A type(c_ptr) which
    use, intrinsic :: iso_c_binding
    implicit none

    type(c_ptr), intent(in), value :: proc
    type(c_ptr), intent(inout) :: ptr

    ptr = proc
end subroutine

subroutine sub_dispatcher_1 (fn, arg1)
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine fn (arg)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: arg
        end subroutine
    end interface

    type(c_ptr), intent(in) :: arg1

    call fn(arg1)
end subroutine

subroutine sub_dispatcher_5 (fn, arg1, arg2, arg3, arg4, arg5)
    use, intrinsic :: iso_c_binding
    implicit none

    interface
        subroutine fn (arg1, arg2, arg3, arg4, arg5)
            use, intrinsic :: iso_c_binding
            implicit none
            type(c_ptr), intent(in), value :: arg1, arg2, arg3, arg4, arg5
        end subroutine
    end interface

    type(c_ptr), intent(in) :: arg1, arg2, arg3, arg4, arg5

    call fn(arg1, arg2, arg3, arg4, arg5)
end subroutine





