module vasp_neci_interface

    use constants, only: dp

contains
    subroutine construct_ijab_one(i, j, a, b, integral)
        implicit none
        integer :: i, j, a, b
        complex(dp) :: integral(1, 1)
        ! Make warnings go away
        i = i; j = j; a = a; b = b; integral = integral
    end subroutine construct_ijab_one
end module vasp_neci_interface
