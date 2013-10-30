interface
    function fn (i, j, k, l, fn2) result(hel) bind(c)
        use constants, only: dp
        use iso_c_hack
        implicit none
        integer, intent(in) :: i, j, k, l
!#ifdef SX
!        type(c_ptr), intent(in) :: fn2
!#else
        type(c_ptr), intent(in), value :: fn2
!#endif
        HElement_t :: hel
    end function
end interface
