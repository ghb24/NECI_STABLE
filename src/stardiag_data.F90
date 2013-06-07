! Copyright (c) 2013, Ali Alavi unless otherwise noted.
! This program is integrated in Molpro with the permission of George Booth and Ali Alavi
 
module StarDiagData
    type star_walker
        ! Det indicates the determinant that the walker is currently at. 
        ! WSign: +ve = .true. , -ve = .false. Initially, dets are HF and 
        ! sign is positive
        integer :: det = 1
        logical :: wSign = .true.
    end type

    interface assignment(=)
        module procedure star_walker_assign
    end interface
    interface operator(.eq.)
        module procedure star_walker_eq
    end interface
    interface operator(.ne.)
        module procedure star_walker_ne
    end interface
    interface operator(.gt.)
        module procedure star_walker_gt
    end interface
    interface operator(.lt.)
        module procedure star_walker_lt
    end interface
contains
    elemental subroutine star_walker_assign (lhs, rhs)
        type(star_walker), intent(out) :: lhs
        type(star_walker), intent(in) :: rhs
        lhs%det = rhs%det
        lhs%wSign = rhs%wSign
    end subroutine
    ! n.b. These comparisons only compare the det term (Which is what we want
    !      to sort by).
    elemental function star_walker_eq (a, b) result(bEq)
        type(star_walker), intent(in) :: a, b
        logical :: bEq
        bEq = a%det .eq. b%det
    end function
    elemental function star_walker_ne (a, b) result(bNEq)
        type(star_walker), intent(in) :: a, b
        logical :: bNEq
        bNEq = a%det .ne. b%det
    end function
    elemental function star_walker_lt (a, b) result(bLt)
        type(star_walker), intent(in) :: a, b
        logical :: bLt
        bLt = a%det .lt. b%det
    end function
    elemental function star_walker_gt (a, b) result(bGt)
        type(star_walker), intent(in) :: a, b
        logical :: bGt
        bGt = a%det .gt. b%det
    end function
end module
