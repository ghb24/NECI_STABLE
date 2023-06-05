#include "macros.h"

module calcrho_mod
    use constants, only: dp
    better_implicit_none
    private
    public :: gethelement, igetexcitlevel, igetexcitlevel_2

contains

    integer pure function icmpdets(d1, d2, nel)
        integer, intent(in) :: d1(1:nel), d2(1:nel), nel
        integer :: i
        do i = 1, nel
            if (d1(i) < d2(i)) then
                icmpdets = -1
                return
            else if (d1(i) > d2(i)) then
                icmpdets = 1
                return
            end if
        end do
        icmpdets = 0
    end function

    integer pure function igetexcitlevel(ni, nj, nel)
        integer, intent(in) :: nI(nEl), nJ(nEl), nEl
        integer :: i, j, ic
        ic = 0
        do i = 1, nel
            do j = 1, nel
                if (ni(i) == nj(j)) then
                    ic = ic + 1
                    exit
                end if
            end do
        end do
        igetexcitlevel = nel - ic
    end function

!This routine is the same as IGetExcitLevel, but for two differences.
!First, this will only work for ordered lists of occupied orbitals
!Secondly, this will only find the number of excitation levels apart up
!to a maximum excitation level of MaxExcit. If the distance away is more than this
!it will simply return MaxExcit+1. i.e. if only want to know if it is a double
!excitation or less, then CALL IGetExcitLevel_2(NI,NJ,NEL,2). This will return
!0, 1, 2, or 3 if the distance is more than a double excitation.
    integer pure function IGetExcitLevel_2(NI, NJ, NEL, MaxExcit)
        integer, intent(in) :: nI(nEl), nJ(nEl), nEl, MaxExcit
        integer :: I, IC, j, MaxIC
        logical :: CommonOrb

        IF (MaxExcit >= NEl) THEN
            MaxIC = NEl - 1
        ELSE
            MaxIC = MaxExcit
        END IF
        I = 1   !Index for electron in NI
        IC = 0  !Number of excitation levels between determinants
        DO WHILE ((IC <= MaxIC) .and. (I <= NEl))
            CommonOrb = .false.
            do j = 1, NEL
            IF (NJ(j) == NI(I)) THEN
                CommonOrb = .true.
                EXIT
            ELSEIF (NJ(j) > NI(I)) THEN
                EXIT
            END IF
            end do
            IF (.not. CommonOrb) THEN
                IC = IC + 1
            END IF
            I = I + 1
        END DO
        IGetExcitLevel_2 = IC
    end function

!.. Very unnecessarily complex, but would be faster for large numbers of
!.. electrons.

    integer pure function igetexcitlevel_(ni, nj, nel)
        integer, intent(in) :: nI(nEl), nJ(nEl), nEl
        integer :: i, j, ic
!.. We only count differences from I to J
        do while (i <= nel .and. j <= nel)
            do while (ni(i) < nj(j) .and. i <= nel)
                i = i + 1
                ic = ic + 1
            end do
            do while (ni(i) > nj(j) .and. i <= nel .and. j <= nel)
                j = j + 1
            end do
            if (ni(i) == nj(j)) then
                i = i + 1
                j = j + 1
            end if
        end do
        ic = ic + (nel + 1 - i)
        igetexcitlevel_ = ic
    end function

!.. GETHELEMENT
!.. Get matrix element of the hamiltonian
    HElement_t(dp) pure FUNCTION GETHELEMENT(II, IJ, HAMIL, LAB, NROW, NDET)
        integer, intent(in) :: ii, ij, lab(*), nrow(nDet), ndet
        HElement_t(dp), intent(in) :: HAMIL(*)
        integer :: i, j, indxrow, imax, k
!.. We only have half of H, so if J<I, return the symmetrical (J,I) element
!.. Or if we have the whole H, it's quicker to look closer to its beginning
        if (ij < ii) then
            i = ij
            j = ii
        else
            i = ii
            j = ij
        end if
        gethelement = 0.0_dp
        indxrow = 1
!.. find the ith row
        do k = 1, i - 1
            indxrow = indxrow + nrow(k)
        end do
        imax = indxrow + nrow(i) - 1
        do k = indxrow, imax
            if (lab(k) > j) return
            if (lab(k) == j) then
                gethelement = hamil(k)
                return
            end if
        end do
    end function

end module calcrho_mod
