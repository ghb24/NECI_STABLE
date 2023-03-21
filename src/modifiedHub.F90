#include "macros.h"

module breathing_Hub
    use SystemData, only: momIndexTable, nBasisMax, G1, breathingCont
    use sym_mod, only: MomPbcSym
    use constants, only: pi, dp

contains

    subroutine setupMomIndexTable()
        implicit none
        integer :: li(3), i, kx, ky, kz, s
        ! fill the table containing the index of a given momentum
        ! such that knowing q allows for direct lookup of the matrix element
        do i = 1, 3
            li(i) = nBasisMax(i, 2) - nBasisMax(i, 1) + 1
        end do
        allocate(momIndexTable(li(1), li(2), li(3), 2), stat=i)
        do s = 1, 2
            do kz = 1, li(3)
                do ky = 1, li(2)
                    do kx = 1, li(1)
                        momIndexTable(kx, ky, kz, s) = kx + (ky - 1) * li(1) + &
                                                       (kz - 1) * li(1) * li(2) + (s - 1) * li(1) * li(2) * li(3)
                    end do
                end do
            end do
        end do
    end subroutine setupMomIndexTable

    elemental function bHubIndexFunction(i, j, k, l) result(bInd)
        ! this function gets the index of the matrix element of the breathing term
        ! corresponding to the four states i,j,k,l

        ! not sure if it is too expensive to do the index generation on the fly
        ! maybe storing all matrix elements ignoring redundancies might be worthwile
        ! even though memory cost will grow significantly (but still no noteworthy memory usage)
        implicit none
        integer, intent(in) :: i, j, k, l
        integer :: bInd, buf(3), li(3), tgt
        ! to get the breathing effect, we need to know the exchanged momentum
        ! therefore, compare the momentum of i to that of k/l with the matching
        ! spin
        unused_var(j)

        unused_var(j)

        if (G1(i)%MS == G1(k)%MS) then
            tgt = k
        else
            tgt = l
        end if
        buf = G1(i)%k - G1(tgt)%k
        call MomPbcSym(buf, nBasisMax)
        bInd = momIndexTable(buf(1), buf(2), buf(3), G1(i)%MS)
    end function bHubIndexFunction

    subroutine setupBreathingCont(prefactor)
        implicit none
        integer :: li(3), i, kx, ky, kz, q(3), s
        real(dp) :: bcBuf, prefactor
        do i = 1, 3
            li(i) = nBasisMax(i, 2) - nBasisMax(i, 1) + 1
        end do
        allocate(breathingCont(li(1) * li(2) * li(3) * 2), stat=s)
        do s = 1, 2
            do kz = 1, li(3)
                do ky = 1, li(2)
                    do kx = 1, li(1)
                        ! momentum is the only input to the matrix elements of the breathing
                        ! contribution
                        q = (/kx, ky, kz/)
                        call MomPbcSym(q, nBasisMax)
                        bcBuf = 0.0_dp
                        do i = 1, 3
                            if (li(i) > 1) bcBuf = bcBuf + cos(2 * pi / li(i) * q(i))
                        end do
                        breathingCont(momIndexTable(kx, ky, kz, s)) = prefactor * bcBuf
                    end do
                end do
            end do
        end do
    end subroutine setupBreathingCont

end module breathing_Hub
