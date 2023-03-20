!.. A stub containing dummy functions in CPMD, which will appear on linking appropriately.
#include "macros.h"
module cpmdstub_mod
    implicit none

contains
    SUBROUTINE INITFINDXI(I, J, K, L, UMATEL, gen2ints)
        use constants, only: dp
        IMPLICIT NONE
        INTEGER I, J, K, L
        real(dp) UMATEL
        logical gen2ints

! Make warnings go away
        i = i; j = j; k = k; l = l; gen2ints = gen2ints
        UMATEL = 0.0_dp
    END

    SUBROUTINE CPMDANTISYMINTEL(G1, UMat2D, HarInt, NStates)
        use constants, only: dp
        use SystemData, only: BasisFN
        use util_mod, only: stop_all
        IMPLICIT NONE
        type(BasisFN) :: G1(*)
        complex(dp) :: HarInt(*)
        integer :: nStates
        HElement_t(dp) :: UMat2D(*)
!eliminate warnings
        G1(1) = G1(1)
        nStates = nStates
        UMat2D(1) = UMat2D(1)
        HarInt(1) = HarInt(1)
        call stop_all("CPMDANTISYMINTEL", "Entering stub routine.")
    END SUBROUTINE

    SUBROUTINE CPMDGETOCC(ISTATE, OCC)
        use constants, only: dp
        IMPLICIT NONE
!         INCLUDE 'elct.inc'
        real(dp) OCC
        INTEGER ISTATE
!         COMMON /ELCT/ F
!         OCC=F(ISTATE)
        OCC = 0.0_dp
        istate = istate
    END

    subroutine CalcHarPIInts(HarInt, nStates)
        use constants, only: dp
        implicit none
        integer nStates
        complex(dp) HarInt(nStates, nStates)

! Eliminate Warnings
        HarInt = HarInt; nStates = nStates
    end subroutine CalcHarPIInts

    subroutine KPntSymInt(i, j, k, l, a, b, c, d)
        implicit none
        integer i, j, k, l, a, b, c, d
! Make warnings go away
        i = i; j = j; k = k; l = l; a = a; b = b; c = c; d = d
    end subroutine KPntSymInt
end module
