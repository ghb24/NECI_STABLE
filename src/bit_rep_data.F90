#include "macros.h"
module bit_rep_data

    use CalcData, only: tUseRealCoeffs
    use constants

    implicit none

    ! Structure of a bit representation:

    ! | 0-NIfD: Det | Yamanouchi | Sign(Re) | Sign(Im) | Flags |
    !
    ! -------
    ! (NIfD + 1) * 64-bits              Orbital rep.
    !  NIfY      * 32-bits              Yamanouchi symbol
    !  1         * 32-bits              Signs (Re)
    ! (1         * 32-bits if needed)   Signs (Im)
    ! (1         * 32-bits if needed)   Flags

    integer :: nIfTot  ! Upper bound of bit representation. In form 0:NIfTot
    integer :: nIfD    ! Final byte representing spatial/spin orbitals

    integer :: nOffY   ! Offset of Yamanouchi symbol. n.b. always on a 64bit
                       ! boundary --> This is in 64-bit integers.
                       ! TODO: ensure Yamanouchi symbols are 32bit.
    integer :: nIfY    ! Number of bytes to represent a Yamanouchi symbol

    integer :: nIfDBO  ! Size used for bit operations (e.g. sorting)

    integer :: nOffFlag   ! Offset of flags. in bytes
    integer :: nIfFlag    ! Number of bytes to contain flags.

    integer :: nOffSgn  ! Offset of signs in integers
    integer :: nIfSgn   ! Number of integers used for signs
    integer :: nIfTotKP ! Upper bound of krylov_vecs.

    integer :: nIfBCast ! Size of data to use in annihilation broadcast

    ! Has the RDM component of the bit representation (For bcast) been inited.
    logical :: bit_rdm_init
    integer :: nOffParent

    ! Somewhere to store the spawning parent coefficient (for funky
    ! initiator thresholds).
    integer :: nOffParentCoeff, nIfParentCoeff

    ! Flags which we can store
    ! RT_M_Merge: Adapted real-time flags
    logical :: tUseFlags

    integer :: flag_counter

    integer, parameter :: flag_deterministic = 0, &
                          flag_determ_parent = 1, &
                          flag_trial = 2, &
                          flag_connected = 3, &
                          flag_has_been_initiator(1) = 4
                          ! RT_M_Merge: These should only be adressed with __REALTIME
                          ! use these unused to mark diagonal "spawns"
#ifdef __PROG_NUMRUNS
    integer, parameter :: flag_initiator(lenof_sign_max) &
                            = (/6, 7, 8, 10, 11, 12, 13, 14, 15, &
                                16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26/), &
                          flag_adi_checked = 27, &
                          flag_static_init(lenof_sign_max) &
                            = (/28, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, &
                                42, 43, 44, 45, 46, 47, 48/), &
                          num_flags = 49
#else
    integer, parameter :: flag_initiator(2) = (/ 6, 7/), &
                          flag_adi_checked = 8, &
                          flag_static_init(2) = (/9, 10/), &
                          num_flags = 11
#endif

contains

    pure function test_flag (ilut, flg) result(bSet)

        ! Test specified flag (0 indexed) in the bit representation.
        !
        ! In:  flg  - Integer index of flag to test
        !      ilut - Bit representation of determinant
        ! Ret: bSet - returns .true. if the flag is set, false otherwise


        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer, intent(in) :: flg
        logical :: bSet

        !Commented out code is for when we need multiple integers for storing flags (unlikely!)
!        ind = NOffFlag + flg / bits_n_int
!        off = mod(flg, bits_n_int)

!        bSet = btest(ilut(ind), off)

        if(tUseFlags) then
            bSet = btest(ilut(NOffFlag), flg)
        else
            bSet = .false.
        endif

    end function test_flag

    pure subroutine extract_sign (ilut, real_sgn)
        integer(n_int), intent(in) :: ilut(0:nIfTot)
        real(dp), intent(out) :: real_sgn(lenof_sign)
        integer(n_int) :: sgn(lenof_sign)

        sgn = iLut(NOffSgn:NOffSgn+lenof_sign-1)
        ! transfer operates elementwise
        real_sgn = transfer(sgn, real_sgn)

    end subroutine extract_sign

end module
