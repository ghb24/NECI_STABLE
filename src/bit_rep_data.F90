module bit_rep_data

    use constants, only: n_int, bits_n_int

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

    ! Flags which we can store
    integer, parameter :: flag_is_initiator(2) = (/0,1/), &
                          flag_parent_initiator(2) = (/0,1/), & ! n.b. the same
                          flag_make_initiator(2) = (/2,3/), &
                          flag_deterministic = 4

    ! Bit mask with all bits unset except for bit 4, corresponding to the
    ! deterministic flag. This is used in ilut_lt and ilut_gt to order the
    ! states such that deterministic states are at the top of the sorted list.
    integer(n_int) :: deterministic_mask

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
        integer :: off, ind

        !Commented out code is for when we need multiple integers for storing flags (unlikely!)
!        ind = NOffFlag + flg / bits_n_int
!        off = mod(flg, bits_n_int)

!        bSet = btest(ilut(ind), off)
        bSet = btest(ilut(NOffFlag), flg)

    end function test_flag

end module
