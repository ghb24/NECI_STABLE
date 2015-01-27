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

    ! Flags which we can store
    logical :: tUseflags
    integer, parameter :: flag_deterministic = 0, &
                          flag_determ_parent = 1, &
                          flag_trial = 2, &
                          flag_connected = 3, &
                          flag_has_been_initiator(1) = 4, &
                          flag_unused1 = 5, &
                          flag_unused2 = 6, &
                          flag_unused3 = 7, &
                          flag_ic0_spawn = 8, &
                          flag_death_done = 9, &
                          flag_negative_sign = 10

#ifdef __PROG_NUMRUNS
    integer, parameter :: flag_initiator(lenof_sign_max) &
                            = (/11, 12, 13, 14, 15, 16, 17, 18, 19, 20, &
                                21, 22, 23, 24, 25, 26, 27, 28, 29, 30/), &
                          num_flags = 31
#else
    integer, parameter :: flag_initiator(2) = (/11,12/), &
                          num_flags = 13
#endif

    ! IMPORTANT
    integer, parameter :: flag_bit_offset = bits_n_int - num_flags
    integer(n_int), parameter :: sign_mask = ishft(not(0_n_int), -num_flags), &
                                 flags_mask = not(sign_mask), &
                                 sign_neg_mask = ibset(sign_mask, &
                                          flag_bit_offset + flag_negative_sign)

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

        bSet = .false.

#if defined(__INT64) && !defined(__PROG_NUMRUNS)
        if ((.not. tUseRealCoeffs) .or. tUseFlags) then
#else
        if (tUseFlags) then
#endif
            bSet = btest(ilut(NOffFlag), flg + flag_bit_offset)
        end if

    end function test_flag

    pure subroutine extract_sign (ilut, real_sgn)
        integer(n_int), intent(in) :: ilut(0:nIfTot)
        real(dp), intent(out) :: real_sgn(lenof_sign)
        integer(n_int) :: sgn(lenof_sign)

#if defined(__INT64) && !defined(__PROG_NUMRUNS)
        if (tUseRealCoeffs) then
            sgn = ilut(NOffSgn:NOffSgn+lenof_sign-1)
            real_sgn = transfer(sgn, real_sgn)
        else
            ! TODO: Should we inline the flag test
            sgn(1) = int(iand(ilut(NOffSgn), sign_mask), sizeof_int)
            if (test_flag(ilut, flag_negative_sign)) sgn(1) = -sgn(1)
            if (lenof_sign == 2) then
                sgn(lenof_sign) = int(ilut(NOffSgn+1), sizeof_int)
            end if
            real_sgn = real(sgn, dp)
        end if
#else
        sgn = iLut(NOffSgn:NOffSgn+lenof_sign-1)
        real_sgn = transfer(sgn, real_sgn)
#endif
    end subroutine extract_sign

end module
