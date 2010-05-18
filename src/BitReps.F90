module bit_reps
    use FciMCData, only: CurrentDets, WalkVecDets, MaxWalkersPart
    use SystemData, only: nel, tCSF, tTruncateCSF, nbasis, csf_trunc_level
    use CalcData, only: tTruncInitiator
    use csf_data, only: csf_yama_bit, csf_test_bit
    use constants, only: lenof_sign, end_n_int, bits_n_int, n_int
    use DetBitOps, only: count_open_orbs
    use bit_rep_data
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

contains

    subroutine allocate_currentdets ()
        
        ! Allocate memory of the correct size for the currentdets array.

        integer :: ierr
        character(*), parameter :: this_routine = 'allocate_currentdets'
        
        allocate (WalkVecDets(0:NIfTot, MaxWalkersPart), stat=ierr)
        if (ierr /= 0) &
            call stop_all (this_routine, "Allocation failed for WalkVecDets")

    end subroutine allocate_currentdets

    subroutine init_bit_rep ()

        ! Set the values of nifd etc.

        character(*), parameter :: this_routine = 'init_bit_rep'

        ! This indicates the upper-bound for the determinants when expressed
        ! in bit-form. This will equal int(nBasis/32).
        ! The actual total length for a determinant in bit form will be
        ! NoIntforDet+1 + nIfY (which is the size of the Yamanouchi Symbol
        nIfD = int(nbasis / bits_n_int)

        ! Could use only 32-bits for this, except that it makes it very
        ! tricky to do do anything like sorting, as the latter 32-bits of the
        ! integer would contain random junk.
        NOffY = NIfD + 1
        if (tCSF) then
            if (tTruncateCSF) then
                NIfY = int(csf_trunc_level / bits_n_int) + 1
            else
                NIfY = int(nel / bits_n_int) + 1
            endif
        else
            NIfY = 0
        endif
        if (NIfY > 1) &
            call stop_all (this_routine, "CSFs with more than bits_n_int &
                          &open-shell electrons are not supported, and are &
                          &probably not a good idea.")

        ! The signs array
        NOffSgn = NOffY + NIfY
        NIfSgn = 1
#ifdef __CMPLX
        NIfSgn = NIfSgn + 1     !TODO: If __INT64, adjust packing into one integer
#endif

        ! The number of integers used for sorting / other bit manipulations
        NIfDBO = NIfD + NIfY

        ! Integers for flags
        if (tTruncInitiator) then
            !If there are other options which require flags, then this criteria must be extended.
            !However, do not increase this value from one, since we should only need max one integer
            !for flags, and this is hardcoded in elsewhere.
            NIfFlag = 1
        else
            NIfFlag = 0
        endif
        NOffFlag = NOffSgn + NIfSgn

        ! The total number of bits_n_int-bit integers used - 1
        NIfTot = NIfD + NIfY + NIfSgn + NIfFlag

        WRITE(6,*) "Setting integer length of determinants as bit-strings to: ", NIfTot + 1
      WRITE(6,*) "Setting integer bit-length of determinants as bit-strings to: ", bits_n_int

         
    end subroutine

    subroutine extract_bit_rep (ilut, nI, sgn, flags)
        
        ! Extract useful terms out of the bit-representation of a walker

        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer, intent(out) :: nI(nel), flags
        integer, dimension(lenof_sign), intent(out) :: sgn

        call decode_bit_det (nI, ilut)

        sgn = iLut(NOffSgn:NOffSgn+lenof_sign-1)
        IF(NOffFlag.eq.1) THEN
            flags = iLut(NOffFlag)
        ELSE
            flags = 0
        ENDIF

    end subroutine extract_bit_rep

    function extract_sign (ilut)
        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer, dimension(lenof_sign), intent(out) :: extract_sign

        extract_sign = iLut(NOffSgn:NOffSgn+lenof_sign-1)
    end subroutine extract_sign

    function extract_flags (iLut)
        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer, intent(out) :: flags

        IF(NOffFlag.eq.1) THEN
            flags = iLut(NOffFlag)
        ELSE
            flags = 0
        ENDIF

    end function extract_flags

    subroutine encode_flags (ilut, flag)

        ! Add new flag information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flag

        iLut(NOffFlag) = flag

    end subroutine encode_flags

    subroutine encode_sign (ilut, sgn)

        ! Add new sign information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, dimension(lenof_sign), intent(in) :: sgn

        iLut(NOffSgn:NOffSgn+NIfSgn-1) = sgn

    end subroutine encode_sign

    subroutine encode_det (ilut, Det)

        ! Add new det information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer(n_int), intent(in) :: Det(0:NIfDBO)

        iLut(0:NIfDBO) = Det

    end subroutine encode_det

    subroutine decode_bit_det (nI, iLut)

        ! This is a routine to take a determinant in bit form and construct 
        ! the natural ordered NEl integer forim of the det.
        ! If CSFs are enabled, transfer the yamanouchi symbol as well.

        integer(n_int), intent(in) :: iLut(0:NIfTot)
        integer, intent(out) :: nI(nel)
        integer :: i, j, elec, pos, nopen
        logical :: bIsCsf

        ! We need to use the CSF decoding routine if CSFs are enable, and we
        ! are below a truncation limit if set.
        bIsCsf = .false.
        if (tCSF) then
            if (tTruncateCSF) then
                nopen = count_open_orbs(iLut)
                if (nopen <= csf_trunc_level) then
                    bIsCsf = .true.
                endif
            else
                bIsCsf = .true.
            endif
        endif

        elec=0
        if (bIsCsf) then
            ! Consider the closed shell electrons first
            do i=0,NIfD
                do j=0,bits_n_int-2,2
                    if (btest(iLut(i),j)) then
                        if (btest(iLut(i),j+1)) then
                            ! An electron pair is in this spatial orbital
                            ! (2 matched spin orbitals)
                            elec = elec + 2
                            nI(elec-1) = (bits_n_int*i) + (j+1)
                            nI(elec) = (bits_n_int*i) + (j+2)
                            if (elec == nel) return
                        endif
                    endif
                enddo
            enddo

            ! Now consider the open shell electrons
            ! TODO: can we move in steps of two, to catch unmatched pairs?
            nopen = 0
            do i=0,NIfD
                do j=0,end_n_int
                    if (btest(iLut(i),j)) then
                        if (.not.btest(iLut(i),xor(j,1))) then
                            elec = elec + 1
                            nI(elec) = (bits_n_int*i) + (j+1)
                            pos = NIfD + 1 + (nopen/bits_n_int)
                            if (btest(iLut(Pos),mod(nopen,bits_n_int))) then
                                nI(elec) = ibset(nI(elec),csf_yama_bit)
                            endif
                            nopen = nopen + 1
                        endif
                    endif
                    if (elec==nel) exit
                enddo
                if (elec==nel) exit
            enddo
            ! If there are any open shell e-, set the csf bit
            nI = ibset(nI, csf_test_bit)
        else
            do i=0,NIfD
                do j=0,end_n_int
                    if(btest(iLut(i),j)) then
                        !An electron is at this orbital
                        elec=elec+1
                        nI(elec)=(i*bits_n_int)+(j+1)
                        if (elec == nel) exit
                    endif
                enddo
                if (elec == nel) exit
            enddo
        endif
    end subroutine decode_bit_det

end module bit_reps
