module bit_reps
    use FciMCData, only: CurrentDets, WalkVecDets, MaxWalkersPart
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

    integer :: nOffF   ! Offset of flags. in bytes
    integer :: nIfF    ! Number of bytes to contain flags.

    integer :: nOffSgn  ! Offset of signs in integers
    integer :: nOffSgnI ! Offset of complex part of signs
    integer :: nIfSgn   ! Number of integers used for signs

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

        ! This indicates the upper-bound for the determinants when expressed
        ! in bit-form. This will equal int(nBasis/32).
        ! The actual total length for a determinant in bit form will be
        ! NoIntforDet+1 + nIfY (which is the size of the Yamanouchi Symbol
        nIfD = int(nbasis / bits_n_int)

        ! Could use only 32-bits for this, except that it makes it very
        ! tricky to do do anything like sorting, as the latter 32-bits of the
        ! integer would contain random junk.
        NOffY = NIfTot + 1
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
        nOffSgn = NOffY + NIfY
        nIfSgn = 1
#ifdef __CMPLX
#ifdef __INT64
        nOffSgnI = nOffSgn
#else
        nOffSgnI = nOffSgn + 1
        nIfSgn = nIfSgn + 1
#endif
#endif

        ! The number of integers used for sorting / other bit manipulations
        NIfDBO = NIfD + NIfY

        ! The total number of bits_n_int-bit integers used - 1
        NIfTot = NIfD + NIfY

    end subroutine

    subroutine extract_bit_rep (ilut, nI, sgn, flags)
        
        ! Extract useful terms out of the bit-representation of a walker

        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer, intent(out) :: nI(nel), sgn, flags
        integer :: off, len

        call decode_bit_det (nI, ilut)

        off = int(nOffSgn / size_n_int)
        sft = mod(nOffSgn / size_n_int)
        sgn = ishft(ilut(off), -8*sft)

        off = int(nOffG / size_n_int)
        sft = mod(nOffF / size_n_int)
        sgn = ishft(ilut(off), -8*sft)
        
    end subroutine extract_bit_rep

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
