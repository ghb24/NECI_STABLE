module bit_reps
    use FciMCData, only: CurrentDets, WalkVecDets, MaxWalkersPart
    use SystemData, only: nel, tCSF, tTruncateCSF, nbasis, csf_trunc_level
    use CalcData, only: tTruncInitiator
    use csf_data, only: csf_yama_bit, csf_test_bit
    use constants, only: lenof_sign, end_n_int, bits_n_int, n_int
    use DetBitOps, only: count_open_orbs
    use bit_rep_data
    use SymExcitDataMod, only: excit_gen_store_type, tBuildOccVirtList, &
                               OrbClassCount, ScratchSize, SymLabelList2, &
                               SymLabelCounts2
    use sym_general_mod, only: ClassCountInd
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


    interface set_flag
        module procedure set_flag_single
        module procedure set_flag_general
    end interface

    ! Which decoding function do we want to use?
    interface decode_bit_det
!        module procedure decode_bit_det_bitwise
        module procedure decode_bit_det_chunks
    end interface

    ! Some (rather nasty) data for the chunkwise decoding
    integer, parameter :: decode_map_arr(0:8,0:255) = reshape(&
        (/0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,2,0,0,0,0,0,0,0,2,1,2,0,0,0,&
        0,0,0,1,3,0,0,0,0,0,0,0,2,1,3,0,0,0,0,0,0,2,2,3,0,0,0,0,0,0,3,1,2,&
        3,0,0,0,0,0,1,4,0,0,0,0,0,0,0,2,1,4,0,0,0,0,0,0,2,2,4,0,0,0,0,0,0,&
        3,1,2,4,0,0,0,0,0,2,3,4,0,0,0,0,0,0,3,1,3,4,0,0,0,0,0,3,2,3,4,0,0,&
        0,0,0,4,1,2,3,4,0,0,0,0,1,5,0,0,0,0,0,0,0,2,1,5,0,0,0,0,0,0,2,2,5,&
        0,0,0,0,0,0,3,1,2,5,0,0,0,0,0,2,3,5,0,0,0,0,0,0,3,1,3,5,0,0,0,0,0,&
        3,2,3,5,0,0,0,0,0,4,1,2,3,5,0,0,0,0,2,4,5,0,0,0,0,0,0,3,1,4,5,0,0,&
        0,0,0,3,2,4,5,0,0,0,0,0,4,1,2,4,5,0,0,0,0,3,3,4,5,0,0,0,0,0,4,1,3,&
        4,5,0,0,0,0,4,2,3,4,5,0,0,0,0,5,1,2,3,4,5,0,0,0,1,6,0,0,0,0,0,0,0,&
        2,1,6,0,0,0,0,0,0,2,2,6,0,0,0,0,0,0,3,1,2,6,0,0,0,0,0,2,3,6,0,0,0,&
        0,0,0,3,1,3,6,0,0,0,0,0,3,2,3,6,0,0,0,0,0,4,1,2,3,6,0,0,0,0,2,4,6,&
        0,0,0,0,0,0,3,1,4,6,0,0,0,0,0,3,2,4,6,0,0,0,0,0,4,1,2,4,6,0,0,0,0,&
        3,3,4,6,0,0,0,0,0,4,1,3,4,6,0,0,0,0,4,2,3,4,6,0,0,0,0,5,1,2,3,4,6,&
        0,0,0,2,5,6,0,0,0,0,0,0,3,1,5,6,0,0,0,0,0,3,2,5,6,0,0,0,0,0,4,1,2,&
        5,6,0,0,0,0,3,3,5,6,0,0,0,0,0,4,1,3,5,6,0,0,0,0,4,2,3,5,6,0,0,0,0,&
        5,1,2,3,5,6,0,0,0,3,4,5,6,0,0,0,0,0,4,1,4,5,6,0,0,0,0,4,2,4,5,6,0,&
        0,0,0,5,1,2,4,5,6,0,0,0,4,3,4,5,6,0,0,0,0,5,1,3,4,5,6,0,0,0,5,2,3,&
        4,5,6,0,0,0,6,1,2,3,4,5,6,0,0,1,7,0,0,0,0,0,0,0,2,1,7,0,0,0,0,0,0,&
        2,2,7,0,0,0,0,0,0,3,1,2,7,0,0,0,0,0,2,3,7,0,0,0,0,0,0,3,1,3,7,0,0,&
        0,0,0,3,2,3,7,0,0,0,0,0,4,1,2,3,7,0,0,0,0,2,4,7,0,0,0,0,0,0,3,1,4,&
        7,0,0,0,0,0,3,2,4,7,0,0,0,0,0,4,1,2,4,7,0,0,0,0,3,3,4,7,0,0,0,0,0,&
        4,1,3,4,7,0,0,0,0,4,2,3,4,7,0,0,0,0,5,1,2,3,4,7,0,0,0,2,5,7,0,0,0,&
        0,0,0,3,1,5,7,0,0,0,0,0,3,2,5,7,0,0,0,0,0,4,1,2,5,7,0,0,0,0,3,3,5,&
        7,0,0,0,0,0,4,1,3,5,7,0,0,0,0,4,2,3,5,7,0,0,0,0,5,1,2,3,5,7,0,0,0,&
        3,4,5,7,0,0,0,0,0,4,1,4,5,7,0,0,0,0,4,2,4,5,7,0,0,0,0,5,1,2,4,5,7,&
        0,0,0,4,3,4,5,7,0,0,0,0,5,1,3,4,5,7,0,0,0,5,2,3,4,5,7,0,0,0,6,1,2,&
        3,4,5,7,0,0,2,6,7,0,0,0,0,0,0,3,1,6,7,0,0,0,0,0,3,2,6,7,0,0,0,0,0,&
        4,1,2,6,7,0,0,0,0,3,3,6,7,0,0,0,0,0,4,1,3,6,7,0,0,0,0,4,2,3,6,7,0,&
        0,0,0,5,1,2,3,6,7,0,0,0,3,4,6,7,0,0,0,0,0,4,1,4,6,7,0,0,0,0,4,2,4,&
        6,7,0,0,0,0,5,1,2,4,6,7,0,0,0,4,3,4,6,7,0,0,0,0,5,1,3,4,6,7,0,0,0,&
        5,2,3,4,6,7,0,0,0,6,1,2,3,4,6,7,0,0,3,5,6,7,0,0,0,0,0,4,1,5,6,7,0,&
        0,0,0,4,2,5,6,7,0,0,0,0,5,1,2,5,6,7,0,0,0,4,3,5,6,7,0,0,0,0,5,1,3,&
        5,6,7,0,0,0,5,2,3,5,6,7,0,0,0,6,1,2,3,5,6,7,0,0,4,4,5,6,7,0,0,0,0,&
        5,1,4,5,6,7,0,0,0,5,2,4,5,6,7,0,0,0,6,1,2,4,5,6,7,0,0,5,3,4,5,6,7,&
        0,0,0,6,1,3,4,5,6,7,0,0,6,2,3,4,5,6,7,0,0,7,1,2,3,4,5,6,7,0,1,8,0,&
        0,0,0,0,0,0,2,1,8,0,0,0,0,0,0,2,2,8,0,0,0,0,0,0,3,1,2,8,0,0,0,0,0,&
        2,3,8,0,0,0,0,0,0,3,1,3,8,0,0,0,0,0,3,2,3,8,0,0,0,0,0,4,1,2,3,8,0,&
        0,0,0,2,4,8,0,0,0,0,0,0,3,1,4,8,0,0,0,0,0,3,2,4,8,0,0,0,0,0,4,1,2,&
        4,8,0,0,0,0,3,3,4,8,0,0,0,0,0,4,1,3,4,8,0,0,0,0,4,2,3,4,8,0,0,0,0,&
        5,1,2,3,4,8,0,0,0,2,5,8,0,0,0,0,0,0,3,1,5,8,0,0,0,0,0,3,2,5,8,0,0,&
        0,0,0,4,1,2,5,8,0,0,0,0,3,3,5,8,0,0,0,0,0,4,1,3,5,8,0,0,0,0,4,2,3,&
        5,8,0,0,0,0,5,1,2,3,5,8,0,0,0,3,4,5,8,0,0,0,0,0,4,1,4,5,8,0,0,0,0,&
        4,2,4,5,8,0,0,0,0,5,1,2,4,5,8,0,0,0,4,3,4,5,8,0,0,0,0,5,1,3,4,5,8,&
        0,0,0,5,2,3,4,5,8,0,0,0,6,1,2,3,4,5,8,0,0,2,6,8,0,0,0,0,0,0,3,1,6,&
        8,0,0,0,0,0,3,2,6,8,0,0,0,0,0,4,1,2,6,8,0,0,0,0,3,3,6,8,0,0,0,0,0,&
        4,1,3,6,8,0,0,0,0,4,2,3,6,8,0,0,0,0,5,1,2,3,6,8,0,0,0,3,4,6,8,0,0,&
        0,0,0,4,1,4,6,8,0,0,0,0,4,2,4,6,8,0,0,0,0,5,1,2,4,6,8,0,0,0,4,3,4,&
        6,8,0,0,0,0,5,1,3,4,6,8,0,0,0,5,2,3,4,6,8,0,0,0,6,1,2,3,4,6,8,0,0,&
        3,5,6,8,0,0,0,0,0,4,1,5,6,8,0,0,0,0,4,2,5,6,8,0,0,0,0,5,1,2,5,6,8,&
        0,0,0,4,3,5,6,8,0,0,0,0,5,1,3,5,6,8,0,0,0,5,2,3,5,6,8,0,0,0,6,1,2,&
        3,5,6,8,0,0,4,4,5,6,8,0,0,0,0,5,1,4,5,6,8,0,0,0,5,2,4,5,6,8,0,0,0,&
        6,1,2,4,5,6,8,0,0,5,3,4,5,6,8,0,0,0,6,1,3,4,5,6,8,0,0,6,2,3,4,5,6,&
        8,0,0,7,1,2,3,4,5,6,8,0,2,7,8,0,0,0,0,0,0,3,1,7,8,0,0,0,0,0,3,2,7,&
        8,0,0,0,0,0,4,1,2,7,8,0,0,0,0,3,3,7,8,0,0,0,0,0,4,1,3,7,8,0,0,0,0,&
        4,2,3,7,8,0,0,0,0,5,1,2,3,7,8,0,0,0,3,4,7,8,0,0,0,0,0,4,1,4,7,8,0,&
        0,0,0,4,2,4,7,8,0,0,0,0,5,1,2,4,7,8,0,0,0,4,3,4,7,8,0,0,0,0,5,1,3,&
        4,7,8,0,0,0,5,2,3,4,7,8,0,0,0,6,1,2,3,4,7,8,0,0,3,5,7,8,0,0,0,0,0,&
        4,1,5,7,8,0,0,0,0,4,2,5,7,8,0,0,0,0,5,1,2,5,7,8,0,0,0,4,3,5,7,8,0,&
        0,0,0,5,1,3,5,7,8,0,0,0,5,2,3,5,7,8,0,0,0,6,1,2,3,5,7,8,0,0,4,4,5,&
        7,8,0,0,0,0,5,1,4,5,7,8,0,0,0,5,2,4,5,7,8,0,0,0,6,1,2,4,5,7,8,0,0,&
        5,3,4,5,7,8,0,0,0,6,1,3,4,5,7,8,0,0,6,2,3,4,5,7,8,0,0,7,1,2,3,4,5,&
        7,8,0,3,6,7,8,0,0,0,0,0,4,1,6,7,8,0,0,0,0,4,2,6,7,8,0,0,0,0,5,1,2,&
        6,7,8,0,0,0,4,3,6,7,8,0,0,0,0,5,1,3,6,7,8,0,0,0,5,2,3,6,7,8,0,0,0,&
        6,1,2,3,6,7,8,0,0,4,4,6,7,8,0,0,0,0,5,1,4,6,7,8,0,0,0,5,2,4,6,7,8,&
        0,0,0,6,1,2,4,6,7,8,0,0,5,3,4,6,7,8,0,0,0,6,1,3,4,6,7,8,0,0,6,2,3,&
        4,6,7,8,0,0,7,1,2,3,4,6,7,8,0,4,5,6,7,8,0,0,0,0,5,1,5,6,7,8,0,0,0,&
        5,2,5,6,7,8,0,0,0,6,1,2,5,6,7,8,0,0,5,3,5,6,7,8,0,0,0,6,1,3,5,6,7,&
        8,0,0,6,2,3,5,6,7,8,0,0,7,1,2,3,5,6,7,8,0,5,4,5,6,7,8,0,0,0,6,1,4,&
        5,6,7,8,0,0,6,2,4,5,6,7,8,0,0,7,1,2,4,5,6,7,8,0,6,3,4,5,6,7,8,0,0,&
        7,1,3,4,5,6,7,8,0,7,2,3,4,5,6,7,8,0,8,1,2,3,4,5,6,7,8/),&
        (/9,256/) )

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
        WRITE(6,*) "Complex walkers in use."
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

        ! N.B. Flags MUST be last!!!!!
        !      If we change this bit, then we need to adjust ilut_lt and 
        !      ilut_gt.

        ! The total number of bits_n_int-bit integers used - 1
        NIfTot = NIfD + NIfY + NIfSgn + NIfFlag

        WRITE(6,*) "Setting integer length of determinants as bit-strings to: ", NIfTot + 1
        WRITE(6,*) "Setting integer bit-length of determinants as bit-strings to: ", bits_n_int
         
    end subroutine

    subroutine extract_bit_rep (ilut, nI, sgn, flags, store)
        
        ! Extract useful terms out of the bit-representation of a walker

        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer, intent(out) :: nI(nel), flags
        integer, dimension(lenof_sign), intent(out) :: sgn
        type(excit_gen_store_type), intent(inout), optional :: store

        if (tBuildOccVirtList .and. present(store)) then
            call decode_bit_det_lists (nI, ilut, store)
        else
            call decode_bit_det (nI, ilut)
        endif

        sgn = iLut(NOffSgn:NOffSgn+lenof_sign-1)
        IF(NifFlag.eq.1) THEN
            flags = iLut(NOffFlag)
        ELSE
            flags = 0
        ENDIF

    end subroutine extract_bit_rep

    pure subroutine extract_sign (ilut,sgn)
        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer, dimension(lenof_sign), intent(out) :: sgn

        sgn = iLut(NOffSgn:NOffSgn+lenof_sign-1)
    end subroutine extract_sign

    function extract_flags (iLut)
        integer(n_int), intent(in) :: ilut(0:nIfTot)
        integer :: extract_flags

        IF(NIfFlag.eq.1) THEN
            extract_flags = iLut(NOffFlag)
        ELSE
            extract_flags = 0
        ENDIF

    end function extract_flags

    function extract_part_sign (ilut, part_type) result(sgn)

        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(in) :: part_type
        integer :: sgn

        sgn = ilut(nOffSgn + part_type - 1)

    end function

    subroutine encode_bit_rep (ilut, Det, sgn, flag)
        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, dimension(lenof_sign), intent(in) :: sgn
        integer(n_int), intent(in) :: Det(0:NIfDBO)
        integer, intent(in) :: flag
        
        iLut(0:NIfDBO) = Det
        iLut(NOffSgn:NOffSgn+NIfSgn-1) = sgn
        IF(NIfFlag.eq.1) THEN
            iLut(NOffFlag) = flag
        ENDIF

    end subroutine encode_bit_rep

    subroutine encode_flags (ilut, flag)

        ! Add new flag information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flag

        iLut(NOffFlag) = flag

    end subroutine encode_flags

    subroutine clear_all_flags (ilut)

        ! Clear all of the flags

        integer(n_int), intent(inout) :: ilut(0:niftot)

        if (NIfFlag > 0) &
            iLut(NOffFlag:NOffFlag+NIfFlag-1) = 0

    end subroutine clear_all_flags

    subroutine encode_sign (ilut, sgn)

        ! Add new sign information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, dimension(lenof_sign), intent(in) :: sgn

        iLut(NOffSgn:NOffSgn+NIfSgn-1) = sgn

    end subroutine encode_sign

    subroutine encode_part_sign (ilut, sgn, part_type)
!This routine is like encode_sign, but it only encodes either the
!real OR imaginary sign of the walker, while leaving the other
!sign untouched. part_type=1 indicates real and part_type=2 imaginary.
!The sign argument is now a simple integer, rather than of dimension(lenof_sign).
        integer(n_int), intent(inout) :: ilut(0:NIfTot)
        integer, intent(in) :: sgn, part_type

        iLut(NOffSgn+part_type-1) = sgn

    end subroutine encode_part_sign
        

    subroutine set_flag_general (ilut, flg, state)

        ! Set or clear the specified flag (0 indexed) according to
        ! the value in state.
        !
        ! In:    flg   - Integer index of flag to set
        !        state - Flag will be set if state is true.
        ! InOut: ilut  - Bit representation of determinant

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flg
        logical, intent(in) :: state

        if (state) then
            call set_flag_single (ilut, flg)
        else
            call clr_flag (ilut, flg)
        endif
    end subroutine set_flag_general

    subroutine set_flag_single (ilut, flg)

        ! Set the specified flag (0 indexed) in the bit representation
        !
        ! In:    flg  - Integer index of flag to set
        ! InOut: ilut - Bit representation of determinant

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flg
!        integer :: off, ind

!        ind = NOffFlag + flg / bits_n_int
!        off = mod(flg, bits_n_int)
!        ilut(ind) = ibset(ilut(ind), off)
        
!This now assumes that we do not have more flags than bits in an integer.
         ilut(NOffFlag) = ibset(ilut(NOffFlag),flg)

    end subroutine set_flag_single

    subroutine copy_flag (ilut_src, ilut_dest, flg)

        ! Copy the selected flag between iluts

        integer(n_int), intent(in) :: ilut_src(0:niftot)
        integer(n_int), intent(inout) :: ilut_dest(0:niftot)
        integer, intent(in) :: flg
        logical :: state

        state = test_flag (ilut_src, flg)
        call set_flag_general (ilut_dest, flg, state)

    end subroutine


    subroutine clr_flag (ilut, flg)

        ! Clear the specified flag (0 indexed) in the bit representation
        !
        ! In:    flg  - Integer index of flag to clear
        ! InOut: ilut - Bit representation of determinant

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer, intent(in) :: flg
!        integer :: off, ind

!        ind = NOffFlag + flg / bits_n_int
!        off = mod(flg, bits_n_int)
!        ilut(ind) = ibclr(ilut(ind), off)

!This now assumes that we do not have more flags than bits in an integer.
        ilut(NOffFlag) = ibclr(ilut(NOffFlag),flg)

    end subroutine clr_flag

    ! function test_flag is in bit_rep_data
    ! This avoids a circular dependence with DetBitOps.

    subroutine encode_det (ilut, Det)

        ! Add new det information to a packaged walker.

        integer(n_int), intent(inout) :: ilut(0:nIfTot)
        integer(n_int), intent(in) :: Det(0:NIfDBO)

        iLut(0:NIfDBO) = Det

    end subroutine encode_det

    subroutine decode_bit_det_lists (nI, iLut, store)

        ! This routine decodes a determinant in bit form and constructs
        ! the natural ordered NEl integer representation of the det.
        !
        ! It also constructs lists of the occupied and unoccupied orbitals
        ! within a symmetry.

        integer(n_int), intent(in) :: iLut(0:niftot)
        integer, intent(out) :: nI(nel)
        type(excit_gen_store_type), intent(inout) :: store
        integer :: i, j, elec, orb, ind, virt(ScratchSize)

        ! Initialise the class counts
        store%ClassCountOcc = 0
        virt = 0

        elec = 0
        do i = 0, NIfD
            do j = 0, end_n_int
                orb = (i * bits_n_int) + (j + 1)
                ind = ClassCountInd(orb)
                if (btest(iLut(i), j)) then
                    !An electron is at this orbital
                    elec = elec + 1
                    nI(elec) = orb

                    ! Update class counts
                    store%ClassCountOcc(ind) = store%ClassCountOcc(ind) + 1

                    ! Store orbital INDEX in list of occ. orbs.
                    store%occ_list(store%ClassCountOcc(ind), ind) = elec

                    if (elec == nel) exit
                else
                    ! Update count
                    virt(ind) = virt(ind)+1

                    ! Store orbital in list of unocc. orbs.
                    store%virt_list(virt(ind), ind) = orb
                endif
            enddo
            if (elec == nel) exit
        enddo

        ! Give final class count
        store%ClassCountUnocc = OrbClassCount - store%ClassCountOcc
        store%tFilled = .true.
        store%scratch3(1) = -1

        ! Fill in the remaineder of the virtuals list
        forall (ind = 1:ScratchSize)
            !if (virt(ind) /= store%ClassCountUnocc(ind)) then
                store%virt_list ( &
                    virt(ind) + 1 : &
                    store%ClassCountUnocc(ind), ind) = &
                SymLabelList2 (&
                    SymLabelCounts2(1, ind) + virt(ind) + &
                        store%ClassCountOcc(ind) : &
                    SymLabelCounts2(1, ind) + OrbClassCount(ind) - 1)
            !endif
        endforall

    end subroutine


    subroutine decode_bit_det_chunks (nI, iLut)

        ! This is a routine to take a determinant in bit form and construct
        ! the natural ordered Nel integer form of the det.
        ! If CSFs are enabled, transfer the Yamanouchi symbol as well.

        integer(n_int), intent(in) :: ilut(0:NIftot)
        integer, intent(out) :: nI(nel)
        integer :: nopen, i, j, k, val, elec, offset, pos
        logical :: bIsCsf

        ! We need to use the CSF decoding routine if CSFs are enabled, and 
        ! we are below a truncation limit if set.
        bIsCsf = .false.
        if (tCSF) then
            if (tTruncateCSF) then
                nopen = count_open_orbs(ilut)
                if (nopen <= csf_trunc_level) then
                    bIsCsf = .true.
                endif
            else
                bIsCsf = .true.
            endif
        endif

        elec = 0
        if (bIsCsf) then
            ! ****************
            ! Currently this just works in the old fashioned way. We aren't
            ! really that worried about CSF efficiency atm.
            ! ****************
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
                        if (.not.btest(iLut(i),ieor(j,1))) then
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
            offset = 0
            do i = 0, NIfD
                do j = 0, bits_n_int - 1, 8
                    val = iand(ishft(ilut(i), -j), Z'FF')
                    do k = 1, decode_map_arr(0, val)
                        elec = elec + 1
                        nI(elec) = offset + decode_map_arr(k, val)
                        if (elec == nel) exit
                    enddo
                    if (elec == nel) exit
                    offset = offset + 8
                enddo
                if (elec == nel) exit
            enddo

        endif

    end subroutine


    subroutine decode_bit_det_bitwise (nI, iLut)

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
                        if (.not.btest(iLut(i),ieor(j,1))) then
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
    end subroutine decode_bit_det_bitwise

end module bit_reps
