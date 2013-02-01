#include "macros.h"

module DetBitOps

    ! A collection of useful operations to perform on the bit-representation
    ! of determinants.

    use Systemdata, only: nel, tCSF, tTruncateCSF, csf_trunc_level
    use CalcData, only: tTruncInitiator
    use bit_rep_data, only: NIfY, NIfTot, NIfD, NOffFlag, NIfFlag, &
                            test_flag, flag_is_initiator,NIfDBO,NOffSgn
    use csf_data, only: iscsf, csf_yama_bit, csf_orbital_mask, csf_test_bit
    use constants, only: n_int,bits_n_int,end_n_int,dp,lenof_sign,sizeof_int

    implicit none

#ifdef __INT64
    ! 10101010 and 01010101 in binary respectively.
!    integer(n_int), parameter :: MaskBeta=Z'5555555555555555'
    integer(n_int), parameter :: MaskBeta=6148914691236517205_n_int
    integer(n_int), parameter :: MaskAlpha=IShft(MaskBeta,1)
#else
    integer(n_int), parameter :: MaskBeta=1431655765_n_int
    integer(n_int), parameter :: MaskAlpha=-1431655766_n_int
#endif

    ! Which count-bits procedure do we use?
    ! http://gurmeetsingh.wordpress.com/2008/08/05/fast-bit-counting-routines/
    ! for a variety of interesting bit counters
    interface CountBits
        !module procedure CountBits_sparse
        !module procedure CountBits_nifty
        module procedure CountBits_elemental
    end interface

    contains

    ! This will count the bits set in a bit-string up to a number nBitsMax, if
    ! provided.
    ! The function will return 0 -> nBitsMax+1
    ! A value of nBitsMax+1 indicates that more bits are set than was expected.
    ! The total number of set bits can exceed nBitsMax+1, however.
    ! Counts bits set in integer array (0:nLast)
    pure integer function CountBits_sparse (iLut, nLast, nBitsMax)
        integer, intent(in), optional :: nBitsMax
        integer, intent(in) :: nLast 
        integer(kind=n_int), intent(in) :: iLut(0:nLast)
        integer(kind=n_int) :: iLutTemp(0:nLast)
        integer :: i, lnBitsMax

        ! By default, allow all the bits to be set
        if (present(nBitsMax)) then
            lnBitsMax = nBitsMax
        else
            lnBitsMax = bits_n_int * (nLast+1)
        endif

        CountBits_sparse = 0
        iLutTemp = iLut
        do i=0,nLast
            do while((iLutTemp(i).ne.0).and.(CountBits_sparse.le.lnBitsMax))
                ! Clear the rightmost set bit
                iLutTemp(i)=IAND(iLutTemp(i),iLutTemp(i)-1)
                CountBits_sparse = CountBits_sparse + 1
            enddo
            if(CountBits_sparse .gt. lnBitsMax) return
        enddo
    end function CountBits_sparse

    ! Try counting using a nifty bit of bitwise arithmetic
    ! See comments for CountBits_sparse and count_set_bits.
    pure integer function Countbits_nifty (iLut, nLast, nBitsMax)
        integer, intent(in), optional :: nBitsMax
        integer, intent(in) :: nLast
        integer(kind=n_int), intent(in) :: iLut(0:nLast)
        integer :: i, lnBitsMax

        ! By default, allow all the bits to be set
        if (present(nBitsMax)) then
            lnBitsMax = nBitsMax
        else
            lnBitsMax = bits_n_int * (nLast+1)
        endif

        CountBits_nifty = 0
        do i=0,nLast
            CountBits_nifty = CountBits_nifty + count_set_bits(iLut(i))
            if (CountBits_nifty .gt. lnBitsMax) then
                CountBits_nifty = lnBitsmax+1
                return
            endif
        enddo       

    end function CountBits_nifty

    ! Using elemental routines rather than an explicit do-loop. Should be
    ! faster.
    pure function CountBits_elemental (iLut, nLast, nBitsMax) result(nbits)
        integer, intent(in), optional :: nBitsMax
        integer, intent(in) :: nLast
        integer(kind=n_int), intent(in) :: iLut(0:nLast)
        integer :: nbits

        nbits = sum(count_set_bits(iLut))
        
        !No advantage to test for this!
!        if (present(nBitsMax)) nbits = min(nBitsmax+1, nbits)
    end function

    ! An elemental routine which will count the number of bits set in one 
    ! (32 bit) integer. We can do similar things for 8bit, 16bit and 64bit.
    ! This makes use of the same counting trick as CountBits_nifty. As nicely
    ! summarised by James:
    !
    ! The general idea is to use a divide and conquer approach.
    ! * Set each 2 bit field to be the sum of the set bits in the two single
    !   bits originally in that field.
    ! * Set each 4 bit field to be the sum of the set bits in the two 2 bit
    !   fields originally in the 4 bit field.
    ! * Set each 8 bit field to be the sum of the set bits in the two 4 bit
    !   fields it contains.
    ! * etc.
    ! Thus we obtain an algorithm like:
    !     x = ( x & 01010101...) + ( (x>>1) & 01010101...)
    !     x = ( x & 00110011...) + ( (x>>2) & 00110011...)
    !     x = ( x & 00001111...) + ( (x>>4) & 00001111...)
    ! etc., where & indicates AND and >> is the shift right operator.
    ! Further optimisations are:
    ! * Any & operations can be omitted where there is no danger that
    ! a field's sum will carry over into the next field.
    ! * The first line can be replaced by:
    !     x = x - ( (x>>1) & 01010101...)
    !   thanks to the population (number of set bits) in an integer
    !   containing p bits being given by:
    !     pop(x) = \sum_{i=0}^{p-1} x/2^i
    ! * Summing 8 bit fields together can be performed via a multiplication
    !   followed by a right shift.
    elemental function count_set_bits (a) result (nbits)
        integer(n_int), intent(in) :: a
        integer :: nbits
        integer(n_int) :: tmp

#ifdef __INT64
        integer(n_int), parameter :: m1 = 6148914691236517205_n_int  !Z'5555555555555555'
        integer(n_int), parameter :: m2 = 3689348814741910323_n_int  !Z'3333333333333333'
        integer(n_int), parameter :: m3 = 1085102592571150095_n_int  !Z'0f0f0f0f0f0f0f0f'
        integer(n_int), parameter :: m4 = 72340172838076673_n_int    !Z'0101010101010101'

        ! For 64 bit integers:
        tmp = a - iand(ishft(a,-1), m1)
        tmp = iand(tmp, m2) + iand(ishft(tmp,-2), m2)
        tmp = iand(tmp, m3) + iand(ishft(tmp,-4), m3)
        nbits = int(ishft(tmp*m4, -56),sizeof_int)
#else
        integer(n_int), parameter :: m1 = 1431655765_n_int    !Z'55555555'
        integer(n_int), parameter :: m2 = 858993459_n_int     !Z'33333333'
        integer(n_int), parameter :: m3 = 252645135_n_int     !Z'0F0F0F0F'
        integer(n_int), parameter :: m4 = 16843009_n_int      !Z'01010101'

        ! For 32 bit integers:
        tmp = a - iand(ishft(a,-1), m1)
        tmp = iand(tmp, m2) + iand(ishft(tmp, -2), m2)
        tmp = iand((tmp+ishft(tmp, -4)), m3) * m4
        nbits = ishft(tmp, -24)
#endif

    end function count_set_bits

    pure integer function count_open_orbs (iLut)
        
        ! Returns the number of unpaired electrons in the determinant.
        !
        ! In:  iLut (0:NIfD) - Source bit det

        integer(kind=n_int), intent(in) :: iLut(0:NIfD)
        integer(kind=n_int), dimension(0:NIfD) :: alpha, beta

        alpha = iand(iLut, MaskAlpha)
        beta = iand(iLut, MaskBeta)
        alpha = ishft(alpha, -1)
        alpha = ieor(alpha, beta)
        
        count_open_orbs = CountBits(alpha, NIfD)
    end function


    pure function FindBitExcitLevel(iLutnI, iLutnJ, maxExLevel) result(IC)

        ! Find the excitation level of one determinant relative to another
        ! given their bit strings (the number of orbitals they differ by)
        !
        ! In:  iLutnI, iLutnJ    - The bit representations
        !      maxExLevel        - An (optional) maximum ex level to consider
        ! Ret: FindBitExcitLevel - The number of orbitals i,j differ by

        integer(kind=n_int), intent(in) :: iLutnI(0:NIfD), iLutnJ(0:NIfD)
        integer, intent(in), optional :: maxExLevel
        integer(kind=n_int) :: tmp(0:NIfD)
        integer :: IC

        ! Obtain a bit string with only the excited orbitals one one det.
        tmp = ieor(iLutnI, iLutnJ)
        tmp = iand(iLutnI, tmp)

        ! Then count them
        ! Since our CountBits routines don't actually make a saving
        ! for counting smaller numbers of bits, no point in even testing
        ! for a maxExLevel!
!        if (present(maxExLevel)) then
!            IC = CountBits(tmp, NIfD, maxExLevel)
!        else
            IC = CountBits(tmp, NIfD)
!        endif

    end function FindBitExcitLevel

    function FindSpatialBitExcitLevel (iLutI, iLutJ, maxExLevel) result(IC)
        
        ! Find the excitation level of one determinant relative to another
        ! given their bit strings, ignoring the spin components of orbitals.
        ! (i.e. the number of spatial orbitals they differ by)
        !
        ! In:  iLutI, iLutJ - The bit representations
        !      maxExLevel   - An (optional) maximum ex level to consider
        ! Ret: IC           - The numbero f orbitals i,j differ by

        integer(kind=n_int), intent(in) :: iLutI(0:NIfD), iLutJ(0:NIfD)
        integer, intent(in), optional :: maxExLevel
        integer :: IC
        integer(kind=n_int), dimension(0:NIfD,2) :: alpha, beta, sing, doub, tmp

        ! Obtain the alphas and betas
        alpha(:,1) = iand(ilutI, MaskAlpha)
        alpha(:,2) = iand(ilutJ, MaskAlpha)
        beta(:,1) = iand(ilutI, MaskBeta)
        beta(:,2) = iand(ilutJ, MaskBeta)

        ! Bit strings separating doubles, and singles shifted to beta pos.
        doub = iand(beta, ishft(alpha, -1))
        doub = ior(doub, ishft(doub, +1))
        sing = ieor(beta, ishft(alpha, -1))

        ! Doubles and singles shifted to betas. Obtain unique orbitals ...
        tmp = ior(doub, sing)
        tmp(:,1) = ieor(tmp(:,1), tmp(:,2))
        tmp(:,1) = iand(tmp(:,1), tmp(:,2))

        ! ... and count them.
        if (present(maxExLevel)) then
            IC = CountBits(tmp(:,1), NIfD, maxExLevel)
        else
            IC = CountBits(tmp(:,1), NIfD)
        endif
    end function FindSpatialBitExcitLevel

    !WARNING - I think this *may* be buggy - use with caution - ghb24 8/6/10
    pure subroutine get_bit_excitmat (ilutI, iLutJ, ex, IC)
        
        ! Obatin the excitation matrix between two determinants from their bit
        ! representation without calculating tSign --> a bit quicker.
        !
        ! In:    iLutI, iLutJ - Bit representations of determinants I,J
        ! InOut: IC           - Specify max IC before bailing, and return
        !                       number of orbital I,J differ by
        ! Out:   ex           - Excitation matrix between I,J

        integer(kind=n_int), intent(in) :: iLutI(0:NIfD), iLutJ(0:NIfD)
        integer, intent(inout)  :: IC
        integer, intent(out), dimension(2,IC) :: ex

        integer(kind=n_int) :: ilut(0:NIfD,2)
        integer :: pos(2), max_ic, i, j, k

        ! Obtain bit representations of I,J containing only unique orbitals
        ilut(:,1) = ieor(ilutI, ilutJ)
        ilut(:,2) = iand(ilutJ, ilut(:,1))
        ilut(:,1) = iand(ilutI, ilut(:,1))

        max_ic = IC
        pos = 0
        IC = 0
        do i=0,NIfD
            do j=0,bits_n_int
                do k=1,2
                    if (pos(k) < max_ic) then
                        if (btest(ilut(i,k), j)) then
                            pos(k) = pos(k) + 1
                            IC = max(IC, pos(k))
                            ex(k, pos(k)) = bits_n_int*i + j + 1
                        endif
                    endif
                enddo
                if (pos(1) >= max_ic .and. pos(2) >= max_ic) return
            enddo
        enddo
    end subroutine
    
    subroutine get_bit_open_unique_ind (iLutI, iLutJ, op_ind, nop, &
                                        tsign_id, nsign, IC)

                                        ! TODO: comment
        ! Obtain the indices of unique open orbitals in I and J. 
        !
        ! In:  ILutI, ILutJ - Bit representations of determinants
        !      IC           - (Max) number of orbitals for I,J to differ by
        ! Out: op_ind       - Array of unique single indices for I,J
        !      nop          - Number of unique singles in each of I,J

        integer(kind=n_int), intent(in) :: iLutI(0:NIfD), iLutJ(0:NIfD)
        integer, intent(in) :: IC
        integer, intent(out) :: op_ind(2*IC, 2), nop(2)
        integer, intent(out) :: tsign_id (2*IC,2), nsign(2)

        integer :: i, j, det, sing_ind(2) 
        integer(kind=n_int) :: ilut(0:NIfD,2), sing(0:NIfD,2)
        integer(kind=n_int) :: alpha(0:NIfD), beta(0:NIfD)

        ! Obtain all the singles in I,J
        alpha = iand(iLutI, MaskAlpha)
        beta = iand(iLutI, MaskBeta)
        alpha = ishft(alpha, -1)
        sing(:,1) = ieor(alpha, beta)

        alpha = iand(iLutJ, MaskAlpha)
        beta = iand(iLutJ, MaskBeta)
        alpha = ishft(alpha, -1)
        sing(:,2) = ieor(alpha, beta)

        ! Obtain bit representations of I,J with only the differing orbitals.
        ilut(:,1) = ieor(sing(:,1), sing(:,2))
        ilut(:,2) = iand(ilut(:,1), sing(:,2))
        ilut(:,1) = iand(ilut(:,1), sing(:,1))

        nop = 0
        sing_ind = 0
        nsign = 0
        do i=0,NIfD
            do j=0,end_n_int
                ! TODO: If CSF, increment in steps of 2.
                do det=1,2
                    if (nop(det) < 2*IC) then
                        ! Update the singles index
                        if (btest(sing(i,det), j)) &
                            sing_ind(det) = sing_ind(det) + 1

                        if (btest(ilut(i,det), j)) then
                            ! If unique single, store its index.
                            nop(det) = nop(det) + 1
                            op_ind(nop(det),det) = sing_ind(det)

                            ! If single comes from a double in the other det,
                            ! then it affects tSign when permuted.
                            ! TODO: tidy and compact
                            if (det == 1) then
                                if (btest(iLutJ(i), ieor(j,1))) then
                                    nsign(1) = nsign(1) + 1
                                    tsign_id(nsign(1),1) = sing_ind(1)
                                endif
                            else if (det == 2) then
                                if (btest(iLutI(i), ieor(j,1))) then
                                    nsign(2) = nsign(2) + 1
                                    tsign_id(nsign(2),2) = sing_ind(2)
                                endif
                            endif
                        endif
                    endif
                enddo

                if (nop(1) >= 2*IC .and. nop(2) >= 2*IC) return
            enddo
        enddo
    end subroutine get_bit_open_unique_ind


    ! This will return true if iLutI is identical to iLutJ and will return 
    ! false otherwise.
    pure function DetBitEQ(iLutI,iLutJ,nLast) result(res)
        integer, intent(in), optional :: nLast
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        logical :: res
        integer :: i, lnLast

        if(iLutI(0).ne.iLutJ(0)) then
            res=.false.
            return
        else
            if (present(nLast)) then
                lnLast = nLast
            else
                lnLast = NIfDBO
            endif

            do i=1,lnLast
                if(iLutI(i).ne.iLutJ(i)) then
                    res=.false.
                    return
                endif
            enddo
        endif
        res=.true.
    end function DetBitEQ

    pure function sign_lt (ilutI, ilutJ) result (bLt)

        ! This is a comparison function between two bit strings of length 
        ! 0:NIfTot, and will return true if absolute value of the sign of 
        ! ilutI is less than ilutJ

        integer(n_int), intent(in) :: iLutI(0:), iLutJ(0:)
        logical :: bLt
        integer :: SignI(lenof_sign), SignJ(lenof_sign)
        real(dp) :: RealSignI(lenof_sign), RealSignJ(lenof_sign)
        real(dp) :: WeightI,WeightJ

        ! Extract the sign. Ensure that we convert to a real, as the integers
        ! themselves mean nothing.
        SignI = int(iLutI(NOffSgn:NOffSgn+lenof_sign-1),sizeof_int)
        SignJ = int(iLutJ(NOffSgn:NOffSgn+lenof_sign-1),sizeof_int)
        RealSignI = transfer(SignI, RealSignI)
        RealSignJ = transfer(SignJ, RealSignJ)

        if(lenof_sign == 1) then
            bLt = abs(RealSignI(1)) < abs(RealSignJ(1))
        else
            WeightI = sqrt(real(RealSignI(1), dp)**2 + &
                           real(RealSignI(lenof_sign), dp)**2)
            WeightJ = sqrt(real(RealSignJ(1), dp)**2 + &
                           real(RealSignJ(lenof_sign), dp)**2)

            bLt = WeightI < WeightJ
        endif
    end function sign_lt

    pure function sign_gt (ilutI, ilutJ) result (bGt)

        ! This is a comparison function between two bit strings of length 
        ! 0:NIfTot, and will return true if the abs sign of ilutI is greater
        ! than ilutJ

        integer(n_int), intent(in) :: iLutI(0:), iLutJ(0:)
        logical :: bGt
        integer :: SignI(lenof_sign), SignJ(lenof_sign)
        real(dp) :: RealSignI(lenof_sign), RealSignJ(lenof_sign)
        real(dp) :: WeightI, WeightJ

        SignI = int(iLutI(NOffSgn:NOffSgn+lenof_sign-1),sizeof_int)
        SignJ = int(iLutJ(NOffSgn:NOffSgn+lenof_sign-1),sizeof_int)
        RealSignI = transfer(SignI, RealSignI)
        RealSignJ = transfer(SignJ, RealSignJ)

        if(lenof_sign == 1) then
            bGt = abs(RealSignI(1)) > abs(RealSignJ(1))
        else
            WeightI = sqrt(real(RealSignI(1), dp)**2 + &
                           real(RealSignI(lenof_sign), dp)**2)
            WeightJ = sqrt(real(RealSignJ(1), dp)**2 + &
                           real(SignJ(lenof_sign), dp)**2)

            bGt = WeightI > WeightJ
        endif
    end function sign_gt

    pure function ilut_lt (ilutI, ilutJ) result (bLt)
!        use util_mod, only: operator(.arrlt.)

        ! A slightly subtler sort than DetBitLt.
        ! Sort the iluts integer by integer, up to the determinant. 
        ! Ignore flag and sign differences.

        integer(n_int), intent(in) :: iLutI(0:), iLutJ(0:)
        integer :: i
        logical :: bLt

        ! Sort by the first item first ...
        do i = 0, NIfDBO    !   NOffFlag - 1
            if (iLutI(i) /= iLutJ(i)) exit
        enddo

        !! Make the comparison
!        if (i >= NOffFlag) then
        if (i > NIfDBO) then
            bLt = .false.
!            if (tTruncInitiator) then
!                !if initiator, sort first by real flag, the imaginary.
!                if (test_flag(ilutI, flag_is_initiator(1)) .and. &
!                    .not. test_flag(ilutJ, flag_is_initiator(1))) then
!                        !I<J if real i is initiator, and j is not
!                        bLt = .true.
!                elseif((test_flag(ilutI, flag_is_initiator(1)).eqv. test_flag(ilutJ, flag_is_initiator(1))) &
!                    .and.(test_flag(ilutI, flag_is_initiator(2))).and..not.test_flag(ilutJ, flag_is_initiator(2))) then
!                        !if real flags the same, I<J if imaginary i is initiator and j is not.
!                        bLt = .true.
!                endif
!
!            endif
        else
            bLt = ilutI(i) < ilutJ(i)
        endif

    end function

    pure function ilut_gt (iLutI, iLutJ) result(bGt)
!        use util_mod, only: operator(.arrgt.)

        ! A slightly subtler sort than DetBitGt.
        ! Sort the iluts integer by integer. If we get to the flags, and they
        ! are occupied, then sort initiators as 'less than' non-initiators
        !
        ! --> Initiators appear earlier in a list than non-initiators

        integer(n_int), intent(in) :: iLutI(0:), iLutJ(0:)
        integer :: i
        logical :: bGt


        !bGt = iLutI .arrgt. iLutJ
        
        ! Sort by the first item first ...
        do i = 0, NIfDBO    !   NOffFlag - 1
            if (ilutI(i) /= iLutJ(i)) exit
        enddo

        ! Make the comparison
!        if (i >= NOffFlag) then
        if (i > NIfDBO) then
            bGt = .false.
!            if (tTruncInitiator) then
!                if (.not. test_flag(ilutI, flag_is_initiator(1)) .and. &
!                    test_flag(ilutJ, flag_is_initiator(1))) then
!                    bGt = .true.
!                elseif ((test_flag(ilutI, flag_is_initiator(1)) .eqv. test_flag(ilutJ, flag_is_initiator(1))) &
!                    .and. (.not. test_flag(ilutI, flag_is_initiator(2)).and.test_flag(ilutJ, flag_is_initiator(2)))) then
!                    !If real flags same, sort by imaginary flags.
!                    bGt = .true.
!                endif
!            endif
        else
            bGt = ilutI(i) > ilutJ(i)
        endif

    end function

    ! This will return true if the determinant has been set to zero, and 
    ! false otherwise.
    pure logical function DetBitZero(iLutI,nLast)
        integer, intent(in), optional :: nLast
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot)
        integer :: i, lnLast
        if(iLutI(0).ne.0) then
            DetBitZero=.false.
            return
        else
            if (present(nLast)) then
                lnLast = nLast
            else
                lnLast = NIftot
            endif
            do i=1,lnLast
                if(iLutI(i).ne.0) then
                    DetBitZero=.false.
                    return
                endif
            enddo
        endif
        DetBitZero=.true.
    end function DetBitZero


    ! This will return 1 if iLutI is "less" than iLutJ, 0 if the determinants
    ! are identical, or -1 if iLutI is "more" than iLutJ
    pure integer function DetBitLT(iLutI,iLutJ,nLast)
        integer, intent(in), optional :: nLast
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer :: i, lnLast

        !First, compare first integers
        IF(iLutI(0).lt.iLutJ(0)) THEN
            DetBitLT=1
        ELSEIF(iLutI(0).eq.iLutJ(0)) THEN
            ! If the integers are the same, then cycle through the rest of 
            ! the integers until we find a difference.
            ! If we don't want to consider all the integers, specify nLast
            if (present(nLast)) then
                lnLast = nLast
            else
                lnLast = NIfDBO
            endif

            do i=1,lnLast
                IF(iLutI(i).lt.iLutJ(i)) THEN
                    DetBitLT=1
                    RETURN
                ELSEIF(iLutI(i).gt.iLutJ(i)) THEN
                    DetBitLT=-1
                    RETURN
                ENDIF
            enddo
            DetBitLT=0
        ELSE
            DetBitLT=-1
        ENDIF
    END FUNCTION DetBitLT

    ! This will return 1 if iLutI is "less" than iLutJ, or -1 if iLutI is 
    ! "more" than iLutJ.  If these are identical, this routine looks at 
    ! iLut2I and iLut2J, and returns 1 if iLut2I is "less" than iLut2J, -1 
    ! if iLut2I is "more than iLut2J, and 0 if these are still identical.
    integer function Det2BitLT(iLutI,iLutJ,iLut2I,iLut2J,nLast)
        integer, intent(in), optional :: nLast
        integer :: i,lnLast
        integer(kind=n_int) :: iLutI(0:NIfTot),iLutJ(0:NIfTot)
        integer(kind=n_int) :: iLut2I(0:NIfTot),iLut2J(0:NIfTot)

        IF(iLutI(0).lt.iLutJ(0)) THEN
            ! First, compare first integers
            Det2BitLT=1
            RETURN
        ELSEIF(iLutI(0).gt.iLutJ(0)) THEN
            Det2BitLT=-1
            RETURN
        ELSEIF(iLutI(0).eq.iLutJ(0)) THEN
            ! If the integers are the same, then cycle through the rest of 
            ! the integers until we find a difference.
            if (present(nLast)) then
                lnLast = nLast
            else
                lnLast = NIfDBO
            endif
            do i=1,lnLast
                IF(iLutI(i).lt.iLutJ(i)) THEN
                    Det2BitLT=1
                    RETURN
                ELSEIF(iLutI(i).gt.iLutJ(i)) THEN
                    Det2BitLT=-1
                    RETURN
                ENDIF
            enddo
            ! If we get through this loop without RETURN-ing, iLutI and iLutJ
            ! are identical, so look to iLut2I and iLut2J            
            IF(iLut2I(0).lt.iLut2J(0)) THEN
                Det2BitLT=1
                RETURN
            ELSEIF(iLut2I(0).gt.iLut2J(0)) THEN
                Det2BitLT=-1
                RETURN
            ELSEIF(iLut2I(0).eq.iLut2J(0)) THEN
                do i=1,lnLast
                    IF(iLut2I(i).lt.iLut2J(i)) THEN
                        Det2BitLT=1
                        RETURN
                    ELSEIF(iLut2I(i).gt.iLut2J(i)) THEN
                        Det2BitLT=-1
                        RETURN
                    ENDIF
                enddo
            ENDIF
        ENDIF
        !If we still have not returned, both determinants are identical. 
        Det2BitLT=0
    END FUNCTION Det2BitLT

    ! This will return 1 if iLutI is "less" than iLutJ, 0 if the determinants 
    ! are identical, or -1 if iLutI is "more" than iLutJ
    ! This particular version checks excitation level initially, then only if
    ! these are the same does it move on to determinants.
    integer function DetExcitBitLT(iLutI,iLutJ,iLutHF,nLast)
        integer, intent(in), optional :: nLast
        integer(kind=n_int), intent(in) :: iLutI(0:NIftot), iLutJ(0:NIfTot)
        integer(kind=n_int), intent(in) :: iLutHF(0:NIfTot)
        integer i, ExcitLevelI, ExcitLevelJ,lnLast
        
        ExcitLevelI = FindBitExcitLevel(iLutI, iLutHF, nel)
        ExcitLevelJ = FindBitExcitLevel(iLutJ, iLutHF, nel)

        ! First order in terms of excitation level.  I.e. if the excitation 
        ! levels are different, we don't care what the determinants are we 
        ! just order in terms of the excitation level.
        IF(ExcitLevelI.lt.ExcitLevelJ) THEN
            DetExcitBitLT=1
            RETURN
        ELSEIF(ExcitLevelI.gt.ExcitLevelJ) THEN
            DetExcitBitLT=-1
            RETURN

        ! If the excitation levels are the same however, we need to look at 
        ! the determinant and order according to this.            
        ELSEIF(ExcitLevelI.eq.ExcitLevelJ) THEN
            ! First, compare first integers
            IF(iLutI(0).lt.iLutJ(0)) THEN
                DetExcitBitLT=1
                RETURN
            ELSEIF(iLutI(0).eq.iLutJ(0)) THEN
                ! If the integers are the same, then cycle through the rest 
                ! of the integers until we find a difference.
                if (present(nLast)) then
                    lnLast = nLast
                else
                    lnLast = NIfDBO
                endif
                do i=1,lnLast
                    IF(iLutI(i).lt.iLutJ(i)) THEN
                        DetExcitBitLT=1
                        RETURN
                    ELSEIF(iLutI(i).gt.iLutJ(i)) THEN
                        DetExcitBitLT=-1
                        RETURN
                    ENDIF
                enddo
            ELSE
                DetExcitBitLT=-1
                RETURN
            ENDIF
            ! If it gets through all this without being returned then the 
            ! two determinants are equal and DetExcitBitLT=0
            DetExcitBitLT=0
        ENDIF
    END FUNCTION DetExcitBitLT

    ! This is a routine to encode a determinant as natural ordered integers
    ! (nI) as a bit string (iLut(0:NIfTot)) where NIfD=INT(nBasis/32)
    ! If this is a csf, the csf is contained afterwards.
    pure subroutine EncodeBitDet(nI,iLut)
        integer, intent(in) :: nI(nel)
        integer(kind=n_int), intent(out) :: iLut(0:NIfTot)
        integer :: i, det, pos, nopen
        logical :: open_shell

        iLut(:)=0
        nopen = 0
        open_shell = .false.
        if (tCSF .and. iscsf (nI)) then
            do i=1,nel
                ! THe first non-paired orbital has yama symbol = 1
                if ((.not. open_shell) .and. &
                    btest(nI(i), csf_yama_bit)) open_shell = .true.

                ! Set the bit in the bit representation
                det = iand(nI(i), csf_orbital_mask)
                iLut((det-1)/bits_n_int) = ibset(iLut((det-1)/bits_n_int),mod(det-1,bits_n_int))

                if (open_shell) then
                    if (btest(nI(i), csf_yama_bit)) then
                        pos = NIfD + 1 + (nopen/bits_n_int)
                        iLut(pos) = ibset(iLut(pos), mod(nopen,bits_n_int))
                    endif
                    nopen = nopen + 1
                endif
            enddo
        else
            do i=1,nel
                pos = (nI(i) - 1) / bits_n_int
                iLut(pos)=ibset(iLut(pos),mod(nI(i)-1,bits_n_int))
            enddo
        endif
    end subroutine EncodeBitDet


    pure function spatial_bit_det (ilut) result(ilut_s)

        ! Convert the spin orbital representation in ilut_s into a spatial
        ! orbital representation, with all singly occupied orbitals in the 
        ! 'beta' position.
        !
        ! In:  ilut   - Spin orbital, bit representation
        ! Out: ilut_s - Spatial orbital, bit representation. Loses all sign
        !               etc. info (i.e. for ints > NIfD --> 0)

        integer(n_int), intent(in) :: ilut(0:NIfTot)
        integer(n_int) :: ilut_s (0:NIfTot)
        integer(n_int), dimension(0:NIfD) :: alpha, beta, a_sft, b_sft

        ! Obtain alpha/beta orbital representations
        alpha = iand(ilut(0:NIfD), MaskAlpha)
        beta = iand(ilut(0:NIfD), MaskBeta)

        ! Shift alphas to beta pos and vice-versa
        a_sft = ishft(alpha, -1)
        b_sft = ishft(beta, +1)

        ! Obtain representation with all singly occupied orbitals in the beta
        ! position, and doubly occupied orbitals doubly occupied
        ilut_s(NIfD+1:NIfTot) = 0
        ilut_s(0:NIfD) = ior(beta, ior(a_sft, iand(b_sft, alpha)))

    end function


    subroutine FindExcitBitDet(iLutnI, iLutnJ, IC, ExcitMat)

        ! This routine will find the bit-representation of an excitation by
        ! constructing the new ilut from the old one and the excitation matrix
        !
        ! In:  iLutnI (0:NIfD) - source bit det
        !      IC              - Excitation level
        !      ExcitMat(2,2)   - Excitation Matrix
        ! Out: iLutnJ (0:NIfD) - New bit det

        integer, intent(in) :: IC, ExcitMat(2,2)
        integer(kind=n_int), intent(in) :: iLutnI (0:NIfTot)
        integer(kind=n_int), intent(inout) :: iLutnJ (0:NIfTot)
        integer :: pos(2,2), bit(2,2), i

        iLutnJ = iLutnI
        if (IC == 0) then
            if (.not.tCSF) then
                call stop_all ("FindExcitBitDet", 'Invalid excitation level')
            endif
        else
            ! Which integer and bit in ilut represent each element?
            pos = (excitmat - 1) / bits_n_int
            bit = mod(excitmat - 1, bits_n_int)

            ! Clear bits for excitation source, and set bits for target
            do i=1,IC
                iLutnJ(pos(1,i)) = ibclr(iLutnJ(pos(1,i)), bit(1,i))
                iLutnJ(pos(2,i)) = ibset(iLutnJ(pos(2,i)), bit(2,i))
            enddo
        endif

    end subroutine FindExcitBitDet

    subroutine shift_det_bit_singles_to_beta (iLut)
        integer(kind=n_int), intent(inout) :: iLut(0:NIfD)
        integer(kind=n_int) :: iA(0:NIfD), iB(0:NIfD)

        ! Extract the betas
        iB = iand(iLut, MaskBeta)
        ! Extract the alphas and shift them into beta positions.
        iA = ishft(iand(iLut, MaskAlpha), -1)
        
        ! Generate the doubles
        iLut = iand(iB, iA)
        iLut = ior(iLut, ishft(iLut, 1))

        ! Generate the singles and include in result
        iLut = ior(iLut, ieor(iA, iB))
    end subroutine

    ! Test if all of the beta singles are in higher numbered orbitals than
    ! the alpha singles.
    logical function is_canonical_ms_order (nI)
        integer, intent(in) :: nI(nel)
        integer(kind=n_int), dimension(0:NIfTot) :: alpha, beta, tmp
        integer :: first_beta_byte, first_beta_bit
        integer i

        call EncodeBitDet(nI, alpha)
        beta = iand(alpha, MaskBeta)
        tmp = iand(alpha, MaskAlpha)

        alpha = iand(tmp, not(ishft(beta,1))) ! Only alpha singles
        beta = iand(beta, not(ishft(tmp,-1)))  ! Only beta singles

        ! Find the first non-zero beta byte
        is_canonical_ms_order = .false.
        do i=0,NIfD
            if (beta(i) /= 0) exit
        enddo
        if (i > NIfD) return
        first_beta_byte = i

        ! Find the last non-zero alpha byte
        do i=NIfD,first_beta_byte,-1
            if (alpha(i) /= 0) exit            
        enddo

        if (i < first_beta_byte) then
            is_canonical_ms_order = .true.
        else
            ! Now we need to consider the bits.
            do i=0,end_n_int
                if (btest(beta(first_beta_byte),i)) exit
            enddo
            first_beta_bit = i

            ! TODO: steps of 2, as alpha/beta even/odd...
            do i=end_n_int,first_beta_bit,-1
                if (btest(alpha(first_beta_byte), i)) exit
            enddo
            if (i < first_beta_bit) is_canonical_ms_order = .true.
        endif
    end function

    pure function TestClosedShellDet (ilut) result(tClosed)

        ! Is the determinant closed shell?

        integer(n_int), intent(in) :: iLut(0:NIfTot)
        integer(n_int) :: alpha(0:NIfD), beta(0:NIfD)
        logical :: tClosed

        ! Separate alphas and betas
        alpha = iand(ilut(0:NIfD), MaskAlpha)
        beta = iand(ilut(0:NIfD), MaskBeta)

        ! Shift and XOR to eliminate paired electrons
        alpha = ieor(beta, ishft(alpha, -1))

        ! Are there any remaining unpaired electrons?
        tClosed = all(alpha == 0)

    end function TestClosedShellDet

    ! Routine to count number of open *SPATIAL* orbitals in a bit-string 
    ! representation of a determinant.
    ! ************************
    ! BROKEN
    ! NOTE: This function name is misleading
    !       It counts the number of unpaired Beta electrons (ignores Alpha)
    !       --> Returns nopen/2 <==> Ms=0
    ! ************************
    pure SUBROUTINE CalcOpenOrbs(iLut,OpenOrbs)
        INTEGER(kind=n_int) :: iLutAlpha(0:NIfD),iLutBeta(0:NIfD)
        integer(n_int), intent(in) :: ilut(0:NIfD)
        integer, intent(out) :: OpenOrbs
        INTEGER :: i
        
        iLutAlpha(:)=0
        iLutBeta(:)=0

        do i=0,NIfD     
                    
            iLutAlpha(i)=IAND(iLut(i),MaskAlpha)    !Seperate the alpha and beta bit strings
            iLutBeta(i)=IAND(iLut(i),MaskBeta)
            iLutAlpha(i)=ISHFT(iLutAlpha(i),-1)     !Shift all alpha bits to the left by one.

            iLutAlpha(i)=NOT(iLutAlpha(i))              ! This NOT means that set bits are now represented by 0s, not 1s
            iLutAlpha(i)=IAND(iLutAlpha(i),iLutBeta(i)) ! Now, only the 1s in the beta string will be counted.

        enddo

        OpenOrbs = CountBits(iLutAlpha,NIfD,NEl)
    END SUBROUTINE CalcOpenOrbs

    function IsAllowedHPHF (ilut, sym_ilut) result (bAllowed)

        ! Is the specified determinant an 'allowed' HPHF function (i.e. can
        ! it be found in the determinant list, or is it the symmetry paired
        ! one?

        integer(n_int), intent(in) :: ilut(0:NIfD)
        integer(n_int) :: ilut_tmp(0:NIfD)
        integer(n_int), intent(out), optional :: sym_ilut(0:NIfD)
        logical :: bAllowed

        if (TestClosedShellDet(ilut)) then
            bAllowed = .true.
        else
            call spin_sym_ilut (ilut, ilut_tmp)
            if (DetBitLt(ilut, ilut_tmp, NIfD) > 0) then
                bAllowed = .false.
            else
                bAllowed = .true.
            endif

            if (present(sym_ilut)) sym_ilut = ilut_tmp
        endif

    end function

    pure subroutine spin_sym_ilut (ilutI, ilutJ)

        ! Generate the spin-coupled determinant of ilutI in ilutJ. Performs
        ! the same operation as FindDetSpinSym rather more concisely.

        integer(n_int), intent(in) :: ilutI(0:NIfD)
        integer(n_int), intent(out) :: ilutJ(0:NIfD)
        integer(n_int) :: ilut_tmp(0:NIfD)

        ilut_tmp = ishft(iand(ilutI, MaskAlpha), -1)
        ilutJ = ishft(iand(ilutI, MaskBeta), +1)
        ilutJ = ior(ilutJ, ilut_tmp)

    end subroutine


    pure function get_single_parity (ilut, src, tgt) result(par)

        ! Find the relative parity of two determinants, where one is ilut
        ! and the other is a single excitation of ilut where orbital src is
        ! swapped with orbital tgt.

        integer, intent(in) :: src, tgt
        integer(n_int), intent(in) :: ilut(0:NIfTot)

        integer(n_int) :: tmp(0:NIfD), mask(0:NIfD)
        integer :: min_orb, max_orb, par, min_int, max_int, cnt
        integer :: min_in_int, max_in_int

        ! We just want to count the orbitals between these two limits.
        min_orb = (min(src, tgt) + 1) - 1
        max_orb = (max(src, tgt) - 1)

        ! Which integers of the bit representation are involved?
        min_int = int(min_orb / bits_n_int)
        max_int = int(max_orb / bits_n_int)

        ! Where in the integer do the revelant bits sit?
        min_in_int = mod(min_orb, bits_n_int)
        max_in_int = mod(max_orb, bits_n_int)

        ! Generate a mask so as to only count the occupied orbitals
        ! between where we started and the end.
        mask(0:min_int-1) = 0
        mask(min_int:max_int) = not(0_n_int)
        mask(max_int+1:NIfD) = 0
        mask(min_int) = &
            iand(mask(min_int), ishft(not(0_n_int), min_in_int))
        mask(max_int) = &
            iand(mask(max_int), not(ishft(not(0_n_int), max_in_int)))

        ! Count the number of occupied orbitals between the source and tgt
        ! orbitals.
        cnt = CountBits(iand(mask, ilut(0:NIfD)), NIfD)

        ! Get the parity from this information.
        if (btest(cnt, 0)) then
            par = -1
        else
            par = 1
        end if

    end function

    function get_double_parity (ilut, src, tgt) result(par)

        ! Find the relative parity of two determinants, where one is ilut
        ! and the other is a single excitation of ilut where orbital src is
        ! swapped with orbital tgt.

        integer, intent(in) :: src(2), tgt(2)
        integer(n_int), intent(in) :: ilut(0:NIfTot)

        integer(n_int) :: tmp(0:NIfD), mask(0:NIfD)
        integer :: min_orb, max_orb, par, min_int, max_int, cnt
        integer :: min_in_int, max_in_int


        if (all(tgt > maxval(src)) .or. all(tgt < minval(src))) then

            ! The source and target orbitals don't overlap
            par = get_single_parity (ilut, src(1), src(2)) * &
                  get_single_parity (ilut, tgt(1), tgt(2))

        !elseif ((minval(src) < minval(tgt)) .eqv. &
        !                                 (maxval(src) > maxval(tgt))) then
        else

            ! All categories of overlapping src and target orbitals are the 
            ! same.
            par = get_single_parity (ilut, minval(src), minval(tgt)) * &
                  get_single_parity (ilut, maxval(src), maxval(tgt))

        end if


    end function


end module

    pure subroutine GetBitExcitation(iLutnI,iLutnJ,Ex,tSign)

        ! A port from hfq. The first of many...
        ! JSS.

        ! In:
        !    iLutnI(basis_length): bit string representation of the Slater
        !        determinant.
        !    iLutnJ(basis_length): bit string representation of the Slater
        !        determinant.
        !    Ex(1,1): contains the maximum excitation level, max_excit, to be
        !        considered.
        ! Out:
        !    Ex(2,max_excit): contains the excitation connected iLutnI to
        !        iLutnJ.  Ex(1,i) is the i-th orbital excited from and Ex(2,i)
        !        is the corresponding orbital excited to.
        !    tSign:
        !        True if an odd number of permutations is required to line up
        !        the determinants.

        use SystemData, only: nel
        use bit_rep_data, only: NIfD
        use DetBitOps, only: CountBits_nifty
        use constants, only: n_int,bits_n_int,end_n_int
        implicit none
        integer(kind=n_int), intent(in) :: iLutnI(0:NIfD), iLutnJ(0:NIfD)
        integer, intent(inout) :: Ex(2,*)
        logical, intent(out) :: tSign
        integer :: i, j, iexcit1, iexcit2, perm, iel1, iel2, max_excit
        logical :: testI, testJ

        tSign=.true.
        max_excit = Ex(1,1)
        Ex(:,:max_excit) = 0

        if (max_excit > 0) then

            iexcit1 = 0
            iexcit2 = 0
            iel1 = 0
            iel2 = 0
            perm = 0

            ! Finding the permutation to align the determinants is non-trivial.
            ! It turns out to be quite easy with bit operations.
            ! The idea is to do a "dumb" permutation where the determinants are
            ! expressed in two sections: orbitals not involved in the excitation
            ! and those that are.  Each section is stored in ascending index
            ! order.
            ! To obtain such ordering requires (for each orbital that is
            ! involved in the excitation) a total of
            ! nel - iel - max_excit + iexcit
            ! where nel is the number of electrons, iel is the position of the
            ! orbital within the list of occupied states in the determinant,
            ! max_excit is the total number of excitations and iexcit is the number
            ! of the "current" orbital involved in excitations.
            ! e.g. Consider (1, 2, 3, 4, 5) -> (1, 3, 5, 6, 7).
            ! (1, 2, 3, 4) goes to (1, 3, 2, 4).
            ! 2 is the first (iexcit=1) orbital found involved in the excitation
            ! and so requires 5 - 2 - 2 + 1 = 2 permutation to shift it to the
            ! first "slot" in the excitation "block" in the list of states.
            ! 4 is the second orbital found and requires 5 - 4 - 2 + 2 = 1
            ! permutation to shift it the end (last "slot" in the excitation
            ! block).
            ! Whilst the resultant number of permutations isn't necessarily the
            ! minimal number for the determinants to align, this is irrelevant
            ! as the Slater--Condon rules only care about whether the number of
            ! permutations are odd or even.

            ! n.b. We don't need to include shift or iexcit in the perm
            !      calculation, as is it symmetric as iexcit reaches the same
            !      maximum value for both src and target iluts
            !shift = nel - max_excit

            do i = 0, NIfD
                if (iLutnI(i) == iLutnJ(i)) cycle
                do j = 0, end_n_int

                    testI = btest(iLutnI(i),j)
                    testJ = btest(iLutnJ(i),j)

                    if (testJ) iel2 = iel2 + 1

                    if (testI) then
                        iel1 = iel1 + 1
                        if (.not.testJ) then
                            ! occupied in iLutnI but not in iLutnJ
                            iexcit1 = iexcit1 + 1
                            Ex(1,iexcit1) = i*bits_n_int+j+1
                            !perm = perm + (shift - iel1 + iexcit1)
                            perm = perm + iel1
                        end if
                    else
                        if (testJ) then
                            ! occupied in iLutnI but not in iLutnJ
                            iexcit2 = iexcit2 + 1
                            Ex(2,iexcit2) = i*bits_n_int+j+1
                            !perm = perm + (shift - iel2 + iexcit2)
                            perm = perm + iel2
                        end if
                    end if
                    if (iexcit1 == max_excit .and. iexcit2 == max_excit) exit
                end do
                if (iexcit1 == max_excit .and. iexcit2 == max_excit) exit
            end do

            ! It seems that this test is faster than btest(perm,0)!
            tSign = mod(perm,2) == 1

            if (iexcit1<max_excit) then
                Ex(:,iexcit1+1) = 0 ! Indicate we've ended the excitation.
            end if

        end if

    end subroutine GetBitExcitation

!This routine will find the largest bit set in a bit-string (i.e. the highest value orbital)
    SUBROUTINE LargestBitSet(iLut,NIfD,LargestOrb)
        use constants, only: bits_n_int,end_n_int,n_int
        IMPLICIT NONE
        INTEGER :: LargestOrb, NIfD,i,j
        INTEGER(KIND=n_int) :: iLut(0:NIfD)

!        do i=NIfD,0,-1
!!Count down through the integers in the bit string.
!!The largest set bit is equal to INT(log_2 (N))
!            IF(iLut(i).ne.0) THEN
!                LargestOrb=NINT(LOG(REAL(iLut(i)+1))*1.4426950408889634)
!                EXIT
!            ENDIF
!        enddo
!        LargestOrb=LargestOrb+(i*32)

        outer: do i=NIfD,0,-1
            do j=end_n_int,0,-1
                IF(BTEST(iLut(i),j)) THEN
                    EXIT outer
                ENDIF
            enddo
        enddo outer
        LargestOrb=(i*bits_n_int)+j+1

    END SUBROUTINE LargestBitSet


!This routine will find the i and a orbitals from a single excitation.
!NOTE! This routine will find i and a, but not distinguish between them. To calculate which one i is,
!you would need to do another XOR with the original orbital and find out which bit this corresponded to.
    SUBROUTINE FindSingleOrbs(iLutnI,iLutnJ,NIfD,Orbs)
        use constants, only: n_int,bits_n_int
        IMPLICIT NONE
        INTEGER :: NIfD,Orbs(2)
        INTEGER(KIND=n_int) :: iLutnI(0:NIfD),iLutnJ(0:NIfD)
        INTEGER(kind=n_int) :: iLutExcited(0:NIfD)

        iLutExcited(:)=IEOR(iLutnI(:),iLutnJ(:))
        CALL LargestBitSet(iLutExcited,NIfD,Orbs(1))
!Found first orbital. Now clear this from the list and search again for the second....
        iLutExcited((Orbs(1)-1)/bits_n_int)=IBCLR(iLutExcited((Orbs(1)-1)/bits_n_int),mod(Orbs(1)-1,bits_n_int))
        CALL LargestBitSet(iLutExcited,NIfD,Orbs(2))

    END SUBROUTINE FindSingleOrbs

