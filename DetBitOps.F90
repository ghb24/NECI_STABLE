!This file contains a load of useful operations to perform on determinants represented as bit-strings.
! Start the process of modularising this bit!!
module DetBitOps
    use Systemdata, only: nel, NIfD, NIfY, NIfTot, tCSF
    implicit none

    ! http://gurmeetsingh.wordpress.com/2008/08/05/fast-bit-counting-routines/
    ! for a variety of interesting bit counters
    interface CountBits
        !module procedure CountBits_sparse
        module procedure CountBits_nifty
    end interface

    contains

    ! This will count the bits set in a bit-string up to a number nBitsMax, if
    ! provided.
    ! The function will return 0 -> nBitsMax+1
    ! A value of nBitsMax+1 indicates that more bits are set than was expected.
    ! The total number of set bits can exceed nBitsMax+1, however.
    ! Counts bits set in integer array (0:nLast)
    integer function CountBits_sparse (iLut, nLast, nBitsMax)
        integer, intent(in), optional :: nBitsMax
        integer, intent(in) :: nLast, iLut(0:nLast)
        integer :: iLutTemp(0:nLast), i, lnBitsMax

        ! By default, allow all the bits to be set
        if (present(nBitsMax)) then
            lnBitsMax = nBitsMax
        else
            lnBitsMax = 32 * (nLast+1)
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
    ! See comments for CountBits_sparse.
    integer function Countbits_nifty (iLut, nLast, nBitsMax)
        integer, intent(in), optional :: nBitsMax
        integer, intent(in) :: nLast, iLut(0:nLast)
        integer ::tmp, i, lnBitsMax

        ! By default, allow all the bits to be set
        if (present(nBitsMax)) then
            lnBitsMax = nBitsMax
        else
            lnBitsMax = 32 * (nLast+1)
        endif

        CountBits_nifty = 0
        do i=0,nLast
            tmp = iLut(i)
            tmp = tmp - iand(ishft(tmp, -1), Z'55555555')
            tmp = iand(tmp, Z'33333333') + iand(ishft(tmp, -2), Z'33333333')
            tmp = iand((tmp+ishft(tmp, -4)), Z'F0F0F0F') * Z'1010101'
            CountBits_nifty = CountBits_nifty + ishft(tmp, -24)

            if (CountBits_nifty .gt. lnBitsMax) then
                CountBits_nifty = lnBitsmax+1
                return
            endif
        enddo       
    end function CountBits_nifty

    integer function count_open_orbs (iLut)
        integer, intent(in) :: iLut(0:NIfD)
        integer, dimension(0:NIfD) :: alpha, beta

        alpha = iand(iLut, Z'AAAAAAAA')
        beta = iand(iLut, Z'55555555')
        alpha = ishft(alpha, -1)
        alpha = ieor(alpha, beta)
        
        count_open_orbs = CountBits(alpha, NIfD)
    end function

    !This will return true if iLutI is identical to iLutJ and will return false otherwise.
    logical function DetBitEQ(iLutI,iLutJ,nLast)
        integer, intent(in), optional :: nLast
        integer, intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer :: i, lnLast

        if(iLutI(0).ne.iLutJ(0)) then
            DetBitEQ=.false.
            return
        else
            if (present(nLast)) then
                lnLast = nLast
            else
                lnLast = NIftot
            endif

            do i=1,lnLast
                if(iLutI(i).ne.iLutJ(i)) then
                    DetBitEQ=.false.
                    return
                endif
            enddo
        endif
        DetBitEQ=.true.
    end function DetBitEQ

    ! This will return 1 if iLutI is "less" than iLutJ, 0 if the determinants
    ! are identical, or -1 if iLutI is "more" than iLutJ
    integer function DetBitLT(iLutI,iLutJ,nLast)
        integer, intent(in), optional :: nLast
        integer, intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
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
                lnLast = NIftot
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
        integer :: iLutI(0:NIfTot),iLutJ(0:NIfTot),i
        integer :: iLut2I(0:NIfTot),iLut2J(0:NIfTot),lnLast

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
                lnLast = NIftot
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
        integer, intent(in) :: iLutI(0:NIftot), iLutJ(0:NIfTot)
        integer, intent(in) :: iLutHF(0:NIfTot)
        integer i, ExcitLevelI, ExcitLevelJ,lnLast
        
        CALL FindBitExcitLevel(iLutI,iLutHF,ExcitLevelI,nel)
        CALL FindBitExcitLevel(iLutJ,iLutHF,ExcitLevelJ,nel)

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
                    lnLast = NIftot
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
    subroutine EncodeBitDet(nI,iLut)
        use csf, only: iscsf, csf_yama_bit, csf_orbital_mask
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: iLut(0:NIfTot)
        integer :: i, det, pos, nopen
        logical :: open_shell

        iLut(:)=0
        nopen = 0
        if (tCSF .and. iscsf (nI)) then
            do i=1,nel
                ! THe first non-paired orbital has yama symbol = 1
                if ((.not. open_shell) .and. &
                    btest(nI(i), csf_yama_bit)) open_shell = .true.

                ! Set the bit in the bit representation
                det = iand(nI(i), csf_orbital_mask)
                iLut((det-1)/32) = ibset(iLut((det-1)/32),mod(det-1,32))

                if (open_shell) then
                    if (btest(nI(i), csf_yama_bit)) then
                        pos = NIfD + 1 + (nopen/32)
                        iLut(pos) = ibset(iLut(pos), mod(nopen,32))
                    endif
                    nopen = nopen + 1
                endif
            enddo
            iLut (NIfD+1) = nopen
        else
            do i=1,nel
                iLut((nI(i)-1)/32)=IBSET(iLut((nI(i)-1)/32),mod(nI(i)-1,32))
            enddo
        endif

    end subroutine EncodeBitDet

    ! This is a routine to take a determinant in bit form and construct 
    ! the natural ordered NEl integer for of the det.
    ! If CSFs are enabled, transefer the yamanouchi symbol as well.
    subroutine DecodeBitDet(nI,iLut)
        use csf, only: csf_yama_bit, csf_test_bit
        integer, intent(in) :: iLut(0:NIfTot)
        integer, intent(out) :: nI(nel)
        integer :: i, j, elec, pos, nopen

        elec=0
        if (tCSF) then
            ! Consider the closed shell electrons first
            do i=0,NIfD
                do j=0,30,2
                    if (btest(iLut(i),j)) then
                        if (btest(iLut(i),j+1)) then
                            ! An electron pair is in this spacial orbital
                            ! (2 matched spin orbitals)
                            elec = elec + 2
                            nI(elec-1) = (32*i) + (j+1)
                            nI(elec) = (32*i) + (j+2)
                            if (elec == nel) return
                        endif
                    endif
                enddo
            enddo

            ! Now consider the open shell electrons
            ! TODO: can we move in steps of two, to catch unmatched pairs?
            nopen = 0
            do i=0,NIfD
                do j=0,31
                    if (btest(iLut(i),j)) then
                        if (.not.btest(iLut(i),xor(j,1))) then
                            elec = elec + 1
                            nI(elec) = (32*i) + (j+1)
                            pos = NIfD + 1 + (nopen/32)
                            if (btest(iLut(Pos),mod(nopen,32))) then
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
                do j=0,31
                    if(btest(iLut(i),j)) then
                        !An electron is at this orbital
                        elec=elec+1
                        nI(elec)=(i*32)+(j+1)
                        if(elec.eq.nel) return
                    endif
                enddo
            enddo
        endif
    end subroutine DecodeBitDet

    subroutine shift_det_bit_singles_to_beta (iLut)
        integer, intent(inout) :: iLut(0:NIfD)
        integer :: iA(0:NIfD), iB(0:NIfD)

        ! Extract the betas
        iB = iand(iLut, Z'55555555')

        ! Extract the alphas and shift them into beta positions.
        iA = ishft(iand(iLut, Z'AAAAAAAA'), -1)
        
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
        integer, dimension(0:NIfTot) :: alpha, beta, tmp
        integer :: first_beta_byte, first_beta_bit
        integer i

        call EncodeBitDet(nI, alpha)
        beta = iand(alpha, Z'55555555')
        tmp = iand(alpha, Z'AAAAAAAA')

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
            do i=0,31
                if (btest(beta(first_beta_byte),i)) exit
            enddo
            first_beta_bit = i

            ! TODO: steps of 2, as alpha/beta even/odd...
            do i=31,first_beta_bit,-1
                if (btest(alpha(first_beta_byte), i)) exit
            enddo
            if (i < first_beta_bit) is_canonical_ms_order = .true.
        endif
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

        use SystemData, only: NIfD, nel
        use DetBitOps, only: CountBits_nifty
        implicit none
        integer, intent(in) :: iLutnI(0:NIfD), iLutnJ(0:NIfD)
        integer, intent(inout) :: Ex(2,*)
        logical, intent(out) :: tSign
        integer :: i, j, iexcit1, iexcit2, perm, iel1, iel2, shift, max_excit
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
            shift = nel - max_excit

            do i = 0, NIfD
                if (iLutnI(i) == iLutnJ(i)) cycle
                do j = 0, 31

                    testI = btest(iLutnI(i),j)
                    testJ = btest(iLutnJ(i),j)

                    if (testJ) iel2 = iel2 + 1

                    if (testI) then
                        iel1 = iel1 + 1
                        if (.not.testJ) then
                            ! occupied in iLutnI but not in iLutnJ
                            iexcit1 = iexcit1 + 1
                            Ex(1,iexcit1) = i*32+j+1
                            perm = perm + (shift - iel1 + iexcit1)
                        end if
                    else
                        if (testJ) then
                            ! occupied in iLutnI but not in iLutnJ
                            iexcit2 = iexcit2 + 1
                            Ex(2,iexcit2) = i*32+j+1
                            perm = perm + (shift - iel2 + iexcit2)
                        end if
                    end if

                end do
            end do

            ! It seems that this test is faster than btest(perm,0)!
            tSign = mod(perm,2) == 1

            if (iexcit1<max_excit) then
                Ex(:,iexcit1+1) = 0 ! Indicate we've ended the excitation.
            end if

        end if

    end subroutine GetBitExcitation


!This routine will find the bit-representation of an excitation by constructing the new iLut from the old one and the excitation matrix.
    SUBROUTINE FindExcitBitDet(iLutnI,iLutnJ,IC,ExcitMat,NIfD)
        IMPLICIT NONE
        INTEGER :: iLutnI(0:NIfD),iLutnJ(0:NIfD),IC,ExcitMat(2,2),NIfD

        iLutnJ(:)=iLutnI(:)
        IF(IC.eq.1) THEN
!Single excitation - clear one bit and set another.
            iLutnJ((ExcitMat(1,1)-1)/32)=IBCLR(iLutnJ((ExcitMat(1,1)-1)/32),mod(ExcitMat(1,1)-1,32))
            iLutnJ((ExcitMat(2,1)-1)/32)=IBSET(iLutnJ((ExcitMat(2,1)-1)/32),mod(ExcitMat(2,1)-1,32))
        ELSE
!Double excitation - clear two bits and set two others.
            iLutnJ((ExcitMat(1,1)-1)/32)=IBCLR(iLutnJ((ExcitMat(1,1)-1)/32),mod(ExcitMat(1,1)-1,32))
            iLutnJ((ExcitMat(2,1)-1)/32)=IBSET(iLutnJ((ExcitMat(2,1)-1)/32),mod(ExcitMat(2,1)-1,32))
            iLutnJ((ExcitMat(1,2)-1)/32)=IBCLR(iLutnJ((ExcitMat(1,2)-1)/32),mod(ExcitMat(1,2)-1,32))
            iLutnJ((ExcitMat(2,2)-1)/32)=IBSET(iLutnJ((ExcitMat(2,2)-1)/32),mod(ExcitMat(2,2)-1,32))
        ENDIF

    END SUBROUTINE FindExcitBitDet

!This function will return true if the determinant is closed shell, or false if not.
    LOGICAL FUNCTION TestClosedShellDet(iLut)
        use systemdata, only: NIfD
        INTEGER :: iLut(0:NIfD),iLutAlpha(0:NIfD),iLutBeta(0:NIfD),MaskAlpha,MaskBeta,i
        
        iLutAlpha(:)=0
        iLutBeta(:)=0
        MaskBeta=1431655765    !This is 1010101... in binary
        MaskAlpha=-1431655766  !This is 0101010... in binary
        TestClosedShellDet=.true.

        do i=0,NIfD

            iLutAlpha(i)=IAND(iLut(i),MaskAlpha)    !Seperate the alpha and beta bit strings
            iLutBeta(i)=IAND(iLut(i),MaskBeta)
            iLutAlpha(i)=ISHFT(iLutAlpha(i),-1)     !Shift all alpha bits to the left by one.
            iLutAlpha(i)=IEOR(iLutAlpha(i),iLutBeta(i)) !Do an XOR on the original beta bits and shifted alpha bits - they should cancel exactly.
            
            IF(iLutAlpha(i).ne.0) THEN
                TestClosedShellDet=.false.  !Det is not closed shell - return
                RETURN
            ENDIF
        enddo

    END FUNCTION TestClosedShellDet

!Routine to count number of open *SPATIAL* orbitals in a bit-string representation of a determinant.
! ************************
! BROKEN
! NOTE: This function name is misleading
!       It counts the number of unpaired Beta electrons (ignores Alpha)
!       --> Returns nopen/2 <==> Ms=0
! ************************
    SUBROUTINE CalcOpenOrbs(iLut,OpenOrbs)
        use systemdata, only: NIfD, nel
        use DetBitOps, only: CountBits
        INTEGER :: iLut(0:NIfD),iLutAlpha(0:NIfD),iLutBeta(0:NIfD),MaskAlpha,MaskBeta,i,OpenOrbs
        
        iLutAlpha(:)=0
        iLutBeta(:)=0
        MaskBeta=1431655765    !This is 1010101... in binary
        MaskAlpha=-1431655766  !This is 0101010... in binary

!        do i=0,NIfD
!
!            iLutAlpha(i)=IAND(iLut(i),MaskAlpha)    !Seperate the alpha and beta bit strings
!            iLutBeta(i)=IAND(iLut(i),MaskBeta)
!            iLutAlpha(i)=ISHFT(iLutAlpha(i),-1)     !Shift all alpha bits to the left by one.
!            iLutAlpha(i)=IEOR(iLutAlpha(i),iLutBeta(i)) !Do an XOR on the original beta bits and shifted alpha bits - only open shell occupied orbitals will remain.
!            
!        enddo
!
!        OpenOrbs = CountBits(iLutAlpha,NIfD,NEl)
!        OpenOrbs=OpenOrbs/2

!Alternatively....use a NOT and an AND to only count half as many set bits

        do i=0,NIfD     
                    
            iLutAlpha(i)=IAND(iLut(i),MaskAlpha)    !Seperate the alpha and beta bit strings
            iLutBeta(i)=IAND(iLut(i),MaskBeta)
            iLutAlpha(i)=ISHFT(iLutAlpha(i),-1)     !Shift all alpha bits to the left by one.

            iLutAlpha(i)=NOT(iLutAlpha(i))              ! This NOT means that set bits are now represented by 0s, not 1s
            iLutAlpha(i)=IAND(iLutAlpha(i),iLutBeta(i)) ! Now, only the 1s in the beta string will be counted.

        enddo

        OpenOrbs = CountBits(iLutAlpha,NIfD,NEl)
    END SUBROUTINE CalcOpenOrbs

!This routine will find the largest bit set in a bit-string (i.e. the highest value orbital)
    SUBROUTINE LargestBitSet(iLut,NIfD,LargestOrb)
        IMPLICIT NONE
        INTEGER :: LargestOrb, iLut(0:NIfD),NIfD,i,j

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
            do j=31,0,-1
                IF(BTEST(iLut(i),j)) THEN
                    EXIT outer
                ENDIF
            enddo
        enddo outer
        LargestOrb=(i*32)+j+1

    END SUBROUTINE LargestBitSet

!This routine will find the excitation level of two determinants in bit strings.
    SUBROUTINE FindBitExcitLevel(iLutnI,iLutnJ,ExcitLevel,MaxExcitLevel)
        use SystemData, only: NIfD
        use DetBitOps, only: CountBits
        IMPLICIT NONE
        INTEGER :: iLutnI(0:NIfD),iLutnJ(0:NIfD),ExcitLevel,MaxExcitLevel
        INTEGER :: iLutExcited(0:NIfD),i,k

!First find a new bit string which just contains the excited orbitals
        iLutExcited(:)=IEOR(iLutnI(:),iLutnJ(:))
        iLutExcited(:)=IAND(iLutExcited(:),iLutnI(:))

!Now, simply count the bits in it...
        ExcitLevel = CountBits(iLutExcited, NIfD, MaxExcitLevel)

    END SUBROUTINE FindBitExcitLevel

!This routine will find the i and a orbitals from a single excitation.
!NOTE! This routine will find i and a, but not distinguish between them. To calculate which one i is,
!you would need to do another XOR with the original orbital and find out which bit this corresponded to.
    SUBROUTINE FindSingleOrbs(iLutnI,iLutnJ,NIfD,Orbs)
        IMPLICIT NONE
        INTEGER :: iLutnI(0:NIfD),iLutnJ(0:NIfD),NIfD,Orbs(2)
        INTEGER :: iLutExcited(0:NIfD)

        iLutExcited(:)=IEOR(iLutnI(:),iLutnJ(:))
        CALL LargestBitSet(iLutExcited,NIfD,Orbs(1))
!Found first orbital. Now clear this from the list and search again for the second....
        iLutExcited((Orbs(1)-1)/32)=IBCLR(iLutExcited((Orbs(1)-1)/32),mod(Orbs(1)-1,32))
        CALL LargestBitSet(iLutExcited,NIfD,Orbs(2))

    END SUBROUTINE FindSingleOrbs

    
! Based on SORTI, SortBitDets sorts determinants as bit strings, and takes the corresponding element from array RB with it (sign)
! RA is the array of determinants of length N to sort
! The RA array elements go from 0:NIfD
! RB is the array of integers to go with the determinant
      SUBROUTINE SortBitDets(N,RA,RB)
      use SystemData, only: NIfTot,NIfDBO
      use DetBitOps, only: DetBitLT
      INTEGER N,I,L,IR,J
      INTEGER RA(0:NIfTot,N)
      INTEGER RB(N)
      INTEGER RRA(0:NIfTot),RRB
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA(:)=RA(:,L)
          RRB=RB(L)
        ELSE
          RRA(:)=RA(:,IR)
          RA(:,IR)=RA(:,1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(:,1)=RRA(:)
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF((DetBitLT(RA(:,J),RA(:,J+1),NIfDBO)).eq.1) J=J+1
          ENDIF
          IF((DetBitLT(RRA(:),RA(:,J),NIfDBO)).eq.1) THEN
            RA(:,I)=RA(:,J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(:,I)=RRA(:)
        RB(I)=RRB

      GO TO 10

      END SUBROUTINE SortBitDets
     
! Based on SortBitDets, this routine sorts determinants firstly by array RA as bit strings, then by the bit strings in array RA2.
! When ordered, the determinants take the corresponding elements from array RB with them (sign).
! RA is the array of determinants of length N to sort first.
! The RA array elements go from 0:NIfD
! RA2 is the array of determinants of length N2 to sort second.
! the RA2 array elements go from 0:NIfD2
! RB is the array of integers to go with the determinant

!        CALL Sort2BitDetsPlus3(MinorValidSpawned,MinorSpawnDets(0:NoIntforDet,1:MinorValidSpawned),NoIntforDet,MinorSpawnParent(0:NoIntforDet,1:MinorValidSpawned),&
!                &NoIntforDet,MinorSpawnSign(MinorValidSpawned))

      SUBROUTINE Sort2BitDetsPlus3(N,RA,RA2,RB)
      use DetBitOps, only: Det2BitLT
      use SystemData, only: NIfTot,NIfDBO
      INTEGER N,I,L,IR,J
      INTEGER RA(0:NIfTot,N),RA2(0:NIfTot,N)
      INTEGER RB(N),RC(N),RD(N)
      INTEGER RRA(0:NIfTot),RRA2(0:NIfTot),RRB
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA(:)=RA(:,L)
          RRA2(:)=RA2(:,L)
          RRB=RB(L)
        ELSE
          RRA(:)=RA(:,IR)
          RA(:,IR)=RA(:,1)
          RRA2(:)=RA2(:,IR)
          RA2(:,IR)=RA2(:,1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(:,1)=RRA(:)
            RA2(:,1)=RRA2(:)
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(Det2BitLT(RA(:,J),RA(:,J+1),RA2(:,J),RA2(:,J+1),NIfDBO).eq.1) J=J+1
          ENDIF
          IF((Det2BitLT(RRA(:),RA(:,J),RRA2(:),RA2(:,J),NIfDBO)).eq.1) THEN
            RA(:,I)=RA(:,J)
            RA2(:,I)=RA2(:,J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(:,I)=RRA(:)
        RA2(:,I)=RRA2(:)
        RB(I)=RRB

      GO TO 10

      END SUBROUTINE Sort2BitDetsPlus3
    
 
! Based on SortBitDets, this sorts determinants as bit strings by excitation level and then determinant and takes the corresponding
! element from array RB with it (sign).
! RA is the array of determinants of length N to sort
! The RA array elements go from 0:NIfD
! RB is the array of integers to go with the determinant
! iLutHF is the HF determinant (in bit string), and NEl is the number of electrons.
      SUBROUTINE SortExcitBitDets(N,RA,RB,iLutHF)
          use DetBitOps, only: DetExcitBitLT
          use SystemData, only: NIftot, nel, NIfDBO
      INTEGER N,I,L,IR,J
      INTEGER RA(0:NIfTot,N)
      INTEGER RB(N)
      INTEGER iLutHF(0:NIfTot)
      INTEGER RRA(0:NIfTot),RRB
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA(:)=RA(:,L)
          RRB=RB(L)
        ELSE
          RRA(:)=RA(:,IR)
          RA(:,IR)=RA(:,1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(:,1)=RRA(:)
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF((DetExcitBitLT(RA(:,J),RA(:,J+1),iLutHF(:),NIfDBO)).eq.1) J=J+1
          ENDIF
          IF((DetExcitBitLT(RRA(:),RA(:,J),iLutHF(:),NIfDBO)).eq.1) THEN
            RA(:,I)=RA(:,J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(:,I)=RRA(:)
        RB(I)=RRB

      GO TO 10

      END SUBROUTINE SortExcitBitDets
    

! Based on SortBitDets, however in SortBitSign, RA is an array of integers (sign) to be sorted in descending order of absolute size, and
! RB is an array of elements going from 0:NIfD (determinants) to be taken with the element of RA.
! RA has length N.
      SUBROUTINE SortBitSign(N,RA,RB)
        use Systemdata, only: NIfTot
      INTEGER N,I,L,IR,J
      INTEGER RA(N)
      INTEGER RB(0:NIfTot,N)
      INTEGER RRA,RRB(0:NIfTot)
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB(:)=RB(:,L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          RRB(:)=RB(:,IR)
          RB(:,IR)=RB(:,1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(:,1)=RRB(:)
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ABS(RA(J)).gt.ABS(RA(J+1))) J=J+1
          ENDIF
          IF(ABS(RRA).gt.ABS(RA(J))) THEN
            RA(I)=RA(J)
            RB(:,I)=RB(:,J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(:,I)=RRB(:)

      GO TO 10

      END SUBROUTINE SortBitSign
   

! Based on SORTI, SortBitDets sorts determinants as bit strings, and takes the corresponding element from array RB with it (sign) and the real array RC (HElems)
! RA is the array of determinants of length N to sort
! The RA array elements go from 0:NIfD
! RB is the array of integers to go with the determinant
      SUBROUTINE SortBitDetswH(N,RA,RB,RC)
        use systemdata, only: NIfTot, nel, NIfDBO
        use DetBitOps, only: DetBitLT
      INTEGER N,I,L,IR,J
      INTEGER RA(0:NIfTot,N)
      INTEGER RB(N)
      REAL*8 RC(N),RRC
      INTEGER RRA(0:NIfTot),RRB
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA(:)=RA(:,L)
          RRB=RB(L)
          RRC=RC(L)
        ELSE
          RRA(:)=RA(:,IR)
          RA(:,IR)=RA(:,1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          RRC=RC(IR)
          RC(IR)=RC(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(:,1)=RRA(:)
            RB(1)=RRB
            RC(1)=RRC
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF((DetBitLT(RA(:,J),RA(:,J+1),NIfDBO)).eq.1) J=J+1
          ENDIF
          IF((DetBitLT(RRA(:),RA(:,J),NIfDBO)).eq.1) THEN
            RA(:,I)=RA(:,J)
            RB(I)=RB(J)
            RC(I)=RC(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(:,I)=RRA(:)
        RB(I)=RRB
        RC(I)=RRC

      GO TO 10

      END SUBROUTINE SortBitDetswH




