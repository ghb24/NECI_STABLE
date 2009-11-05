!This file contains a load of useful operations to perform on determinants represented as bit-strings.

!This will return true if iLutI is identical to iLutJ and will return false otherwise.
    LOGICAL FUNCTION DetBitEQ(iLutI,iLutJ,NIfD)
        IMPLICIT NONE
        INTEGER :: iLutI(0:NIfD),iLutJ(0:NIfD),NIfD,i

        IF(iLutI(0).ne.iLutJ(0)) THEN
            DetBitEQ=.false.
            RETURN
        ELSE
            do i=1,NIfD
                IF(iLutI(i).ne.iLutJ(i)) THEN
                    DetBitEQ=.false.
                    RETURN
                ENDIF
            enddo
        ENDIF
        DetBitEQ=.true.

    END FUNCTION DetBitEQ

!if(Det2BitLT(MinorStarDets(0:NIfD,nup),DetCurr(0:NIfD),NIfD,MinorStarParent(0:NIfD,nup),DetCurr2(0:NIfD),NIfD).eq.1) then 

!This will return 1 if iLutI is "less" than iLutJ, or -1 if iLutI is "more" than iLutJ.  If these are identical, this routine looks at 
!iLut2I and iLut2J, and returns 1 if iLut2I is "less" than iLut2J, -1 if iLut2I is "more than iLut2J, and 0 if these are still identical.
    INTEGER FUNCTION Det2BitLT(iLutI,iLutJ,NIfD,iLut2I,iLut2J,NIfD2)
        IMPLICIT NONE
        INTEGER :: iLutI(0:NIfD),iLutJ(0:NIfD),NIfD,i
        INTEGER :: iLut2I(0:NIfD2),iLut2J(0:NIfD2),NIfD2

        IF(iLutI(0).lt.iLutJ(0)) THEN
!First, compare first integers
            Det2BitLT=1
            RETURN
        ELSEIF(iLutI(0).gt.iLutJ(0)) THEN
            Det2BitLT=-1
            RETURN
        ELSEIF(iLutI(0).eq.iLutJ(0)) THEN
!If the integers are the same, then cycle through the rest of the integers until we find a difference.
            do i=1,NIfD
                IF(iLutI(i).lt.iLutJ(i)) THEN
                    Det2BitLT=1
                    RETURN
                ELSEIF(iLutI(i).gt.iLutJ(i)) THEN
                    Det2BitLT=-1
                    RETURN
                ENDIF
            enddo
!If we get through this loop without RETURN-ing, iLutI and iLutJ are identical, so look to iLut2I and iLut2J            
            IF(iLut2I(0).lt.iLut2J(0)) THEN
                Det2BitLT=1
                RETURN
            ELSEIF(iLut2I(0).gt.iLut2J(0)) THEN
                Det2BitLT=-1
                RETURN
            ELSEIF(iLut2I(0).eq.iLut2J(0)) THEN
                do i=1,NIfD2
                    IF(iLut2I(i).lt.iLut2J(i)) THEN
                        Det2BitLT=1
                        RETURN
                    ELSEIF(iLut2I(i).gt.iLut2J(i)) THEN
                        Det2BitLT=-1
                        RETURN
                    ENDIF
                enddo
            ENDIF
!If we still have not returned, both determinants are identical.               
        ENDIF
        Det2BitLT=0

    END FUNCTION Det2BitLT


!This will return 1 if iLutI is "less" than iLutJ, 0 if the determinants are identical, or -1 if iLutI is "more" than iLutJ
    INTEGER FUNCTION DetBitLT(iLutI,iLutJ,NIfD)
        IMPLICIT NONE
        INTEGER :: iLutI(0:NIfD),iLutJ(0:NIfD),NIfD,i

        IF(iLutI(0).lt.iLutJ(0)) THEN
!First, compare first integers
            DetBitLT=1
        ELSEIF(iLutI(0).eq.iLutJ(0)) THEN
!If the integers are the same, then cycle through the rest of the integers until we find a difference.
            do i=1,NIfD
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


!This will return 1 if iLutI is "less" than iLutJ, 0 if the determinants are identical, or -1 if iLutI is "more" than iLutJ
!This particular version checks excitation level initially, then only if these are the same does it move on to determinants.
    INTEGER FUNCTION DetExcitBitLT(iLutI,iLutJ,iLutHF,NIfD,NEl)
        IMPLICIT NONE
        INTEGER :: iLutI(0:NIfD),iLutJ(0:NIfD),iLutHF(0:NIfD),NIfD,i,ExcitLevelI,ExcitLevelJ,NEl
        
        CALL FindBitExcitLevel(iLutI,iLutHF,NIfD,ExcitLevelI,NEl)
        CALL FindBitExcitLevel(iLutJ,iLutHF,NIfD,ExcitLevelJ,NEl)

! First order in terms of excitation level.  I.e. if the excitation levels are different, we don't care what the determinants
! are we just order in terms of the excitation level.
        IF(ExcitLevelI.lt.ExcitLevelJ) THEN
            DetExcitBitLT=1
            RETURN
        ELSEIF(ExcitLevelI.gt.ExcitLevelJ) THEN
            DetExcitBitLT=-1
            RETURN

! If the excitation levels are the same however, we need to look at the determinant and order according to this.            
        ELSEIF(ExcitLevelI.eq.ExcitLevelJ) THEN
        
            IF(iLutI(0).lt.iLutJ(0)) THEN
! First, compare first integers
                DetExcitBitLT=1
                RETURN
            ELSEIF(iLutI(0).eq.iLutJ(0)) THEN
! If the integers are the same, then cycle through the rest of the integers until we find a difference.
                do i=1,NIfD
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
! If it gets through all this without being returned then the two determinants are equal and DetExcitBitLT=0
            DetExcitBitLT=0

        ENDIF

    END FUNCTION DetExcitBitLT

! Start the process of modularising this bit!!
module DetBitOps
    use Systemdata, only: nel, NIfD, NIfY, NIfTot, tCSF
    implicit none
    contains
    ! This is a routine to encode a determinant as natural ordered integers
    ! (nI) as a bit string (iLut(0:NIfTot)) where NIfD=INT(nBasis/32)
    ! If this is a csf, the csf is contained afterwards.
    subroutine EncodeBitDet(nI,iLut)
        use csf, only: iscsf, csf_yama_bit, csf_orbital_mask
        implicit none
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: iLut(0:NIfTot)
        integer :: i, det, pos, nopen
        logical :: open_shell

        iLut(:)=0
        nopen = 0
        if (tCSF .and. iscsf (nI)) then
            do i=1,nel
                ! The first non-paired orbital has yama symbol = 1
                if ((.not. open_shell) .and. &
                    btest(nI(i), csf_yama_bit)) open_shell = .true.

                ! Set the bit in the bit representation
                det = iand(nI(i), csf_orbital_mask)
                iLut((det-1)/32) = ibset(iLut((det-1)/32),mod(det-1,32))

                if (open_shell) then
                    if (btest(nI(i), csf_yama_bit)) then
                        pos = NIfD + 2 + (nopen/32)
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
        use SystemData, only: nel, NIfD, NIfY, NIfTot, tCSF
        use csf, only: csf_yama_bit, csf_test_bit
        implicit none
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
                            pos = NIfD + 2 + (nopen/32)
                            if (btest(iLut(Pos),mod(nopen,32))) then
                                nI(elec) = ibset(nI(elec),csf_yama_bit)
                            endif
                            nopen = nopen + 1
                        endif
                    endif
                    if ((elec==nel) .or. (nopen==iLut(NIfD+1))) exit
                enddo
                if ((elec==nel) .or. (nopen==iLut(NIfD+1))) exit
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
end module

    SUBROUTINE GetBitExcitation(iLutnI,iLutnJ,NIfD,NEl,Ex,tSign)
        IMPLICIT NONE
        INTEGER :: iLutnI(0:NIfD),iLutnJ(0:NIfD),Ex(2,*),i,BitExcitMat(0:NIfD),BitCommonOrbs(0:NIfD),FromInd,ToInd,j
        INTEGER :: nIElec,nJElec,Diff,NIfD,NEl,iMaxSize
        LOGICAL :: tSign

        tSign=.false.
        ToInd=1
        FromInd=1
        nIElec=0
        nJElec=0
        iMaxSize=EX(1,1) 

        do i=0,NIfD
            BitExcitMat(i)=IEOR(iLutnI(i),iLutnJ(i))
            BitCommonOrbs(i)=IAND(iLutnI(i),iLutnJ(i))
            
            IF((BitExcitMat(i).eq.0).and.(BitCommonOrbs(i).eq.0)) CYCLE

!            WRITE(6,*) "BitExcitMat: "
!            do j=0,31
!                IF(BTEST(BitExcitMat(i),j)) THEN
!                    WRITE(6,"(A)",advance='no') " 1 "
!                ELSE
!                    WRITE(6,"(A)",advance='no') " 0 "
!                ENDIF
!            enddo
!            WRITE(6,*) "*** iLutnI: "
!            do j=0,31
!                IF(BTEST(iLutnI(i),j)) THEN
!                    WRITE(6,"(A)",advance='no') " 1 "
!                ELSE
!                    WRITE(6,"(A)",advance='no') " 0 "
!                ENDIF
!            enddo


            do j=0,31
!First, find the parity of the excitations.
!we need to search the common orbitals, in order to check the parity shift.

                IF(BTEST(BitCommonOrbs(i),j)) THEN
                    IF(nIElec.ne.nJElec) THEN
                        Diff=abs(nIElec-nJElec)
                        IF(mod(Diff,2).eq.1) THEN
                            tSign=.not.tSign
                        ENDIF
                    ENDIF
                    BitCommonOrbs(i)=IBCLR(BitCommonOrbs(i),j)
!                ELSE
!                    IF(BTEST(iLutnI(i),j)) nIElec=nIElec+1
!                    IF(BTEST(iLutnJ(i),j)) nJElec=nJElec+1
                ENDIF

!Now, create the Ex matrix
                IF(BTEST(BitExcitMat(i),j)) THEN
!This is either a 'from' orbital, or a 'to' orbital - find out which.
                    IF(BTEST(iLutnI(i),j)) THEN
!This is one of the orbitals that an electron is excited *from*
                        nIElec=nIElec+1
                        Ex(1,FromInd)=(i*32)+(j+1)
                        FromInd=FromInd+1
                        BitExcitMat(i)=IBCLR(BitExcitMat(i),j)  !Remove the bit
                        IF((BitExcitMat(i).eq.0).and.(BitCommonOrbs(i).eq.0)) EXIT
                    ELSE
!This is one of the orbitals that an electron is excited *to*
                        nJElec=nJElec+1
                        Ex(2,ToInd)=(i*32)+(j+1)
                        ToInd=ToInd+1
                        BitExcitMat(i)=IBCLR(BitExcitMat(i),j)  !Remove the bit
                        IF((BitExcitMat(i).eq.0).and.(BitCommonOrbs(i).eq.0)) EXIT
                    ENDIF
                ENDIF

            enddo
        enddo
        if (ToInd.LE.iMaxSize) THEN
           EX(1,ToInd)=0  !Indicate that we've ended the excit
           EX(2,ToInd)=0  !Indicate that we've ended the excit
        ENDIF
!        IF(ToInd.ne.FromInd) CALL Stop_All("GetBitExcitation","Error in constructing excitation matrix")

    END SUBROUTINE GetBitExcitation


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
    LOGICAL FUNCTION TestClosedShellDet(iLut,NIfD)
        INTEGER :: iLut(0:NIfD),iLutAlpha(0:NIfD),iLutBeta(0:NIfD),MaskAlpha,MaskBeta,i,NIfD
        
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
    SUBROUTINE CalcOpenOrbs(iLut,NIfD,NEl,OpenOrbs)
        INTEGER :: iLut(0:NIfD),iLutAlpha(0:NIfD),iLutBeta(0:NIfD),MaskAlpha,MaskBeta,i,NIfD,NEl,OpenOrbs
        
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
!        CALL CountBits(iLutAlpha,NIfD,OpenOrbs,NEl)
!        OpenOrbs=OpenOrbs/2

!Alternatively....use a NOT and an AND to only count half as many set bits

        do i=0,NIfD     
                    
            iLutAlpha(i)=IAND(iLut(i),MaskAlpha)    !Seperate the alpha and beta bit strings
            iLutBeta(i)=IAND(iLut(i),MaskBeta)
            iLutAlpha(i)=ISHFT(iLutAlpha(i),-1)     !Shift all alpha bits to the left by one.

            iLutAlpha(i)=NOT(iLutAlpha(i))              ! This NOT means that set bits are now represented by 0s, not 1s
            iLutAlpha(i)=IAND(iLutAlpha(i),iLutBeta(i)) ! Now, only the 1s in the beta string will be counted.

        enddo

        CALL CountBits(iLutAlpha,NIfD,OpenOrbs,NEl)


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

!This will count the bits set in a bit-string up to a number nBitsMax.
!NBits will return 0 -> nBitsMax+1
    SUBROUTINE CountBits(iLut,NIfD,nBits,nBitsMax)
        IMPLICIT NONE
        INTEGER :: iLut(0:NIfD),NIfD,i,nBitsMax,nBits,iLutTemp(0:NIfD)

        nBits=0
        iLutTemp(:)=iLut(:)
        do i=0,NIfD
            do while((iLutTemp(i).ne.0).and.(nBits.le.nBitsMax))
                iLutTemp(i)=IAND(iLutTemp(i),iLutTemp(i)-1) !This clears the rightmost set bit
                nBits=nBits+1
            enddo
            IF(nBits.gt.nBitsMax) RETURN
        enddo

    END SUBROUTINE CountBits

!This routine will find the excitation level of two determinants in bit strings.
    SUBROUTINE FindBitExcitLevel(iLutnI,iLutnJ,NIfD,ExcitLevel,MaxExcitLevel)
        IMPLICIT NONE
        INTEGER :: iLutnI(0:NIfD),iLutnJ(0:NIfD),ExcitLevel,NIfD,MaxExcitLevel
        INTEGER :: iLutExcited(0:NIfD),i,k

!First find a new bit string which just contains the excited orbitals
        iLutExcited(:)=IEOR(iLutnI(:),iLutnJ(:))
        iLutExcited(:)=IAND(iLutExcited(:),iLutnI(:))

!Now, simply count the bits in it...
        CALL CountBits(iLutExcited,NIfD,ExcitLevel,MaxExcitLevel)

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
      SUBROUTINE SortBitDets(N,RA,NIfD,RB)
      INTEGER N,NIfD,I,L,IR,J
      INTEGER RA(0:NIfD,N)
      INTEGER RB(N)
      INTEGER RRA(0:NIfD),RRB
      INTEGER DetBitLT
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA(0:NIfD)=RA(0:NIfD,L)
          RRB=RB(L)
        ELSE
          RRA(0:NIfD)=RA(0:NIfD,IR)
          RA(0:NIfD,IR)=RA(0:NIfD,1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(0:NIfD,1)=RRA(0:NIfD)
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF((DetBitLT(RA(0:NIfD,J),RA(0:NIfD,J+1),NIfD)).eq.1) J=J+1
          ENDIF
          IF((DetBitLT(RRA(0:NIfD),RA(0:NIfD,J),NIfD)).eq.1) THEN
            RA(0:NIfD,I)=RA(0:NIfD,J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(0:NIfD,I)=RRA(0:NIfD)
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

      SUBROUTINE Sort2BitDetsPlus3(N,RA,NIfD,RA2,NIfD2,RB)
      INTEGER N,NIfD,NIfD2,I,L,IR,J
      INTEGER RA(0:NIfD,N),RA2(0:NIfD2,N)
      INTEGER RB(N),RC(N),RD(N)
      INTEGER RRA(0:NIfD),RRA2(0:NIfD2),RRB
      INTEGER Det2BitLT
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA(0:NIfD)=RA(0:NIfD,L)
          RRA2(0:NIfD2)=RA2(0:NIfD2,L)
          RRB=RB(L)
        ELSE
          RRA(0:NIfD)=RA(0:NIfD,IR)
          RA(0:NIfD,IR)=RA(0:NIfD,1)
          RRA2(0:NIfD2)=RA2(0:NIfD2,IR)
          RA2(0:NIfD2,IR)=RA2(0:NIfD2,1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(0:NIfD,1)=RRA(0:NIfD)
            RA2(0:NIfD2,1)=RRA2(0:NIfD2)
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(Det2BitLT(RA(0:NIfD,J),RA(0:NIfD,J+1),NIfD,RA2(0:NIfD2,J),RA2(0:NIfD2,J+1),NIfD2).eq.1) J=J+1
          ENDIF
          IF((Det2BitLT(RRA(0:NIfD),RA(0:NIfD,J),NIfD,RRA2(0:NIfD2),RA2(0:NIfD2,J),NIfD2)).eq.1) THEN
            RA(0:NIfD,I)=RA(0:NIfD,J)
            RA2(0:NIfD2,I)=RA2(0:NIfD2,J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(0:NIfD,I)=RRA(0:NIfD)
        RA2(0:NIfD2,I)=RRA2(0:NIfD2)
        RB(I)=RRB

      GO TO 10

      END SUBROUTINE Sort2BitDetsPlus3
    
 
! Based on SortBitDets, this sorts determinants as bit strings by excitation level and then determinant and takes the corresponding
! element from array RB with it (sign).
! RA is the array of determinants of length N to sort
! The RA array elements go from 0:NIfD
! RB is the array of integers to go with the determinant
! iLutHF is the HF determinant (in bit string), and NEl is the number of electrons.
      SUBROUTINE SortExcitBitDets(N,RA,NIfD,RB,iLutHF,NEl)
      INTEGER N,NIfD,I,L,IR,J,NEl
      INTEGER RA(0:NIfD,N)
      INTEGER RB(N)
      INTEGER iLutHF(0:NIfD)
      INTEGER RRA(0:NIfD),RRB
      INTEGER DetExcitBitLT
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA(0:NIfD)=RA(0:NIfD,L)
          RRB=RB(L)
        ELSE
          RRA(0:NIfD)=RA(0:NIfD,IR)
          RA(0:NIfD,IR)=RA(0:NIfD,1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(0:NIfD,1)=RRA(0:NIfD)
            RB(1)=RRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF((DetExcitBitLT(RA(0:NIfD,J),RA(0:NIfD,J+1),iLutHF(0:NIfD),NIfD,NEl)).eq.1) J=J+1
          ENDIF
          IF((DetExcitBitLT(RRA(0:NIfD),RA(0:NIfD,J),iLutHF(0:NIfD),NIfD,NEl)).eq.1) THEN
            RA(0:NIfD,I)=RA(0:NIfD,J)
            RB(I)=RB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(0:NIfD,I)=RRA(0:NIfD)
        RB(I)=RRB

      GO TO 10

      END SUBROUTINE SortExcitBitDets
    

! Based on SortBitDets, however in SortBitSign, RA is an array of integers (sign) to be sorted in descending order of absolute size, and
! RB is an array of elements going from 0:NIfD (determinants) to be taken with the element of RA.
! RA has length N.
      SUBROUTINE SortBitSign(N,RA,NIfD,RB)
      INTEGER N,NIfD,I,L,IR,J
      INTEGER RA(N)
      INTEGER RB(0:NIfD,N)
      INTEGER RRA,RRB(0:NIfD)
      INTEGER DetBitLT
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          RRB(0:NIfD)=RB(0:NIfD,L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          RRB(0:NIfD)=RB(0:NIfD,IR)
          RB(0:NIfD,IR)=RB(0:NIfD,1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RB(0:NIfD,1)=RRB(0:NIfD)
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
            RB(0:NIfD,I)=RB(0:NIfD,J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        RB(0:NIfD,I)=RRB(0:NIfD)

      GO TO 10

      END SUBROUTINE SortBitSign
   

! Based on SORTI, SortBitDets sorts determinants as bit strings, and takes the corresponding element from array RB with it (sign) and the real array RC (HElems)
! RA is the array of determinants of length N to sort
! The RA array elements go from 0:NIfD
! RB is the array of integers to go with the determinant
      SUBROUTINE SortBitDetswH(N,RA,NIfD,RB,RC)
      INTEGER N,NIfD,I,L,IR,J
      INTEGER RA(0:NIfD,N)
      INTEGER RB(N)
      REAL*8 RC(N),RRC
      INTEGER RRA(0:NIfD),RRB
      INTEGER DetBitLT
 
      IF(N.LE.1) RETURN
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA(0:NIfD)=RA(0:NIfD,L)
          RRB=RB(L)
          RRC=RC(L)
        ELSE
          RRA(0:NIfD)=RA(0:NIfD,IR)
          RA(0:NIfD,IR)=RA(0:NIfD,1)
          RRB=RB(IR)
          RB(IR)=RB(1)
          RRC=RC(IR)
          RC(IR)=RC(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(0:NIfD,1)=RRA(0:NIfD)
            RB(1)=RRB
            RC(1)=RRC
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF((DetBitLT(RA(0:NIfD,J),RA(0:NIfD,J+1),NIfD)).eq.1) J=J+1
          ENDIF
          IF((DetBitLT(RRA(0:NIfD),RA(0:NIfD,J),NIfD)).eq.1) THEN
            RA(0:NIfD,I)=RA(0:NIfD,J)
            RB(I)=RB(J)
            RC(I)=RC(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(0:NIfD,I)=RRA(0:NIfD)
        RB(I)=RRB
        RC(I)=RRC

      GO TO 10

      END SUBROUTINE SortBitDetswH




