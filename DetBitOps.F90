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
            RETURN
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
        ELSE
            DetBitLT=-1
            RETURN
        ENDIF
        DetBitLT=0

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


!This is a routine to encode a determinant as natural ordered integers (nI) as a bit string (iLut(0:NIfD))
!where NIfD=INT(nBasis/32)
    SUBROUTINE EncodeBitDet(nI,iLut,NEl,NIfD)
        IMPLICIT NONE
        INTEGER :: nI(NEl),iLut(0:NIfD),NoIntforDet,i,NEl,NIfD

        iLut(:)=0
        do i=1,NEl
            iLut((nI(i)-1)/32)=IBSET(iLut((nI(i)-1)/32),mod(nI(i)-1,32))
        enddo

    END SUBROUTINE EncodeBitDet

!This is a routine to take a determinant in bit form and construct the natural ordered NEl integer for of the det.
    SUBROUTINE DecodeBitDet(nI,iLut,NEl,NIfD)
        IMPLICIT NONE
        INTEGER :: nI(NEl),iLut(0:NIfD),NIfD,i,j,NEl,Elec

        Elec=0
        nI(:)=0
        do i=0,NIfD
            do j=0,31
                IF(BTEST(iLut(i),j)) THEN
!An electron is at this orbital
                    nI(Elec+1)=(i*32)+(j+1)
                    Elec=Elec+1
                    IF(Elec.eq.NEl) RETURN
                ENDIF
            enddo
        enddo

    END SUBROUTINE DecodeBitDet

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

    SUBROUTINE CalcOpenOrbs(iLut,NIfD,NEl,OpenOrbs)
        INTEGER :: iLut(0:NIfD),iLutAlpha(0:NIfD),iLutBeta(0:NIfD),MaskAlpha,MaskBeta,i,NIfD,NEl,OpenOrbs
        
        iLutAlpha(:)=0
        iLutBeta(:)=0
        MaskBeta=1431655765    !This is 1010101... in binary
        MaskAlpha=-1431655766  !This is 0101010... in binary

        do i=0,NIfD

            iLutAlpha(i)=IAND(iLut(i),MaskAlpha)    !Seperate the alpha and beta bit strings
            iLutBeta(i)=IAND(iLut(i),MaskBeta)
            iLutAlpha(i)=ISHFT(iLutAlpha(i),-1)     !Shift all alpha bits to the left by one.
            iLutAlpha(i)=IEOR(iLutAlpha(i),iLutBeta(i)) !Do an XOR on the original beta bits and shifted alpha bits - only open shell occupied orbitals will remain.
            
        enddo

        CALL CountBits(iLutAlpha,NIfD,OpenOrbs,NEl)
        OpenOrbs=OpenOrbs/2

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




