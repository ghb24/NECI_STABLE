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
                IF(iLutI(0).ne.iLutJ(0)) THEN
                    DetBitEQ=.false.
                    RETURN
                ENDIF
            enddo
        ENDIF
        DetBitEQ=.true.

    END FUNCTION DetBitEQ

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
            IF((DetBitLT(RA(0:NIfD,J),RA(0:NIfD,J+1),NIfD)).eq.-1) J=J+1
          ENDIF
          IF((DetBitLT(RRA(0:NIfD),RA(0:NIfD,J),NIfD)).eq.-1) THEN
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




