PROGRAM ConvertPOPSFILE
      !Convert a SymDETS file to a POPSFILE
      IMPLICIT NONE
      INTEGER NIfD,nBasis,ExcitLev,nDet,i,Num,TotWalkers,TotParts,Parts,WriteParts,j,k,DetBitLT,junk,nOrbFull,nElecFroz,Index,NIfFull
      REAL*8 :: Norm,Weight,TotWeight,rat,r,FracPart
      INTEGER , ALLOCATABLE :: Det(:),Dets(:,:),Signs(:),Mapping(:),FullFrozOrbs(:)
      LOGICAL :: tHPHF,IsAllowedHPHFDet,IsOpenShell,tInit,tMapping

      junk=0

      WRITE(6,*) "How many basis functions?..."
      READ(*,*) nBasis

      WRITE(6,*) "How many particles do you want to write out?..."
      READ(*,*) TotWalkers

      NIfD=nBasis/32
      WRITE(6,*) "Number of integers needed to encode determinants: ",NIfD+1
      ALLOCATE(Det(0:NIfD))

      WRITE(6,*) "Do you want this POPSFILE to be suitable for use in an HPHF calculation? (T/F) "
      READ(*,*) tHPHF
      WRITE(6,*) "Is this POPSFILE to be used for an i-FCIQMC calculation (T), or a full spawning run (F)?"
      READ(*,*) tInit

      WRITE(6,*) "Do you want to use an orbital mapping to a full system?..."
      READ(*,*) tMapping

      IF(tMapping) THEN
          WRITE(6,*) "How many orbitals in the full space?..."
          READ(*,*) nOrbFull

          WRITE(6,*) "How many electrons have been frozen?..."
          READ(*,*) nElecFroz

          NIfFull=nOrbFull/32

          ALLOCATE(FullFrozOrbs(nElecFroz))

          IF(nElecFroz.ne.0) THEN
              OPEN(9,FILE='CoreFrozOrbs',STATUS='old',action='read')

              do i=1,nElecFroz
                  READ(9,*) FullFrozOrbs(i)
              enddo
              CLOSE(9)
          ENDIF
          ALLOCATE(Mapping(nBasis))
          OPEN(9,FILE='Mapping',STATUS='old',action='read')
          do i=1,nBasis
              READ(9,*) Index,Mapping(Index)   !Orbital index in frozen system , Orbital index in full system.
          enddo
          CLOSE(9)

      ELSE
          NIfFull=NIfD

      ENDIF


      OPEN(8,FILE='SymDETS',STATUS='old',action='read')

      Weight=0.D0
      Norm=0.D0
      TotWeight=0.D0

      do while(.true.)

101     READ(8,*,END=199) nDet,ExcitLev,Det(0:NIfD),Weight

          IF(tHPHF) THEN
              IF(IsAllowedHPHFDet(Det(0:NIfD),NIfD)) THEN
                  IF(IsOpenShell(Det(0:NIfD),NIfD)) THEN
                      Norm=Norm+Weight*Weight*2.D0
                      TotWeight=TotWeight+abs(Weight)*(sqrt(2.D0))
                  ELSE
                      Norm=Norm+Weight*Weight
                      TotWeight=TotWeight+abs(Weight)
                  ENDIF
              ENDIF
          ELSE
              TotWeight=TotWeight+abs(Weight)
              Norm=Norm+Weight**2
          ENDIF

      enddo

199   CONTINUE
      REWIND(8)
      WRITE(6,*) "Total determinants in SymDETS file: ",nDet
      WRITE(6,*) "Normalisation for wavefunction: ",Norm
      WRITE(6,*) "Total weight of amplitudes: ",TotWeight

      open(13,FILE='TEMP',STATUS='unknown')

      call Random_Seed()

      TotParts=0
      WriteParts=0
      do i=1,nDet

          READ(8,*) Num,ExcitLev,Det(0:NIfD),Weight
          IF(tHPHF) THEN
              IF(IsAllowedHPHFDet(Det(0:NIfD),NIfD)) THEN
                  IF(IsOpenShell(Det(0:NIfD),NIfD)) Weight=Weight*(sqrt(2.D0))
              ELSE
                  CYCLE
              ENDIF
          ENDIF 

          rat=REAL(TotWalkers,8)*abs(Weight)/TotWeight
          Parts=INT(rat)
          FracPart=rat-REAL(Parts,8)

          call random_number(r)
          IF(rat.gt.r) THEN
              Parts=Parts+1
          ENDIF

          IF(abs(Parts).gt.0) THEN
              IF(Weight.lt.0.D0) THEN
                  Parts=-Parts
              ENDIF

              WRITE(13,*) i,Det(0:NIfD),Parts
              TotParts=TotParts+abs(Parts)

              WriteParts=WriteParts+1
          ENDIF

      enddo

      WRITE(6,*) "Total number of particles stochastically created in POPSFILE = ", TotParts
      WRITE(6,*) "Total number of determinants occupied in POPSFILE = ", WriteParts
      CLOSE(8)

      ALLOCATE(Dets(0:NIfFull,WriteParts))
      ALLOCATE(Signs(WriteParts))

      REWIND(13)
      do i=1,WriteParts
          READ(13,*) j,Det(0:NIfD),Signs(i)
          IF(tMapping) THEN
              CALL ConvertDetToNewBasis(Det(0:NIfD),Dets(0:NIfFull,i),NIfD,NIfFull,nElecFroz,FullFrozOrbs,Mapping,nBasis)
          ELSE
              Dets(0:NIfFull,i)=Det(0:NIfD)
          ENDIF
      enddo
      close(13,status='delete')


      WRITE(6,*) "Sorting list of determinants..."
      CALL SortBitDets(WriteParts,Dets,Signs,NIfFull)

!      do i=2,WriteParts
!          WRITE(6,*) DetBitLT(Dets(0:NIfD,i),Dets(0:NIfD,i-1),NIfD)
!      enddo

      WRITE(6,*) "Writing to POPSFILE..."

      OPEN(17,FILE="POPSFILE",Status='replace')
      WRITE(17,*) WriteParts,"   TOTWALKERS (all nodes)"
      WRITE(17,*) 0.D0,"   DIAG SHIFT"
      WRITE(17,*) 0,"   SUMNOATHF (all nodes)"
      WRITE(17,*) 0.D0,"   SUMENUM ( \sum<D0|H|Psi> - all nodes)"
      WRITE(17,*) 0,"   PREVIOUS CYCLES"
      do i=1,WriteParts
          do k=0,NIfFull
              WRITE(17,"(I20)",advance='no') Dets(k,i)
          enddo
          IF(tInit) WRITE(17,"(I20)",advance='no') junk
          WRITE(17,*) Signs(i)
      enddo

      CLOSE(17)

END PROGRAM ConvertPOPSFILE

    LOGICAL FUNCTION IsOpenShell(iLut,NIfD)
        IMPLICIT NONE
        INTEGER :: iLut(0:NIfD),iLutAlpha(0:NIfD),iLutBeta(0:NIfD),MaskAlpha,MaskBeta,i,NIfD

        iLutAlpha(:)=0
        iLutBeta(:)=0
        MaskBeta=1431655765    !This is 1010101... in binary
        MaskAlpha=-1431655766  !This is 0101010... in binary
        IsOpenShell=.false.

        do i=0,NIfD

            iLutAlpha(i)=IAND(iLut(i),MaskAlpha)    !Seperate the alpha and beta bit strings
            iLutBeta(i)=IAND(iLut(i),MaskBeta)
            iLutAlpha(i)=ISHFT(iLutAlpha(i),-1)     !Shift all alpha bits to the left by one.
            iLutAlpha(i)=IEOR(iLutAlpha(i),iLutBeta(i)) !Do an XOR on the original beta bits and shifted alpha bits - they should cancel exactly.

            IF(iLutAlpha(i).ne.0) THEN
                IsOpenShell=.true.  !Det is open shell - return
                RETURN
            ENDIF
        enddo

    END FUNCTION IsOpenShell


    LOGICAL FUNCTION IsAllowedHPHFDet(iLutnI,NIfTot)
        INTEGER :: iLutSym(0:NIfTot),iLutnI(0:NIfTot),iLutTemp(0:NIfTot),i

        CALL FindExcitBitDetSym(iLutnI,iLutSym,NIfTot)

        ! iLutnI is 'less' than iLutSym, so iLutSym is the determinant with 
        ! the first open-shell = alpha. Swap them around.
        ! Only count up to NIfD to avoid Yamanouchi symbol etc.
        i=DetBitLT(iLutnI,iLutSym,NIfTot)

        IF(i.eq.1) THEN
            IsAllowedHPHFDet=.false.

!            iLutTemp(:)=iLutnI(:)
!            iLutnI(:)=iLutSym(:)
!            iLutSym(:)=iLutTemp(:)
!            tSwapped=.true.
!        ELSEIF(i.eq.0) THEN
!            CALL Stop_All("ReturnAlphaOpenDet","Shouldn't have closed shell determinants in here")
        ELSE
            IsAllowedHPHFDet=.true.
!            tSwapped=.false.
        ENDIF

    END FUNCTION IsAllowedHPHFDet

!In closed-shell systems with equal number of alpha and beta strings, the amplitude of a determinant in the final CI wavefunction is the same
!when the alpha and beta electrons are swapped (for S=0, see Helgakker for more details). It will sometimes be necessary to find this other
!determinant when spawning. This routine will find the bit-representation of an excitation by constructing the symmetric iLut from the its
!symmetric partner, also in bit form.
    SUBROUTINE FindExcitBitDetSym(iLut,iLutSym,NIfTot)
        IMPLICIT NONE
        INTEGER :: iLut(0:NIfTot),iLutSym(0:NIfTot)
        INTEGER :: iLutAlpha(0:NIfTot),iLutBeta(0:NIfTot),MaskAlpha,MaskBeta,i,NIfTot

!        WRITE(6,*) "******"
        iLutSym(:)=0
        iLutAlpha(:)=0
        iLutBeta(:)=0
        MaskBeta=1431655765    !This is 1010101... in binary
        MaskAlpha=-1431655766  !This is 0101010... in binary

!        WRITE(6,*) "MaskAlpha: "
!        do i=0,31
!            IF(BTEST(MaskAlpha,i)) THEN
!                WRITE(6,"(I3)",advance='no') 1
!            ELSE
!                WRITE(6,"(I3)",advance='no') 0
!            ENDIF
!        enddo
!        WRITE(6,*) ""
!        WRITE(6,*) "MaskBeta: "
!        do i=0,31
!            IF(BTEST(MaskBeta,i)) THEN
!                WRITE(6,"(I3)",advance='no') 1
!            ELSE
!                WRITE(6,"(I3)",advance='no') 0
!            ENDIF
!        enddo
!        WRITE(6,*) ""

        do i=0,NIfTot

            iLutAlpha(i)=IAND(iLut(i),MaskAlpha)    !Seperate the alpha and beta bit strings
            iLutBeta(i)=IAND(iLut(i),MaskBeta)

            iLutAlpha(i)=ISHFT(iLutAlpha(i),-1)  !Shift all alpha bits to the left by one.
            iLutBeta(i)=ISHFT(iLutBeta(i),1)   !Shift all beta bits to the right by one.

            iLutSym(i)=IOR(iLutAlpha(i),iLutBeta(i))    !Combine the bit strings to give the final bit representation.

!            WRITE(6,*) "ILut: "
!            do j=0,31
!                IF(BTEST(iLut(i),j)) THEN
!                    WRITE(6,"(I3)",advance='no') 1
!                ELSE
!                    WRITE(6,"(I3)",advance='no') 0
!                ENDIF
!            enddo
!            WRITE(6,*) ""
!            WRITE(6,*) "iLutSym: "
!            do j=0,31
!                IF(BTEST(iLutSym(i),j)) THEN
!                    WRITE(6,"(I3)",advance='no') 1
!                ELSE
!                    WRITE(6,"(I3)",advance='no') 0
!                ENDIF
!            enddo
!            WRITE(6,*) ""

        enddo

    END SUBROUTINE FindExcitBitDetSym




    ! This will return 1 if iLutI is "less" than iLutJ, 0 if the determinants
    ! are identical, or -1 if iLutI is "more" than iLutJ
    integer function DetBitLT(iLutI,iLutJ,NIfTot)
        integer, intent(in) :: NIfTot
        integer, intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer :: i, lnLast

        !First, compare first integers
        IF(iLutI(0).lt.iLutJ(0)) THEN
            DetBitLT=1
        ELSEIF(iLutI(0).eq.iLutJ(0)) THEN
            lnLast = NIftot

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


! Based on SORTI, SortBitDets sorts determinants as bit strings, and takes the corresponding element from array RB with it (sign)
! RA is the array of determinants of length N to sort
! The RA array elements go from 0:NIfD
! RB is the array of integers to go with the determinant
      SUBROUTINE SortBitDets(N,RA,RB,NIfTot)
      INTEGER N,I,L,IR,J
      INTEGER RA(0:NIfTot,N)
      INTEGER RB(N)
      INTEGER RRA(0:NIfTot),RRB
      INTEGER :: DetBitLT
 
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
            IF((DetBitLT(RA(:,J),RA(:,J+1),NIfTot)).eq.1) J=J+1
          ENDIF
          IF((DetBitLT(RRA(:),RA(:,J),NIfTot)).eq.1) THEN
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


      SUBROUTINE ConvertDetToNewBasis(iLut,iLutNew,NIfD,NIfFull,nElecFroz,FullFrozOrbs,Mapping,nBasis)
        IMPLICIT NONE
        INTEGER :: iLut(0:NIfD),iLutNew(0:NIfFull),NIfD,NIfFull,nElecFroz,FullFrozOrbs(nElecFroz),Mapping(nBasis),nBasis
        INTEGER :: i,j,elec,pos,OccOrb

            iLutNew(0:NIfFull)=0

            do i=0,NIfD
                do j=0,31
                    if(btest(iLut(i),j)) then
                        !An electron is at this orbital
                        elec=elec+1
                        OccOrb=(i*32)+(j+1)

                        pos = (Mapping(OccOrb)-1)/32
                        iLutNew(pos)=ibset(iLutNew(pos),mod(Mapping(OccOrb)-1,32))

                    endif
                enddo
            enddo

            do i=1,nElecFroz
                pos = (FullFrozOrbs(i) - 1) / 32
                iLutNew(pos)=ibset(iLutNew(pos),mod(FullFrozOrbs(i) - 1,32))
            enddo

      END SUBROUTINE ConvertDetToNewBasis

