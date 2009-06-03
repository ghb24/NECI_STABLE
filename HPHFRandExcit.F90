MODULE HPHFRandExcitMod
!Half-projected HF wavefunctions are a linear combination of two HF determinants, where all alphas -> betas and betas -> alpha to create the pair.
!In closed-shell systems, these two determinants have the same FCI amplitude, and so it is easier to treat them as a pair.
!The probability of creating this HPHF where both are from pairs of spin-coupled determinants (i & j -> a & b):
![ P(i->a) + P(i->b) + P(j->a) + P(j->b) ]/2
!We therefore need to find the excitation matrix between the determinant which wasn't excited and the determinant which was created.

    use SystemData, only: nEl,tMerTwist,NIfD
    use SymData, only: nSymLabels
    use mt95 , only : genrand_real2
    use GenRandSymExcitNUMod , only : GenRandSymExcitScratchNU,ConstructClassCounts,CalcNonUniPGen 
    IMPLICIT NONE

    contains

!nI will always be the determinant with the first open-shell having an alpha spin-orbital occupied.
    SUBROUTINE GenRandHPHFExcit(nI,iLutnI,nJ,iLutnJ,pDoub,exFlag,pGen)
        INTEGER :: nI(NEl),iLutnI(0:NIfD),iLutnJ(0:NIfD),nJ(NEl),exFlag,ExcitMat(2,2),IC
        INTEGER :: iLutnJ2(0:NIfD),nI2(NEl),nJ2(NEl),Ex2(2,2),ExcitLevel,iLutnI2(0:NIfD)
        REAL*8 :: pDoub,pGen,r,pGen2
        INTEGER :: ClassCount2(2,0:nSymLabels-1),ClassCount3(2,0:nSymLabels-1)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1),ClassCountUnocc3(2,0:nSymLabels-1)
        LOGICAL :: tGenClassCountnI,tGenClassCountnI2,TestClosedShellDet,tParity,tSign

        tGenClassCountnI=.false.
        tGenClassCountnI2=.false.

!Test is nI is a closed-shell determinant
        IF(TestClosedShellDet(iLutnI,NIfD)) THEN
!If determinant is closed shell, then all probabilities are the same, so P=2*Prob since both spins are equally likely to be generated (as long as generates open shell HPHF).
!Just need to return the right spin.

            CALL GenRandSymExcitScratchNU(nI,iLutnI,nJ,pDoub,IC,ExcitMat,tParity,exFlag,pGen,ClassCount2,ClassCountUnocc2,tGenClassCountnI)
            
!Create bit representation of excitation - iLutnJ
            CALL FindExcitBitDet(iLutnI,iLutnJ,IC,ExcitMat,NIfD)

            IF(IC.eq.2) THEN
                IF(.not.TestClosedShellDet(iLutnJ,NIfD)) THEN
                    pGen=pGen*2.D0
                    CALL ReturnAlphaOpenDet(nJ,iLutnJ,iLutnJ2,.true.)
                ELSE
!Excitation is closed shell: Closed shell -> Closed Shell
                    RETURN
                ENDIF
            ELSE
!Excitation is definitely open-shell
                pGen=pGen*2.D0
                CALL ReturnAlphaOpenDet(nJ,iLutnJ,iLutnJ2,.true.)
            ENDIF

            RETURN
        ENDIF

!If det is open-shell we choose one of the determinants with 50% chance to create an excitation from.
        IF(tMerTwist) THEN
            CALL genrand_real2(r)
        ELSE
            CALL RANLUX(r,1)
        ENDIF
!This will find the full ordered form for nI2 and its bit representation. (Is this always needed?)
        CALL FindDetSpinSym(nI,nI2,NEl)
        CALL FindExcitBitDetSym(iLutnI,iLutnI2)

        IF(r.lt.0.D5) THEN
!Excite from nJ from nI
            CALL GenRandSymExcitScratchNU(nI,iLutnI,nJ,pDoub,IC,ExcitMat,tParity,exFlag,pGen,ClassCount2,ClassCountUnocc2,tGenClassCountnI)

!Find Bit-representation of excitation.
            CALL FindExcitBitDet(iLutnI,iLutnJ,IC,ExcitMat,NIfD)
            IF(TestClosedShellDet(iLutnJ,NIfD)) THEN
!Excitation created is a closed shell determinant. Both determinants are connected to it, and crucially with the same probability. This means that the final pGen is unchanged.
                RETURN
            ENDIF

!We may have been able to excite from nI2 to this determinant. see if it in connected.
            CALL FindBitExcitLevel(iLutnI2,iLutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
!                CALL DecodeBitDet(nI2,iLutnI2,NIfD)
                Ex2(1,1)=ExcitLevel
                CALL GetExcitation(nI2,nJ,NEl,Ex2,tSign)
                tGenClassCountnI2=.true.
                CALL ConstructClassCounts(nI2,ClassCount3,ClassCountUnocc3)
                CALL CalcNonUniPGen(Ex2,ExcitLevel,ClassCount3,ClassCountUnocc3,pDoub,pGen2)
                pGen=pGen+pGen2
            ENDIF

        ELSE
!Excite from the spin-pair of nI (called nI2)

!            CALL DecodeBitDet(nI2,iLutnI2,NEl,NIfD)
!            CALL FindDetSpinSym(nI,nI2,NEl)
            CALL GenRandSymExcitScratchNU(nI2,iLutnI2,nJ,pDoub,IC,ExcitMat,tParity,exFlag,pGen,ClassCount3,ClassCountUnocc3,tGenClassCountnI2)

!Find Bit-representation of excitation.
            CALL FindExcitBitDet(iLutnI2,iLutnJ,IC,ExcitMat,NIfD)
            IF(TestClosedShellDet(iLutnJ,NIfD)) THEN
!Excitation created is a closed shell determinant. Both determinants are connected to it, and crucially with the same probability. This means that the final pGen is unchanged.
                RETURN
            ENDIF

!We know we have gone from open-shell HPHF to open-shell HPHF. We need all four pGens.
!We have nI2 -> nJ. Find nI -> nJ. First, we need to know whether it is connected or not.
            CALL FindBitExcitLevel(iLutnI,iLutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                Ex2(1,1)=ExcitLevel
                CALL GetExcitation(nI,nJ,NEl,Ex2,tSign)
!We need to calculate the new classcount arrays for the original determinant passed in.
                tGenClassCountnI=.true.
                CALL ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2)
                CALL CalcNonUniPGen(Ex2,ExcitLevel,ClassCount2,ClassCountUnocc2,pDoub,pGen2)
                pGen=pGen+pGen2
            ENDIF

        ENDIF

!We also need to look at how we *could* have excited to the spin-coupled determinant of nJ.
        CALL FindExcitBitDetSym(iLutnJ,iLutnJ2)
        CALL FindDetSpinSym(nJ,nJ2,NEl)

!Firstly, nI2 -> nJ2
        CALL FindBitExcitLevel(iLutnI2,iLutnJ2,NIfD,ExcitLevel,2)
        IF(ExcitLevel.le.2) THEN
            Ex2(1,1)=ExcitLevel
            CALL GetExcitation(nI2,nJ2,NEl,Ex2,tSign)
            IF(.not.tGenClassCountnI2) THEN
!                tGenClassCountnI2=.true.
                CALL ConstructClassCounts(nI2,ClassCount3,ClassCountUnocc3)
            ENDIF
            CALL CalcNonUniPGen(Ex2,ExcitLevel,ClassCount3,ClassCountUnocc3,pDoub,pGen2)
            pGen=pGen+pGen2
        ENDIF

!Finally, nI -> nJ2
        CALL FindBitExcitLevel(iLutnI,iLutnJ2,NIfD,ExcitLevel,2)
        IF(ExcitLevel.le.2) THEN
            Ex2(1,1)=ExcitLevel
            CALL GetExcitation(nI,nJ2,NEl,Ex2,tSign)
            IF(.not.tGenClassCountnI) THEN
!                tGenClassCountnI=.true.
                CALL ConstructClassCounts(nI,ClassCount2,ClassCountUnocc2)
            ENDIF
            CALL CalcNonUniPGen(Ex2,ExcitLevel,ClassCount2,ClassCountUnocc2,pDoub,pGen2)
            pGen=pGen+pGen2
        ENDIF

        pGen=pGen/2.D0  !Normalize pGens.

!Excitation is open-shell. We need to find the correct spin-pair to send back. 
        CALL ReturnAlphaOpenDet(nJ,iLutnJ,iLutnJ2,.false.)

    END SUBROUTINE GenRandHPHFExcit

!This routine will take a determinant, and create the determinant whose first open-shell spatial orbital contains an alpha electron.
!If the first open-shell electron is a beta orbital, then the balue of the bit-string will be smaller. We are interested in returning
!the larger of the open-shell bit strings since this will correspond to the first open-shell electron being an alpha.
!iLutnI (nI) is returned as this determinant, with iLutSym (nJ) being the other.
!If tCalciLutSym is false, iLutSym will be calculated from iLutnI. Otherwise, it won't.
    SUBROUTINE ReturnAlphaOpenDet(nI,iLutnI,iLutSym,tCalciLutSym)
        INTEGER :: iLutSym(0:NIfD),nI(NEl),iLutnI(0:NIfD),nJ(NEl),iLutTemp(0:NIfD),i
        LOGICAL :: tCalciLutSym,DetBitLT

        IF(tCalciLutSym) THEN
            CALL FindExcitBitDetSym(iLutnI,iLutSym)
        ENDIF

        i=DetBitLT(iLutnI,iLutSym,NIfD)
        IF(i.eq.1) THEN
!iLutnI is 'less' than iLutSym, so iLutSym is the determinant with the first open-shell = alpha. Swap them around.
            iLutTemp(:)=iLutnI(:)
            iLutnI(:)=iLutSym(:)
            iLutSym(:)=iLutTemp(:)
            CALL FindDetSpinSym(nI,nJ,NEl)
            nI(:)=nJ(:)
        ELSEIF(i.eq.0) THEN
            CALL Stop_All("ReturnAlphaOpenDet","Shouldn't have closed shell determinants in here")
        ENDIF

    END SUBROUTINE ReturnAlphaOpenDet
        

!There will be a quicker way to do this without needing the sort.
!This create the spin-coupled determinant of nI in nJ in natural ordered form.
    SUBROUTINE FindDetSpinSym(nI,nJ,NEl)
        INTEGER :: nI(NEl),nJ(NEl),NEl,i

        do i=1,NEl
            IF(mod(nI(i),2).eq.0) THEN
!electron is an alpha - change it to a beta (remove one)
                nJ(i)=nI(i)-1
            ELSE
                nJ(i)=nI(i)+1
            ENDIF
        enddo
        CALL NECI_SORTI(NEl,nJ)

    END SUBROUTINE FindDetSpinSym


!In closed-shell systems with equal number of alpha and beta strings, the amplitude of a determinant in the final CI wavefunction is the same
!when the alpha and beta electrons are swapped (for S=0, see Helgakker for more details). It will sometimes be necessary to find this other
!determinant when spawning. This routine will find the bit-representation of an excitation by constructing the symmetric iLut from the its
!symmetric partner, also in bit form.
    SUBROUTINE FindExcitBitDetSym(iLut,iLutSym)
        IMPLICIT NONE
        INTEGER :: iLut(0:NIfD),iLutSym(0:NIfD)
        INTEGER :: iLutAlpha(0:NIfD),iLutBeta(0:NIfD),MaskAlpha,MaskBeta,i

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

        do i=0,NIfD

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


END MODULE HPHFRandExcitMod


