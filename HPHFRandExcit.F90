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
!    SAVE
!    INTEGER :: Count=0

    contains

!nI will always be the determinant with the first open-shell having an alpha spin-orbital occupied.
    SUBROUTINE GenRandHPHFExcit(nI,iLutnI,nJ,iLutnJ,pDoub,exFlag,pGen)
        INTEGER :: nI(NEl),iLutnI(0:NIfD),iLutnJ(0:NIfD),nJ(NEl),exFlag,ExcitMat(2,2),IC
        INTEGER :: iLutnJ2(0:NIfD),nI2(NEl),nJ2(NEl),Ex2(2,2),ExcitLevel,iLutnI2(0:NIfD)
        REAL*8 :: pDoub,pGen,r,pGen2
        INTEGER :: ClassCount2(2,0:nSymLabels-1),ClassCount3(2,0:nSymLabels-1)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1),ClassCountUnocc3(2,0:nSymLabels-1)
        LOGICAL :: tGenClassCountnI,tGenClassCountnI2,TestClosedShellDet,tParity,tSign

!        Count=Count+1
!        WRITE(6,*) "COUNT: ",Count
!        CALL FLUSH(6)
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
!Excite to nJ from nI
            CALL GenRandSymExcitScratchNU(nI,iLutnI,nJ,pDoub,IC,ExcitMat,tParity,exFlag,pGen,ClassCount2,ClassCountUnocc2,tGenClassCountnI)

!Find Bit-representation of excitation.
            CALL FindExcitBitDet(iLutnI,iLutnJ,IC,ExcitMat,NIfD)
            IF(TestClosedShellDet(iLutnJ,NIfD)) THEN
!Excitation created is a closed shell determinant. Both determinants are connected to it, and crucially with the same probability. This means that the final pGen is unchanged.
                RETURN
            ENDIF

!We may have been able to excite from nI2 to this determinant. see if it in connected.
            CALL FindBitExcitLevel(iLutnI2,iLutnJ,NIfD,ExcitLevel,2)
            IF((ExcitLevel.le.2).and.(ExcitLevel.ne.0)) THEN
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
            IF((ExcitLevel.le.2).and.(ExcitLevel.ne.0)) THEN
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
        IF((ExcitLevel.le.2).and.(ExcitLevel.ne.0)) THEN
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
        IF((ExcitLevel.le.2).and.(ExcitLevel.ne.0)) THEN
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

!This routine will be faster than the one above - it only generates excitations from one of the determinants in the HPHF function to be excited from.
!This relies on the fact that both determinants in the HPHF function to be excited from will always be connected to all excited HPHF functions.
!nI will always need to be a unique choice of determinant within each HPHF function, and then we never actually need to refer to its spin-coupled partner.
    SUBROUTINE GenRandHPHFExcit2Scratch(nI,iLutnI,nJ,iLutnJ,pDoub,exFlag,pGen,ClassCount2,ClassCountUnocc2,tGenClassCountnI)
        INTEGER :: nI(NEl),iLutnI(0:NIfD),iLutnJ(0:NIfD),nJ(NEl),exFlag,IC,ExcitMat(2,2)
        INTEGER :: iLutnJ2(0:NIfD),nJ2(NEl),Ex2(2,2),ExcitLevel
        REAL*8 :: pDoub,pGen,pGen2
        INTEGER :: ClassCount2(2,0:nSymLabels-1)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)
        LOGICAL :: tGenClassCountnI,TestClosedShellDet,tSign

!        Count=Count+1
!        WRITE(6,*) "COUNT: ",Count
!        CALL FLUSH(6)

!Create excitation of uniquely chosen determinant in this HPHF function.
        CALL GenRandSymExcitScratchNU(nI,iLutnI,nJ,pDoub,IC,ExcitMat,tSign,exFlag,pGen,ClassCount2,ClassCountUnocc2,tGenClassCountnI)
!Create bit representation of excitation - iLutnJ
        CALL FindExcitBitDet(iLutnI,iLutnJ,IC,ExcitMat,NIfD)
            
!Test!
!        CALL CalcNonUniPGen(ExcitMat,IC,ClassCount2,ClassCountUnocc2,pDoub,pGen2)
!        IF(abs(pGen-pGen2).gt.1.D-7) THEN
!            WRITE(6,*) "*******, PGens Incorrect"
!            CALL Stop_All("ouvbou","OUBOU")
!        ENDIF

!        IF(Count.eq.4) THEN
!            WRITE(6,*) "***",nI(:),iLutnJ(:)
!        ENDIF

        IF(TestClosedShellDet(iLutnJ,NIfD)) THEN
!There is only one way which we could have generated the excitation nJ since it has no spin-partner. Also, we will always return the 'correct' version.
            RETURN
        ELSE
!Open shell excitation - could we have generated the spin-coupled determinant instead?

            CALL FindExcitBitDetSym(iLutnJ,iLutnJ2)
!            IF(Count.eq.4) THEN
!                WRITE(6,*) "***",nI(:),iLutnJ(:),iLutnJ2(:)
!            ENDIF
            CALL FindBitExcitLevel(iLutnJ2,iLutnI,NIfD,ExcitLevel,2)

            IF((ExcitLevel.eq.2).or.(ExcitLevel.eq.1)) THEN

                Ex2(1,1)=ExcitLevel
                CALL DecodeBitDet(nJ2,iLutnJ2,NEl,NIfD)
!                IF(ExcitLevel.eq.1) THEN
!                    WRITE(6,*) "SINGLE EXCITATION",nJ2(:),nI(:)
!                ENDIF

                CALL GetExcitation(nI,nJ2,NEl,Ex2,tSign) !This could be done more efficiently... !***!
                CALL CalcNonUniPGen(Ex2,ExcitLevel,ClassCount2,ClassCountUnocc2,pDoub,pGen2)    !Check this is done efficiently... !***!
!                IF(abs(pGen-pGen2).gt.1.D-7) THEN
!!We cannot guarentee that the pGens are going to be the same - in fact, generally, they wont be.
!                    WRITE(6,*) nI(:),nJ2(:)
!                    WRITE(6,*) "*************  PGens NOT equal",pGen,pGen2
!                ENDIF
                pGen=pGen+pGen2
                CALL ReturnAlphaOpenDet(nJ,iLutnJ,iLutnJ2,.false.)  !Here, we actually know nJ, so wouldn't need to regenerate it...
                RETURN
            ENDIF

            CALL ReturnAlphaOpenDet(nJ,iLutnJ,iLutnJ2,.false.)

        ENDIF

    END SUBROUTINE GenRandHPHFExcit2Scratch

!This routine will take a determinant, and create the determinant whose final open-shell spatial orbital contains an alpha electron.
!If the final open-shell electron is a beta orbital, then the balue of the bit-string will be smaller. We are interested in returning
!the larger of the open-shell bit strings since this will correspond to the final open-shell electron being an alpha.
!This rationalization may well break down when it comes to the negative bit (32), however, this may not matter, since all we really
!need is a unique description of a HPHF...?
!iLutnI (nI) is returned as this determinant, with iLutSym (nJ) being the other.
!If tCalciLutSym is false, iLutSym will be calculated from iLutnI. Otherwise, it won't.
    SUBROUTINE ReturnAlphaOpenDet(nI,iLutnI,iLutSym,tCalciLutSym)
        INTEGER :: iLutSym(0:NIfD),nI(NEl),iLutnI(0:NIfD),nJ(NEl),iLutTemp(0:NIfD),i
        INTEGER :: DetBitLT
        LOGICAL :: tCalciLutSym

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
        INTEGER :: nTemp(NEl),nI(NEl),nJ(NEl),NEl,i

        do i=1,NEl
            IF(mod(nI(i),2).eq.0) THEN
!electron is an alpha - change it to a beta (remove one)
!However, we only want to do this if the electron before it is not the beta in the same spatial orbital
                IF((i.eq.1).or.((nI(i)-1).ne.nI(i-1))) THEN
                    nJ(i)=nI(i)-1
                ELSE
                    nJ(i)=nI(i)
                ENDIF
            ELSE
                IF((i.eq.NEl).or.((nI(i)+1).ne.nI(i+1))) THEN
                    nJ(i)=nI(i)+1
                ELSE
                    nJ(i)=nI(i)
                ENDIF
            ENDIF
        enddo

!        nTemp(:)=nJ(:)
!        CALL NECI_SORTI(NEl,nTemp)
!        do i=1,NEl
!            IF(nTemp(i).ne.nJ(i)) THEN
!                STOP 'Massive Error'
!            ENDIF
!        enddo

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

    
!This routine will take a HPHF nI, and find Iterations number of excitations of it. It will then histogram these, summing in 1/pGen for every occurance of
!the excitation. This means that all excitations should be 0 or 1 after enough iterations. It will then count the excitations and compare the number to the
!number of excitations generated using the full enumeration excitation generation.
    SUBROUTINE TestGenRandHPHFExcit(nI,Iterations,pDoub)
        Use SystemData , only : NEl,nBasis,G1,nBasisMax,NIfD
        IMPLICIT NONE
        INTEGER :: ClassCount2(2,0:nSymLabels-1),nIX(NEl)
        INTEGER :: ClassCountUnocc2(2,0:nSymLabels-1)
        INTEGER :: i,Iterations,nI(NEl),nJ(NEl),DetConn,nI2(NEl),DetConn2,iUniqueHPHF,iUniqueBeta,PartInd,ierr,iExcit
        REAL*8 :: pDoub,pGen
        LOGICAL :: Unique,TestClosedShellDet,DetBitEQ,Die,tGenClassCountnI
        INTEGER :: iLutnI(0:NIfD),iLutnJ(0:NIfD),iLutnI2(0:NIfD),iLutSym(0:NIfD)
        INTEGER , ALLOCATABLE :: ConnsAlpha(:,:),ConnsBeta(:,:),ExcitGen(:),UniqueHPHFList(:,:)
        REAL*8 , ALLOCATABLE :: Weights(:)
        INTEGER :: iMaxExcit,nStore(6),nExcitMemLen,j,k,l

        CALL EncodeBitDet(nI,iLutnI,NEl,NIfD)
        CALL FindDetSpinSym(nI,nI2,NEl)
        CALL EncodeBitDet(nI2,iLutnI2,NEl,NIfD)
        IF(TestClosedShellDet(iLutnI,NIfD)) THEN
            IF(.not.DetBitEQ(iLutnI,iLutnI2,NIfD)) THEN
                CALL Stop_All("TestGenRandHPHFExcit","Closed shell determinant entered, but alpha and betas different...")
            ENDIF
        ENDIF
        WRITE(6,*) "nI: ",nI(:)
        WRITE(6,*) ""
        WRITE(6,*) "nISym: ",nI2(:)
        WRITE(6,*) ""
        WRITE(6,*) "iLutnI: ",iLutnI(:)
        WRITE(6,*) "iLutnISymL ",iLutnI2(:)
        WRITE(6,*) "***"
        WRITE(6,*) Iterations,pDoub
!        WRITE(6,*) "nSymLabels: ",nSymLabels
        CALL FLUSH(6)

!First, we need to enumerate all possible HPHF wavefunctions from each spin-pair of determinants.
!These need to be stored in an array
!Find the number of symmetry allowed excitations there should be by looking at the full excitation generator.
!Setup excit generators for this determinant
        iMaxExcit=0
        nStore(1:6)=0
        CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
        ALLOCATE(EXCITGEN(nExcitMemLen),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
        EXCITGEN(:)=0
        CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,EXCITGEN,nJ,iMaxExcit,0,nStore,3)
        CALL GetSymExcitCount(EXCITGEN,DetConn)
        WRITE(6,*) "Alpha determinant has ",DetConn," total excitations:"
        ALLOCATE(ConnsAlpha(0:NIfD,DetConn))
        i=1
        lp2: do while(.true.)
            CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.false.,EXCITGEN,nJ,iExcit,0,nStore,3)
            IF(nJ(1).eq.0) exit lp2
            CALL EncodeBitDet(nJ,iLutnJ,NEl,NIfD)
            IF(.not.TestClosedShellDet(iLutnJ,NIfD)) THEN
                CALL ReturnAlphaOpenDet(nJ,iLutnJ,iLutSym,.true.)
            ENDIF
!            WRITE(6,"(4I4,A,I4,A,I13)") nJ(:), " *** ", iExcit, " *** ", iLutnJ(:)
            ConnsAlpha(0:NIfD,i)=iLutnJ(:)
            i=i+1
        enddo lp2

!Now we also need to store the excitations from the other spin-coupled determinant.
        iMaxExcit=0
        nStore(1:6)=0
        DEALLOCATE(EXCITGEN)
        
        CALL GenSymExcitIt2(nI2,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
        ALLOCATE(EXCITGEN(nExcitMemLen),stat=ierr)
        IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
        EXCITGEN(:)=0
        CALL GenSymExcitIt2(nI2,NEl,G1,nBasis,nBasisMax,.TRUE.,EXCITGEN,nJ,iMaxExcit,0,nStore,3)
        CALL GetSymExcitCount(EXCITGEN,DetConn2)
        WRITE(6,*) "Beta determinant has ",DetConn2," total excitations"
        ALLOCATE(ConnsBeta(0:NIfD,DetConn2))
        i=1
        lp: do while(.true.)
            CALL GenSymExcitIt2(nI2,NEl,G1,nBasis,nBasisMax,.false.,EXCITGEN,nJ,iExcit,0,nStore,3)
            IF(nJ(1).eq.0) exit lp
            CALL EncodeBitDet(nJ,iLutnJ,NEl,NIfD)
            IF(.not.TestClosedShellDet(iLutnJ,NIfD)) THEN
                CALL ReturnAlphaOpenDet(nJ,iLutnJ,iLutSym,.true.)
            ENDIF
!            WRITE(6,"(4I4,A,I4,A,I13)") nJ(:), " *** ",iExcit," *** ",iLutnJ(:)
            ConnsBeta(0:NIfD,i)=iLutnJ(:)
            i=i+1
        enddo lp
        DEALLOCATE(EXCITGEN)

!Now we need to find how many HPHF functions there are.
        iUniqueHPHF=0
        do j=1,DetConn
!Run though all HPHF in the first array
            Unique=.true.
            do k=j-1,1,-1
!Run backwards through the array to see if this HPHF has come before
                IF(DetBitEQ(ConnsAlpha(0:NIfD,k),ConnsAlpha(0:NIfD,j),NIfD)) THEN
!This HPHF has already been counted before...
                    Unique=.false.
                    EXIT
                ENDIF
            enddo
            IF(Unique) THEN
!Unique HPHF found, count it
                iUniqueHPHF=iUniqueHPHF+1
            ENDIF
        enddo

        iUniqueBeta=0

!Now look through all the excitations for the spin-coupled determinant from the original HPHF...
        do j=1,DetConn2
!Run though all excitations in the first array, *and* up to where we are in the second array
            Unique=.true.
            do k=1,DetConn
                IF(DetBitEQ(ConnsAlpha(0:NIfD,k),ConnsBeta(0:NIfD,j),NIfD)) THEN
                    Unique=.false.
                    EXIT
                ENDIF
            enddo
            IF(Unique) THEN
!Need to search backwards through the entries we've already looked at in this array...
                do k=j-1,1,-1
                    IF(DetBitEQ(ConnsBeta(0:NIfD,k),ConnsBeta(0:NIfD,j),NIfD)) THEN
                        Unique=.false.
                        EXIT
                    ENDIF
                enddo
            ENDIF
            IF(Unique) THEN
                iUniqueHPHF=iUniqueHPHF+1
                iUniqueBeta=iUniqueBeta+1
            ENDIF
        enddo

        WRITE(6,*) "There are ",iUniqueHPHF," unique HPHF wavefunctions from the HPHF given."
        
        WRITE(6,*) "There are ",iUniqueBeta," unique HPHF wavefunctions from the spin-coupled determinant, which are not in a alpha version."
        IF(iUniqueBeta.ne.0) THEN
            WRITE(6,*) "HPHF from beta, but not from alpha!"
            CALL FLUSH(6)
            STOP
        ENDIF

        ALLOCATE(UniqueHPHFList(0:NIfD,iUniqueHPHF))
        UniqueHPHFList(:,:)=0
!Now fill the list of HPHF Excitations.
        iUniqueHPHF=0
        do j=1,DetConn
!Run though all HPHF in the first array
            Unique=.true.
            do k=j-1,1,-1
!Run backwards through the array to see if this HPHF has come before
                IF(DetBitEQ(ConnsAlpha(0:NIfD,k),ConnsAlpha(0:NIfD,j),NIfD)) THEN
!This HPHF has already been counted before...
                    Unique=.false.
                    EXIT
                ENDIF
            enddo
            IF(Unique) THEN
!Unique HPHF found, count it
                iUniqueHPHF=iUniqueHPHF+1
                UniqueHPHFList(:,iUniqueHPHF)=ConnsAlpha(0:NIfD,j)
            ENDIF
        enddo

!Now look through all the excitations for the spin-coupled determinant from the original HPHF...
        do j=1,DetConn2
!Run though all excitations in the first array, *and* up to where we are in the second array
            Unique=.true.
            do k=1,DetConn
                IF(DetBitEQ(ConnsAlpha(0:NIfD,k),ConnsBeta(0:NIfD,j),NIfD)) THEN
                    Unique=.false.
                    EXIT
                ENDIF
            enddo
            IF(Unique) THEN
!Need to search backwards through the entries we've already looked at in this array...
                do k=j-1,1,-1
                    IF(DetBitEQ(ConnsBeta(0:NIfD,k),ConnsBeta(0:NIfD,j),NIfD)) THEN
                        Unique=.false.
                        EXIT
                    ENDIF
                enddo
            ENDIF
            IF(Unique) THEN
                iUniqueHPHF=iUniqueHPHF+1
                UniqueHPHFList(:,iUniqueHPHF)=ConnsBeta(0:NIfD,j)
            ENDIF
        enddo

!Now sort the list, so that it can be easily binary searched.
        ALLOCATE(ExcitGen(iUniqueHPHF))
        CALL SortBitDets(iUniqueHPHF,UniqueHPHFList(:,1:iUniqueHPHF),NIfD,ExcitGen)
        DEALLOCATE(ExcitGen)

        WRITE(6,*) "Unique HPHF wavefunctions are: "
        do i=1,iUniqueHPHF
            WRITE(6,*) UniqueHPHFList(0:NIfD,i)
        enddo

        ALLOCATE(Weights(iUniqueHPHF))
        Weights(:)=0.D0
        tGenClassCountnI=.false.

        do i=1,Iterations

            IF(mod(i,10000).eq.0) WRITE(6,"(A,I10)") "Iteration: ",i

            CALL GenRandHPHFExcit(nI,iLutnI,nJ,iLutnJ,pDoub,3,pGen)
            CALL GenRandHPHFExcit2Scratch(nI,iLutnI,nJ,iLutnJ,pDoub,3,pGen,ClassCount2,ClassCountUnocc2,tGenClassCountnI)
!            CALL GenRandSymExcitNU(nI,iLut,nJ,pDoub,IC,ExcitMat,TParity,exFlag,pGen)

!Search through the list of HPHF wavefunctions to find slot.
            CALL BinSearchListHPHF(iLutnJ,UniqueHPHFList(0:NIfd,1:iUniqueHPHF),iUniqueHPHF,1,iUniqueHPHF,PartInd,Unique)

            IF(.not.Unique) THEN
                CALL Stop_All("TestGenRandHPHFExcit","Cannot find excitation in list of allowed excitations")
            ENDIF

            Weights(PartInd)=Weights(PartInd)+(1.D0/pGen)
             
!Check excitation
!            CALL IsSymAllowedExcit(nI,nJ,IC,ExcitMat,SymAllowed)

        enddo
        
        OPEN(8,FILE="PGenHist",STATUS="UNKNOWN")

!normalise excitation probabilities
        Die=.false.
        do i=1,iUniqueHPHF
            Weights(i)=Weights(i)/real(Iterations,8)
            IF(abs(Weights(i)-1.D0).gt.0.1) THEN
                WRITE(6,*) "Error here!"
                Die=.true.
            ENDIF
            WRITE(6,*) i,UniqueHPHFList(0:NIfD,i),Weights(i)
            CALL DecodeBitDet(nIX,UniqueHPHFList(0:NIfD,i),NEl,NIfD)
            WRITE(6,*) nIX(:)
            WRITE(8,*) i,UniqueHPHFList(0:NIfD,i),Weights(i)
        enddo

        CLOSE(8)
        IF(Die) THEN
            CALL Stop_All("IUB","TestFail")
        ENDIF

    END SUBROUTINE TestGenRandHPHFExcit

    SUBROUTINE BinSearchListHPHF(iLut,List,Length,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: iLut(0:NIfD),MinInd,MaxInd,PartInd
        INTEGER :: List(0:NIfD,Length),Length
        INTEGER :: i,j,N,Comp,DetBitLT
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL FLUSH(6)
        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(List(:,MaxInd),iLut(:),NIfD)
            IF(Comp.eq.0) THEN
                tSuccess=.true.
                PartInd=MaxInd
                RETURN
            ELSE
                tSuccess=.false.
                PartInd=MinInd
            ENDIF
        ENDIF
        do while(j-i.gt.0)  !End when the upper and lower bound are the same.
            N=(i+j)/2       !Find the midpoint of the two indices
!            WRITE(6,*) i,j,n

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp=DetBitLT(List(:,N),iLut(:),NIfD)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for. Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
                    Comp=DetBitLT(List(:,i+1),iLut(:),NIfD)
                    IF(Comp.eq.0) THEN
                        tSuccess=.true.
                        PartInd=i+1
                        RETURN
                    ELSEIF(Comp.eq.1) THEN
!final entry is less than the one we want.
                        tSuccess=.false.
                        PartInd=i+1
                        RETURN
                    ELSE
                        tSuccess=.false.
                        PartInd=i
                        RETURN
                    ENDIF

                ELSEIF(i.eq.MinInd) THEN
                    tSuccess=.false.
                    PartInd=i
                    RETURN
                ELSE
                    i=j
                ENDIF


            ELSEIF(Comp.eq.-1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j=N
            ELSE
!We have failed - exit loop
                i=j
            ENDIF

        enddo

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess=.false.
        PartInd=MAX(MinInd,i-1)

    END SUBROUTINE BinSearchListHPHF

END MODULE HPHFRandExcitMod


