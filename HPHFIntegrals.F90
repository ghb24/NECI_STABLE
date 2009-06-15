!Find the HElement between two half-projected hartree-fock determinants (different ones).
!nI and nJ have to be uniquely chosen, so that their spin-coupled determinant will not arise.
SUBROUTINE HPHFGetOffDiagHElement(nI,nJ,MatEl)
    Use HElem
    Use SystemData , only : NIfD,NEl,nBasisMax,G1,nBasis,Brr
    use SystemData, only : ECore,ALat,NMSH
    use IntegralsData, only : UMat,FCK,NMAX
    use HPHFRandExcitMod , only : FindDetSpinSym,FindExcitBitDetSym
    IMPLICIT NONE
    INTEGER :: iLutnI(0:NIfD),iLutnJ(0:NIfD),nI(NEl),nI2(NEl),nJ(NEl),nJ2(NEl),iLutnI2(0:NIfD),iLutnJ2(0:NIfD)
    INTEGER :: ExcitLevel,OpenOrbsI,OpenOrbsJ,Ex(2,NEl)
    TYPE(HElement) :: MatEl,MatEl2
    LOGICAL :: TestClosedShellDet,tSymmetricInts,DetBitEQ

    MatEl=HElement(0.D0)

    CALL EncodeBitDet(nI,iLutnI,NEl,NIfD)
    CALL EncodeBitDet(nJ,iLutnJ,NEl,NIfD)
    IF(DetBitEQ(iLutnI,iLutnJ,NIfD)) THEN
!Do not allow an 'off-diagonal' matrix element. The problem is that the HPHF excitation generator can generate the same HPHF function. We do not want to allow spawns here.
        RETURN
    ENDIF

    IF(TestClosedShellDet(iLutnI,NIfD)) THEN
        IF(TestClosedShellDet(iLutnJ,NIfD)) THEN
!Closed Shell -> Closed Shell. Both alpha and beta of the same orbital have been moved to the same new orbital. The matrix element is the same as before.
!            WRITE(6,*) "Closed Shell -> Closed Shell"
!            CALL FLUSH(6)
            CALL FindBitExcitLevel(iLutnI,iLutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                IF(ExcitLevel.le.0) THEN
                    WRITE(6,*) iLutnI(:),iLutnJ(:),ExcitLevel
                    WRITE(6,*) "***",nI(:)
                    WRITE(6,*) "***",nJ(:)
                    CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart9")
                ENDIF
                CALL SltCnd(nEl,nBasisMax,nBasis,nI,nJ,G1,nEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl)
            ENDIF
            RETURN
        ELSE
!Closed Shell -> Open Shell. <X|H|Y>= 1/sqrt(2) [Hia + Hib], or with a minus if iLutnJ has an odd number of open orbitals.
!            WRITE(6,*) "Closed Shell -> Open Shell"
!            CALL FLUSH(6)
            CALL CalcOpenOrbs(iLutnJ,NIfD,NEl,OpenOrbsJ)

            CALL FindDetSpinSym(nJ,nJ2,NEl)
            CALL FindExcitBitDetSym(iLutnJ,iLutnJ2)

!First, find <nI|H|nJ>
            CALL FindBitExcitLevel(iLutnI,iLutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart8")
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
            ENDIF
!Now, find <nI|H|nJ2>
            CALL FindBitExcitLevel(iLutnI,iLutnJ2,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart7")
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                IF(mod(OpenOrbsJ,2).eq.0) THEN
!Open-shell HPHF is a symmetric combination - therefore, add contribution.
                    MatEl=MatEl+MatEl2
                ELSE
                    MatEl=MatEl-MatEl2
                ENDIF
            ENDIF
            MatEl%v=MatEl%v/SQRT(2.D0)
            RETURN
        ENDIF
    ELSE
!Initial HPHF is open-shell. Find the spin-coupled determinant.
        CALL FindDetSpinSym(nI,nI2,NEl)
        CALL FindExcitBitDetSym(iLutnI,iLutnI2)

!Original HPHF is antisymmetric if OpenOrbs is odd, or symmetric if its even.
        CALL CalcOpenOrbs(iLutnI,NIfD,NEl,OpenOrbsI)

        IF(TestClosedShellDet(iLutnJ,NIfD)) THEN
!            WRITE(6,*) "Open Shell -> Closed Shell"
!            CALL FLUSH(6)
!OpenShell -> Closed Shell. I am pretty sure that if one of the determinants is connected, then the other is connected with the same IC (+matrix element?) Test this later.
            CALL FindBitExcitLevel(iLutnI,iLutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart6")
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
            ENDIF
            CALL FindBitExcitLevel(iLutnI2,ilutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart5")
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI2,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                IF(mod(OpenOrbsI,2).eq.0) THEN
!Open-shell HPHF is a symmetric combination - add integral contribution.
                    MatEl=MatEl+MatEl2
                ELSE
                    MatEl=MatEl-MatEl2
                ENDIF
            ENDIF
            MatEl%v=MatEl%v/SQRT(2.D0)
            RETURN

        ELSE
!OpenShell -> Open Shell. Find the spin pair of nJ.
!            WRITE(6,*) "Open Shell -> Open Shell"
!            CALL FLUSH(6)
            CALL FindDetSpinSym(nJ,nJ2,NEl)
            CALL FindExcitBitDetSym(iLutnJ,iLutnJ2)
        
!We need to find out whether the nJ HPHF wavefunction is symmetric or antisymmetric. This is dependant on the number of open shell orbitals.
            CALL CalcOpenOrbs(iLutnJ,NIfD,NEl,OpenOrbsJ)
            IF((mod(OpenOrbsI,2).eq.0).and.(mod(OpenOrbsJ,2).eq.0)) THEN
                tSymmetricInts=.true.
            ELSE
                tSymmetricInts=.false.
            ENDIF

!Matrix element is 1/2 [Hia + Hib + Hja + Hjb] when both HPHF functions are symmetric
            CALL FindBitExcitLevel(iLutnI,iLutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                IF(ExcitLevel.le.0) THEN
                    WRITE(6,*) ExcitLevel,iLutnI,iLutnJ
                    CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart4")
                ENDIF
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
            ENDIF
            CALL FindBitExcitLevel(iLutnI2,ilutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart3")
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI2,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                IF(tSymmetricInts.or.((mod(OpenOrbsI,2).eq.0).and.(mod(OpenOrbsJ,2).eq.1))) THEN
!Positive is Sym -> sym or Sym -> AntiSym
                    MatEl=MatEl+MatEl2
                ELSE
                    MatEl=MatEl-MatEl2
                ENDIF
            ENDIF
            CALL FindBitExcitLevel(iLutnI,iLutnJ2,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart2")
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                IF(tSymmetricInts.or.((mod(OpenOrbsI,2).eq.1).and.(mod(OpenOrbsJ,2).eq.0))) THEN
!Positive is sym -> sym or antisym -> sym
                    MatEl=MatEl+MatEl2
                ELSE
                    MatEl=MatEl-MatEl2
                ENDIF
            ENDIF
            CALL FindBitExcitLevel(iLutnI2,ilutnJ2,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart1")
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI2,nJ2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                IF(tSymmetricInts.or.((mod(OpenOrbsI,2).eq.1).and.(mod(OpenOrbsJ,2).eq.1))) THEN
!Positive is sym -> sym or antisym -> antisym
                    MatEl=MatEl+MatEl2
                ELSE
                    MatEl=MatEl-MatEl2
                ENDIF
            ENDIF
            MatEl%v=MatEl%v/2.D0
        ENDIF

    ENDIF

END SUBROUTINE HPHFGetOffDiagHElement


SUBROUTINE HPHFGetDiagHElement(nI,MatEl)
    Use HElem
    Use SystemData , only : NIfD,NEl,nBasisMax,G1,nBasis,Brr
    use SystemData, only : ECore,ALat,NMSH
    use IntegralsData, only : UMat,FCK,NMAX
    use HPHFRandExcitMod , only : FindDetSpinSym,FindExcitBitDetSym
    IMPLICIT NONE
    INTEGER :: nI(NEl),nI2(NEl),ExcitLevel,OpenOrbs
    INTEGER :: iLutnI(0:NIfD),iLutnI2(0:NIfD)
    TYPE(HElement) :: MatEl,MatEl2
    LOGICAL :: TestClosedShellDet

    MatEl=HElement(0.D0)

    CALL EncodeBitDet(nI,iLutnI,NEl,NIfD)
    IF(TestClosedShellDet(iLutnI,NIfD)) THEN
        CALL SltCnd(nEl,nBasisMax,nBasis,nI,nI,G1,nEl,NMSH,FCK,NMAX,ALAT,UMat,MatEl)
        MatEl%v=MatEl%v+ECore
        RETURN
    ELSE
        CALL CalcOpenOrbs(iLutnI,NIfD,NEl,OpenOrbs)
!Open Shell Determinant. Find the spin pair
        CALL FindDetSpinSym(nI,nI2,NEl)
        CALL FindExcitBitDetSym(iLutnI,iLutnI2)

        MatEl2=HElement(0.D0)
!<X|H|X> = 1/2 [ <i|H|i> + <j|H|j> ] + <i|H|j> where i and j are the two spin-coupled dets which make up X
!In the case of the antisymmetric pair, the cross term is subtracted.
        CALL SltCnd(nEl,nBasisMax,nBasis,nI,nI,G1,nEl,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
        MatEl=MatEl+MatEl2
!        WRITE(6,*) MatEl2%v,ECore
        MatEl2=HElement(0.D0)

        CALL SltCnd(nEl,nBasisMax,nBasis,nI2,nI2,G1,nEl,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
        MatEl=MatEl+MatEl2
!        WRITE(6,*) MatEl2%v,ECore
        MatEl2=HElement(0.D0)

        MatEl%v=MatEl%v/2.D0

!See if they are connected for the 'cross' term
        CALL FindBitExcitLevel(iLutnI,iLutnI2,NIfD,ExcitLevel,2)
        IF(ExcitLevel.le.2) THEN
            IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetDiagHElement","Determinants are a forbidden excitation level apart")
            CALL SltCnd(NEl,nBasisMax,nBasis,nI,nI2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
            IF(mod(OpenOrbs,2).eq.1) THEN
!Subtract cross terms if determinant is antisymmetric.
                MatEl=MatEl-MatEl2
            ELSE
                MatEl=MatEl+MatEl2
            ENDIF
!            WRITE(6,*) MatEl2
        ENDIF
        MatEl%v=MatEl%v+ECore
    ENDIF

END SUBROUTINE HPHFGetDiagHElement
