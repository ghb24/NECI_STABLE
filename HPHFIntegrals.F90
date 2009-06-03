!Find the HElement between two half-projected hartree-fock determinants (different ones).
SUBROUTINE HPHFGetOffDiagHElement(nI,nJ,MatEl)
    Use HElem
    Use SystemData , only : NIfD,NEl,nBasisMax,G1,nBasis,Brr
    use SystemData, only : ECore,ALat,NMSH
    use IntegralsData, only : UMat,FCK,NMAX
    use HPHFRandExcitMod , only : FindDetSpinSym,FindExcitBitDetSym
    IMPLICIT NONE
    INTEGER :: iLutnI(0:NIfD),iLutnJ(0:NIfD),nI(NEl),nI2(NEl),nJ(NEl),nJ2(NEl),iLutnI2(0:NIfD),iLutnJ2(0:NIfD)
    INTEGER :: ExcitLevel
    TYPE(HElement) :: MatEl,MatEl2
    LOGICAL :: TestClosedShellDet

    MatEl=HElement(0.D0)

    CALL EncodeBitDet(nI,iLutnI,NEl,NIfD)
    CALL EncodeBitDet(nJ,iLutnJ,NEl,NIfD)

    IF(TestClosedShellDet(iLutnI,NIfD)) THEN
        IF(TestClosedShellDet(iLutnJ)) THEN
!Closed Shell -> Closed Shell. Both alpha and beta of the same orbital have been moved to the same new orbital. The matrix element is the same as before.
            CALL SltCnd(nEl,nBasisMax,nBasis,nI,nJ,G1,nEl-2,NMSH,FCK,NMAX,ALAT,UMat,MatEl)
            RETURN
        ELSE
!Closed Shell -> Open Shell. <X|H|Y>= 1/sqrt(2) [Hia + Hib]
            CALL FindDetSpinSym(nJ,nJ2,NEl)
            CALL FindExcitBitDetSym(iLutnJ,iLutnJ2,NIfD)

!First, find <nI|H|nJ>
            CALL FindBitExcitLevel(iLutnI,iLutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
            ENDIF
!Now, find <nI|H|nJ2>
            CALL FindBitExcitLevel(iLutnI,iLutnJ2,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
            ENDIF
            MatEl%v=MatEl%v/SQRT(2.D0)
            RETURN
        ENDIF
    ELSE
!Initial HPHF is open-shell. Find the spin-coupled determinant.
        CALL FindDetSpinSym(nI,nI2,NEl)
        CALL FindExcitBitDetSym(iLutnI,iLutnI2,NIfD)

        IF(TestClosedShellDet(iLutnJ)) THEN
!OpenShell -> Closed Shell. I am pretty sure that if one of the determinants is connected, then the other is connected with the same IC (+matrix element?) Test this later.
            CALL FindBitExcitLevel(iLutnI,iLutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
            ENDIF
            CALL FindBitExcitLevel(iLutnI2,ilutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI2,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
            ENDIF
            MatEl%v=MatEl%v/SQRT(2.D0)
            RETURN
        ELSE
!OpenShell -> Open Shell. Find the spin pair of nJ.
            CALL FindDetSpinSym(nJ,nJ2,NEl)
            CALL FindExcitBitDetSym(iLutnJ,iLutnJ2,NIfD)

!Matrix element is 1/2 [Hia + Hib + Hja + Hjb]
            CALL FindBitExcitLevel(iLutnI,iLutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
            ENDIF
            CALL FindBitExcitLevel(iLutnI2,ilutnJ,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI2,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
            ENDIF
            CALL FindBitExcitLevel(iLutnI,iLutnJ2,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
            ENDIF
            CALL FindBitExcitLevel(iLutnI2,ilutnJ2,NIfD,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                MatEl2=HElement(0.D0)
                CALL SltCnd(NEl,nBasisMax,nBasis,nI2,nJ2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                MatEl=MatEl+MatEl2
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
    INTEGER :: nI(NEl),nI2(NEl),ExcitLevel
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
!Open Shell Determinant. Find the spin pair
        CALL FindDetSpinSym(nI,nI2,NEl)
        CALL FindExcitBitDetSym(iLutnI,iLutnI2,NIfD)

        MatEl2=HElement(0.D0)
!<X|H|X> = 1/2 [ <i|H|i> + <j|H|j> ] + <i|H|j> where i and j are the two spin-coupled dets which make up X
        CALL SltCnd(nEl,nBasisMax,nBasis,nI,nI,G1,nEl,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
        MatEl=MatEl+MatEl2
        MatEl2=HElement(0.D0)

        CALL SltCnd(nEl,nBasisMax,nBasis,nI2,nI2,G1,nEl,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
        MatEl=MatEl+MatEl2
        MatEl2=HElement(0.D0)

        MatEl%v=MatEl%v/2.D0

!See if they are connected for the 'cross' term
        CALL FindBitExcitLevel(iLutnI,iLutnI2,NIfD,ExcitLevel,2)
        IF(ExcitLevel.le.2) THEN
            CALL SltCnd(NEl,nBasisMax,nBasis,nI,nI2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
            MatEl=MatEl+MatEl2
        ENDIF
        MatEl%v=MatEl%v+ECore
    ENDIF

END SUBROUTINE HPHFGetDiagHElement
