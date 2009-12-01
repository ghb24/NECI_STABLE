!Find the HElement between two half-projected hartree-fock determinants (different ones).
!nI and nJ have to be uniquely chosen, so that their spin-coupled determinant will not arise.
SUBROUTINE HPHFGetOffDiagHElement(nI,nJ,iLutnI,iLutnJ,MatEl)
    Use HElem
    Use SystemData , only : NEl,nBasisMax,G1,nBasis,Brr,NIftot,NIfDBO
    use SystemData, only : ECore,ALat,NMSH
    use IntegralsData, only : UMat,FCK,NMAX
    use HPHFRandExcitMod , only : FindDetSpinSym,FindExcitBitDetSym
    use DetBitOps, only: DetBitEQ
    IMPLICIT NONE
    INTEGER :: iLutnI(0:NIfTot),iLutnJ(0:NIfTot),nI(NEl),nI2(NEl),nJ(NEl),nJ2(NEl),iLutnI2(0:NIfTot),iLutnJ2(0:NIfTot)
    INTEGER :: ExcitLevel,OpenOrbsI,OpenOrbsJ,Ex(2,2)
    TYPE(HElement) :: MatEl,MatEl2
    LOGICAL :: TestClosedShellDet,tSymmetricInts,tSign

    MatEl%v=0.D0

!    CALL EncodeBitDet(nI,iLutnI)
!    CALL EncodeBitDet(nJ,iLutnJ)
    IF(DetBitEQ(iLutnI,iLutnJ,NIfDBO)) THEN
!Do not allow an 'off-diagonal' matrix element. The problem is that the HPHF excitation generator can generate the same HPHF function. We do not want to allow spawns here.
        RETURN
    ENDIF
    

    IF(TestClosedShellDet(iLutnI)) THEN
        IF(TestClosedShellDet(iLutnJ)) THEN
!Closed Shell -> Closed Shell. Both alpha and beta of the same orbital have been moved to the same new orbital. The matrix element is the same as before.
!            WRITE(6,*) "Closed Shell -> Closed Shell"
!            CALL FLUSH(6)
            CALL FindBitExcitLevel(iLutnI,iLutnJ,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                Ex(1,1)=ExcitLevel
                CALL GetExcitation(nI,nJ,NEl,Ex,tSign)
!                CALL GetBitExcitation(iLutnI,iLutnJ,Ex,tSign)
!                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart9")
                CALL SltCndExcit2(nEl,nBasisMax,nBasis,nI,nJ,G1,nEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl,Ex,tSign)
!                CALL SltCnd(nEl,nBasisMax,nBasis,nI,nJ,G1,nEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl)
            ENDIF
!            WRITE(6,*) "1 ",MatEl%v
            RETURN
        ELSE
!Closed Shell -> Open Shell. <X|H|Y>= 1/sqrt(2) [Hia + Hib], or with a minus if iLutnJ has an odd number of open orbitals.
!            WRITE(6,*) "Closed Shell -> Open Shell"
!            CALL FLUSH(6)
!First, find <nI|H|nJ>
            CALL FindBitExcitLevel(iLutnI,iLutnJ,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
                Ex(1,1)=ExcitLevel
                CALL GetExcitation(nI,nJ,NEl,Ex,tSign)
!                CALL GetBitExcitation(iLutnI,iLutnJ,Ex,tSign)
                CALL SltCndExcit2(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl,Ex,tSign)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl,Ex,tSign)
                MatEl%v=SQRT(2.D0)*MatEl%v
            ENDIF


!Below is the old way of doing this, however, the integrals are the same (for the same reason as os -> cs)
!            CALL CalcOpenOrbs(iLutnJ,OpenOrbsJ)
!
!            CALL FindDetSpinSym(nJ,nJ2,NEl)
!            CALL FindExcitBitDetSym(iLutnJ,iLutnJ2)
!
!!First, find <nI|H|nJ>
!            CALL FindBitExcitLevel(iLutnI,iLutnJ,ExcitLevel,2)
!            IF(ExcitLevel.le.2) THEN
!                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart8")
!                MatEl2=HElement(0.D0)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!                MatEl=MatEl+MatEl2
!                WRITE(6,*) 1,MatEl2%v,ExcitLevel
!            ENDIF
!!Now, find <nI|H|nJ2>
!            CALL FindBitExcitLevel(iLutnI,iLutnJ2,ExcitLevel,2)
!            IF(ExcitLevel.le.2) THEN
!                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart7")
!                MatEl2=HElement(0.D0)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!                IF(mod(OpenOrbsJ,2).eq.0) THEN
!!Open-shell HPHF is a symmetric combination - therefore, add contribution.
!                    MatEl=MatEl+MatEl2
!                    WRITE(6,*) 2,MatEl2%v,ExcitLevel
!                ELSE
!                    MatEl=MatEl-MatEl2
!                    WRITE(6,*) 2,-1.D0*MatEl2%v,ExcitLevel
!                ENDIF
!            ENDIF
!            MatEl%v=MatEl%v/SQRT(2.D0)
!            WRITE(6,*) "2 ",MatEl%v
            RETURN
        ENDIF
    ELSE
        IF(TestClosedShellDet(iLutnJ)) THEN
!            WRITE(6,*) "Open Shell -> Closed Shell"
!            CALL FLUSH(6)
!            WRITE(6,*) "***"
!OpenShell -> Closed Shell. I am pretty sure that if one of the determinants is connected, then the other is connected with the same IC (+matrix element?) Test this later.
            CALL FindBitExcitLevel(iLutnI,iLutnJ,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
!                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart6")
                Ex(1,1)=ExcitLevel
                CALL GetExcitation(nI,nJ,NEl,Ex,tSign)
!                CALL GetBitExcitation(iLutnI,iLutnJ,Ex,tSign)
                CALL SltCndExcit2(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl,Ex,tSign)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl)
!                WRITE(6,*) 1,MatEl2%v,Excitlevel
!                WRITE(6,*) 2,MatEl2%v,Excitlevel
                MatEl%v=SQRT(2.D0)*MatEl%v
!                WRITE(6,*) "3 ",MatEl%v,nI(:),nJ(:)
            ENDIF


!This is the old way of doing it, however, both integrals have the same absolute value, and the sign from the sym/anti sym cancels the sign change of the two determinants.
!            CALL FindBitExcitLevel(iLutnI,iLutnJ,ExcitLevel,2)
!            IF(ExcitLevel.le.2) THEN
!                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart6")
!                MatEl2=HElement(0.D0)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!                WRITE(6,*) 1,MatEl2%v,Excitlevel
!                MatEl=MatEl+MatEl2
!            ENDIF
!            CALL FindBitExcitLevel(iLutnI2,ilutnJ,ExcitLevel,2)
!            IF(ExcitLevel.le.2) THEN
!                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart5")
!                MatEl2=HElement(0.D0)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI2,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!                IF(mod(OpenOrbsI,2).eq.0) THEN
!!Open-shell HPHF is a symmetric combination - add integral contribution.
!                    MatEl=MatEl+MatEl2
!                    WRITE(6,*) 2,MatEl2%v,ExcitLevel
!                ELSE
!                    MatEl=MatEl-MatEl2
!                    WRITE(6,*) 2,-MatEl2%v,ExcitLevel
!                ENDIF
!            ENDIF
!            MatEl%v=MatEl%v/SQRT(2.D0)
            RETURN

        ELSE
!OpenShell -> Open Shell. Find the spin pair of nJ.
!Initial HPHF is open-shell. Find the spin-coupled determinant.
            CALL FindExcitBitDetSym(iLutnI,iLutnI2)

!            WRITE(6,*) "Open Shell -> Open Shell"
!            CALL FLUSH(6)
!            WRITE(6,*) "***"

            CALL FindBitExcitLevel(iLutnI,iLutnJ,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
!                IF(ExcitLevel.le.0) THEN
!                    WRITE(6,*) ExcitLevel,iLutnI,iLutnJ
!                    CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart4")
!                ENDIF
!                MatEl2%v=0.D0
                Ex(1,1)=ExcitLevel
                CALL GetExcitation(nI,nJ,NEl,Ex,tSign)
!                CALL GetBitExcitation(iLutnI,iLutnJ,Ex,tSign)
                CALL SltCndExcit2(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl,Ex,tSign)
!                WRITE(6,*) "MatEl Old: ",MatEl
!                CALL FLUSH(6)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl)
!                WRITE(6,*) 1,REAL(MatEl2%v,8),ExcitLevel
!                WRITE(6,*) 4,REAL(MatEl2%v,8),ExcitLevel
!                MatEl=MatEl+MatEl2
               
            ENDIF
            
            CALL FindBitExcitLevel(iLutnI2,ilutnJ,ExcitLevel,2)
            IF(ExcitLevel.le.2) THEN
!We need to find out whether the nJ HPHF wavefunction is symmetric or antisymmetric. This is dependant on the number of open shell orbitals.
                CALL FindDetSpinSym(nI,nI2,NEl)
                CALL CalcOpenOrbs(iLutnJ,OpenOrbsJ)
!Original HPHF is antisymmetric if OpenOrbs is odd, or symmetric if its even.
                CALL CalcOpenOrbs(iLutnI,OpenOrbsI)
!                CALL GetBitExcitation(iLutnI2,iLutnJ,Ex,tSign)
                Ex(1,1)=ExcitLevel
                CALL GetExcitation(nI2,nJ,NEl,Ex,tSign)

!                IF((mod(OpenOrbsI,2).eq.0).and.(mod(OpenOrbsJ,2).eq.0)) THEN
!                    tSymmetricInts=.true.
!                ELSE
!                    tSymmetricInts=.false.
!                ENDIF
!                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart3")
                MatEl2%v=0.D0
                CALL SltCndExcit2(NEl,nBasisMax,nBasis,nI2,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2,Ex,tSign)
!                WRITE(6,*) "MatEl2 Old: ",MatEl2
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI2,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
                IF(((mod(OpenOrbsI,2).eq.0).and.(mod(OpenOrbsJ,2).eq.0)).or.((mod(OpenOrbsI,2).eq.0).and.(mod(OpenOrbsJ,2).eq.1))) THEN
!                    WRITE(6,*) 2,REAL(MatEl2%v,8),ExcitLevel
!                    WRITE(6,*) 3,REAL(MatEl2%v,8),ExcitLevel
                    MatEl=MatEl+MatEl2
                ELSE
!                    WRITE(6,*) 2,-1.D0*REAL(MatEl2%v,8),ExcitLevel
!                    WRITE(6,*) 3,-1.D0*REAL(MatEl2%v,8),ExcitLevel
                    MatEl=MatEl-MatEl2
                ENDIF
            ENDIF
        
!Below is the old method for calculating these integrals.

!            CALL FindDetSpinSym(nJ,nJ2,NEl)
!            CALL FindExcitBitDetSym(iLutnJ,iLutnJ2)
!            CALL CalcOpenOrbs(iLutnJ,OpenOrbsJ)
!            IF((mod(OpenOrbsI,2).eq.0).and.(mod(OpenOrbsJ,2).eq.0)) THEN
!                tSymmetricInts=.true.
!            ELSE
!                tSymmetricInts=.false.
!            ENDIF
!
!
!
!!Matrix element is 1/2 [Hia + Hib + Hja + Hjb] when both HPHF functions are symmetric
!            CALL FindBitExcitLevel(iLutnI,iLutnJ,ExcitLevel,2)
!            IF(ExcitLevel.le.2) THEN
!                IF(ExcitLevel.le.0) THEN
!                    WRITE(6,*) ExcitLevel,iLutnI,iLutnJ
!                    CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart4")
!                ENDIF
!                MatEl2=HElement(0.D0)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!                WRITE(6,*) 1,REAL(MatEl2%v,8),ExcitLevel
!                MatEl=MatEl+MatEl2
!            ENDIF
!            CALL FindBitExcitLevel(iLutnI2,ilutnJ,ExcitLevel,2)
!            IF(ExcitLevel.le.2) THEN
!                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart3")
!                MatEl2=HElement(0.D0)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI2,nJ,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!                IF(tSymmetricInts.or.((mod(OpenOrbsI,2).eq.0).and.(mod(OpenOrbsJ,2).eq.1))) THEN
!!Positive is Sym -> sym or Sym -> AntiSym
!                    MatEl=MatEl+MatEl2
!                    WRITE(6,*) 2,REAL(MatEl2%v,8),ExcitLevel
!                ELSE
!                    MatEl=MatEl-MatEl2
!                    WRITE(6,*) 2,-1.D0*REAL(MatEl2%v,8),ExcitLevel
!                ENDIF
!            ENDIF
!            CALL FindBitExcitLevel(iLutnI,iLutnJ2,ExcitLevel,2)
!            IF(ExcitLevel.le.2) THEN
!                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart2")
!                MatEl2=HElement(0.D0)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI,nJ2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!                IF(tSymmetricInts.or.((mod(OpenOrbsI,2).eq.1).and.(mod(OpenOrbsJ,2).eq.0))) THEN
!!Positive is sym -> sym or antisym -> sym
!                    WRITE(6,*) 3,REAL(MatEl2%v,8),ExcitLevel
!                    MatEl=MatEl+MatEl2
!                ELSE
!                    WRITE(6,*) 3,-1.D0*REAL(MatEl2%v,8),ExcitLevel
!                    MatEl=MatEl-MatEl2
!                ENDIF
!            ENDIF
!            CALL FindBitExcitLevel(iLutnI2,ilutnJ2,ExcitLevel,2)
!            IF(ExcitLevel.le.2) THEN
!                IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetOffDiagHElement","Determinants are a forbidden excitation level apart1")
!                MatEl2=HElement(0.D0)
!                CALL SltCnd(NEl,nBasisMax,nBasis,nI2,nJ2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!                IF(tSymmetricInts.or.((mod(OpenOrbsI,2).eq.1).and.(mod(OpenOrbsJ,2).eq.1))) THEN
!!Positive is sym -> sym or antisym -> antisym
!                    WRITE(6,*) 4,REAL(MatEl2%v,8),ExcitLevel
!                    MatEl=MatEl+MatEl2
!                ELSE
!                    WRITE(6,*) 4,-1.D0*REAL(MatEl2%v,8),ExcitLevel
!                    MatEl=MatEl-MatEl2
!                ENDIF
!            ENDIF
!            MatEl%v=MatEl%v/2.D0
!            WRITE(6,*) MatEl%v
!            WRITE(6,*) "4 ",MatEl%v
        ENDIF

    ENDIF

END SUBROUTINE HPHFGetOffDiagHElement


SUBROUTINE HPHFGetDiagHElement(nI,iLutnI,MatEl)
    Use HElem
    Use SystemData , only : NEl,nBasisMax,G1,nBasis,Brr
    use SystemData, only : ECore,ALat,NMSH, NIfTot
    use IntegralsData, only : UMat,FCK,NMAX
    use HPHFRandExcitMod , only : FindDetSpinSym,FindExcitBitDetSym
    IMPLICIT NONE
    INTEGER :: nI(NEl),nI2(NEl),ExcitLevel,OpenOrbs
    INTEGER :: iLutnI(0:NIfTot),iLutnI2(0:NIfTot)
    TYPE(HElement) :: MatEl,MatEl2
    LOGICAL :: TestClosedShellDet

    MatEl%v=0.D0

!    CALL EncodeBitDet(nI,iLutnI)
    IF(TestClosedShellDet(iLutnI)) THEN
        CALL SltCnd(nEl,nBasisMax,nBasis,nI,nI,G1,nEl,NMSH,FCK,NMAX,ALAT,UMat,MatEl)
        MatEl%v=MatEl%v+ECore
        RETURN
    ELSE
!<i|H|i> = <j|H|j>, so no need to calculate both.
        MatEl2%v=0.D0
!<X|H|X> = 1/2 [ <i|H|i> + <j|H|j> ] + <i|H|j> where i and j are the two spin-coupled dets which make up X
!In the case of the antisymmetric pair, the cross term is subtracted.
        CALL SltCnd(nEl,nBasisMax,nBasis,nI,nI,G1,nEl,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
        MatEl%v=MatEl%v+MatEl2%v

!See if there is a cross-term
        CALL FindExcitBitDetSym(iLutnI,iLutnI2)
        CALL FindBitExcitLevel(iLutnI,iLutnI2,ExcitLevel,2)
        IF(ExcitLevel.le.2) THEN
            MatEl2%v=0.D0
            CALL CalcOpenOrbs(iLutnI,OpenOrbs)
            CALL FindDetSpinSym(nI,nI2,NEl)
!            IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetDiagHElement","Determinants are a forbidden excitation level apart")
            CALL SltCnd(NEl,nBasisMax,nBasis,nI,nI2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
            IF(mod(OpenOrbs,2).eq.1) THEN
!Subtract cross terms if determinant is antisymmetric.
                MatEl=MatEl-MatEl2
            ELSE
                MatEl=MatEl+MatEl2
            ENDIF
        ENDIF
        MatEl%v=MatEl%v+ECore


!Below is the old way of doing it - extra effort...
!        CALL CalcOpenOrbs(iLutnI,OpenOrbs)
!!Open Shell Determinant. Find the spin pair
!        CALL FindDetSpinSym(nI,nI2,NEl)
!        CALL FindExcitBitDetSym(iLutnI,iLutnI2)
!
!        MatEl2=HElement(0.D0)
!!<X|H|X> = 1/2 [ <i|H|i> + <j|H|j> ] + <i|H|j> where i and j are the two spin-coupled dets which make up X
!!In the case of the antisymmetric pair, the cross term is subtracted.
!        CALL SltCnd(nEl,nBasisMax,nBasis,nI,nI,G1,nEl,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!        MatEl=MatEl+MatEl2
!        WRITE(6,*) 1,MatEl2%v
!!        WRITE(6,*) MatEl2%v,ECore
!        MatEl2=HElement(0.D0)
!
!        CALL SltCnd(nEl,nBasisMax,nBasis,nI2,nI2,G1,nEl,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!        MatEl=MatEl+MatEl2
!        WRITE(6,*) 2,MatEl2%v
!!        WRITE(6,*) MatEl2%v,ECore
!        MatEl2=HElement(0.D0)
!
!        MatEl%v=MatEl%v/2.D0
!
!!See if they are connected for the 'cross' term
!        CALL FindBitExcitLevel(iLutnI,iLutnI2,ExcitLevel,2)
!        IF(ExcitLevel.le.2) THEN
!            IF(ExcitLevel.le.0) CALL Stop_All("HPHFGetDiagHElement","Determinants are a forbidden excitation level apart")
!            CALL SltCnd(NEl,nBasisMax,nBasis,nI,nI2,G1,NEl-ExcitLevel,NMSH,FCK,NMAX,ALAT,UMat,MatEl2)
!            IF(mod(OpenOrbs,2).eq.1) THEN
!!Subtract cross terms if determinant is antisymmetric.
!                MatEl=MatEl-MatEl2
!                WRITE(6,*) 3,MatEl2%v
!            ELSE
!                MatEl=MatEl+MatEl2
!                WRITE(6,*) 4,MatEl2%v
!            ENDIF
!!            WRITE(6,*) MatEl2
!        ENDIF
!        MatEl%v=MatEl%v+ECore
    ENDIF

END SUBROUTINE HPHFGetDiagHElement
