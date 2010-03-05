module hphf_integrals
    use HElem
    use SystemData, only: NEl, nBasisMax, G1, nBasis, Brr, NIftot, NIfDBO, &
                          tHub, ECore, ALat, NMSH, NIftot
    use IntegralsData, only: UMat,FCK,NMAX
    use HPHFRandExcitMod, only: FindDetSpinSym, FindExcitBitDetSym
    use DetBitOps, only: DetBitEQ, FindExcitBitDet, FindBitExcitLevel
    use sltcnd_mod, only: sltcnd, sltcnd_excit
    use IntegralsData, only: UMat,FCK,NMAX
    implicit none

    
    contains
    function hphf_off_diag_helement (nI, nJ, iLutnI, iLutnJ) result(hel)

        ! Find the HELement between two half-projected hartree-fock 
        ! determinants (different ones). NI and nJ have to be uniquely 
        ! chosen, so that their spin-coupled determinant will not arise.
        !
        ! In:  nI, nJ         - Determinants to consider
        !      iLutnI, iLutnJ - Bit representations of I,J
        ! Ret: hel          - The calculated matrix element

        integer, intent(in) :: nI(nel), nJ(nel)
        integer, intent(in) :: iLutnI(0:NIfTot), iLutnJ(0:NIfTot)
        type(HElement) :: hel

        integer :: nI2(nel), iLutnI2(0:NIfTot)
        integer :: ExcitLevel, OpenOrbsI, OpenOrbsJ, Ex(2,2)
        type(HElement) :: MatEl2
        logical :: TestClosedShellDet, tSymmetricInts, tSign

        if (DetBitEQ(iLutnI, iLutnJ, NIfDBO)) then
            ! Do not allow an 'off-diagonal' matrix element. The problem is 
            ! that the HPHF excitation generator can generate the same HPHF 
            ! function. We do not want to allow spawns here.
            hel = HElement(0)
            return
        endif

        hel = sltcnd (nI, nJ, iLutnI, iLutnJ)
        if (TestClosedShellDet(iLutnI)) then
            if (.not. TestClosedShellDet(iLutnJ)) then
                ! Closed shell --> Open shell, <X|H|Y> = 1/sqrt(2) [Hia + Hib]
                ! or with minus if iLutnJ has an odd number of spin orbitals.
                ! OTHERWISE Closed shell -> closed shell. Both the alpha and 
                ! beta of the same orbital have been moved to the same new 
                ! orbital. The matrix element is the same as before.
                hel = hel * HELement(sqrt(2.d0))
            endif
        else
            if (TestClosedShellDet(iLutnJ)) then
                ! Open shell -> Closed shell. I am pretty sure that if one of
                ! the determinants is connected, then the other is connected 
                ! with the same IC (+matrix element?). Test this later.
                hel = hel * HElement(sqrt(2.d0))
            else
                ! Open shell -> Open shell. Find the spin pair of nJ.
                call FindExcitBitDetSym(iLutnI, iLutnI2)
                ExcitLevel = FindBitExcitLevel(iLutnI2, ilutnJ, 2)

                if (ExcitLevel.le.2) then
                    ! We need to find out whether the nJ HPHF wavefunction is 
                    ! symmetric or antisymmetric. This is dependant on the 
                    ! number of open shell orbitals.
                    call FindDetSpinSym(nI, nI2, nel)
                    call CalcOpenOrbs(iLutnJ, OpenOrbsJ)

                    ! Original HPHF is antispmmetric if OpenOrbs is odd, or
                    ! symmetric if it is even.
                    call CalcOpenOrbs(iLutnI,OpenOrbsI)
                    Ex(1,1)=ExcitLevel
                    call GetBitExcitation(iLutnI2,iLutnJ,Ex,tSign)

                    MatEl2 = sltcnd_excit (nI2, nJ, ExcitLevel, Ex, tSign)

                    if (((mod(OpenOrbsI,2) == 0).and.(mod(OpenOrbsJ,2) == 0))&
                        .or. ((mod(OpenOrbsI,2) == 0) .and. &
                              (mod(OpenOrbsJ,2) == 1))) then
                        hel = hel + MatEl2
                    else
                        hel = hel - MatEl2
                    endif
                endif
            
            endif
        endif
    end function hphf_off_diag_helement


    function hphf_diag_helement (nI, iLutnI) result(hel)

        ! Find the diagonal HElment for a half-projected hartree-fock 
        ! determinant.
        !
        ! In:  nI      - Determinant to consider
        !      iLutnI  - Bit representation of I
        ! Ret: hel   - The calculated matrix element

        integer, intent(in) :: nI(nel), iLutnI(0:NIfTot)
        type(HElement) :: hel

        integer :: nI2(nel), iLutnI2(0:NIfTot)
        integer :: ExcitLevel, OpenOrbs
        type(HElement) :: MatEl2
        logical :: TestClosedShellDet

        hel = sltcnd_excit (nI, nI, 0)
        if (.not. TestClosedShellDet(iLutnI)) then
            ! <i|H|i> = <j|H|j>, so no need to calculate both.
            ! <X|H|X> = 1/2 [ <i|H|i> + <j|H|j> ] + <i|H|j> where i and j are
            ! the two spin-coupled dets which make up X. In the case of the 
            ! antisymmetric pair, the cross term is subtracted.

            ! See if there is a cross-term
            call FindExcitBitDetSym(iLutnI, iLutnI2)
            ExcitLevel = FindBitExcitLevel(iLutnI, iLutnI2, 2)
            if (ExcitLevel.le.2) then
                call CalcOpenOrbs (iLutnI, OpenOrbs)
                call FindDetSpinSym (nI, nI2, nel)
                MatEl2 = sltcnd (nI, nI2, iLutnI, iLutnI2)

                if (mod(OpenOrbs,2).eq.1) then
                    ! Subtract cross terms if determinant is antisymmetric.
                    hel = hel - MatEl2
                else
                    hel = hel + MatEl2
                endif
            endif
        endif

        hel = hel + HELement(ECore)
    end function hphf_diag_helement
end module
