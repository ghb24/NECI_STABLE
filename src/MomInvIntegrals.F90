module MI_integrals
    use constants, only: dp,n_int,sizeof_int
    use MomInv, only: IsBitMomSelfInv, InvertMomDet, ReturnMomAllowedDet, InvertMomBitDet 
    use MomInv, only: IsMomSelfInv, CalcMomAllowedBitDet
    use SystemData, only: NEl, G1, nBasis, tAntisym_MI, Ecore, modk_offdiag
    use DetBitOps, only: DetBitEQ, FindExcitBitDet, FindBitExcitLevel
    use sltcnd_mod, only: sltcnd, sltcnd_excit
    use bit_reps, only: NIfD, NIfTot, NIfDBO, decode_bit_det
    implicit none
    
    
    interface MI_off_diag_helement
        module procedure MI_off_diag_helement_norm
        module procedure MI_off_diag_helement_spawn
    end interface

    contains

    function MI_spawn_sign (nI, nJ, iLutI, iLutJ, ic, ex, &
                                  parity, HElGen) result (hel)
        integer, intent(in) :: nI(nel), nJ(nel), ic, ex(2,2)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: parity
        integer :: iUnused
        HElement_t :: hel
        HElement_t, intent(in) :: HElGen

        hel=HElGen

        ! Avoid warnings
        iUnused = IC; iUnused = ex(1,1); iUnused = nI(1); iUnused = nJ(1)
        iUnused = int(iLutI(0),sizeof_int); iUnused = int(iLutJ(0),sizeof_int)
        iUnused = parity

    end function MI_spawn_sign

    function MI_off_diag_helement_spawn (nI, nJ, iLutI, iLutJ, ic, ex, &
                                           parity, HElGen) result (hel)
        integer, intent(in) :: nI(nel), nJ(nel), ic, ex(2,2)
        integer(kind=n_int), intent(in) :: iLutI(0:NIfTot), iLutJ(0:NIfTot)
        integer, intent(in) :: parity
        integer :: iUnused
        HElement_t :: hel
        HElement_t , intent(in) :: HElGen

        hel = MI_off_diag_helement_norm (nI, nJ, iLutI, iLutJ)

        if (IC /= 0 .and. modk_offdiag) &
            hel = -abs(hel)

        ! Avoid warnings
        iUnused = IC; iUnused = ex(1,1); iUnused = parity

    end function MI_off_diag_helement_spawn

    function MI_off_diag_helement_norm (nI, nJ, iLutnI, iLutnJ) result(hel)

        ! Find the  between two momentum coupled            
        ! determinants (different ones). NI and nJ have to be uniquely 
        ! chosen, so that their spin-coupled determinant will not arise.
        !
        ! In:  nI, nJ         - Determinants to consider
        !      iLutnI, iLutnJ - Bit representations of I,J
        ! Ret: hel          - The calculated matrix element

        integer, intent(in) :: nI(nel), nJ(nel)
        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot), iLutnJ(0:NIfTot)
        HElement_t :: hel

        integer :: nI2(nel)
        integer(kind=n_int) :: iLutnI2(0:NIfTot)
        integer :: ExcitLevel, Ex(2,2)
        HElement_t :: MatEl2
        integer :: parity

        if (DetBitEQ(iLutnI, iLutnJ, NIfDBO)) then
            ! Do not allow a 'diagonal' matrix element. The problem is 
            ! that the HPHF excitation generator can generate the same HPHF 
            ! function. We do not want to allow spawns here.
            hel = (0.0_dp)
            return
        endif
        
        hel = sltcnd (nI, iLutnI, iLutnJ)
        if (IsMomSelfInv(nI,iLutnI)) then
            if(tAntisym_MI) then
                !Self-inverse MI functions should have no couplings if dets antisymmetric
                hel = 0.0_dp
            elseif (.not. IsMomSelfInv(nJ,iLutnJ)) then
                ! SelfInv MI --> Mom-paired MI, <X|H|Y> = 1/sqrt(2) [Hia + Hib]
                ! OTHERWISE Closed shell -> closed shell. Both the alpha and 
                ! beta of the same orbital have been moved to the same new 
                ! orbital. The matrix element is the same as before.
                hel = hel * (sqrt(2.0_dp))
            endif
        else
            if (IsMomSelfInv(nJ,iLutnJ)) then
                if(tAntisym_MI) then
                    !Self-inverse MI functions should have no couplings if dets antisymmetric
                    hel = 0.0_dp
                else
                    ! Mom-paired MI -> SelfInv MI. If one of
                    ! the determinants is connected, then the other is connected 
                    ! with the same IC & matrix element
                    hel = hel * sqrt(2.0_dp)
                endif
            else
                ! Mom-paired MI -> Mom-paired MI. Find the momentum coupled det of nJ.
                call InvertMomBitDet(iLutnI, iLutnI2, nI)
!                call FindExcitBitDetSym(iLutnI, iLutnI2)
                ExcitLevel = FindBitExcitLevel(iLutnI2, ilutnJ, 2)

                !Cross-term present
                if (ExcitLevel.le.2) then
                    
                    !Calculate cross-term matrix element
                    Ex(1,1)=ExcitLevel
                    call GetBitExcitation(iLutnI2,iLutnJ,Ex,parity)
                    call decode_bit_det(nI2,iLutnI2)

                    MatEl2 = sltcnd_excit (nI2, ExcitLevel, Ex, parity)

                    if(tAntisym_MI) then
                        hel = hel - MatEl2
                    else
                        hel = hel + MatEl2
                    endif
                endif
            endif
        endif
    end function MI_off_diag_helement_norm


    function MI_diag_helement(nI,iLutnI) result(hel)

        ! Finds the diagonal HElement for a momentum inverted
        ! function (potentially two momentum coupled determinants)
        integer, intent(in) :: nI(nel) 
        integer(kind=n_int), intent(in) :: iLutnI(0:NIfTot)
        HElement_t :: hel

        integer :: nI2(nel)
        integer(kind=n_int) :: iLutnI2(0:NIfTot)
        integer :: ExcitLevel
        HElement_t :: MatEl2

        hel = sltcnd_excit (nI,0)
        if(.not.IsMomSelfInv(nI,iLutnI)) then
            !Two determinants in function.
            ! See if there is a cross-term
            call InvertMomBitDet(iLutnI,iLutnI2,nI)
            ExcitLevel = FindBitExcitLevel(iLutnI, iLutnI2, 2)
            if (ExcitLevel.le.2) then
                MatEl2 = sltcnd (nI,  iLutnI, iLutnI2)

                if(tAntisym_MI) then
                    !Subtract cross terms
                    hel = hel - MatEl2
                else
                    hel = hel + MatEl2
                endif
            endif
        endif

        hel = hel + ECore

    end function MI_diag_helement

end module MI_integrals
