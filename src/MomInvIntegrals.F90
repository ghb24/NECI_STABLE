module MI_integrals
    use constants, only: dp,n_int
    use MomInv, only: IsBitMomSelfInv, InvertMomDet, ReturnMomAllowedDet, InvertMomBitDet 
    use MomInv, only: IsMomSelfInv, CalcMomAllowedBitDet
    use SystemData, only: NEl, G1, nBasis, tAntisym_MI
    use DetBitOps, only: DetBitEQ, FindExcitBitDet, FindBitExcitLevel
    use sltcnd_mod, only: sltcnd, sltcnd_excit
    use bit_reps, only: NIfD, NIfTot, NIfDBO, decode_bit_det

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
        if(.not.IsMomSelfInv(nI,iLutnI) then
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
