#include "macros.h"

MODULE ISKRandExcit 
!ISK (Inversion-Symmetry K-point) wavefunctions are a linear combinations of two determinants, 
!where all orbitals are exchanged with 
!the equivalent orbital where the k-points of all occupied orbitals have been switched.
!In simple notation we will consider an excitation where determinants i and j are in the original ISK,
!and determinants a and b in the excited ISK.
!As with HPHF wavefunctions, a large simplification occurs when it is realised that P(i->a) = P(j->b),
!therefore all excitation are connected to the determinantal excitations of just one of the constituent determinants.
!This means that only one constituent determinant will be considered in the space.

    use SystemData, only: nel, tCSF, Alat, G1, nbasis, nbasismax, nmsh, arr
    use IntegralsData, only: UMat, fck, nMax
    use SymData, only: nSymLabels
    use dSFMT_interface, only : genrand_real2_dSFMT
    use GenRandSymExcitNUMod, only: gen_rand_excit, ConstructClassCounts, &
                                    CalcNonUniPGen, ScratchSize 
    use DetBitOps, only: DetBitLT, DetBitEQ, FindExcitBitDet, &
                         FindBitExcitLevel,MaskAlpha,MaskBeta
    use FciMCData, only: pDoubles
    use constants, only: dp,n_int,bits_n_int
    use HElem
    use sltcnd_mod, only: sltcnd_excit
    use bit_reps, only: NIfD, NIfDBO, NIfTot
    use sort_mod
    IMPLICIT NONE

    contains

    subroutine gen_ISK_excit (nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                               tParity, pGen, HEl, tFilled, ClassCount2, &
                               ClassCountUnocc2, scratch)
        use FciMCData, only: tGenMatHEl

        integer, intent(in) :: nI(nel) 
        integer(kind=n_int), intent(in) :: iLutnI(0:niftot)
        integer, intent(in) :: exFlag
        integer, intent(out) :: nJ(nel)
        integer(kind=n_int), intent(out) :: iLutnJ(0:niftot)
        integer, intent(out) :: IC, ExcitMat(2,2)
        integer, intent(inout) :: ClassCount2(ScratchSize)
        integer, intent(inout) :: ClasscountUnocc2(ScratchSize)
        integer, intent(inout) :: scratch(ScratchSize) ! Not used
        logical, intent(out) :: tParity ! Not used
        logical, intent(inout) :: tFilled
        real*8, intent(out) :: pGen
        HElement_t, intent(out) :: HEl
        character(*), parameter :: this_routine='gen_ISK_excit'

        !First, generate a random excitation from the determinant which is given in the argument
        call gen_rand_excit (nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                             tSignOrig, pGen, HEl, tFilled, Classcount2, &
                             ClassCountUnocc2, scratch)
         if(IsNullDet(nJ)) return    !No excitation created

!Create bit representation of excitation - iLutnJ
        call FindExcitBitDet(iLutnI,iLutnJ,IC,ExcitMat)

        if(is_self_inverse(nJ,iLutnJ)) then
!There is only one way which we could have generated the excitation nJ since it has no determinant partner.
!Also, we will always return the 'correct' determinant.
            if(tGenMatHEl) then
                !Generate matrix element of excitation.
                call stop_all(this_routine,"tGenMatHEl not yet implemented")
            endif
        else
!This excitation *does* have a inverse determinant. Could we have generated this instead?
!Find the inverse
            call returned_invsym(nJ,nJSym,iLutnJ,iLutnJSym,.true.,tSwapped)

!Calculate whether the 'cross-term' is non-zero. i.e. is the original determinant connected to the 
!inverse determinant of the excitation.
            if(tSwapped) then
                !We have swapped the excitation, so that iLutnJ now contains the cross determinant
                call ISK_cross_det_conn(iLutnI,iLutnJ,tSame_ISK,tCrossConnected,CrossIC,CrossEx,tSignCross)
            else
                call ISK_cross_det_conn(iLutnI,iLutnJSym,tSame_ISK,tCrossConnected,CrossIC,CrossEx,tSignCross)
            endif

            if(tSame_ISK) then
                !We have created the same ISK - return null ISK
                nJ(1)=0
                if(tGenHEl) HEl=0.D0
            elseif(tCrossConnected) then
                !The cross-term is connected. Calculate the probability that we created this det instead (in same ISK)
                CALL CalcNonUniPGen(nI,CrossEx,CrossIC,ClassCount2,ClassCountUnocc2,pDoubles,pGen2)
                pGen=pGen+pGen2

                if(tGenMatHEl) then
                    !calculate matrix element to open shell excitation.
                endif
            else
                !The cross-term is not connected. Therefore, there was only one way we could have generated this ISK.
                !nI *CANNOT* be a self-inverse ISK here, otherwise we should have had a cross-term.
                !Check this, but remove check once we know it is working.
                if(is_self_inverse(nI,iLutnI)) call stop_all(this_routine,"Dodgy logic - this should not be self-inv")
                
                if(tGenMatEl) then
                endif
            endif
        endif

    end subroutine gen_ISK_excit


    !determine whether two determinants are connected.
    !Return whether the ISKs are actually the same
    !If the cross term is found to exist, return the excitation matrix and parity of the excitation
    !TODO: Optimisation - do we need to calculate the excitation matrix, or can we just pass it in?
    pure subroutine ISK_cross_det_conn(iLutnI,iLutnJSym,tSame_ISK,tcross_conn,CrossIC,ExCross,tExSign)
        integer(n_int), intent(in) :: iLutnI(0:NIfTot),iLutnJSym(0:NIfTot)
        logical, intent(out) :: tcross_conn,tSame_ISK,tExSign
        integer, intent(out) :: ExCross(2,2),CrossIC
        integer :: SymB,ijSymProd

        tSame_ISK=.false.   !Is the ISK generated actually the one we started with?

        !Comment out these lines when we are sure above call is working!
        CrossIC = FindBitExcitLevel(iLutnI, iLutnJSym, 2)

        if(CrossIC.eq.0) then
            !We have generated the same ISK!!
            !The inverse of the excitation is the determinant that we started with!
            !This does not want to be allowed.
            tSame_ISK=.true.
            return
        elseif(CrossIC.gt.2) then
            tcross_conn=.false.
            return
        endif

        !cross-term is allowed by excitation level. Is it allowed by momentum conservation?
        !calculate excitation matrix
        ExCross(1,1)=CrossIC
        call GetBitExcitation(iLutnI,iLutnJSym,ExCross,tExSign)
        !i = ex(1,1)
        !j = ex(1,2)
        if(CrossIC.eq.1) then
            if(SpinOrbSymLabel(ExCross(1,1)).ne.SpinOrbSymLabel(ExCross(2,1))) then
                !symmetry forbidden
                tcross_conn=.false.
            else
                tcross_conn=.true.
            endif
        else
            !Excitlevel *must* equal 2
            ijSymProd=SymTableLabels(SpinOrbSymLabel(ExCross(1,1)),SpinOrbSymLabel(ExCross(1,2)))
            SymB=SymTableLabels(ijSymProd,SymInvLabel(ExCross(2,1)))
            if(SymB.ne.SpinOrbSymLabel(ExCross(2,2))) then
                !symmetry forbidden
                tcross_conn=.false.
            else
                tcross_conn=.true.
            endif
        endif

    end subroutine ISK_cross_det_conn

!Function to determine whether a determinant is its own self inverse when inverted.
!If it is (returns true), then there is no other determinant in the ISK function.
!Input: the determinant to test, in both natural ordered, and bit representations.
!TODO: This could be sped up by a factor of two, since we only actually need to 
!search through half of the electrons, since if it is symmetric, then the other 
!half should be inverses of ones that we've already tested!
    pure function is_self_inverse(nI,iLut) result(tSelfInv)
        integer, intent(in) :: nI(Nel)
        intent(n_int), intent(in) :: iLut(0:niftot)
        logical, intent(out) :: tSelfInv
        integer :: pos

        tSelfInv=.true.

        do i=1,NEl
            !run through the electrons, find inverse and test whether it is in the bit representation of the determinant.
            pos = KPntInvSymOrb(nI(i))-1

            if(.not. (btest(iLut(pos/bits_n_int),mod(pos,bits_n_int)))) then
                !the inverse orbital is not found in the determinant.
                !therefore the determinant can not be a self inverse.
                tSelfInv=.false.
                return
            endif
        enddo

    end function is_self_inverse

!This routine will take a determinant (nI and nISym) and return in the same place the correct
!determinant out of the two inverse determinants which make up the ISK. If this is actually its
!partner determinant, then tSwapped is returned as true, and nISym and iLutnISym are returned
!as the original determinant.
!The flag tCalcnISym let the routine know whether to calculate the symmetry
!partner of nI initially, or whether it is passed in correctly.
!Whether one determinant or its partner is returned as nI is purely based on which bit representation is 'largest'.
    subroutine returned_invsym(nI,nISym,iLutnI,iLutnISym,tCalcnISym,tSwapped)
        integer, intent(inout) :: nI(NEl),nISym(NEl)
        integer :: nTemp(NEl)
        intent(n_int) , intent(inout) :: iLutnI(0:NIfTot),iLutnISym(0:NIfTot)
        intent(n_int) :: iLutTemp(0:NIfTot)
        logical, intent(in) :: tCalcnISym
        logical, intent(out) :: tSwapped

        if(tCalcnISym) then
            call find_invsym_det(nI,nISym,iLutnISym)
        endif

        ! iLutnI is 'less' than iLutSym, so iLutSym is the determinant with 
        ! the first open-shell = alpha. Swap them around.
        i=DetBitLT(iLutnI,iLutSym,NIfD)
        if(i.eq.1) then
            iLutTemp(:)=iLutnI(:)
            iLutnI(:)=iLutnISym(:)
            iLutnISym(:)=iLutTemp(:)
            nTemp(:)=nI(:)
            nI(:)=nISym(:)
            nISym(:)=nTemp(:)
            tSwapped=.true.
        elseif(i.eq.0) then
            CALL Stop_All("ReturnAlphaOpenDet","Shouldn't have self-inverse determinants in here")
        else
            tSwapped=.false.
        endif

    end subroutine returned_invsym

    subroutine find_invsym_det(nI,nISym,iLutnISym)
        integer , intent(in) :: nI(NEl)
        integer(n_int) , intent(out) :: iLutnISym(0:NIfTot)
        integer , intent(out) :: nISym(NEl)
        integer :: orb,pos

        iLutnISym(:)=0

        do i=1,NEl
            orb=KPntInvSymOrb(nI(i))
            nISym(i)=orb
            pos=(orb-1)/bits_n_int
            iLutnISym(pos)=ibset(iLutnISym(pos),mod(orb-1,bits_n_int))
        enddo

    end subroutine find_invsym_det



END MODULE ISKRandExcit 
