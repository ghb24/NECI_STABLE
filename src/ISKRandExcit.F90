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
    use constants, only: dp,n_int
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

        !First, generate a random excitation from the determinant which is given in the argument
        call gen_rand_excit (nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                             tSignOrig, pGen, HEl, tFilled, Classcount2, &
                             ClassCountUnocc2, scratch)
         if(IsNullDet(nJ)) return    !No excitation created

!Create bit representation of excitation - iLutnJ
        call FindExcitBitDet(iLutnI,iLutnJ,IC,ExcitMat)


    end subroutine gen_ISK_excit

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

END MODULE ISKRandExcit 
