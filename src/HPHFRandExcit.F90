#include "macros.h"

MODULE HPHFRandExcitMod
!Half-projected HF wavefunctions are a linear combination of two HF determinants,
!where all alphas -> betas and betas -> alpha to create the pair.
!In closed-shell systems, these two determinants have the same FCI amplitude, and so it is easier to treat them as a pair.
!The probability of creating this HPHF where both are from pairs of spin-coupled determinants (i & j -> a & b):
![ P(i->a) + P(i->b) + P(j->a) + P(j->b) ]/2
!We therefore need to find the excitation matrix between the determinant which wasn't
!excited and the determinant which was created.

    use SystemData, only: nel, tCSF, Alat, G1, nbasis, nbasismax, nmsh, arr, &
                          tOddS_HPHF, modk_offdiag, tGen_4ind_weighted, &
                          tGen_4ind_reverse, tLatticeGens, tGen_4ind_2, tHUB, &
                          tUEG, tUEGNewGenerator, t_pchb_excitgen
    use IntegralsData, only: UMat, fck, nMax
    use SymData, only: nSymLabels
    use dSFMT_interface, only : genrand_real2_dSFMT
    use GenRandSymExcitNUMod, only: gen_rand_excit, calc_pgen_symrandexcit2, &
                                    ScratchSize, CalcPGenLattice
    use excit_gens_int_weighted, only: gen_excit_4ind_weighted, &
                                       gen_excit_4ind_reverse, &
                                       calc_pgen_4ind_weighted, &
                                       calc_pgen_4ind_reverse
    use DetBitOps, only: DetBitLT, DetBitEQ, FindExcitBitDet, &
                         FindBitExcitLevel, MaskAlpha, MaskBeta, &
                         TestClosedShellDet, CalcOpenOrbs, IsAllowedHPHF, &
                         DetBitEQ
    use FciMCData, only: pDoubles, excit_gen_store_type, ilutRef
    use constants, only: dp,n_int, EPS
    use sltcnd_mod, only: sltcnd_excit
    use bit_reps, only: NIfD, NIfDBO, NIfTot
    use SymExcitDataMod, only: excit_gen_store_type
    use excit_gen_5, only: calc_pgen_4ind_weighted2, gen_excit_4ind_weighted2
    use pchb_excitgen, only: calc_pgen_pchb, gen_rand_excit_pchb
    use sort_mod
    use HElem
    use CalcData, only: t_matele_cutoff, matele_cutoff, t_back_spawn, t_back_spawn_flex
    use back_spawn_excit_gen, only: gen_excit_back_spawn, calc_pgen_back_spawn, &
                                    gen_excit_back_spawn_ueg, calc_pgen_back_spawn_ueg, &
                                    calc_pgen_back_spawn_hubbard, gen_excit_back_spawn_hubbard, &
                                    gen_excit_back_spawn_ueg_new, calc_pgen_back_spawn_ueg_new
    IMPLICIT NONE
!    SAVE
!    INTEGER :: Count=0

    contains

!Calculate probability of exciting from HPHF nI to HPHF nJ
!It is imperative that when using this routine, the 'correct' determinant is sent in
!i.e. the unique determinant representation of the two HPHF functions. This is because
!the classcount arrays will be different for the two determinants.
!tSameFunc will be returned as true if the two HPHF functions are the same
    subroutine CalcPGenHPHF (nI,iLutnI,nJ,iLutnJ,ex,ClassCount,ClassCountUnocc,pDoubles,pGen,tSameFunc)
        integer, intent(in) :: nI(nel)
        integer(kind=n_int), intent(in) :: iLutnI(0:niftot),iLutnJ(0:niftot)
        integer, intent(in) :: ClassCount(ScratchSize),ClassCountUnocc(ScratchSize)
        integer, intent(in) :: nJ(nel),ex(2,2)
        real(dp), intent(in) :: pDoubles
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tSameFunc
        logical :: tSign,tSwapped
        real(dp) :: pGen2
        integer :: ic
        integer :: Ex2(2,2),nJ_loc(nel),nJ2(nel)
        integer(kind=n_int) :: iLutnJ_loc(0:niftot),iLutnJ2(0:niftot)

        tSameFunc = .false.
        pGen = 0.0_dp

        IF(TestClosedShellDet(iLutnJ)) THEN
            !nJ is CS, therefore, only one way of generating it.
            ic = FindBitExcitLevel(iLutnI, iLutnJ, 2)
            if(ic.eq.0) then
                tSameFunc=.true.
                return
            endif
            if(ic.le.2) then
                call CalcNonUniPGen(nI,ilutnI,ex,ic,ClassCount,ClassCountUnocc,pDoubles,pGen)
            endif
        else
            !nJ is openshell. Add the probabilities of generating each pair (if both connected)
            nJ_loc = nJ
            iLutnJ_loc = iLutnJ
            CALL ReturnAlphaOpenDet(nJ_loc,nJ2,iLutnJ_loc,iLutnJ2,.true.,.true.,tSwapped)

            !First find nI -> nJ
            ic = FindBitExcitLevel(iLutnI, iLutnJ_loc, 2)
            if(ic.eq.0) then
                tSameFunc=.true.
                return
            endif
            if(ic.le.2) then
                if(.not.tSwapped) then
                    !ex is correct for this excitation
                    call CalcNonUnipGen(nI,ilutnI,ex,ic,ClassCount,ClassCountUnocc,pDoubles,pGen)
                else
                    Ex2(1,1)=ic
                    call GetBitExcitation(iLutnI,iLutnJ_loc,Ex2,tSign)
                    call CalcNonUnipGen(nI,ilutnI,Ex2,ic,ClassCount,ClassCountUnocc,pDoubles,pGen)
                endif
            endif

            !Now consider nI -> nJ2 and add the probabilities
            ic = FindBitExcitLevel(iLutnI, iLutnJ2, 2)
            if(ic.eq.0) then
                tSameFunc=.true.
                return
            endif
            if(ic.le.2) then
                if(tSwapped) then
                    !ex is correct for this excitation
                    call CalcNonUnipGen(nI, ilutnI,ex,ic,ClassCount,ClassCountUnocc,pDoubles,pGen2)
                else
                    Ex2(1,1)=ic
                    call GetBitExcitation(iLutnI,iLutnJ2,Ex2,tSign)
                    call CalcNonUnipGen(nI, ilutnI,Ex2,ic,ClassCount,ClassCountUnocc,pDoubles,pGen2)
                endif
                pGen = pGen + pGen2
            endif
        endif

    end subroutine CalcPGenHPHF

    subroutine gen_hphf_excit (nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                               tParity, pGen, HEl, store, part_type)

        use FciMCData, only: tGenMatHEl

        ! Generate an HPHF excitation using only one of the determinants in
        ! the source HPHF function.
        !
        ! --> Relies on both determinants in the HPHF function being connected
        !     to all excited HPHF functions.
        ! --> nI will always need to be a unique choice of determinant within
        !     each HPHF function, and then we nevver need to (explicitly)
        !     refering to its spin-coupled partner.
        ! --> If tGenMatEl is true, the Hamiltonian matrix element between the
        !     two determinants will be calculated, and returned in Hel.

        integer, intent(in) :: nI(nel)
        integer(kind=n_int), intent(in) :: iLutnI(0:niftot)
        integer, intent(in) :: exFlag
        integer, intent(out) :: nJ(nel)
        integer(kind=n_int), intent(out) :: iLutnJ(0:niftot)
        integer, intent(out) :: IC, ExcitMat(2,2)
        logical, intent(out) :: tParity ! Not used
        real(dp), intent(out) :: pGen
        HElement_t(dp), intent(out) :: HEl
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = "gen_hphf_excit"

        integer(kind=n_int) :: iLutnJ2(0:niftot)
        integer :: openOrbsI, openOrbsJ, nJ2(nel), ex2(2,2), excitLevel
        real(dp) :: pGen2
        HElement_t(dp) :: MatEl, MatEl2
        logical :: tSign, tSignOrig
        logical :: tSwapped

        ! Avoid warnings
        tParity = .false.

        ! [W.D] this whole hphf should be optimized.. and cleaned up
        ! because it is a mess really..
        ! Generate a normal excitation.

        if (t_back_spawn .or. t_back_spawn_flex) then
            if (tUEGNewGenerator .and. tLatticeGens) then
                call gen_excit_back_spawn_ueg_new(nI, ilutnI, nJ, ilutnJ, exFlag, ic, &
                                          ExcitMat, tSignOrig, pgen, Hel, store, part_type)
            else if (tUEG .and. tLatticeGens) then
                call gen_excit_back_spawn_ueg(nI, ilutnI, nJ, ilutnJ, exFlag, ic, &
                                          ExcitMat, tSignOrig, pgen, Hel, store, part_type)
            else if (tHUB .and. tLatticeGens) then
                call gen_excit_back_spawn_hubbard (nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat,&
                                 tSignOrig, pGen, HEl, store, part_type)
            else
                call gen_excit_back_spawn(nI, ilutnI, nJ, ilutnJ, exFlag, ic, &
                                          ExcitMat, tSignOrig, pgen, Hel, store, part_type)
            end if

        else if (tGen_4ind_weighted) then
            call gen_excit_4ind_weighted (nI, ilutnI, nJ, ilutnJ, exFlag, ic, &
                                          ExcitMat, tSignOrig, pGen, Hel,&
                                          store)
        else if (tGen_4ind_reverse) then
            call gen_excit_4ind_reverse (nI, ilutnI, nJ, ilutnJ, exFlag, ic, &
                                          ExcitMat, tSignOrig, pGen, Hel,&
                                          store)
        else if (tGen_4ind_2) then
            call gen_excit_4ind_weighted2(nI, ilutnI, nJ, ilutnJ, exFlag, ic, &
                                          ExcitMat, tSignOrig, pGen, Hel, &
                                          store)
        else if (t_pchb_excitgen) then
            call gen_rand_excit_pchb(nI, ilutnI, nJ, iLutnJ, exFlag, IC, ExcitMat,&
                                 tSignOrig, pGen, HEl, store)
        else
            call gen_rand_excit (nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat,&
                                 tSignOrig, pGen, HEl, store)
        end if

!        Count=Count+1
!        WRITE(6,*) "COUNT: ",Count
!        CALL neci_flush(6)

        ! Create excitation of uniquely chosen determinant in this HPHF
        ! function.
        IF(IsNullDet(nJ)) RETURN

        ! Create bit representation of excitation - iLutnJ.
        ! n.b. 4ind_weighted does this already.
        if (.not. (tGen_4ind_weighted .or. tGen_4ind_reverse .or. tGen_4ind_2)) &
            CALL FindExcitBitDet(iLutnI,iLutnJ,IC,ExcitMat)

!Test!
!        CALL CalcNonUniPGen(ExcitMat,IC,ClassCount2,ClassCountUnocc2,pDoub,pGen2)
!        IF(abs(pGen-pGen2).gt.1.0e-7_dp_dp) THEN
!            WRITE(6,*) "*******, PGens Incorrect"
!            CALL Stop_All("ouvbou","OUBOU")
!        ENDIF

        IF(TestClosedShellDet(iLutnJ)) THEN
!There is only one way which we could have generated the excitation nJ since it has
!no spin-partner. Also, we will always return the 'correct' version.
            IF(tGenMatHEl) THEN
!Generate matrix element -> HPHF to closed shell det.
                IF(TestClosedShellDet(iLutnI)) THEN
                    !Closed shell -> Closed Shell
                    if(tOddS_HPHF) then
                        call stop_all("gen_hphf_excit","Should not be at closed shell det with Odd S")
                    else
                        HEl = sltcnd_excit (nI, IC, ExcitMat, tSignOrig)
                    endif
                ELSE
                    !Open shell -> Closed Shell
                    if(tOddS_HPHF) then
                        !Odd S States cannot have CS components
                        HEl=0.0_dp
                    else
                        MatEl = sltcnd_excit (nI, IC, ExcitMat, tSignOrig)
                        HEl=MatEl*SQRT(2.0_dp)
                    endif
                ENDIF
                if (IC /= 0 .and. modk_offdiag) hel = -abs(hel)
            ENDIF
        ELSE
!Open shell excitation - could we have generated the spin-coupled determinant instead?

!Find the open shell version.
            CALL ReturnAlphaOpenDet(nJ,nJ2,iLutnJ,iLutnJ2,.true.,.true.,tSwapped)

!            CALL FindExcitBitDetSym(iLutnJ,iLutnJ2)
!Try and find if spin-coupled determinant from excitation is attached.
            IF(tSwapped) THEN
                ExcitLevel = FindBitExcitLevel(iLutnI, iLutnJ, 2)
            ELSE
                ExcitLevel = FindBitExcitLevel(iLutnI, iLutnJ2, 2)
            ENDIF

            IF((ExcitLevel.eq.2).or.(ExcitLevel.eq.1)) THEN     !This is if we have all determinants in the two HPHFs connected...

                Ex2(1,1)=ExcitLevel
!                CALL DecodeBitDet(nJ2,iLutnJ2)     !This could be done better !***!

                IF(tSwapped) THEN
!                    CALL GetExcitation(nI,nJ,NEl,Ex2,tSign) !This could be done more efficiently... !***!
                    CALL GetBitExcitation(iLutnI,iLutnJ,Ex2,tSign)
                ELSE
!                    CALL GetExcitation(nI,nJ2,NEl,Ex2,tSign)
                    CALL GetBitExcitation(iLutnI,iLutnJ2,Ex2,tSign)
                ENDIF
                CALL CalcNonUniPGen(nI, ilutnI, Ex2, ExcitLevel, &
                                    store%ClassCountOcc, &
                                    store%ClassCountUnocc, pDoubles, pGen2, part_type)
!!We cannot guarentee that the pGens are going to be the same - in fact, generally, they wont be.
                pGen=pGen+pGen2

                IF(tGenMatHEl) THEN
!Generate matrix element to open shell excitation
                    IF(TestClosedShellDet(iLutnI)) THEN    !Closed shell -> Open shell : Want to sum in SQRT(2)* Hij
                        if(tOddS_HPHF) then
                            !Cannot have CS components
                            HEl=0.0_dp
                        else
                            IF(tSwapped) THEN
                                MatEl = sltcnd_excit (nI, IC, Ex2, tSign)
                            ELSE
                                MatEl = sltcnd_excit (nI, IC, ExcitMat, &
                                                      tSignOrig)
                            ENDIF
                            HEl=MatEl*SQRT(2.0_dp)
                        endif
                    ELSE     !Open shell -> Open shell

!First find nI -> nJ. If nJ has swapped, then this will be different.
                        IF(tSwapped) THEN
                            MatEl = sltcnd_excit (nI, ExcitLevel, Ex2, &
                                                  tSign)
                        ELSE
                            MatEl = sltcnd_excit (nI, IC, ExcitMat, &
                                                  tSignOrig)
                        ENDIF

                        !now nI2 -> nJ (modelled as nI -> nJ2 with appropriate sign modifications)

!                        CALL FindDetSpinSym(nI,nI2,NEl)
!                        CALL FindExcitBitDetSym(iLutnI,iLutnI2)
!                        ExcitLevel2 = FindBitExcitLevel(iLutnI2, iLutnJ, 2)

                        IF((ExcitLevel.eq.2).or.(ExcitLevel.eq.1)) THEN

                            CALL CalcOpenOrbs(iLutnJ,OpenOrbsJ)
                            CALL CalcOpenOrbs(iLutnI,OpenOrbsI)

                            IF(tSwapped) THEN
                                IF((OpenOrbsJ+OpenOrbsI).eq.3) tSignOrig=.not.tSignOrig
 !I.e. J odd and I even or vice versa, but since these can only be at max quads, then they can only have 1/2 open orbs

                                MatEl2 = sltcnd_excit (nI, IC, ExcitMat, &
                                                       tSignOrig)
                            ELSE
!I.e. J odd and I even or vice versa, but since these can only be at max quads, then they can only have 1/2 open orbs
                                IF((OpenOrbsJ+OpenOrbsI).eq.3) tSign=.not.tSign

                                MatEl2 = sltcnd_excit (nI,  ExcitLevel, &
                                                       Ex2, tSign)
                            ENDIF

                            IF(tOddS_HPHF) THEN
!again, since these can only be at max quads, then they can only have 1/2 open orbs...
                                IF(OpenOrbsI.eq.2) THEN
                                    MatEl=MatEl-MatEl2
                                ELSE
                                    MatEl=MatEl+MatEl2
                                ENDIF
                            ELSE
!again, since these can only be at max quads, then they can only have 1/2 open orbs...
                                IF(OpenOrbsI.eq.2) THEN
                                    MatEl=MatEl+MatEl2
                                ELSE
                                    MatEl=MatEl-MatEl2
                                ENDIF
                            ENDIF
!                            WRITE(6,*) "MatEl2 NEW: ",MatEl2
                        ENDIF
                        HEl=MatEl

                    ENDIF   !Endif from open/closed shell det
                    if (IC /= 0 .and. modk_offdiag) hel = -abs(hel)

                ENDIF   !Endif want to generate matrix element

!                CALL ReturnAlphaOpenDet(nJ,nJ2,iLutnJ,iLutnJ2,.false.,.false.,tSwapped)
!Here, we actually know nJ, so don't need to regenerate it...

            ELSEIF(ExcitLevel.eq.0) THEN
!We have generated the same HPHF. MatEl wants to be zero.
                nJ(1)=0
                IF(tGenMatHEl) THEN
                    HEl=0.0_dp
                ENDIF

            ELSE    !Open-shell to Open-shell, but with no cross-connection.

!                CALL ReturnAlphaOpenDet(nJ,nJ2,iLutnJ,iLutnJ2,.false.,.true.,tSwapped)

                IF(tGenMatHEl) THEN
!iLutnI MUST be open-shell here, since otherwise it would have been connected to
!iLutnJ2. Also, we know the cross connection (i.e. MatEl2 = 0)
                    IF(tSwapped) THEN
                        CALL CalcOpenOrbs(iLutnJ,OpenOrbsJ)
!                        CALL CalcOpenOrbs(iLutnI,OpenOrbsI)
!     IF(((mod(OpenOrbsI,2).eq.1).and.(mod(OpenOrbsJ,2).eq.1)).or.((mod(OpenOrbsI,2).eq.0).and.(mod(OpenOrbsJ,2).eq.1))) THEN
                        IF(tOddS_HPHF) then
                            IF(mod(OpenOrbsJ,2).eq.0) THEN
    !                            WRITE(6,*) "Swapped parity"
                                tSignOrig=.not.tSignOrig
                            ENDIF
                        ELSE
                            IF(mod(OpenOrbsJ,2).eq.1) THEN
    !                            WRITE(6,*) "Swapped parity"
                                tSignOrig=.not.tSignOrig
                            ENDIF
                        ENDIF
                        MatEl = sltcnd_excit(nI,  IC, ExcitMat, tSignOrig)
                    ELSE
                        MatEl = sltcnd_excit (nI, IC, ExcitMat, tSignOrig)
                    ENDIF

                    HEl=MatEl
                    if (IC /= 0 .and. modk_offdiag) hel = -abs(hel)

                ENDIF


            ENDIF

        ENDIF

        ! [W.D.]
        ! i should also abort here already if the matrix element
        ! if below a threshold to optimize the calculation
!         if (abs(Hel) < EPS) then
!             nJ(1) = 0
!             pgen = 0.0_dp
!             Hel = 0.0_dp
!             return
!         end if
!
!         if (t_matele_cutoff) then
!             if (abs(Hel) < matele_cutoff) then
!                 Hel = 0.0_dp
!                 nJ(1) = 0
!                 pgen = 0.0_dp
!                 return
!             end if
!         end if

!        CALL HPHFGetOffDiagHElement(nI,nJ,iLutnI,iLutnJ,MatEl2)
!        IF((MatEl2-MatEl).gt.1.0e-7_dp_dp) THEN
!            WRITE(6,*) MatEl2,MatEl
!            CALL Stop_All("ikb","Error in getting correct HEl - 2")
!        ENDIF
!        if(TestClosedShellDet(iLutnJ)) then
!            write(6,*) "Trying to excite TO closed shell det - why!?"
!            write(6,*) nI
!            write(6,*) "***"
!            write(6,*) iLutnI
!            write(6,*) "HEl: ",HEl
!            if(HEl.ne.0.0_dp) call stop_all("gen_hphf_excit","WHY?!")
!            call neci_flush(6)
!        endif

    end subroutine

!This routine will take a determinant, and create the determinant whose final open-shell
!spatial orbital contains an alpha electron.
!If the final open-shell electron is a beta orbital, then the balue of the bit-string
!will be smaller. We are interested in returning
!the larger of the open-shell bit strings since this will correspond to the final
!open-shell electron being an alpha.
!This rationalization may well break down when it comes to the negative bit (32),
!however, this may not matter, since all we really
!need is a unique description of a HPHF...?
!iLutnI (nI) is returned as this determinant, with iLutSym (nJ) being the other.
!If tCalciLutSym is false, iLutSym will be calculated from iLutnI. Otherwise, it won't.
    SUBROUTINE ReturnAlphaOpenDet(nI,nJ,iLutnI,iLutSym,tCalciLutSym,tCalcnISym,tSwapped)
        INTEGER(KIND=n_int), intent(inout) :: iLutSym(0:NIfTot),iLutnI(0:NIfTot)
        integer(kind=n_int) :: iLutTemp(0:NIfTot)
        INTEGER :: i,nTemp(NEl)
        integer, intent(inout) :: nJ(NEl),nI(NEl)
        LOGICAL, intent(in) :: tCalciLutSym,tCalcnISym
        logical, intent(out) :: tSwapped

        IF(tCalciLutSym) THEN
            CALL FindExcitBitDetSym(iLutnI,iLutSym)
        ENDIF
        IF(tCalcnISym) THEN
            CALL FindDetSpinSym(nI,nJ,NEl)
        ENDIF

        ! iLutnI is 'less' than iLutSym, so iLutSym is the determinant with
        ! the first open-shell = alpha. Swap them around.
        ! Only count up to NIfD to avoid Yamanouchi symbol etc.
        i=DetBitLT(iLutnI, iLutSym, NIfD)
        IF(i.eq.1) THEN
            iLutTemp(:)=iLutnI(:)
            iLutnI(:)=iLutSym(:)
            iLutSym(:)=iLutTemp(:)
!            CALL FindDetSpinSym(nI,nJ,NEl)
            nTemp(:)=nI(:)
            nI(:)=nJ(:)
            nJ(:)=nTemp(:)
            tSwapped=.true.
        ELSEIF(i.eq.0) THEN
            CALL Stop_All("ReturnAlphaOpenDet","Shouldn't have closed shell determinants in here")
        ELSE
            tSwapped=.false.
        ENDIF

    END SUBROUTINE ReturnAlphaOpenDet


!This create the spin-coupled determinant of nI in nJ in natural ordered form.
    PURE SUBROUTINE FindDetSpinSym(nI,nJ,NEl)
        INTEGER, intent(in) :: NEl,nI(NEl)
        integer, intent(out) :: nJ(NEl)
        integer :: i

        do i = 1, nel

            ! If electron is an alpha electron, change it to a beta (unless
            ! it is part of a closed pair of electrons).
            if (is_alpha(nI(i))) then
                if (i == 1) then
                    nJ(i) = nI(i) - 1
                elseif(get_beta(nI(i)) /= nI(i-1)) then
                    nJ(i) = nI(i) - 1
                else
                    nJ(i) = nI(i)
                endif
            ! vice-versa for beta.
            else
                if (i == nel) then
                    nJ(i) = nI(i) + 1
                elseif(get_alpha(nI(i)) /= nI(i+1)) then
                    nJ(i) = nI(i) + 1
                else
                    nJ(i) = nI(i)
                endif
            endif
        enddo

!        nTemp(:)=nJ(:)
!        CALL NECI_SORTI(NEl,nTemp)
!        do i=1,NEl
!            IF(nTemp(i).ne.nJ(i)) THEN
!                STOP 'Massive Error'
!            ENDIF
!        enddo

    END SUBROUTINE FindDetSpinSym

!In closed-shell systems with equal number of alpha and beta strings, the amplitude of a
!determinant in the final CI wavefunction is the same
!when the alpha and beta electrons are swapped (for S=0, see Helgakker for more details).
!It will sometimes be necessary to find this other
!determinant when spawning. This routine will find the bit-representation of an excitation
!by constructing the symmetric iLut from the its
!symmetric partner, also in bit form.
    PURE SUBROUTINE FindExcitBitDetSym(iLut,iLutSym)
        IMPLICIT NONE
        INTEGER(KIND=n_int) , intent(in) :: iLut(0:NIfTot)
        INTEGER(KIND=n_int) , intent(out) :: iLutSym(0:NIfTot)
        INTEGER(KIND=n_int) :: iLutAlpha(0:NIfTot),iLutBeta(0:NIfTot)
        INTEGER :: i

        iLutSym(:)=0
        iLutAlpha(:)=0
        iLutBeta(:)=0

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


!!This routine will take a HPHF nI, and find Iterations number of excitations of it.
!It will then histogram these, summing in 1/pGen for every occurance of
!!the excitation. This means that all excitations should be 0 or 1 after enough iterations.
!It will then count the excitations and compare the number to the
!!number of excitations generated using the full enumeration excitation generation.
!    SUBROUTINE TestGenRandHPHFExcit(nI,Iterations,pDoub)
!        Use SystemData , only : NEl,nBasis,G1,nBasisMax
!        use DetBitOps, only: EncodeBitDet
!        use bit_reps, only: decode_bit_det
!        use GenRandSymExcitNuMod, only: scratchsize
!        use FciMCData, only: tGenMatHEl
!        use util_mod, only: get_free_unit
!        IMPLICIT NONE
!        INTEGER :: ClassCount2(ScratchSize),nIX(NEl)
!        INTEGER :: ClassCountUnocc2(ScratchSize)
!        INTEGER :: i,Iterations,nI(NEl),nJ(NEl),DetConn,nI2(NEl),nJ2(NEl),DetConn2,iUniqueHPHF,iUniqueBeta,PartInd,ierr,iExcit
!        real(dp) :: pDoub,pGen
!        LOGICAL :: Unique,TestClosedShellDet,Die,tGenClassCountnI,tSwapped
!        INTEGER(KIND=n_int) :: iLutnI(0:NIfTot),iLutnJ(0:NIfTot),iLutnI2(0:NIfTot),iLutSym(0:NIfTot)
!        INTEGER(KIND=n_int), ALLOCATABLE :: ConnsAlpha(:,:),ConnsBeta(:,:),UniqueHPHFList(:,:)
!        INTEGER , ALLOCATABLE :: ExcitGen(:)
!        real(dp) , ALLOCATABLE :: Weights(:)
!        INTEGER :: iMaxExcit,nStore(6),nExcitMemLen,j,k,l, iunit
!        integer :: icunused, exunused_var(2,2), scratch3(scratchsize)
!        integer :: icunused, exunused_var(2,2), scratch3(scratchsize)
!        logical :: tParityunused, tTmp
!
!        CALL EncodeBitDet(nI,iLutnI)
!        CALL FindDetSpinSym(nI,nI2,NEl)
!        CALL EncodeBitDet(nI2,iLutnI2)
!        IF(TestClosedShellDet(iLutnI)) THEN
!            IF(.not.DetBitEQ(iLutnI,iLutnI2,NIfDBO)) THEN
!                CALL Stop_All("TestGenRandHPHFExcit","Closed shell determinant entered, but alpha and betas different...")
!            ENDIF
!        ENDIF
!        WRITE(6,*) "nI: ",nI(:)
!        WRITE(6,*) ""
!        WRITE(6,*) "nISym: ",nI2(:)
!        WRITE(6,*) ""
!        WRITE(6,*) "iLutnI: ",iLutnI(:)
!        WRITE(6,*) "iLutnISymL ",iLutnI2(:)
!        WRITE(6,*) "***"
!        WRITE(6,*) Iterations,pDoub
!!        WRITE(6,*) "nSymLabels: ",nSymLabels
!        CALL neci_flush(6)
!
!!First, we need to enumerate all possible HPHF wavefunctions from each spin-pair of determinants.
!!These need to be stored in an array
!!Find the number of symmetry allowed excitations there should be by looking at the full excitation generator.
!!Setup excit generators for this determinant
!        iMaxExcit=0
!        nStore(1:6)=0
!        CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
!        ALLOCATE(EXCITGEN(nExcitMemLen),stat=ierr)
!        IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
!        EXCITGEN(:)=0
!        CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.TRUE.,EXCITGEN,nJ,iMaxExcit,0,nStore,3)
!        CALL GetSymExcitCount(EXCITGEN,DetConn)
!        WRITE(6,*) "Alpha determinant has ",DetConn," total excitations:"
!        ALLOCATE(ConnsAlpha(0:NIfTot,DetConn))
!        i=1
!        lp2: do while(.true.)
!            CALL GenSymExcitIt2(nI,NEl,G1,nBasis,nBasisMax,.false.,EXCITGEN,nJ,iExcit,0,nStore,3)
!            IF(IsNullDet(nJ)) exit lp2
!            CALL EncodeBitDet(nJ,iLutnJ)
!            IF(.not.TestClosedShellDet(iLutnJ)) THEN
!                CALL ReturnAlphaOpenDet(nJ,nJ2,iLutnJ,iLutSym,.true.,.true.,tSwapped)
!            ENDIF
!!            WRITE(6,"(4I4,A,I4,A,I13)") nJ(:), " *** ", iExcit, " *** ", iLutnJ(:)
!            ConnsAlpha(0:NIfD,i)=iLutnJ(:)
!            i=i+1
!        enddo lp2
!
!!Now we also need to store the excitations from the other spin-coupled determinant.
!        iMaxExcit=0
!        nStore(1:6)=0
!        DEALLOCATE(EXCITGEN)
!
!        CALL GenSymExcitIt2(nI2,NEl,G1,nBasis,nBasisMax,.TRUE.,nExcitMemLen,nJ,iMaxExcit,0,nStore,3)
!        ALLOCATE(EXCITGEN(nExcitMemLen),stat=ierr)
!        IF(ierr.ne.0) CALL Stop_All("SetupExcitGen","Problem allocating excitation generator")
!        EXCITGEN(:)=0
!        CALL GenSymExcitIt2(nI2,NEl,G1,nBasis,nBasisMax,.TRUE.,EXCITGEN,nJ,iMaxExcit,0,nStore,3)
!        CALL GetSymExcitCount(EXCITGEN,DetConn2)
!        WRITE(6,*) "Beta determinant has ",DetConn2," total excitations"
!        ALLOCATE(ConnsBeta(0:NIfTot,DetConn2))
!        i=1
!        lp: do while(.true.)
!            CALL GenSymExcitIt2(nI2,NEl,G1,nBasis,nBasisMax,.false.,EXCITGEN,nJ,iExcit,0,nStore,3)
!            IF(IsNullDet(nJ)) exit lp
!            CALL EncodeBitDet(nJ,iLutnJ)
!            IF(.not.TestClosedShellDet(iLutnJ)) THEN
!                CALL ReturnAlphaOpenDet(nJ,nJ2,iLutnJ,iLutSym,.true.,.true.,tSwapped)
!            ENDIF
!!            WRITE(6,"(4I4,A,I4,A,I13)") nJ(:), " *** ",iExcit," *** ",iLutnJ(:)
!            ConnsBeta(:,i)=iLutnJ(:)
!            i=i+1
!        enddo lp
!        DEALLOCATE(EXCITGEN)
!
!!Now we need to find how many HPHF functions there are.
!        iUniqueHPHF=0
!        do j=1,DetConn
!!Run though all HPHF in the first array
!            Unique=.true.
!            do k=j-1,1,-1
!!Run backwards through the array to see if this HPHF has come before
!                IF(DetBitEQ(ConnsAlpha(0:NIfTot,k),ConnsAlpha(0:NIfTot,j),NIfDBO)) THEN
!!This HPHF has already been counted before...
!                    Unique=.false.
!                    EXIT
!                ENDIF
!            enddo
!            IF(Unique) THEN
!!Unique HPHF found, count it
!                iUniqueHPHF=iUniqueHPHF+1
!            ENDIF
!        enddo
!
!        iUniqueBeta=0
!
!!Now look through all the excitations for the spin-coupled determinant from the original HPHF...
!        do j=1,DetConn2
!!Run though all excitations in the first array, *and* up to where we are in the second array
!            Unique=.true.
!            do k=1,DetConn
!                IF(DetBitEQ(ConnsAlpha(:,k),ConnsBeta(:,j),NIfDBO)) THEN
!                    Unique=.false.
!                    EXIT
!                ENDIF
!            enddo
!            IF(Unique) THEN
!!Need to search backwards through the entries we've already looked at in this array...
!                do k=j-1,1,-1
!                    IF(DetBitEQ(ConnsBeta(0:NIfTot,k),ConnsBeta(0:NIfTot,j),NIfDBO)) THEN
!                        Unique=.false.
!                        EXIT
!                    ENDIF
!                enddo
!            ENDIF
!            IF(Unique) THEN
!                iUniqueHPHF=iUniqueHPHF+1
!                iUniqueBeta=iUniqueBeta+1
!            ENDIF
!        enddo
!
!        WRITE(6,*) "There are ",iUniqueHPHF," unique HPHF wavefunctions from the HPHF given."
!
!        WRITE(6,*) "There are ",iUniqueBeta," unique HPHF wavefunctions from the spin-coupled
!determinant, which are not in a alpha version."
!        IF(iUniqueBeta.ne.0) THEN
!            WRITE(6,*) "HPHF from beta, but not from alpha!"
!            CALL neci_flush(6)
!            STOP
!        ENDIF
!
!        ALLOCATE(UniqueHPHFList(0:NIfTot,iUniqueHPHF))
!        UniqueHPHFList(:,:)=0
!!Now fill the list of HPHF Excitations.
!        iUniqueHPHF=0
!        do j=1,DetConn
!!Run though all HPHF in the first array
!            Unique=.true.
!            do k=j-1,1,-1
!!Run backwards through the array to see if this HPHF has come before
!                IF(DetBitEQ(ConnsAlpha(0:NIfTot,k),ConnsAlpha(0:NIfTot,j),NIfDBO)) THEN
!!This HPHF has already been counted before...
!                    Unique=.false.
!                    EXIT
!                ENDIF
!            enddo
!            IF(Unique) THEN
!!Unique HPHF found, count it
!                iUniqueHPHF=iUniqueHPHF+1
!                UniqueHPHFList(:,iUniqueHPHF)=ConnsAlpha(:,j)
!            ENDIF
!        enddo
!
!!Now look through all the excitations for the spin-coupled determinant from the original HPHF...
!        do j=1,DetConn2
!!Run though all excitations in the first array, *and* up to where we are in the second array
!            Unique=.true.
!            do k=1,DetConn
!                IF(DetBitEQ(ConnsAlpha(:,k),ConnsBeta(:,j),NIfDBO)) THEN
!                    Unique=.false.
!                    EXIT
!                ENDIF
!            enddo
!            IF(Unique) THEN
!!Need to search backwards through the entries we've already looked at in this array...
!                do k=j-1,1,-1
!                    IF(DetBitEQ(ConnsBeta(:,k),ConnsBeta(:,j),NIfDBO)) THEN
!                        Unique=.false.
!                        EXIT
!                    ENDIF
!                enddo
!            ENDIF
!            IF(Unique) THEN
!                iUniqueHPHF=iUniqueHPHF+1
!                UniqueHPHFList(:,iUniqueHPHF)=ConnsBeta(0:NIfD,j)
!            ENDIF
!        enddo
!
!!Now sort the list, so that it can be easily binary searched.
!        ALLOCATE(ExcitGen(iUniqueHPHF))
!        call sort (UniqueHPHFList(:,1:iUniqueHPHF), ExcitGen)
!        DEALLOCATE(ExcitGen)
!
!        WRITE(6,*) "Unique HPHF wavefunctions are: "
!        do i=1,iUniqueHPHF
!            WRITE(6,*) UniqueHPHFList(0:NIfTot,i)
!        enddo
!
!        ALLOCATE(Weights(iUniqueHPHF))
!        Weights(:)=0.0_dp
!        tGenClassCountnI=.false.
!
!        do i=1,Iterations
!
!            IF(mod(i,10000).eq.0) WRITE(6,"(A,I10)") "Iteration: ",i
!
!            CALL GenRandHPHFExcit(nI,iLutnI,nJ,iLutnJ,pDoub,3,pGen)
!            tTmp = tGenMatHEl
!            tGenMatHel = .false.
!            call gen_hphf_excit (nI, iLutnI, nJ, iLutnJ, 3, icunused, &
!                                 exunused, tparityunused, pGen, &
!                                 tGenClassCountnI, ClassCount2, &
!                                 ClassCountUnocc2, scratch3)
!            tGenMatHel = tTmp
!!            CALL GenRandSymExcitNU(nI,iLut,nJ,pDoub,IC,ExcitMat,TParity,exFlag,pGen)
!
!!Search through the list of HPHF wavefunctions to find slot.
!            CALL BinSearchListHPHF(iLutnJ,UniqueHPHFList(:,1:iUniqueHPHF),iUniqueHPHF,1,iUniqueHPHF,PartInd,Unique)
!
!            IF(.not.Unique) THEN
!                CALL Stop_All("TestGenRandHPHFExcit","Cannot find excitation in list of allowed excitations")
!            ENDIF
!
!            Weights(PartInd)=Weights(PartInd)+(1.0_dp/pGen)
!
!!Check excitation
!!            CALL IsSymAllowedExcit(nI,nJ,IC,ExcitMat)
!
!        enddo
!
!        iunit = get_free_unit()
!        OPEN(iunit,FILE="PGenHist",STATUS="UNKNOWN")
!
!!normalise excitation probabilities
!        Die=.false.
!        do i=1,iUniqueHPHF
!            Weights(i)=Weights(i)/real(Iterations,8)
!            IF(abs(Weights(i)-1.0_dp).gt.0.1) THEN
!                WRITE(6,*) "Error here!"
!                Die=.true.
!            ENDIF
!            WRITE(6,*) i,UniqueHPHFList(0:NIfTot,i),Weights(i)
!            call decode_bit_det (nIX, UniqueHPHFList(0:NIfTot,i))
!            WRITE(6,*) nIX(:)
!            WRITE(iunit,"(I8,G25.10)",advance='no') i,Weights(i)
!            do j=0,NIfTot-1
!                WRITE(iunit,"(I16)",advance='no') UniqueHPHFList(j,i)
!            enddo
!            WRITE(iunit,"(I16)") UniqueHPHFList(NIfTot,i)
!        enddo
!
!        CLOSE(iunit)
!        IF(Die) THEN
!            CALL Stop_All("IUB","TestFail")
!        ENDIF
!
!    END SUBROUTINE TestGenRandHPHFExcit

    SUBROUTINE BinSearchListHPHF(iLut,List,Length,MinInd,MaxInd,PartInd,tSuccess)
        INTEGER :: Length,MinInd,MaxInd,PartInd
        INTEGER(KIND=n_int) :: iLut(0:NIfTot),List(0:NIfTot,Length)
        INTEGER :: i,j,N,Comp
        LOGICAL :: tSuccess

!        WRITE(6,*) "Binary searching between ",MinInd, " and ",MaxInd
!        CALL neci_flush(6)
        i=MinInd
        j=MaxInd
        IF(i-j.eq.0) THEN
            Comp=DetBitLT(List(:,MaxInd),iLut(:),NIfDBO)
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
            Comp=DetBitLT(List(:,N),iLut(:),NIfDBO)

            IF(Comp.eq.0) THEN
!Praise the lord, we've found it!
                tSuccess=.true.
                PartInd=N
                RETURN
            ELSEIF((Comp.eq.1).and.(i.ne.N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for.
!Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i=N
            ELSEIF(i.eq.N) THEN


                IF(i.eq.MaxInd-1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
                    Comp=DetBitLT(List(:,i+1),iLut(:),NIfDBO)
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



    subroutine CalcNonUniPGen(nI, ilutI, ex, ic, ClassCount2, &
                              ClassCountUnocc2, pDoub, pGen, part_type)

        ! This routine will calculate the PGen between two connected
        ! determinants, nI and nJ which are IC excitations of each other, using
        ! the unbiased scheme.
        !
        ! Only the excitation matrix is needed (1,*) are the i,j orbs, and
        ! (2,*) are the a,b orbs. This is the prob of generating nJ FROM nI,
        ! not the other way round. Passed in is also the ClassCount2 arrays for
        ! nI, and the probability of picking a double.
        !
        ! A word of warning: The routine does not check that the determinants
        ! are indeed connected, and may well return a non-zero probability even
        ! if they arent. Therefore, make sure that they are at most double
        ! excitations of each other.
        !
        ! nI is the determinant from which the excitation comes from.

        use bit_reps, only: get_initiator_flag
        use bit_rep_data, only: test_flag

        integer, intent(in) :: nI(nel), ex(2,2), ic
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ClassCount2(ScratchSize)
        integer, intent(in) :: ClassCountUnocc2(ScratchSize)
        real(dp), intent(in) :: pDoub
        real(dp), intent(out) :: pGen
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'CalcNonUniPGen'

        integer :: temp_part_type

        ! We need to consider which of the excitation generators are in use,
        ! and call the correct routine in each case.
        ASSERT(.not. (tCSF)) ! .or. tSpinProjDets

        pgen = 0.0_dp

        ! i have to make sure to catch all this call to this function correctly
        if (present(part_type)) then
            temp_part_type = part_type
        else
            temp_part_type = 1
        end if

        ! does it help to avoid recalculating for the reference?
        ! do i need to  check if it is actually a non-initiator?
        ! i guess i do.. or i go the unnecessary way of checking again in
        ! the called back-spawn functions
        if ((t_back_spawn .or. t_back_spawn_flex) .and. &
            (.not. DetBitEq(ilutI,ilutRef(:,temp_part_type),nifdbo)) .and. &
            (.not. test_flag(ilutI, get_initiator_flag(temp_part_type)))) then
            ! i just realised this also has to be done for the hubbard
            ! and the ueg model.. -> create those functions!
            if (tHUB .and. tLatticeGens) then
                pgen = calc_pgen_back_spawn_hubbard(nI, ilutI, ex, ic, temp_part_type)
            else if (tUEGNewGenerator .and. tLatticeGens) then
                pgen = calc_pgen_back_spawn_ueg_new(nI, ilutI, ex, ic, temp_part_type)
            else if (tUEG .and. tLatticeGens) then
                pgen = calc_pgen_back_spawn_ueg(ilutI, ex, ic, temp_part_type)
            else
                pgen = calc_pgen_back_spawn(nI, ilutI, ex, ic, temp_part_type)
            end if

        ! this if construct is not well setup.. this can fail..
        else
            if (tLatticeGens) then
                if (ic == 2) then
                    call CalcPGenLattice (ex, pGen)
                else
                    pGen = 0
                end if
            else if (tGen_4ind_2) then
                pgen = calc_pgen_4ind_weighted2(nI, ilutI, ex, ic)

            else if (tGen_4ind_weighted) then
                pgen = calc_pgen_4ind_weighted (nI, ilutI, ex, ic, &
                                                ClassCountUnocc2)
            else if (tGen_4ind_reverse) then
                pgen = calc_pgen_4ind_reverse (nI, ilutI, ex, ic)
            else if (t_pchb_excitgen) then
                pgen = calc_pgen_pchb(nI, ex, ic, ClassCount2, ClassCountUnocc2)
            else
                ! Here we assume that the normal excitation generators in
                ! symrandexcit2.F90 are being used.
                call calc_pgen_symrandexcit2 (nI, ex, ic, ClassCount2, &
                                              ClassCountUnocc2, pDoub, pGen)
            end if
        end if

    end subroutine



END MODULE HPHFRandExcitMod


