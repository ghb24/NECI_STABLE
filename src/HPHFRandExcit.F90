#include "macros.h"

MODULE HPHFRandExcitMod
!Half-projected HF wavefunctions are a linear combination of two HF determinants,
!where all alphas -> betas and betas -> alpha to create the pair.
!In closed-shell systems, these two determinants have the same FCI amplitude, and so it is easier to treat them as a pair.
!The probability of creating this HPHF where both are from pairs of spin-coupled determinants (i & j -> a & b):
![ P(i->a) + P(i->b) + P(j->a) + P(j->b) ]/2
!We therefore need to find the excitation matrix between the determinant which wasn't
!excited and the determinant which was created.

    use SystemData, only: nel, Alat, G1, nbasis, nbasismax, nmsh, arr, &
                          tOddS_HPHF, modk_offdiag, tGen_4ind_weighted, &
                          tGen_4ind_reverse, tLatticeGens, tGen_4ind_2, tHUB, &
                          tUEG, tUEGNewGenerator, t_new_real_space_hubbard, &
                          t_tJ_model, t_heisenberg_model, t_lattice_model, &
                          t_k_space_hubbard, t_3_body_excits, t_uniform_excits, &
                          t_trans_corr_hop, t_spin_dependent_transcorr, &
                          t_pchb_excitgen, t_mol_3_body, t_ueg_3_body, tGUGA, &
                          t_pcpp_excitgen, max_ex_level, t_guga_pchb

    use IntegralsData, only: UMat, fck, nMax

    use SymData, only: nSymLabels

    use dSFMT_interface, only: genrand_real2_dSFMT

    use GenRandSymExcitNUMod, only: gen_rand_excit, calc_pgen_symrandexcit2, &
                                    ScratchSize, CalcPGenLattice, construct_class_counts
    use excit_gens_int_weighted, only: gen_excit_4ind_weighted, &
                                       gen_excit_4ind_reverse, &
                                       calc_pgen_4ind_weighted, &
                                       calc_pgen_4ind_reverse
    use tc_three_body_excitgen, only: gen_excit_mol_tc, calc_pgen_triple
    use DetBitOps, only: DetBitLT, DetBitEQ, FindExcitBitDet, &
                         FindBitExcitLevel, MaskAlpha, MaskBeta, &
                         TestClosedShellDet, CalcOpenOrbs, IsAllowedHPHF, &
                         DetBitEQ, GetBitExcitation

    use FciMCData, only: pDoubles, ilutRef

    use constants, only: dp, n_int, EPS, maxExcit

    use sltcnd_mod, only: dyn_sltcnd_excit_old

    use bit_reps, only: NIfD, NIfTot

    use SymExcitDataMod, only: excit_gen_store_type

    use excit_gen_5, only: calc_pgen_4ind_weighted2, gen_excit_4ind_weighted2

    use pcpp_excitgen, only: calc_pgen_pcpp, gen_rand_excit_pcpp, create_elec_map

    use sort_mod

    use HElem

    use CalcData, only: t_matele_cutoff, matele_cutoff, t_back_spawn, t_back_spawn_flex

    use back_spawn_excit_gen, only: gen_excit_back_spawn, calc_pgen_back_spawn, &
                                    gen_excit_back_spawn_ueg, calc_pgen_back_spawn_ueg, &
                                    calc_pgen_back_spawn_hubbard, gen_excit_back_spawn_hubbard, &
                                    gen_excit_back_spawn_ueg_new, calc_pgen_back_spawn_ueg_new

    use real_space_hubbard, only: gen_excit_rs_hubbard, calc_pgen_rs_hubbard, &
                                  gen_excit_rs_hubbard_transcorr, &
                                  calc_pgen_rs_hubbard_transcorr, &
                                  gen_excit_rs_hubbard_transcorr_uniform, &
                                  calc_pgen_rs_hubbard_transcorr_uniform, &
                                  gen_excit_rs_hubbard_spin_dependent_transcorr, &
                                  calc_pgen_rs_hubbard_spin_dependent_transcorr

    use tJ_model, only: gen_excit_tj_model, gen_excit_heisenberg_model, &
                        calc_pgen_tJ_model, calc_pgen_heisenberg_model

    use lattice_mod, only: get_helement_lattice

    use k_space_hubbard, only: gen_excit_k_space_hub, calc_pgen_k_space_hubbard, &
                               gen_excit_uniform_k_space_hub

    use exc_gen_classes, only: current_exc_generator

    use procedure_pointers, only: generate_excitation_t

    use guga_pchb_excitgen, only: calc_pgen_guga_pchb

    use guga_bitRepOps, only: current_csf_i

    use util_mod, only: stop_all

    IMPLICIT NONE

    procedure(generate_excitation_t), pointer :: exc_generator_for_HPHF => null()

contains

!Calculate probability of exciting from HPHF nI to HPHF nJ
!It is imperative that when using this routine, the 'correct' determinant is sent in
!i.e. the unique determinant representation of the two HPHF functions. This is because
!the classcount arrays will be different for the two determinants.
!tSameFunc will be returned as true if the two HPHF functions are the same
    subroutine CalcPGenHPHF(nI, iLutnI, nJ, iLutnJ, ex, ClassCount, ClassCountUnocc, pDoubles, pGen, tSameFunc)
        integer, intent(in) :: nI(nel)
        integer(kind=n_int), intent(in) :: iLutnI(0:niftot), iLutnJ(0:niftot)
        integer, intent(in) :: ClassCount(ScratchSize), ClassCountUnocc(ScratchSize)
        integer, intent(in) :: nJ(nel), ex(2, maxExcit)
        real(dp), intent(in) :: pDoubles
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tSameFunc
        logical :: tSign, tSwapped
        real(dp) :: pGen2
        integer :: ic
        integer :: Ex2(2, maxExcit), nJ_loc(nel), nJ2(nel)
        integer(kind=n_int) :: iLutnJ_loc(0:niftot), iLutnJ2(0:niftot)
#ifdef DEBUG_
        character(*), parameter :: this_routine = "CalcPGenHPHF"
#endif
        tSameFunc = .false.
        pGen = 0.0_dp

        IF(TestClosedShellDet(iLutnJ)) THEN
            !nJ is CS, therefore, only one way of generating it.
            ic = FindBitExcitLevel(iLutnI, iLutnJ, 2)
            if(ic == 0) then
                tSameFunc = .true.
                return
            end if
            if(ic <= max_ex_level) then
                call CalcNonUniPGen(nI, ilutnI, ex, ic, ClassCount, ClassCountUnocc, pDoubles, pGen)
            end if
        else
            !nJ is openshell. Add the probabilities of generating each pair (if both connected)
            nJ_loc = nJ
            iLutnJ_loc = iLutnJ
            CALL ReturnAlphaOpenDet(nJ_loc, nJ2, iLutnJ_loc, iLutnJ2, .true., .true., tSwapped)

            !First find nI -> nJ
            ic = FindBitExcitLevel(iLutnI, iLutnJ_loc, 2)
            if(ic == 0) then
                tSameFunc = .true.
                return
            end if
            if(ic <= max_ex_level) then
                if(.not. tSwapped) then
                    !ex is correct for this excitation
                    call CalcNonUnipGen(nI, ilutnI, ex, ic, ClassCount, ClassCountUnocc, pDoubles, pGen)
                else
                    Ex2(1, 1) = ic
                    call GetBitExcitation(iLutnI, iLutnJ_loc, Ex2, tSign)
                    call CalcNonUnipGen(nI, ilutnI, Ex2, ic, ClassCount, ClassCountUnocc, pDoubles, pGen)
                end if
            end if

            !Now consider nI -> nJ2 and add the probabilities
            ic = FindBitExcitLevel(iLutnI, iLutnJ2, 2)
            if(ic == 0) then
                tSameFunc = .true.
                return
            end if
            if(ic <= max_ex_level) then
                if(tSwapped) then
                    !ex is correct for this excitation
                    call CalcNonUnipGen(nI, ilutnI, ex, ic, ClassCount, ClassCountUnocc, pDoubles, pGen2)
                else
                    Ex2(1, 1) = ic
                    call GetBitExcitation(iLutnI, iLutnJ2, Ex2, tSign)
                    call CalcNonUnipGen(nI, ilutnI, Ex2, ic, ClassCount, ClassCountUnocc, pDoubles, pGen2)
                end if
                pGen = pGen + pGen2
            end if
        end if

    end subroutine CalcPGenHPHF

    subroutine gen_hphf_excit(nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat, &
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
        integer, intent(out) :: IC, ExcitMat(2, maxExcit)
        logical, intent(out) :: tParity ! Not used
        real(dp), intent(out) :: pGen
        HElement_t(dp), intent(out) :: HEl
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = "gen_hphf_excit"

        integer(kind=n_int) :: iLutnJ2(0:niftot)
        integer :: openOrbsI, openOrbsJ, nJ2(nel), ex2(2, maxExcit), excitLevel
        real(dp) :: pGen2
        HElement_t(dp) :: MatEl, MatEl2
        logical :: tSign, tSignOrig
        logical :: tSwapped
        integer :: temp_ex(2, maxExcit)

#ifdef WARNING_WORKAROUND_
        tParity = .false.
#endif

        call exc_generator_for_HPHF(nI, iLutnI, nJ, iLutnJ, exFlag, IC, ExcitMat, &
                                    tSignOrig, pGen, HEl, store, part_type)
        ! Create excitation of uniquely chosen determinant in this HPHF
        ! function.
        IF(IsNullDet(nJ)) RETURN

        ! Create bit representation of excitation - iLutnJ.
        ! n.b. 4ind_weighted does this already.
        if (.not. any([tGen_4ind_weighted, tGen_4ind_reverse, tGen_4ind_2])) then
            CALL FindExcitBitDet(iLutnI, iLutnJ, IC, ExcitMat)
        end if

        IF(TestClosedShellDet(iLutnJ)) THEN
!There is only one way which we could have generated the excitation nJ since it has
!no spin-partner. Also, we will always return the 'correct' version.
            IF(tGenMatHEl) THEN
!Generate matrix element -> HPHF to closed shell det.
                IF(TestClosedShellDet(iLutnI)) THEN
                    !Closed shell -> Closed Shell
                    if(tOddS_HPHF) then
                        call stop_all("gen_hphf_excit", "Should not be at closed shell det with Odd S")
                    else
                        ! [W.D. 30.10.2017]
                        ! have to change here to use the real-space hubbard
                        ! routines.. thats why this whole HPHF should be
                        ! reworked.
                        ! [W.D. 13.11.2017]
                        ! somehow i reintroduced a bug in the HPHF + hubbard
                        ! implementation with "fixes" in here -> check that!
                        temp_ex(1, :) = ExcitMat(2, :)
                        temp_ex(2, :) = ExcitMat(1, :)
                        if(t_lattice_model) then
                            if(t_k_space_hubbard .or. &
                               (t_new_real_space_hubbard .and. t_trans_corr_hop)) then
                                hel = get_helement_lattice(nJ, ic, temp_ex, tSignOrig)
                            else
                                call Stop_All(this_routine, &
                                              "no closed shell to closed shell possible in real-space lattice models!")
                            end if
                        else
                            HEl = dyn_sltcnd_excit_old(nJ, IC, temp_ex, tSignOrig)
                        end if
                    end if
                ELSE
                    !Open shell -> Closed Shell
                    if(tOddS_HPHF) then
                        !Odd S States cannot have CS components
                        HEl = 0.0_dp
                    else
                        temp_ex(1, :) = ExcitMat(2, :)
                        temp_ex(2, :) = ExcitMat(1, :)
                        if(t_lattice_model) then
                            Matel = get_helement_lattice(nJ, ic, temp_ex, tSignOrig)
                        else
                            MatEl = dyn_sltcnd_excit_old(nJ, IC, temp_ex, tSignOrig)
                        end if
                        HEl = MatEl * SQRT(2.0_dp)
                    end if
                end if
                if(IC /= 0 .and. modk_offdiag) hel = -abs(hel)
            end if
        ELSE
!Open shell excitation - could we have generated the spin-coupled determinant instead?

!Find the open shell version.
            CALL ReturnAlphaOpenDet(nJ, nJ2, iLutnJ, iLutnJ2, .true., .true., tSwapped)

!Try and find if spin-coupled determinant from excitation is attached.
            IF(tSwapped) THEN
                ExcitLevel = FindBitExcitLevel(iLutnI, iLutnJ, 2)
            ELSE
                ExcitLevel = FindBitExcitLevel(iLutnI, iLutnJ2, 2)
            end if

            IF(ExcitLevel <= max_ex_level .and. ExcitLevel > 0) THEN     !This is if we have all determinants in the two HPHFs connected...

                Ex2(1, 1) = ExcitLevel

                IF(tSwapped) THEN
                    CALL GetBitExcitation(iLutnI, iLutnJ, Ex2, tSign)
                ELSE
                    CALL GetBitExcitation(iLutnI, iLutnJ2, Ex2, tSign)
                end if
                ! As we are passing store%ClassCountOcc/Unocc, we have to make sure they are
                ! set
                if(.not. store%tFilled) then
                    CALL construct_class_counts(nI, store%ClassCountOcc, &
                        store%ClassCountUnocc)
                    if(t_pcpp_excitgen) store%elec_map = create_elec_map(ilutnI)
                    store%tFilled = .true.
                end if

                CALL CalcNonUniPGen(nI, ilutnI, Ex2, ExcitLevel, &
                                    store%ClassCountOcc, &
                                    store%ClassCountUnocc, pDoubles, pGen2, part_type)

!!We cannot guarentee that the pGens are going to be the same - in fact, generally, they wont be.
                pGen = pGen + pGen2

                IF(tGenMatHEl) THEN
!Generate matrix element to open shell excitation
                    IF(TestClosedShellDet(iLutnI)) THEN    !Closed shell -> Open shell : Want to sum in SQRT(2)* Hij
                        if(tOddS_HPHF) then
                            !Cannot have CS components
                            HEl = 0.0_dp
                        else
                            if(tSwapped) then
                                temp_ex(1, :) = ex2(2, :)
                                temp_ex(2, :) = ex2(1, :)
                                if(t_lattice_model) then
                                    ASSERT(.not. t_heisenberg_model)
                                    MatEl = get_helement_lattice(nJ, ic, temp_ex, tSign)
                                else
                                    MatEl = dyn_sltcnd_excit_old(nJ, IC, temp_ex, tSign)
                                end if
                            else
                                temp_ex(1, :) = ExcitMat(2, :)
                                temp_ex(2, :) = ExcitMat(1, :)
                                if(t_lattice_model) then
                                    ASSERT(.not. t_heisenberg_model)
                                    MatEl = get_helement_lattice(nJ, ic, temp_ex, tSignOrig)
                                else
                                    MatEl = dyn_sltcnd_excit_old(nJ, IC, temp_ex, tSignOrig)
                                end if
                            end if
                            HEl = MatEl * SQRT(2.0_dp)
                        end if
                    ELSE     !Open shell -> Open shell

!First find nI -> nJ. If nJ has swapped, then this will be different.
                        if(tSwapped) then
                            temp_ex(1, :) = ex2(2, :)
                            temp_ex(2, :) = ex2(1, :)
                            if(t_lattice_model) then
                                MatEl = get_helement_lattice(nJ, ExcitLevel, temp_ex, tSign)
                            else
                                MatEl = dyn_sltcnd_excit_old(nJ, ExcitLevel, temp_ex, tSign)
                            end if
                        else
                            temp_ex(1, :) = ExcitMat(2, :)
                            temp_ex(2, :) = ExcitMat(1, :)
                            if(t_lattice_model) then
                                MatEl = get_helement_lattice(nJ, ic, temp_ex, tSignOrig)
                            else
                                MatEl = dyn_sltcnd_excit_old(nJ, IC, temp_ex, tSignOrig)
                            end if

                        end if

                        !now nI2 -> nJ (modelled as nI -> nJ2 with appropriate sign modifications)

                        IF((ExcitLevel == 2) .or. (ExcitLevel == 1)) THEN

                            CALL CalcOpenOrbs(iLutnJ, OpenOrbsJ)
                            CALL CalcOpenOrbs(iLutnI, OpenOrbsI)

                            IF(tSwapped) THEN
                                IF((OpenOrbsJ + OpenOrbsI) == 3) tSignOrig = .not. tSignOrig
                                !I.e. J odd and I even or vice versa, but since these can only be at max quads, then they can only have 1/2 open orbs
                                temp_ex(1, :) = ExcitMat(2, :)
                                temp_ex(2, :) = ExcitMat(1, :)

                                if(t_lattice_model) then
                                    ! here i want to get nI -> nJ2
                                    ! when it was swapped the original
                                    ! ExcitMat and tSignOrig are associated
                                    ! with the excitation
                                    MatEl2 = get_helement_lattice(nJ2, ic, temp_ex, tSignOrig)
                                else
                                    MatEl2 = dyn_sltcnd_excit_old(nJ2, IC, temp_ex, tSignOrig)
                                end if
                            ELSE
!I.e. J odd and I even or vice versa, but since these can only be at max quads, then they can only have 1/2 open orbs
                                IF((OpenOrbsJ + OpenOrbsI) == 3) tSign = .not. tSign
                                temp_ex(1, :) = ex2(2, :)
                                temp_ex(2, :) = ex2(1, :)

                                if(t_lattice_model) then
                                    ! if they were not swapped Ex2 and tSign
                                    ! are associated with nI -> nJ2
                                    MatEl2 = get_helement_lattice(nJ2, ExcitLevel, temp_ex, tSign)
                                else
                                    MatEl2 = dyn_sltcnd_excit_old(nJ2, ExcitLevel, temp_ex, tSign)
                                end if
                            end if

                            IF(tOddS_HPHF) THEN
!again, since these can only be at max quads, then they can only have 1/2 open orbs...
                                IF(OpenOrbsI == 2) THEN
                                    MatEl = MatEl - MatEl2
                                ELSE
                                    MatEl = MatEl + MatEl2
                                end if
                            ELSE
!again, since these can only be at max quads, then they can only have 1/2 open orbs...
                                IF(OpenOrbsI == 2) THEN
                                    MatEl = MatEl + MatEl2
                                ELSE
                                    MatEl = MatEl - MatEl2
                                end if
                            end if
                        end if
                        HEl = MatEl

                    end if   !Endif from open/closed shell det
                    if(IC /= 0 .and. modk_offdiag) hel = -abs(hel)

                end if   !Endif want to generate matrix element

!Here, we actually know nJ, so don't need to regenerate it...

            else if(ExcitLevel == 0) THEN
!We have generated the same HPHF. MatEl wants to be zero.
                nJ(1) = 0
                IF(tGenMatHEl) THEN
                    HEl = 0.0_dp
                end if

            ELSE    !Open-shell to Open-shell, but with no cross-connection.

                IF(tGenMatHEl) THEN
!iLutnI MUST be open-shell here, since otherwise it would have been connected to
!iLutnJ2. Also, we know the cross connection (i.e. MatEl2 = 0)
                    ! WD: Here I am not 100% sure if I always take ExcitMat..
                    temp_ex(1, :) = ExcitMat(2, :)
                    temp_ex(2, :) = ExcitMat(1, :)
                    IF(tSwapped) THEN
                        CALL CalcOpenOrbs(iLutnJ, OpenOrbsJ)
                        IF(tOddS_HPHF) then
                            IF(mod(OpenOrbsJ, 2) == 0) THEN
                                tSignOrig = .not. tSignOrig
                            end if
                        ELSE
                            IF(mod(OpenOrbsJ, 2) == 1) THEN
                                tSignOrig = .not. tSignOrig
                            end if
                        end if

                        if(t_lattice_model) then
                            MatEl = get_helement_lattice(nJ2, ic, temp_ex, tSignOrig)
                        else
                            MatEl = dyn_sltcnd_excit_old(nJ2, IC, temp_ex, tSignOrig)
                        end if
                    ELSE
                        if(t_lattice_model) then
                            MatEl = get_helement_lattice(nJ, ic, temp_ex, tSignOrig)
                        else
                            MatEl = dyn_sltcnd_excit_old(nJ, IC, temp_ex, tSignOrig)
                        end if
                    end if

                    HEl = MatEl
                    if(IC /= 0 .and. modk_offdiag) hel = -abs(hel)

                end if

            end if

        end if

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
    SUBROUTINE ReturnAlphaOpenDet(nI, nJ, iLutnI, iLutSym, tCalciLutSym, tCalcnISym, tSwapped)
        INTEGER(KIND=n_int), intent(inout) :: iLutSym(0:NIfTot), iLutnI(0:NIfTot)
        integer(kind=n_int) :: iLutTemp(0:NIfTot)
        INTEGER :: i, nTemp(NEl)
        integer, intent(inout) :: nJ(NEl), nI(NEl)
        LOGICAL, intent(in) :: tCalciLutSym, tCalcnISym
        logical, intent(out) :: tSwapped

        IF(tCalciLutSym) THEN
            CALL FindExcitBitDetSym(iLutnI, iLutSym)
        end if
        IF(tCalcnISym) THEN
            CALL FindDetSpinSym(nI, nJ, NEl)
        end if

        ! iLutnI is 'less' than iLutSym, so iLutSym is the determinant with
        ! the first open-shell = alpha. Swap them around.
        ! Only count up to NIfD to avoid Yamanouchi symbol etc.
        i = DetBitLT(iLutnI, iLutSym, NIfD)
        IF(i == 1) THEN
            iLutTemp(:) = iLutnI(:)
            iLutnI(:) = iLutSym(:)
            iLutSym(:) = iLutTemp(:)
!            CALL FindDetSpinSym(nI,nJ,NEl)
            nTemp(:) = nI(:)
            nI(:) = nJ(:)
            nJ(:) = nTemp(:)
            tSwapped = .true.
        else if(i == 0) THEN
            CALL Stop_All("ReturnAlphaOpenDet", "Shouldn't have closed shell determinants in here")
        ELSE
            tSwapped = .false.
        end if

    END SUBROUTINE ReturnAlphaOpenDet

!This create the spin-coupled determinant of nI in nJ in natural ordered form.
    PURE SUBROUTINE FindDetSpinSym(nI, nJ, NEl)
        INTEGER, intent(in) :: NEl, nI(NEl)
        integer, intent(out) :: nJ(NEl)
        integer :: i


        ! for debug compilation treat first entry seperately
        if (is_alpha(nI(1))) then
            nJ(1) = nI(1) - 1
        else
            if (get_alpha(nI(1)) /= nI(2)) then
                nJ(1) = nI(1) + 1
            else
                nJ(1) = nI(1)
            end if
        end if

        do i = 2, nel
            ! If electron is an alpha electron, change it to a beta (unless
            ! it is part of a closed pair of electrons).
            if(is_alpha(nI(i))) then
                if(i == 1) then
                    nJ(i) = nI(i) - 1
                else if(get_beta(nI(i)) /= nI(i - 1)) then
                    nJ(i) = nI(i) - 1
                else
                    nJ(i) = nI(i)
                end if
                ! vice-versa for beta.
            else
                if(i == nel) then
                    nJ(i) = nI(i) + 1
                else if(get_alpha(nI(i)) /= nI(i + 1)) then
                    nJ(i) = nI(i) + 1
                else
                    nJ(i) = nI(i)
                end if
            end if
        end do

    END SUBROUTINE FindDetSpinSym

!In closed-shell systems with equal number of alpha and beta strings, the amplitude of a
!determinant in the final CI wavefunction is the same
!when the alpha and beta electrons are swapped (for S=0, see Helgakker for more details).
!It will sometimes be necessary to find this other
!determinant when spawning. This routine will find the bit-representation of an excitation
!by constructing the symmetric iLut from the its
!symmetric partner, also in bit form.
    PURE SUBROUTINE FindExcitBitDetSym(iLut, iLutSym)
        IMPLICIT NONE
        INTEGER(KIND=n_int), intent(in) :: iLut(0:NIfTot)
        INTEGER(KIND=n_int), intent(out) :: iLutSym(0:NIfTot)
        INTEGER(KIND=n_int) :: iLutAlpha(0:NIfTot), iLutBeta(0:NIfTot)
        INTEGER :: i

        iLutSym(:) = 0
        iLutAlpha(:) = 0
        iLutBeta(:) = 0

        do i = 0, NIfD

            iLutAlpha(i) = IAND(iLut(i), MaskAlpha)    !Seperate the alpha and beta bit strings
            iLutBeta(i) = IAND(iLut(i), MaskBeta)

            iLutAlpha(i) = ISHFT(iLutAlpha(i), -1)  !Shift all alpha bits to the left by one.
            iLutBeta(i) = ISHFT(iLutBeta(i), 1)   !Shift all beta bits to the right by one.

            iLutSym(i) = IOR(iLutAlpha(i), iLutBeta(i))    !Combine the bit strings to give the final bit representation.

        end do

    END SUBROUTINE FindExcitBitDetSym

!!This routine will take a HPHF nI, and find Iterations number of excitations of it.
!It will then histogram these, summing in 1/pGen for every occurance of
!!the excitation. This means that all excitations should be 0 or 1 after enough iterations.
!It will then count the excitations and compare the number to the
!!number of excitations generated using the full enumeration excitation generation.

    SUBROUTINE BinSearchListHPHF(iLut, List, Length, MinInd, MaxInd, PartInd, tSuccess)
        INTEGER :: Length, MinInd, MaxInd, PartInd
        INTEGER(KIND=n_int) :: iLut(0:NIfTot), List(0:NIfTot, Length)
        INTEGER :: i, j, N, Comp
        LOGICAL :: tSuccess

        i = MinInd
        j = MaxInd
        IF(i - j == 0) THEN
            Comp = DetBitLT(List(:, MaxInd), iLut(:), nifd)
            IF(Comp == 0) THEN
                tSuccess = .true.
                PartInd = MaxInd
                RETURN
            ELSE
                tSuccess = .false.
                PartInd = MinInd
            end if
        end if
        do while(j - i > 0)  !End when the upper and lower bound are the same.
            N = (i + j) / 2       !Find the midpoint of the two indices

!Comp is 1 if CyrrebtDets(N) is "less" than iLut, and -1 if it is more or 0 if they are the same
            Comp = DetBitLT(List(:, N), iLut(:), nifd)

            IF(Comp == 0) THEN
!Praise the lord, we've found it!
                tSuccess = .true.
                PartInd = N
                RETURN
            else if((Comp == 1) .and. (i /= N)) THEN
!The value of the determinant at N is LESS than the determinant we're looking for.
!Therefore, move the lower bound of the search up to N.
!However, if the lower bound is already equal to N then the two bounds are consecutive and we have failed...
                i = N
            else if(i == N) THEN

                IF(i == MaxInd - 1) THEN
!This deals with the case where we are interested in the final/first entry in the list. Check the final entry of the list and leave
!We need to check the last index.
                    Comp = DetBitLT(List(:, i + 1), iLut(:), nifd)
                    IF(Comp == 0) THEN
                        tSuccess = .true.
                        PartInd = i + 1
                        RETURN
                    else if(Comp == 1) THEN
!final entry is less than the one we want.
                        tSuccess = .false.
                        PartInd = i + 1
                        RETURN
                    ELSE
                        tSuccess = .false.
                        PartInd = i
                        RETURN
                    end if

                else if(i == MinInd) THEN
                    tSuccess = .false.
                    PartInd = i
                    RETURN
                ELSE
                    i = j
                end if

            else if(Comp == -1) THEN
!The value of the determinant at N is MORE than the determinant we're looking for. Move the upper bound of the search down to N.
                j = N
            ELSE
!We have failed - exit loop
                i = j
            end if

        end do

!If we have failed, then we want to find the index that is one less than where the particle would have been.
        tSuccess = .false.
        PartInd = MAX(MinInd, i - 1)

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
        use bit_rep_data, only: test_flag, IlutBits

        integer, intent(in) :: nI(nel), ex(2, maxExcit), ic
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ClassCount2(ScratchSize)
        integer, intent(in) :: ClassCountUnocc2(ScratchSize)
        real(dp), intent(in) :: pDoub
        real(dp), intent(out) :: pGen
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = 'CalcNonUniPGen'

        integer :: temp_part_type
        integer(n_int) :: ilutJ(0:IlutBits%len_tot)

        ! We need to consider which of the excitation generators are in use,
        ! and call the correct routine in each case.

        pgen = 0.0_dp

        ! i have to make sure to catch all this call to this function correctly
        if(present(part_type)) then
            temp_part_type = part_type
        else
            temp_part_type = 1
        end if

        ! does it help to avoid recalculating for the reference?
        ! do i need to  check if it is actually a non-initiator?
        ! i guess i do.. or i go the unnecessary way of checking again in
        ! the called back-spawn functions
        if(ic == 3) then
            pgen = calc_pgen_triple(nI, ex)
        else if((t_back_spawn .or. t_back_spawn_flex) .and. &
                (.not. DetBitEq(ilutI, ilutRef(:, temp_part_type), nifd)) .and. &
                (.not. test_flag(ilutI, get_initiator_flag(temp_part_type)))) then
            ! i just realised this also has to be done for the hubbard
            ! and the ueg model.. -> create those functions!
            if(tHUB .and. tLatticeGens) then
                pgen = calc_pgen_back_spawn_hubbard(nI, ilutI, ex, ic, temp_part_type)
            else if(tUEGNewGenerator .and. tLatticeGens) then
                pgen = calc_pgen_back_spawn_ueg_new(nI, ilutI, ex, ic, temp_part_type)
            else if(tUEG .and. tLatticeGens) then
                pgen = calc_pgen_back_spawn_ueg(ilutI, ex, ic, temp_part_type)
            else
                pgen = calc_pgen_back_spawn(nI, ilutI, ex, ic, temp_part_type)
            end if
        else if(tGen_4ind_2) then
            pgen = calc_pgen_4ind_weighted2(nI, ilutI, ex, ic)
        else if(tGen_4ind_weighted) then
            pgen = calc_pgen_4ind_weighted(nI, ilutI, ex, ic, &
                                           ClassCountUnocc2)
        else if(tGen_4ind_reverse) then
            pgen = calc_pgen_4ind_reverse(nI, ilutI, ex, ic)

            ! this if construct is not well setup.. this can fail..
        else
            if(tLatticeGens) then
                if(ic == 2) then
                    call CalcPGenLattice(ex, pGen)
                else
                    pGen = 0
                end if
            else if(tGen_4ind_2) then
                pgen = calc_pgen_4ind_weighted2(nI, ilutI, ex, ic)

            else if(tGen_4ind_weighted) then
                pgen = calc_pgen_4ind_weighted(nI, ilutI, ex, ic, &
                                               ClassCountUnocc2)
            else if(tGen_4ind_reverse) then
                pgen = calc_pgen_4ind_reverse(nI, ilutI, ex, ic)

            else if(t_new_real_space_hubbard) then
                if(t_trans_corr_hop) then
                    if(t_uniform_excits) then
                        pgen = calc_pgen_rs_hubbard_transcorr_uniform(ex, ic)
                    else
                        pgen = calc_pgen_rs_hubbard_transcorr(nI, ilutI, ex, ic)
                    end if
                else if(t_spin_dependent_transcorr) then
                    pgen = calc_pgen_rs_hubbard_spin_dependent_transcorr(nI, ilutI, ex, ic)
                else
                    pgen = calc_pgen_rs_hubbard(ilutI, ex, ic)
                end if

            else if(t_tJ_model) then
                pgen = calc_pgen_tJ_model(ilutI, ex, ic)

            else if(t_heisenberg_model) then
                pgen = calc_pgen_heisenberg_model(ilutI, ex, ic)

            else if(t_k_space_hubbard) then
                ! change with Kais uniform excitgen implementation
                if(t_uniform_excits) then
                    if(ic == 2) then
                        call CalcPGenLattice(ex, pgen)
                    else
                        pgen = 0.0_dp
                    end if
                else
                    pgen = calc_pgen_k_space_hubbard(nI, ilutI, ex, ic)
                end if
            else if (t_guga_pchb) then
                pgen = calc_pgen_guga_pchb(ilutI, current_csf_i, ilutJ)
            else if(t_pcpp_excitgen) then
                pgen = calc_pgen_pcpp(ilutI, ex, ic)
            else if (allocated(current_exc_generator)) then
                pgen = current_exc_generator%get_pgen(nI, ilutI, ex, ic, ClassCount2, ClassCountUnocc2)
            else
                ! Here we assume that the normal excitation generators in
                ! symrandexcit2.F90 are being used.
                call calc_pgen_symrandexcit2(nI, ex, ic, ClassCount2, &
                                             ClassCountUnocc2, pDoub, pGen)
            end if
        end if

    end subroutine

END MODULE HPHFRandExcitMod
