#include "macros.h"

module back_spawn_excit_gen

    use constants, only: dp, n_int, EPS, bits_n_int, maxExcit, stdout
    use SystemData, only: nel, G1, nbasis, tHPHF, NMAXX, NMAXY, NMAXZ, &
                          tOrbECutoff, OrbECutoff, nOccBeta, nOccAlpha, ElecPairs, &
                          tHPHF
    use bit_rep_data, only: niftot
    use SymExcitDataMod, only: excit_gen_store_type, SpinOrbSymLabel, &
                               kPointToBasisFn
    use bit_reps, only: test_flag, get_initiator_flag
    use FciMCData, only: pSingles, projedet, pDoubles
    use dSFMT_interface, only: genrand_real2_dSFMT
    use excit_gens_int_weighted, only: gen_single_4ind_ex, select_orb_sing, &
                                       pick_weighted_elecs, select_orb, &
                                       pgen_select_orb, pgen_weighted_elecs, pgen_single_4ind
    use excit_gen_5, only: gen_double_4ind_ex2, pick_a_orb, pgen_select_a_orb, &
                           calc_pgen_4ind_weighted2
    use CalcData, only: t_back_spawn_flex, t_back_spawn_occ_virt, t_back_spawn, &
                        occ_virt_level, t_back_spawn_flex_option, t_back_spawn_option
    use GenRandSymExcitNUMod, only: ClassCountInd, RandExcitSymLabelProd, &
                                    CreateExcitLattice, CalcPGenLattice
    use back_spawn, only: check_electron_location, pick_virtual_electrons_double, &
                          pick_occupied_orbital_single, pick_virtual_electron_single, &
                          pick_occupied_orbital, pick_second_occupied_orbital, &
                          is_in_ref, pick_occupied_orbital_ueg, &
                          pick_virtual_electrons_double_hubbard, pick_occupied_orbital_hubbard, &
                          is_allowed_ueg_k_vector
    use get_excit, only: make_single, make_double
    use Determinants, only: write_det, get_helement
    use ueg_excit_gens, only: gen_double_ueg, create_ab_list_ueg, pick_uniform_elecs, &
                              calc_pgen_ueg
    use util_mod, only: operator(.div.)

    use util_mod_numerical, only: binary_search_first_ge

    use lattice_models_utils, only: make_ilutJ, get_orb_from_kpoints, get_ispn

    use excit_gens_int_weighted, only: get_paired_cc_ind

#ifdef DEBUG_
    use SystemData, only: tNoFailAb
#endif

    implicit none

contains

    subroutine gen_excit_back_spawn_ueg_new(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                            ExcitMat, tParity, pgen, HelGen, store, part_type)
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ic, ExcitMat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = "gen_excit_back_spawn_ueg_new"

        integer :: iUnused

#ifdef DEBUG_
        real(dp) :: pgen2
        HElement_t(dp) :: temp_hel
#endif

        HelGen = 0.0_dp
        iUnused = exFlag
        iUnused = store%nopen

        ic = 2

        ! do i want to implement both old and new back spawn??
        ! no!
        if ((t_back_spawn_flex .or. t_back_spawn) .and. &
            (.not. test_flag(ilutI, get_initiator_flag(part_type)))) then

            call gen_double_back_spawn_ueg_new(nI, ilutI, part_type, nJ, ilutJ, tParity, &
                                               ExcitMat, pgen)

#ifdef DEBUG_
            if (.not. IsNullDet(nJ)) then
                pgen2 = calc_pgen_back_spawn_ueg_new(nI, ilutI, ExcitMat, ic, part_type)
                if (abs(pgen - pgen2) > 1.0e-6_dp) then
                    if (tHPHF) then
                        print *, "due to circular dependence, no matrix element calc possible!"
                        temp_hel = 0.0_dp
                    else
                        temp_hel = get_helement(nI, nJ, ilutI, ilutJ)
                    end if

                    write(stdout, *) 'Calculated and actual pgens differ. for non-initiator'
                    write(stdout, *) 'This will break HPHF calculations'
                    write(stdout, *) 'reference det: '
                    call write_det(6, projedet(:, part_type_to_run(part_type)), .true.)
                    call write_det(6, nI, .false.)
                    write(stdout, '(" --> ")', advance='no')
                    call write_det(6, nJ, .true.)
                    write(stdout, *) 'Excitation matrix: ', ExcitMat(1, 1:ic), '-->', &
                        ExcitMat(2, 1:ic)
                    write(stdout, *) 'Generated pGen:  ', pgen
                    write(stdout, *) 'Calculated pGen: ', pgen2
                    write(stdout, *) 'matrix element: ', temp_hel
                    call stop_all(this_routine, "Invalid pGen")
                end if
            end if
#endif
        else

            call gen_double_ueg(nI, ilutI, nJ, ilutJ, tParity, ExcitMat, pgen)

#ifdef DEBUG_
            if (.not. IsNullDet(nJ)) then
                pgen2 = calc_pgen_ueg(ilutI, ExcitMat, ic)
                if (abs(pgen - pgen2) > 1.0e-6_dp) then
                    if (tHPHF) then
                        print *, "due to circular dependence, no matrix element calc possible!"
!                         temp_hel = hphf_off_diag_helement(nI,nJ,ilutI,ilutJ)
                        temp_hel = 0.0_dp
                    else
                        temp_hel = get_helement(nI, nJ, ilutI, ilutJ)
                    end if

                    write(stdout, *) 'Calculated and actual pgens differ. for non-initiator'
                    write(stdout, *) 'This will break HPHF calculations'
                    call write_det(6, nI, .false.)
                    write(stdout, '(" --> ")', advance='no')
                    call write_det(6, nJ, .true.)
                    write(stdout, *) 'Excitation matrix: ', ExcitMat(1, 1:ic), '-->', &
                        ExcitMat(2, 1:ic)
                    write(stdout, *) 'Generated pGen:  ', pgen
                    write(stdout, *) 'Calculated pGen: ', pgen2
                    write(stdout, *) 'matrix element: ', temp_hel
                    call stop_all(this_routine, "Invalid pGen")
                end if
            end if
#endif
        end if

    end subroutine gen_excit_back_spawn_ueg_new

    subroutine gen_double_back_spawn_ueg_new(nI, ilutI, part_type, nJ, ilutJ, tPar, ex, pgen)
        integer, intent(in) :: nI(nel), part_type
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tPar
        real(dp), intent(out) :: pgen

        integer :: elecs(2), src(2), ispn, sum_ml, loc, orb_a, orb_b
        real(dp) :: p_elec, p_orb, dummy, cum_arr(nBasis), cum_sum

        logical :: t_temp_back_spawn

        t_temp_back_spawn = .false.

        if (t_back_spawn) then
            call pick_virtual_electrons_double(nI, part_type, elecs, src, ispn, &
                                               sum_ml, p_elec)

            loc = -1
        else
            call pick_uniform_elecs(elecs, p_elec)
            src = nI(elecs)

            ispn = get_ispn(src)

            loc = check_electron_location(src, 2, part_type)
        end if

        if (elecs(1) == 0) then
            nJ(1) = 0
            pgen = 0.0_dp
            return
        end if

        if ((loc == 2) .or. (loc == 1 .and. occ_virt_level /= -1) .or. &
            (loc == 0 .and. occ_virt_level >= 1)) then

            t_temp_back_spawn = .true.
            call pick_occupied_orbital_ueg(ilutI, src, iSpn, part_type, p_orb, &
                                           dummy, orb_a)
            ! it can happen that there are no valid orbitals
            if (orb_a == 0) then
                nJ(1) = 0
                pgen = 0.0_dp
                return
            end if

            ! i think in this case i have to multiply if both are in the
            ! reference or?

        else

            ! i could refactor that in a smaller function:
            call create_ab_list_ueg(ilutI, src, cum_arr, cum_sum)

            if (cum_sum < EPS) then
                nJ(1) = 0
                pgen = 0.0_dp
                return
            end if

            orb_a = binary_search_first_ge(cum_arr, genrand_real2_dsfmt() * cum_sum)

            if (orb_a == 1) then
                p_orb = cum_arr(1) / cum_sum
            else
                p_orb = (cum_arr(orb_a) - cum_arr(orb_a - 1)) / cum_sum
            end if

        end if

        orb_b = get_orb_from_kpoints(src(1), src(2), orb_a)

        ! thats it i guess..
        ! but i might have to be careful:
        ! in the back-spawn mechanism any order of the orbitals is possible..
        ! in the new ueg by NB not.. there orb_b > orb_a is strictly
        ! enforced.. hm.. how should we handle that?
        ! wait a minute.. since when is this called with the bare
        ! eleci and elecj?? yes it is apparently..
        call make_double(nI, nJ, elecs(1), elecs(2), orb_a, orb_b, ex, tpar)

        ilutJ = make_ilutJ(ilutI, ex, 2)

        pgen = p_elec * p_orb

        if (t_temp_back_spawn .and. is_in_ref(orb_b, part_type)) then
            pgen = 2.0_dp * pgen
        end if

    end subroutine gen_double_back_spawn_ueg_new

    function calc_pgen_back_spawn_ueg_new(nI, ilutI, ex, ic, part_type) result(pgen)
        ! i also need immmidiately a calc_pgen function!
        integer, intent(in) :: nI(nel), ex(2, 2), ic, part_type
        integer(n_int), intent(in) :: ilutI(0:niftot)
        real(dp) :: pgen

        integer :: elecs(2), src(2), ispn, sum_ml, dummy_src(2), dumm_iSpn
        integer :: orb_a, tgt(2), loc
        real(dp) :: p_elec, p_orb, dummy, cum_arr(nBasis), cum_sum

        if (ic /= 2) then
            pgen = 0.0_dp
            return
        end if

        ! and maybe i should also enable that i call this outside of
        ! knowledge of the initator status..
        if (test_flag(ilutI, get_initiator_flag(part_type))) then
            pgen = calc_pgen_ueg(ilutI, ex, ic)
        else

            src = get_src(ex)
            tgt = get_tgt(ex)
            ispn = get_ispn(src)

            if (t_back_spawn) then
                call pick_virtual_electrons_double(nI, part_type, elecs, dummy_src, &
                                                   dumm_iSpn, sum_ml, p_elec, .true.)
                loc = -1
            else
                p_elec = 1.0_dp / real(ElecPairs, dp)
                loc = check_electron_location(src, 2, part_type)
            end if

            if ((loc == 2) .or. (loc == 1 .and. occ_virt_level /= -1) .or. &
                (loc == 0 .and. occ_virt_level >= 1)) then
                ! argh.. wait a minute.. i have to ensure that i only do that
                ! for back-spawn flex!

                call pick_occupied_orbital_ueg(ilutI, src, iSpn, part_type, p_orb, &
                                               dummy, orb_a, .true.)

                ! do i need to multiply if both are in the reference?
                ! i guess so.. since then i could have picked it in both
                ! orders..
                if (is_in_ref(tgt(1), part_type) .and. is_in_ref(tgt(2), part_type)) then
                    p_orb = 2.0_dp * p_orb
                end if
            else

                ! i could refactor that in a smaller function:
                call create_ab_list_ueg(ilutI, src, cum_arr, cum_sum)

                tgt = get_tgt(ex)

                orb_a = tgt(1)

                if (orb_a == 1) then
                    p_orb = cum_arr(orb_a) / cum_sum
                else
                    p_orb = (cum_arr(orb_a) - cum_arr(orb_a - 1)) / cum_sum
                end if
            end if

            pgen = p_orb * p_elec

        end if

    end function calc_pgen_back_spawn_ueg_new

    subroutine gen_excit_back_spawn_hubbard(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                            ExcitMat, tParity, pgen, HelGen, store, part_type)
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ic, ExcitMat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = "gen_excit_back_spawn_hubbard"

        integer :: iUnused
#ifdef DEBUG_
        real(dp) :: pgen2
        HElement_t(dp) :: temp_hel
#endif
        ! so now pack everything necessary for a ueg excitation generator.

        ! why is ilutJ not created here? do it!
        HelGen = 0.0_dp
        iUnused = exFlag
        iUnused = store%nopen

        ic = 2

        ! this function gets pointed to if tUEG, t_back_spawn_flex and tLatticeGens
        ! is set
        ! BUT: again: we also need to take care of the HPHF keyword..
        ! which is a mess

        ! implement it for now without this tNoFailAb flag
        ASSERT(.not. tNoFailAb)

        if ((t_back_spawn .or. t_back_spawn_flex) .and. .not. &
            test_flag(ilutI, get_initiator_flag(part_type))) then

            call gen_double_back_spawn_hubbard(nI, ilutI, nJ, ilutJ, tParity, ExcitMat, &
                                               pgen)

#ifdef DEBUG_
            if (.not. IsNullDet(nJ)) then
                pgen2 = calc_pgen_back_spawn_hubbard(nI, ilutI, ExcitMat, ic, part_type)
                if (abs(pgen - pgen2) > 1.0e-6_dp) then
                    if (tHPHF) then
                        print *, "due to circular dependence, no matrix element calc possible!"
!                         temp_hel = hphf_off_diag_helement(nI,nJ,ilutI,ilutJ)
                        temp_hel = 0.0_dp
                    else
                        temp_hel = get_helement(nI, nJ, ilutI, ilutJ)
                    end if

                    write(stdout, *) 'Calculated and actual pgens differ. for non-initiator'
                    write(stdout, *) 'This will break HPHF calculations'
                    call write_det(6, nI, .false.)
                    write(stdout, '(" --> ")', advance='no')
                    call write_det(6, nJ, .true.)
                    write(stdout, *) 'Excitation matrix: ', ExcitMat(1, 1:ic), '-->', &
                        ExcitMat(2, 1:ic)
                    write(stdout, *) 'Generated pGen:  ', pgen
                    write(stdout, *) 'Calculated pGen: ', pgen2
                    write(stdout, *) 'matrix element: ', temp_hel
                    call stop_all(this_routine, "Invalid pGen")
                end if
            end if
#endif
        else
            ! do i want to rewrite the old on or just reuse? for now reuse:
            call CreateExcitLattice(nI, ilutI, nJ, tParity, ExcitMat, pgen, part_type)

            if (.not. IsNullDet(nJ)) ilutJ = make_ilutJ(ilutI, ExcitMat, ic)

#ifdef DEBUG_
            if (.not. IsNullDet(nJ)) then
                call CalcPGenLattice(ExcitMat, pgen2)
                if (abs(pgen - pgen2) > 1.0e-6_dp) then
                    if (tHPHF) then
                        print *, "due to circular dependence, no matrix element calc possible!"
!                         temp_hel = hphf_off_diag_helement(nI,nJ,ilutI,ilutJ)
                        temp_hel = 0.0_dp
                    else
                        temp_hel = get_helement(nI, nJ, ilutI, ilutJ)
                    end if

                    write(stdout, *) 'Calculated and actual pgens differ. for non-initiator'
                    write(stdout, *) 'This will break HPHF calculations'
                    call write_det(6, nI, .false.)
                    write(stdout, '(" --> ")', advance='no')
                    call write_det(6, nJ, .true.)
                    write(stdout, *) 'Excitation matrix: ', ExcitMat(1, 1:ic), '-->', &
                        ExcitMat(2, 1:ic)
                    write(stdout, *) 'Generated pGen:  ', pgen
                    write(stdout, *) 'Calculated pGen: ', pgen2
                    write(stdout, *) 'matrix element: ', temp_hel
                    call stop_all(this_routine, "Invalid pGen")
                end if
            end if
#endif
        end if

    end subroutine gen_excit_back_spawn_hubbard

    subroutine gen_double_back_spawn_hubbard(nI, ilutI, nJ, ilutJ, tParity, ExcitMat, &
                                             pgen, part_type)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ExcitMat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        integer, optional :: part_type
        character(*), parameter :: this_routine = "gen_double_back_spawn_hubbard"

        integer :: elec_i, elec_j, iSpn, orb_b, src(2), elecs(2), &
                   loc, temp_part_type, orb_a
        real(dp) :: x, pAIJ, dummy, mult, pgen_elec
        logical :: tAllowedExcit, t_temp_back_spawn

        ! damn.. remember if i initialize stuff above it implicitly assumes
        ! the (save) attribut!
        t_temp_back_spawn = .false.

        if (present(part_type)) then
            temp_part_type = part_type
        else
            temp_part_type = 1
        end if

        ! i know here that back-spawn is on and this is a non-inititator:
        ! so for the beginning implementation we pick the two electrons
        ! randomly
        ! i should make this to a function below: maybe with the hubbard
        ! flag as additional input to provide us with two independent
        ! electrons in any order possible!
        ! maybe for now.. i still want to retain the original
        ! back-spawn functionality
        if (t_back_spawn) then
            call pick_virtual_electrons_double_hubbard(nI, temp_part_type, elecs, src, ispn, &
                                                       pgen_elec)
            ! check if enogh electrons are in the virtual
            if (elecs(1) == 0) then
                nJ(1) = 0
                pgen = 0.0_dp
                return
            end if

            elec_i = elecs(1)
            elec_j = elecs(2)

            loc = -1
        else

            do
                elec_i = 1 + int(genrand_real2_dsfmt() * nel)
                do
                    elec_j = 1 + int(genrand_real2_dsfmt() * nel)
                    if (elec_j /= elec_i) exit
                end do
                ispn = get_ispn([nI(elec_i), nI(elec_j)])

                ! i need opposite spin-excitations in the hubbard model
                if (iSpn == 2) exit
            end do

            ! do i need 2* here?
            pgen_elec = 1.0_dp / real(nOccBeta * nOccAlpha, dp)

            src = nI([elec_i, elec_j])

            loc = check_electron_location(src, 2, temp_part_type)
        end if

        ! i need to call nI(elecs)
        ! this test should be fine since, if back_spawn is active loc should
        ! always be 0 so we are not risking of picking occupied orbitals
        ! below..

        ! wait a minute.. thats incorrect or?
        ! if we have both in the occupied manifold we want to restrict
        ! our orbital choice..
        ! i can decide based on the occ_virt_level or?
        if ((loc == 2) .or. (loc == 1 .and. occ_virt_level /= -1) .or. &
            (loc == 0 .and. occ_virt_level >= 1)) then
            t_temp_back_spawn = .true.
            ! i think i need a new one for the ueg also.. since there is not
            ! such a spin restriction as in the hubbard model
            ! i think just using the "normal" occupied picker should be fine
            ! how does this compile even??
            ! no.. write a new one without the spin-restriction
            call pick_occupied_orbital_hubbard(ilutI, temp_part_type, pAIJ, orb_a)
            ! i should take k-space symmetry into accound in picking
            ! orb a

            ! i guess we can have no possible excitations..
            if (orb_a == 0) then
                nJ(1) = 0
                pgen = 0.0_dp
                return
            end if
        else
            ! otherwise pick freely..
            ! in the hubbard case we always have iSpn = 2
            do
                orb_a = int(nBasis * genrand_real2_dsfmt()) + 1

                ! i guess it is more efficient to check if the momentum is
                ! correctly conserved here..

                if (IsNotOcc(ilutI, orb_a)) exit
            end do

            paij = 1.0_dp / real(nBasis - nel, dp)

        end if

        ! i guess this is not necessary, but anyway..
        IF ((orb_a < 0) .or. (orb_a > nBasis)) THEN
            CALL Stop_All("CreateDoubExcitLattice", "Incorrect basis function generated")
        end if

!         ! Is kb allowed by the size of the space?
!         ! Currently only applies when NMAXX etc. are set by the CELL keyword
!         tAllowedExcit = is_allowed_ueg_k_vector(src(1), src(2), orb_a)
!
!         IF(.not.tAllowedExcit) THEN
!             nJ(1)=0
!             RETURN
!         end if

        orb_b = get_orb_from_kpoints(src(1), src(2), orb_a)

        IF (orb_b <= 0 .or. orb_a == orb_b) THEN
            nJ(1) = 0
            pgen = 0.0_dp
            RETURN
        end if

        ! Is b occupied?
        if (IsOcc(ilutI, orb_b)) then
            nj(1) = 0
            pgen = 0.0_dp
            return
        end if

        ! Find the new determinant
        call make_double(nI, nJ, elec_i, elec_j, orb_a, &
                         orb_b, ExcitMat, tParity)

        ilutJ = make_ilutJ(ilutI, ExcitMat, 2)

        ! we knoe it is back-spawn + hubbard remember! so paij is already above
        !calculate generation probabilities
        ! note, p(b|ij)=p(a|ij) for this system
        ! if b is also in the occupied manifold in the case of back_spawn_flex
        mult = 2.0_dp
        if (t_temp_back_spawn .and. (.not. is_in_ref(orb_b, temp_part_type))) then
            mult = 1.0_dp
        end if

        pgen = mult * pgen_elec * pAIJ

    end subroutine gen_double_back_spawn_hubbard

    function calc_pgen_back_spawn_hubbard(nI, ilutI, ex, ic, part_type) result(pgen)
        integer, intent(in) :: nI(nel), ex(2, 2), ic, part_type
        integer(n_int), intent(in) :: ilutI(0:niftot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = "calc_pgen_back_spawn"

        integer :: d_elecs(2), d_src(2), d_ispn, src(2), loc, tgt(2), d_orb
        real(dp) :: pgen_elec, paij, mult

        ! the hubbard pgen recalculation is pretty similar to the ueg one
        ! below..
        if (ic /= 2) then
            pgen = 0.0_dp
            return
        end if

        if (test_flag(ilutI, get_initiator_flag(part_type))) then

            ! can i use the already provided routine also for hubbard models?
            ! yes!
            call CalcPGenLattice(ex, pgen)
        else
            ! in the hubbard case it can also be that we pick the electrons
            ! with the old back-spawn method
            src = get_src(ex)

            if (t_back_spawn) then

                call pick_virtual_electrons_double_hubbard(nI, part_type, d_elecs, d_src, &
                                                           d_ispn, pgen_elec, .true.)

                loc = -1
            else

                ! otherwise:
                pgen_elec = 1.0_dp / real(nOccBeta * nOccAlpha, dp)

                loc = check_electron_location(src, 2, part_type)
            end if

            ! otherwise i have to calculate stuff
            if ((loc == 2) .or. (loc == 1 .and. occ_virt_level /= -1) .or. &
                (loc == 0 .and. occ_virt_level >= 1)) then
                ! i only do that in the back-spawn flex case..

                call pick_occupied_orbital_hubbard(ilutI, part_type, paij, &
                                                   d_orb, .true.)

                tgt = get_tgt(ex)

                ! ok i just realized.. the excitation matrix gets always
                ! sorted.. thats not good for my implementation..
                ! atleast that means that i cant assume that the second
                ! picked orbital is always at position tgt(2)
                ! is it enough to check if both are in the reference?
                ! because when i am here, i know that atleast one is
                ! definetly in the reference.. and if the second also is
                ! i could have picked the orbitals the other way around..
                if (is_in_ref(tgt(1), part_type) .and. is_in_ref(tgt(2), part_type)) then
                    mult = 2.0_dp
                else
                    mult = 1.0_dp
                end if

            else
                ! otherwise we can pick the orbital (a) freely
                paij = 1.0_dp / real(nbasis - nel, dp)
                mult = 2.0_dp

            end if

            pgen = mult * paij * pgen_elec

        end if

    end function calc_pgen_back_spawn_hubbard

    ! to disentangle the mess of excitation generators also implement an new
    ! one for the UEG and the hubbard model seperately
    subroutine gen_excit_back_spawn_ueg(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                        ExcitMat, tParity, pgen, HelGen, store, part_type)
        ! if back-spawn and ueg is turned on point to this excitation
        ! generator! check if we hit all the relevant parts in the code though
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ic, ExcitMat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = "gen_excit_back_spawn_ueg"

        integer :: iUnused
#ifdef DEBUG_
        real(dp) :: pgen2
        HElement_t(dp) :: temp_hel
#endif
        ! so now pack everything necessary for a ueg excitation generator.

        HElGen = 0.0_dp
        iUnused = exFlag
        iUnused = store%nopen

        ic = 2

        ! this function gets pointed to if tUEG, t_back_spawn_flex and tLatticeGens
        ! is set
        ! BUT: again: we also need to take care of the HPHF keyword..
        ! which is a mess

        ! implement it for now without this tNoFailAb flag
        ASSERT(.not. tNoFailAb)

        if (t_back_spawn_flex .and. .not. test_flag(ilutI, get_initiator_flag(part_type))) then

            call gen_double_back_spawn_ueg(nI, ilutI, nJ, ilutJ, tParity, ExcitMat, &
                                           pgen)

#ifdef DEBUG_
            if (.not. IsNullDet(nJ)) then
                pgen2 = calc_pgen_back_spawn_ueg(ilutI, ExcitMat, ic, part_type)
                if (abs(pgen - pgen2) > 1.0e-6_dp) then
                    if (tHPHF) then
                        print *, "due to circular dependence, no matrix element calc possible!"
!                         temp_hel = hphf_off_diag_helement(nI,nJ,ilutI,ilutJ)
                        temp_hel = 0.0_dp
                    else
                        temp_hel = get_helement(nI, nJ, ilutI, ilutJ)
                    end if

                    write(stdout, *) 'Calculated and actual pgens differ. for non-initiator'
                    write(stdout, *) 'This will break HPHF calculations'
                    call write_det(6, nI, .false.)
                    write(stdout, '(" --> ")', advance='no')
                    call write_det(6, nJ, .true.)
                    write(stdout, *) 'Excitation matrix: ', ExcitMat(1, 1:ic), '-->', &
                        ExcitMat(2, 1:ic)
                    write(stdout, *) 'Generated pGen:  ', pgen
                    write(stdout, *) 'Calculated pGen: ', pgen2
                    write(stdout, *) 'matrix element: ', temp_hel
                    call stop_all(this_routine, "Invalid pGen")
                end if
            end if
#endif
        else
            ! do i want to rewrite the old on or just reuse? for now reuse:
            call CreateExcitLattice(nI, ilutI, nJ, tParity, ExcitMat, pgen, part_type)

            if (.not. IsNullDet(nJ)) ilutJ = make_ilutJ(ilutI, ExcitMat, ic)

#ifdef DEBUG_
            if (.not. IsNullDet(nJ)) then
                call CalcPGenLattice(ExcitMat, pgen2)
                if (abs(pgen - pgen2) > 1.0e-6_dp) then
                    if (tHPHF) then
                        print *, "due to circular dependence, no matrix element calc possible!"
!                         temp_hel = hphf_off_diag_helement(nI,nJ,ilutI,ilutJ)
                        temp_hel = 0.0_dp
                    else
                        temp_hel = get_helement(nI, nJ, ilutI, ilutJ)
                    end if

                    write(stdout, *) 'Calculated and actual pgens differ. for non-initiator'
                    write(stdout, *) 'This will break HPHF calculations'
                    call write_det(6, nI, .false.)
                    write(stdout, '(" --> ")', advance='no')
                    call write_det(6, nJ, .true.)
                    write(stdout, *) 'Excitation matrix: ', ExcitMat(1, 1:ic), '-->', &
                        ExcitMat(2, 1:ic)
                    write(stdout, *) 'Generated pGen:  ', pgen
                    write(stdout, *) 'Calculated pGen: ', pgen2
                    write(stdout, *) 'matrix element: ', temp_hel
                    call stop_all(this_routine, "Invalid pGen")
                end if
            end if
#endif

        end if

    end subroutine gen_excit_back_spawn_ueg

    subroutine gen_double_back_spawn_ueg(nI, ilutI, nJ, ilutJ, tParity, ExcitMat, &
                                         pgen, part_type)
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ExcitMat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        integer, optional :: part_type
        character(*), parameter :: this_routine = "gen_double_back_spawn_ueg"

        integer :: elec_i, elec_j, iSpn, orb_b, src(2), &
                   loc, temp_part_type, ki(3), kj(3), ka(3), kb(3), kb_ms, TestEnergyB, &
                   iSpinIndex, orb_a
        real(dp) :: x, pAIJ, dummy, mult
        logical :: tAllowedExcit, t_temp_back_spawn

        t_temp_back_spawn = .false.

        if (present(part_type)) then
            temp_part_type = part_type
        else
            temp_part_type = 1
        end if

        ! i know here that back-spawn is on and this is a non-inititator:
        ! so for the beginning implementation we pick the two electrons
        ! randomly
        ! i should make this to a function below: maybe with the hubbard
        ! flag as additional input to provide us with two independent
        ! electrons in any order possible!
        ! i do not need two do loops in the case of the UEG
        elec_i = 1 + int(genrand_real2_dsfmt() * nel)
        do
            elec_j = 1 + int(genrand_real2_dsfmt() * nel)
            if (elec_j /= elec_i) exit
        end do

        src = nI([elec_i, elec_j])
        iSpn = get_ispn(src)

        ! i need to call nI(elecs)
        loc = check_electron_location(src, 2, temp_part_type)

        ! wait a minute.. thats incorrect or?
        ! if we have both in the occupied manifold we want to restrict
        ! our orbital choice..
        ! i can decide based on the occ_virt_level or?
        if ((loc == 2) .or. (loc == 1 .and. occ_virt_level /= -1) .or. &
            (loc == 0 .and. occ_virt_level >= 1)) then
            t_temp_back_spawn = .true.
            ! i think i need a new one for the ueg also.. since there is not
            ! such a spin restriction as in the hubbard model
            ! i think just using the "normal" occupied picker should be fine
            ! how does this compile even??
            ! no.. write a new one without the spin-restriction
            call pick_occupied_orbital_ueg(ilutI, src, ispn, temp_part_type, pAIJ, &
                                           dummy, orb_a)

            if (orb_a == 0) then
                nJ(1) = 0
                pgen = 0.0_dp
                return
            end if

        else
            ! otherwise pick freely..
            do
                x = genrand_real2_dSFMT()
                if (iSpn == 2) then
                    orb_a = int(nBasis * x) + 1

                    pAIJ = 1.0_dp / real(nBasis - nel, dp)
                else
                    orb_a = 2 * (INT(nBasis / 2 * x) + 1) & ! 2*(a number between 1 and nbasis/2) gives the alpha spin
                    & - (1 - (iSpn / 3)) ! alpha numbered even, iSpn/3 returns 1 for alpha/alpha, 0 for beta/beta
                    if (iSpn == 1) then
                        paij = 1.0_dp / (nbasis / 2 - noccbeta)
                    else !(ispn = 3)..
                        paij = 1.0_dp / (nbasis / 2 - noccalpha)
                    end if
                end if
                if (IsNotOcc(ilutI, orb_a)) exit
            end do

        end if

        ! i guess this is not necessary, but anyway..
        IF ((orb_a < 0) .or. (orb_a > nBasis)) THEN
            CALL Stop_All("CreateDoubExcitLattice", "Incorrect basis function generated")
        end if

        tAllowedExcit = is_allowed_ueg_k_vector(src(1), src(2), orb_a)

        IF (.not. tAllowedExcit) THEN
            nJ(1) = 0
            RETURN
        end if

        orb_b = get_orb_from_kpoints(src(1), src(2), orb_a)

        IF (orb_b == -1 .or. orb_a == orb_b) THEN
            nJ(1) = 0
            pgen = 0.0_dp
            RETURN
        end if

        ! Is b occupied?
        if (.not. tAllowedExcit .or. IsOcc(ilutI, orb_b)) then
            nj(1) = 0
            pgen = 0.0_dp
            return
        end if

        ! Find the new determinant
        call make_double(nI, nJ, elec_i, elec_j, orb_a, &
                         orb_b, ExcitMat, tParity)

        ilutJ = make_ilutJ(ilutI, ExcitMat, 2)

        ! we knoe it is back-spawn + ueg remember! so paij is already above
        !calculate generation probabilities
        ! note, p(b|ij)=p(a|ij) for this system
        ! i have to be careful.. only if b is also in the occupied manifold
        ! if we forced the pick.. then it is * 2
        mult = 2.0_dp
        if (t_temp_back_spawn .and. (.not. is_in_ref(orb_b, part_type))) then
            mult = 1.0_dp
        end if

        pgen = 2.0_dp * mult * paij / real(nel * (nel - 1), dp)

    end subroutine gen_double_back_spawn_ueg

    ! also write a wrapper-like routine for an excitation generator if
    ! back-spawn is activated.. to not mess up all the old functions too much.
    subroutine gen_excit_back_spawn(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                    ExcitMat, tParity, pgen, HelGen, store, part_type)
        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nel), ic, ExcitMat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = "gen_excit_back_spawn"

#ifdef DEBUG_
        HElement_t(dp) :: temp_hel
        real(dp) :: pgen2
#endif
        unused_var(exFlag); unused_var(store)

        HElGen = 0.0_dp
        ! check the non-initiator criteria beforehand
        ! i also have to consider that back-spawn gets turned on later on
        ! so i have to check if back-spawn is active already or not..

        if ((t_back_spawn_flex .or. t_back_spawn) .and. &
            .not. test_flag(ilutI, get_initiator_flag(part_type))) then

            ! otherwise use custom made ones
            if (genrand_real2_dSFMT() < pSingles) then

                ic = 1
                call gen_single_back_spawn(nI, ilutI, part_type, nJ, ilutJ, ExcitMat, &
                                           tParity, pgen)
                pgen = pgen * pSingles

            else

                ic = 2
                call gen_double_back_spawn(nI, ilutI, part_type, nJ, ilutJ, ExcitMat, &
                                           tParity, pgen)
                pgen = pgen * pDoubles

            end if

#ifdef DEBUG_
            if (.not. IsNullDet(nJ)) then
                pgen2 = calc_pgen_back_spawn(nI, ilutI, ExcitMat, ic, part_type)
                if (abs(pgen - pgen2) > 1.0e-6_dp) then
                    if (tHPHF) then
                        print *, "due to circular dependence, no matrix element calc possible!"
!                         temp_hel = hphf_off_diag_helement(nI,nJ,ilutI,ilutJ)
                        temp_hel = 0.0_dp
                    else
                        temp_hel = get_helement(nI, nJ, ilutI, ilutJ)
                    end if

                    write(stdout, *) 'Calculated and actual pgens differ. for non-initiator'
                    write(stdout, *) 'This will break HPHF calculations'
                    write(stdout, *) "reference determinant: "
                    call write_det(6, projedet(:, part_type_to_run(part_type)), .true.)
                    call write_det(6, nI, .false.)
                    write(stdout, '(" --> ")', advance='no')
                    call write_det(6, nJ, .true.)
                    write(stdout, *) 'Excitation matrix: ', ExcitMat(1, 1:ic), '-->', &
                        ExcitMat(2, 1:ic)
                    write(stdout, *) 'Generated pGen:  ', pgen
                    write(stdout, *) 'Calculated pGen: ', pgen2
                    write(stdout, *) 'matrix element: ', temp_hel
                    call stop_all(this_routine, "Invalid pGen")
                end if
            end if
#endif
        else

            ! do the "normal" excitation type if it is an initiator
            if (genrand_real2_dSFMT() < pSingles) then

                ic = 1
                call gen_single_4ind_ex(nI, ilutI, nJ, ilutJ, ExcitMat, &
                                        tParity, pGen)
                pgen = pgen * pSingles

            else

                ic = 2
                call gen_double_4ind_ex2(nI, ilutI, nJ, ilutJ, ExcitMat, &
                                         tParity, pGen)
                pgen = pgen * pDoubles

            end if
#ifdef DEBUG_
            if (.not. IsNullDet(nJ)) then
                pgen2 = calc_pgen_4ind_weighted2(nI, ilutI, ExcitMat, ic)
                if (abs(pgen - pgen2) > 1.0e-6_dp) then
                    if (tHPHF) then
                        print *, "due to circular dependence, no matrix element calc possible!"
!                         temp_hel = hphf_off_diag_helement(nI,nJ,ilutI,ilutJ)
                        temp_hel = 0.0_dp
                    else
                        temp_hel = get_helement(nI, nJ, ilutI, ilutJ)
                    end if

                    write(stdout, *) 'Calculated and actual pgens differ. initiator'
                    write(stdout, *) 'This will break HPHF calculations'
                    call write_det(6, nI, .false.)
                    write(stdout, '(" --> ")', advance='no')
                    call write_det(6, nJ, .true.)
                    write(stdout, *) 'Excitation matrix: ', ExcitMat(1, 1:ic), '-->', &
                        ExcitMat(2, 1:ic)
                    write(stdout, *) 'Generated pGen:  ', pgen
                    write(stdout, *) 'Calculated pGen: ', pgen2
                    write(stdout, *) 'matrix element: ', temp_hel
                    call stop_all(this_routine, "Invalid pGen")
                end if
            end if
#endif
        end if
    end subroutine gen_excit_back_spawn

    subroutine gen_single_back_spawn(nI, ilutI, part_type, nJ, ilutJ, ex, tPar, pgen)
        ! specialised single excitation routine for the back-spawn method
        ! for the moment i still have to decide, which back-spawn method is
        ! in use..
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(in) :: part_type
        integer, intent(out) :: nJ(nel), ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tPar
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "gen_single_back_spawn"

        integer :: elec, src, cc_index, loc, tgt
        real(dp) :: pgen_elec

        ! depending on the method we pick electrons accordingly
        if (t_back_spawn_flex) then
            elec = 1 + floor(genrand_real2_dSFMT() * nel)

            loc = check_electron_location([nI(elec), 0], 1, part_type)

            pgen_elec = 1.0_dp / real(nel, dp)
        else
            call pick_virtual_electron_single(nI, part_type, elec, pgen_elec)
            loc = -1

        end if

        src = nI(elec)

        cc_index = ClassCountInd(get_spin(src), SpinOrbSymLabel(src), &
                                 G1(src)%Ml)

        ! i have to make the logic easier here at some point..
        ! we just need to do some extensive testing and decide on one
        ! back-spawn method and abandon the rest..
        ! i hope this logic is correct: otherwise check on the bottom
        if ((t_back_spawn_occ_virt) .or. (t_back_spawn_flex .and. ( &
                                          (loc == 2 .and. occ_virt_level /= -1) .or. occ_virt_level == 2))) then

            call pick_occupied_orbital_single(ilutI, src, cc_index, part_type, pgen, tgt)

        else

            tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen)

        end if
!
!         if (t_back_spawn_occ_virt) then
!             call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt)
!
!         else if (t_back_spawn_flex) then
!
!             if (loc == 2 .and. occ_virt_level /= -1) then
!                 call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt)
!
!             else
!                 if (occ_virt_level == 2) then
!                     call pick_occupied_orbital_single(nI, src, cc_index, pgen, tgt)
!                 else
!                     tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen)
!                 end if
!             end if
!         else
!             tgt = select_orb_sing(nI, ilutI, src, cc_index, pgen)
!         end if

        if (tgt == 0) then
            nJ(1) = 0
            pgen = 0.0_dp
            return
        end if

        call make_single(nI, nJ, elec, tgt, ex, tPar)

        ilutJ = make_ilutJ(ilutI, ex, 1)

        pgen = pgen * pgen_elec

    end subroutine gen_single_back_spawn

    subroutine gen_double_back_spawn(nI, ilutI, part_type, nJ, ilutJ, ex, tPar, pgen)
        ! the double excitation routine for the back-spawn method
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: part_type
        integer, intent(out) :: nJ(nel), ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NIfTot)
        logical, intent(out) :: tpar
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "gen_double_back_spawn"

        integer :: elecs(2), sym_product, ispn, sum_ml, src(2), loc, cc_a, cc_b, &
                   orbs(2)
        real(dp) :: int_cpt(2), cum_sum(2), sum_pair(2), cpt_pair(2), &
                    cum_arr(nbasis)

        if (t_back_spawn_flex) then
            ! Pick the electrons in a weighted fashion
            call pick_weighted_elecs(nI, elecs, src, sym_product, ispn, sum_ml, &
                                     pgen)

            loc = check_electron_location(src, 2, part_type)

        else
            call pick_virtual_electrons_double(nI, part_type, elecs, src, ispn, &
                                               sum_ml, pgen)

            ! here it could be that there are no valid electrons..
            if (elecs(1) == 0) then
                nJ(1) = 0
                pgen = 0.0_dp
                return
            end if

            ! TODO!! thats the problem here! there are circular dependencies
            ! when I call that from here.. so i guess i need to setup the
            ! back_spawn excitation routines in a seperate file
            sym_product = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
                                                SpinOrbSymLabel(src(2)))

            ! indicate no flex option is used
            loc = -1

        end if

        ! at some point i have to make this whole logic easier here..
        ! in the end i think i have to remove some possibilities and stick to
        ! one back-spawn method.
        if ((t_back_spawn_occ_virt) .or. (t_back_spawn_flex .and. ( &
                                          (loc == 1 .and. occ_virt_level /= -1) .or. (loc == 2) .or. &
                                          (loc == 0 .and. occ_virt_level >= 1)))) then

            call pick_occupied_orbital(ilutI, src, ispn, part_type, int_cpt(1), cum_sum(1), &
                                       orbs(1))

        else

            orbs(1) = pick_a_orb(ilutI, src, iSpn, int_cpt(1), cum_sum(1), cum_arr)

        end if

        if (orbs(1) /= 0) then
            cc_a = ClasSCountInd(orbs(1))
            cc_b = get_paired_cc_ind(cc_a, sym_product, sum_ml, iSpn)

            if (t_back_spawn_flex .and. ((loc == 2 .and. occ_virt_level /= -1) .or. &
                                         (occ_virt_level == 2))) then

                call pick_second_occupied_orbital(ilutI, cc_b, orbs(1), ispn, &
                                                  part_type, int_cpt(2), cum_sum(2), orbs(2))

            else

                orbs(2) = select_orb(ilutI, src, cc_b, orbs(1), int_cpt(2), &
                                     cum_sum(2))
            end if

            ASSERT((.not. (is_beta(orbs(2)) .and. .not. is_beta(orbs(1)))))

        end if

        if (any(orbs == 0)) then
            nJ(1) = 0
            pgen = 0.0_dp
            return
        end if

        ! can i exit right away if this happens??
        ! i am pretty sure this means
        if (any(cum_sum < EPS)) then
            cum_sum = 1.0_dp
            int_cpt = 0.0_dp
            if (.not. same_spin(orbs(1), orbs(2))) then
                nJ(1) = 0
                pgen = 0.0_dp
                return
            end if
        end if

        ! only on parallel excitations.. and symmetric exciation generator is
        ! turned off for now in the back-spawning
        if (same_spin(orbs(1), orbs(2))) then
            if (t_back_spawn_occ_virt .or. (t_back_spawn_flex .and. ( &
                                            (loc == 1 .and. (occ_virt_level == 0 .or. occ_virt_level == 1)) &
                                            .or. (loc == 2 .and. occ_virt_level == -1) .or. &
                                            (loc == 0 .and. occ_virt_level == 1)))) then

                if (is_in_ref(orbs(2), part_type)) then
                    ! if (b) is also in the occupied manifold i could have
                    ! picked the other way around..
                    ! with the same uniform probability:
                    cpt_pair(1) = int_cpt(1)
                    sum_pair(1) = cum_sum(1)
                    ! and then (a) would have been picked according to the
                    ! "normal" procedure
                    call pgen_select_orb(ilutI, src, orbs(2), orbs(1), &
                                         cpt_pair(2), sum_pair(2))
                else
                    ! if (b) is not in the occupied this does not work or?
                    ! since i am forcing (a) to be in the occupied..
                    ! so remove this pgen:
                    cpt_pair = 0.0_dp
                    sum_pair = 1.0_dp
                end if

            else if (t_back_spawn_flex .and. ((loc == 0 .and. occ_virt_level == 2) &
                                              .or. (loc == 1 .and. occ_virt_level == 2) .or. &
                                              (loc == 2 .and. occ_virt_level /= -1))) then

                ! we have picked both in the occupied manifold
                cpt_pair = int_cpt
                sum_pair = cum_sum

            else

                ! otherwise "normal"
                call pgen_select_a_orb(ilutI, src, orbs(2), iSpn, cpt_pair(1), &
                                       sum_pair(1), cum_arr, .false.)
                call pgen_select_orb(ilutI, src, orbs(2), orbs(1), &
                                     cpt_pair(2), sum_pair(2))

            end if

        else

            cpt_pair = 0.0_dp
            sum_Pair = 1.0_dp

        end if

        if (any(sum_pair < EPS)) then
            cpt_pair = 0.0_dp
            sum_pair = 1.0_dp
            ! can something like that slip through:
            if (any(cum_sum < EPS)) then
                pgen = 0.0_dp
                nJ(1) = 0
                return
            end if
        end if

        pgen = pgen * (product(int_cpt) / product(cum_sum) + &
                       product(cpt_pair) / product(sum_pair))

        ! And generate the actual excitation
        call make_double(nI, nJ, elecs(1), elecs(2), orbs(1), orbs(2), &
                         ex, tpar)

        ilutJ = make_ilutJ(ilutI, ex, 2)

    end subroutine gen_double_back_spawn

    function calc_pgen_back_spawn_ueg(ilutI, ex, ic, part_type) result(pgen)
        integer, intent(in) :: ex(2, 2), ic, part_type
        integer(n_int), intent(in) :: ilutI(0:niftot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = "calc_pgen_back_spawn_ueg"

        real(dp) :: cum_sum, pAIJ
        integer :: dummy_orb, ispn, loc, src(2), tgt(2)

        ! i have to write a generation probability calculator for the hphf
        ! implementation this should only be called if it is a non-init
        ! and back-spawn is on.. otherwise just use the already provided
        ! ones
        if (ic /= 2) then
            pgen = 0.0_dp
            return
        end if

        if (test_flag(ilutI, get_initiator_flag(part_type))) then

            ! i have to do some symmetry setup beforehand..
            ! or i do it by hand to avoid the unnecessary overhead..
            ! nope.. i actually only need:
            ! although i should not land here i guess..
            ! this functionality i could actually unit-test.. damn..
            call CalcPGenLattice(ex, pgen)
        else
            ! do the back-spawn pgen..
            ! elec were picked randomly(but in both orders!)
            ! and then we have to check the electron location, to know if
            ! there was a restriction on the first orbitals
            ! the electrons are in the first column or?
            src = get_src(ex)
            loc = check_electron_location(src, 2, part_type)

            ! and with the implementation right now it is only restricted
            ! if both were in the occup
            ! we extend this also in the ueg for occ_virt_level
            if ((loc == 2) .or. (loc == 1 .and. occ_virt_level /= -1) .or. &
                (loc == 0 .and. occ_virt_level >= 1)) then

                ! i could just call  the orbital picking routine..
                ! although this is inefficient..
                ! and i have to see if it would have been possible the
                ! other way around too..
                ! i need ispn here..
                ispn = get_ispn(src)
                call pick_occupied_orbital_ueg(ilutI, src, ispn, part_type, pAIJ, cum_sum, &
                                               dummy_orb, .true.)

                ! i think i can't just use 2*p(a|ij) since this assumption
                ! comes from p(b|ij) = p(a|ij) which is not true by
                ! default if we exclude a to the occupied manifold..
                ! i actually have to check if b is also in the
                ! occupied
                ! only mult by 2 if b is also in the occupied manifolf and
                ! thus it would have been possible to pick it the other
                ! way around
                pgen = 2.0_dp * pAIJ / real(nel * (nel - 1), dp)
                ! i should check the hole location too..
                tgt = get_tgt(ex)
                ! is tgt(1) always the first picked? yes it is!
                if (is_in_ref(tgt(1), part_type) .and. is_in_ref(tgt(2), part_type)) then
                    ! in this case we can recalc p(b|ij)
                    pgen = 2.0_dp * pgen
                end if

            else
                ! otherwise the orbital (a) is free
                ! so it is the number of available orbitals, depending on
                ! the spin
                ! here i can just use the formulas for the standard ueg..
                ! or just use the above provided routine
                call CalcPGenLattice(ex, pgen)
            end if
        end if

    end function calc_pgen_back_spawn_ueg

    function calc_pgen_back_spawn(nI, ilutI, ex, ic, part_type) result(pgen)
        ! to use HPHF keyword and also the test if the pgens are correct
        ! or just to be able and to be sure and more save i need a way to
        ! recalculate the generation probability also for a back-spawn
        ! excitation from a non-initiator determinant
        ! so although thats a hassle implement it, otherwise i cannot be
        ! quite sure about the method
        integer, intent(in) :: nI(nel), ex(2, 2), ic, part_type
        integer(n_int), intent(in) :: ilutI(0:niftot)
        real(dp) :: pgen
        character(*), parameter :: this_routine = "calc_pgen_back_spawn"

        integer :: dummy, ssrc, stgt, cc_index, src(2), tgt(2), dummy_elecs(2), &
                   dummy_orbs(2), ispn, loc, sum_ml, sym_prod, cc_a, cc_b, &
                   dummy_ispn, dummy_sum_ml
        real(dp) :: elec_pgen, orb_pgen, int_cpt(2), cum_sum(2), cpt_pair(2), &
                    sum_pair(2), cum_arr(nbasis)
        logical :: t_gen_list, t_in_ref, t_par
        ! i should only call this function when i am on a non-inititator or?
        ! in the HPHF framework i could end up here also on a non-initiator
        ! so here i should then go to 4ind-2 pgen calculator if it is a
        ! inititator

        if (test_flag(ilutI, get_initiator_flag(part_type))) then

            pgen = calc_pgen_4ind_weighted2(nI, ilutI, ex, ic)

        else
            ! in the hphf framework we definetly only end up here if the
            ! back-spawn method is already turned on.. but is this
            ! ensured in the other places the function might get called?
            ! decision: if we call this function we want to get the pgen in
            ! the back-spawn method, so ensure outside if it is turned on
            ! or not.

            if (ic == 1) then

                ! depending on the implementation
                ! here since i want to calculate the real pgen if back_spawn
                ! is or would be on, i should work with the
                ! _option keywords..
                ssrc = ex(1, 1)
                stgt = ex(2, 1)

                if (t_back_spawn_flex_option) then
                    elec_pgen = 1.0_dp / real(nel, dp)

                    loc = check_electron_location([ssrc, 0], 1, part_type)

                else
                    ! it is slow anyway.. so just do the dummy implementation
                    ! as Ali mentioned in the HPHF implementation it is not
                    ! that common to have a doubly connected HPHF det

                    call pick_virtual_electron_single(nI, part_type, dummy, elec_pgen, .true.)

                    loc = -1

                end if

                if ((t_back_spawn_occ_virt) .or. (t_back_spawn_flex .and. &
                                                  ((loc == 2 .and. occ_virt_level /= -1) .or. occ_virt_level == 2))) then

                    ! the symmetry of the orbital is known after all
                    cc_index = ClassCountInd(get_spin(stgt), SpinOrbSymLabel(stgt), &
                                             G1(stgt)%Ml)

                    ! reuse the routine..
                    call pick_occupied_orbital_single(ilutI, ssrc, cc_index, &
                                                      part_type, orb_pgen, dummy, .true.)

                else

                    ! reuse the other pgen single function
                    ! but i have to multiply out the electron pgen
                    orb_pgen = pgen_single_4ind(nI, ilutI, ssrc, stgt) * &
                               real(nel, dp)

                end if

                pgen = pSingles * elec_pgen * orb_pgen

            else if (ic == 2) then

                src = get_src(ex)

                ispn = get_ispn(src)

                sum_ml = sum(G1(src)%ml)

                sym_prod = RandExcitSymLabelProd(SpinOrbSymLabel(src(1)), &
                                                 SpinOrbSymLabel(src(2)))

                ! i have to be careful with the orbitals..
                ! because i guess those get ordered.. and i have to check if
                ! it is possible to have picked the orbitals in either
                ! order..
                ! ok.. and there is the restriction to pick a beta orbital
                ! in src(1) first if it is a anti-parallel excitation
                ! ok this means atleast src(1) is definetly the first picked
                ! orbital.. is this also the case in HPHF?? argh i hate that
                ! stuff
                ! do this testing here once
                ! todo: i have to fix this since here i do not know anymore
                ! what was orb a and orb b since ex is sorted..
                tgt = get_tgt(ex)
                t_in_ref = (is_in_ref(tgt(1), part_type) .and. is_in_ref(tgt(2), part_type))
                t_par = (is_beta(tgt(1)) .eqv. is_beta(tgt(2)))

                ! now i know.. for opposite spin excitations the beta orbital
                ! of tgt was the first picked one! so it would be best to
                ! switch it back, so i can correctly recalculate the
                ! pgen.. for parallel spin excitations the order does not
                ! matter
                if (.not. t_par) then
                    if (.not. is_beta(tgt(1))) then
                        tgt = [tgt(2), tgt(1)]
                    end if
                end if

                if (t_back_spawn_flex) then

                    elec_pgen = pgen_weighted_elecs(nI, src)

                    loc = check_electron_location(src, 2, part_type)

                else

                    call pick_virtual_electrons_double(nI, part_type, dummy_elecs, &
                                                       dummy_orbs, dummy_ispn, dummy_sum_ml, elec_pgen, .true.)

                    loc = -1

                end if
                ! for some back_spawn_flex it can happen that we have no
                ! restrictions on the orbitals.. but thats rare i guess..
                if (t_back_spawn_occ_virt .or. (t_back_spawn_flex .and. &
                                                ((loc == 1 .and. occ_virt_level /= -1) .or. loc == 2 .or. &
                                                 (loc == 0 .and. occ_virt_level >= 1)))) then

                    ! in this case it matters which orbital was picked
                    ! first..
                    ! in this case we can be sure that atleast on of the
                    ! orbitals is in the reference..
                    if (t_par) then
                        if (.not. is_in_ref(tgt(1), part_type)) then
                            ! then it is incorrectly ordered.. reorder!
                            tgt = [tgt(2), tgt(1)]
                        end if
                    end if

                    call pick_occupied_orbital(ilutI, src, ispn, &
                                               part_type, int_cpt(1), cum_sum(1), dummy_orbs(1), .true.)

                    ! i can atleast do some stuff for picking it the other
                    ! way or?
                    ! because if it is a parallel spin-excitations and orbital
                    ! (b) is in the reference the probs p(b|ij) are the same
                    if (t_par .and. t_in_ref) then
                        cpt_pair(1) = int_cpt(1)
                        sum_pair(1) = cum_sum(1)

                    else
                        cpt_pair = 0.0_dp
                        sum_pair = 1.0_dp
                        ! maybe i could set a flag that it does not have to
                        ! be recalced otherwise..
                    end if

                    ! here i should go on to calculate p(b|aij) since in the
                    ! other case we have everything..

                    if (t_back_spawn_flex .and. ( &
                        (loc == 2 .and. occ_virt_level /= -1) .or. occ_virt_level == 2)) then

                        ! this is the case where orbital (b) is also restricted

                        cc_a = ClassCountInd(tgt(1))
                        cc_b = get_paired_cc_ind(cc_a, sym_prod, sum_ml, ispn)
                        call pick_second_occupied_orbital(ilutI, cc_b, tgt(1), &
                                                          ispn, part_type, int_cpt(2), cum_sum(2), dummy_orbs(2), .true.)

                        ! and ofc both probs are the same if the spin-fits
                        if (t_par) then
                            cpt_pair = int_cpt
                            sum_pair = cum_sum
                        else
                            cpt_pair = 0.0_dp
                            sum_pair = 1.0_dp
                        end if

                    else
                        ! the order should not matter or??
                        ! or does it? and i always force a beta orbital
                        ! to be picked first.. hm..
                        ! the order does matter! and the excitation matrix
                        ! is always ordered at this point but during picking
                        ! the orbitals it is not!
                        ! so we have to be more careful here!
                        call pgen_select_orb(ilutI, src, tgt(1), tgt(2), int_cpt(2), &
                                             cum_sum(2))

                        ! in this case (a) was restricted but (b) was not
                        ! so check if the other way would have been
                        ! possible
                        if (t_par .and. t_in_ref) then
                            ! p(b|ij) already set above
                            call pgen_select_orb(ilutI, src, tgt(2), tgt(1), &
                                                 cpt_pair(2), sum_pair(2))

                        else
                            cpt_pair = 0.0_dp
                            sum_pair = 1.0_dp
                        end if

                    end if

                else

                    call pgen_select_a_orb(ilutI, src, tgt(1), ispn, int_cpt(1), &
                                           cum_sum(1), cum_arr, .true.)

                    ! but i guess in this case i can already cover all the
                    ! pgen..
                    if (int_cpt(1) > EPS) then
                        call pgen_select_orb(ilutI, src, tgt(1), tgt(2), &
                                             int_cpt(2), cum_sum(2))

                        t_gen_list = .false.
                    else
                        t_gen_list = .false.
                        int_cpt = 0.0_dp
                        cum_sum = 1.0_dp
                    end if
                    ! otherwise i can pick orbital (a) freely.. which also
                    ! means that i definetly can pick (b) freely and do it
                    ! the other way around..

                    call pgen_select_a_orb(ilutI, src, tgt(2), ispn, cpt_pair(1), &
                                           sum_pair(1), cum_arr, t_gen_list)

                    if (cpt_pair(1) > EPS) then
                        call pgen_select_orb(ilutI, src, tgt(2), tgt(1), &
                                             cpt_pair(2), sum_pair(2))
                    else
                        cpt_pair = 0.0_dp
                        sum_pair = 1.0_dp
                    end if
                end if

                ! how do we handle this now??
                ! i only want to handle these cases, which i have not
                ! recalculated above already. especially for the p(b|ij)
                ! probability..
                ! duh.. i have to actually do something with the above
                ! calculated pgens..

                ! now i have to figure p(ab), p(ba) stuff.. and implement this
                ! above.. annoying..

                if (any(cum_sum < EPS)) then
                    int_cpt = 0.0_dp
                    cum_sum = 1.0_dp
                end if
                if (any(sum_pair < EPS)) then
                    cpt_pair = 0.0_dp
                    sum_pair = 1.0_dp
                end if

                pgen = pDoubles * elec_pgen * (product(int_cpt) / product(cum_sum) + &
                                               product(cpt_pair) / product(sum_pair))

            else

                pgen = 0.0_dp

            end if
        end if

    end function calc_pgen_back_spawn
end module back_spawn_excit_gen
