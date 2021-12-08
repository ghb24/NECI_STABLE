#include "macros.h"
! GUGA excitations module:
! contains as much excitation related functionality as possible
module guga_excitations
    ! modules
    use CalcData, only: t_trunc_guga_pgen, t_trunc_guga_matel, &
                        trunc_guga_pgen, trunc_guga_matel, t_trunc_guga_pgen_noninits, &
                        t_matele_cutoff, matele_cutoff

    use SystemData, only: nEl, nBasis, ElecPairs, G1, nSpatOrbs, &
                          tGen_guga_weighted, nBasisMax, treal, &
                          t_approx_exchange, t_approx_exchange_noninits, &
                          is_init_guga, t_heisenberg_model, t_tJ_model, t_mixed_hubbard, &
                          t_guga_pchb, thub

    use constants, only: dp, n_int, bits_n_int, Root2, &
                         EPS, bni_, bn2_, stdout, int_rdm

    use bit_reps, only: decode_bit_det

    use bit_rep_data, only: GugaBits, niftot, nifguga, nifd

    use DetBitOps, only: count_open_orbs, return_ms

    use guga_data, only: ExcitationInformation_t, getSingleMatrixElement, &
                         getDoubleMatrixElement, getMixedFullStop, &
                         WeightData_t, orbitalIndex, &
                         tag_excitations, tag_tmp_excits, &
                         excit_type, gen_type, excit_names

    use guga_bitRepOps, only: isProperCSF_ilut, calcB_vector_ilut, getDeltaB, &
                              setDeltaB, count_open_orbs_ij, &
                              encode_matrix_element, update_matrix_element, &
                              extract_matrix_element, write_det_guga, &
                              write_guga_list, add_guga_lists, count_alpha_orbs_ij, &
                              count_beta_orbs_ij, findFirstSwitch, findLastSwitch, &
                              calcStepvector, find_switches, &
                              calcOcc_vector_int, &
                              extract_h_element, encode_stochastic_rdm_info, &
                              get_preceeding_opposites, &
                              CSF_Info_t, csf_ref

    use guga_matrixElements, only: calc_guga_matrix_element, calc_mixed_contr_integral, &
                calcremainingswitches_excitinfo_double, calcremainingswitches_excitinfo_single, &
                calcstartprob, calcstayingprob, endfx, endgx, &
                init_fullstartweight, init_singleweight

    use OneEInts, only: GetTMatEl

    use procedure_pointers, only: get_umat_el

    use dSFMT_interface, only: genrand_real2_dSFMT

    use FciMCData, only: ilutRef, tFillingStochRDMOnFly

    use util_mod, only: binary_search_first_ge, abs_l1, operator(.isclose.), &
                        operator(.div.), near_zero, stop_all, neci_flush

    use GenRandSymExcitNUMod, only: RandExcitSymLabelProd

    use SymExcitDataMod, only: SpinOrbSymLabel, OrbClassCount, SymLabelCounts2, &
                               KPointToBasisFn, sym_label_list_spat

    use sym_general_mod, only: ClassCountInd

    use excit_gens_int_weighted, only: get_paired_cc_ind

    use umatcache, only: gtID

    use guga_procedure_pointers, only: &
                calc_orbital_pgen_contr, &
                calc_orbital_pgen_contrib_start, calc_orbital_pgen_contrib_end

    use Umatcache, only: UMat2D

    use timing_neci, only: timer, set_timer, halt_timer, get_total_time

    use guga_types, only: WeightObj_t

    use MemoryManager, only: LogMemAlloc, LogMemDealloc

    use sym_mod, only: MomPbcSym

    use guga_bitRepOps, only: contract_2_rdm_ind

    use lattice_models_utils, only: create_all_open_shell_dets



    use guga_matrixElements, only: &
        calc_mixed_start_contr_sym, calc_mixed_end_contr_sym, &
        calc_mixed_contr_sym



    use guga_matrixElements, only: &
        get_forced_zero_double, getminus_double, getminus_semistart, &
        getplus_double, getplus_semistart, init_doubleweight, &
        init_semistartweight

    better_implicit_none

    private
    public :: global_excitinfo, print_excitinfo, &
              assign_excitinfo_values_double, assign_excitinfo_values_single, &
              actHamiltonian, calcdoubleexcitation_withweight, &
              calcnonoverlapdouble, calcsingleoverlaplowering, calcsingleoverlapraising, &
              calcsingleoverlapmixed, calcdoublelowering, calcdoubleraising, &
              calcdoublel2r, calcdoubler2l, calcfullstoplowering, calcfullstopraising, &
              calcfullstopl2r, calcfullstopr2l, calcfullstartlowering, &
              calcfullstartraising, calcfullstartl2r, calcfullstartr2l, &
              calcfullstartfullstopalike, calcfullstartfullstopmixed, &
              calcremainingswitches_excitinfo_double, checkcompatibility, &
              createsinglestart, singleupdate, singleend, init_singleweight, &
              calcremainingswitches_excitinfo_single, excitationIdentifier, &
              pickorbs_sym_uniform_ueg_single, pickorbs_sym_uniform_ueg_double, &
              pickorbs_sym_uniform_mol_single, pickorbs_sym_uniform_mol_double, &
              calc_orbital_pgen_contr_ueg, calc_orbital_pgen_contr_mol, &
              calc_mixed_x2x_ueg, &
              calc_mixed_end_contr_sym, &
              calc_off_diag_guga_ref_direct, &
              pickorbs_real_hubbard_single, pickorbs_real_hubbard_double, &
              excitationIdentifier_single, excitationIdentifier_double, &
              init_doubleWeight, init_semiStartWeight, init_fullStartWeight, &
              calcremainingswitches_single, calcallexcitations_single, &
              calcallexcitations_double, &
              createstochasticstart_single, &
              pickrandomorb_scalar, pickrandomorb_forced, pickrandomorb_vector, &
              pickrandomorb_restricted, singlestochasticupdate, &
              singlestochasticend, &
              calcfullstartfullstopmixedstochastic, mixedfullstartstochastic, &
              doubleupdatestochastic, calcfullstartraisingstochastic, &
              calcfullstartloweringstochastic, calcfullstopraisingstochastic, &
              calcfullstoploweringstochastic, calcsingleoverlapmixedstochastic, &
              mixedfullstopstochastic, calcloweringsemistartstochastic, &
              calcraisingsemistartstochastic, calcloweringsemistopstochastic, &
              calcraisingsemistopstochastic, calcfullstartl2r_stochastic, &
              calcfullstartr2l_stochastic, calcfullstopl2r_stochastic, &
              calcfullstopr2l_stochastic, calcdoubleloweringstochastic, &
              calcdoubleraisingstochastic, calcdoublel2r2l_stochastic, &
              calcdoubler2l2r_stochastic, calcdoublel2r_stochastic, &
              calcdoubler2l_stochastic, calcallexcitations, &
              pick_elec_pair_uniform_guga, get_guga_integral_contrib, &
              calc_pgen_mol_guga_single, get_excit_level_from_excitInfo, &
              get_guga_integral_contrib_spat, calc_orbital_pgen_contrib_start_def, &
              calc_orbital_pgen_contrib_end_def, create_hamiltonian_guga, &
              csf_to_sds_ilut, csf_vector_to_sds, checkCompatibility_single

    ! use a global excitationInformation type variable to store information
    ! about the last generated excitation to analyze matrix elements and
    ! pgens
    type(ExcitationInformation_t) :: global_excitInfo

    abstract interface
        function calc_pgen_general(csf_i, i) result(pgen)
            import :: dp, CSF_Info_t
            type(CSF_Info_t), intent(in) :: csf_i
            integer, intent(in) :: i
            real(dp) :: pgen
        end function calc_pgen_general
    end interface

    interface excitationIdentifier
        module procedure excitationIdentifier_single
        module procedure excitationIdentifier_double
    end interface excitationIdentifier

    interface calcAllExcitations
        module procedure calcAllExcitations_single
        module procedure calcAllExcitations_double
        module procedure calcAllExcitations_excitInfo_single
    end interface calcAllExcitations

contains

    subroutine csf_vector_to_sds(csfs, csf_coeffs, sds, sd_coeffs, ms)
        integer(n_int), intent(in) :: csfs(:,:)
        real(dp), intent(in) :: csf_coeffs(:)
        real(dp), intent(in), optional :: ms
        integer(n_int), intent(out), allocatable :: sds(:,:)
        real(dp), intent(out), allocatable :: sd_coeffs(:)

        real(dp) :: ms_
        integer :: n_sds, spin, n_tot, i
        integer(n_int), allocatable :: temp_all(:,:), temp_sds(:,:)
        real(dp), allocatable :: temp_coeffs(:)


        spin = abs(return_ms(csfs(:,1)))
        def_default(ms_, ms, spin/2.)

        n_sds = 2 ** nSpatorbs
        allocate(temp_all(0:GugaBits%len_tot,n_sds), source = 0_n_int)

        n_tot = 0

        do i = 1, size(csfs,2)
            call csf_to_sds_ilut(csfs(:,i), temp_sds, temp_coeffs, ms_, csf_coeffs(i))
            call add_guga_lists(n_tot, size(temp_sds,2), temp_all, temp_sds)
        end do

        allocate(sds(0:GugaBits%len_tot, n_tot), source = temp_all(:,1:n_tot))
        allocate(sd_coeffs(n_tot), source = 0.0_dp)

        do i = 1, n_tot
            sd_coeffs(i) = extract_matrix_element(sds(:,i), 1)
        end do

    end subroutine csf_vector_to_sds

    subroutine csf_to_sds_ilut(csf, sds, weights, ms, coeff)
        integer(n_int), intent(in) :: csf(0:GugaBits%len_tot)
        real(dp), intent(in), optional :: ms
        real(dp), intent(in), optional :: coeff
        integer(n_int), intent(out), allocatable :: sds(:,:)
        real(dp), intent(out), allocatable :: weights(:)
        character(*), parameter :: this_routine = "csf_to_sds_ilut"

        real(dp) :: ms_
        integer :: spin, n_alpha, n_beta
        integer :: i, j, step(nSpatorbs), delta_k
        integer :: nI(nel)
        real(dp) :: bVec(nSpatorbs), aVec(nSpatorbs), lambda_k, x, coeff_
        integer(n_int), allocatable :: all_sds(:,:)
        real(dp), allocatable :: all_weights(:)

        def_default(coeff_, coeff, 1.0_dp)


        ! only works for heisenberg model for now..
        if (any(calcOcc_vector_int(csf) /= 1)) then
            call stop_all(this_routine, "only implemented for heisenberg for now")
        end if
        if (near_zero(coeff_)) then
            allocate(sds(0:GugaBits%len_tot, 0))
            allocate(weights(0))
            return
        end if
        spin = abs(return_ms(csf))

        def_default(ms_, ms, spin/2.)

        ! construct SDs by attaching (N/2 + ms) up spins and (N/2 - ms) down spins
        n_alpha = nint(nSpatorbs / 2. + ms_)
        n_beta = nint(nSpatorbs / 2. - ms_)

        all_sds = create_all_open_shell_dets(nSpatorbs, n_beta, n_alpha)

        allocate(all_weights(size(all_sds,2)), source = 0.0_dp)

        step = calcStepvector(csf(0:GugaBits%len_orb))
        bVec = calcB_vector_ilut(csf(0:GugaBits%len_orb))
        aVec = real(([(i, i = 1,nSpatorbs)] - bVec),dp) / 2.0

        do i = 1, size(all_sds,2)
            x = 1.0_dp

            call decode_bit_det(nI, all_sds(:,i))

            do j = 1, nSpatorbs
                if (.not. near_zero(x)) then
                    delta_k = mod(get_spin(nI(j)),2)
                    lambda_k = get_preceeding_opposites(nI, j)

                    select case (step(j))

                    case (1)
                        x = x * sqrt((aVec(j) + bVec(j) - lambda_k)/bVec(j))

                    case (2)
                        x = x * (-1.0_dp)**(bVec(j)+delta_k) * &
                            sqrt((lambda_k - aVec(j) + 1.0_dp)/(bVec(j)+2.0_dp))

                    case (3)

                        x = x * (-1.0_dp)**bVec(j)

                    end select

                end if
            end do

            all_weights(i) =  x

        end do


        j = 1
        do i = 1, size(all_sds,2)
            if (.not. near_zero(all_weights(i))) then
                all_sds(:,j) = all_sds(:,i)
                all_weights(j) = all_weights(i)
                j = j + 1
            end if
        end do

        ! i might have to normalize
        all_weights(1:j-1) = coeff_ * all_weights(1:j-1) / sqrt(sum(all_weights(1:j-1)**2))

        allocate(sds(0:GugaBits%len_tot, j - 1), source = 0_n_int)
        allocate(weights(j-1), source = all_weights(1:j-1))

        do i = 1, j - 1
            sds(0:0, i) = all_sds(:, i)
            call encode_matrix_element(sds(:,i), all_weights(i), 1)
        end do

    end subroutine csf_to_sds_ilut


    function calc_off_diag_guga_ref_direct(ilut, csf_i, run, exlevel) result(hel)
        integer(n_int), intent(in) :: ilut(0:niftot)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in), optional :: run
        integer, intent(out), optional :: exlevel
        HElement_t(dp) :: hel

        integer :: run_
        type(ExcitationInformation_t) :: excitInfo

        def_default(run_, run, 1)
        call calc_guga_matrix_element(ilut, csf_i, ilutRef(0:niftot, run_), csf_ref(run_), excitInfo, hel, .true.)

        if (present(exlevel)) then
            if (excitInfo%valid) then
                exlevel = merge(1, 2, excitInfo%typ == excit_type%single)
            else
                ! non-valid > 3 excit
                exlevel = nel
            end if
        end if
    end function calc_off_diag_guga_ref_direct

    function calc_guga_mat_wrapper(ilutI, csf_i, ilutJ, csf_j) result(mat_ele)
        integer(n_int), intent(in) :: ilutI(0:niftot), ilutJ(0:niftot)
        type(CSF_Info_t), intent(in) :: csf_i, csf_j
        HElement_t(dp) :: mat_ele

        type(ExcitationInformation_t) :: excitInfo

        call calc_guga_matrix_element(ilutI, csf_i, ilutJ, csf_j, excitInfo, mat_ele, t_hamil = .true.)

    end function calc_guga_mat_wrapper

    function create_hamiltonian_guga(ilut_list) result(hamil)
        integer(n_int), intent(in) :: ilut_list(:,:)
        HElement_t(dp) :: hamil(size(ilut_list,2), size(ilut_list,2))

        type(CSF_Info_t) :: csf_i, csf_j
        integer :: i, j

        do i = 1, size(ilut_list,2)
            csf_i = CSF_Info_t(ilut_list(:, i))
            do j = 1, size(ilut_list,2)
                csf_j = CSF_Info_t(ilut_list(:, j))
                hamil(i,j) = calc_guga_mat_wrapper(ilut_list(:, j), csf_j, ilut_list(:,i), csf_i)
            end do
        end do

    end function create_hamiltonian_guga


    ! write up all the specific stochastic excitation routines

    subroutine calcFullStartFullStopMixedStochastic(ilut, csf_i, excitInfo, t, pgen, &
                                                posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcFullStartFullStopMixedStochastic"

        type(WeightObj_t) :: weights
        real(dp) ::  branch_pgen, temp_pgen, above_cpt, below_cpt, rdm_mat, p_orig
        HElement_t(dp) :: integral
        integer :: iOrb, i, j, k, l, typ

        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(.not. isThree(ilut, excitInfo%fullStart))
        ASSERT(.not. isZero(ilut, excitInfo%fullEnd))
        ASSERT(.not. isThree(ilut, excitInfo%fullEnd))

        if (present(opt_weight)) then
            weights = opt_weight
        else
            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. &
                    (.not. is_init_guga))) then
                weights = init_forced_end_exchange_weight(csf_i, excitInfo%fullEnd)
            else
                weights = init_doubleWeight(csf_i, excitInfo%fullEnd)
            end if
        end if

        if (t_approx_exchange .or. (t_approx_exchange_noninits .and. &
                (.not. is_init_guga))) then
            call forced_mixed_start(ilut, csf_i, excitInfo, t, branch_pgen)
        else
            call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)
        end if

        ! set pgen to 0 in case of early exit
        pgen = 0.0_dp

        ! in this excitation there has to be a switch somewhere in the double
        ! overlap region, so only x1 matrix element counts essentially
        ! and this can be 0 at the full-start already -> so check that here
        ! and in case abort
        ! should i do that in the mixedFullStartStochastic routine

        ! do that x1 matrix element in the routine and only check probWeight here
        check_abort_excit(branch_pgen, t)

        temp_pgen = 1.0_dp

        do iOrb = excitInfo%fullStart + 1, excitInfo%fullEnd - 1
            call doubleUpdateStochastic(ilut, csf_i, iOrb, excitInfo, weights, negSwitches, &
                                        posSwitches, t, branch_pgen)

            ! zero x1 - elements can also happen in the double update
            if (near_zero(extract_matrix_element(t, 2)) .or. near_zero(branch_pgen)) then
                t = 0_n_int
                return
            end if
        end do

        call mixedFullStopStochastic(ilut, csf_i, excitInfo, t)

        ! check if there was a change in the stepvector in the double
        ! overlap region
        if (.not. near_zero(extract_matrix_element(t, 1))) then
            t = 0_n_int
            return
        end if

        ! if we do RDMs also store the x0 and x1 coupling coeffs
        ! and I need to do it before the routines below since excitInfo
        ! gets changed there
        if (tFillingStochRDMOnFly) then
            ! i need to unbias against the total pgen later on in the
            ! RDM sampling otherwise the rdm-bias factor is not correct!
            ! encode the necessary information in the rdm-matele!
            i = excitInfo%i
            j = excitInfo%j
            k = excitInfo%k
            l = excitInfo%l
            typ = excitInfo%typ
            rdm_mat = extract_matrix_element(t, 2)
            call calc_orbital_pgen_contr(csf_i, [2 * i, 2 * j], above_cpt, &
                                         below_cpt)
            p_orig = (above_cpt + below_cpt) * branch_pgen
            if (.not. (t_heisenberg_model .or. t_tJ_model)) then
                p_orig = p_orig / real(ElecPairs, dp)
            end if
        end if

        global_excitInfo = excitInfo

        if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then

            if (getDeltaB(t) == 0) then
                t = 0_n_int
                return
            end if

            call calc_mixed_contr_integral(ilut, csf_i, t, excitInfo%fullStart, &
                excitInfo%fullEnd, integral)

            pgen = branch_pgen

            ! just to be save that a switch always happens at the end
            ! print that out for now
        else
            call calc_mixed_contr_sym(ilut, t, csf_i, excitInfo, pgen, integral)
        end if

        if (near_zero(integral)) then
            t = 0_n_int
            pgen = 0.0_dp
            return
        end if

        if (tFillingStochRDMOnFly) then
            if (.not. near_zero(p_orig)) then
                call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                           contract_2_rdm_ind(i, j, k, l, excit_lvl = 2, &
                           excit_typ=typ), x0 = 0.0_dp, x1 = rdm_mat * pgen / p_orig)
            end if
        end if

        call encode_matrix_element(t, 0.0_dp, 2)
        call encode_matrix_element(t, integral, 1)

    end subroutine calcFullStartFullStopMixedStochastic

    subroutine calc_orbital_pgen_contr_mol(csf_i, occ_orbs, cpt_a, cpt_b)
        ! calculates the cumulatice probability list for different
        ! full-start -> full-stop mixed excitations, used in the recalculation of
        ! contrbuting pgens from different picked orbitals
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cpt_a, cpt_b

        integer :: i, j, orb
        real(dp) :: cum_sum, cpt_ab, cpt_ba, ba_sum, ab_sum

        ! given 2 already picked electrons, this routine creates a list of
        ! p(a)*p(b|a) probabilities to pick the already determined holes
        ! is this so easy??
        ! what do i know? i know there is a possible switch between i and j
        ! or otherwise those indices would not have been taken, since an
        ! excitation was already created
        ! i know electrons i and j are from different orbitals and singly
        ! occupied
        ! have to loop once over all orbitals to create list for p(x)
        ! most def, but that should be easy since there are not so many
        ! restrictions, but to determine p(b|a) i have to take into account
        ! the symmetries and the guga restrictions..
        ! but i need the whole information to get cum_sum correct, for the
        ! given electrons i and j..
        ! maybe split up the loop to different sectors, where i know more info

        ! UPDATE: calculate it more directly! also for UEG models!
        ! just output the nececcary p(a)*p(b|a) and vv. and not a whole
        ! list!

        ! chanve that routine now, to use the general pre-generated cum-list
        ! for the current CSF
        i = gtID(minval(occ_orbs))
        j = gtID(maxval(occ_orbs))

        cum_sum = 0.0_dp

        if (tGen_guga_weighted) then
            do orb = 1, i - 1
                ! calc. the p(a)
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)

                end if

            end do
        end if

        ! deal with orb (i) in specific way:
        cpt_a = get_guga_integral_contrib(occ_orbs, i, -1)

        cum_sum = cum_sum + cpt_a

        ! also get p(b|a)
        ! did i get that the wrong way around??
        call pgen_select_orb_guga_mol(csf_i, occ_orbs, i, j, cpt_ba, ba_sum, i, .true.)

        if (tGen_guga_weighted) then
            do orb = i + 1, j - 1
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)
                end if
            end do
        end if

        ! deal with j also speciallly
        cpt_b = get_guga_integral_contrib(occ_orbs, j, -1)

        cum_sum = cum_sum + cpt_b

        ! and get p(a|b)
        call pgen_select_orb_guga_mol(csf_i, occ_orbs, j, i, cpt_ab, ab_sum, -j, .true.)

        ! and deal with rest:

        if (tGen_guga_weighted) then
            do orb = j + 1, nSpatOrbs
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)
                end if
            end do
        end if

        if (.not. tGen_guga_weighted) then
            cum_sum = csf_i%cum_list(nSpatOrbs)
        end if

        if (near_zero(cum_sum) .or. near_zero(ab_sum) .or. near_zero(ba_sum)) then
            cpt_a = 0.0_dp
            cpt_b = 0.0_dp
        else
            ! and get hopefully correct final values:
            cpt_a = cpt_a / cum_sum * cpt_ba / ba_sum
            cpt_b = cpt_b / cum_sum * cpt_ab / ab_sum
        end if

    end subroutine calc_orbital_pgen_contr_mol

    subroutine calc_orbital_pgen_contr_ueg(csf_i, occ_orbs, above_cpt, below_cpt)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: above_cpt, below_cpt

        real(dp) :: cum_sum, cum_arr(nSpatOrbs)
        type(ExcitationInformation_t) :: tmp_excitInfo(nSpatOrbs)
        integer :: tmp_orbArr(nSpatOrbs)
        integer :: i, j

        call gen_ab_cum_list_1_1(csf_i, occ_orbs, cum_arr, tmp_excitInfo, tmp_orbArr)
        cum_sum = cum_arr(nSpatOrbs)

        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))

        ! get orbital probability at beginning: question is can they be
        ! 0, and can they be independently 0, and if yes how is the
        ! influence...
        above_cpt = (cum_arr(j) - cum_arr(j - 1)) / cum_sum

        if (i == 1) then
            below_cpt = cum_arr(1) / cum_sum
        else
            below_cpt = (cum_arr(i) - cum_arr(i - 1)) / cum_sum
        end if

    end subroutine calc_orbital_pgen_contr_ueg

    subroutine calcDoubleR2L_stochastic(ilut, csf_i, excitInfo, t, branch_pgen, &
                                        posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcDoubleR2L_stochastic"

        integer :: iOrb, start2, ende1, ende2, start1, switch
        type(WeightObj_t) :: weights
        real(dp) :: temp_pgen
        HElement_t(dp) :: integral

        ASSERT(.not. isThree(ilut, excitInfo%fullStart))
        ASSERT(.not. isZero(ilut, excitInfo%secondStart))
        ASSERT(.not. isZero(ilut, excitInfo%firstEnd))
        ASSERT(.not. isThree(ilut, excitInfo%fullEnd))

        start1 = excitInfo%fullStart
        start2 = excitInfo%secondStart
        ende1 = excitInfo%firstEnd
        ende2 = excitInfo%fullEnd

        ! : create correct weights:
        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_fullDoubleWeight(csf_i, start2, ende1, ende2, negSwitches(start2), &
                                            negSwitches(ende1), posSwitches(start2), posSwitches(ende1), &
                                            csf_i%B_real(start2), csf_i%B_real(ende1))
        end if

        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = start1 + 1, start2 - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)

            branch_pgen = branch_pgen * temp_pgen
            ! check validity
            check_abort_excit(branch_pgen, t)

        end do

        ! change weights... maybe need both single and double type weights
        ! then do lowering semi start
        weights = weights%ptr

        call calcLoweringSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                             posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = start2 + 1, ende1 - 1
            call doubleUpdateStochastic(ilut, csf_i, iOrb, excitInfo, weights, negSwitches, &
                                        posSwitches, t, branch_pgen)
            ! check validity
            check_abort_excit(branch_pgen, t)

        end do

        ! then update weights and and to lowering semi-stop
        weights = weights%ptr

        call calcRaisingSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                           posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        excitInfo%currentGen = excitInfo%lastGen

        do iOrb = ende1 + 1, ende2 - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)

            branch_pgen = branch_pgen * temp_pgen
            ! check validity
            check_abort_excit(branch_pgen, t)

        end do

        ! and finally to end step
        call singleStochasticEnd(csf_i, excitInfo, t)

        ! if we do RDMs also store the x0 and x1 coupling coeffs
        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                   excit_lvl=2, excit_typ=excitInfo%typ), &
                                    x0=extract_matrix_element(t, 1), &
                                    x1=extract_matrix_element(t, 2))
        end if

        ! for the additional contributing integrals:
        ! i have to consider that there might be a non-overlap excitation,
        ! which leads to the same excitation if there is no change in the
        ! stepvector in the overlap region!
        ! where i have to parts of the matrix elements, but which can be
        ! expressed in terms of the already calulated one i think..
        ! atleast in the case if R2L and L2R, as both the bottom and top
        ! parts are the same and only the semi-start ans semi-stop have to
        ! be modified..
        ! and also the x1 element has to be dismissed i guess..
        ! so i should change how i deal with the x1 elements and maybe
        ! keep them out here in these functions..

        ! determine if a switch happended:
        switch = findFirstSwitch(ilut, t, start2 + 1, ende1)
        ! i have to rethink the matrix element contribution by the
        ! equivalent non-overlap excitations, since its possible that
        ! a +-2 excitations stays on the +-2 branch in the double overlap
        ! region all the time and then also can be produced by a non-overlap
        ! but then the matrix element is different then the calculation used
        ! right now..
        ! no!!! sinnce the excitations leading to such a csf wouldn be valid
        ! single excitations-.. since the deltaB values at the end would not
        ! match

        if (switch > 0) then
            ! a switch happened and only mixed overlap contributes
            ! after determining how to deal with different x0 and x1 parts
            ! wait: if a switch happened i know that the x0-element is 0
            ! and only the x1-element of the overlap region counts!
            integral = extract_matrix_element(t, 2) * (get_umat_el(start1, ende2, ende1, start2) + &
                                                       get_umat_el(ende2, start1, start2, ende1)) / 2.0_dp

            if (near_zero(integral)) then
                branch_pgen = 0.0_dp
                t = 0
            else
                call encode_matrix_element(t, 0.0_dp, 2)
                call encode_matrix_element(t, integral, 1)
            end if

        else
            ! no switch happened: so i have to think about the
            ! non-overlap contribution
            ! in the R2L and L2R case it is really easy if i keep the
            ! x0 and x1 elements seperate up until here.
            ! since the types of generators are the same and only the x0
            ! element(which does not get changed in the overlap region)
            ! are nececarry for the non-overlap excitation , it only gets
            ! a -t^2 = -1/2 factor...
            ! so to get the non-overlap version i need to multiply by -2.0!
            integral = (-extract_matrix_element(t, 1) * (get_umat_el(start1, ende2, start2, ende1) + &
                                                     get_umat_el(ende2, start1, ende1, start2)) * 2.0_dp + (extract_matrix_element(t, 1) + &
                                                              extract_matrix_element(t, 2)) * (get_umat_el(start1, ende2, ende1, start2) + &
                                                                                        get_umat_el(ende2, start1, start2, ende1))) / 2.0_dp

            if (near_zero(integral)) then
                branch_pgen = 0.0_dp
                t = 0_n_int
            else
                call encode_matrix_element(t, 0.0_dp, 2)
                call encode_matrix_element(t, integral, 1)
            end if
        end if
    end subroutine calcDoubleR2L_stochastic

    subroutine calcDoubleL2R_stochastic(ilut, csf_i, excitInfo, t, branch_pgen, &
                                        posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcDoubleL2R_stochastic"

        integer :: iOrb, start2, ende1, ende2, start1, switch
        type(WeightObj_t) :: weights
        real(dp) :: temp_pgen
        HElement_t(dp) :: integral

        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(.not. isThree(ilut, excitInfo%secondStart))
        ASSERT(.not. isThree(ilut, excitInfo%firstEnd))
        ASSERT(.not. isZero(ilut, excitInfo%fullEnd))

        start1 = excitInfo%fullStart
        start2 = excitInfo%secondStart
        ende1 = excitInfo%firstEnd
        ende2 = excitInfo%fullEnd

        ! : create correct weights:
        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_fullDoubleWeight(csf_i, start2, ende1, ende2, negSwitches(start2), &
                                            negSwitches(ende1), posSwitches(start2), posSwitches(ende1), &
                                            csf_i%B_real(start2), csf_i%B_real(ende1))
        end if

        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = start1 + 1, start2 - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            branch_pgen = branch_pgen * temp_pgen

            ! check validity
            check_abort_excit(branch_pgen, t)

        end do

        ! change weights... maybe need both single and double type weights
        ! then do lowering semi start
        weights = weights%ptr

        call calcRaisingSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                            posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = start2 + 1, ende1 - 1
            call doubleUpdateStochastic(ilut, csf_i, iOrb, excitInfo, weights, negSwitches, &
                                        posSwitches, t, branch_pgen)
            ! check validity
            check_abort_excit(branch_pgen, t)

        end do

        ! then update weights and and to lowering semi-stop
        weights = weights%ptr

        call calcLoweringSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                            posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        excitInfo%currentGen = excitInfo%lastGen

        do iOrb = ende1 + 1, ende2 - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            branch_pgen = branch_pgen * temp_pgen
            ! check validity
            check_abort_excit(branch_pgen, t)

        end do

        ! and finally to end step
        call singleStochasticEnd(csf_i, excitInfo, t)

        ! if we do RDMs also store the x0 and x1 coupling coeffs
        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                    contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                           excit_lvl=2, excit_typ=excitInfo%typ), &
                                            x0=extract_matrix_element(t, 1), &
                                            x1=extract_matrix_element(t, 2))
        end if

        ! for the additional contributing integrals:
        ! i have to consider that there might be a non-overlap excitation,
        ! which leads to the same excitation if there is no change in the
        ! stepvector in the overlap region!
        ! where i have to parts of the matrix elements, but which can be
        ! expressed in terms of the already calulated one i think..
        ! atleast in the case if R2L and L2R, as both the bottom and top
        ! parts are the same and only the semi-start ans semi-stop have to
        ! be modified..
        ! and also the x1 element has to be dismissed i guess..
        ! so i should change how i deal with the x1 elements and maybe
        ! keep them out here in these functions..
        switch = findFirstSwitch(ilut, t, start2 + 1, ende1)
        if (switch > 0) then
            integral = extract_matrix_element(t, 2) * (get_umat_el(ende1, start2, start1, ende2) + &
                                                       get_umat_el(start2, ende1, ende2, start1)) / 2.0_dp

            if (near_zero(integral)) then
                branch_pgen = 0.0_dp
                t = 0
            else
                call encode_matrix_element(t, 0.0_dp, 2)
                call encode_matrix_element(t, integral, 1)
            end if
        else
            integral = (-extract_matrix_element(t, 1) * (get_umat_el(start2, ende1, start1, ende2) + &
                                                     get_umat_el(ende1, start2, ende2, start1)) * 2.0_dp + (extract_matrix_element(t, 1) + &
                                                              extract_matrix_element(t, 2)) * (get_umat_el(ende1, start2, start1, ende2) + &
                                                                                        get_umat_el(start2, ende1, ende2, start1))) / 2.0_dp

            if (near_zero(integral)) then
                branch_pgen = 0.0_dp
                t = 0_n_int
            else
                call encode_matrix_element(t, 0.0_dp, 2)
                call encode_matrix_element(t, integral, 1)
            end if
        end if
    end subroutine calcDoubleL2R_stochastic

    subroutine calcDoubleL2R2L_stochastic(ilut, csf_i, excitInfo, t, branch_pgen, &
                                          posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcDoubleL2R2L_stochastic"

        integer :: iOrb, start2, ende1, ende2, start1, switch
        type(WeightObj_t) :: weights
        real(dp) :: temp_pgen
        HElement_t(dp) :: integral

        ! have to create this additional routine to more efficiently
        ! incorporate the integral contributions, since it is vastly different
        ! when involving mixed generators, but thats still todo!
        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(.not. isThree(ilut, excitInfo%secondStart))
        ASSERT(.not. isZero(ilut, excitInfo%firstEnd))
        ASSERT(.not. isThree(ilut, excitInfo%fullEnd))

        start1 = excitInfo%fullStart
        start2 = excitInfo%secondStart
        ende1 = excitInfo%firstEnd
        ende2 = excitInfo%fullEnd

        if (present(opt_weight)) then
            weights = opt_weight
        else
            ! : create correct weights:
            weights = init_fullDoubleWeight(csf_i, start2, ende1, ende2, negSwitches(start2), &
                                            negSwitches(ende1), posSwitches(start2), posSwitches(ende1), &
                                            csf_i%B_real(start2), csf_i%B_real(ende1))
        end if

        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)

        ! check validity (defined in macros.h)
        check_abort_excit(branch_pgen, t)

        do iOrb = start1 + 1, start2 - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            ! check validity
            branch_pgen = branch_pgen * temp_pgen

            check_abort_excit(branch_pgen, t)

        end do

        ! change weights... maybe need both single and double type weights
        ! then do lowering semi start
        weights = weights%ptr

        call calcRaisingSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                            posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = start2 + 1, ende1 - 1
            call doubleUpdateStochastic(ilut, csf_i, iOrb, excitInfo, weights, negSwitches, &
                                        posSwitches, t, branch_pgen)
            ! check validity
            check_abort_excit(branch_pgen, t)

        end do

        ! then update weights and and to lowering semi-stop
        weights = weights%ptr

        call calcRaisingSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                           posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = ende1 + 1, ende2 - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            ! check validity
            branch_pgen = branch_pgen * temp_pgen

            check_abort_excit(branch_pgen, t)

        end do

        ! and finally to end step
        call singleStochasticEnd(csf_i, excitInfo, t)

        ! if we do RDMs also store the x0 and x1 coupling coeffs
        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                           excit_lvl=2, excit_typ=excitInfo%typ), &
                                            x0=extract_matrix_element(t, 1), &
                                            x1=extract_matrix_element(t, 2))
        end if

        ! todo: think about the additional integral contributions and the
        ! relative sign of different order influences...
        ! and not sure yet if in this case i can use this function generally
        ! for R, RR, RR, R
        ! and L -> RL -> RL -> L
        ! since the integral/order influences are different maybe... todo!
        ! see above too!
        ! but yeah matrix elements definitly depends on excitation here.
        ! if there is no change in the overlap region, there is also the
        ! non-overlap contribution! otherwise not. maybe set some flag to
        ! indicate if such a stepvector change happened or not.

        ! update: on the matrix elements..
        ! i have to consider the non-overlap excitation, which can lead to
        ! the same excitation if there is no change in the stepvector
        ! in the overlap region
        ! the problem with the L2R2L and R2L2R funcitons is that the
        ! generator type changes for the non-overlap excitation so the
        ! matrix elements have to be changed more than in the R2L and L2R case

        ! if a switch happend -> also no additional contribution
        switch = findFirstSwitch(ilut, t, start2 + 1, ende1)

        if (switch > 0) then
            integral = extract_matrix_element(t, 2) * (get_umat_el(ende2, start2, start1, ende1) + &
                                                       get_umat_el(start2, ende2, ende1, start1)) / 2.0_dp

            if (near_zero(integral)) then
                branch_pgen = 0.0_dp
                t = 0_n_int
            else
                call encode_matrix_element(t, 0.0_dp, 2)
                call encode_matrix_element(t, integral, 1)
            end if

        else
            ! the no switch happened i have to get the additional contributions
            ! by the non-overlap version, but now the generator type at the
            ! end is different as the on-the fly calculated ...
            ! so the x0 matrix element changes by more (or even recalculate?)
            ! the bottom contribution, stays the same -sqrt(2) to cancel the
            ! -t
            ! the contribution at the semi-stop stays the same +t
            ! so those to get to -2.0_dp
            ! and bullshit that generator changes!! is also the same
            ! so it is exactly the same as in the R2L and L2R cases!! phew

            ! have also check if the non-overlap matrix elements is 0...
            ! this can happen unfortunately
            integral = (-extract_matrix_element(t, 1) * (get_umat_el(start2, ende2, start1, ende1) + &
                                                     get_umat_el(ende2, start2, ende1, start1)) * 2.0_dp + (extract_matrix_element(t, 1) + &
                                                              extract_matrix_element(t, 2)) * (get_umat_el(ende2, start2, start1, ende1) + &
                                                                                        get_umat_el(start2, ende2, ende1, start1))) / 2.0_dp

            if (near_zero(integral)) then
                branch_pgen = 0.0_dp
                t = 0_n_int
            else
                call encode_matrix_element(t, 0.0_dp, 2)
                call encode_matrix_element(t, integral, 1)
            end if
        end if

    end subroutine calcDoubleL2R2L_stochastic

    subroutine calcDoubleRaisingStochastic(ilut, csf_i, excitInfo, t, branch_pgen, &
                                           posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcDoubleRaisingStochastic"

        integer :: iOrb, start2, ende1, ende2, start1
        type(WeightObj_t) :: weights
        real(dp) :: temp_pgen
        HElement_t(dp) :: integral

        ASSERT(.not. isThree(ilut, excitInfo%fullStart))
        ASSERT(.not. isThree(ilut, excitInfo%secondStart))
        ASSERT(.not. isZero(ilut, excitInfo%firstEnd))
        ASSERT(.not. isZero(ilut, excitInfo%fullEnd))

        start1 = excitInfo%fullStart
        start2 = excitInfo%secondStart
        ende1 = excitInfo%firstEnd
        ende2 = excitInfo%fullEnd

        ! : create correct weights:
        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_fullDoubleWeight(csf_i, start2, ende1, ende2, negSwitches(start2), &
                                            negSwitches(ende1), posSwitches(start2), posSwitches(ende1), &
                                            csf_i%B_real(start2), csf_i%B_real(ende1))
        end if

        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = start1 + 1, start2 - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            ! check validity
            branch_pgen = branch_pgen * temp_pgen

            check_abort_excit(branch_pgen, t)

        end do

        ! change weights... maybe need both single and double type weights
        ! then do lowering semi start
        ! just point to the next weight:
        weights = weights%ptr

        call calcRaisingSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                            posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = start2 + 1, ende1 - 1
            call doubleUpdateStochastic(ilut, csf_i, iOrb, excitInfo, weights, negSwitches, &
                                        posSwitches, t, branch_pgen)
            ! check validity

            check_abort_excit(branch_pgen, t)

        end do

        ! then update weights and and to lowering semi-stop
        weights = weights%ptr

        call calcRaisingSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                           posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = ende1 + 1, ende2 - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            ! check validity
            branch_pgen = branch_pgen * temp_pgen

            check_abort_excit(branch_pgen, t)

        end do

        ! and finally to end step
        call singleStochasticEnd(csf_i, excitInfo, t)

        ! if we do RDMs also store the x0 and x1 coupling coeffs
        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                       excit_lvl=2, excit_typ=excitInfo%typ), &
                                            x0=extract_matrix_element(t, 1), &
                                            x1=extract_matrix_element(t, 2))
        end if

        ! todo: think about the additional integral contributions and the
        ! relative sign of different order influences...
        ! and not sure yet if in this case i can use this function generally
        ! for R, RR, RR, R
        ! and L -> RL -> RL -> L
        ! since the integral/order influences are different maybe... todo!

        ! think more about the sign influence in these kind of excitations..
        ! according to my notes, if i keep the x0 and x1 elements seperate
        ! until out here i can express it neatly in form of:
        ! x0(U1 + U2) + x1(U1 - U2)
        ! but not quite sure about that anymore... hopefully yes!
        integral = (extract_matrix_element(t, 1) * (get_umat_el(start1, start2, ende1, ende2) + &
                                                   get_umat_el(start2, start1, ende2, ende1) + get_umat_el(start1, start2, ende2, ende1) + &
                                                    get_umat_el(start2, start1, ende1, ende2)) + excitInfo%order * excitInfo%order1 * &
                    extract_matrix_element(t, 2) * ( &
                    get_umat_el(start1, start2, ende1, ende2) + get_umat_el(start2, start1, ende2, ende1) - &
                    get_umat_el(start1, start2, ende2, ende1) - get_umat_el(start2, start1, ende1, ende2))) / 2.0_dp

        if (near_zero(integral)) then
            branch_pgen = 0.0_dp
            t = 0_n_int
        else
            call encode_matrix_element(t, 0.0_dp, 2)
            call encode_matrix_element(t, integral, 1)
        end if

    end subroutine calcDoubleRaisingStochastic

    subroutine calcDoubleR2L2R_stochastic(ilut, csf_i, excitInfo, t, branch_pgen, &
                                          posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcDoubleR2L2R_stochastic"

        integer :: iOrb, switch
        type(WeightObj_t) :: weights
        real(dp) :: temp_pgen
        HElement_t(dp) :: integral

        associate (i => excitInfo%i, j => excitInfo%j, k => excitInfo%k, &
                   l => excitInfo%l, start1 => excitInfo%fullstart, &
                   start2 => excitInfo%secondStart, ende1 => excitInfo%firstEnd, &
                   ende2 => excitInfo%fullEnd, typ => excitInfo%typ)

            ASSERT(.not. isThree(ilut, start1))
            ASSERT(.not. isZero(ilut, start2))
            ASSERT(.not. isThree(ilut, ende1))
            ASSERT(.not. isZero(ilut, ende2))

            if (present(opt_weight)) then
                weights = opt_weight
            else
                ! : create correct weights:
                weights = init_fullDoubleWeight(csf_i, start2, ende1, ende2, negSwitches(start2), &
                                                negSwitches(ende1), posSwitches(start2), posSwitches(ende1), &
                                                csf_i%B_real(start2), csf_i%B_real(ende1))
            end if

            call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                              negSwitches, t, branch_pgen)

            ! check validity
            check_abort_excit(branch_pgen, t)

            do iOrb = start1 + 1, start2 - 1
                call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                            negSwitches, t, temp_pgen)
                ! check validity
                branch_pgen = branch_pgen * temp_pgen

                check_abort_excit(branch_pgen, t)

            end do

            ! change weights... maybe need both single and double type weights
            ! then do lowering semi start
            weights = weights%ptr

            call calcLoweringSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                                 posSwitches, t, branch_pgen)

            ! check validity
            check_abort_excit(branch_pgen, t)

            do iOrb = start2 + 1, ende1 - 1
                call doubleUpdateStochastic(ilut, csf_i, iOrb, excitInfo, weights, negSwitches, &
                                            posSwitches, t, branch_pgen)
                ! check validity

                check_abort_excit(branch_pgen, t)

            end do

            ! then update weights and and to lowering semi-stop
            weights = weights%ptr

            call calcLoweringSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                                posSwitches, t, branch_pgen)

            ! check validity
            check_abort_excit(branch_pgen, t)

            do iOrb = ende1 + 1, ende2 - 1
                call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                            negSwitches, t, temp_pgen)

                branch_pgen = branch_pgen * temp_pgen
                ! check validity
                check_abort_excit(branch_pgen, t)

            end do

            ! and finally to end step
            call singleStochasticEnd(csf_i, excitInfo, t)

            ! if we do RDMs also store the x0 and x1 coupling coeffs
            if (tFillingStochRDMOnFly) then
                call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                        contract_2_rdm_ind(i, j, k, l, excit_lvl=2, excit_typ=typ), &
                                                x0=extract_matrix_element(t, 1), &
                                                x1=extract_matrix_element(t, 2))
            end if

            ! todo: think about the additional integral contributions and the
            ! relative sign of different order influences...
            ! and not sure yet if in this case i can use this function generally
            ! for L -> LL -> LL -> L
            ! and R -> RL -> RL -> R
            ! since the integral/order influences are different maybe... todo!

            ! update: on the matrix elements..
            ! i have to consider the non-overlap excitation, which can lead to
            ! the same excitation if there is no change in the stepvector
            ! in the overlap region
            ! the problem with the L2R2L and R2L2R funcitons is that the
            ! generator type changes for the non-overlap excitation so the
            ! matrix elements have to be changed more than in the R2L and L2R case

            switch = findFirstSwitch(ilut, t, start2 + 1, ende1)
            if (switch > 0) then
                integral = extract_matrix_element(t, 2) * (get_umat_el(start1, ende1, ende2, start2) + &
                                                           get_umat_el(ende1, start1, start2, ende2)) / 2.0_dp

                if (near_zero(integral)) then
                    branch_pgen = 0.0_dp
                    t = 0_n_int
                else
                    call encode_matrix_element(t, 0.0_dp, 2)
                    call encode_matrix_element(t, integral, 1)
                end if
            else
                integral = (-extract_matrix_element(t, 1) * (get_umat_el(start1, ende1, start2, ende2) + &
                                                     get_umat_el(ende1, start1, ende2, start2)) * 2.0_dp + (extract_matrix_element(t, 1) + &
                                                              extract_matrix_element(t, 2)) * (get_umat_el(start1, ende1, ende2, start2) + &
                                                                                        get_umat_el(ende1, start1, start2, ende2))) / 2.0_dp

                if (near_zero(integral)) then
                    branch_pgen = 0.0_dp
                    t = 0_n_int
                else
                    call encode_matrix_element(t, 0.0_dp, 2)
                    call encode_matrix_element(t, integral, 1)
                end if
            end if

        end associate

    end subroutine calcDoubleR2L2R_stochastic

    subroutine calcDoubleLoweringStochastic(ilut, csf_i, excitInfo, t, branch_pgen, &
                                            posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcDoubleLoweringStochastic"

        integer :: iOrb, start2, ende1, ende2, start1
        type(WeightObj_t) :: weights
        real(dp) :: temp_pgen
        HElement_t(dp) :: integral

        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(.not. isZero(ilut, excitInfo%secondStart))
        ASSERT(.not. isThree(ilut, excitInfo%firstEnd))
        ASSERT(.not. isThree(ilut, excitInfo%fullEnd))

        start1 = excitInfo%fullStart
        start2 = excitInfo%secondStart
        ende1 = excitInfo%firstEnd
        ende2 = excitInfo%fullEnd

        ! : create correct weights:
        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_fullDoubleWeight(csf_i, start2, ende1, ende2, negSwitches(start2), &
                                            negSwitches(ende1), posSwitches(start2), posSwitches(ende1), &
                                            csf_i%B_real(start2), csf_i%B_real(ende1))
        end if

        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = start1 + 1, start2 - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            ! check validity
            branch_pgen = branch_pgen * temp_pgen

            check_abort_excit(branch_pgen, t)
        end do

        ! change weights... maybe need both single and double type weights
        ! then do lowering semi start
        ! can i just do:
        weights = weights%ptr

        ! branch_pgen gets update insde the routine!
        call calcLoweringSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                             posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = start2 + 1, ende1 - 1
            ! branch_pgen gets updated inside update routine
            call doubleUpdateStochastic(ilut, csf_i, iOrb, excitInfo, weights, negSwitches, &
                                        posSwitches, t, branch_pgen)
            ! here only need to have probweight, since i cant only check x1 element
            ! check validity

            check_abort_excit(branch_pgen, t)
        end do

        ! then update weights and and to lowering semi-stop
        weights = weights%ptr

        ! branch_pgen gets updated inside funciton
        call calcLoweringSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                            posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        do iOrb = ende1 + 1, ende2 - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            ! check validity
            branch_pgen = branch_pgen * temp_pgen

            check_abort_excit(branch_pgen, t)
        end do

        ! and finally to end step
        call singleStochasticEnd(csf_i, excitInfo, t)

        ! if we do RDMs also store the x0 and x1 coupling coeffs
        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                       excit_lvl=2, excit_typ=excitInfo%typ), &
                                            x0=extract_matrix_element(t, 1), &
                                            x1=extract_matrix_element(t, 2))
        end if

        ! todo: think about the additional integral contributions and the
        ! relative sign of different order influences...
        ! and not sure yet if in this case i can use this function generally
        ! for L -> LL -> LL -> L
        ! and R -> RL -> RL -> R
        ! since the integral/order influences are different maybe... todo!

        ! think more about the sign influence in these kind of excitations..
        ! according to my notes, if i keep the x0 and x1 elements seperate
        ! until out here i can express it neatly in form of:
        ! x0(U1 + U2) + x1(U1 - U2)
        ! but not quite sure about that anymore... hopefully yes!
        ! try that for now!

        integral = (extract_matrix_element(t, 1) * (get_umat_el(ende1, ende2, start1, start2) + &
                                                   get_umat_el(ende2, ende1, start2, start1) + get_umat_el(ende2, ende1, start1, start2) + &
                                                    get_umat_el(ende1, ende2, start2, start1)) + excitInfo%order * excitInfo%order1 * &
                    extract_matrix_element(t, 2) * ( &
                    get_umat_el(ende1, ende2, start1, start2) + get_umat_el(ende2, ende1, start2, start1) - &
                    get_umat_el(ende2, ende1, start1, start2) - get_umat_el(ende1, ende2, start2, start1))) / 2.0_dp

        if (near_zero(integral)) then
            branch_pgen = 0.0_dp
            t = 0_n_int
        else
            call encode_matrix_element(t, 0.0_dp, 2)
            call encode_matrix_element(t, integral, 1)
        end if

    end subroutine calcDoubleLoweringStochastic

    subroutine calcFullStopL2R_stochastic(ilut, csf_i, excitInfo, t, pgen, &
                                          posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight

        type(WeightObj_t) :: weights
        integer :: st, se, en, i, j, k, l, elecInd, holeInd
        real(dp) :: branch_pgen, &
                    temp_pgen, rdm_mat, p_orig, orb_pgen
        HElement_t(dp) :: integral

        st = excitInfo%fullStart
        se = excitInfo%secondStart
        en = excitInfo%fullEnd

        ! init weights
        if (present(opt_weight)) then
            weights = opt_weight
        else
            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                ! the weights should be the only necessary change to force
                ! a switch at the end, as the other branches get 0 weight..
                weights = init_forced_end_semistart_weight(csf_i, se, en, negSwitches(se), &
                                                           posSwitches(se), csf_i%B_real(se))

            else
                weights = init_semiStartWeight(csf_i, se, en, negSwitches(se), &
                                               posSwitches(se), csf_i%B_real(se))
            end if
        end if

        ! create st
        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)

        ! in case of early access the pgen should be set to 0
        pgen = 0.0_dp

        ! check validity
        check_abort_excit(branch_pgen, t)

        do i = st + 1, se - 1
            call singleStochasticUpdate(ilut, csf_i, i, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            branch_pgen = branch_pgen * temp_pgen
            ! check validity

            check_abort_excit(branch_pgen, t)
        end do

        ! do the specific se-st
        ! try the new reusing of the weights object..
        weights = weights%ptr

        call calcRaisingSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                            posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        ! do the specific double update to ensure a switch
        ! although switch can also happen at end only...
        ! actually that would be, in the full-stop case, temporary measure...
        ! but would unjust favor certain types of excitations..
        do i = se + 1, en - 1
            call doubleUpdateStochastic(ilut, csf_i, i, excitInfo, &
                                        weights, negSwitches, posSwitches, t, branch_pgen)

            if (near_zero(extract_matrix_element(t, 2)) .or. near_zero(branch_pgen)) then
                t = 0_n_int
                return
            end if
        end do

        call mixedFullStopStochastic(ilut, csf_i, excitInfo, t)

        ! check if matrix element is non-zero and if a switch happened
        if (.not. near_zero(extract_matrix_element(t, 1))) then
            t = 0_n_int
            branch_pgen = 0.0_dp
            return
        end if

        if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
            if (getDeltaB(t) == 0) then
                t = 0_n_int
                branch_pgen = 0.0_dp
                return
            end if
        end if

        if (near_zero(extract_matrix_element(t, 2))) then
            branch_pgen = 0.0_dp
            t = 0_n_int
            return
        end if

        ! if we do RDMs also store the x0 and x1 coupling coeffs
        ! and I need to do it before the routines below since excitInfo
        ! gets changed there
        if (tFillingStochRDMOnFly) then
            ! i need to unbias against the total pgen later on in the
            ! RDM sampling otherwise the rdm-bias factor is not correct!
            ! encode the necessary information in the rdm-matele!
            i = excitInfo%i
            j = excitInfo%j
            k = excitInfo%k
            l = excitInfo%l
            elecInd = st
            holeInd = se
            rdm_mat = extract_matrix_element(t, 2)
            call calc_orbital_pgen_contrib_end(&
                    csf_i, [2 * elecInd, 2 * en], holeInd, orb_pgen)
            p_orig = orb_pgen * branch_pgen / real(ElecPairs, dp)
            if (csf_i%stepvector(elecInd) == 3) p_orig = p_orig * 2.0_dp
        end if

        call encode_matrix_element(t, extract_matrix_element(t, 2), 1)

        ! actually I should provide a new routine, which "just"
        ! calculates the matrix element contribution and not
        ! the modified pgen, as the spatial orbitals are now fixed
        ! this could be done in the initialisation, where i just
        ! point to a new function, which only calculates the
        ! matrix element contribution

        global_excitInfo = excitInfo

        if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
            call calc_mixed_end_contr_approx(t, csf_i, excitInfo, integral)
            pgen = branch_pgen

        else
            call calc_mixed_end_contr_sym(ilut, csf_i, t, excitInfo, branch_pgen, pgen, &
                                          integral)
        end if

        if (tFillingStochRDMOnFly) then
            if (.not. near_zero(p_orig)) then
                call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                        contract_2_rdm_ind(i, j, k, l, excit_lvl=2, &
                                       excit_typ=excitInfo%typ), x0=0.0_dp, &
                                                x1=rdm_mat * pgen / p_orig)
            end if
        end if

        call encode_matrix_element(t, 0.0_dp, 2)
        call update_matrix_element(t, (get_umat_el(en, se, st, en) + &
                                       get_umat_el(se, en, en, st)) / 2.0_dp + integral, 1)

    end subroutine calcFullStopL2R_stochastic

    subroutine calc_mixed_end_contr_approx(t, csf_i, excitInfo, integral)
        ! for the approx. mixed end contribution i "just" need to
        ! calculate the correct matrix element influences
        integer(n_int), intent(in) :: t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        HElement_t(dp), intent(out) :: integral
        character(*), parameter :: this_routine = "calc_mixed_end_contr_approx"

        integer :: st, se, en, elecInd, holeInd, step, sw, i
        real(dp) :: top_cont, mat_ele, stay_mat, end_mat
        logical :: above_flag

        ! do as much stuff as possible beforehand
        st = excitInfo%fullStart
        se = excitInfo%secondStart
        en = excitInfo%fullEnd
        if (excitInfo%typ == excit_type%fullstop_L_to_R) then
            elecInd = st
            holeInd = se
        else if (excitInfo%typ == excit_type%fullstop_R_to_L) then
            elecInd = se
            holeInd = st
        else
            call stop_all(this_routine, "should not be here!")
        end if

        integral = h_cast(0.0_dp)

        step = csf_i%stepvector(en)

        ! i am sure the last switch happens at the full-stop!
        sw = en

        if (en < nSpatOrbs) then
            select case (step)
            case (1)
                if (isOne(t, en)) then
                    top_cont = -Root2 * sqrt((csf_i%B_real(en) + 2.0_dp) / &
                                             csf_i%B_real(en))

                else
                    top_cont = -Root2 / sqrt(csf_i%B_real(en) * (csf_i%B_real(en) + 2.0_dp))

                end if
            case (2)
                if (isOne(t, en)) then
                    top_cont = -Root2 / sqrt(csf_i%B_real(en) * (csf_i%B_real(en) + 2.0_dp))

                else
                    top_cont = Root2 * sqrt(csf_i%B_real(en) / &
                                            (csf_i%B_real(en) + 2.0_dp))
                end if

            case default
                call stop_all(this_routine, "wrong stepvalues!")

            end select

            if (.not. near_zero(top_cont)) then

                above_flag = .false.
                mat_ele = 1.0_dp

                do i = en + 1, nSpatOrbs
                    if (csf_i%Occ_int(i) /= 1) cycle

                    ! then check if thats the last step
                    if (csf_i%stepvector(i) == 2 .and. csf_i%B_int(i) == 0) then
                        above_flag = .true.
                    end if

                    ! in the other routine i check if the orbital pgen
                    ! is 0 for the above orbitals.. do I need to do that
                    ! also here?? or is this implicit if the matrix
                    ! element will be 0??

                    step = csf_i%stepvector(i)

                    call getDoubleMatrixElement(step, step, 0, gen_type%L, gen_type%R, csf_i%B_real(i), &
                                                1.0_dp, x1_element=stay_mat)

                    call getMixedFullStop(step, step, 0, csf_i%B_real(i), &
                                          x1_element=end_mat)

                    ! this check should never be true, but just to be sure
                    if (near_zero(stay_mat)) above_flag = .true.

                    if (.not. near_zero(end_mat)) then
                        integral = integral + end_mat * mat_ele * &
                                   (get_umat_el(i, holeInd, elecInd, i) + &
                                    get_umat_el(holeInd, i, i, elecInd)) / 2.0_dp
                    end if

                    if (above_flag) exit

                    ! otherwise update your running matrix element vars
                    mat_ele = mat_ele * stay_mat

                end do

                integral = integral * top_cont
            end if
        end if

    end subroutine calc_mixed_end_contr_approx

    subroutine calcFullStopR2L_stochastic(ilut, csf_i, excitInfo, t, pgen, &
                                          posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight

        type(WeightObj_t) :: weights
        integer :: st, se, en, i, j, k, l, elecInd, holeInd
        real(dp) :: branch_pgen, &
                    temp_pgen, p_orig, rdm_mat, orb_pgen
        HElement_t(dp) :: integral

        st = excitInfo%fullStart
        se = excitInfo%secondStart
        en = excitInfo%fullEnd

        ! init weights
        if (present(opt_weight)) then
            weights = opt_weight
        else
            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                weights = init_forced_end_semistart_weight(csf_i, se, en, negSwitches(se), &
                                                           posSwitches(se), csf_i%B_real(se))
            else
                weights = init_semiStartWeight(csf_i, se, en, negSwitches(se), &
                                               posSwitches(se), csf_i%B_real(se))
            end if
        end if

        ! create start
        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)

        ! in case of early exit pgen should be set to 0
        pgen = 0.0_dp
        ! check validity

        check_abort_excit(branch_pgen, t)

        do i = st + 1, se - 1
            call singleStochasticUpdate(ilut, csf_i, i, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            branch_pgen = branch_pgen * temp_pgen
            ! check validity

            check_abort_excit(branch_pgen, t)
        end do

        ! do the specific semi-start
        weights = weights%ptr

        call calcLoweringSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                             posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        ! do the specific double update to ensure a switch
        ! although switch can also happen at end only...
        ! actually that would be, in the full-stop case, temporary measure...
        ! but would unjust favor certain types of excitations..
        do i = se + 1, en - 1
            call doubleUpdateStochastic(ilut, csf_i, i, excitInfo, &
                                        weights, negSwitches, posSwitches, t, branch_pgen)
            ! also here there has to be a switch at some point so check x1
            if (near_zero(extract_matrix_element(t, 2)) .or. near_zero(branch_pgen)) then
                t = 0_n_int
                return
            end if
        end do

        call mixedFullStopStochastic(ilut, csf_i, excitInfo, t)

        ! check if matrix element is non-zero and if a switch happened
        if (.not. near_zero(extract_matrix_element(t, 1))) then
            t = 0_n_int
            return
        end if
        if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
            ! make it crude for now, that we only check if the delta B value
            ! is non-zero at the end, otherwise abort this spawn..
            if (getDeltaB(t) == 0) then
                t = 0_n_int
                return
            end if
        end if

        if (near_zero(extract_matrix_element(t, 2))) then
            t = 0_n_int
            return
        end if

        ! if we do RDMs also store the x0 and x1 coupling coeffs
        ! and I need to do it before the routines below since excitInfo
        ! gets changed there
        if (tFillingStochRDMOnFly) then
            ! i need to unbias against the total pgen later on in the
            ! RDM sampling otherwise the rdm-bias factor is not correct!
            ! encode the necessary information in the rdm-matele!
            i = excitInfo%i
            j = excitInfo%j
            k = excitInfo%k
            l = excitInfo%l
            elecInd = se
            holeInd = st
            rdm_mat = extract_matrix_element(t, 2)
            call calc_orbital_pgen_contrib_end(&
                    csf_i, [2 * elecInd, 2 * en], holeInd, orb_pgen)
            p_orig = orb_pgen * branch_pgen / real(ElecPairs, dp)
            if (csf_i%stepvector(elecInd) == 3) p_orig = p_orig * 2.0_dp

        end if

        ! the x1-element is still encoded in the second entry..
        ! move it to the first elemen
        call encode_matrix_element(t, extract_matrix_element(t, 2), 1)

        global_excitInfo = excitInfo

        if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then

            call calc_mixed_end_contr_approx(t, csf_i, excitInfo, integral)
            pgen = branch_pgen

        else
            call calc_mixed_end_contr_sym(ilut, csf_i, t, excitInfo, branch_pgen, pgen, integral)
        end if

        if (tFillingStochRDMOnFly) then
            if (.not. near_zero(p_orig)) then
                call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                        contract_2_rdm_ind(i, j, k, l, excit_lvl=2, &
                                           excit_typ=excitInfo%typ), x0=0.0_dp, &
                                                x1=rdm_mat * pgen / p_orig)
            end if
        end if

        call encode_matrix_element(t, 0.0_dp, 2)
        call update_matrix_element(t, (get_umat_el(en, st, se, en) + &
                                       get_umat_el(st, en, en, se)) / 2.0_dp + integral, 1)

    end subroutine calcFullStopR2L_stochastic

    subroutine doubleUpdateStochastic(ilut, csf_i, s, excitInfo, weights, negSwitches, &
                                      posSwitches, t, probWeight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: s
        type(ExcitationInformation_t), intent(in) :: excitInfo
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: negSwitches(nSpatOrbs), posSwitches(nSpatOrbs)
        integer(n_int), intent(inout) :: t(0:nifguga)
        real(dp), intent(inout) :: probWeight
        character(*), parameter :: this_routine = "doubleUpdateStochastic"

        real(dp) :: minusWeight, plusWeight, zeroWeight
        integer :: gen1, gen2, deltaB
        real(dp) :: bVal, tempWeight_0, tempWeight_1, order

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(s > 0 .and. s <= nSpatOrbs)

        if (csf_i%Occ_int(s) /= 1) then
            ! no change in stepvector or matrix element in this case
            return
        end if

        ! and for more readibility extract certain values:
        gen1 = excitInfo%gen1
        gen2 = excitInfo%gen2
        bVal = csf_i%B_real(s)
        ! stupid! only need at order at semistarts and semistops and not for
        ! the overlap region
        order = 1.0_dp

        deltaB = getDeltaB(t)

        ! new idea: make the combination stepvalue + deltaB !
        ! this give me 6 distinct integer quantities which i can choose
        ! from in a select case statement!

        select case (csf_i%stepvector(s) + deltaB)
            ! depending on the deltaB value different possibs
        case (3)
            ! d=1 + b=2 = 3
            ! only staying
            call getDoubleMatrixElement(1, 1, deltaB, gen1, gen2, bVal, &
                                        order, x1_element=tempWeight_1)

            tempWeight_0 = 0.0_dp

        case (-1)
            ! d=1 + b=-2 = -1
            ! both 0 and -2 branches are possible
            minusWeight = weights%proc%minus(negSwitches(s), &
                                             bVal, weights%dat)
            zeroWeight = weights%proc%zero(negSwitches(s), &
                                           posSwitches(s), bVal, weights%dat)

            if (near_zero(minusWeight + zeroWeight)) then
                probWeight = 0.0_dp
                t = 0
                return
            end if

            minusWeight = calcStayingProb(minusWeight, zeroWeight, bVal)

            if (genrand_real2_dSFMT() < minusWeight) then
                ! stay on -2 branch

                call getDoubleMatrixElement(1, 1, deltaB, gen1, gen2, &
                                            bVal, order, x1_element=tempWeight_1)

                tempWeight_0 = 0.0_dp

                probWeight = probWeight * minusWeight

            else
                ! switch to 0 branch
                ! 1 -> 2
                clr_orb(t, 2 * s - 1)
                set_orb(t, 2 * s)

                call setDeltaB(0, t)

                call getDoubleMatrixElement(2, 1, deltaB, gen1, gen2, &
                                            bVal, order, x1_element=tempWeight_1)

                tempWeight_0 = 0.0_dp

                probWeight = probWeight * (1.0_dp - minusWeight)

            end if

        case (1)
            ! d=1 + b=0 = 1
            ! 0 branch arrives -> have to check b value
            ! hm... is this a bug?? if d = 1, b cant be 0, its atleast
            ! 1.. and then a switch to +2 cant happen..
            ! but that should have been dealt with the weights below
            ! probably.. so thats why it didnt matter probably..
            if (csf_i%B_int(s) == 1) then
                ! only staying branch
                call getDoubleMatrixElement(1, 1, deltaB, gen1, gen2, bVal, &
                                            order, tempWeight_0, tempWeight_1)

            else
                ! both 0 and +2 branch are possible
                plusWeight = weights%proc%plus(posSwitches(s), &
                                               bVal, weights%dat)
                zeroWeight = weights%proc%zero(negSwitches(s), &
                                               posSwitches(s), bVal, weights%dat)

                if (near_zero(plusWeight + zeroWeight)) then
                    probWeight = 0.0_dp
                    t = 0
                end if

                zeroWeight = calcStayingProb(zeroWeight, plusWeight, bVal)

                if (genrand_real2_dSFMT() < zeroWeight) then
                    ! stay on 0 branch
                    call getDoubleMatrixElement(1, 1, deltaB, gen1, gen2, &
                                                bVal, order, tempWeight_0, tempWeight_1)

                    probWeight = probWeight * zeroWeight

                else
                    ! switch to +2 branch
                    ! 1 -> 2
                    clr_orb(t, 2 * s - 1)
                    set_orb(t, 2 * s)

                    call setDeltaB(2, t)

                    call getDoubleMatrixElement(2, 1, deltaB, gen1, gen2, &
                                                bVal, order, x1_element=tempWeight_1)

                    tempWeight_0 = 0.0_dp

                    probWeight = probWeight * (1.0_dp - zeroWeight)
                end if
            end if
        case (0)
            ! d=2 + b=-2 : 0
            ! only staying
            call getDoubleMatrixElement(2, 2, deltaB, gen1, gen2, bVal, &
                                        order, x1_element=tempWeight_1)

            tempWeight_0 = 0.0_dp

        case (2)
            ! d=2 + b=0 : 2
            ! always -2 and 0 branching possible
            minusWeight = weights%proc%minus(negSwitches(s), &
                                             bVal, weights%dat)
            zeroWeight = weights%proc%zero(negSwitches(s), &
                                           posSwitches(s), bVal, weights%dat)

            ! here should be the only case where b*w_0 + w_- = 0 ->
            ! have to avoid divison by 0 then
            ! but in this case only the stayin on 0 branch is valid
            ! since it should be always > 0 and the -branch has 0 weight
            if (near_zero(bVal * zeroWeight + minusWeight)) then
                zeroWeight = 1.0_dp
            else
                zeroWeight = calcStayingProb(zeroWeight, minusWeight, bVal)
            end if

            if (near_zero(zeroWeight + minusWeight)) then
                probWeight = 0.0_dp
                t = 0
                return
            end if

            if (genrand_real2_dSFMT() < zeroWeight) then
                ! stay on 0 branch
                call getDoubleMatrixElement(2, 2, deltaB, gen1, gen2, bVal, &
                                            order, tempWeight_0, tempWeight_1)

                probWeight = probWeight * zeroWeight

            else
                ! switch to -2 branch
                ! 2 -> 1
                set_orb(t, 2 * s - 1)
                clr_orb(t, 2 * s)

                call setDeltaB(-2, t)

                call getDoubleMatrixElement(1, 2, deltaB, gen1, gen2, &
                                            bVal, order, x1_element=tempWeight_1)

                tempWeight_0 = 0.0_dp

                probWeight = probWeight * (1.0_dp - zeroWeight)

            end if

        case (4)
            ! d=2 + b=2 : 4

            ! have to check b value if branching is possible
            if (csf_i%B_int(s) < 2) then

                ! only switch possible
                ! 2 -> 1
                set_orb(t, 2 * s - 1)
                clr_orb(t, 2 * s)

                call setDeltaB(0, t)

                call getDoubleMatrixElement(1, 2, deltaB, gen1, gen2, &
                                            bVal, order, x1_element=tempWeight_1)

                tempWeight_0 = 0.0_dp

            else
                ! both 0 and +2 branches possible
                plusWeight = weights%proc%plus(posSwitches(s), &
                                               bVal, weights%dat)
                zeroWeight = weights%proc%zero(negSwitches(s), &
                                               posSwitches(s), bVal, weights%dat)

                if (near_zero(plusWeight + zeroWeight)) then
                    probWeight = 0.0_dp
                    t = 0
                    return
                end if

                plusWeight = calcStayingProb(plusWeight, zeroWeight, bVal)

                if (genrand_real2_dSFMT() < plusWeight) then
                    ! stay on +2 branch
                    call getDoubleMatrixElement(2, 2, deltaB, gen1, gen2, &
                                                bVal, order, x1_element=tempWeight_1)

                    tempWeight_0 = 0.0_dp

                    probWeight = probWeight * plusWeight

                else
                    ! switch to 0 branch
                    ! 2 -> 1
                    set_orb(t, 2 * s - 1)
                    clr_orb(t, 2 * s)

                    call setDeltaB(0, t)
                    call getDoubleMatrixElement(1, 2, deltaB, gen1, gen2, &
                                                bVal, order, x1_element=tempWeight_1)

                    tempWeight_0 = 0.0_dp

                    probWeight = probWeight * (1.0_dp - plusWeight)

                end if
            end if
        end select

        call update_matrix_element(t, tempWeight_0, 1)
        call update_matrix_element(t, tempWeight_1, 2)

        if (near_zero(tempWeight_0) .and. near_zero(tempWeight_1)) then
            probWeight = 0.0_dp
            t = 0
            return
        end if

        if (t_trunc_guga_pgen .or. &
            (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
            if (probWeight < trunc_guga_pgen) then
                probWeight = 0.0_dp
                t = 0_n_int
            end if
        end if

        if (t_trunc_guga_matel) then
            if (abs(extract_matrix_element(t, 1)) < trunc_guga_matel .and. &
                abs(extract_matrix_element(t, 2)) < trunc_guga_matel) then
                probWeight = 0.0_dp
                t = 0_n_int
            end if
        end if
    end subroutine doubleUpdateStochastic

    subroutine mixedFullStopStochastic(ilut, csf_i, excitInfo, t)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(inout) :: t(0:nifguga)
        character(*), parameter :: this_routine = "mixedFullStopStochastic"

        integer :: ende, deltaB
        real(dp) :: bVal, tempWeight_0, tempWeight_1

        ASSERT(.not. isThree(ilut, excitInfo%fullEnd))
        ASSERT(.not. isZero(ilut, excitInfo%fullEnd))
        ! no 3 allowed at end or else it would be single-excitation-like

        ende = excitInfo%fullEnd
        bVal = csf_i%B_real(ende)

        deltaB = getDeltaB(t)

        ! NOTE: also remember for mixed full-stops there has to be a switch at
        ! some point of the double overlap region or else it is single-excitation
        ! like... so only consider x1 matrix element here..

        ! combine deltaB and stepvalue info here to reduce if statements
        ! asserts dont work anymore with new select case statements
        ! do it out here:
#ifdef DEBUG_
        if (csf_i%stepvector(ende) == 1) then
            ASSERT(deltaB /= 2)
        end if
        if (csf_i%stepvector(ende) == 2) then
            ASSERT(deltaB /= -2)
        end if
#endif

        select case (deltaB + csf_i%stepvector(ende))
        case (1)
            ! d=1 + b=0 : 1
            ! ! +2 branch not allowed here
            ! not sure if i can access only the x1 element down there..
            call getMixedFullStop(1, 1, 0, bVal, tempWeight_0, tempWeight_1)

        case (-1)
            ! d=1 + b=-2 : -1
            ! deltaB = -2
            ! switch 1 -> 2
            set_orb(t, 2 * ende)
            clr_orb(t, 2 * ende - 1)

            ! matrix element is 1 in this case
            tempWeight_1 = 1.0_dp
            tempWeight_0 = 0.0_dp

        case (2)
            ! d=2 + b=0 : 2

            call getMixedFullStop(2, 2, 0, bVal, tempWeight_0, tempWeight_1)

        case (4)
            ! d=2 + b=2 : 4
            ! deltab = 2
            ! switch 2 -> 1
            set_orb(t, 2 * ende - 1)
            clr_orb(t, 2 * ende)

            tempWeight_1 = 1.0_dp
            tempWeight_0 = 0.0_dp
        end select

        ! thats kind of stupid what i did here...
        ! check if x0-matrix element is non-zero which is an indicator that
        ! no switch happened in the double overlap region

        ! just for completion reasons still update the x0 matrix element here
        ! although it should be 0 anyway..
        call update_matrix_element(t, tempWeight_0, 1)
        call update_matrix_element(t, tempWeight_1, 2)

        ! todo... have to write some excitation cancellation function...
        ! to get back to the start of an excitation or somehow ensure that
        ! switch happens...
    end subroutine mixedFullStopStochastic

    subroutine calcRaisingSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                              posSwitches, t, probWeight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: negSwitches(nSpatOrbs), posSwitches(nSpatOrbs)
        integer(n_int), intent(inout) :: t(0:nifguga)
        real(dp), intent(inout) :: probWeight
        character(*), parameter :: this_routine = "calcRaisingSemiStartStochastic"

        integer :: se, deltaB
        real(dp) :: tempWeight_0, tempWeight_1, minusWeight, &
                    plusWeight, zeroWeight, bVal

        ASSERT(.not. isThree(ilut, excitInfo%secondStart))

        ! i can be sure that there is no 3 or 0 at the fullEnd, or otherwise
        ! this would be single-excitation like.
        se = excitInfo%secondStart
        bVal = csf_i%B_real(se)

        deltaB = getDeltaB(t)

        ! do non-choosing possibs first
        ! why does this cause a segfault on compilation with gfortran??
        ! do some debugging:

        ! fix for gfortran compilation for some reasono
        ! i can probably fix it when i finally get to this point in
        ! test running

        select case (csf_i%stepvector(se))
        case (1)
            ! 1 -> 3
            set_orb(t, 2 * se)

            call setDeltaB(deltaB + 1, t)

            call getDoubleMatrixElement(3, 1, deltaB, excitInfo%gen1, &
                                        excitInfo%gen2, bVal, excitInfo%order, &
                                        tempWeight_0, tempWeight_1)

        case (2)
            ! 2 -> 3
            set_orb(t, 2 * se - 1)

            call setDeltaB(deltaB - 1, t)

            call getDoubleMatrixElement(3, 2, deltaB, excitInfo%gen1, &
                                        excitInfo%gen2, bVal, excitInfo%order, &
                                        tempWeight_0, tempWeight_1)

        case (0)
            ! 0:
            ! for arriving -1 branch branching is always possible
            if (deltaB == -1) then
                ! here the choice is between 0 and -2 branch
                minusWeight = weights%proc%minus(negSwitches(se), bVal, weights%dat)
                zeroWeight = weights%proc%zero(negSwitches(se), posSwitches(se), &
                                               bVal, weights%dat)

                if (near_zero(minusWeight + zeroWeight)) then
                    probWeight = 0.0_dp
                    t = 0
                    return
                end if

                ! cant directly cant assign it to probWeight
                zeroWeight = calcStartProb(zeroWeight, minusWeight)

                if (genrand_real2_dSFMT() < zeroWeight) then
                    ! go to 0 branch
                    ! 0 -> 2
                    set_orb(t, 2 * se)

                    call setDeltaB(0, t)

                    call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, excitInfo%order, &
                                                tempWeight_0, tempWeight_1)

                    probWeight = probWeight * zeroWeight

                else
                    ! to to -2 branch
                    ! 0 -> 1
                    set_orb(t, 2 * se - 1)

                    call setDeltaB(-2, t)

                    call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, excitInfo%order, &
                                                x1_element=tempWeight_1)

                    tempWeight_0 = 0.0_dp

                    probWeight = probWeight * (1.0_dp - zeroWeight)
                end if
            else
                ! +1 branch arriving -> have to check b values
                ! UPDATE: include b value check into probWeight calculation
                ! todo
                if (csf_i%B_int(se) < 2) then
                    ! only 0 branch possible
                    ! todo: in this forced cases due to the b value, have to
                    ! think about, how that influences the probWeight...
                    ! 0 -> 1
                    set_orb(t, 2 * se - 1)

                    call setDeltaB(0, t)

                    call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, excitInfo%order, &
                                                tempWeight_0, tempWeight_1)

                else
                    ! need probs
                    plusWeight = weights%proc%plus(posSwitches(se), &
                                                   bVal, weights%dat)
                    zeroWeight = weights%proc%zero(negSwitches(se), posSwitches(se), &
                                                   bVal, weights%dat)

                    if (near_zero(plusWeight + zeroWeight)) then
                        probWeight = 0.0_dp
                        t = 0
                        return
                    end if

                    ! cant directly cant assign it to probWeight
                    zeroWeight = calcStartProb(zeroWeight, plusWeight)

                    if (genrand_real2_dSFMT() < zeroWeight) then
                        ! do 0 branch
                        ! 0 -> 1
                        set_orb(t, 2 * se - 1)

                        call setDeltaB(0, t)

                        call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, &
                                                    tempWeight_0, tempWeight_1)

                        probWeight = probWeight * zeroWeight

                    else
                        ! do +2 branch
                        ! 0 -> 2
                        set_orb(t, 2 * se)

                        call setDeltaB(2, t)

                        call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, &
                                                    x1_element=tempWeight_1)

                        tempWeight_0 = 0.0_dp

                        probWeight = probWeight * (1.0_dp - zeroWeight)
                    end if
                end if
            end if
        end select

        call encode_matrix_element(t, extract_matrix_element(t, 1) * tempWeight_1, 2)
        call update_matrix_element(t, tempWeight_0, 1)

        if (near_zero(tempWeight_0) .and. near_zero(tempWeight_1)) then
            probWeight = 0.0_dp
            t = 0
        end if

    end subroutine calcRaisingSemiStartStochastic

    subroutine calcLoweringSemiStartStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                               posSwitches, t, probWeight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: negSwitches(nSpatOrbs), posSwitches(nSpatOrbs)
        integer(n_int), intent(inout) :: t(0:nifguga)
        real(dp), intent(inout) :: probWeight
        character(*), parameter :: this_routine = "calcLoweringSemiStartStochastic"

        integer :: se, deltaB
        real(dp) :: tempWeight_0, tempWeight_1, minusWeight, &
                    plusWeight, zeroWeight, bVal

        ASSERT(.not. isZero(ilut, excitInfo%secondStart))

        ! i can be sure that there is no 3 or 0 at the fullEnd, or otherwise
        ! this would be single-excitation like.

        se = excitInfo%secondStart
        bVal = csf_i%B_real(se)

        deltaB = getDeltaB(t)

        ! do non-choosing possibs first
        ! why does this cause a segfault on compilation with gfortran??
        ! do some debugging:
        ! same gfortran compilex issue fix as above

        select case (csf_i%stepvector(se))
        case (1)
            ! 1 -> 0
            clr_orb(t, 2 * se - 1)

            call setDeltaB(deltaB + 1, t)

            call getDoubleMatrixElement(0, 1, deltaB, excitInfo%gen1, &
                                        excitInfo%gen2, bVal, excitInfo%order, &
                                        tempWeight_0, tempWeight_1)

        case (2)
            ! 2 -> 0
            clr_orb(t, 2 * se)

            call setDeltaB(deltaB - 1, t)

            call getDoubleMatrixElement(0, 2, deltaB, excitInfo%gen1, &
                                        excitInfo%gen2, bVal, excitInfo%order, &
                                        tempWeight_0, tempWeight_1)

        case (3)
            ! 3:
            ! for arriving -1 branch branching is always possible
            if (deltaB == -1) then
                ! here the choice is between 0 and -2 branch
                minusWeight = weights%proc%minus(negSwitches(se), &
                                                 bVal, weights%dat)
                zeroWeight = weights%proc%zero(negSwitches(se), posSwitches(se), &
                                               bVal, weights%dat)

                if (near_zero(minusWeight + zeroWeight)) then
                    probWeight = 0.0_dp
                    t = 0
                    return
                end if
                ! cant directly cant assign it to probWeight
                zeroWeight = calcStartProb(zeroWeight, minusWeight)

                if (genrand_real2_dSFMT() < zeroWeight) then
                    ! go to 0 branch
                    ! 3 -> 2
                    clr_orb(t, 2 * se - 1)

                    call setDeltaB(0, t)

                    call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, excitInfo%order, &
                                                tempWeight_0, tempWeight_1)

                    probWeight = probWeight * zeroWeight

                else
                    ! to to -2 branch
                    ! 3 -> 1
                    clr_orb(t, 2 * se)

                    call setDeltaB(-2, t)

                    call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, excitInfo%order, &
                                                x1_element=tempWeight_1)

                    tempWeight_0 = 0.0_dp

                    probWeight = probWeight * (1.0_dp - zeroWeight)
                end if
            else
                ! +1 branch arriving -> have to check b values
                if (csf_i%B_int(se) < 2) then
                    ! only 0 branch possible
                    ! todo: in this forced cases due to the b value, have to
                    ! think about, how that influences the probWeight...
                    ! 3 -> 1
                    clr_orb(t, 2 * se)

                    call setDeltaB(0, t)

                    call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, excitInfo%order, &
                                                tempWeight_0, tempWeight_1)

                else
                    ! need probs
                    plusWeight = weights%proc%plus(posSwitches(se), &
                                                   bVal, weights%dat)
                    zeroWeight = weights%proc%zero(negSwitches(se), posSwitches(se), &
                                                   bVal, weights%dat)

                    if (near_zero(plusWeight + zeroWeight)) then
                        probWeight = 0.0_dp
                        t = 0
                        return
                    end if

                    ! cant directly cant assign it to probWeight
                    zeroWeight = calcStartProb(zeroWeight, plusWeight)

                    if (genrand_real2_dSFMT() < zeroWeight) then
                        ! do 0 branch
                        ! 3 -> 1
                        clr_orb(t, 2 * se)

                        call setDeltaB(0, t)

                        call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, &
                                                    tempWeight_0, tempWeight_1)

                        probWeight = probWeight * zeroWeight

                    else
                        ! do +2 branch
                        ! 3 -> 2
                        clr_orb(t, 2 * se - 1)

                        call setDeltaB(2, t)

                        call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, &
                                                    x1_element=tempWeight_1)

                        tempWeight_0 = 0.0_dp

                        probWeight = probWeight * (1.0_dp - zeroWeight)
                    end if
                end if
            end if
        end select

        call encode_matrix_element(t, extract_matrix_element(t, 1) * tempWeight_1, 2)
        call update_matrix_element(t, tempWeight_0, 1)

        if (near_zero(tempWeight_0) .and. near_zero(tempWeight_1)) then
            probWeight = 0.0_dp
            t = 0
        end if

    end subroutine calcLoweringSemiStartStochastic

    subroutine calcRaisingSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                             posSwitches, t, probWeight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: negSwitches(nSpatOrbs), posSwitches(nSpatOrbs)
        integer(n_int), intent(inout) :: t(0:nifguga)
        real(dp), intent(inout) :: probWeight
        character(*), parameter :: this_routine = "calcRaisingSemiStopStochastic"

        integer :: semi, deltaB
        real(dp) :: tempWeight_0, tempWeight_1, minusWeight, &
                    plusWeight, bVal

        ASSERT(.not. isZero(ilut, excitInfo%firstEnd))

        semi = excitInfo%firstEnd
        bVal = csf_i%B_real(semi)

        ! in the stochastic case i am sure that at there is no 3 or 0 at the
        ! full start... so definetly all deltaB branches can arrive here
        ! first deal with the forced ones
        deltaB = getDeltaB(t)

        select case (csf_i%stepvector(semi))
        case (1)
            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                if (getDeltaB(t) == 2) then
                    t = 0_n_int
                    probWeight = 0.0_dp
                    return
                end if
            end if

            ASSERT(getDeltaB(t) /= 2)

            ! do the 1 -> 0 switch
            clr_orb(t, 2 * semi - 1)

            call setDeltaB(deltaB + 1, t)

            call getDoubleMatrixElement(0, 1, deltaB, excitInfo%gen1, &
                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                        tempWeight_0, tempWeight_1)

        case (2)
            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                if (getDeltaB(t) == -2) then
                    t = 0_n_int
                    probWeight = 0.0_dp
                    return
                end if
            end if

            ASSERT(getDeltaB(t) /= -2)

            ! do 2 -> 0
            clr_orb(t, 2 * semi)

            call setDeltaB(deltaB - 1, t)

            call getDoubleMatrixElement(0, 2, deltaB, excitInfo%gen1, &
                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                        tempWeight_0, tempWeight_1)

        case (3)
            ! its a 3 -> have to check b if a 0 branch arrives
            select case (deltaB)
            case (-2)
                ! do 3 -> 2
                clr_orb(t, 2 * semi - 1)

                call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                            excitInfo%gen2, bVal, excitInfo%order1, &
                                            x1_element=tempWeight_1)

                tempWeight_0 = 0.0_dp

                call setDeltaB(-1, t)

            case (2)
                ! do 3 -> 1
                clr_orb(t, 2 * semi)

                call setDeltaB(+1, t)

                call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                            excitInfo%gen2, bVal, excitInfo%order1, &
                                            x1_element=tempWeight_1)

                tempWeight_0 = 0.0_dp

            case (0)
                ! deltaB = 0 branch arrives -> check b
                if (csf_i%B_int(semi) == 0) then
                    ! only 3 -> 1 possble
                    clr_orb(t, 2 * semi)

                    call setDeltaB(-1, t)

                    call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, excitInfo%order1, &
                                                tempWeight_0, tempWeight_1)

                else
                    ! finally have to check the probs
                    minusWeight = weights%proc%minus(negSwitches(semi), &
                                                     bVal, weights%dat)
                    plusWeight = weights%proc%plus(posSwitches(semi), &
                                                   bVal, weights%dat)

                    if (near_zero(minusWeight + plusWeight)) then
                        probWeight = 0.0_dp
                        t = 0
                        return
                    end if

                    ! here i cant directly save it to probWeight ...
                    minusWeight = calcStartProb(minusWeight, plusWeight)

                    if (genrand_real2_dSFMT() < minusWeight) then
                        ! do -1 branch:
                        ! set 3 -> 1
                        clr_orb(t, 2 * semi)

                        call setDeltaB(-1, t)

                        call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order1, &
                                                    tempWeight_0, tempWeight_1)

                        probWeight = probWeight * minusWeight

                    else
                        ! do +1 branch
                        ! set 3 -> 2
                        clr_orb(t, 2 * semi - 1)

                        call setDeltaB(1, t)

                        call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order1, &
                                                    tempWeight_0, tempWeight_1)

                        probWeight = probWeight * (1 - minusWeight)
                    end if
                end if
            end select
        end select

        ! if its a mixed full-start raising into lowering -> check i a
        ! switch happened in the double overlap region
        if (excitInfo%typ == excit_type%fullstart_R_to_L) then
            ! this is indicated by a non-zero x0-matrix element
            if (.not. near_zero(extract_matrix_element(t, 1) * tempWeight_0)) then
                probWeight = 0.0_dp
                t = 0
                return
            end if
        end if

        ! do not combine the x0 and x1 elements at this point, because its
        ! easier to deal with the contributing integrals if we keep them
        ! seperate until the end!

        if (near_zero(extract_matrix_element(t, 1) * tempWeight_0) .and. &
            near_zero(extract_matrix_element(t, 2) * tempWeight_1)) then
            probWeight = 0.0_dp
            t = 0_n_int
            return
        end if

        call update_matrix_element(t, tempWeight_0, 1)
        call update_matrix_element(t, tempWeight_1, 2)

    end subroutine calcRaisingSemiStopStochastic

    subroutine calcLoweringSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                              posSwitches, t, probWeight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: negSwitches(nSpatOrbs), posSwitches(nSpatOrbs)
        integer(n_int), intent(inout) :: t(0:nifguga)
        real(dp), intent(inout) :: probWeight
        character(*), parameter :: this_routine = "calcLoweringSemiStopStochastic"

        integer :: semi, deltaB
        real(dp) :: tempWeight_0, tempWeight_1, minusWeight, &
                    plusWeight, bVal

        ASSERT(.not. isThree(ilut, excitInfo%firstEnd))

        semi = excitInfo%firstEnd
        bVal = csf_i%B_real(semi)

        ! in the stochastic case i am sure that at there is no 3 or 0 at the
        ! full start... so definetly all deltaB branches can arrive here
        ! first deal with the forced ones
        deltaB = getDeltaB(t)

        select case (csf_i%stepvector(semi))
        case (1)
            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                if (getDeltaB(t) == 2) then
                    t = 0
                    return
                end if
            end if

            ASSERT(getDeltaB(t) /= 2)

            ! do the 1 -> 3 switch
            set_orb(t, 2 * semi)

            call setDeltaB(deltaB + 1, t)

            call getDoubleMatrixElement(3, 1, deltaB, excitInfo%gen1, &
                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                        tempWeight_0, tempWeight_1)

        case (2)
            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                if (getDeltaB(t) == -2) then
                    t = 0
                    return
                end if
            end if

            ASSERT(getDeltaB(t) /= -2)

            ! do 2 -> 3
            set_orb(t, 2 * semi - 1)

            call setDeltaB(deltaB - 1, t)

            call getDoubleMatrixElement(3, 2, deltaB, excitInfo%gen1, &
                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                        tempWeight_0, tempWeight_1)

        case (0)
            ! its a 0 -> have to check b if a 0 branch arrives
            select case (deltaB)
            case (-2)
                ! do 0 -> 2
                set_orb(t, 2 * semi)

                call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                            excitInfo%gen2, bVal, excitInfo%order1, &
                                            x1_element=tempWeight_1)

                tempWeight_0 = 0.0_dp

                call setDeltaB(-1, t)

            case (2)
                ! do 0 -> 1
                set_orb(t, 2 * semi - 1)

                call setDeltaB(+1, t)

                call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                            excitInfo%gen2, bVal, excitInfo%order1, &
                                            x1_element=tempWeight_1)

                tempWeight_0 = 0.0_dp

            case (0)
                ! deltaB=0 branch arrives -> check b
                if (csf_i%B_int(semi) == 0) then
                    ! only 0 -> 1 possble
                    set_orb(t, 2 * semi - 1)

                    call setDeltaB(-1, t)

                    call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, excitInfo%order1, &
                                                tempWeight_0, tempWeight_1)

                else
                    ! finally have to check the probs
                    minusWeight = weights%proc%minus(negSwitches(semi), &
                                                     bVal, weights%dat)
                    plusWeight = weights%proc%plus(posSwitches(semi), &
                                                   bVal, weights%dat)

                    if (near_zero(minusWeight + plusWeight)) then
                        probWeight = 0.0_dp
                        t = 0
                        return
                    end if

                    ! here i cant directly save it to probWeight ...
                    minusWeight = calcStartProb(minusWeight, plusWeight)

                    if (genrand_real2_dSFMT() < minusWeight) then
                        ! do -1 branch:
                        ! set 0 -> 1
                        set_orb(t, 2 * semi - 1)

                        call setDeltaB(-1, t)

                        call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order1, &
                                                    tempWeight_0, tempWeight_1)

                        probWeight = probWeight * minusWeight

                    else
                        ! do +1 branch
                        ! set 0 -> 2
                        set_orb(t, 2 * semi)

                        call setDeltaB(1, t)

                        call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order1, &
                                                    tempWeight_0, tempWeight_1)

                        probWeight = probWeight * (1.0_dp - minusWeight)
                    end if
                end if
            end select
        end select

        ! for mixed fullstart check if no switch happened in the double
        ! overlap region, indicated by a non-zero-x0 matrix element
        if (excitInfo%typ == excit_type%fullStart_L_to_R) then
            if (.not. near_zero(extract_matrix_element(t, 1) * tempWeight_0)) then
                probWeight = 0.0_dp
                t = 0
                return
            end if
        end if

        if (near_zero(extract_matrix_element(t, 1) * tempWeight_0) .and. &
            near_zero(extract_matrix_element(t, 2) * tempWeight_1)) then
            probWeight = 0.0_dp
            t = 0
            return
        end if

        call update_matrix_element(t, tempWeight_0, 1)
        call update_matrix_element(t, tempWeight_1, 2)

    end subroutine calcLoweringSemiStopStochastic

    subroutine calcFullStartR2L_stochastic(ilut, csf_i, excitInfo, t, pgen, &
                                           posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcFullStartR2L_stochastic"

        integer :: i, st, en, se, gen, j, k, l, holeInd, elecInd
        type(WeightObj_t) :: weights
        real(dp) :: branch_pgen, temp_pgen, &
                    p_orig, rdm_mat, orb_pgen
        HElement_t(dp) :: integral
        procedure(calc_pgen_general), pointer :: calc_pgen_yix_start

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(.not. isThree(ilut, excitInfo%fullStart))

        ! create the fullStart
        st = excitInfo%fullStart
        en = excitInfo%fullEnd
        se = excitInfo%firstEnd
        gen = excitInfo%lastGen

        ! create correct weights:
        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_fullStartWeight(csf_i, se, en, negSwitches(se), &
                                           posSwitches(se), csf_i%B_real(se))
        end if

        if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
            ! todo the switch
            call forced_mixed_start(ilut, csf_i, excitInfo, t, branch_pgen)
        else
            call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)
        end if

        ! set pgen to 0 in case of early exit
        pgen = 0.0_dp

        ! check validity
        check_abort_excit(branch_pgen, t)

        ! then for the overlap region i need a double update routine, which
        ! somehow garuantees a switch happens at some point, to avoid
        ! single excitations

        do i = st + 1, se - 1
            call doubleUpdateStochastic(ilut, csf_i, i, excitInfo, &
                                        weights, negSwitches, posSwitches, t, branch_pgen)

            ! check validity
            if (near_zero(extract_matrix_element(t, 2)) .or. near_zero(branch_pgen)) then
                t = 0_n_int
                return
            end if
        end do

        ! then deal with specific semi-stop
        ! and update weights here
        ! i also could use the fact that the single weights are already
        ! initialized within the fullstart weights or?
        ! and then use smth like
        weights = weights%ptr

        call calcRaisingSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                           posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        excitInfo%currentGen = excitInfo%lastGen

        do i = se + 1, en - 1
            call singleStochasticUpdate(ilut, csf_i, i, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            branch_pgen = branch_pgen * temp_pgen

            if (t_trunc_guga_pgen .or. &
                (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                if (branch_pgen < trunc_guga_pgen) then
                    t = 0_n_int
                    return
                end if
            end if

            ! check validity
            check_abort_excit(branch_pgen, t)

        end do

        call singleStochasticEnd(csf_i, excitInfo, t)

        ! if we do RDMs also store the x0 and x1 coupling coeffs
        ! and I need to do it before the routines below since excitInfo
        ! gets changed there
        if (tFillingStochRDMOnFly) then
            ! i need to unbias against the total pgen later on in the
            ! RDM sampling otherwise the rdm-bias factor is not correct!
            ! encode the necessary information in the rdm-matele!
            i = excitInfo%i
            j = excitInfo%j
            k = excitInfo%k
            l = excitInfo%l
            elecInd = se
            holeInd = en
            rdm_mat = extract_matrix_element(t, 2)
            call calc_orbital_pgen_contrib_start(&
                    csf_i, [2 * st, 2 * elecInd], holeInd, orb_pgen)
            p_orig = orb_pgen * branch_pgen / real(ElecPairs, dp)
            if (csf_i%stepvector(elecInd) == 3) p_orig = p_orig * 2.0_dp
        end if

        call encode_matrix_element(t, extract_matrix_element(t, 1) + &
                                   extract_matrix_element(t, 2), 1)
        call encode_matrix_element(t, 0.0_dp, 2)

        ! update: since i can formulate everything in terms of the already
        ! calculated matrix element i can abort here if it is zero
        if (near_zero(extract_matrix_element(t, 1))) then
            t = 0_n_int
            return
        end if

        global_excitInfo = excitInfo

        if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
            call calc_mixed_start_contr_approx(t, csf_i, excitInfo, integral)
            pgen = branch_pgen

        else
            call calc_mixed_start_contr_sym(ilut, csf_i, t, excitInfo, branch_pgen, pgen, &
                                            integral)
        end if

        if (tFillingStochRDMOnFly) then
            if (.not. near_zero(p_orig)) then
                call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                                    contract_2_rdm_ind(i, j, k, l, excit_lvl=2, &
                                       excit_typ=excitInfo%typ), x0=0.0_dp, &
                                                x1=rdm_mat * pgen / p_orig)
            end if
        end if

        ! and finally update the matrix element with all contributions
        call update_matrix_element(t, (get_umat_el(st, en, se, st) + &
                                       get_umat_el(en, st, st, se)) / 2.0_dp + integral, 1)

    end subroutine calcFullStartR2L_stochastic

    subroutine calc_mixed_start_contr_approx(t, csf_i, excitInfo, integral)
        integer(n_int), intent(in) :: t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        HElement_t(dp), intent(out) :: integral
        character(*), parameter :: this_routine = "calc_mixed_start_contr_approx"

        integer :: st, se, en, elecInd, holeInd, sw, step, i
        real(dp) :: bot_cont, mat_ele, stay_mat, start_mat
        logical :: below_flag

        st = excitInfo%fullStart
        se = excitInfo%firstEnd
        en = excitInfo%fullEnd
        ! depending on the type of excitaiton, calculation of orbital pgens
        ! change
        if (excitInfo%typ == excit_type%fullStart_L_to_R) then
            elecInd = en
            holeInd = se
        else if (excitInfo%typ == excit_type%fullstart_R_to_L) then
            elecInd = se
            holeInd = en
        else
            call stop_all(this_routine, "should not be here!")
        end if

        ! I am sure first switch is at full-start
        sw = st

        ! what can i precalculate beforehand?
        step = csf_i%stepvector(st)

        integral = h_cast(0.0_dp)

        if (step == 1) then

            ASSERT(isTwo(t, st))

            bot_cont = -sqrt(2.0_dp / ((csf_i%B_real(st) - 1.0_dp) * &
                                       (csf_i%B_real(st) + 1.0_dp)))

        else

            ASSERT(isOne(t, st))

            bot_cont = -sqrt(2.0_dp / ((csf_i%B_real(st) + 1.0_dp) * &
                                       (csf_i%B_real(st) + 3.0_dp)))

        end if

        if (.not. near_zero(bot_cont)) then

            mat_ele = 1.0_dp
            below_flag = .false.

            do i = st - 1, 1, -1
                if (csf_i%Occ_int(i) /= 1) cycle

                ! then check if thats the last stepvalue to consider
                if (csf_i%stepvector(i) == 1 .and. csf_i%B_int(i) == 1) then
                    below_flag = .true.
                end if

                ! then deal with the matrix element and branching probabilities
                step = csf_i%stepvector(i)

                ! get both start and staying matrix elements -> and update
                ! matrix element contributions on the fly to avoid second loop!
                call getDoubleMatrixElement(step, step, -1, gen_type%R, gen_type%L, csf_i%B_real(i), &
                                            1.0_dp, x1_element=start_mat)

                call getDoubleMatrixElement(step, step, 0, gen_type%R, gen_type%L, csf_i%B_real(i), &
                                            1.0_dp, x1_element=stay_mat)

                if (near_zero(stay_mat)) below_flag = .true.
                ! "normally" matrix element shouldnt be 0 anymore... still check

                if (.not. near_zero(start_mat)) then
                    integral = integral + start_mat * mat_ele * (get_umat_el(i, holeInd, elecInd, i) &
                                                                 + get_umat_el(holeInd, i, i, elecInd)) / 2.0_dp
                end if

                ! also update matrix element on the fly
                mat_ele = stay_mat * mat_ele

            end do
        end if

        ! maybe I need to deal with the pgens here, if they are not
        ! correctly considered outside..

    end subroutine calc_mixed_start_contr_approx

    subroutine calcFullStartL2R_stochastic(ilut, csf_i, excitInfo, t, pgen, &
                                           posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcFullStartL2R_stochastic"

        integer :: i, st, en, se, gen, j, k, l, elecInd, holeInd
        type(WeightObj_t) :: weights
        real(dp) :: branch_pgen, temp_pgen, p_orig, orb_pgen, rdm_mat
        HElement_t(dp) :: integral
        procedure(calc_pgen_general), pointer :: calc_pgen_yix_start

        ! the integral contributions to this mixed excitaiton full starts are
        ! more involved then i have imagined. since all the orbitals from below
        ! the first stepvector change(which has to happen at some point of the
        ! excitation, or else it is single-like,(which is also still unclear
        ! how i should implement that)) contribute to the same excitation.
        ! similar to the full-start full-stop mixed generator case.
        ! to ensure a stepvector change i could write some kind of restricted
        ! double excitiation update function, which accumulates switch
        ! probability with decreasing switch possibilities...
        ! and then during the excitation creation i have to save the first
        ! stepvector change and depending on that calculate the remaining
        ! integral contribution. todo
        ! and same with the fullstartr2l, and fullstoprl2/l2r
        ! for both full starts and full stop routines, the contributions are
        ! the same atleast.
        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(.not. isThree(ilut, excitInfo%fullStart))

        ! create the fullStart
        st = excitInfo%fullStart
        en = excitInfo%fullEnd
        se = excitInfo%firstEnd
        gen = excitInfo%lastGen

        ! create correct weights:
        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_fullStartWeight(csf_i, se, en, negSwitches(se), &
                                           posSwitches(se), csf_i%B_real(se))
        end if

        ! in the case of the approximate exchange excitations I need to
        ! force a switch at the beginning
        if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
            ! do the switch
            call forced_mixed_start(ilut, csf_i, excitInfo, t, branch_pgen)

        else
            call mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)
        end if

        ! in case of early exit pgen should be set to 0
        pgen = 0.0_dp

        check_abort_excit(branch_pgen, t)

        ! in the mixed fullstart case there has to be a switch at some point
        ! in the double overlap region, or else it would be a single-like
        ! excitation -> so check if x1 matrix element is 0

        ! then for the overlap region i need a double update routine, which
        ! somehow garuantees a switch happens at some point, to avoid
        ! single excitations
        do i = st + 1, se - 1
            call doubleUpdateStochastic(ilut, csf_i, i, excitInfo, &
                                        weights, negSwitches, posSwitches, t, branch_pgen)

            ! to keep it general, i cant only check weights in doubleUpdate
            if (near_zero(extract_matrix_element(t, 2)) .or. near_zero(branch_pgen)) then
                t = 0_n_int
                return
            end if
        end do

        ! then deal with specific semi-stop
        ! and update weights here
        weights = weights%ptr

        call calcLoweringSemiStopStochastic(ilut, csf_i, excitInfo, weights, negSwitches, &
                                            posSwitches, t, branch_pgen)

        ! check validity
        check_abort_excit(branch_pgen, t)

        excitInfo%currentGen = excitInfo%lastGen

        do i = se + 1, en - 1
            call singleStochasticUpdate(ilut, csf_i, i, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            branch_pgen = branch_pgen * temp_pgen

            if (t_trunc_guga_pgen .or. &
                (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                if (branch_pgen < trunc_guga_pgen) then
                    t = 0_n_int
                    return
                end if
            end if

            ! check validity
            check_abort_excit(branch_pgen, t)
        end do

        call singleStochasticEnd(csf_i, excitInfo, t)

        ! if we do RDMs also store the x0 and x1 coupling coeffs
        ! and I need to do it before the routines below since excitInfo
        ! gets changed there
        if (tFillingStochRDMOnFly) then
            ! i need to unbias against the total pgen later on in the
            ! RDM sampling otherwise the rdm-bias factor is not correct!
            ! encode the necessary information in the rdm-matele!
            i = excitInfo%i
            j = excitInfo%j
            k = excitInfo%k
            l = excitInfo%l
            elecInd = en
            holeInd = se
            rdm_mat = extract_matrix_element(t, 2)
            call calc_orbital_pgen_contrib_start(&
                    csf_i, [2 * st, 2 * elecInd], holeInd, orb_pgen)
            p_orig = orb_pgen * branch_pgen / real(ElecPairs, dp)
            if (csf_i%stepvector(elecInd) == 3) p_orig = p_orig * 2.0_dp
        end if

        ! put everything in first entry
        call encode_matrix_element(t, extract_matrix_element(t, 1) + &
                                   extract_matrix_element(t, 2), 1)
        call encode_matrix_element(t, 0.0_dp, 2)

        ! now have to check if there has be a change in the double overlap
        ! region, or else its just a single-like excitation:
        ! todo write a routine for that
        ! switchFlag = checkForSwitch(ilut, t, st, se-1)

        ! if a switch happened i have to calculated the additional contributing
        ! matrix elements from all indices below the fullstart
        ! can i formulate that in terms of the already calulated matrix
        ! element, as in the diagonal matrix element calculation?
        ! probably... but todo!

        ! update: since i can formulate everything in terms of the already
        ! calculated matrix element i can abort here if it is zero
        if (near_zero(extract_matrix_element(t, 1))) then
            t = 0_n_int
            return
        end if

        global_excitInfo = excitInfo

        if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
            call calc_mixed_start_contr_approx(t, csf_i, excitInfo, integral)
            pgen = branch_pgen
        else
            call calc_mixed_start_contr_sym(ilut, csf_i, t, excitInfo, branch_pgen, pgen, integral)
        end if

        if (tFillingStochRDMOnFly) then
            if (.not. near_zero(p_orig)) then
                call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                            contract_2_rdm_ind(i, j, k, l, excit_lvl=2, &
                                       excit_typ=excitInfo%typ), x0=0.0_dp, &
                                                x1=rdm_mat * pgen / p_orig)
            end if
        end if

        ! and finally update the matrix element with all contributions
        call update_matrix_element(t, (get_umat_el(st, se, en, st) + &
                                       get_umat_el(se, st, st, en)) / 2.0_dp + integral, 1)

    end subroutine calcFullStartL2R_stochastic

    subroutine calc_mixed_x2x_ueg(ilut, csf_i, t, excitInfo, branch_pgen, pgen, &
                                  integral, rdm_ind, rdm_mat)
        integer(n_int), intent(in) :: ilut(0:nifguga), t(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        real(dp), intent(inout) :: branch_pgen
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: integral
        integer(int_rdm), intent(out), allocatable, optional :: rdm_ind(:)
        real(dp), intent(out), allocatable, optional :: rdm_mat(:)
        character(*), parameter :: this_routine = "calc_mixed_x2x_ueg"

        pgen = 0.0_dp
        integral = 0.0_dp
        unused_var(ilut); unused_var(t); unused_var(excitInfo); unused_var(branch_pgen);
        unused_var(csf_i)
        if (present(rdm_ind)) then
            allocate(rdm_ind(0), source=0_int_rdm)
        end if
        if (present(rdm_mat)) then
            allocate(rdm_mat(0), source=0.0_dp)
        end if
        call stop_all(this_routine, &
                      "in Hubbard/UEG calculations with full k-point symmetry, this excitation shouldnt be reached!")
    end subroutine calc_mixed_x2x_ueg

    subroutine calc_orbital_pgen_contrib_end_def(csf_i, occ_orbs, orb_a, orb_pgen)
        ! write a combined function for both r2l and l2r since its only
        ! one difference -> only one if condition to adjust for both!
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen

        integer :: i, j, orb
        real(dp) :: cum_sum, cpt_a, cpt_b, cpt_ba, cpt_ab, ba_sum, ab_sum

        ! electron indices
        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))

        cum_sum = 0.0_dp
        if (tGen_guga_weighted) then
            do orb = 1, orb_a - 1
                ! calc. the p(a)
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)

                end if

            end do
        end if

        ! deal with orb (a) in specific way:
        cpt_a = get_guga_integral_contrib(occ_orbs, orb_a, -1)

        cum_sum = cum_sum + cpt_a

        ! also get p(b|a)
        ! depending if its a r2l or l2r full-stop:
        if (i < orb_a) then
            ! its a L2R -> so no restrictions
            call pgen_select_orb_guga_mol(csf_i, occ_orbs, orb_a, j, cpt_ba, ba_sum)
        else
            ! its a R2L so orbital i is off-limits
            call pgen_select_orb_guga_mol(csf_i, occ_orbs, orb_a, j, cpt_ba, ba_sum, i)
        end if

        ! change to the fullstart into fullstop: loop until orbital j for the
        ! fullstop implementattion
        if (tGen_guga_weighted) then
            do orb = orb_a + 1, j - 1
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)
                end if
            end do
        end if

        ! deal with j also speciallly
        cpt_b = get_guga_integral_contrib(occ_orbs, j, -1)

        cum_sum = cum_sum + cpt_b

        ! and get p(a|b)
        ! only orbitals below j are allowed!
        call pgen_select_orb_guga_mol(csf_i, occ_orbs, j, orb_a, cpt_ab, ab_sum, -j, .true.)

        ! and deal with rest:
        if (tGen_guga_weighted) then
            do orb = j + 1, nSpatOrbs
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)
                end if
            end do
        end if

        if (.not. tGen_guga_weighted) then
            cum_sum = csf_i%cum_list(nSpatOrbs)
        end if

        if (near_zero(cum_sum) .or. near_zero(ab_sum) .or. near_zero(ba_sum)) then
            orb_pgen = 0.0_dp
        else
            ! and get hopefully correct final values:
            cpt_a = cpt_a / cum_sum * cpt_ba / ba_sum
            cpt_b = cpt_b / cum_sum * cpt_ab / ab_sum

            ! and add them up to the final orbital pgen
            orb_pgen = cpt_a + cpt_b
        end if

    end subroutine calc_orbital_pgen_contrib_end_def

    subroutine calc_orbital_pgen_contrib_start_def(csf_i, occ_orbs, orb_a, orb_pgen)
        ! write a combined function for both r2l and l2r since its only
        ! one difference -> only one if condition to adjust for both!
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a
        real(dp), intent(out) :: orb_pgen

        integer :: i, j, orb
        real(dp) :: cum_sum, cpt_a, cpt_b, cpt_ba, cpt_ab, ba_sum, ab_sum

        ! electron indices
        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))

        ! damn... i need the probability of the elec-pair picking too or?

        cum_sum = 0.0_dp
        if (tGen_guga_weighted) then
            do orb = 1, i - 1
                ! calc. the p(a)
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)
                end if

            end do
        end if

        ! deal with orb (i) in specific way:
        cpt_a = get_guga_integral_contrib(occ_orbs, i, -1)

        cum_sum = cum_sum + cpt_a

        ! also get p(b|a)
        ! did i mix up i and j here below.. what do i assume picked already
        ! here? if its really p(b|a) i should switch up i and j..
        ! hm.. the inputted j is the holeInd i which is fixed.. i gets looped
        ! over on the outside and is assumed picked first or 2nd here??
        ! taking i and j here is wrong! i is the open orbital, but j
        ! is the already picked electron! it has to be orb_a here or?
        call pgen_select_orb_guga_mol(csf_i, occ_orbs, i, orb_a, cpt_ba, ba_sum, i, .true.)

        ! change to the fullstart into fullstop: loop until orbital a
        if (tGen_guga_weighted) then
            do orb = i + 1, orb_a - 1
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)
                end if
            end do
        end if

        ! deal with a also speciallly
        cpt_b = get_guga_integral_contrib(occ_orbs, orb_a, -1)

        cum_sum = cum_sum + cpt_b

        ! and get p(a|b)
        ! here the only difference between r2l and l2r fullstarts come into
        ! play!
        if (orb_a > j) then
            ! then orb_j is off-limits
            call pgen_select_orb_guga_mol(csf_i, occ_orbs, orb_a, i, cpt_ab, ab_sum, j)
        else
            ! in this case there is no restriction guga-wise..
            call pgen_select_orb_guga_mol(csf_i, occ_orbs, orb_a, i, cpt_ab, ab_sum)
        end if

        ! and deal with rest:
        if (tGen_guga_weighted) then
            do orb = orb_a + 1, nSpatOrbs
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)
                end if
            end do
        end if

        if (.not. tGen_guga_weighted) then
            cum_sum = csf_i%cum_list(nSpatOrbs)
        end if

        if (near_zero(cum_sum) .or. near_zero(ab_sum) .or. near_zero(ba_sum)) then
            orb_pgen = 0.0_dp
        else
            ! and get hopefully correct final values:
            cpt_a = cpt_a / cum_sum * cpt_ba / ba_sum
            cpt_b = cpt_b / cum_sum * cpt_ab / ab_sum

            ! and add them up to the final orbital pgen
            orb_pgen = cpt_a + cpt_b
        end if

    end subroutine calc_orbital_pgen_contrib_start_def

    subroutine forced_mixed_start(ilut, csf_i, excitInfo, t, probWeight)
        ! NOTE: mixed full-start matrix elements are stores in the same row
        ! as delta B = -1 ones -> so access getDoubleMatrixElement with
        ! db = -1 below!
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: probWeight
        character(*), parameter :: this_routine = "forced_mixed_start"

        real(dp) :: tempWeight, bVal, tempWeight_1
        integer :: st

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(.not. isThree(ilut, excitInfo%fullStart))

        ! cant be 0 or zero matrix element, but cant be 3 either or only
        ! deltaB=0 branch would have non-zero matrix element and that would
        ! to a single-excitation-like DE

        st = excitInfo%fullStart
        bVal = csf_i%B_real(st)

        t = ilut

        select case (csf_i%stepvector(st))
        case (1)
            if (csf_i%B_int(st) < 2) then
                ! actually not possible in this case, as only 0 branch valid..
                probWeight = 0.0_dp
                t = 0_n_int
                return
            end if

            ! otherwise switch 1 -> 2
            clr_one(t, st)
            set_two(t, st)

            call getDoubleMatrixElement(2, 1, -1, gen_type%L, gen_type%R, &
                                        bVal, 1.0_dp, x1_element=tempWeight_1)

            ! x0 matrix element:
            tempWeight = 0.0_dp

            call setDeltaB(2, t)

        case (2)
            ! choose -2 branch 2 -> 1
            clr_two(t, st)
            set_one(t, st)

            call getDoubleMatrixElement(1, 2, -1, gen_type%L, gen_type%R, &
                                        bVal, 1.0_dp, x1_element=tempWeight_1)

            tempWeight = 0.0_dp

            call setDeltaB(-2, t)
#ifdef DEBUG_
        case default
            call stop_all(this_routine, "wrong stepvalue!")
#endif
        end select

        call encode_matrix_element(t, tempWeight_1, 2)
        call encode_matrix_element(t, tempWeight, 1)

        if (near_zero(abs(tempWeight) + abs(tempWeight_1))) then
            probWeight = 0.0_dp
            t = 0_n_int
            return
        end if

        ! and since there is no choice: branch_pgen is 1
        probWeight = 1.0_dp

    end subroutine forced_mixed_start

    subroutine mixedFullStartStochastic(ilut, csf_i, excitInfo, weights, posSwitches, &
                                        negSwitches, t, probWeight)
        ! NOTE: mixed full-start matrix elements are stores in the same row
        ! as delta B = -1 ones -> so access getDoubleMatrixElement with
        ! db = -1 below!
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: probWeight
        character(*), parameter :: this_routine = "mixedFullStartStochastic"

        real(dp) :: minusWeight, plusWeight, zeroWeight, tempWeight, bVal, &
                    tempWeight_1
        integer :: st

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(.not. isThree(ilut, excitInfo%fullStart))

        ! cant be 0 or zero matrix element, but cant be 3 either or only
        ! deltaB=0 branch would have non-zero matrix element and that would
        ! to a single-excitation-like DE

        st = excitInfo%fullStart
        bVal = csf_i%B_real(st)

        t = ilut

        ! another note on the matrix elements... for a mixed full-start
        ! i know that the branch has to switch to +-2 at some point to
        ! yield a non-single excitation -> so i could just ignore the x0
        ! matrix element, as it definetly gets 0 at some point...
        ! but i also then have to write a specific double update function,
        ! which also deals with this kind of "restricted" or better enforced
        ! double excitation regime!

        ! no do not do that, as i have to check probweights inbetween if
        ! excitation yields non-zero matrix element
        select case (csf_i%stepvector(st))
        case (1)
            if (csf_i%B_int(st) < 2) then
                ! only 0-branch start possible, which always has weight > 0
                ! no change in stepvector, so just calc matrix element
                call getDoubleMatrixElement(1, 1, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, tempWeight, tempWeight_1)

                call setDeltaB(0, t)

                probWeight = 1.0_dp
            else
                ! both branches possible
                plusWeight = weights%proc%plus(posSwitches(st), &
                                               bVal, weights%dat)
                zeroWeight = weights%proc%zero(negSwitches(st), &
                                               posSwitches(st), bVal, weights%dat)

                probWeight = calcStartProb(zeroWeight, plusWeight)

                if (genrand_real2_dSFMT() < probWeight) then
                    ! choose 0 branch

                    call getDoubleMatrixElement(1, 1, -1, gen_type%L, gen_type%R, &
                                                bVal, 1.0_dp, tempWeight, tempWeight_1)

                    call setDeltaB(0, t)

                else
                    ! +2 branch:
                    ! switch 1 -> 2
                    clr_orb(t, 2 * st - 1)
                    set_orb(t, 2 * st)

                    call getDoubleMatrixElement(2, 1, -1, gen_type%L, gen_type%R, &
                                                bVal, 1.0_dp, x1_element=tempWeight_1)

                    tempWeight = 0.0_dp

                    call setDeltaB(2, t)

                    probWeight = 1.0_dp - probWeight
                end if

            end if

        case (2)
            ! here always both branches are possible
            minusWeight = weights%proc%minus(negSwitches(st), &
                                             bVal, weights%dat)
            zeroWeight = weights%proc%zero(negSwitches(st), &
                                           posSwitches(st), bVal, weights%dat)

            probWeight = calcStartProb(zeroWeight, minusWeight)

            if (genrand_real2_dSFMT() < probWeight) then
                ! choose 0 branch
                call getDoubleMatrixElement(2, 2, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, tempWeight, tempWeight_1)

                call setDeltaB(0, t)

            else
                ! choose -2 branch:
                ! change 2 -> 1
                set_orb(t, 2 * st - 1)
                clr_orb(t, 2 * st)

                call getDoubleMatrixElement(1, 2, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, x1_element=tempWeight_1)

                tempWeight = 0.0_dp

                call setDeltaB(-2, t)

                probWeight = 1.0_dp - probWeight

            end if
#ifdef DEBUG_
        case default
            call stop_all(this_routine, "wrong stepvalue!")
#endif
        end select

        call encode_matrix_element(t, tempWeight_1, 2)
        call encode_matrix_element(t, tempWeight, 1)

        ! since i only need this routine in excitations, where there has to be
        ! a switch in the double overlap part i could check if it is 0 and then
        ! set probWeight to 0 and return
        if (near_zero(abs(tempWeight) + abs(tempWeight_1))) then
            probWeight = 0.0_dp
            t = 0
            return
        end if

    end subroutine mixedFullStartStochastic

    subroutine calcSingleOverlapMixedStochastic(ilut, csf_i, excitInfo, t, branch_pgen, &
                                                posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcSingleOverlapMixedStochastic"

        type(WeightObj_t) :: weights
        real(dp) :: tempWeight, bVal, temp_pgen
        HElement_t(dp) :: umat
        integer :: iOrb, deltaB

        ! first check the umat element, although if picked correctly it should
        ! be non-zero anyway, there are only 2 symmetric contributions to this
        ! so the 1/2 in the 2-electron part of the hamiltonian get cancelled
        ! and actually have to check type of excitation here...
        ! have to settle on what the multiple picked orbitals is... so
        ! i would not need an if statement here to determine the type of gens
        !todo: can make that without if statement here by using correct ijkl

        if (excitInfo%typ == excit_type%single_overlap_L_to_R) then

            umat = (get_umat_el(excitInfo%firstEnd, excitInfo%secondStart, &
                                excitInfo%fullStart, excitInfo%fullEnd) + &
                    get_umat_el(excitInfo%secondStart, excitInfo%firstEnd, &
                                excitInfo%fullEnd, excitInfo%fullStart)) / 2.0_dp

        else if (excitInfo%typ == excit_type%single_overlap_R_to_L) then

            umat = (get_umat_el(excitInfo%fullStart, excitInfo%fullEnd, &
                                excitInfo%firstEnd, excitInfo%secondStart) + &
                    get_umat_el(excitInfo%fullEnd, excitInfo%fullStart, &
                                excitInfo%secondStart, excitInfo%firstEnd)) / 2.0_dp
        else
            call stop_all(this_routine, "shouldnt be here!")
        end if

        ! todo : correct sum of contributing 2-body integrals and correct
        ! indexing!
        if (near_zero(umat)) then
            branch_pgen = 0.0_dp
            t = 0_n_int
            return
        end if

        ! in the mixed single overlap case its just like a regular single
        ! excitation except the special change in stepvector at the
        ! single overlap site!
        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_singleWeight(csf_i, excitInfo%fullEnd)
        end if

        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)

        ! check if weights were 0
        check_abort_excit(branch_pgen, t)

        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            ! check and update weights
            branch_pgen = branch_pgen * temp_pgen

            check_abort_excit(branch_pgen, t)
        end do

        iOrb = excitInfo%secondStart
        bVal = csf_i%B_real(iOrb)

        deltaB = getDeltaB(t)

        if (excitInfo%firstGen == gen_type%L) then
            ! lowering gen ends here
            ASSERT(isZero(ilut, iOrb))

            ! have to get double excitation matrix elements in here..
            call getDoubleMatrixElement(3, 0, deltaB, excitInfo%firstGen, &
                                        excitInfo%lastGen, bVal, 1.0_dp, tempWeight)

            ! change 0 -> 3
            set_orb(t, 2 * iOrb)
            set_orb(t, 2 * iOrb - 1)

        else
            ! raising gen ends here
            ASSERT(isThree(ilut, iOrb))

            call getDoubleMatrixElement(0, 3, deltaB, excitInfo%firstGen, &
                                        excitInfo%lastGen, bVal, 1.0_dp, tempWeight)

            ! change 3 -> 0
            clr_orb(t, 2 * iOrb)
            clr_orb(t, 2 * iOrb - 1)

        end if
        excitInfo%currentGen = excitInfo%lastGen

        do iOrb = excitInfo%secondStart + 1, excitInfo%fullEnd - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            ! check and update weights
            branch_pgen = branch_pgen * temp_pgen

            check_abort_excit(branch_pgen, t)
        end do

        call singleStochasticEnd(csf_i, excitInfo, t)

        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                                            contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                                               excit_lvl=2, excit_typ=excitInfo%typ), x1=0.0_dp, &
                                            x0=extract_matrix_element(t, 1) * tempWeight)
        end if

        ! for efficiency only encode umat here
        call encode_matrix_element(t, 0.0_dp, 2)
        call update_matrix_element(t, tempWeight * umat, 1)

    end subroutine calcSingleOverlapMixedStochastic

    subroutine calcFullstopRaisingStochastic(ilut, csf_i, excitInfo, t, branch_pgen, &
                                             posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcFullstopRaisingStochastic"

        real(dp) :: nOpen, tempWeight, bVal, temp_pgen
        HElement_t(dp) :: umat
        type(WeightObj_t) :: weights
        integer :: iOrb, deltaB

        ASSERT(isThree(ilut, excitInfo%fullEnd))
        ASSERT(isProperCSF_ilut(ilut))

        ! could check U matrix element here.. although in the final version
        ! of the orbital picker, it probably should be accessed already
        ! for the cauchy-schwarz criteria
        umat = (get_umat_el(excitInfo%fullStart, excitInfo%secondStart, &
                            excitInfo%fullEnd, excitInfo%fullEnd) + get_umat_el( &
                excitInfo%secondStart, excitInfo%fullStart, excitInfo%fullEnd, &
                excitInfo%fullEnd)) / 2.0_dp

        ! todo: correct sum and indexing..
        if (near_zero(umat)) then
            branch_pgen = 0.0_dp
            t = 0_n_int
            return
        end if

        ! create weight object here
        ! i think i only need single excitations weights here, since
        ! the semi stop in this case is like an end step...
        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_singleWeight(csf_i, excitInfo%secondStart)
        end if

        ! i only need normal single stochastic start then..
        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)

        ! matrix element cannot be zero but both branch weights might have been
        check_abort_excit(branch_pgen, t)

        ! then stochastc single update
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            branch_pgen = branch_pgen * temp_pgen
            ! check if branch weights turned 0
            check_abort_excit(branch_pgen, t)
        end do

        ! write it for specific lowering semi-start
        iOrb = excitInfo%secondStart
        bVal = csf_i%B_real(iOrb)

        ! hope that the weighs correctly assert that only fitting deltaB
        ! values arrive here, to ensure a deltaB=0 branch in the DE overlap
        deltaB = getDeltaB(t)
        select case (csf_i%stepvector(iOrb))
        case (1)
            ASSERT(deltaB == -1)

            ! change 1 -> 3
            set_orb(t, 2 * iOrb)

            call getDoubleMatrixElement(3, 1, deltaB, gen_type%R, gen_type%R, &
                                        bVal, excitInfo%order, tempWeight)

            call setDeltaB(0, t)

        case (2)
            ASSERT(deltaB == 1)

            ! change 2 -> 3
            set_orb(t, 2 * iOrb - 1)

            call getDoubleMatrixElement(3, 2, deltaB, gen_type%R, gen_type%R, &
                                        bVal, 1.0_dp, tempWeight)

            call setDeltaB(0, t)

        case (0)
            ASSERT(isZero(ilut, iOrb))

            if (deltaB == -1) then

                ! change 0->2
                set_orb(t, 2 * iOrb)

                call getDoubleMatrixElement(2, 0, deltaB, gen_type%R, gen_type%R, &
                                            bVal, 1.0_dp, tempWeight)

            else
                ! change 0->1
                set_orb(t, 2 * iOrb - 1)

                call getDoubleMatrixElement(1, 0, deltaB, gen_type%R, gen_type%R, &
                                            bVal, 1.0_dp, tempWeight)

            end if

            call setDeltaB(0, t)

#ifdef DEBUG_
        case (3)
            call stop_all(this_routine, "wrong stepvalue!")
#endif
        end select

        ! only deltaB = 0 branch, so only number of open orbitals important
        ! for matrix element in overlap region

        nOpen = (-1.0_dp)**real(&
                    count_open_orbs_ij(csf_i, excitInfo%secondStart + 1, excitInfo%fullEnd - 1), &
                dp)

        iOrb = excitInfo%fullEnd

        ! change 3 -> 0
        clr_orb(t, 2 * iOrb)
        clr_orb(t, 2 * iOrb - 1)

        call update_matrix_element(t, tempWeight * nOpen * Root2, 1)

        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                                            contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                                               excit_lvl=2, excit_typ=excitInfo%typ), x1=0.0_dp, &
                                            x0=extract_matrix_element(t, 1))
        end if

        ! also just encode all matrix element contributions here
        call encode_matrix_element(t, 0.0_dp, 2)
        call update_matrix_element(t, umat, 1)

    end subroutine calcFullstopRaisingStochastic

    subroutine calcFullstopLoweringStochastic(ilut, csf_i, excitInfo, t, branch_pgen, &
                                              posSwitches, negSwitches, opt_weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcFullstopLoweringStochastic"

        real(dp) :: nOpen, tempWeight, bVal, temp_pgen
        HElement_t(dp) :: umat
        type(WeightObj_t) :: weights
        integer :: iOrb, deltaB

        ASSERT(isZero(ilut, excitInfo%fullEnd))
        ASSERT(isProperCSF_ilut(ilut))

        ! could check U matrix element here.. although in the final version
        ! of the orbital picker, it probably should be accessed already
        ! for the cauchy-schwarz criteria
        umat = (get_umat_el(excitInfo%fullEnd, excitInfo%fullEnd, &
                            excitInfo%fullStart, excitInfo%secondStart) + &
                get_umat_el(excitInfo%fullEnd, excitInfo%fullEnd, &
                            excitInfo%secondStart, excitInfo%fullStart)) / 2.0_dp

        ! todo correct combination of sum terms and correct indexcing
        if (near_zero(umat)) then
            branch_pgen = 0.0_dp
            t = 0_n_int
            return
        end if

        ! create weight object here
        ! i think i only need single excitations weights here, since
        ! the semi stop in this case is like an end step...
        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_singleWeight(csf_i, excitInfo%secondStart)
        end if

        ! i only need normal single stochastic start then..
        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, t, branch_pgen)

        ! matrix elements cant be 0, but maybe both branch weights were...
        check_abort_excit(branch_pgen, t)

        ! then stochastc single update
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, posSwitches, &
                                        negSwitches, t, temp_pgen)
            branch_pgen = branch_pgen * temp_pgen
            ! matrix elements cant be 0 but branch weight might be..
            check_abort_excit(branch_pgen, t)
        end do

        ! write it for specific lowering semi-start
        iOrb = excitInfo%secondStart
        bVal = csf_i%B_real(iOrb)

        ! hope that the weighs correctly assert that only fitting deltaB
        ! values arrive here, to ensure a deltaB=0 branch in the DE overlap
        deltaB = getDeltaB(t)

        select case (csf_i%stepvector(iOrb))
        case (1)
            ASSERT(getDeltaB(t) == -1)

            ! change 1 -> 0
            clr_orb(t, 2 * iOrb - 1)

            call getDoubleMatrixElement(0, 1, deltaB, gen_type%L, gen_type%L, &
                                        bVal, excitInfo%order, tempWeight)

            call setDeltaB(0, t)

        case (2)
            ASSERT(getDeltaB(t) == 1)

            ! change 2 -> 0
            clr_orb(t, 2 * iOrb)

            call getDoubleMatrixElement(0, 2, deltaB, gen_type%L, gen_type%L, &
                                        bVal, 1.0_dp, tempWeight)

            call setDeltaB(0, t)

        case (3)
            ASSERT(isThree(ilut, iOrb))

            if (deltaB == -1) then

                ! change 3->2
                clr_orb(t, 2 * iOrb - 1)

                call getDoubleMatrixElement(2, 3, deltaB, gen_type%L, gen_type%L, &
                                            bVal, 1.0_dp, tempWeight)

            else
                ! change 3->1
                clr_orb(t, 2 * iOrb)

                call getDoubleMatrixElement(1, 3, deltaB, gen_type%L, gen_type%L, &
                                            bVal, 1.0_dp, tempWeight)

            end if

            call setDeltaB(0, t)

#ifdef DEBUG_
        case (0)
            call stop_all(this_routine, "wrong stepvalue!")
#endif
        end select

        ! only deltaB = 1 branch, so only number of open orbitals important
        ! for matrix element in overlap region

        nOpen = (-1.0_dp)**real(count_open_orbs_ij(csf_i, excitInfo%secondStart + 1, &
                                                   excitInfo%fullEnd - 1), dp)

        iOrb = excitInfo%fullEnd

        ! change 0 -> 3
        set_orb(t, 2 * iOrb)
        set_orb(t, 2 * iOrb - 1)

        call update_matrix_element(t, tempWeight * nOpen * Root2, 1)

        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                                            contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                                               excit_lvl=2, excit_typ=excitInfo%typ), x1=0.0_dp, &
                                            x0=extract_matrix_element(t, 1))
        end if

        ! update all matrix element contributions at once
        call encode_matrix_element(t, 0.0_dp, 2)
        call update_matrix_element(t, umat, 1)

    end subroutine calcFullstopLoweringStochastic

    subroutine calcFullStartLoweringStochastic(ilut, csf_i, excitInfo, t, branch_pgen, &
                                               posSwitches, negSwitches, opt_weight)
        ! in this case there is no ambiguity in the matrix elements, as they
        ! are uniquely determined and thus can be efficiently calculated on
        ! the fly
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcFullStartLoweringStochastic"

        real(dp) :: tempWeight, minusWeight, plusWeight, nOpen, bVal, temp_pgen
        HElement_t(dp) :: umat
        integer :: start, ende, semi, gen, iOrb, deltaB
        type(WeightObj_t) :: weights

        ASSERT(isThree(ilut, excitInfo%fullStart))
        ASSERT(isProperCSF_ilut(ilut))

        ! create the fullStart
        start = excitInfo%fullStart
        ende = excitInfo%fullEnd
        semi = excitInfo%firstEnd
        gen = excitInfo%firstGen
        bVal = csf_i%B_real(semi)

        ! could check U matrix element here.. although in the final version
        ! of the orbital picker, it probably should be accessed already
        ! for the cauchy-schwarz criteria

        ! in fullstart lowering the x1-elements are always zero -> so no
        ! minus sign in the integral contribution -> so just add the two
        ! relevant elements!
        umat = (get_umat_el(ende, semi, start, start) + &
                get_umat_el(semi, ende, start, start)) / 2.0_dp

        ! todo! correct combination of umat sum terms.. and correct indexing
        if (near_zero(umat)) then
            branch_pgen = 0.0_dp
            t = 0_n_int
            return
        end if

        ! set t
        t = ilut

        ! set 3->0
        clr_orb(t, 2 * start)
        clr_orb(t, 2 * start - 1)

        ! double overlap region influence only determined by number of open
        ! orbitals
        nOpen = real(count_open_orbs_ij(csf_i, start, semi - 1), dp)

        deltaB = getDeltaB(t)

        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_singleWeight(csf_i, ende)
        end if
        ! could encode matrix element here, but do it later for more efficiency

        tempWeight = 0.0_dp

        select case (csf_i%stepvector(semi))
        case (0)
            if (csf_i%B_int(semi) > 0) then

                minusWeight = weights%proc%minus(negSwitches(semi), bVal, weights%dat)
                plusWeight = weights%proc%plus(posSwitches(semi), bVal, weights%dat)

                if (near_zero(minusWeight + plusWeight)) then
                    branch_pgen = 0.0_dp
                    t = 0_n_int
                    return
                end if

                branch_pgen = calcStartProb(minusWeight, plusWeight)

                if (genrand_real2_dSFMT() < branch_pgen) then
                    ! do -1 branch
                    ! do 0->1 first -1 branch first
                    set_orb(t, 2 * semi - 1)

                    ! order does not matter since only x0 matrix element
                    call getDoubleMatrixElement(1, 0, deltaB, gen_type%L, gen_type%L, bVal, &
                                                1.0_dp, tempWeight)

                    call setDeltaB(-1, t)

                else
                    ! do +1 branch
                    ! then do 0->2: +1 branch(change from already set 1
                    set_orb(t, 2 * semi)

                    call getDoubleMatrixElement(2, 0, deltaB, gen_type%L, gen_type%L, bVal, &
                                                1.0_dp, tempWeight)

                    call setDeltaB(1, t)
                    branch_pgen = 1.0_dp - branch_pgen

                end if
            else
                ! only 0->1, -1 branch possible
                set_orb(t, 2 * semi - 1)

                ! order does not matter since only x0 matrix element
                call getDoubleMatrixElement(1, 0, deltaB, gen_type%L, gen_type%L, bVal, &
                                            1.0_dp, tempWeight)

                call setDeltaB(-1, t)

                branch_pgen = 1.0_dp
            end if
        case (1)
            ! only one excitation possible, which also has to have
            ! non-zero weight or otherwise i wouldnt even be here
            ! and reuse already provided t
            ! 1 -> 3
            set_orb(t, 2 * semi)

            call setDeltaB(1, t)
            ! still get the matrix elements here to deal with this
            ! generally
            call getDoubleMatrixElement(3, 1, 0, gen_type%L, gen_type%L, bVal, &
                                        1.0_dp, tempWeight)

            branch_pgen = 1.0_dp

        case (2)
            ! here we have
            ! 2 -> 3
            set_orb(t, 2 * semi - 1)

            call setDeltaB(-1, t)

            call getDoubleMatrixElement(3, 2, 0, gen_type%L, gen_type%L, bVal, &
                                        1.0_dp, tempWeight)

            branch_pgen = 1.0_dp

            ! encode fullstart contribution and pseudo overlap region here

        end select

        call encode_matrix_element(t, 0.0_dp, 2)
        call encode_matrix_element(t, Root2 * tempWeight * (-1.0_dp)**nOpen, 1)

        ! and then we have to do just a regular single excitation
        do iOrb = semi + 1, ende - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, &
                                        posSwitches, negSwitches, t, temp_pgen)
            ! matrix element cant be 0 but maybe both branch weights were..
            branch_pgen = branch_pgen * temp_pgen
            check_abort_excit(branch_pgen, t)
        end do

        call singleStochasticEnd(csf_i, excitInfo, t)

        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                                            contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                                               excit_lvl=2, excit_typ=excitInfo%typ), x1=0.0_dp, &
                                            x0=extract_matrix_element(t, 1))
        end if

        call update_matrix_element(t, umat, 1)

    end subroutine calcFullStartLoweringStochastic

    subroutine calcFullStartRaisingStochastic(ilut, csf_i, excitInfo, t, &
                                              branch_pgen, posSwitches, negSwitches, opt_weight)
        ! in this case there is no ambiguity in the matrix elements, as they
        ! are uniquely determined and thus can be efficiently calculated on
        ! the fly
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out) :: t(0:nifguga)
        real(dp), intent(out) :: branch_pgen
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in), optional :: opt_weight
        character(*), parameter :: this_routine = "calcFullStartRaisingStochastic"

        real(dp) :: tempWeight, minusWeight, plusWeight, nOpen, bVal, temp_pgen
        HElement_t(dp) :: umat
        integer :: start, ende, semi, gen, iOrb, deltaB
        type(WeightObj_t) :: weights

        ASSERT(isZero(ilut, excitInfo%fullStart))
        ASSERT(isProperCSF_ilut(ilut))

        ! create the fullStart
        start = excitInfo%fullStart
        ende = excitInfo%fullEnd
        semi = excitInfo%firstEnd
        gen = excitInfo%firstGen
        bVal = csf_i%B_real(semi)

        ! could check U matrix element here.. although in the final version
        ! of the orbital picker, it probably should be accessed already
        ! for the cauchy-schwarz criteria
        ! but anyway need it for the final matrix element
        ! here in this case only this and the version with both index pairs
        ! exchange correspond to the reached excitation -> so the 1/2 in fron
        ! of the 2 electron part cancels in the hamiltonian.
        umat = (get_umat_el(start, start, semi, ende) + &
                get_umat_el(start, start, ende, semi)) / 2.0_dp

        ! but its not only this matrix element isnt it? have to also take the
        ! one with exchanged indices into account, but where only the type
        ! excitation determines the sign between the 2...
        ! also have to correctly index the routine with phycisist notation..
        ! check if i can deduce that by only the order entries of excitInfo..

        ! better do a excitation abortion here !TODO!! change that to the
        ! correct sum of umats
        if (near_zero(umat)) then
            branch_pgen = 0.0_dp
            t = 0_n_int
            return
        end if

        ! set t
        t = ilut

        ! set 0->3
        set_orb(t, 2 * start)
        set_orb(t, 2 * start - 1)

        ! double overlap region influence only determined by number of open
        ! orbitals
        nOpen = real(count_open_orbs_ij(csf_i, start, semi - 1), dp)

        deltaB = getDeltaB(t)
        if (present(opt_weight)) then
            weights = opt_weight
        else
            weights = init_singleWeight(csf_i, ende)
        end if

        ! could encode matrix element here, but do it later for more efficiency

        select case (csf_i%stepvector(semi))
        case (3)
            if (csf_i%B_int(semi) > 0) then
                minusWeight = weights%proc%minus(negSwitches(semi), bVal, weights%dat)
                plusWeight = weights%proc%plus(posSwitches(semi), bVal, weights%dat)

                ! if both branches are zero, i have to abort the excitation
                ! altough that shouldnt happen...
                if (near_zero(minusWeight + plusWeight)) then
                    branch_pgen = 0.0_dp
                    t = 0_n_int
                    return
                end if

                branch_pgen = calcStartProb(minusWeight, plusWeight)

                if (genrand_real2_dSFMT() < branch_pgen) then
                    ! do -1 branch
                    ! do 3->1 first -1 branch first
                    clr_orb(t, 2 * semi)

                    ! order does not matter since only x0 matrix element
                    call getDoubleMatrixElement(1, 3, deltaB, gen_type%R, gen_type%R, bVal, &
                                                1.0_dp, tempWeight)

                    call setDeltaB(-1, t)

                else
                    ! do +1 branch
                    ! then do 3->2: +1 branch(change from already set 1
                    clr_orb(t, 2 * semi - 1)

                    call getDoubleMatrixElement(2, 3, deltaB, gen_type%R, gen_type%R, bVal, &
                                                1.0_dp, tempWeight)

                    call setDeltaB(1, t)
                    branch_pgen = 1.0_dp - branch_pgen

                end if
            else
                ! only 3->1, -1 branch possible
                clr_orb(t, 2 * semi)

                ! order does not matter since only x0 matrix element
                call getDoubleMatrixElement(1, 3, deltaB, gen_type%R, gen_type%R, bVal, &
                                            1.0_dp, tempWeight)

                call setDeltaB(-1, t)

                branch_pgen = 1.0_dp
            end if
        case (1)
            ! only one excitation possible, which also has to have
            ! non-zero weight or otherwise i wouldnt even be here
            ! and reuse already provided t
            ! 1 -> 0
            clr_orb(t, 2 * semi - 1)

            call getDoubleMatrixElement(0, 1, 0, gen_type%R, gen_type%R, bVal, &
                                        1.0_dp, tempWeight)

            call setDeltaB(1, t)

            branch_pgen = 1.0_dp

        case (2)

            ! here we have
            ! 2 -> 0
            clr_orb(t, 2 * semi)

            call getDoubleMatrixElement(0, 2, 0, gen_type%R, gen_type%R, bVal, &
                                        1.0_dp, tempWeight)

            call setDeltaB(-1, t)

            branch_pgen = 1.0_dp

            ! encode fullstart contribution and pseudo overlap region here
            ! too in one go. one additional -1 due to semistop

        end select

        call encode_matrix_element(t, 0.0_dp, 2)
        call encode_matrix_element(t, Root2 * tempWeight * (-1.0_dp)**nOpen, 1)

        ! x0 matrix elements cant be 0 at a RR semistop
        ! and then we have to do just a regular single excitation
        do iOrb = semi + 1, ende - 1
            call singleStochasticUpdate(ilut, csf_i, iOrb, excitInfo, weights, &
                                        posSwitches, negSwitches, t, temp_pgen)
            branch_pgen = branch_pgen * temp_pgen
            ! in the single overlap regions matrix elements cant actually
            ! become 0! -> so no check for matrix element is necessary here
            ! but i could be set to zero, due to both branches having
            ! zero probabilistic weight, which also should not happen, but just
            ! to be sure! this gets indicated by probWeight = 0 and t = 0
            check_abort_excit(branch_pgen, t)
        end do

        call singleStochasticEnd(csf_i, excitInfo, t)

        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                                            contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                                               excit_lvl=2, excit_typ=excitInfo%typ), x1=0.0_dp, &
                                            x0=extract_matrix_element(t, 1))
        end if

        call update_matrix_element(t, umat, 1)

    end subroutine calcFullStartRaisingStochastic

    subroutine singleStochasticEnd(csf_i, excitInfo, t)
        ! routine to end a stochastic excitation for a single generator
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(inout) :: t(0:nifguga)
        character(*), parameter :: this_routine = "singleStochasticEnd"

        integer :: ende, gen, deltaB
        real(dp) :: bVal, tempWeight

        ende = excitInfo%fullEnd
        deltaB = getDeltaB(t)
        tempWeight = extract_matrix_element(t, 1)
        bVal = csf_i%B_real(ende)
        gen = excitInfo%currentGen
        ! should be really similar to full excitation routine, except all the
        ! clean up stuff
        select case (csf_i%stepvector(ende))
        case (0)
            ! if it is zero it implies i have a lowering generator, otherwise
            ! i wouldnt even be here.
            if (deltaB == -1) then
                ! set alpha bit 0 -> 2
                set_orb(t, 2 * ende)

                tempWeight = getSingleMatrixElement(2, 0, deltaB, gen, bVal)

            else
                ! set beta bit 0 -> 1
                set_orb(t, 2 * ende - 1)

                tempWeight = getSingleMatrixElement(1, 0, deltaB, gen, bVal)

            end if

        case (3)
            ! if d = 3, it implies a raising generator or otherwise we
            ! would not be here ..
            if (deltaB == -1) then
                ! clear beta bit
                clr_orb(t, 2 * ende - 1)

                tempWeight = getSingleMatrixElement(2, 3, deltaB, gen, bVal)

            else
                ! clear alpha bit
                clr_orb(t, 2 * ende)

                tempWeight = getSingleMatrixElement(1, 3, deltaB, gen, bVal)

            end if

        case (1)

            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                if (deltaB /= -1) then
                    t = 0
                    return
                end if
            end if

            ASSERT(deltaB == -1)

            if (gen == gen_type%R) then
                clr_orb(t, 2 * ende - 1)

                tempWeight = getSingleMatrixElement(0, 1, deltaB, gen, bVal)

            else ! lowering generator

                ! for lowering gen: 1 -> 3: set alpha bit
                set_orb(t, 2 * ende)

                tempWeight = getSingleMatrixElement(3, 1, deltaB, gen, bVal)

            end if

        case (2)
            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                if (deltaB /= 1) then
                    t = 0
                    return
                end if
            end if
            ASSERT(deltaB == 1)

            if (gen == gen_type%R) then
                ! for raising: 2 -> 0: clear alpha bit
                clr_orb(t, 2 * ende)

                tempWeight = getSingleMatrixElement(0, 2, deltaB, gen, bVal)

            else ! lowering gen
                ! for lowering gen: 2 -> 3 : set beta bit
                set_orb(t, 2 * ende - 1)

                tempWeight = getSingleMatrixElement(3, 2, deltaB, gen, bVal)
            end if
        end select

        call update_matrix_element(t, tempWeight, 1)
        call update_matrix_element(t, tempWeight, 2)

        ! do excitaiton abortion
        if (near_zero(extract_matrix_element(t, 1)) .and. &
            near_zero(extract_matrix_element(t, 2))) then
            t = 0_n_int
        end if

    end subroutine singleStochasticEnd

    subroutine singleStochasticUpdate(ilut, csf_i, s, excitInfo, weights, posSwitches, &
                                      negSwitches, t, probWeight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: s
        type(ExcitationInformation_t), intent(in) :: excitInfo
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer(n_int), intent(out) :: t(0:nifguga) ! to use macros, have to use short names..
        real(dp), intent(out) :: probWeight
        character(*), parameter :: this_routine = "singleStochasticUpdate"

        real(dp) :: tempWeight, bVal, minusWeight, plusWeight
        integer :: gen, deltaB, step
        ! is also very similar to full single update
        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(s > 0 .and. s <= nSpatOrbs)

        ! change! have to update the x1 matrix elements here too,
        ! to keep them seperated to use the individually for the contributing
        ! integrals

        ! set standard probWeight
        probWeight = 1.0_dp

        select case (csf_i%stepvector(s))
        case (0)
            ! do nothin actually.. not even change matrix elements
            return

        case (3)
            ! only change matrix element to negative one
            call update_matrix_element(t, -1.0_dp, 1)
            call update_matrix_element(t, -1.0_dp, 2)
            return

        end select

        ! to generally use this function i need to define a current generator...
        gen = excitInfo%currentGen
        bVal = csf_i%B_real(s)

        ! only need weights in certain cases..
        deltaB = getDeltaB(t)

        ! is it possible to only check for possible switches and else handle
        ! it in one go? try!
        select case (csf_i%stepvector(s) + deltaB)
        case (0)
            ! d=1 + b=-1 : 0
            ! no bValue restrictions here, so that should be handlebar
            plusWeight = weights%proc%plus(posSwitches(s), bVal, weights%dat)
            minusWeight = weights%proc%minus(negSwitches(s), bVal, weights%dat)

            ! here do a check if not both weights are 0...
            if (near_zero(plusWeight + minusWeight)) then
                probWeight = 0.0_dp
                t = 0
                return
            end if

            ! calc. staying probabiliy
            probWeight = calcStayingProb(minusWeight, plusWeight, bVal)

            if (genrand_real2_dSFMT() < probWeight) then
                ! stay on -1 branch
                tempWeight = getSingleMatrixElement(1, 1, deltaB, gen, bVal)

            else
                ! switch to +1 branch 1 -> 2
                set_orb(t, 2 * s)
                clr_orb(t, 2 * s - 1)

                call setDeltaB(1, t)

                tempWeight = getSingleMatrixElement(2, 1, deltaB, gen, bVal)

                probWeight = 1.0_dp - probWeight
            end if

        case (3)
            ! d=2 + b=1 : 3
            ! do i need bValue check here?
            ! probably... to distinguish forced switches
            if (csf_i%B_int(s) > 0) then
                ! both excitations possible
                plusWeight = weights%proc%plus(posSwitches(s), bVal, weights%dat)
                minusWeight = weights%proc%minus(negSwitches(s), bVal, weights%dat)
                ! here do a check if not both weights are 0...
                if (near_zero(plusWeight + minusWeight)) then
                    probWeight = 0.0_dp
                    t = 0
                    return
                end if

                probWeight = calcStayingProb(plusWeight, minusWeight, bVal)

                if (genrand_real2_dSFMT() < probWeight) then
                    ASSERT(plusWeight > 0.0_dp)
                    ! stay on +1 branch
                    tempWeight = getSingleMatrixElement(2, 2, deltaB, gen, bVal)

                else

                    ASSERT(minusWeight > 0.0_dp)
                    ! do the 2 -> 1 switch
                    clr_orb(t, 2 * s)
                    set_orb(t, 2 * s - 1)

                    ! change deltaB
                    call setDeltaB(-1, t)

                    ! and calc matrix elements
                    tempWeight = getSingleMatrixElement(1, 2, deltaB, gen, bVal)

                    probWeight = 1.0_dp - probWeight

                end if

            else
                ! forced switch
#ifdef DEBUG_
                minusWeight = weights%proc%minus(negSwitches(s), bVal, weights%dat)
                if (.not. minusWeight > 0.0_dp) then
                    print *, "+", plusWeight
                    print *, "-", minusWeight
                    print *, negSwitches
                    print *, "overlap:", excitInfo%overlap

                    call write_det_guga(stdout, ilut)
                    call write_det_guga(stdout, t)
                    call print_excitInfo(excitInfo)
                end if

                ASSERT(minusWeight > 0.0_dp)
#endif
                ! do 2 -> 1forced switch
                clr_orb(t, 2 * s)
                set_orb(t, 2 * s - 1)

                ! change deltaB
                call setDeltaB(-1, t)

                ! and calc matrix elements
                tempWeight = getSingleMatrixElement(1, 2, deltaB, gen, bVal)
            end if
        case (1, 2)
            ! d=1 + b=1:  2
            ! d=2 + b=-1: 1

            ! just staying possibility...

            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                plusWeight = weights%proc%plus(posSwitches(s), bVal, weights%dat)
                minusWeight = weights%proc%minus(negSwitches(s), bVal, weights%dat)

                if ((deltaB == -1 .and. near_zero(minusWeight)) &
                    .or. (deltaB == 1 .and. near_zero(plusWeight))) then
                    t = 0_n_int
                    probWeight = 0.0_dp
                    return
                end if
            end if

#ifdef DEBUG_
            ! for sanity check for weights here too just to be sure to not get
            ! into non compatible excitations
            plusWeight = weights%proc%plus(posSwitches(s), bVal, weights%dat)
            minusWeight = weights%proc%minus(negSwitches(s), bVal, weights%dat)
            if (deltaB == -1) then
                ASSERT(minusWeight > 0.0_dp)
            end if
            if (deltaB == +1) then
                ASSERT(plusWeight > 0.0_dp)
            end if
#endif

            ! can i efficiently code that up?
            ! deltaB stays the same.., stepvector stays the same.. essentiall
            ! only matrix element changes.. -> need deltaB, generators and stepvalue
            step = csf_i%stepvector(s)

            ! to get correct bValue use it for next spatial orbital..
            tempWeight = getSingleMatrixElement(step, step, deltaB, gen, bVal)

        end select

        call update_matrix_element(t, tempWeight, 1)
        call update_matrix_element(t, tempWeight, 2)
        if (t_trunc_guga_matel) then
            if (abs(extract_matrix_element(t, 1)) < trunc_guga_matel) then
                probWeight = 0.0_dp
                t = 0_n_int
            end if
        end if

    end subroutine singleStochasticUpdate

    subroutine createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                            negSwitches, t, probWeight)
        ! create a stochastic start for a single generator
        ! essentially exactly the same as in the full excitation calculation
        ! except that at a 0/3 start we have to do it stochastically ...
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        type(WeightObj_t), intent(in) :: weights
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer(n_int), intent(out) :: t(0:nifguga) ! to use macros, have to use short names..
        real(dp), intent(out) :: probWeight
        character(*), parameter :: this_routine = "createStochasticStart_single"

        real(dp) :: tempWeight, minusWeight, plusWeight, bVal
        integer :: st, gen

        ! and also the possibility of the excitation should be ensured due to
        ! the correct choosing of excitation indices..
#ifdef DEBUG_
        ! also assert we are not calling it for a weight gen. accidently
        ASSERT(excitInfo%currentGen /= 0)
        ASSERT(isProperCSF_ilut(ilut))
        ! also check if calculated b vector really fits to ilut
        ASSERT(all(csf_i%B_real .isclose. calcB_vector_ilut(ilut(0:nifd))))
        if (excitInfo%currentGen == gen_type%R) then
            ASSERT(.not. isThree(ilut, excitInfo%fullStart))
        else if (excitInfo%currentGen == gen_type%L) then
            ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        end if
        ! also have checked if atleast on branch way can lead to an excitaiton
#endif

        ! for more oversight
        st = excitInfo%fullStart
        gen = excitInfo%currentGen
        bVal = csf_i%B_real(st)

        t = ilut

        select case (csf_i%stepvector(st))
        case (1)
            ! set corresponding orbital to 0 or 3 depending on generator type
            if (gen == gen_type%R) then ! raising gen case
                ! set the alpha orbital also to 1 to make d=1 -> d'=3
                set_orb(t, 2 * st)

                ! does it work like that:
                ! would have all necessary values and if statements to directly
                ! calculated matrix element here... maybe do it here after all
                ! to be more efficient...
                tempWeight = getSingleMatrixElement(3, 1, +1, gen, bVal)

            else ! lowering gen case
                ! clear beta orbital to make d=1 -> d'=0
                clr_orb(t, 2 * st - 1)

                tempWeight = getSingleMatrixElement(0, 1, +1, gen, bVal)

            end if
            ! set deltaB:
            ! for both cases deltaB = +1, set that as signed weight
            ! update: use previously unused flags to encode deltaB
            call setDeltaB(1, t)

            ! have to set probabilistic weight to 1, since only 1 possibility
            probWeight = 1.0_dp

        case (2)
            if (gen == gen_type%R) then
                set_orb(t, 2 * st - 1)

                ! matrix elements
                tempWeight = getSingleMatrixElement(3, 2, -1, gen, bVal)

            else
                clr_orb(t, 2 * st)

                ! matrix elements
                tempWeight = getSingleMatrixElement(0, 2, -1, gen, bVal)

            end if
            call setDeltaB(-1, t)
            probWeight = 1.0_dp

        case (0, 3)
            ! if the deltaB = -1 weights is > 0 this one always works
            ! write weight calculation function! is a overkill here,
            ! but makes it easier to use afterwards with stochastic excitaions

            ! use new weight objects here.
            minusWeight = weights%proc%minus(negSwitches(st), bVal, weights%dat)
            plusWeight = weights%proc%plus(posSwitches(st), bVal, weights%dat)

            ! if both weights an b value allow both excitations, have to choose
            ! one stochastically

            ! also have to consider, that i might not actually can assure
            ! possible excitations, since im not yet sure if i check for
            ! possible switches... -> check that todo
            ! for now assume i checked for that, with checkCompatibility() or
            ! something similar..
            ! do not need the previous if structs since the weights
            ! already care for the facts if some excitation is incompatible..
            ! I only have to avoid a completely improper excitation...
            ! if b == 0 and minusWeight == 0 -> no excitation possible todo!

            if (near_zero(plusWeight + minusWeight) .or. &
                near_zero(bVal + minusWeight)) then
                probWeight = 0.0_dp
                t = 0
                return
            end if

            ! calc starting weight to go for the minusBranch
            probWeight = calcStartProb(minusWeight, plusWeight)

            if (genrand_real2_dSFMT() < probWeight) then
                ! do the -1 branch
                if (gen == gen_type%R) then
                    ASSERT(isZero(ilut, st))
                    ! set beta bit
                    set_orb(t, 2 * st - 1)

                    ! matrix element
                    tempWeight = getSingleMatrixElement(1, 0, -1, gen, bVal)

                else
                    ASSERT(isThree(ilut, st))
                    ! clear alpha bit
                    clr_orb(t, 2 * st)

                    ! matrix element
                    tempWeight = getSingleMatrixElement(1, 3, -1, gen, bVal)

                end if

                call setDeltaB(-1, t)

            else
                ! do the +1 branch
                if (gen == gen_type%R) then
                    ASSERT(isZero(ilut, st))
                    ! set alpha bit
                    set_orb(t, 2 * st)

                    ! matrix element
                    tempWeight = getSingleMatrixElement(2, 0, +1, gen, bVal)

                else
                    ASSERT(isThree(ilut, st))
                    ! cleat beta bit
                    clr_orb(t, 2 * st - 1)

                    ! matrix element
                    tempWeight = getSingleMatrixElement(2, 3, +1, gen, bVal)
                end if
                call setDeltaB(+1, t)
                ! update the probabilistic weight with the actual use one
                probWeight = 1.0_dp - probWeight
            end if
        end select

        call encode_matrix_element(t, tempWeight, 1)

    end subroutine createStochasticStart_single


    elemental function init_forced_end_exchange_weight(csf_i, sOrb) result(forced_double)
        ! obj has the same structure as the semi-start weight, reuse them!
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb
        type(WeightObj_t) :: forced_double
        character(*), parameter :: this_routine = "init_forced_end_exchange_weight"

        ASSERT(sOrb > 0 .and. sOrb <= nSpatOrbs)

        forced_double%dat%F = endFx(csf_i, sOrb)
        forced_double%dat%G = endGx(csf_i, sOrb)

        forced_double%proc%minus => getMinus_double
        forced_double%proc%plus => getPlus_double
        forced_double%proc%zero => get_forced_zero_double

        forced_double%initialized = .true.
    end function init_forced_end_exchange_weight

    function init_forced_end_semistart_weight(csf_i, sOrb, pOrb, negSwitches, posSwitches, bVal) &
        result(forced_semistart)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb, pOrb
        real(dp), intent(in) :: negSwitches, posSwitches, bVal
        type(WeightObj_t) :: forced_semistart

        ! TODO(@Oskar): I think this is invalid code
        type(WeightObj_t), target, save :: double

        forced_semistart%dat%F = endFx(csf_i, sorb)
        forced_semistart%dat%G = endGx(csf_i, sorb)

        double = init_forced_end_exchange_weight(csf_i, porb)

        forced_semistart%ptr => double

        forced_semistart%dat%minus = double%proc%minus(negSwitches, bVal, double%dat)
        forced_semistart%dat%plus = double%proc%plus(posSwitches, bVal, double%dat)
        forced_semistart%dat%zero = double%proc%zero(negSwitches, posSwitches, bVal, double%dat)

        forced_semistart%proc%minus => getMinus_semiStart
        forced_semistart%proc%plus => getPlus_semiStart

        forced_semistart%initialized = .true.
    end function init_forced_end_semistart_weight

    function init_fullDoubleWeight(csf_i, sOrb, pOrb, oOrb, negSwitches1, &
                                   negSwitches2, posSwitches1, posSwitches2, bVal1, bVal2) result(fullDouble)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb, pOrb, oOrb
        real(dp), intent(in) :: negSwitches1, negSwitches2, posSwitches1, &
                                posSwitches2, bVal1, bVal2
        type(WeightObj_t) :: fullDouble
        character(*), parameter :: this_routine = "init_fullDoubleWeight"

        type(WeightObj_t), target, save :: fullStart

        ASSERT(sOrb > 0 .and. sOrb <= nSpatOrbs)
        ASSERT(pOrb > 0 .and. pOrb <= nSpatOrbs)
        ASSERT(oOrb > 0 .and. oOrb <= nSpatOrbs)
        ASSERT(negSwitches1 >= 0.0_dp)
        ASSERT(negSwitches2 >= 0.0_dp)
        ASSERT(posSwitches1 >= 0.0_dp)
        ASSERT(posSwitches2 >= 0.0_dp)

        fullDouble%dat%F = endFx(csf_i, sOrb)
        fullDouble%dat%G = endGx(csf_i, sOrb)

        fullStart = init_fullStartWeight(csf_i, pOrb, oOrb, negSwitches2, &
                                         posSwitches2, bVal2)

        fullDouble%ptr => fullStart

        fullDouble%dat%minus = fullStart%proc%minus(negSwitches1, bVal1, fullStart%dat)
        fullDouble%dat%plus = fullStart%proc%plus(posSwitches1, bVal1, fullStart%dat)
        fullDouble%dat%zero = fullStart%proc%zero(negSwitches1, posSwitches1, &
                                                  bVal1, fullStart%dat)

        fullDouble%proc%minus => getMinus_semiStart
        fullDouble%proc%plus => getPlus_semiStart

        fullDouble%initialized = .true.

    end function init_fullDoubleWeight

    function init_singleOverlapRaising(csf_i, sOrb, pOrb, negSwitches, posSwitches, &
                                       bVal) result(singleRaising)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb, pOrb
        real(dp), intent(in) :: negSwitches, posSwitches, bVal
        type(WeightObj_t) :: singleRaising
        character(*), parameter :: this_routine = "init_singleOverlapRaising"

        type(WeightObj_t) :: single

        ASSERT(sOrb > 0 .and. sOrb <= nSpatOrbs)
        ASSERT(pOrb > 0 .and. pOrb <= nSpatOrbs)
        ASSERT(negSwitches >= 0.0_dp)
        ASSERT(posSwitches >= 0.0_dp)

        ! misuse zero data element as this L-function in this case!
        singleRaising%dat%F = endFx(csf_i, sOrb)
        singleRaising%dat%G = endGx(csf_i, sOrb)
        singleRaising%dat%zero = endLx(csf_i, sOrb)

        single = init_singleWeight(csf_i, pOrb)

        singleRaising%dat%minus = single%proc%minus(negSwitches, bVal, single%dat)
        singleRaising%dat%plus = single%proc%plus(posSwitches, bVal, single%dat)

        singleRaising%proc%minus => getMinus_overlapRaising
        singleRaising%proc%plus => getPlus_overlapRaising

        singleRaising%initialized = .true.

    end function init_singleOverlapRaising

    function getMinus_overlapRaising(nSwitches, bVal, dat) result(minusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: dat
        real(dp) :: minusWeight
        character(*), parameter :: this_routine = "getMinus_overlapRaising"

        ASSERT(nSwitches >= 0.0_dp)

        if (near_zero(bVal)) then
            !minusWeight = dat%F*(dat%minus + (1.0_dp - dat%G)*dat%plus) + &
            !    nSwitches*dat%G*(dat%zero*dat%plus + (1.0_dp - dat%F)*dat%minus)
            minusWeight = dat%F * (dat%minus + (1.0_dp - dat%G) * dat%plus) + &
                          nSwitches * dat%G * (dat%plus + (1.0_dp - dat%F) * dat%minus)

        else
            ! minusWeight = dat%F*(dat%minus + (1.0_dp - dat%G)*dat%plus) + &
            !     nSwitches*dat%G/bVal*(dat%zero*dat%plus + (1.0_dp - dat%F)*dat%minus)
            minusWeight = dat%F * (dat%minus + (1.0_dp - dat%G) * dat%plus) + &
                          nSwitches * dat%G / bVal * (dat%plus + (1.0_dp - dat%F) * dat%minus)
        end if

        ASSERT(minusWeight >= 0.0_dp)

    end function getMinus_overlapRaising

    function getPlus_overlapRaising(nSwitches, bVal, dat) result(plusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: dat
        real(dp) :: plusWeight
        character(*), parameter :: this_routine = "getPlus_overlapRaising"

        ASSERT(nSwitches >= 0.0_dp)

        if (near_zero(bVal)) then
            plusWeight = 0.0_dp
        else
            !plusWeight = dat%G*(dat%zero*dat%plus + (1.0_dp - dat%F)*dat%minus) + &
            !    nSwitches*dat%F/bVal * (dat%minus + (1.0_dp - dat%G)*dat%plus)
            plusWeight = dat%G * (dat%plus + (1.0_dp - dat%F) * dat%minus) + &
                         nSwitches * dat%F / bVal * (dat%minus + (1.0_dp - dat%G) * dat%plus)
        end if

        ASSERT(plusWeight >= 0.0_dp)

    end function getPlus_overlapRaising

    function init_singleOverlapLowering(csf_i, sOrb, pOrb, negSwitches, &
                                        posSwitches, bVal) result(singleLowering)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb, pOrb
        real(dp), intent(in) :: negSwitches, posSwitches, bVal
        type(WeightObj_t) :: singleLowering
        character(*), parameter :: this_routine = "init_singleOverlapLowering"

        type(WeightObj_t), target, save :: single

        ASSERT(sOrb > 0 .and. sOrb <= nSpatOrbs)
        ASSERT(pOrb > 0 .and. pOrb <= nSpatOrbs)
        ASSERT(negSwitches >= 0.0_dp)
        ASSERT(posSwitches >= 0.0_dp)

        ! misuse zero data element as this L-function in this case!
        singleLowering%dat%F = endFx(csf_i, sOrb)
        singleLowering%dat%G = endGx(csf_i, sOrb)

        single = init_singleWeight(csf_i, pOrb)

        singleLowering%ptr => single

        singleLowering%dat%minus = single%proc%minus(negSwitches, bVal, single%dat)
        singleLowering%dat%plus = single%proc%plus(posSwitches, bVal, single%dat)

        singleLowering%proc%minus => getMinus_overlapLowering
        singleLowering%proc%plus => getPlus_overlapLowering

    end function init_singleOverlapLowering

    function getMinus_overlapLowering(nSwitches, bVal, dat) result(minusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: dat
        real(dp) :: minusWeight
        character(*), parameter :: this_routine = "getMinus_overlapLowering"

        ASSERT(nSwitches >= 0.0_dp)

        ! if at the current spot b is 0, only the -1 branch is technically
        ! possible and the +1 branch always is forbidden. so i essentially
        ! just have to check if the -1 branch has a non-zero weight, and
        ! set the +1 branch to 0 probability
        ! if non-zero, wird die -1 wahrscheinlichkeit dann eh zu eins
        ! normalisiert.
        if (near_zero(bVal)) then
            minusWeight = dat%G * dat%minus + dat%F * (1.0_dp - dat%G) * dat%plus + &
                          nSwitches * (dat%F * dat%plus + dat%G * (1.0_dp - dat%F) * dat%minus)

        else
            minusWeight = dat%G * dat%minus + dat%F * (1.0_dp - dat%G) * dat%plus + &
                          nSwitches / bVal * (dat%F * dat%plus + dat%G * (1.0_dp - dat%F) * dat%minus)

        end if

        ASSERT(minusWeight >= 0.0_dp)

    end function getMinus_overlapLowering

    function getPlus_overlapLowering(nSwitches, bVal, dat) result(plusWeight)
        real(dp), intent(in) :: nSwitches, bVal
        type(WeightData_t), intent(in) :: dat
        real(dp) :: plusWeight
        character(*), parameter :: this_routine = "getPlus_overlapLowering"

        ASSERT(nSwitches >= 0.0_dp)

        ! if b == 0, set the plus weight to zero, to handle that case
        if (near_zero(bVal)) then
            plusWeight = 0.0_dp
        else

            plusWeight = dat%F * dat%plus + dat%G * (1.0_dp - dat%F) * dat%minus + &
                         nSwitches / bVal * (dat%G * dat%minus + dat%F * (1.0_dp - dat%G) * dat%plus)

        end if
        ASSERT(plusWeight >= 0.0_dp)

    end function getPlus_overlapLowering

    subroutine actHamiltonian(ilut, csf_i, excitations, nTot, t_singles_only, &
            t_print_time, t_full)
        ! subroutine to calculate the action of the full Hamiltonian on a
        ! a single CSF given in ilut bit representation and outputs a list
        ! of excitations also in ilut format, where the exact matrix element
        ! is stored, where usually the signed walker weight of an occupied
        ! determinant is stored
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nTot
        logical, intent(in), optional :: t_singles_only, t_print_time, t_full
        character(*), parameter :: this_routine = "actHamiltonian"
        type(timer), save :: proc_timer

        integer(n_int), allocatable :: tmp_all_excits(:, :)
        integer :: i, j, k, l, nExcits, nMax, ierr
        integer(n_int), allocatable :: tempExcits(:,:)
        logical :: t_singles_only_, t_print_time_, t_full_
        integer :: n
        real(dp) :: cmp

        def_default(t_singles_only_, t_singles_only, .false.)
        def_default(t_print_time_, t_print_time, .false.)
        def_default(t_full_, t_full, .true.)

        if (.not. present(t_singles_only)) then
            t_singles_only_ = .false.
        else
            t_singles_only_ = t_singles_only
        end if

        ASSERT(isProperCSF_ilut(ilut))

        ! finally add some timing routines to this as it gets quite heavy
        proc_timer%timer_name = this_routine
        call set_timer(proc_timer)

        ! essentially have to calculate the diagonal matrix element of ilut
        ! and additionally loop over all possible single excitaiton (i,j)
        ! and double excitation (i,j,k,l) index combinations and for them
        ! calculate all possible excitations and matrix elements

        ! have to somehow store all the already created excitations,...
        ! so I have to initialize an array somehow... how exactly am i going
        ! to do that... the ilut list combinig routine says there should be no
        ! repetitions in the list, but if i initialize a huge array to zero
        ! and fill that up, does that cause problems? -> ask simon/ali

        ! and how big does this array have to be?
        ! for each excitation 2**nOpenOrbs and there are at worst
        ! (nOrbitals/2)**4 excitation, but not all with the full excitation
        ! range.. but for now initiate it with this worse than worst case
        ! hm it seems i can give a number of determinants input to
        ! routine add_ilut_lists, so maybe i constrain the number of
        ! processed states with that
        nMax = 6 + 4 * (nSpatOrbs)**4 * (count_open_orbs(ilut) + 1)
        ! todo: get an exact maximum formula for that which holds in all cases!
        ! because that number above is WAYY too big!
        ! also for excitations where the indices are known!
        ! and fot those pescy non-overlap excitations!

        allocate(tmp_all_excits(0:nifguga, nMax), stat=ierr)
        call LogMemAlloc('tmp_all_excits', (nifguga + 1) * nMax, 8, this_routine, tag_tmp_excits)

        ! maybe have to set to zero here then, or initialize somehow different

        nTot = 0

        ! single excitations:
        ! does it help to not calculate stuff if not necessary..
        ! probably yes..
        if (.not. ((thub .and. .not. treal) .or. t_heisenberg_model)) then
        do i = 1, nSpatOrbs
            do j = 1, nSpatOrbs
                if (i == j) cycle ! do not calc. diagonal elements here
                ! need integral contributions too...
                ! maybe only loop over non-zero index combinations or
                ! otherwise a lot of effort wasted ...
                ! and i might have to loop over created excitations here too
                ! to multiply with the integrals... better to do that
                ! within the excitation creation!

                ! have to think where and when to allocate excitations
                ! or maybe call the routine below with a dummy variable.
                ! and store calculated ones in a different variable here..
                call calcAllExcitations(ilut, csf_i, i, j, tempExcits, nExcits, t_full_)
#ifdef DEBUG_
                do n = 1, nExcits
                    ASSERT(isProperCSF_ilut(tempExcits(:, n), .true.))
                end do
#endif

                ! merge created lists.. think how that will be implemented
                if (nExcits > 0) then
                    call add_guga_lists(nTot, nExcits, tmp_all_excits, tempExcits)
                    ! nTot gets automatically updated too
                end if

            end do
        end do
        end if

        ! double excitations
        ! do it really primitive for now. -> make it more elaborate later
        if (.not. t_singles_only_) then
        if (.not. (thub .and. treal)) then
        do i = 1, nSpatOrbs
            do j = 1, nSpatOrbs
                do k = 1, nSpatOrbs
                    do l = 1, nSpatOrbs
                        if (i == j .and. k == l) cycle

                        ! if (t_heisenberg_model) then
                        !     if (.not. (i == l .and. j == k)) cycle
                        ! end if
                        ! not only i=j and k=l index combinations can lead to
                        ! diagonal contributions but also i = l and k = j
                        ! mixed fullstart-> fullstop excitations can have a
                        ! diagonal term -> not only can, but always will have
                        ! so i have to exclude this term in the exact excitation
                        ! calculation when states get added up, since otherwise
                        ! i would count these terms to often, as they are
                        ! already taken into account in the diagonal term
                        ! calculation

                        ! integral contributions
                        call calcAllExcitations(ilut, csf_i, i, j, k, l, tempExcits, &
                                                nExcits, t_full_)

#ifdef DEBUG_
                        do n = 1, nExcits
                            if (.not. isProperCSF_ilut(tempExcits(:, n), .true.)) then
                                call write_guga_list(stdout, tempExcits(:, 1:nExcits))
                                print *, "ijkl:", i, j, k, l
                            end if
                            ASSERT(isProperCSF_ilut(tempExcits(:, n), .true.))
                        end do
#endif

                        if (i == l .and. k == j) then
                            ! exclude the possible diagonal term
                            ! i assume for now that this term is always there
                            ! atleast the x0-matrix element is never zero for
                            ! these diagonal excitations -> so check if it
                            ! stays in the same position: yes it stays in
                            ! first position always..
                            ! for non-open orbitals i
                            ! NO assert since i exclude excitations where
                            ! x1 is always 0
                            ! or otherwise somthings wrong... and then only
                            ! update the list if another excitation is
                            ! encountered
                            if (nExcits > 1) then

                                call add_guga_lists(nTot, nExcits - 1, tmp_all_excits, &
                                                    tempExcits(:, 2:))

                            end if
                        else
                            if (nExcits > 0) then
                                call add_guga_lists(nTot, nExcits, tmp_all_excits, &
                                                    tempExcits)
                            end if

                        end if

                    end do
                end do
            end do
        end do
        end if
        end if

        ! do an additional check if matrix element is 0

        j = 1
        if (t_matele_cutoff) then
            cmp = matele_cutoff
        else
            cmp = EPS
        end if

        do i = 1, nTot
            if (abs(extract_h_element(tmp_all_excits(:, i))) < cmp) cycle

            tmp_all_excits(:, j) = tmp_all_excits(:, i)

            j = j + 1

        end do

        nTot = j - 1

        allocate(excitations(0:nifguga, nTot), stat=ierr)
        ! hm to log that does not make so much sense.. since it gets called
        ! more than once and is only a temporary array..
        call LogMemAlloc('excitations', nTot, 8, this_routine, tag_excitations)

        excitations = tmp_all_excits(:, 1:nTot)

        deallocate(tmp_all_excits)
        call LogMemDealloc(this_routine, tag_tmp_excits)

        ! allocate the outputted excitations here since otherwise i store
        ! a way too big array than necessary

        ! i probably should resize the excitaions array.. since usually nMax
        ! is way bigger then nTot..

        call halt_timer(proc_timer)
        if (t_print_time_) then
            write(stdout, *) " Exact Hamiltonian application done! "
            write(stdout, *) " Elapsed time: ", get_total_time(proc_timer)
            call neci_flush(stdout)
        end if

    end subroutine actHamiltonian

    subroutine calcAllExcitations_excitInfo_single(ilut, csf_i, excitInfo, posSwitches, &
               negSwitches, tmatFlag, excitations, nExcits)
        ! excitation calculation if excitInfo is already calculated
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        logical, intent(in) :: tmatFlag
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits

        HElement_t(dp) :: tmat
        integer :: iOrb
        type(WeightObj_t) :: weights
        integer(n_int), allocatable :: tempExcits(:, :)

        ! to do combine (i,j) version and this version to one

        ! if tmatFlag is false do not check or include the tmat element
        ! for double excitations which essentially only involve sinlge excits

        ! i think this is always called with the flag as false... hm,
        ! maybe redundant and remove.
        ! when is this called actually? because this influences if i can use
        ! the csf_i%stepvector info and stuff..
        ! i am pretty sure that in all the calls to this routine
        ! the necessary quantities are known!
        if (tmatFlag) then
            tmat = getTmatEl(2 * excitInfo%i, 2 * excitInfo%j)
            if (near_zero(tmat)) then
                nExcits = 0
                allocate(excitations(0, 0))
                return
            end if
        else
            tmat = 1.0_dp
        end if

        ! also do not need to check if excitation is compatible, since this
        ! has already been done
        weights = init_singleWeight(csf_i, excitInfo%fullEnd)

        ! have to give probabilistic weight object as input, to deal
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, &
                               negSwitches, weights, tempExcits, nExcits)

        ! to not call getTmatEl again in createSingleStart loop over
        ! the atmost two excitations here and multiply with tmat
        do iOrb = excitInfo%fullStart + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, &
                              negSwitches, weights, tempExcits, nExcits)
        end do

        call singleEnd(ilut, csf_i, excitInfo, tempExcits, &
                       nExcits, excitations)

    end subroutine calcAllExcitations_excitInfo_single

    subroutine calcAllExcitations_single(ilut, csf_i, i, j, excitations, nExcits, t_full)
        ! function to calculate all possible single excitation for a CSF
        ! given in (ilut) format and indices (i,j). used to calculated
        ! H|D> to calculate the projected energy.
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: i, j
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        logical, intent(in), optional :: t_full
        character(*), parameter :: this_routine = "calcAllExcitations_single_ilut"

        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: tempExcits(:, :)
        integer :: ierr, iOrb, iEx, st
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs), &
                    plusWeight, minusWeight
        HElement_t(dp) :: tmat
        type(WeightObj_t) :: weights
        logical :: t_full_
        def_default(t_full_, t_full, .true.)

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        ! set nExcits to zero by default and only change when excitations
        ! happen
        nExcits = 0

        ! first check the single particle matrix element, it it is zero leave
        ! have index it with spin orbitals: assume non-UHF basis
        if (t_full_) then
            tmat = GetTMatEl(2 * i, 2 * j)
        else
            tmat = h_cast(1.0_dp)
        end if

        if (near_zero(tmat)) then
            allocate(excitations(0, 0), stat=ierr)
            return
        end if

        ! first determine excitation info:
        excitInfo = excitationIdentifier(i, j)

        ! determine the main restrictions of stepvalues depending on generator
        ! type -> no raising start at d = 3, no lowering start at d = 0 etc
        if (excitInfo%gen1 /= 0) then
            if (csf_i%stepvector(i) == 3 .or. csf_i%stepvector(j) == 0) then
                allocate(excitations(0, 0), stat=ierr)
                return
            end if

            if (t_tJ_model) then
                ! restrict hops to singly occupied orbitals in t-J model
                if (csf_i%Occ_int(i) == 1) then
                    allocate(excitations(0, 0), stat=ierr)
                    return
                end if
            end if
        end if

        ! also check necessary ending restrictions and probability weights
        ! to check if no excitation is possible
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitches, negSwitches)

        ! change it here to also use the functions involving
        ! the WeightObj_t objects.. to efficiently determine
        ! if excitations have to aborted
        ! and to make the whole code more general
        weights = init_singleWeight(csf_i, excitInfo%fullEnd)
        plusWeight = weights%proc%plus(posSwitches(excitInfo%fullStart), &
                                       csf_i%B_real(excitInfo%fullStart), weights%dat)
        minusWeight = weights%proc%minus(negSwitches(excitInfo%fullStart), &
                                         csf_i%B_real(excitInfo%fullStart), weights%dat)

        ! calc total weight functions of all possible branches
        ! if that is zero no excitation possible..
        ! do not need exact weight, but only knowledge if it is zero...
        ! so I do not yet need the bVector

        ! now we have to really calculate the excitations
        ! but also use the stochastic weight functions to not
        ! calculate unnecessary excitations...
        ! maybe here allocate arrays with the maximum number of
        ! possible excitations, and fill it up step after step
        ! 2^|i-j| is more then the maximum number of excitations..
        ! this could be a problem with this kind of implementation
        ! actually... -> since for every possible number of
        ! i,j index pair( and for double even i,j,k,l) this gets
        ! out of hand for big systems. maybe i have to switch back
        ! to determine the kind of excitations between to arbitrarily
        ! given CSFs and then calc. the overlap between them...
        ! wait and talk to Ali about that...

        ! should calculate bVector beforehand, and maybe store whole
        ! b vector for all excitation calculation (i,j) for a given
        ! CSF, since always the same and needed...

        ! when i use the remaining switches and end restrictions
        ! correctly i should not need to delete already created
        ! excitations. -> so i can fill up the excitation list step
        ! by step. -> and maybe use top most value as indication of
        ! delta b value... so i dont need two lists or some more
        ! advanced data-structure... but when i do it ilut format
        ! maybe i can use the already provided flag structure which
        ! is usually used for the sign of the walkers on a given CSF.
        ! have to figure out how to access and effectively adress
        ! them

        ! do it again in this kind of fashion:

        st = excitInfo%fullStart
        ! check compatibility of chosen indices

        if ((csf_i%stepvector(st) == 1 .and. near_zero(plusWeight)) .or. &
            (csf_i%stepvector(st) == 2 .and. near_zero(minusWeight)) .or. &
            near_zero(minusWeight + plusWeight)) then
            allocate(excitations(0, 0), stat=ierr)
            return
        end if

        ! have to give probabilistic weight object as input, to deal
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, &
                               negSwitches, weights, tempExcits, nExcits)

        ! to not call getTmatEl again in createSingleStart loop over
        ! the atmost two excitations here and multiply with tmat

        do iOrb = excitInfo%fullStart + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, &
                              negSwitches, weights, tempExcits, nExcits)
        end do

        call singleEnd(ilut, csf_i, excitInfo, tempExcits, &
                       nExcits, excitations)

        ! encode IC = 1 in the deltB information of the GUGA
        ! excitation to handle it in the remaining NECI code
        ! correctly
        do iEx = 1, nExcits
            call encode_matrix_element(excitations(:, iEx), 0.0_dp, 2)
            call update_matrix_element(excitations(:, iEx), tmat, 1)
            call setDeltaB(1, excitations(:, iEx))
        end do

    end subroutine calcAllExcitations_single

    subroutine doubleUpdate(ilut, csf_i, sO, excitInfo, weights, tempExcits, nExcits, &
                            negSwitches, posSwitches)
        ! for double excitation weights are usually always already calculated
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sO
        type(ExcitationInformation_t), intent(in) :: excitInfo
        type(WeightObj_t), intent(in) :: weights
        integer(n_int), intent(inout) :: tempExcits(:, :)
        integer, intent(inout) :: nExcits
        real(dp), intent(in) :: negSwitches(nSpatOrbs), posSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "doubleUpdate"

        real(dp) :: minusWeight, plusWeight, zeroWeight
        integer :: gen1, gen2, iEx, deltaB
        real(dp) :: bVal, tempWeight_0, tempWeight_1, order
        integer(n_int) :: t(0:nifguga), s(0:nifguga)

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(sO > 0 .and. sO <= nSpatOrbs)

        if (csf_i%Occ_int(sO) /= 1) then
            ! in this case no change in stepvector or matrix element
            return
        end if

        ! and for more readibility extract certain values:
        gen1 = excitInfo%gen1
        gen2 = excitInfo%gen2
        bVal = csf_i%B_real(sO)
        ! stupid!! order only needed at semistarts and semistops! so just set
        ! it to 1.0 here
        order = 1.0_dp

        ! else extract the relevant weights:
        minusWeight = weights%proc%minus(negSwitches(sO), bVal, weights%dat)
        plusWeight = weights%proc%plus(posSwitches(sO), bVal, weights%dat)
        zeroWeight = weights%proc%zero(negSwitches(sO), posSwitches(sO), bVal, &
                                       weights%dat)

        ! try new implementation with checking all 3 weights...
        ! ===================================================================
        if (csf_i%stepvector(sO) == 1) then
            if (minusWeight > 0.0_dp .and. plusWeight > 0.0_dp .and. zeroWeight > 0.0_dp) then
                ! all excitations are possible in this case.
                ! first do the on track
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    ! need it two times for matrix elements
                    s = t
                    deltaB = getDeltaB(t)

                    ! just calc matrix element for 1 -> 1 on track
                    call getDoubleMatrixElement(1, 1, deltaB, gen1, gen2, &
                                                bVal, order, tempWeight_0, tempWeight_1)

                    call update_matrix_element(t, tempWeight_0, 1)

                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t

                    ! if deltaB != 2 do :
                    if (deltaB /= 2) then
                        ! then do the 1 -> 2 switch
                        clr_orb(s, 2 * sO - 1)
                        set_orb(s, 2 * sO)

                        ! change deltaB
                        call setDeltaB(deltaB + 2, s)

                        ! calc matrix elements, only x1 matrix elements here
                        call getDoubleMatrixElement(2, 1, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                        call encode_matrix_element(s, 0.0_dp, 1)
                        call update_matrix_element(s, tempWeight_1, 2)

                        nExcits = nExcits + 1

                        tempExcits(:, nExcits) = s

                    end if
                end do

            else if (near_zero(plusWeight) .and. (.not. near_zero(minusWeight)) &
                     .and. (.not. near_zero(zeroWeight))) then
                ! all purely + branches are not possible
                ! do the on track possibles first
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    s = t
                    deltaB = getDeltaB(t)
                    ! i think i never should have a +2 here then now...
                    ASSERT(deltaB /= 2)

                    call getDoubleMatrixElement(1, 1, deltaB, gen1, gen2, &
                                                bVal, order, tempWeight_0, tempWeight_1)

                    call update_matrix_element(t, tempWeight_0, 1)
                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t

                    ! and do the 1->2 to 0 branch in case of dB = -2
                    if (deltaB == -2) then
                        ! 1 -> 2 switch
                        clr_orb(s, 2 * sO - 1)
                        set_orb(s, 2 * sO)

                        ! change deltaB
                        call setDeltaB(deltaB + 2, s)

                        ! calc matrix elements, only x1 matrix elements here
                        call getDoubleMatrixElement(2, 1, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                        call encode_matrix_element(s, 0.0_dp, 1)
                        call update_matrix_element(s, tempWeight_1, 2)

                        nExcits = nExcits + 1

                        tempExcits(:, nExcits) = s
                    end if
                end do

            else if (near_zero(minusWeight) .and. (.not. near_zero(plusWeight)) &
                     .and. (.not. near_zero(zeroWeight))) then

                ! -2 branch not possible.. not sure if i ever come here with
                ! dB = -2 then.. -> check
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    s = t
                    deltaB = getDeltaB(t)
                    ! do the on-track specifically if a -2 comes
                    if (deltaB == -2) then
                        clr_orb(t, 2 * sO - 1)
                        set_orb(t, 2 * sO)

                        call setDeltaB(0, t)

                        call getDoubleMatrixElement(2, 1, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                        call encode_matrix_element(t, 0.0_dp, 1)
                        call update_matrix_element(t, tempWeight_1, 2)

                        tempExcits(:, iEx) = t

                    else
                        ! elsewise do the staying possib
                        call getDoubleMatrixElement(1, 1, deltaB, gen1, gen2, &
                                                    bVal, order, tempWeight_0, tempWeight_1)

                        call update_matrix_element(t, tempWeight_0, 1)
                        call update_matrix_element(t, tempWeight_1, 2)

                        tempExcits(:, iEx) = t

                        ! and if dB = 0, do the possible +2 switch too
                        if (deltaB == 0) then
                            clr_orb(s, 2 * sO - 1)
                            set_orb(s, 2 * sO)

                            call setDeltaB(2, s)

                            call getDoubleMatrixElement(2, 1, deltaB, gen1, gen2, &
                                                        bVal, order, x1_element=tempWeight_1)

                            call encode_matrix_element(s, 0.0_dp, 1)
                            call update_matrix_element(s, tempWeight_1, 2)

                            nExcits = nExcits + 1

                            tempExcits(:, nExcits) = s
                        end if
                    end if
                end do
            else if (near_zero(minusWeight) .and. near_zero(plusWeight) &
                     .and. (.not. near_zero(zeroWeight))) then
                ! only the 0 branches are possible
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)

                    ASSERT(deltaB /= +2)

                    if (deltaB == -2) then
                        ! do the 1->2 0 swicth
                        clr_orb(t, 2 * sO - 1)
                        set_orb(t, 2 * sO)

                        call setDeltaB(0, t)

                        call getDoubleMatrixElement(2, 1, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)
                        tempWeight_0 = 0.0_dp

                    else
                        ! if 0 comes just get the matrix elements
                        call getDoubleMatrixElement(1, 1, deltaB, gen1, gen2, &
                                                    bVal, order, tempWeight_0, tempWeight_1)

                    end if
                    ! and update matrix elements
                    call update_matrix_element(t, tempWeight_0, 1)
                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t
                end do
            else if ((.not. near_zero(plusWeight)) .and. (.not. near_zero(minusWeight)) &
                     .and. near_zero(zeroWeight)) then
                ! the continuing 0-branches are not allowed
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)

                    if (deltaB == 0) then
                        ! do the switch i a 0 branch comes
                        clr_orb(t, 2 * sO - 1)
                        set_orb(t, 2 * sO)

                        call setDeltaB(2, t)

                        call getDoubleMatrixElement(2, 1, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                    else
                        ! otherwise just get the matrix elemet
                        ! x0-always zero in this case
                        call getDoubleMatrixElement(1, 1, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)
                    end if

                    call encode_matrix_element(t, 0.0_dp, 1)
                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t
                end do
            else if (near_zero(zeroWeight) .and. near_zero(plusWeight) &
                     .and. (.not. near_zero(minusWeight))) then
                ! only -2 staying branch is possible... so i should just be here
                ! with a -2
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    ASSERT(getDeltaB(t) == -2)

                    call getDoubleMatrixElement(1, 1, -2, gen1, gen2, &
                                                bVal, order, x1_element=tempWeight_1)

                    ! x0 element should already be 0
                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t
                end do
            else if (near_zero(zeroWeight) .and. near_zero(zeroWeight) &
                     .and. (.not. near_zero(plusWeight))) then
                ! only the +2 cont. branches are possible
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)

                    ASSERT(deltaB /= -2)

                    if (deltaB == 0) then
                        ! do the switch to +2
                        clr_orb(t, 2 * sO - 1)
                        set_orb(t, 2 * sO)

                        call setDeltaB(2, t)

                        call getDoubleMatrixElement(2, 1, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                    else
                        call getDoubleMatrixElement(1, 1, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)
                    end if
                    call encode_matrix_element(t, 0.0_dp, 1)
                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t
                end do
            else
                ! should not be here... all weights 0#
                call stop_all(this_routine, "should not be here! all 3 weights 0 in double update at sO=1")
            end if

            ! also do the possibilities at a step=2
        else if (csf_i%stepvector(sO) == 2) then
            if (minusWeight > 0.0_dp .and. plusWeight > 0.0_dp .and. zeroWeight > 0.0_dp) then
                ! everything possible!
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    s = t

                    deltaB = getDeltaB(t)

                    ! do on track first
                    call getDoubleMatrixElement(2, 2, deltaB, gen1, gen2, &
                                                bVal, order, tempWeight_0, tempWeight_1)

                    call update_matrix_element(t, tempWeight_0, 1)
                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t

                    ! if /= -2 2 -> 1 branching also possible
                    if (deltaB /= -2) then
                        clr_orb(s, 2 * sO)
                        set_orb(s, 2 * sO - 1)

                        call setDeltaB(deltaB - 2, s)

                        call getDoubleMatrixElement(1, 2, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                        call encode_matrix_element(s, 0.0_dp, 1)
                        call update_matrix_element(s, tempWeight_1, 2)

                        nExcits = nExcits + 1
                        tempExcits(:, nExcits) = s
                    end if
                end do
            else if (near_zero(plusWeight) .and. (.not. near_zero(minusWeight)) &
                     .and. (.not. near_zero(zeroWeight))) then
                ! +2 cont. not possible
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)
                    if (deltaB == -2) then
                        call getDoubleMatrixElement(2, 2, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                        call update_matrix_element(t, tempWeight_1, 2)

                        tempExcits(:, iEx) = t

                    else if (deltaB == 2) then

                        clr_orb(t, 2 * sO)
                        set_orb(t, 2 * sO - 1)

                        call setDeltaB(0, t)

                        call getDoubleMatrixElement(1, 2, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                        call update_matrix_element(t, tempWeight_1, 2)

                        tempExcits(:, iEx) = t

                    else
                        s = t
                        call getDoubleMatrixElement(2, 2, deltaB, gen1, gen2, &
                                                    bVal, order, tempWeight_0, tempWeight_1)

                        call update_matrix_element(t, tempWeight_0, 1)
                        call update_matrix_element(t, tempWeight_1, 2)

                        tempExcits(:, iEx) = t

                        ! then do switch
                        clr_orb(s, 2 * sO)
                        set_orb(s, 2 * sO - 1)

                        call setDeltaB(-2, s)

                        call getDoubleMatrixElement(1, 2, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                        call encode_matrix_element(s, 0.0_dp, 1)
                        call update_matrix_element(s, tempWeight_1, 2)

                        nExcits = nExcits + 1

                        tempExcits(:, nExcits) = s

                    end if
                end do
            else if (near_zero(minusWeight) .and. (.not. near_zero(plusWeight)) &
                     .and. (.not. near_zero(zeroWeight))) then
                ! cont. -2 branches not possible
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)

                    s = t
                    ASSERT(deltaB /= -2)

                    ! do staying first

                    call getDoubleMatrixElement(2, 2, deltaB, gen1, gen2, &
                                                bVal, order, tempWeight_0, tempWeight_1)

                    call update_matrix_element(t, tempWeight_0, 1)
                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t

                    ! then do possible switch
                    if (deltaB == 2) then
                        clr_orb(s, 2 * sO)
                        set_orb(s, 2 * sO - 1)

                        call setDeltaB(0, s)

                        call getDoubleMatrixElement(1, 2, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                        call update_matrix_element(s, tempWeight_1, 2)

                        nExcits = nExcits + 1

                        tempExcits(:, nExcits) = s

                    end if
                end do
            else if (near_zero(minusWeight) .and. near_zero(plusWeight) &
                     .and. (.not. near_zero(zeroWeight))) then
                ! only cont. 0 branch valid
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)
                    ASSERT(deltaB /= -2)

                    if (deltaB == 2) then
                        ! do switch
                        clr_orb(t, 2 * sO)
                        set_orb(t, 2 * sO - 1)

                        call setDeltaB(0, t)

                        call getDoubleMatrixElement(1, 2, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)
                        tempWeight_0 = 0.0_dp
                    else
                        call getDoubleMatrixElement(2, 2, deltaB, gen1, gen2, &
                                                    bVal, order, tempWeight_0, tempWeight_1)

                    end if

                    call update_matrix_element(t, tempWeight_0, 1)
                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t
                end do
            else if ((.not. near_zero(plusWeight)) .and. (.not. near_zero(minusWeight)) &
                     .and. near_zero(zeroWeight)) then
                ! no 0 branch valid
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)
                    ! only switch if 0 branch arrives
                    if (deltaB == 0) then
                        ! do swicht
                        clr_orb(t, 2 * sO)
                        set_orb(t, 2 * sO - 1)

                        call setDeltaB(-2, t)

                        call getDoubleMatrixElement(1, 2, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)
                        call encode_matrix_element(t, 0.0_dp, 1)
                    else
                        ! staying is possible

                        call getDoubleMatrixElement(2, 2, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                    end if
                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t
                end do
            else if (near_zero(zeroWeight) .and. near_zero(plusWeight) &
                     .and. (.not. near_zero(minusWeight))) then
                ! only -2 branches valis
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)
                    ASSERT(deltaB /= 2)
                    if (deltaB == 0) then
                        clr_orb(t, 2 * sO)
                        set_orb(t, 2 * sO - 1)

                        call setDeltaB(-2, t)

                        call getDoubleMatrixElement(1, 2, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)
                        call encode_matrix_element(t, 0.0_dp, 1)
                    else
                        call getDoubleMatrixElement(2, 2, deltaB, gen1, gen2, &
                                                    bVal, order, x1_element=tempWeight_1)

                    end if
                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t
                end do
            else if (near_zero(zeroWeight) .and. near_zero(minusWeight) &
                     .and. (.not. near_zero(plusWeight))) then
                ! only +2 staying branch is possible -> assert that only that
                ! comes
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    ASSERT(getDeltaB(t) == 2)

                    call getDoubleMatrixElement(2, 2, 2, gen1, gen2, &
                                                bVal, order, x1_element=tempWeight_1)

                    call update_matrix_element(t, tempWeight_1, 2)

                    tempExcits(:, iEx) = t
                end do
            else
                ! should not be here.. all 3 weights 0
                call stop_all(this_routine, "should not be here. all 3 weights 0 at s=2 in!")
            end if
        end if

    end subroutine doubleUpdate

    subroutine singleUpdate(ilut, csf_i, sOrb, excitInfo, posSwitches, negSwitches, &
                            weightObj, tempExcits, nExcits)
        ! update function for calculation of all single excitations of a
        ! given CSF ilut. this implementation takes a general
        ! probabilistic weight obj as input, to calculate the specific
        ! probabilisitc weights, determining the different needed ones for
        ! specific excitations
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t) :: weightObj
        integer(n_int), intent(inout) :: tempExcits(:, :)
        integer, intent(inout) :: nExcits
        character(*), parameter :: this_routine = "singleUpdate"

        integer(n_int) :: t(0:nifguga)
        integer :: iEx, deltaB, gen
        real(dp) :: plusWeight, minusWeight, tempWeight, bVal

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(sOrb > 0 .and. sOrb < nSpatOrbs)

        select case (csf_i%stepvector(sOrb))
!         if (csf_i%stepvector(sOrb) == 0) then
!         if (isZero(ilut, sOrb)) then
            ! do nothin actually.. not even change matrix elements
        case (0)

            return

!         else if (csf_i%stepvector(sOrb) == 3) then
!         else if (isThree(ilut, sOrb)) then
            ! only change matrix element to negative one
        case (3)

            do iEx = 1, nExcits
                call update_matrix_element(tempExcits(:, iEx), -1.0_dp, 1)
            end do
            return

        end select
!         end if

        ! to generally use this function i need to define a current generator...
        gen = excitInfo%currentGen
        bVal = csf_i%B_real(sOrb)
        ! if 1 or 2 as stepvalue need probWeight
        ! here the change to the regular excitations is made:
        ! smth like: have yet to determine, when weight object gets initialized
        ! but have to check b-value here to avoid division by 0
        ! or include that in the probWeight calculation?
        ! b can only be 0 after a 2. which is only a switch possib. for +1
        ! in the single excitation atleast. and if i set the plus weight to
        ! 0 and the -1 weight depending on the switch possibilities and the
        ! end value i could somehow include that?!
        ! have the bValue restriction now included in the weight calc.
        ! atleast for pure single excitations... have to additionally do that
        ! for the more complicated excitation types too
        plusWeight = weightObj%proc%plus(posSwitches(sOrb), bVal, weightObj%dat)
        minusWeight = weightObj%proc%minus(negSwitches(sOrb), bVal, weightObj%dat)

        ! have to do some sort of abort_excitations funciton here too
        ASSERT(.not. near_zero(plusWeight + minusWeight))

        if (csf_i%stepvector(sOrb) == 1) then
            ! if its a deltaB = -1 this is a switch possib, indepentent
            ! of the b-vector..
            ! have to loop over already created excitations

            ! maybe check if weight is positive before loop, as it is always
            ! the same and then do an according update..
            ! is positive weight is 0, negative weight has to be >0
            ! or else we wouldn be here -> no switches just update -1 branches
            if (near_zero(plusWeight)) then
                ! no switches lead to  a nonzero excitation, just update
                ! matrix element and stay on track
                do iEx = 1, nExcits
                    ! need deltaB value
#ifdef DEBUG_
                    deltaB = getDeltaB(tempExcits(:, iEx))
                    ! check if a -1 is encountert
                    ASSERT(deltaB == -1)
#endif

                    call update_matrix_element(tempExcits(:, iEx), &
                                               getSingleMatrixElement(1, 1, -1, gen, bVal), 1)

                end do
                ! when negative weight is 0, positiv weight has to be > 0
                ! so update positive branches and switch negative ones
            else if (near_zero(minusWeight)) then
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)

                    if (deltaB == 1) then
                        ! only encode new matrix element in this case
                        call update_matrix_element(t, &
                                                   getSingleMatrixElement(1, 1, deltaB, gen, bVal), 1)

                    else
                        ! switch also 1 -> 2
                        clr_orb(t, 2 * sOrb - 1)
                        set_orb(t, 2 * sOrb)

                        call setDeltaB(1, t)

                        call update_matrix_element(t, &
                                                   getSingleMatrixElement(2, 1, deltaB, gen, bVal), 1)

                    end if
                    tempExcits(:, iEx) = t
                end do
                ! both weights are positiv, staying on track and switching possib
            else
                ! check if deltaB=-1
                do iEx = 1, nExcits
                    ! staying on track possible for all excitations. calculate
                    deltaB = getDeltaB(tempExcits(:, iEx))

                    tempWeight = extract_matrix_element(tempExcits(:, iEx), 1)

                    call encode_matrix_element(tempExcits(:, iEx), tempWeight * &
                                               getSingleMatrixElement(1, 1, deltaB, gen, bVal), 1)

                    if (deltaB == -1) then
                        ! for deltaB=-1 branch a 1 is a switch possibility, if
                        ! weight of positive branch is non-zero.
                        ! new excitations get appended at the end of tempExcits
                        nExcits = nExcits + 1
                        t = tempExcits(:, iEx)
                        ! have to change corresponding bits
                        clr_orb(t, 2 * sOrb - 1)
                        set_orb(t, 2 * sOrb)

                        ! update matrix element and deltaB flag on new branch
                        call setDeltaB(1, t)
                        ! forgot the decision with which deltaB i should access
                        call encode_matrix_element(t, tempWeight * &
                                                   getSingleMatrixElement(2, 1, deltaB, gen, bVal), 1)

                        ! and fill into list with correct deltaB
                        tempExcits(:, nExcits) = t

                    end if
                end do
            end if

        else
            ! if it is a 2, it gets a bit more complicated, with the forced
            ! switches and non-possible stayings...
            ! if it is a forcesd switch but the negative weight is zero
            ! i should set the matrix element to zero and move on
            ! if the weight is not zero but its still a forced switch
            ! i should update it on the spot, and not add any additional
            ! excitations. but if switch is possible and weights is
            ! bigger then zero do it as above

            ! if b < 1 there is always a forced switch for +1 branch,
            ! but this situation is not covered by my probWeights...
            ! hence here i have to additionally check and it could lead to
            ! no excitations at all.

            ! did some bullshit thinking here....

            ! can change down here something now, since i changed the bValue
            ! inclusion in the weight functions..

            if (near_zero(plusWeight)) then
                ! update on spot and switch
                ! in this case all +1 branches HAVE to switch, but leave them
                ! in same position
                do iEx = 1, nExcits

                    t = tempExcits(:, iEx)

                    if (getDeltaB(tempExcits(:, iEx)) == 1) then

                        ! change stepvectors and branch

                        set_orb(t, 2 * sOrb - 1)
                        clr_orb(t, 2 * sOrb)

                        call setDeltaB(-1, t)

                        tempWeight = getSingleMatrixElement(1, 2, 1, gen, bVal)

                    else
                        ! else just update the matrix elements :
                        tempWeight = getSingleMatrixElement(2, 2, -1, gen, bVal)

                    end if

                    call update_matrix_element(t, tempWeight, 1)

                    tempExcits(:, iEx) = t
                end do

            else if (near_zero(minusWeight)) then
                ! update on spot stay
                ! in this case staying on branch is possible for +1 and a
                ! -1 branch would have a zero weight, so just update matrix
                ! but no additional excitations...

                do iEx = 1, nExcits
#ifdef DEBUG_
                    deltaB = getDeltaB(tempExcits(:, iEx))
!                     ASSERT(deltaB == +1)
#endif

                    call update_matrix_element(tempExcits(:, iEx), &
                                               getSingleMatrixElement(2, 2, 1, gen, bVal), 1)

                end do

            else if (minusWeight > 0.0_dp .and. plusWeight > 0.0_dp) then
                ! update stay on spot and add additional excitations
                do iEx = 1, nExcits
                    ! first update matrix elements from staying excitations
                    deltaB = getDeltaB(tempExcits(:, iEx))

                    tempWeight = extract_matrix_element(tempExcits(:, iEx), 1)

                    call update_matrix_element(tempExcits(:, iEx), &
                                               getSingleMatrixElement(2, 2, deltaB, gen, bVal), 1)

                    if (deltaB == 1) then
                        ! for deltaB +1 branches switch possibility
                        nExcits = nExcits + 1
                        ! change bits and branch and store at end
                        t = tempExcits(:, iEx)

                        set_orb(t, 2 * sOrb - 1)
                        clr_orb(t, 2 * sOrb)

                        call setDeltaB(-1, t)

                        call encode_matrix_element(t, tempWeight * &
                                                   getSingleMatrixElement(1, 2, deltaB, gen, bVal), 1)

                        tempExcits(:, nExcits) = t

                    end if
                end do
            else
                ! something went wrong in this case
                call stop_all(this_routine, &
                              "something went wrong in the excitation generation, shouldnt be here!")
            end if
        end if

    end subroutine singleUpdate

    subroutine singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)
        ! end function to calculate all single excitations of a given CSF
        ! ilut
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(inout), allocatable :: tempExcits(:, :)
        integer, intent(inout) :: nExcits
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        character(*), parameter :: this_routine = "singleEnd"
        integer :: iEx, cnt, ende, ierr, gen, deltaB
        integer(n_int) :: t(0:nifguga)
        real(dp) :: bVal, tempWeight

        ! with the correct and consistent use of the probabilistic weight
        ! functions during the excitation creation, and the assertment that
        ! the indices and provided CSF are compatible to provide possible
        ! excitation, no wrong value deltaB branches should arrive here..

        ! nevertheless do some asserts in debug mode

#ifdef DEBUG_
        ASSERT(excitInfo%currentGen /= 0)
        ASSERT(isProperCSF_ilut(ilut))
        if (excitInfo%currentGen == gen_type%R) then
            ASSERT(.not. isZero(ilut, excitInfo%fullEnd))
        else if (excitInfo%currentGen == gen_type%L) then
            ASSERT(.not. isThree(ilut, excitInfo%fullEnd))
        end if
#endif

        ! for easier readability
        ende = excitInfo%fullEnd
        ! hm... not sure how to treat the end b value... essentially need
        ! it from one orbital more than -> maybe jsut add or distract one
        ! if 2/1 step
        bVal = csf_i%B_real(ende)
        gen = excitInfo%currentGen

        ! although I do not think I have to check if deltaB values fit,
        ! still do it for now and check if its really always correct...
        select case (csf_i%stepvector(ende))
        case (0)
            ! if it is zero it implies i have a lowering generator, otherwise
            ! i wouldnt even be here.
            do iEx = 1, nExcits
                ! here depending on deltab b set to 1 or 2
                t = tempExcits(:, iEx)

                deltaB = getDeltaB(t)

                if (deltaB == -1) then
                    ! set alpha bit
                    set_orb(t, 2 * ende)
                    tempWeight = getSingleMatrixElement(2, 0, deltaB, gen, bVal)

                else
                    ! set beta bit
                    set_orb(t, 2 * ende - 1)

                    tempWeight = getSingleMatrixElement(1, 0, deltaB, gen, bVal)

                end if
                ! and fill it in with matrix element and maybe even store
                ! the calculated matrix element in the signed weight
                ! entry here... or update weights() -> design decision todo
                call update_matrix_element(t, tempWeight, 1)
                tempExcits(:, iEx) = t

            end do

        case (3)
            ! if d = 3, it implies a raising generator or otherwise we
            ! would not be here ..
            do iEx = 1, nExcits
                t = tempExcits(:, iEx)
                deltaB = getDeltaB(t)

                if (deltaB == -1) then
                    ! clear beta bit
                    clr_orb(t, 2 * ende - 1)

                    tempWeight = getSingleMatrixElement(2, 3, deltaB, gen, bVal)

                else
                    ! clear alpha bit
                    clr_orb(t, 2 * ende)

                    tempWeight = getSingleMatrixElement(1, 3, deltaB, gen, bVal)

                end if

                call update_matrix_element(t, tempWeight, 1)

                tempExcits(:, iEx) = t

            end do

        case (1)
            ! d=1 needs a deltaB -1 branch, so delete +1 excitations if they
            ! happen to get here. alhough that shouldnt happen with the correct
            ! use of probabilistic weight functions.
            cnt = 1

            ! have to check generator type in this case, as both are possible
            ! and should do that outside of the loop as always the same...
            if (gen == gen_type%R) then
                do iEx = 1, nExcits
                    deltaB = getDeltaB(tempExcits(:, cnt))
                    if (deltaB == 1) then
                        ! insert last excitation in current place and delete one
                        ! excitation then cycle
                        tempExcits(:, cnt) = tempExcits(:, nExcits)
                        nExcits = nExcits - 1
                        ! also make noise for now to indicate smth didnt work
                        ! as expected
                        cycle
                    end if
                    ! otherwise deal normally with them
                    t = tempExcits(:, cnt)
                    ! for raising gen 1->0 : clear beta bit
                    clr_orb(t, 2 * ende - 1)

                    ! matrix element is just 1 so dont do anything

                    tempExcits(:, cnt) = t
                    cnt = cnt + 1
                end do
            else ! lowering generator
                do iEx = 1, nExcits
                    deltaB = getDeltaB(tempExcits(:, cnt))

                    if (deltaB == 1) then
                        tempExcits(:, cnt) = tempExcits(:, nExcits)
                        nExcits = nExcits - 1
                        cycle
                    end if
                    t = tempExcits(:, cnt)
                    ! for lowering gen: 1 -> 3: set alpha bit
                    set_orb(t, 2 * ende)

                    tempWeight = getSingleMatrixElement(3, 1, deltaB, gen, bVal)
                    call update_matrix_element(t, tempWeight, 1)

                    tempExcits(:, cnt) = t

                    cnt = cnt + 1
                end do
            end if

        case (2)
            cnt = 1

            if (gen == gen_type%R) then
                do iEx = 1, nExcits
                    deltaB = getDeltaB(tempExcits(:, cnt))
                    if (deltaB == -1) then
                        tempExcits(:, cnt) = tempExcits(:, nExcits)
                        nExcits = nExcits - 1
                        cycle
                    end if
                    t = tempExcits(:, cnt)
                    ! for raising: 2 -> 0: clear alpha bit
                    clr_orb(t, 2 * ende)

                    ! matrix element is just 1 -> so dont do anything
                    tempExcits(:, cnt) = t

                    cnt = cnt + 1
                end do
            else ! lowering gen
                do iEx = 1, nExcits
                    deltaB = getDeltaB(tempExcits(:, cnt))
                    if (deltaB == -1) then
                        tempExcits(:, cnt) = tempExcits(:, nExcits)
                        nExcits = nExcits - 1
                        cycle
                    end if
                    t = tempExcits(:, cnt)
                    ! for lowering gen: 2 -> 3 : set beta bit
                    set_orb(t, 2 * ende - 1)

                    tempWeight = getSingleMatrixElement(3, 2, deltaB, gen, bVal)
                    call update_matrix_element(t, tempWeight, 1)

                    tempExcits(:, cnt) = t

                    cnt = cnt + 1
                end do
            end if
        end select
        ! now have to cut tempExcits to only necessary
        ! and used entries and output it!

        cnt = 1
        ! fill up the excitations and store matrix element
        do iEx = 1, nExcits
            ! maybe check here again to not have excitations with zero matrix
            ! elements
            if (near_zero(extract_matrix_element(tempExcits(:, iEx), 1))) cycle

            tempExcits(:, cnt) = tempExcits(:, iEx)

            cnt = cnt + 1

        end do
        ! also update nExcits one final time if something gets cut
        nExcits = cnt - 1

        allocate(excitations(0:nifguga, nExcits), stat=ierr)

        excitations = tempExcits(:, 1:nExcits)

        ! and get rid of container variables
        deallocate(tempExcits)

    end subroutine singleEnd

    subroutine createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                                 weightObj, tempExcits, nExcits)
        ! subroutine to create full single excitations starts for ilut
        ! allocates the necessary arrays and fills up first excitations,
        ! depending on stepvalue in ilut, the corresponding b value, and
        ! probabilitstic weight functions(to exclude excitations eventually
        ! leading to zero weight)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        type(WeightObj_t), intent(in) :: weightObj
        integer(n_int), intent(out), allocatable :: tempExcits(:, :)
        integer, intent(out) :: nExcits
        character(*), parameter :: this_routine = "createSingleStart"
        integer :: ierr, nmax, st, gen
        integer(n_int) :: t(0:nifguga)
        real(dp) :: minusWeight, plusWeight, tempWeight, bVal

        ! have already asserted that start and end values of stepvector and
        ! generator type are consistent to allow for an excitation.
        ! maybe still assert here in debug mode atleast
#ifdef DEBUG_
        ! also assert we are not calling it for a weight gen. accidently
        ASSERT(excitInfo%currentGen /= 0)
        ASSERT(isProperCSF_ilut(ilut))
        ! also check if calculated b vector really fits to ilut
        if (excitInfo%currentGen == gen_type%R) then
            ASSERT(.not. isThree(ilut, excitInfo%fullStart))
        else if (excitInfo%currentGen == gen_type%L) then
            ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        end if
        ! also have checked if atleast on branch way can lead to an excitaiton
#endif

        ! for more oversight
        st = excitInfo%fullStart
        gen = excitInfo%currentGen
        bVal = csf_i%B_real(st)

        ! first have to allocate both arrays for the determinant and weight list
        ! worse then worst case for single excitations are 2^|i-j| excitations
        ! for a given CSF -> for now use that to allocate the arrays first
        ! update only need the number of open orbitals between i and j, and
        ! some additional room if 0/3 at start
        ! use already provided open orbital counting function.
        ! nMax = 2**(ende - start)
        nMax = 4 + 4 * 2**count_open_orbs_ij(csf_i, st, excitInfo%fullEnd)
        allocate(tempExcits(0:nifguga, nMax), stat=ierr)

        ! create start depending on stepvalue of ilut at start, b value,
        ! which is already calculated at the start of the excitation call and
        ! probabilistic weights
        ! for now just scratch it up in an inefficient way.. optimize later
        ! copy ilut into first row(or column?) of tempExcits
        tempExcits(:, 1) = ilut

        ! i think determinants get initiated with zero matrix element, but not
        ! sure ... anyway set matrix element to 1

        ! maybe need temporary ilut storage
        t = ilut
        nExcits = 0
        select case (csf_i%stepvector(st))
        case (1)
            ! set corresponding orbital to 0 or 3 depending on generator type
            if (gen == gen_type%R) then ! raising gen case
                ! set the alpha orbital also to 1 to make d=1 -> d'=3
                set_orb(t, 2 * st)

                ! does it work like that:
                ! would have all necessary values and if statements to directly
                ! calculated matrix element here... maybe do it here after all
                ! to be more efficient...
                tempWeight = getSingleMatrixElement(3, 1, +1, gen, bVal)

            else ! lowering gen case
                ! clear beta orbital to make d=1 -> d'=0
                clr_orb(t, 2 * st - 1)

                tempWeight = getSingleMatrixElement(0, 1, +1, gen, bVal)

            end if
            ! save number of already stored excitations
            nExcits = nExcits + 1
            ! store matrix element
            call encode_matrix_element(t, tempWeight, 1)
            ! set deltaB:
            ! for both cases deltaB = +1, set that as signed weight
            ! update: use previously unused flags to encode deltaB
            call setDeltaB(1, t)
            ! and store it in list:
            tempExcits(:, nExcits) = t

        case (2)
            if (gen == gen_type%R) then
                set_orb(t, 2 * st - 1)

                ! matrix elements
                tempWeight = getSingleMatrixElement(3, 2, -1, gen, bVal)

            else
                clr_orb(t, 2 * st)

                ! matrix elements
                tempWeight = getSingleMatrixElement(0, 2, -1, gen, bVal)

            end if
            nExcits = nExcits + 1
            call encode_matrix_element(t, tempWeight, 1)
            call setDeltaB(-1, t)
            tempExcits(:, nExcits) = t

        case (0, 3)
            ! if the deltaB = -1 weights is > 0 this one always works
            ! write weight calculation function! is a overkill here,
            ! but makes it easier to use afterwards with stochastic excitaions

            ! b value restriction now implemented ind the weights...
            ! -> plus weight is 0 if bValue == 0

            ! use new weight objects here.
            minusWeight = weightObj%proc%minus(negSwitches(st), bVal, weightObj%dat)
            plusWeight = weightObj%proc%plus(posSwitches(st), bVal, weightObj%dat)

            ! still have to do some abort excitation routine if both weights
            ! are 0
            ASSERT(.not. near_zero(minusWeight + plusWeight))

            if (minusWeight > 0.0_dp) then
                if (gen == gen_type%R) then
                    ! set beta bit
                    set_orb(t, 2 * st - 1)

                    ! matrix element
                    tempWeight = getSingleMatrixElement(1, 0, -1, gen, bVal)

                else
                    ! clear alpha bit
                    clr_orb(t, 2 * st)

                    ! matrix element
                    tempWeight = getSingleMatrixElement(1, 3, -1, gen, bVal)

                end if
                nExcits = nExcits + 1
                call encode_matrix_element(t, tempWeight, 1)
                call setDeltaB(-1, t)
                tempExcits(:, nExcits) = t

            end if

            ! if plusWeight and the bValue allows it also a deltaB=+1 branch
            ! dont need bVal check anymore, since already implemented in
            ! weight calc.
            if (plusWeight > 0.0_dp) then

                ! reset t
                t = ilut

                if (gen == gen_type%R) then
                    ! set alpha bit
                    set_orb(t, 2 * st)

                    ! matrix element
                    tempWeight = getSingleMatrixElement(2, 0, +1, gen, bVal)

                else
                    ! cleat beta bit
                    clr_orb(t, 2 * st - 1)

                    ! matrix element
                    tempWeight = getSingleMatrixElement(2, 3, +1, gen, bVal)
                end if
                ! update list and number
                nExcits = nExcits + 1
                call encode_matrix_element(t, tempWeight, 1)
                call setDeltaB(+1, t)
                tempExcits(:, nExcits) = t

            end if
        end select

    end subroutine createSingleStart

    elemental function endLx(csf_i, ind) result(ret)
        ! special function to determine if branching at the single overlap
        ! site of a RR excitation is possible
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: ind
        real(dp) :: ret

        if (csf_i%stepvector(ind) == 2 .and. csf_i%B_int(ind) > 1) then
            ret = 1.0_dp
        else
            ret = 0.0_dp
        end if
    end function endLx

    subroutine calcAllExcitations_double(ilut, csf_i, i, j, k, l, excitations, nExcits,&
            t_full)
        ! function to calculate all possible double excitation for a CSF
        ! given in (ilut) format and indices (i,j,k,l).
        ! used to calculate the action of the Hamiltonian H|D> to calculate
        ! the projected energy
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: i, j, k, l
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        logical, intent(in), optional :: t_full
        character(*), parameter :: this_routine = "calcAllExcitations_double"

        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        HElement_t(dp) :: umat
        integer :: ierr, n, exlevel
        type(ExcitationInformation_t) :: excitInfo
        logical :: compFlag, t_full_
        def_default(t_full_, t_full, .true.)

        ! if called with k = l = 0 -> call single version of function
        if (k == 0 .and. l == 0) then
            call calcAllExcitations(ilut, csf_i, i, j, excitations, nExcits, t_full_)
            return
        end if

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)
        ASSERT(k > 0 .and. k <= nSpatOrbs)
        ASSERT(l > 0 .and. l <= nSpatOrbs)
        ASSERT(isProperCSF_ilut(ilut))

        ! default:
        nExcits = 0

        ! first check two-particle integral
        if (t_full_) then
            umat = get_umat_el(i, k, j, l)
        else
            umat = h_cast(1.0_dp)
        end if

        if (near_zero(umat)) then
            allocate(excitations(0, 0), stat=ierr)
            return
        end if

        ! otherwise get excitation information
        excitInfo = excitationIdentifier(i, j, k, l)

        ! screw it. for now write a function which checks if indices and ilut
        ! are compatible, and not initiate the excitation right away
        ! but check cases again
        ! in checkCompatibility the number of switches is already
        ! calulated to check if the probabilistic weights fit... maybe but
        ! that out and reuse.. to not waste any effort.
        call checkCompatibility(csf_i, excitInfo, compFlag, posSwitches, negSwitches)

        if (.not. compFlag) then
            allocate(excitations(0, 0), stat=ierr)
            return
        end if

        ! with the more involved excitation identification, i should write
        ! really specific double excitation creators

        ! maybe I should adapt this here also for other type of lattice
        ! models with restricted 2-body interaction..
        if (t_mixed_hubbard) then
            select case (excitInfo%typ)
            case (excit_type%single, &
                  excit_type%raising, &
                  excit_type%lowering, &
                  excit_type%single_overlap_lowering, &
                  excit_type%single_overlap_raising)

                allocate(excitations(0, 0), stat=ierr)
                return
            end select
        end if

        if (t_heisenberg_model) then
            if (excitInfo%typ /= excit_type%fullstart_stop_mixed) then
                allocate(excitations(0,0), stat = ierr)
                return
            end if
        end if

        select case (excitInfo%typ)
        case (excit_type%single)
            ! shouldnt be here.. onyl single excits and full weight gens
            allocate(excitations(0, 0), stat=ierr)
            return

        case (excit_type%raising) ! weight + lowering gen.
            ! can be treated almost like a single excitation
            ! essentially the same, except if d(w) == 3 in the excitaton regime
            call calcDoubleExcitation_withWeight(ilut, csf_i, excitInfo, excitations, &
                                                 nExcits, posSwitches, negSwitches)

            exlevel = 1

        case (excit_type%lowering) ! weight + raising gen
            call calcDoubleExcitation_withWeight(ilut, csf_i, excitInfo, excitations, &
                                                 nExcits, posSwitches, negSwitches)

            exlevel = 1

        case (excit_type%non_overlap) ! non overlap
            call calcNonOverlapDouble(ilut, csf_i, excitInfo, excitations, nExcits, &
                                      posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%single_overlap_lowering) ! single overlap two lowering
            ! how can i efficiently adress that?
            ! can i write that efficiently in one function or do i need more?
            ! probably need more... i already determined
            call calcSingleOverlapLowering(ilut, csf_i, excitInfo, excitations, nExcits, &
                                           posSwitches, negSwitches)

            exlevel = 1

        case (excit_type%single_overlap_raising) ! single overlap raising
            call calcSingleOverlapRaising(ilut, csf_i, excitInfo, excitations, nExcits, &
                                          posSwitches, negSwitches)

            exlevel = 1

        case (excit_type%single_overlap_L_to_R) ! single overlap lowering into raising
            call calcSingleOverlapMixed(ilut, csf_i, excitInfo, excitations, nExcits, &
                                        posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%single_overlap_R_to_L) ! single overlap raising into lowering
            call calcSingleOverlapMixed(ilut, csf_i, excitInfo, excitations, nExcits, &
                                        posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_lowering) ! normal double overlap two lowering
            call calcDoubleLowering(ilut, csf_i, excitInfo, excitations, nExcits, &
                                    posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_raising) ! normal double overlap two raising
            call calcDoubleRaising(ilut, csf_i, excitInfo, excitations, nExcits, &
                                   posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_L_to_R_to_L) ! lowering into raising into lowering
            call calcDoubleRaising(ilut, csf_i, excitInfo, excitations, nExcits, &
                                   posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_R_to_L_to_R) ! raising into lowering into raising
            call calcDoubleLowering(ilut, csf_i, excitInfo, excitations, nExcits, &
                                    posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_L_to_R) ! lowering into raising double
            call calcDoubleL2R(ilut, csf_i, excitInfo, excitations, nExcits, &
                               posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%double_R_to_L) ! raising into lowering double
            call calcDoubleR2L(ilut, csf_i, excitInfo, excitations, nExcits, &
                               posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%fullstop_lowering) ! full stop 2 lowering
            ! can i write a function for both alike generator combinations
            ! i think i can
            call calcFullstopLowering(ilut, csf_i, excitInfo, excitations, nExcits, &
                                      posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%fullstop_raising) ! full stop 2 raising
            call calcFullstopRaising(ilut, csf_i, excitInfo, excitations, nExcits, &
                                     posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%fullstop_L_to_R) ! full stop lowering into raising
            call calcFullStopL2R(ilut, csf_i, excitInfo, excitations, nExcits, &
                                 posSwitches, negSwitches)

            ! in this case there is also the possibility for one single-like
            ! excitation if there is no change in the double overlap region!
            ! todo! how to fix that? is that so important? its only max. 1
            exlevel = 2

        case (excit_type%fullstop_R_to_L) ! full stop raising into lowering
            call calcFullStopR2L(ilut, csf_i, excitInfo, excitations, nExcits, &
                                 posSwitches, negSwitches)

            ! same as for 16
            exlevel = 2

        case (excit_type%fullstart_lowering) ! full start 2 lowering
            call calcFullStartLowering(ilut, csf_i, excitInfo, excitations, nExcits, &
                                       posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%fullstart_raising) ! full start 2 raising
            call calcFulLStartRaising(ilut, csf_i, excitInfo, excitations, nExcits, &
                                      posSwitches, negSwitches)

            exlevel = 2

        case (excit_type%fullStart_L_to_R) ! full start lowering into raising
            call calcFullStartL2R(ilut, csf_i, excitInfo, excitations, nExcits, &
                                  posSwitches, negSwitches)

            ! same as for 16
            exlevel = 2

        case (excit_type%fullstart_R_to_L) ! full start raising into lowering
            call calcFullStartR2L(ilut, csf_i, excitInfo, excitations, nExcits, &
                                  posSwitches, negSwitches)

            ! same as for 16
            exlevel = 2

        case (excit_type%fullstart_stop_alike) ! full start into full stop alike
            call calcFullStartFullStopAlike(ilut, csf_i, excitInfo, excitations)
            nExcits = 1

            exlevel = 2

        case (excit_type%fullstart_stop_mixed) ! full start into full stop mixed
            call calcFullStartFullStopMixed(ilut, csf_i, excitInfo, excitations, nExcits, &
                                            posSwitches, negSwitches)

            ! same as for 16
            exlevel = 2

        end select

        ! for checking purposes encode umat/2 here to compare with sandeeps
        ! results for now..
        do n = 1, nExcits

            call encode_matrix_element(excitations(:, n), 0.0_dp, 2)
            if (t_full_) then
                call update_matrix_element(excitations(:, n), umat / 2.0_dp, 1)
            else
                call encode_matrix_element(excitations(:, n), 1.0_dp, 1)
            end if
            ! and also use the deltaB value for finished excitations to
            ! indicate the level of excitation IC for the remaining NECI code
            call setDeltaB(exlevel, excitations(:, n))

        end do

    end subroutine calcAllExcitations_double

    subroutine calcDoubleR2L(ilut, csf_i, excitInfo, excitations, nExcits, &
                             posSwitches, negSwitches)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcDoubleR2L"

        integer :: iOrb, start1, start2, ende1, ende2
        type(WeightObj_t) :: weights
        real(dp) :: plusWeight, minusWeight, zeroWeight
        integer(n_int), allocatable :: tempExcits(:, :)
        !todo asserts

        ASSERT(.not. isThree(ilut, excitInfo%fullStart))
        ASSERT(.not. isZero(ilut, excitInfo%secondStart))
        ASSERT(.not. isZero(ilut, excitInfo%firstEnd))
        ASSERT(.not. isThree(ilut, excitInfo%fullEnd))

        start1 = excitInfo%fullStart
        start2 = excitInfo%secondStart
        ende1 = excitInfo%firstEnd
        ende2 = excitInfo%fullEnd

        ! todo: create correct weights:
        weights = init_fullDoubleWeight(csf_i, start2, ende1, ende2, negSwitches(start2), &
                                        negSwitches(ende1), posSwitches(start2), posSwitches(ende1), &
                                        csf_i%B_real(start2), csf_i%B_real(ende1))

        excitInfo%currentGen = excitInfo%firstGen

        ! then do single start:
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                               weights, tempExcits, nExcits)

        ! and single update until semi start
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! change weights... maybe need both single and double type weights
        ! maybe semistart is wrong here..
        weights = weights%ptr

        minusWeight = weights%proc%minus(negSwitches(start2), csf_i%B_real(start2), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(start2), csf_i%B_real(start2), weights%dat)
        zeroWeight = weights%proc%zero(negSwitches(start2), posSwitches(start2), &
                                       csf_i%B_real(start2), weights%dat)

        ! then do lowering semi start
        call calcLoweringSemiStart(ilut, csf_i, excitInfo, &
                                   tempExcits, nExcits, plusWeight, minusWeight, zeroWeight)

        ! then do double excitation over double excitation region
        do iOrb = excitInfo%secondStart + 1, excitInfo%firstEnd - 1
            call doubleUpdate(ilut, csf_i, iOrb, excitInfo, weights, tempExcits, nExcits, &
                              negSwitches, posSwitches)
        end do

        ! update weights again:
        weights = weights%ptr
        minusWeight = weights%proc%minus(negSwitches(ende1), csf_i%B_real(ende1), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(ende1), csf_i%B_real(ende1), weights%dat)

        ! then do lowering semi stop
        call calcRaisingSemiStop(ilut, csf_i, excitInfo, tempExcits, nExcits, plusWeight, &
                                 minusWeight)

        ! have to set used generators correctly
        excitInfo%currentGen = excitInfo%lastGen

        ! and then do final single region again
        do iOrb = excitInfo%firstEnd + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! and finally end step
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

        ! that should be it...

    end subroutine calcDoubleR2L

    subroutine calcDoubleL2R(ilut, csf_i, excitInfo, excitations, nExcits, &
                             posSwitches, negSwitches)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcDoubleL2R"

        integer :: iOrb, start1, start2, ende1, ende2
        type(WeightObj_t) :: weights
        real(dp) :: plusWeight, minusWeight, zeroWeight
        integer(n_int), allocatable :: tempExcits(:, :)
        !todo asserts

        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(.not. isThree(ilut, excitInfo%secondStart))
        ASSERT(.not. isThree(ilut, excitInfo%firstEnd))
        ASSERT(.not. isZero(ilut, excitInfo%fullEnd))

        start1 = excitInfo%fullStart
        start2 = excitInfo%secondStart
        ende1 = excitInfo%firstEnd
        ende2 = excitInfo%fullEnd

        !  create correct weights:
        weights = init_fullDoubleWeight(csf_i, start2, ende1, ende2, negSwitches(start2), &
                                        negSwitches(ende1), posSwitches(start2), posSwitches(ende1), &
                                        csf_i%B_real(start2), csf_i%B_real(ende1))

        excitInfo%currentGen = excitInfo%firstGen

        ! then do single start:
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                               weights, tempExcits, nExcits)

        ! and single update until semi start
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! change weights... maybe need both single and double type weights
        ! then do lowering semi start
        weights = weights%ptr

        minusWeight = weights%proc%minus(negSwitches(start2), csf_i%B_real(start2), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(start2), csf_i%B_real(start2), weights%dat)
        zeroWeight = weights%proc%zero(negSwitches(start2), posSwitches(start2), &
                                       csf_i%B_real(start2), weights%dat)

        call calcRaisingSemiStart(ilut, csf_i, excitInfo, &
                                  tempExcits, nExcits, plusWeight, minusWeight, zeroWeight)

        ! then do double excitation over double excitation region
        do iOrb = excitInfo%secondStart + 1, excitInfo%firstEnd - 1
            call doubleUpdate(ilut, csf_i, iOrb, excitInfo, weights, tempExcits, nExcits, &
                              negSwitches, posSwitches)

        end do

        ! update weights again: todo
        weights = weights%ptr
        minusWeight = weights%proc%minus(negSwitches(ende1), csf_i%B_real(ende1), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(ende1), csf_i%B_real(ende1), weights%dat)

        ! then do lowering semi stop
        call calcLoweringSemiStop(ilut, csf_i, excitInfo, tempExcits, nExcits, plusWeight, &
                                  minusWeight)

        ! have to set used generators correctly
        excitInfo%currentGen = excitInfo%lastGen

        ! and then do final single region again
        do iOrb = excitInfo%firstEnd + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! and finally end step
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

        ! that should be it...

    end subroutine calcDoubleL2R

    subroutine calcDoubleRaising(ilut, csf_i, excitInfo, excitations, nExcits, &
                                 posSwitches, negSwitches)
        ! this function can deal with 2 raising and also the mixed L->R->L
        ! case since the called functions are the same

        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcDoubleRaising"

        integer :: iOrb, start2, ende1, ende2
        type(WeightObj_t) :: weights
        real(dp) :: plusWeight, minusWeight, zeroWeight
        integer(n_int), allocatable :: tempExcits(:, :)
        !todo asserts

#ifdef DEBUG_
        if (excitInfo%gen1 == excitInfo%gen2) then
            ! two raising
            ASSERT(.not. isThree(ilut, excitInfo%fullStart))
            ASSERT(.not. isThree(ilut, excitInfo%secondStart))
            ASSERT(.not. isZero(ilut, excitInfo%firstEnd))
            ASSERT(.not. isZero(ilut, excitInfo%fullEnd))
        else
            ASSERT(.not. isZero(ilut, excitInfo%fullStart))
            ASSERT(.not. isThree(ilut, excitInfo%secondStart))
            ASSERT(.not. isZero(ilut, excitInfo%firstEnd))
            ASSERT(.not. isThree(ilut, excitInfo%fullEnd))
        end if
#endif

        start2 = excitInfo%secondStart
        ende1 = excitInfo%firstEnd
        ende2 = excitInfo%fullEnd

        ! create correct weights:
        weights = init_fullDoubleWeight(csf_i, start2, ende1, ende2, negSwitches(start2), &
                                        negSwitches(ende1), posSwitches(start2), posSwitches(ende1), &
                                        csf_i%B_real(start2), csf_i%B_real(ende1))

        excitInfo%currentGen = excitInfo%firstGen
        ! then do single start:
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                               weights, tempExcits, nExcits)

        ! and single update until semi start
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        weights = weights%ptr

        minusWeight = weights%proc%minus(negSwitches(start2), csf_i%B_real(start2), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(start2), csf_i%B_real(start2), weights%dat)
        zeroWeight = weights%proc%zero(negSwitches(start2), posSwitches(start2), &
                                       csf_i%B_real(start2), weights%dat)

        ! change weights... maybe need both single and double type weights
        ! then do lowering semi start
        call calcRaisingSemiStart(ilut, csf_i, excitInfo, &
                                  tempExcits, nExcits, plusWeight, minusWeight, zeroWeight)

        ! then do double excitation over double excitation region
        do iOrb = excitInfo%secondStart + 1, excitInfo%firstEnd - 1
            call doubleUpdate(ilut, csf_i, iOrb, excitInfo, weights, tempExcits, nExcits, &
                              negSwitches, posSwitches)

        end do

        ! update weights again:
        weights = weights%ptr
        minusWeight = weights%proc%minus(negSwitches(ende1), csf_i%B_real(ende1), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(ende1), csf_i%B_real(ende1), weights%dat)

        ! then do lowering semi stop
        call calcRaisingSemiStop(ilut, csf_i, excitInfo, tempExcits, nExcits, plusWeight, &
                                 minusWeight)

        ! have to set used generators correctly (dont actually have to do it
        ! here since they dont change

        ! and then do final single region again
        do iOrb = excitInfo%firstEnd + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! and finally end step
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

        ! that should be it...

    end subroutine calcDoubleRaising

    subroutine calcDoubleLowering(ilut, csf_i, excitInfo, excitations, nExcits, &
                                  posSwitches, negSwitches)
        ! this function can deal with 2 lowering and the mixed R->L-R
        ! case, since the called functions are the same
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcDoubleLowering"

        integer :: iOrb, start2, ende1, ende2
        type(WeightObj_t) :: weights
        real(dp) :: plusWeight, minusWeight, zeroWeight
        integer(n_int), allocatable :: tempExcits(:, :)
        !todo asserts

#ifdef DEBUG_
        if (excitInfo%gen1 == excitInfo%gen2) then
            ASSERT(.not. isZero(ilut, excitInfo%fullStart))
            ASSERT(.not. isZero(ilut, excitInfo%secondStart))
            ASSERT(.not. isThree(ilut, excitInfo%firstEnd))
            ASSERT(.not. isThree(ilut, excitInfo%fullEnd))
        else
            ASSERT(.not. isThree(ilut, excitInfo%fullStart))
            ASSERT(.not. isZero(ilut, excitInfo%secondStart))
            ASSERT(.not. isThree(ilut, excitInfo%firstEnd))
            ASSERT(.not. isZero(ilut, excitInfo%fullEnd))
        end if
#endif

        start2 = excitInfo%secondStart
        ende1 = excitInfo%firstEnd
        ende2 = excitInfo%fullEnd

        ! : create correct weights:
        weights = init_fullDoubleWeight(csf_i, start2, ende1, ende2, negSwitches(start2), &
                                        negSwitches(ende1), posSwitches(start2), posSwitches(ende1), &
                                        csf_i%B_real(start2), csf_i%B_real(ende1))

        excitInfo%currentGen = excitInfo%firstGen
        ! then do single start:
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                               weights, tempExcits, nExcits)

        ! and single update until semi start
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! change weights... maybe need both single and double type weights
        ! then do lowering semi start
        weights = weights%ptr

        minusWeight = weights%proc%minus(negSwitches(start2), csf_i%B_real(start2), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(start2), csf_i%B_real(start2), weights%dat)
        zeroWeight = weights%proc%zero(negSwitches(start2), posSwitches(start2), &
                                       csf_i%B_real(start2), weights%dat)

        call calcLoweringSemiStart(ilut, csf_i, excitInfo, &
                                   tempExcits, nExcits, plusWeight, minusWeight, zeroWeight)

        ! then do double excitation over double excitation region
        do iOrb = excitInfo%secondStart + 1, excitInfo%firstEnd - 1
            call doubleUpdate(ilut, csf_i, iOrb, excitInfo, weights, tempExcits, nExcits, &
                              negSwitches, posSwitches)
        end do

        ! update weights again:
        weights = weights%ptr
        minusWeight = weights%proc%minus(negSwitches(ende1), csf_i%B_real(ende1), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(ende1), csf_i%B_real(ende1), weights%dat)

        ! then do lowering semi stop
        call calcLoweringSemiStop(ilut, csf_i, excitInfo, tempExcits, nExcits, plusWeight, &
                                  minusWeight)

        ! have to set the used generators correctly to handle more versions

        ! and then do final single region again
        do iOrb = excitInfo%firstEnd + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! and finally end step
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

    end subroutine calcDoubleLowering

    subroutine calcFullStartFullStopAlike(ilut, csf_i, excitInfo, excitations)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        character(*), parameter :: this_routine = "calcFullStartFullStopAlike"

        integer :: ierr
        integer(n_int) :: t(0:nifguga)
        real(dp) :: nOpen
        ! full start full stop alike is pretty easy. actually there is only
        ! one possible excitation, and the matrix element sign, just depends
        ! on the number of open orbitals in the excitation range,
        ! so just count that and check if excitation is possible.

        ! assert again just to be save
        ASSERT(isZero(ilut, excitInfo%i))
        ASSERT(isThree(ilut, excitInfo%j))

        ! is only one excitation possible
        allocate(excitations(0:nifguga, 1), stat=ierr)

        ! where everything is the same as in ilut, except at full start and stop
        t = ilut

        ! change start/end depending on type of excitation
        set_orb(t, 2 * excitInfo%i)
        set_orb(t, 2 * excitInfo%i - 1)

        clr_orb(t, 2 * excitInfo%j)
        clr_orb(t, 2 * excitInfo%j - 1)

        ! matrix element deüends only on the number of open orbitals in the
        ! excitaiton region
        nOpen = real(count_open_orbs_ij(csf_i, excitInfo%fullStart, excitInfo%fullEnd), dp)

        ! update! the sum over two-particle integrals involves a 1/2, which
        ! does not get compensated here by

        call encode_matrix_element(t, 0.0_dp, 2)
        call encode_matrix_element(t, 2.0_dp * (-1.0_dp)**nOpen, 1)

        if (tFillingStochRDMOnFly) then
            call encode_stochastic_rdm_info(GugaBits, t, rdm_ind= &
                                            contract_2_rdm_ind(excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l, &
                                                               excit_lvl=2, excit_typ=excitInfo%typ), x1=0.0_dp, &
                                            x0=extract_matrix_element(t, 1))
        end if

        excitations(:, 1) = t

    end subroutine calcFullStartFullStopAlike

    subroutine calcFullStartFullStopMixed(ilut, csf_i, excitInfo, excitations, nExcits, &
                                          posSwitches, negSwitches)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(n_int), allocatable :: tempExcits(:, :)
        real(dp) :: plusWeight, minusWeight, zeroWeight
        type(WeightObj_t) :: weights
        integer :: iOrb

        ! if 3 at full start or full end, this is then a diagonal element
        ! actually, and should thus be already treated by diagonal matrix
        ! calculator. -> should i take that into account for
        ! compatibility check and excitation type determination?
        ! probably yes! -> so assume such an excitation is not coming

        associate(st => excitInfo%fullStart, en => excitInfo%fullEnd)

            if (csf_i%stepvector(st) == 0 .or. csf_i%stepvector(en) == 0) then
                nExcits = 0
                allocate(excitations(0, 0))
                return
            end if

            if (csf_i%stepvector(st) == 3 .or. csf_i%stepvector(en) == 3) then
                nExcits = 1
                allocate(excitations(0:nifguga, nExcits))
                excitations(:, 1) = ilut
                call encode_matrix_element(excitations(:, 1), 0.0_dp, 2)
                call encode_matrix_element(excitations(:, 1), &
                                           -real(csf_i%Occ_int(st) * csf_i%Occ_int(en), dp) / 2.0_dp, 1)

                return
            end if

            ! so only have to deal with 1, or 2 at start and end -> write this
            ! function specifically for these cases.

            ! can i just use already implemented fullStart? i think
            weights = init_doubleWeight(csf_i, en)
            plusWeight = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            minusWeight = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)
            zeroWeight = weights%proc%zero(negSwitches(st), posSwitches(st), &
                                           csf_i%B_real(st), weights%dat)

            ! then call it
            call mixedFullStart(ilut, csf_i, excitInfo, plusWeight, minusWeight, zeroWeight, tempExcits, &
                                nExcits)

            ! and just do double update for the excitation region
            do iOrb = st + 1, en - 1
                call doubleUpdate(ilut, csf_i, iOrb, excitInfo, weights, tempExcits, nExcits, &
                                  negSwitches, posSwitches)
            end do

            ! and then to already implemented mixed end
            call mixedFullStop(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

        end associate

    end subroutine calcFullStartFullStopMixed

    subroutine calcFullStartR2L(ilut, csf_i, excitInfo, excitations, nExcits, &
                                posSwitches, negSwitches, t_no_singles_opt)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        logical, intent(in), optional :: t_no_singles_opt
        character(*), parameter :: this_routine = "calcFullStartR2L"

        integer :: ierr, iOrb, start, ende, semi, gen, start2
        type(WeightObj_t) :: weights
        real(dp) :: minusWeight, plusWeight, zeroWeight
        integer(n_int), allocatable :: tempExcits(:, :)
        logical :: t_no_singles

        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(isProperCSF_ilut(ilut))

        ! create the fullStart
        start = excitInfo%fullStart
        ende = excitInfo%fullEnd
        semi = excitInfo%firstEnd
        gen = excitInfo%firstGen

        if (present(t_no_singles_opt)) then
            t_no_singles = t_no_singles_opt
        else
            t_no_singles = .false.
        end if

        if (t_no_singles .and. csf_i%stepvector(start) == 3) then
            nExcits = 0
            allocate(excitations(0, 0), stat=ierr)
            return
        end if

        ! set up weights
        start2 = excitInfo%secondStart

        if (t_mixed_hubbard) then
            if (csf_i%stepvector(start) == 3) then
                nExcits = 0
                allocate(excitations(0, 0), stat=ierr)
                return
            end if
        end if

        ! create correct weights:
        weights = init_fullStartWeight(csf_i, semi, ende, negSwitches(semi), &
                                       posSwitches(semi), csf_i%B_real(semi))

        minusWeight = weights%proc%minus(negSwitches(start), csf_i%B_real(start), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(start), csf_i%B_real(start), weights%dat)
        zeroWeight = weights%proc%zero(negSwitches(start), posSwitches(start), &
                                       csf_i%B_real(start), weights%dat)

        ! check if first value is 3, so only 0 branch is compatible
        call mixedFullStart(ilut, csf_i, excitInfo, plusWeight, minusWeight, zeroWeight, tempExcits, &
                            nExcits)

        ! then do pseudo double until semi stop
        ! should check for LR(3) start here, have to do nothing if a 3 at
        ! the full start since all matrix elements are one..

        if (csf_i%stepvector(start) /= 3) then
            do iOrb = start + 1, semi - 1
                call doubleUpdate(ilut, csf_i, iOrb, excitInfo, weights, tempExcits, nExcits, &
                                  negSwitches, posSwitches)
            end do
        end if

        ! then deal with the specific semi-stop here
        ! but update weights here..
        ! then reset weights !
        ! do i only need single weight here?
        weights = weights%ptr

        plusWeight = weights%proc%plus(posSwitches(semi), csf_i%B_real(semi), &
                                       weights%dat)
        minusWeight = weights%proc%minus(negSwitches(semi), csf_i%B_real(semi), &
                                         weights%dat)

        call calcRaisingSemiStop(ilut, csf_i, excitInfo, tempExcits, nExcits, plusWeight, &
                                 minusWeight, t_no_singles)

        excitInfo%currentGen = excitInfo%lastGen
        ! and continue on with single excitation region
        do iOrb = semi + 1, ende - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! and normal single end
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

    end subroutine calcFullStartR2L

    subroutine calcFullStartL2R(ilut, csf_i, excitInfo, excitations, nExcits, &
                                posSwitches, negSwitches, t_no_singles_opt)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        logical, intent(in), optional :: t_no_singles_opt
        character(*), parameter :: this_routine = "calcFullStartL2R"

        integer :: ierr, iOrb, start, ende, semi, gen
        real(dp) :: minusWeight, plusWeight, zeroWeight
        type(WeightObj_t) :: weights
        integer(n_int), allocatable :: tempExcits(:, :)
        logical :: t_no_singles

        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(isProperCSF_ilut(ilut))

        ! create the fullStart
        start = excitInfo%fullStart
        ende = excitInfo%fullEnd
        semi = excitInfo%firstEnd
        gen = excitInfo%firstGen

        if (present(t_no_singles_opt)) then
            t_no_singles = t_no_singles_opt
        else
            t_no_singles = .false.
        end if

        if (t_no_singles .and. csf_i%stepvector(start) == 3) then
            nExcits = 0
            allocate(excitations(0, 0), stat=ierr)
            return
        end if

        ! create correct weights:
        weights = init_fullStartWeight(csf_i, semi, ende, negSwitches(semi), &
                                       posSwitches(semi), csf_i%B_real(semi))

        minusWeight = weights%proc%minus(negSwitches(start), csf_i%B_real(start), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(start), csf_i%B_real(start), weights%dat)
        zeroWeight = weights%proc%zero(negSwitches(start), posSwitches(start), &
                                       csf_i%B_real(start), weights%dat)

        if (t_mixed_hubbard) then
            if (csf_i%stepvector(start) == 3) then
                nExcits = 0
                allocate(excitations(0, 0), stat=ierr)
                return
            end if
        end if

        ! check if first value is 3, so only 0 branch is compatible
        call mixedFullStart(ilut, csf_i, excitInfo, plusWeight, minusWeight, zeroWeight, tempExcits, &
                            nExcits)

        ! then do pseudo double until semi stop
        ! should check for LR(3) start here, have to do nothing if a 3 at
        ! the full start since all matrix elements are one..
        if (csf_i%stepvector(start) /= 3) then
            do iOrb = start + 1, semi - 1
                call doubleUpdate(ilut, csf_i, iOrb, excitInfo, weights, tempExcits, nExcits, &
                                  negSwitches, posSwitches)
            end do
        end if

        if (t_mixed_hubbard) then
            ! TODO abort x0 branch..
        end if

        ! then deal with the specific semi-stop here
        ! but todo update weights here..
        weights = weights%ptr
        plusWeight = weights%proc%plus(posSwitches(semi), csf_i%B_real(semi), &
                                       weights%dat)
        minusWeight = weights%proc%minus(negSwitches(semi), csf_i%B_real(semi), &
                                         weights%dat)

        call calcLoweringSemiStop(ilut, csf_i, excitInfo, tempExcits, nExcits, plusWeight, &
                                  minusWeight, t_no_singles)

        ! then reset weights todo!
        excitInfo%currentGen = excitInfo%lastGen

        ! and continue on with single excitation region
        do iOrb = semi + 1, ende - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! and normal single end
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

    end subroutine calcFullStartL2R

    subroutine calcRaisingSemiStop(ilut, csf_i, excitInfo, tempExcits, nExcits, plusWeight, &
                                   minusWeight, t_no_singles_opt)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), intent(inout), allocatable :: tempExcits(:, :)
        integer, intent(inout) :: nExcits
        real(dp), intent(in) :: plusWeight, minusWeight
        logical, intent(in), optional :: t_no_singles_opt
        character(*), parameter :: this_routine = "calcRaisingSemiStop"

        integer :: se, iEx, deltaB, st, ss
        integer(n_int) :: t(0:nifguga), s(0:nifguga)
        real(dp) :: tempWeight, tempWeight_0, tempWeight_1, bVal
        logical :: t_no_singles

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(.not. isZero(ilut, excitInfo%firstEnd))
        ASSERT(plusWeight >= 0.0_dp)
        ASSERT(minusWeight >= 0.0_dp)

        se = excitInfo%firstEnd
        bVal = csf_i%B_real(se)
        st = excitInfo%fullStart
        ss = excitInfo%secondStart

        if (present(t_no_singles_opt)) then
            t_no_singles = t_no_singles_opt
        else
            t_no_singles = .false.
        end if

        ! do it similar to semi-starts and check if its a pseudo excit
        ! TODO!! have to make additional if here, to check i comes from
        ! a full start, because i also want to use it for normal double
        ! excitations, and there it doesnt matter what the stepvector value
        ! at the fullstart is!!!
        ! but for a double raising fullstart eg. there cant be a 3 at the
        ! fullstart... -> how stupid

        if (csf_i%stepvector(st) == 3 .and. st == ss) then
            ! only 0 branches in this case
            ! first do the non-branching possibs
            select case (csf_i%stepvector(se))
            case (1)
                ! 1 -> 0 switch
                do iEx = 1, nExcits

                    t = tempExcits(:, iEx)

                    ASSERT(getDeltaB(t) == 0)
                    ! also need to assert weight is non-zero as it should be

                    ASSERT(plusWeight > 0.0_dp)

                    clr_orb(t, 2 * se - 1)

                    call getDoubleMatrixElement(0, 1, 0, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                    call update_matrix_element(t, tempWeight, 1)
                    call encode_matrix_element(t, 0.0_dp, 2)

                    call setDeltaB(1, t)

                    tempExcits(:, iEx) = t

                end do

            case (2)
                ! 2 -> 0 switch

                if (near_zero(minusWeight)) then
                    nExcits = 0
                    tempExcits = 0
                    return
                end if
                do iEx = 1, nExcits

                    t = tempExcits(:, iEx)

                    ASSERT(getDeltaB(t) == 0)
                    ASSERT(minusWeight > 0.0_dp)

                    clr_orb(t, 2 * se)

                    call getDoubleMatrixElement(0, 2, 0, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                    call update_matrix_element(t, tempWeight, 1)
                    call encode_matrix_element(t, 0.0_dp, 2)

                    call setDeltaB(-1, t)

                    tempExcits(:, iEx) = t

                end do

            case (3)
                ! three case -> branching possibilities according to b an weights
                if (csf_i%B_int(se) > 0 .and. plusWeight > 0.0_dp &
                    .and. minusWeight > 0.0_dp) then
                    ! both excitations are possible then
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        ASSERT(getDeltaB(t) == 0)

                        ! have to store t weight since i need it twice
                        tempWeight_1 = extract_matrix_element(t, 1)

                        ! do 3->1 first
                        clr_orb(t, 2 * se)

                        call setDeltaB(-1, t)

                        call getDoubleMatrixElement(1, 3, 0, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                        call update_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                        ! then change to 3->2 branch
                        set_orb(t, 2 * se)
                        clr_orb(t, 2 * se - 1)

                        call setDeltaB(+1, t)

                        call getDoubleMatrixElement(2, 3, 0, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                        call encode_matrix_element(t, tempWeight * tempWeight_1, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        nExcits = nExcits + 1

                        tempExcits(:, nExcits) = t

                    end do

                else if (csf_i%B_int(se) == 0 .or. near_zero(plusWeight)) then
                    ! only -1 branch possibloe
                    if (near_zero(minusWeight)) then
                        nExcits = 0
                        tempExcits = 0
                        return
                    end if

                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        ASSERT(getDeltaB(t) == 0)
                        ASSERT(minusWeight > 0.0_dp)

                        ! 3 -> 1
                        clr_orb(t, 2 * se)

                        call setDeltaB(-1, t)

                        call getDoubleMatrixElement(1, 3, 0, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                        call update_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                    end do

                else if (near_zero(minusWeight) .and. csf_i%B_int(se) > 0) then
                    ! only +1 branches possible
                    if (near_zero(plusWeight)) then
                        nExcits = 0
                        tempExcits = 0
                        return
                    end if

                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        ASSERT(getDeltaB(t) == 0)
                        ASSERT(plusWeight > 0.0_dp)

                        ! 3 -> 2
                        clr_orb(t, 2 * se - 1)

                        call setDeltaB(+1, t)

                        call getDoubleMatrixElement(2, 3, 0, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                        call update_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                    end do

                else if (near_zero(minusWeight) .and. csf_i%B_int(se) == 0) then
                    ! in this case no excitaiton is possible due to b value todo
                    call stop_all(this_routine, "implement cancelled excitations")

                else
                    ! shouldnt be here...
                    call stop_all(this_routine, "somethin went wrong! shouldnt be here!")
                end if
            end select

        else
            ! also +2, and -2 branches arriving possible!
            ! again the non-branching values first
            select case (csf_i%stepvector(se))
            case (1)
                ! 1 -> 0 for 0 and -2 branch
                ! have to check different weight combinations
                ! although excitations should only get there if the weight
                ! fits.... -> check that

                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)

                    ! check if weights fit.
                    ! do 1 -> 0
                    clr_orb(t, 2 * se - 1)

                    call setDeltaB(deltaB + 1, t)

                    call getDoubleMatrixElement(0, 1, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, &
                                                excitInfo%order1, tempWeight_0, tempWeight_1)

                    ! after semi-stop i only need the sum of the matrix
                    ! elements

                    if (t_no_singles) then
                        if (.not. near_zero(extract_matrix_element(t, 1))) then
                            ASSERT(DeltaB == 0)
                            call encode_matrix_element(t, 0.0_dp, 1)
                            call encode_matrix_element(t, 0.0_dp, 2)

                        end if
                    end if

                    tempWeight = extract_matrix_element(t, 1) * tempWeight_0 + &
                                 extract_matrix_element(t, 2) * tempWeight_1

                    call encode_matrix_element(t, tempWeight, 1)
                    call encode_matrix_element(t, 0.0_dp, 2)

                    tempExcits(:, iEx) = t

                end do

            case (2)
                ! 2 -> 0 for 0 and +2 branches
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)

                    ! do 2 -> 0
                    clr_orb(t, 2 * se)

                    call setDeltaB(deltaB - 1, t)

                    call getDoubleMatrixElement(0, 2, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, &
                                                excitInfo%order1, tempWeight_0, tempWeight_1)

                    if (t_no_singles) then
                        if (.not. near_zero(extract_matrix_element(t, 1))) then
                            ASSERT(DeltaB == 0)
                            call encode_matrix_element(t, 0.0_dp, 1)
                            call encode_matrix_element(t, 0.0_dp, 2)

                        end if
                    end if

                    tempWeight = extract_matrix_element(t, 1) * tempWeight_0 + &
                                 extract_matrix_element(t, 2) * tempWeight_1

                    call encode_matrix_element(t, tempWeight, 1)
                    call encode_matrix_element(t, 0.0_dp, 2)

                    tempExcits(:, iEx) = t

                end do

            case (3)
                ! 3 case -> have to check weights more thourougly
                ! when a certain +-2 excitation comes to a 0 semi-stop
                ! the ongoing weights should allow an excitation or
                ! otherwise the excitation shouldnt even have been created!
                ! for the 0 branch arriving i have to check if a branching
                ! is possible.. and have to do that outside of the do-loops
                ! to be more efficient
                if (csf_i%B_int(se) > 0 .and. plusWeight > 0.0_dp &
                    .and. minusWeight > 0.0_dp) then
                    ! all excitations for 0 branch possible
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (t_no_singles) then
                            if (.not. near_zero(extract_matrix_element(t, 1))) then
                                ASSERT(DeltaB == 0)
                                call encode_matrix_element(t, 0.0_dp, 1)
                                call encode_matrix_element(t, 0.0_dp, 2)

                                call encode_matrix_element(tempExcits(:, iex), 0.0_dp, 1)
                                call encode_matrix_element(tempExcits(:, iex), 0.0_dp, 2)
                            end if
                        end if

                        if (deltaB == 0) then
                            ! do 3 -> 1 first
                            clr_orb(t, 2 * se)

                            call setDeltaB(-1, t)

                            call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        tempWeight_0, tempWeight_1)

                            tempWeight = extract_matrix_element(t, 1) * tempWeight_0 + &
                                         extract_matrix_element(t, 2) * tempWeight_1

                            call encode_matrix_element(t, tempWeight, 1)
                            call encode_matrix_element(t, 0.0_dp, 2)

                            ! need fresh temp. ilut for matrix elements
                            s = tempExcits(:, iEx)

                            tempExcits(:, iEx) = t

                            ! then do 3->2 also
                            clr_orb(s, 2 * se - 1)

                            call setDeltaB(+1, s)

                            call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        tempWeight_0, tempWeight_1)

                            tempWeight = extract_matrix_element(s, 1) * tempWeight_0 + &
                                         extract_matrix_element(s, 2) * tempWeight_1

                            call encode_matrix_element(s, tempWeight, 1)
                            call encode_matrix_element(s, 0.0_dp, 2)

                            nExcits = nExcits + 1

                            tempExcits(:, nExcits) = s

                        else if (deltaB == -2) then
                            ! only 3 -> 2 branch
                            clr_orb(t, 2 * se - 1)

                            call setDeltaB(-1, t)

                            ! not sure anymore how matrix elements get
                            ! adressed. maybe with incoming/outgoing deltaB
                            ! todo! check that!!!
                            ! x0 matrix element is always 0 in this case
                            ! only take out x1 element
                            call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        x1_element=tempWeight_1)

                            tempWeight_1 = extract_matrix_element(t, 2) * tempWeight_1
                            call encode_matrix_element(t, tempWeight_1, 1)
                            call encode_matrix_element(t, 0.0_dp, 2)

                            tempExcits(:, iEx) = t

                        else
                            ! only 3 -> 1 branch
                            clr_orb(t, 2 * se)

                            call setDeltaB(deltaB - 1, t)

                            call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        x1_element=tempWeight_1)

                            tempWeight_1 = extract_matrix_element(t, 2) * tempWeight_1
                            call encode_matrix_element(t, tempWeight_1, 1)
                            call encode_matrix_element(t, 0.0_dp, 2)

                            tempExcits(:, iEx) = t

                        end if
                    end do

                else if (csf_i%B_int(se) == 0 .or. near_zero(plusWeight)) then
                    ! only -1 branch when 0 branch arrives... the switch from
                    ! +2 -> +1 branch shouldnt be affected, since i wouldn not
                    ! arrive at semi.stop if 0 weight, and if b value would
                    ! be so low it also wouldn be possible to have this kind
                    ! of excitation at this point
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (t_no_singles) then
                            if (.not. near_zero(extract_matrix_element(t, 1))) then
                                ASSERT(DeltaB == 0)
                                call encode_matrix_element(t, 0.0_dp, 1)
                                call encode_matrix_element(t, 0.0_dp, 2)

                            end if
                        end if

                        ASSERT(deltaB /= 2)

                        if (deltaB == 0) then
                            ! do 3 -> 1 switch
                            clr_orb(t, 2 * se)

                            call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        tempWeight_0, tempWeight_1)

                            tempWeight = extract_matrix_element(t, 1) * tempWeight_0 + &
                                         extract_matrix_element(t, 2) * tempWeight_1

                        else
                            ! -2 branch

                            ! do 3->2
                            clr_orb(t, 2 * se - 1)

                            call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        x1_element=tempWeight_1)

                            tempWeight = extract_matrix_element(t, 2) * tempWeight_1

                        end if

                        call setDeltaB(-1, t)

                        call encode_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                    end do

                else if (csf_i%B_int(se) > 0 .and. near_zero(minusWeight)) then
                    ! only +1 branch possible afterwards
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (t_no_singles) then
                            if (.not. near_zero(extract_matrix_element(t, 1))) then
                                ASSERT(DeltaB == 0)
                                call encode_matrix_element(t, 0.0_dp, 1)
                                call encode_matrix_element(t, 0.0_dp, 2)

                            end if
                        end if

                        ASSERT(deltaB /= -2)

                        if (deltaB == 0) then
                            ! do 3->2 +1 branch
                            clr_orb(t, 2 * se - 1)

                            call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        tempWeight_0, tempWeight_1)

                            tempWeight = extract_matrix_element(t, 1) * tempWeight_0 + &
                                         extract_matrix_element(t, 2) * tempWeight_1

                        else
                            ! +2 branch
                            ! do 3 -> 1
                            clr_orb(t, 2 * se)

                            call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        x1_element=tempWeight_1)

                            tempWeight = extract_matrix_element(t, 2) * tempWeight_1

                        end if

                        call setDeltaB(+1, t)

                        call encode_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                    end do

                else if (csf_i%B_int(se) == 0 .and. near_zero(plusWeight)) then
                    ! broken excitation due to b value restriction
                    ! todo how to deal with that ...
                    call stop_all(this_routine, "broken excitation due to b value. todo!")

                else
                    ! shouldnt be here
                    call stop_all(this_routine, "something went wrong. shouldnt be here!")
                end if

            end select
        end if

    end subroutine calcRaisingSemiStop

    subroutine calcLoweringSemiStop(ilut, csf_i, excitInfo, tempExcits, nExcits, plusWeight, &
                                    minusWeight, t_no_singles_opt)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), intent(inout), allocatable :: tempExcits(:, :)
        integer, intent(inout) :: nExcits
        real(dp), intent(in) :: plusWeight, minusWeight
        logical, intent(in), optional :: t_no_singles_opt
        character(*), parameter :: this_routine = "calcLoweringSemiStop"

        integer :: se, iEx, deltaB, st, ss
        integer(n_int) :: t(0:nifguga), s(0:nifguga)
        real(dp) :: tempWeight, tempWeight_0, tempWeight_1, bVal
        logical :: t_no_singles

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(.not. isThree(ilut, excitInfo%firstEnd))
        ASSERT(plusWeight >= 0.0_dp)
        ASSERT(minusWeight >= 0.0_dp)

        se = excitInfo%firstEnd
        bVal = csf_i%B_real(se)
        st = excitInfo%fullStart
        ss = excitInfo%secondStart

        if (present(t_no_singles_opt)) then
            t_no_singles = t_no_singles_opt
        else
            t_no_singles = .false.
        end if

        ! do it similar to semi-starts and check if its a pseudo excit
        ! TODO!! have to make additional if here, to check i comes from
        ! a full start, because i also want to use it for normal double
        ! excitations, and there it doesnt matter what the stepvector value
        ! at the fullstart is!!!
        if (csf_i%stepvector(st) == 3 .and. ss == st) then
            ! only 0 branches in this case
            ! first do the non-branching possibs
            select case (csf_i%stepvector(se))
            case (1)
                ! 1 -> 3 switch
                if (near_zero(plusWeight)) then
                    tempExcits = 0
                    nExcits = 0
                    return
                end if

                do iEx = 1, nExcits

                    t = tempExcits(:, iEx)

                    ASSERT(getDeltaB(t) == 0)
                    ! also need to assert weight is non-zero as it should be
                    ! for mixed full-starts or general double excitations
                    ! i cant really guaruantee it comes here with
                    ! a plusWeight > 0, since 0 branch is always a
                    ! possibility -> so i have to remove excitations if the
                    ! weight is 0

                    set_orb(t, 2 * se)

                    call getDoubleMatrixElement(3, 1, 0, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                    call update_matrix_element(t, tempWeight, 1)
                    call encode_matrix_element(t, 0.0_dp, 2)

                    call setDeltaB(1, t)

                    tempExcits(:, iEx) = t

                end do

            case (2)
                ! 2 -> 3 switch
                if (near_zero(minusWeight)) then
                    tempExcits = 0
                    nExcits = 0
                    return
                end if

                do iEx = 1, nExcits

                    t = tempExcits(:, iEx)

                    ASSERT(getDeltaB(t) == 0)

                    set_orb(t, 2 * se - 1)

                    call getDoubleMatrixElement(3, 2, 0, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                    call update_matrix_element(t, tempWeight, 1)
                    call encode_matrix_element(t, 0.0_dp, 2)

                    call setDeltaB(-1, t)

                    tempExcits(:, iEx) = t

                end do

            case (0)
                ! zero case -> branching possibilities according to b an weights
                if (csf_i%B_int(se) > 0 .and. plusWeight > 0.0_dp &
                    .and. minusWeight > 0.0_dp) then
                    ! both excitations are possible then
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        ASSERT(getDeltaB(t) == 0)

                        ! have to store t weight since i need it twice
                        tempWeight_1 = extract_matrix_element(t, 1)

                        ! do 0->1 first
                        set_orb(t, 2 * se - 1)

                        call setDeltaB(-1, t)

                        call getDoubleMatrixElement(1, 0, 0, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                        call update_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                        ! then change to 0->2 branch
                        set_orb(t, 2 * se)
                        clr_orb(t, 2 * se - 1)

                        call setDeltaB(+1, t)

                        call getDoubleMatrixElement(2, 0, 0, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                        call encode_matrix_element(t, tempWeight * tempWeight_1, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        nExcits = nExcits + 1

                        tempExcits(:, nExcits) = t

                    end do

                else if (csf_i%B_int(se) == 0 .or. near_zero(plusWeight)) then
                    ! only -1 branch possible
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        ASSERT(getDeltaB(t) == 0)
                        ASSERT(minusWeight > 0.0_dp)

                        ! 0 -> 1
                        set_orb(t, 2 * se - 1)

                        call setDeltaB(-1, t)

                        call getDoubleMatrixElement(1, 0, 0, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                        call update_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                    end do

                else if (near_zero(minusWeight) .and. csf_i%B_int(se) > 0) then
                    ! only +1 branches possible
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        ASSERT(getDeltaB(t) == 0)
                        ASSERT(plusWeight > 0.0_dp)

                        ! 0 -> 2
                        set_orb(t, 2 * se)

                        call setDeltaB(+1, t)

                        call getDoubleMatrixElement(2, 0, 0, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, 1.0_dp, tempWeight)

                        call update_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                    end do

                else if (near_zero(minusWeight) .and. csf_i%B_int(se) == 0) then
                    ! in this case no excitaiton is possible due to b value todo
                    call stop_all(this_routine, "implement cancelled excitations")

                else
                    ! shouldnt be here...
                    call stop_all(this_routine, "somethin went wrong! shouldnt be here!")
                end if
            end select

        else
            ! also +2, and -2 branches arriving possible!
            ! again the non-branching values first
            select case (csf_i%stepvector(se))
            case (1)
                ! 1 -> 3 for 0 and -2 branch
                ! have to check different weight combinations
                ! although excitations should only get there if the weight
                ! fits.... -> check that

                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)

                    if (t_no_singles) then
                        if (.not. near_zero(extract_matrix_element(t, 1))) then
                            ASSERT(DeltaB == 0)
                            call encode_matrix_element(t, 0.0_dp, 1)
                            call encode_matrix_element(t, 0.0_dp, 2)
                        end if
                    end if

                    ! check if weights fit.
                    ! do not need actually, since no choice in excitations...
                    ! and if dealt correctly until now it should be fine.
                    ! hopefully atleast

                    ! do 1 -> 3
                    set_orb(t, 2 * se)

                    call setDeltaB(deltaB + 1, t)

                    call getDoubleMatrixElement(3, 1, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, &
                                                excitInfo%order1, tempWeight_0, tempWeight_1)

                    ! after semi-stop i only need the sum of the matrix
                    ! elements

                    tempWeight = extract_matrix_element(t, 1) * tempWeight_0 + &
                                 extract_matrix_element(t, 2) * tempWeight_1

                    call encode_matrix_element(t, tempWeight, 1)
                    call encode_matrix_element(t, 0.0_dp, 2)

                    tempExcits(:, iEx) = t

                end do

            case (2)
                ! 2 -> 3 for 0 and +2 branches
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)
                    deltaB = getDeltaB(t)

                    if (t_no_singles) then
                        if (.not. near_zero(extract_matrix_element(t, 1))) then
                            ASSERT(DeltaB == 0)
                            call encode_matrix_element(t, 0.0_dp, 1)
                            call encode_matrix_element(t, 0.0_dp, 2)
                        end if
                    end if

                    ! do 2 -> 3
                    set_orb(t, 2 * se - 1)

                    call setDeltaB(deltaB - 1, t)

                    call getDoubleMatrixElement(3, 2, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, &
                                                excitInfo%order1, tempWeight_0, tempWeight_1)

                    tempWeight = extract_matrix_element(t, 1) * tempWeight_0 + &
                                 extract_matrix_element(t, 2) * tempWeight_1

                    call encode_matrix_element(t, tempWeight, 1)
                    call encode_matrix_element(t, 0.0_dp, 2)

                    tempExcits(:, iEx) = t

                end do

            case (0)
                ! 0 case -> have to check weights more thourougly
                ! when a certain +-2 excitation comes to a 0 semi-stop
                ! the ongoing weights should allow an excitation or
                ! otherwise the excitation shouldnt even have been created!
                ! for the 0 branch arriving i have to check if a branching
                ! is possible.. and have to do that outside of the do-loops
                ! to be more efficient
                if (csf_i%B_int(se) > 0 .and. plusWeight > 0.0_dp &
                    .and. minusWeight > 0.0_dp) then
                    ! all excitations for 0 branch possible
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (t_no_singles) then
                            if (.not. near_zero(extract_matrix_element(t, 1))) then
                                ASSERT(DeltaB == 0)
                                call encode_matrix_element(t, 0.0_dp, 1)
                                call encode_matrix_element(t, 0.0_dp, 2)

                                call encode_matrix_element(tempExcits(:, iex), 0.0_dp, 1)
                                call encode_matrix_element(tempExcits(:, iex), 0.0_dp, 2)

                            end if
                        end if

                        if (deltaB == 0) then
                            ! do 0 -> 1 first
                            set_orb(t, 2 * se - 1)

                            call setDeltaB(-1, t)

                            call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        tempWeight_0, tempWeight_1)

                            tempWeight = extract_matrix_element(t, 1) * tempWeight_0 + &
                                         extract_matrix_element(t, 2) * tempWeight_1

                            call encode_matrix_element(t, tempWeight, 1)
                            call encode_matrix_element(t, 0.0_dp, 2)

                            ! need fresh temp. ilut for matrix elements
                            s = tempExcits(:, iEx)

                            tempExcits(:, iEx) = t

                            ! then do 0->2 also
                            set_orb(s, 2 * se)

                            call setDeltaB(+1, s)

                            call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        tempWeight_0, tempWeight_1)

                            tempWeight = extract_matrix_element(s, 1) * tempWeight_0 + &
                                         extract_matrix_element(s, 2) * tempWeight_1

                            call encode_matrix_element(s, tempWeight, 1)
                            call encode_matrix_element(s, 0.0_dp, 2)

                            nExcits = nExcits + 1

                            tempExcits(:, nExcits) = s

                        else if (deltaB == -2) then
                            ! only 0 -> 2 branch
                            set_orb(t, 2 * se)

                            call setDeltaB(-1, t)

                            ! not sure anymore how matrix elements get
                            ! adressed. maybe with incoming/outgoing deltaB
                            ! todo! check that!!!
                            ! x0 matrix element is always 0 in this case
                            ! only take out x1 element
                            call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        x1_element=tempWeight_1)

                            tempWeight_1 = extract_matrix_element(t, 2) * tempWeight_1
                            call encode_matrix_element(t, tempWeight_1, 1)
                            call encode_matrix_element(t, 0.0_dp, 2)

                            tempExcits(:, iEx) = t

                        else
                            ! only 0 -> 1 branch
                            set_orb(t, 2 * se - 1)

                            call setDeltaB(deltaB - 1, t)

                            call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        x1_element=tempWeight_1)

                            tempWeight_1 = extract_matrix_element(t, 2) * tempWeight_1
                            call encode_matrix_element(t, tempWeight_1, 1)
                            call encode_matrix_element(t, 0.0_dp, 2)

                            tempExcits(:, iEx) = t

                        end if
                    end do

                else if (csf_i%B_int(se) == 0 .or. near_zero(plusWeight)) then
                    ! only -1 branch when 0 branch arrives... the switch from
                    ! +2 -> +1 branch shouldnt be affected, since i wouldn not
                    ! arrive at semi.stop if 0 weight, and if b value would
                    ! be so low it also wouldn be possible to have this kind
                    ! of excitation at this point
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (t_no_singles) then
                            if (.not. near_zero(extract_matrix_element(t, 1))) then
                                ASSERT(DeltaB == 0)
                                call encode_matrix_element(t, 0.0_dp, 1)
                                call encode_matrix_element(t, 0.0_dp, 2)

                            end if
                        end if

                        ASSERT(deltaB /= 2)
                        if (deltaB == 0) then
                            ! do 0 -> 1 switch
                            set_orb(t, 2 * se - 1)

                            call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        tempWeight_0, tempWeight_1)

                            tempWeight = extract_matrix_element(t, 1) * tempWeight_0 + &
                                         extract_matrix_element(t, 2) * tempWeight_1

                        else
                            ! -2 branch

                            ! do 0->2
                            set_orb(t, 2 * se)

                            call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        x1_element=tempWeight_1)

                            tempWeight = extract_matrix_element(t, 2) * tempWeight_1

                        end if

                        call setDeltaB(-1, t)

                        call encode_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                    end do

                else if (csf_i%B_int(se) > 0 .and. near_zero(minusWeight)) then
                    ! only +1 branch possible afterwards
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (t_no_singles) then
                            if (.not. near_zero(extract_matrix_element(t, 1))) then
                                ASSERT(DeltaB == 0)
                                call encode_matrix_element(t, 0.0_dp, 1)
                                call encode_matrix_element(t, 0.0_dp, 2)

                            end if
                        end if

                        ASSERT(deltaB /= -2)

                        if (deltaB == 0) then
                            ! do 0->2 +1 branch
                            set_orb(t, 2 * se)

                            call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        tempWeight_0, tempWeight_1)

                            tempWeight = extract_matrix_element(t, 1) * tempWeight_0 + &
                                         extract_matrix_element(t, 2) * tempWeight_1

                        else
                            ! +2 branch
                            ! do 0 -> 1
                            set_orb(t, 2 * se - 1)

                            call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order1, &
                                                        x1_element=tempWeight_1)

                            tempWeight = extract_matrix_element(t, 2) * tempWeight_1

                        end if

                        call setDeltaB(+1, t)

                        call encode_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                    end do

                else if (csf_i%B_int(se) == 0 .and. near_zero(plusWeight)) then
                    ! broken excitation due to b value restriction
                    ! todo how to deal with that ...
                    call stop_all(this_routine, "broken excitation due to b value. todo!")

                else
                    ! shouldnt be here
                    call stop_all(this_routine, "something went wrong. shouldnt be here!")
                end if

            end select
        end if

    end subroutine calcLoweringSemiStop

    subroutine mixedFullStart(ilut, csf_i, excitInfo, plusWeight, minusWeight, &
                              zeroWeight, tempExcits, nExcits)
        ! remember full-start matrix element are stored in the same row
        ! as deltaB = -1 mixed ones... so access matrix element below with
        ! deltaB = -1 !!
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp), intent(in) :: plusWeight, minusWeight, zeroWeight
        integer(n_int), intent(out), allocatable :: tempExcits(:, :)
        integer, intent(out) :: nExcits
        character(*), parameter :: this_routine = "mixedFullStart"

        integer :: st, nmax, ierr
        integer(n_int) :: t(0:nifguga)
        real(dp) :: tempWeight, tempWeight_1, bVal

        ASSERT(.not. isZero(ilut, excitInfo%fullStart))
        ASSERT(isProperCSF_ilut(ilut))

        ! depending on the stepvector create certain starting branches

        st = excitInfo%fullStart
        bVal = csf_i%B_real(st)

        ! determine worst case amount of excitations:
        nMax = 2 + 2**count_open_orbs_ij(csf_i, st, excitInfo%fullEnd)
        allocate(tempExcits(0:nifguga, nMax), stat=ierr)

        ! assert that at least one of the weights is non-zero
        ! otherwise checkCompatibility() has done smth wrong

        t = ilut

        select case (csf_i%stepvector(st))
        case (3)
            ! only deltaB 0 branch possible, and even no change in stepvector
            call encode_matrix_element(t, -Root2, 1)
            call encode_matrix_element(t, 0.0_dp, 2)

            call setDeltaB(0, t)
            tempExcits(:, 1) = t
            nExcits = 1

            ! do the mixed fullstart new with all weights contributed for
        case (1)
            ASSERT(zeroWeight + plusWeight > 0.0_dp)
            if (zeroWeight > 0.0_dp .and. plusWeight > 0.0_dp) then

                ! depending on weights, and b value maybe 2 excitations possible
                ! zero weight always > 0, so only check other
                ! both excitations possible
                ! first do 1->1
                call setDeltaB(0, t)

                call getDoubleMatrixElement(1, 1, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, tempWeight, tempWeight_1)

                call encode_matrix_element(t, tempWeight, 1)
                call encode_matrix_element(t, tempWeight_1, 2)

                tempExcits(:, 1) = t

                ! the  change it to the 1->2 start

                clr_orb(t, 2 * st - 1)
                set_orb(t, 2 * st)

                call setDeltaB(2, t)

                call getDoubleMatrixElement(2, 1, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, x1_element=tempWeight_1)

                call encode_matrix_element(t, 0.0_dp, 1)
                call encode_matrix_element(t, tempWeight_1, 2)

                tempExcits(:, 2) = t

                nExcits = 2

            else if (near_zero(plusWeight)) then
                ! only 0 branch possible
                call setDeltaB(0, t)

                call getDoubleMatrixElement(1, 1, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, tempWeight, tempWeight_1)

                call encode_matrix_element(t, tempWeight, 1)
                call encode_matrix_element(t, tempWeight_1, 2)

                tempExcits(:, 1) = t

                nExcits = 1

            else if (near_zero(zeroWeight)) then
                ! only the switch to the +2 branch valid

                clr_orb(t, 2 * st - 1)
                set_orb(t, 2 * st)

                call setDeltaB(2, t)

                call getDoubleMatrixElement(2, 1, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, x1_element=tempWeight_1)

                call encode_matrix_element(t, 0.0_dp, 1)
                call encode_matrix_element(t, tempWeight_1, 2)

                tempExcits(:, 1) = t

                nExcits = 1

            else
                ! something went wrong probably... 0 branch always
                ! possible remember, so only 0 branch not possible..
                call stop_all(this_routine, "something went wrong. should not be here!")
            end if

        case (2)
            ! here b value is never a problem so i just have to check for
            ! minusWeight
            ASSERT(minusWeight + zeroWeight > 0.0_dp)

            if (near_zero(minusWeight)) then
                ! only 0 branch possible
                call setDeltaB(0, t)

                call getDoubleMatrixElement(2, 2, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, tempWeight, tempWeight_1)

                call encode_matrix_element(t, tempWeight, 1)
                call encode_matrix_element(t, tempWeight_1, 2)

                tempExcits(:, 1) = t

                nExcits = 1

            else if (minusWeight > 0.0_dp .and. zeroWeight > 0.0_dp) then
                ! both branches possible, ! do 0 first
                call setDeltaB(0, t)

                call getDoubleMatrixElement(2, 2, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, tempWeight, tempWeight_1)

                call encode_matrix_element(t, tempWeight, 1)
                call encode_matrix_element(t, tempWeight_1, 2)

                tempExcits(:, 1) = t

                ! then change that to 2->1 start
                set_orb(t, 2 * st - 1)
                clr_orb(t, 2 * st)

                call setDeltaB(-2, t)

                call getDoubleMatrixElement(1, 2, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, x1_element=tempWeight_1)

                call encode_matrix_element(t, 0.0_dp, 1)
                call encode_matrix_element(t, tempWeight_1, 2)

                tempExcits(:, 2) = t

                nExcits = 2

            else if (near_zero(zeroWeight)) then
                ! only -2 start possible

                ! then change that to 2->1 start
                set_orb(t, 2 * st - 1)
                clr_orb(t, 2 * st)

                call setDeltaB(-2, t)

                call getDoubleMatrixElement(1, 2, -1, gen_type%L, gen_type%R, &
                                            bVal, 1.0_dp, x1_element=tempWeight_1)

                call encode_matrix_element(t, 0.0_dp, 1)
                call encode_matrix_element(t, tempWeight_1, 2)

                tempExcits(:, 1) = t

                nExcits = 1

            else
                ! smth went wrong
                call stop_all(this_routine, "both starting branch weights are 0")

            end if
        end select

    end subroutine mixedFullStart

    subroutine calcFullStartRaising(ilut, csf_i, excitInfo, excitations, nExcits, &
                                    posSwitches, negSwitches)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcFullStartRaising"

        integer :: nMax, ierr, iOrb, start, ende, semi, gen
        integer(n_int) :: t(0:nifguga)
        real(dp) :: tempWeight, minusWeight, plusWeight, bVal, nOpen
        type(WeightObj_t) :: weights
        integer(n_int), allocatable :: tempExcits(:, :)

        ASSERT(isZero(ilut, excitInfo%fullStart))
        ASSERT(isProperCSF_ilut(ilut))

        ! create the fullStart
        start = excitInfo%fullStart
        ende = excitInfo%fullEnd
        semi = excitInfo%firstEnd
        gen = excitInfo%lastGen
        bVal = csf_i%B_real(semi)
        ! first have to allocate both arrays for the determinant and weight list
        ! worse then worst case for single excitations are 2^|i-j| excitations
        ! for a given CSF -> for now use that to allocate the arrays first
        ! update only need the number of open orbitals between i and j, and
        ! some additional room if 0/3 at start
        ! use already provided open orbital counting function.
        ! nMax = 2**(ende - start)
        nMax = 2 + 2**count_open_orbs_ij(csf_i, start, ende)
        allocate(tempExcits(0:nifguga, nMax), stat=ierr)

        t = ilut

        ! have to make the switch and store the first matrix element in it too
        ! additionally also already calculate the sign coming from the
        ! pseudo double excitation which only depends on the number of open
        ! orbitals in the overlap region
        nOpen = real(count_open_orbs_ij(csf_i, start, semi - 1), dp)

        ! set 0->3
        set_orb(t, 2 * start)
        set_orb(t, 2 * start - 1)

        ! could encode matrix element here, but do it later for more efficiency
        ! call encode_part_sign(t,Root2 * (-1.0)**nMax, 1)

        ! also dont need to store it here, also later
        ! tempExcits(:,1) = t

        ! since only the 0 branch is taken, where there are now changes in
        ! the stepvector and the matrix element is just a sign i can directly
        ! go to the lowering semi stop
        ! i could write it specifically here, and i should for efficiency

        ! have to calc. weights here, which are only the normel single
        ! excitation weights
        weights = init_singleWeight(csf_i, ende)
        minusWeight = weights%proc%minus(negSwitches(semi), csf_i%B_real(semi), weights%dat)
        plusWeight = weights%proc%plus(posSwitches(semi), csf_i%B_real(semi), weights%dat)

        ASSERT(.not. isZero(ilut, semi))
        ASSERT(minusWeight + plusWeight > 0.0_dp)

        call encode_matrix_element(t, 0.0_dp, 2)

        select case (csf_i%stepvector(semi))
        case (3)
            ! have to check a few things with 3 semi stop

            ! did something wrong with matrix elements before here.. the 3
            ! semi-stop has an aditional sign...but all matrix elements
            ! the same..
            call encode_matrix_element(t, (-1.0_dp)**(nOpen + 1.0_dp), 1)

            ! bullshit...
            if (minusWeight > 0.0_dp .and. plusWeight > 0.0_dp) then
                ! both +1 and -1 excitations possible
                ! do 3->1 first -1 branch first
                clr_orb(t, 2 * semi)

                call setDeltaB(-1, t)

                tempExcits(:, 1) = t

                ! then do 3->2: +1 branch(change from already set 1
                clr_orb(t, 2 * semi - 1)
                set_orb(t, 2 * semi)

                call setDeltaB(1, t)

                tempExcits(:, 2) = t

                nExcits = 2

            else if (near_zero(plusWeight)) then
                ! only -1 branch possible
                ! do 3->1 first -1 branch first
                clr_orb(t, 2 * semi)

                call setDeltaB(-1, t)

                tempExcits(:, 1) = t

                nExcits = 1
            else if (near_zero(minusWeight)) then
                ! only +1 branch possible
                ! then do 3->2: +1 branch
                clr_orb(t, 2 * semi - 1)

                call getDoubleMatrixElement(2, 3, 0, gen_type%R, gen_type%R, bVal, 1.0_dp, tempWeight)

                call encode_matrix_element(t, Root2 * tempWeight * ((-1.0_dp)**nOpen), 1)

                call setDeltaB(1, t)

                tempExcits(:, 1) = t

                nExcits = 1

            else
                ! something went wrong, or i have to cancel exciation.. todo how
                call stop_all(this_routine, "something went wrong. shouldnt be here!")

            end if

        case (1)
            ! only one excitation possible, which also has to have
            ! non-zero weight or otherwise i wouldnt even be here
            ! and reuse already provided t
            ! 1 -> 0
            clr_orb(t, 2 * semi - 1)

            call getDoubleMatrixElement(0, 1, 0, gen_type%R, gen_type%R, bVal, 1.0_dp, tempWeight)

            call setDeltaB(1, t)

            ! encode fullstart contribution and pseudo overlap region here
            ! too in one go. one additional -1 due to semistop

            ! no Root2 here, since it cancels with t from matEle,
            call encode_matrix_element(t, Root2 * tempWeight * ((-1.0_dp)**(nOpen)), 1)

            tempExcits(:, 1) = t

            nExcits = 1

        case (2)
            ! here we have
            ! 2 -> 0
            clr_orb(t, 2 * semi)

            call getDoubleMatrixElement(0, 2, 0, gen_type%R, gen_type%R, bVal, 1.0_dp, tempWeight)

            call setDeltaB(-1, t)

            ! encode fullstart contribution and pseudo overlap region here
            ! too in one go. one additional -1 due to semistop

            ! no Root2 here, since it cancels with t from matEle,
            call encode_matrix_element(t, Root2 * tempWeight * ((-1.0_dp)**(nOpen)), 1)

            tempExcits(:, 1) = t

            nExcits = 1

        end select

        ! and then we have to do just a regular single excitation
        excitInfo%currentGen = excitInfo%lastGen

        do iOrb = semi + 1, ende - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! and do a regular end step
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

        ! that should be it...

    end subroutine calcFullStartRaising

    subroutine calcFullStartLowering(ilut, csf_i, excitInfo, excitations, nExcits, &
                                     posSwitches, negSwitches)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcFullStartLowering"

        integer :: nMax, iOrb, start, ende, semi, gen
        integer(n_int) :: t(0:nifguga)
        real(dp) :: tempWeight, minusWeight, plusWeight, bVal, nOpen
        type(WeightObj_t) :: weights
        integer(n_int), allocatable :: tempExcits(:, :)

        ASSERT(isThree(ilut, excitInfo%fullStart))
        ASSERT(isProperCSF_ilut(ilut))

        ! create the fullStart
        start = excitInfo%fullStart
        ende = excitInfo%fullEnd
        semi = excitInfo%firstEnd
        gen = excitInfo%firstGen
        bVal = csf_i%B_real(semi)
        ! first have to allocate both arrays for the determinant and weight list
        ! worse then worst case for single excitations are 2^|i-j| excitations
        ! for a given CSF -> for now use that to allocate the arrays first
        ! update only need the number of open orbitals between i and j, and
        ! some additional room if 0/3 at start
        ! use already provided open orbital counting function.
        ! nMax = 2**(ende - start)
        nMax = 2 + 2**count_open_orbs_ij(csf_i, start, ende)
        allocate(tempExcits(0:nifguga, nMax))

        t = ilut

        ! have to make the switch and store the first matrix element in it too
        ! additionally also already calculate the sign coming from the
        ! pseudo double excitation which only depends on the number of open
        ! orbitals in the overlap region
        ! just also count the semi here to take that additional -sign into account
        nOpen = real(count_open_orbs_ij(csf_i, start, semi), dp)

        ! set 3->0
        clr_orb(t, 2 * start)
        clr_orb(t, 2 * start - 1)

        ! since only the 0 branch is taken, where there are now changes in
        ! the stepvector and the matrix element is just a sign i can directly
        ! go to the lowering semi stop
        ! i could write it specifically here, and i should for efficiency

        ! have to calc. weights here, which are only the normel single
        ! excitation weights
        weights = init_singleWeight(csf_i, ende)
        minusWeight = weights%proc%minus(negSwitches(semi), bVal, weights%dat)
        plusWeight = weights%proc%plus(posSwitches(semi), bVal, weights%dat)

        ASSERT(.not. isThree(ilut, semi))
        ASSERT(minusWeight + plusWeight > 0.0_dp)

        call encode_matrix_element(t, 0.0_dp, 2)

        select case (csf_i%stepvector(semi))
        case (0)
            ! have to check a few things with 0 start
            if (minusWeight > 0.0_dp .and. plusWeight > 0.0_dp) then
                ! both +1 and -1 excitations possible
                ! do 0->1 first -1 branch first
                set_orb(t, 2 * semi - 1)

                call getDoubleMatrixElement(1, 0, 0, gen_type%L, gen_type%L, bVal, 1.0_dp, tempWeight)

                call encode_matrix_element(t, Root2 * tempWeight * ((-1.0_dp)**nOpen), 1)

                call setDeltaB(-1, t)

                tempExcits(:, 1) = t

                ! then do 0->2: +1 branch
                clr_orb(t, 2 * semi - 1)
                set_orb(t, 2 * semi)

                call getDoubleMatrixElement(2, 0, 0, gen_type%L, gen_type%L, bVal, 1.0_dp, tempWeight)

                call encode_matrix_element(t, Root2 * tempWeight * ((-1.0_dp)**nOpen), 1)

                call setDeltaB(1, t)

                tempExcits(:, 2) = t

                nExcits = 2

            else if (near_zero(plusWeight)) then
                ! only -1 branch possible
                ! do 0->1 first -1 branch first
                set_orb(t, 2 * semi - 1)

                call getDoubleMatrixElement(1, 0, 0, gen_type%L, gen_type%L, bVal, 1.0_dp, tempWeight)

                call encode_matrix_element(t, Root2 * tempWeight * ((-1.0_dp)**nOpen), 1)

                call setDeltaB(-1, t)

                tempExcits(:, 1) = t

                nExcits = 1
            else if (near_zero(minusWeight)) then
                ! only +1 branch possible

                ! then do 0->2: +1 branch
                set_orb(t, 2 * semi)

                call getDoubleMatrixElement(2, 0, 0, gen_type%L, gen_type%L, bVal, 1.0_dp, tempWeight)

                call encode_matrix_element(t, Root2 * tempWeight * ((-1.0_dp)**nOpen), 1)

                call setDeltaB(1, t)

                tempExcits(:, 1) = t

                nExcits = 1

            else
                ! something went wrong, or i have to cancel exciation.. todo how
                call stop_all(this_routine, "something went wrong. shouldnt be here!")

            end if

        case (1)
            ! only one excitation possible, which also has to have
            ! non-zero weight or otherwise i wouldnt even be here
            ! and reuse already provided t
            ! 1 -> 3
            set_orb(t, 2 * semi)

            call setDeltaB(1, t)

            ! encode fullstart contribution and pseudo overlap region here
            ! too in one go. one additional -1 due to semistop
            call encode_matrix_element(t, ((-1.0_dp)**nOpen), 1)

            tempExcits(:, 1) = t

            nExcits = 1

        case (2)
            ! here we have
            ! 2 -> 3
            set_orb(t, 2 * semi - 1)

            call setDeltaB(-1, t)

            ! encode fullstart contribution and pseudo overlap region here
            ! too in one go. one additional -1 due to semistop
            call encode_matrix_element(t, ((-1.0_dp)**nOpen), 1)

            tempExcits(:, 1) = t

            nExcits = 1

        end select

        ! and then we have to do just a regular single excitation
        excitInfo%currentGen = excitInfo%lastGen

        do iOrb = semi + 1, ende - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! and do a regular end step
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

        ! that should be it...

    end subroutine calcFullStartLowering

    subroutine calcFullStopR2L(ilut, csf_i, excitInfo, excitations, nExcits, &
                               posSwitches, negSwitches, t_no_singles_opt)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        logical, intent(in), optional :: t_no_singles_opt

        ! not sure if single or double weight is necessary here...
        type(WeightObj_t) :: weights
        integer(n_int), allocatable :: tempExcits(:, :)
        integer(n_int) :: t(0:nifguga)
        integer :: iOrb, iEx, cnt, st, se, en, ierr
        real(dp) :: minusWeight, plusWeight, zeroWeight
        logical :: t_no_singles

        st = excitInfo%fullStart
        se = excitInfo%secondStart
        en = excitInfo%fullEnd

        if (present(t_no_singles_opt)) then
            t_no_singles = t_no_singles_opt
        else
            t_no_singles = .false.
        end if

        if (t_no_singles .and. csf_i%stepvector(en) == 3) then
            allocate(excitations(0, 0), stat=ierr)
            nExcits = 0
            return
        end if

        ! init weights
        weights = init_semiStartWeight(csf_i, se, en, negSwitches(se), &
                                       posSwitches(se), csf_i%B_real(se))

        excitInfo%currentGen = excitInfo%firstGen

        ! create start
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                               weights, tempExcits, nExcits)

        ! loop until semi-start
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! can i write a general raisingSemiStart function, to reuse in
        ! other types of excitations?
        ! have to init specific prob. weights
        weights = weights%ptr
        ! depending if there is a 3 at the fullend there can be no switching
        ! possible
        if (csf_i%stepvector(en) == 3) then
            zeroWeight = weights%proc%zero(0.0_dp, 0.0_dp, csf_i%B_real(se), weights%dat)
            ! is plus and minus weight zero then? i think so
            minusWeight = 0.0_dp
            plusWeight = 0.0_dp

            if (t_mixed_hubbard) then
                nExcits = 0
                allocate(excitations(0, 0), stat=ierr)
                return
            end if
        else
            minusWeight = weights%proc%minus(negSwitches(se), csf_i%B_real(se), weights%dat)
            plusWeight = weights%proc%plus(posSwitches(se), csf_i%B_real(se), weights%dat)
            zeroWeight = weights%proc%zero(negSwitches(se), posSwitches(se), &
                                           csf_i%B_real(se), weights%dat)
        end if
        ! do i have to give them posSwitches and negSwitches or could I
        ! just put in actual weight values?
        call calcLoweringSemiStart(ilut, csf_i, excitInfo, tempExcits, nExcits, &
                                   plusWeight, minusWeight, zeroWeight)

        if (csf_i%stepvector(en) == 3) then
            ! in mixed generator case there is no sign coming from the
            ! overlap region, so only finish up the excitation, by just
            ! giving it the correct matrix element

            ! and have to check if non-zero and store in final
            ! excitations list

            cnt = 1
            do iEx = 1, nExcits
                t = tempExcits(:, iEx)

                if (near_zero(extract_matrix_element(t, 1))) cycle

                ! also no change in stepvector in this case
                call update_matrix_element(t, Root2, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                tempExcits(:, cnt) = t
                cnt = cnt + 1
            end do

            nExcits = cnt - 1

            allocate(excitations(0:nifguga, nExcits), stat=ierr)
            excitations = tempExcits(:, 1:nExcits)

            deallocate(tempExcits)

        else
            ! then i have to do a proper double excitation until the end
            ! with the new weights also
            do iOrb = excitInfo%secondStart + 1, excitInfo%fullEnd - 1
                call doubleUpdate(ilut, csf_i, iOrb, excitInfo, weights, tempExcits, nExcits, &
                                  negSwitches, posSwitches)
            end do

            ! and then do a mixed fullstop. also write a general function for that
            call mixedFullStop(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations, &
                               t_no_singles)

        end if

    end subroutine calcFullStopR2L

    subroutine calcFullStopL2R(ilut, csf_i, excitInfo, excitations, nExcits, &
                               posSwitches, negSwitches, t_no_singles_opt)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        logical, intent(in), optional :: t_no_singles_opt

        ! not sure if single or double weight is necessary here...
        type(WeightObj_t) :: weights
        integer(n_int), allocatable :: tempExcits(:, :)
        integer(n_int) :: t(0:nifguga)
        integer :: iOrb, iEx, cnt, st, se, en, ierr
        real(dp) :: minusWeight, plusWeight, zeroWeight
        logical :: t_no_singles

        st = excitInfo%fullStart
        se = excitInfo%secondStart
        en = excitInfo%fullEnd

        if (present(t_no_singles_opt)) then
            t_no_singles = t_no_singles_opt
        else
            t_no_singles = .false.
        end if

        if (t_no_singles .and. csf_i%stepvector(en) == 3) then
            allocate(excitations(0, 0), stat=ierr)
            nExcits = 0
            return
        end if

        ! init weights
        weights = init_semiStartWeight(csf_i, se, en, negSwitches(se), &
                                       posSwitches(se), csf_i%B_real(se))

        excitInfo%currentGen = excitInfo%firstGen
        ! create start
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                               weights, tempExcits, nExcits)

        ! loop until semi-start
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! can i write a general raisingSemiStart function, to reuse in
        ! other types of excitations?
        ! have to init specific prob. weights todo
        ! do i have to give them posSwitches and negSwitches or could I
        ! just put in actual weight values?
        weights = weights%ptr
        if (csf_i%stepvector(en) == 3) then
            minusWeight = 0.0_dp
            plusWeight = 0.0_dp
            zeroWeight = weights%proc%zero(0.0_dp, 0.0_dp, csf_i%B_real(se), weights%dat)

            if (t_mixed_hubbard) then
                allocate(excitations(0, 0), stat=ierr)
                nExcits = 0
                return
            end if

        else
            minusWeight = weights%proc%minus(negSwitches(se), csf_i%B_real(se), weights%dat)
            plusWeight = weights%proc%plus(posSwitches(se), csf_i%B_real(se), weights%dat)
            zeroWeight = weights%proc%zero(negSwitches(se), posSwitches(se), &
                                           csf_i%B_real(se), weights%dat)
        end if

        call calcRaisingSemiStart(ilut, csf_i, excitInfo, tempExcits, nExcits, &
                                  plusWeight, minusWeight, zeroWeight)

        if (csf_i%stepvector(en) == 3) then
            ! in mixed generator case there is no sign coming from the
            ! overlap region, so only finish up the excitation, by just
            ! giving it the correct matrix element

            ! and have to check if non-zero and store in final
            ! excitations list

            cnt = 1
            do iEx = 1, nExcits
                t = tempExcits(:, iEx)

                if (near_zero(extract_matrix_element(t, 1))) cycle

                ! also no change in stepvector in this case
                call update_matrix_element(t, Root2, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                tempExcits(:, cnt) = t
                cnt = cnt + 1
            end do

            nExcits = cnt - 1

            allocate(excitations(0:nifguga, nExcits), stat=ierr)
            excitations = tempExcits(:, 1:nExcits)

            deallocate(tempExcits)

        else
            ! then i have to do a proper double excitation until the end
            ! with the new weights also
            do iOrb = excitInfo%secondStart + 1, excitInfo%fullEnd - 1
                call doubleUpdate(ilut, csf_i, iOrb, excitInfo, weights, tempExcits, nExcits, &
                                  negSwitches, posSwitches)
            end do

            ! and then do a mixed fullstop. also write a general function for that
            call mixedFullStop(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations, &
                               t_no_singles)

        end if

    end subroutine calcFullStopL2R

    subroutine mixedFullStop(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations, &
                             t_no_singles_opt)
        ! full stop routine for mixed generators
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(inout), allocatable :: tempExcits(:, :)
        integer, intent(inout) :: nExcits
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        logical, intent(in), optional :: t_no_singles_opt
        character(*), parameter :: this_routine = "mixedFullStop"

        integer :: iEx, ende, cnt, ierr, deltaB
        real(dp) :: bVal, tempWeight_0, tempWeight_1
        integer(n_int) :: t(0:nifguga)
        logical :: t_no_singles

        ASSERT(.not. isZero(ilut, excitInfo%fullEnd))
        ! not sure if i should check if ilut is not 3... since i could handle
        ! that differently

        if (present(t_no_singles_opt)) then
            t_no_singles = t_no_singles_opt
        else
            t_no_singles = .false.
        end if

        ende = excitInfo%fullEnd
        bVal = csf_i%B_real(ende)

        select case (csf_i%stepvector(ende))
        case (1)
            ! ending for -2 and 0 branches
            do iEx = 1, nExcits

                t = tempExcits(:, iEx)
                deltaB = getDeltaB(t)

                if (t_no_singles) then
                    ! check if the x0 is greater than 0, which indicates
                    ! purely delta-b = 0 route.., which should be disregarded
                    ! if we do not want to take singles into account.
                    if (.not. near_zero(extract_matrix_element(t, 1))) then
                        ! just to be sure this should also mean:
                        ASSERT(DeltaB == 0)

                        ! encode 0 matrix element, which then the loop at
                        ! the bottom should take care of and throw it out
                        call encode_matrix_element(t, 0.0_dp, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                        cycle

                    end if
                end if

                ! a +2 branch shouldnt arrive here
                ASSERT(deltaB /= 2)

                if (deltaB == 0) then
                    ! no change in stepvector so only calc matrix element
                    call getMixedFullStop(1, 1, 0, bVal, tempWeight_0, tempWeight_1)

                    ! finalize added up matrix elements in final step
                    tempWeight_0 = extract_matrix_element(t, 1) * tempWeight_0 + &
                                   extract_matrix_element(t, 2) * tempWeight_1

                else
                    ! -2 branch: change 1->2
                    set_orb(t, 2 * ende)
                    clr_orb(t, 2 * ende - 1)

                    ! matrix element in this case: only the x1 element is 1,
                    ! so just take the x1_element from the ilut
                    tempWeight_0 = extract_matrix_element(t, 2)

                end if
                ! and encode that into real part
                call encode_matrix_element(t, tempWeight_0, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                tempExcits(:, iEx) = t

            end do

        case (2)
            ! ending for 0 and +2 branches
            do iEx = 1, nExcits
                t = tempExcits(:, iEx)
                deltaB = getDeltaB(t)

                if (t_no_singles) then
                    ! check if the x0 is greater than 0, which indicates
                    ! purely delta-b = 0 route.., which should be disregarded
                    ! if we do not want to take singles into account.
                    if (.not. near_zero(extract_matrix_element(t, 1))) then
                        ! just to be sure this should also mean:
                        ASSERT(DeltaB == 0)

                        ! encode 0 matrix element, which then the loop at
                        ! the bottom should take care of and throw it out
                        call encode_matrix_element(t, 0.0_dp, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        tempExcits(:, iEx) = t

                        cycle

                    end if
                end if

                ASSERT(deltaB /= -2)

                if (deltaB == 0) then
                    call getMixedFullStop(2, 2, 0, bVal, tempWeight_0, tempWeight_1)

                    tempWeight_0 = extract_matrix_element(t, 1) * tempWeight_0 + &
                                   extract_matrix_element(t, 2) * tempWeight_1

                else
                    ! set 2-> 1 for +2 branch
                    set_orb(t, 2 * ende - 1)
                    clr_orb(t, 2 * ende)

                    tempWeight_0 = extract_matrix_element(t, 2)

                end if
                call encode_matrix_element(t, tempWeight_0, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                tempExcits(:, iEx) = t

            end do

        case (3)
            ! not sure if i ever access this function with 3 at end but
            ! nevertheless implement it for now...
            ! here only the 0 branches arrive, and whatever happend to switch
            ! at some point in the ovelap region has 0 matrix element
            do iEx = 1, nExcits
                t = tempExcits(:, iEx)
                deltaB = getDeltaB(t)

                ASSERT(deltaB == 0)

                call update_matrix_element(t, Root2, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                tempExcits(:, iEx) = t
            end do
        end select

        ! check again if there are maybe zero matrix elements
        cnt = 1
        do iEx = 1, nExcits
            if (near_zero(extract_matrix_element(tempExcits(:, iEx), 1))) cycle

            if (t_mixed_hubbard) then
                ! for mixed hubbard (and maybe other lattice systems,
                ! i do not want single excitations here (i guess.)
                if (.not. near_zero(extract_matrix_element(tempExcits(:, iEx), 1))) cycle
            end if

            tempExcits(:, cnt) = tempExcits(:, iEx)

            cnt = cnt + 1
        end do

        nExcits = cnt - 1

        ! do finishing up stuff..
        allocate(excitations(0:nifguga, nExcits), stat=ierr)

        excitations = tempExcits(:, 1:nExcits)

        deallocate(tempExcits)

    end subroutine mixedFullStop

    subroutine calcLoweringSemiStart(ilut, csf_i, excitInfo, tempExcits, nExcits, &
                                     plusWeight, minusWeight, zeroWeight)
        ! try at creating a reusable raising generator semi start
        ! just realize, for double excitations i have to save 2 matrix elements
        ! so hopefully a way is to use the imaginary matrix element storage
        ! for these kind of excitations! -> ask simon
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(inout) :: tempExcits(:, :)
        integer, intent(inout) :: nExcits
        real(dp), intent(in) :: plusWeight, minusWeight, zeroWeight
        character(*), parameter :: this_routine = "calcLoweringSemiStart"

        integer :: iEx, deltaB, s, cnt, en, fe
        real(dp) :: tempWeight_0, tempWeight_1, tempWeight, bVal
        integer(n_int) :: t(0:nifguga)

        ! at semi start with alike generator full stops only the deltaB=0
        ! are compatible, i hope this is correctly accounted by the weights..
        s = excitInfo%secondStart
        bVal = csf_i%B_real(s)
        en = excitInfo%fullEnd
        fe = excitInfo%firstEnd

        ! now if there is a LR(3) end it also restricts these type of
        ! exitations to pseudo doubles -> make distinction
        ! the second if is to only apply this to fullstop cases but also be
        ! able to use it fore general double excitations
        if (csf_i%stepvector(en) == 3 .and. fe == en) then
            ! here only swiches to 0 branch are possible
            select case (csf_i%stepvector(s))
            case (1)
                ! only -1 branches can lead to 0 branch
                ! not yet asssured by weights that only -1 branch is
                ! arriving -> exclude wrong excitations
                cnt = 0
                do iEx = 1, nExcits
                    if (getDeltaB(tempExcits(:, iEx)) == -1) then
                        t = tempExcits(:, iEx)

                        ! change 1 -> 0
                        clr_orb(t, 2 * s - 1)

                        call getDoubleMatrixElement(0, 1, -1, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight)

                        call update_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        call setDeltaB(0, t)

                        cnt = cnt + 1
                        tempExcits(:, cnt) = t
                    end if
                end do
                nExcits = cnt

            case (2)
                ! only +1 can lead to 0 branch
                ! not yet correctly accounted in weights, so an arriving +1
                ! branch is ensured... todo
                cnt = 0

                do iEx = 1, nExcits
                    if (getDeltaB(tempExcits(:, iEx)) == +1) then

                        t = tempExcits(:, iEx)

                        ! change 2-> 0
                        clr_orb(t, 2 * s)

                        call getDoubleMatrixElement(0, 2, 1, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight)

                        call update_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        call setDeltaB(0, t)

                        cnt = cnt + 1
                        tempExcits(:, cnt) = t
                    end if

                end do
                nExcits = cnt

            case (3)
                ! has to be 3 in this generator combination.
                ASSERT(isThree(ilut, s))

                ! switches to 0 branch always possible
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)

                    deltaB = getDeltaB(t)

                    if (deltaB == -1) then

                        ! change 3 -> 2
                        clr_orb(t, 2 * s - 1)

                        call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight)

                    else
                        ! change 3->1
                        clr_orb(t, 2 * s)

                        call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight)

                    end if

                    call update_matrix_element(t, tempWeight, 1)
                    call encode_matrix_element(t, 0.0_dp, 2)

                    call setDeltaB(0, t)

                    tempExcits(:, iEx) = t

                end do
            end select

        else

            select case (csf_i%stepvector(s))
            case (1)
                ! always works if weights are fitting, check possibs.

                ! think i can do it generally for both... since its the same
                ! operations only... only need to check if atleast on of the
                ! relevant weights is nonzero

                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)

                    deltaB = getDeltaB(t)

                    ! change 1 -> 0
                    clr_orb(t, 2 * s - 1)

                    call setDeltaB(deltaB + 1, t)

                    call getDoubleMatrixElement(0, 1, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, &
                                                excitInfo%order, tempWeight_0, tempWeight_1)

                    call encode_matrix_element(t, tempWeight_1 * &
                                               extract_matrix_element(t, 1), 2)

                    call update_matrix_element(t, tempWeight_0, 1)

                    tempExcits(:, iEx) = t
                end do

            case (2)
                ! should also work generally, since if a weight is zero the
                ! non compatible branch shouldnt even arrive here

                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)

                    deltaB = getDeltaB(t)

                    ! change 2 -> 0
                    clr_orb(t, 2 * s)

                    call setDeltaB(deltaB - 1, t)

                    call getDoubleMatrixElement(0, 2, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, &
                                                excitInfo%order, tempWeight_0, tempWeight_1)

                    call encode_matrix_element(t, tempWeight_1 * &
                                               extract_matrix_element(t, 1), 2)

                    call update_matrix_element(t, tempWeight_0, 1)

                    tempExcits(:, iEx) = t
                end do

            case (3)
                ! has to be 3 in lowering case

                ! update: for a fulldouble excitation the 0-weight
                ! is not always > 0 dependending on the semistop and
                ! fullend values!
                ! so do it new!
                if ((.not. near_zero(minusWeight)) .and. (.not. near_zero(plusWeight)) &
                    .and. (.not. near_zero(zeroWeight))) then
                    ! all branches possible
                    ! how to most efficiently do all that...
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        tempWeight = extract_matrix_element(t, 1)

                        deltaB = getDeltaB(t)

                        ! first do 3 -> 1
                        clr_orb(t, 2 * s)
                        call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, &
                                                    tempWeight_0, tempWeight_1)

                        call update_matrix_element(t, tempWeight_0, 1)
                        call encode_matrix_element(t, tempWeight * tempWeight_1, 2)

                        call setDeltaB(deltaB - 1, t)

                        tempExcits(:, iEx) = t

                        ! then do 3->2
                        nExcits = nExcits + 1
                        clr_orb(t, 2 * s - 1)
                        set_orb(t, 2 * s)

                        call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, &
                                                    tempWeight_0, tempWeight_1)

                        call encode_matrix_element(t, tempWeight_0 * tempWeight, 1)
                        call encode_matrix_element(t, tempWeight_1 * tempWeight, 2)

                        call setDeltaB(deltaB + 1, t)

                        tempExcits(:, nExcits) = t
                    end do

                else if ((.not. near_zero(minusWeight)) .and. (.not. near_zero(plusWeight)) &
                         .and. near_zero(zeroWeight)) then
                    ! cont. 0 branches not possible
                    ! here a arriving -1 branch can only become a -2 branch
                    ! and a +1 branch a +2 branch
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (deltaB == -1) then
                            ! then there is a 3 -> 1 switch
                            clr_orb(t, 2 * s)

                            call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order, x1_element=tempWeight_1)

                            call setDeltaB(deltaB - 1, t)

                        else
                            ! there is a 3 -> 2 switch
                            clr_orb(t, 2 * s - 1)

                            call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order, x1_element=tempWeight_1)

                            call setDeltaB(deltaB + 1, t)
                        end if

                        call encode_matrix_element(t, extract_matrix_element(t, 1) * &
                                                   tempWeight_1, 2)
                        call encode_matrix_element(t, 0.0_dp, 1)

                        tempExcits(:, iEx) = t
                    end do

                else if ((.not. near_zero(minusWeight)) .and. near_zero(plusWeight) &
                         .and. (.not. near_zero(zeroWeight))) then
                    ! cont. + branches not possible
                    ! so for an arriving -1 branch both -2 and 0 branch are
                    ! possible but for a +1 branch only the 0 branch
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        ! make all possible 3->1 excitation first
                        deltaB = getDeltaB(t)

                        clr_orb(t, 2 * s)

                        ! todo: see how deltaB access correctly works...
                        ! have to use new deltaB here i Think!!!
                        call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight_0, tempWeight_1)

                        ! need unchanged matrix element for branching
                        tempWeight = extract_matrix_element(t, 1)

                        call encode_matrix_element(t, tempWeight * &
                                                   tempWeight_1, 2)
                        call update_matrix_element(t, tempWeight_0, 1)

                        call setDeltaB(deltaB - 1, t)

                        tempExcits(:, iEx) = t

                        ! and for +1 branch there is also a switch

                        if (deltaB == -1) then
                            ! have to switch from previuosly switched 1 -> 2
                            set_orb(t, 2 * s)
                            clr_orb(t, 2 * s - 1)

                            call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, &
                                                        excitInfo%order, tempWeight_0, tempWeight_1)

                            call encode_matrix_element(t, tempWeight_0 * &
                                                       tempWeight, 1)
                            call encode_matrix_element(t, tempWeight_1 * &
                                                       tempWeight, 2)

                            call setDeltaB(0, t)

                            nExcits = nExcits + 1
                            tempExcits(:, nExcits) = t
                        end if
                    end do

                else if ((.not. near_zero(minusWeight)) .and. near_zero(plusWeight) &
                         .and. near_zero(zeroWeight)) then
                    ! only - branch possible
                    ! so only a -1 arriving branch can cont. to a -2 branch
                    ! hopefully the weights coming before handle that correctly
                    ! so no +1 branch arrives here -> check!
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (deltaB == +1) then
                            call print_indices(excitInfo)
                            call write_guga_list(stdout, tempExcits(:, 1:nExcits))
                            ASSERT(deltaB /= 1)
                        end if

                        ! do the 3 -> 1 switch to -2 branch
                        clr_orb(t, 2 * s)

                        call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bval, excitInfo%order, x1_element=tempWeight_1)

                        call setDeltaB(deltaB - 1, t)

                        call encode_matrix_element(t, extract_matrix_element(t, 1) * &
                                                   tempWeight_1, 2)
                        call encode_matrix_element(t, 0.0_dp, 1)

                        tempExcits(:, iEx) = t
                    end do

                else if (near_zero(minusWeight) .and. (.not. near_zero(plusWeight)) &
                         .and. (.not. near_zero(zeroWeight))) then
                    ! cont. - branch not possible
                    ! so both options are possible for an incoming +1 branch
                    ! and only the 0 branch for the -1
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        ! make all possible 3->2 excitation first
                        deltaB = getDeltaB(t)

                        clr_orb(t, 2 * s - 1)

                        ! todo: see how deltaB access correctly works...
                        ! have to use new deltaB here i Think!!!
                        call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight_0, tempWeight_1)

                        ! need unchanged matrix element for branching
                        tempWeight = extract_matrix_element(t, 1)

                        call encode_matrix_element(t, tempWeight_0 * &
                                                   tempWeight, 1)
                        call encode_matrix_element(t, tempWeight_1 * &
                                                   tempWeight, 2)

                        call setDeltaB(deltaB + 1, t)

                        tempExcits(:, iEx) = t

                        ! and for +1 branch there is also a switch

                        if (deltaB == +1) then
                            ! have to switch from previuosly switched 2 -> 1
                            clr_orb(t, 2 * s)
                            set_orb(t, 2 * s - 1)

                            call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, &
                                                        excitInfo%order, tempWeight_0, tempWeight_1)

                            call encode_matrix_element(t, tempWeight_0 * &
                                                       tempWeight, 1)
                            call encode_matrix_element(t, tempWeight_1 * &
                                                       tempWeight, 2)

                            call setDeltaB(0, t)

                            nExcits = nExcits + 1
                            tempExcits(:, nExcits) = t
                        end if
                    end do

                else if (near_zero(minusWeight) .and. (.not. near_zero(plusWeight)) &
                         .and. near_zero(zeroWeight)) then
                    ! only + possible
                    ! only an arriving +1 branch can go on.. so check and
                    ! hope there is no -1 branch arriving
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (deltaB == -1) then
                            call print_indices(excitInfo)
                            call write_guga_list(stdout, tempExcits(:, 1:nExcits))
                            ASSERT(deltaB /= -1)
                        end if

                        ! do the 3 -> 2 switch to the +2 branch
                        clr_orb(t, 2 * s - 1)

                        call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, x1_element=tempWeight_1)

                        call setDeltaB(deltaB + 1, t)

                        call encode_matrix_element(t, extract_matrix_element(t, 1) * &
                                                   tempWeight_1, 2)
                        call encode_matrix_element(t, 0.0_dp, 1)

                        tempExcits(:, iEx) = t
                    end do

                else if (near_zero(minusWeight) .and. near_zero(plusWeight) &
                         .and. (.not. near_zero(zeroWeight))) then
                    ! only 0 branch possible
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        deltaB = getDeltaB(t)

                        if (deltaB == -1) then
                            ! -1 branch 3->2
                            clr_orb(t, 2 * s - 1)

                            call getDoubleMatrixElement(2, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order, tempWeight_0, &
                                                        tempWeight_1)

                        else
                            ! +1 branch: 3->1
                            clr_orb(t, 2 * s)

                            call getDoubleMatrixElement(1, 3, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order, tempWeight_0, &
                                                        tempWeight_1)

                        end if

                        call setDeltaB(0, t)

                        call encode_matrix_element(t, tempWeight_1 * &
                                                   extract_matrix_element(t, 1), 2)
                        call update_matrix_element(t, tempWeight_0, 1)

                        tempExcits(:, iEx) = t
                    end do

                else
                    ! something went wrong -> no branches possible
                    call stop_all(this_routine, " no branches possible in raising semistart!")
                end if

            end select
        end if

    end subroutine calcLoweringSemiStart

    subroutine calcRaisingSemiStart(ilut, csf_i, excitInfo, tempExcits, nExcits, &
                                    plusWeight, minusWeight, zeroWeight)
        ! try at creating a reusable raising generator semi start
        ! just realize, for double excitations i have to save 2 matrix elements
        ! so hopefully a way is to use the imaginary matrix element storage
        ! for these kind of excitations! -> ask simon
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(inout) :: tempExcits(:, :)
        integer, intent(inout) :: nExcits
        real(dp), intent(in) :: plusWeight, minusWeight, zeroWeight
        character(*), parameter :: this_routine = "calcRaisingSemiStart"

        integer :: iEx, deltaB, s, cnt, en, fe
        real(dp) :: tempWeight_0, tempWeight_1, tempWeight, bVal
        integer(n_int) :: t(0:nifguga)

        ! at semi start with alike generator full stops only the deltaB=0
        ! are compatible, i hope this is correctly accounted by the weights..
        s = excitInfo%secondStart
        bVal = csf_i%B_real(s)
        fe = excitInfo%firstEnd
        en = excitInfo%fullEnd

        ! now if there is a LR(3) end it also restricts these type of
        ! exitations to pseudo doubles -> make distinction
        ! to use it for general double excitation i also have to check if
        ! i am dealing with a full stop. otherwise stepvalue doesnt matter
        ! at end
        ! could avoid that is Three(end) now since included in weights...
        if (csf_i%stepvector(en) == 3 .and. en == fe) then
            ! here only swiches to 0 branch are possible
            select case (csf_i%stepvector(s))
            case (1)
                ! only -1 branches can lead to 0 branch
                ! not yet assured correctly by weights that only -1 branches
                ! arrive...
                cnt = 0
                do iEx = 1, nExcits
                    if (getDeltaB(tempExcits(:, iEx)) == -1) then
                        t = tempExcits(:, iEx)

                        ! change 1 -> 3
                        set_orb(t, 2 * s)

                        call getDoubleMatrixElement(3, 1, -1, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight)

                        call update_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        call setDeltaB(0, t)

                        cnt = cnt + 1
                        tempExcits(:, cnt) = t
                    end if
                end do
                nExcits = cnt

            case (2)
                ! only +1 can lead to 0 branch
                cnt = 0
                do iEx = 1, nExcits
                    if (getDeltaB(tempExcits(:, iEx)) == +1) then

                        t = tempExcits(:, iEx)

                        ! change 2-> 3
                        set_orb(t, 2 * s - 1)

                        call getDoubleMatrixElement(3, 2, 1, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight)

                        call update_matrix_element(t, tempWeight, 1)
                        call encode_matrix_element(t, 0.0_dp, 2)

                        call setDeltaB(0, t)

                        cnt = cnt + 1
                        tempExcits(:, cnt) = t
                    end if
                end do
                nExcits = cnt

            case (0)
                ! has to be 0 in this generator combination.
                ASSERT(isZero(ilut, s))

                ! switches to 0 branch always possible
                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)

                    deltaB = getDeltaB(t)

                    if (deltaB == -1) then

                        ! change 0->2
                        set_orb(t, 2 * s)

                        call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight)

                    else
                        ! change 0->1
                        set_orb(t, 2 * s - 1)

                        call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight)

                    end if

                    call update_matrix_element(t, tempWeight, 1)
                    call encode_matrix_element(t, 0.0_dp, 2)

                    call setDeltaB(0, t)

                    tempExcits(:, iEx) = t

                end do
            end select

        else

            select case (csf_i%stepvector(s))
            case (1)
                ! always works if weights are fitting, check possibs.

                ! think i can do it generally for both... since its the same
                ! operations only... only need to check if atleast on of the
                ! relevant weights is nonzero
                ! since no choice here do not have to check...
                ! hope this implementation is correct ...

                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)

                    deltaB = getDeltaB(t)

                    ! change 1 -> 3
                    set_orb(t, 2 * s)

                    call setDeltaB(deltaB + 1, t)

                    call getDoubleMatrixElement(3, 1, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, &
                                                excitInfo%order, tempWeight_0, tempWeight_1)

                    call encode_matrix_element(t, tempWeight_1 * &
                                               extract_matrix_element(t, 1), 2)

                    call update_matrix_element(t, tempWeight_0, 1)

                    tempExcits(:, iEx) = t
                end do

            case (2)
                ! should also work generally, since if a weight is zero the
                ! non compatible branch shouldnt even arrive here

                do iEx = 1, nExcits
                    t = tempExcits(:, iEx)

                    deltaB = getDeltaB(t)

                    ! change 2 -> 3
                    set_orb(t, 2 * s - 1)

                    call setDeltaB(deltaB - 1, t)

                    call getDoubleMatrixElement(3, 2, deltaB, excitInfo%gen1, &
                                                excitInfo%gen2, bVal, &
                                                excitInfo%order, tempWeight_0, tempWeight_1)

                    call encode_matrix_element(t, tempWeight_1 * &
                                               extract_matrix_element(t, 1), 2)

                    call update_matrix_element(t, tempWeight_0, 1)

                    tempExcits(:, iEx) = t
                end do

            case (0)
                ! has to be 0 in raising case

                ! update: for a fulldouble excitation the 0-weight
                ! is not always > 0 dependending on the semistop and
                ! fullend values!
                ! so do it new!
                if ((.not. near_zero(minusWeight)) .and. (.not. near_zero(plusWeight)) &
                    .and. (.not. near_zero(zeroWeight))) then
                    ! all branches possible
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        tempWeight = extract_matrix_element(t, 1)

                        deltaB = getDeltaB(t)

                        ! first do 0 -> 1
                        set_orb(t, 2 * s - 1)
                        call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, &
                                                    tempWeight_0, tempWeight_1)

                        call encode_matrix_element(t, tempWeight_0 * &
                                                   tempWeight, 1)
                        call encode_matrix_element(t, tempWeight_1 * &
                                                   tempWeight, 2)

                        call setDeltaB(deltaB - 1, t)

                        tempExcits(:, iEx) = t

                        ! then do 0->2
                        nExcits = nExcits + 1
                        clr_orb(t, 2 * s - 1)
                        set_orb(t, 2 * s)

                        call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, &
                                                    tempWeight_0, tempWeight_1)

                        call encode_matrix_element(t, tempWeight_0 * &
                                                   tempWeight, 1)
                        call encode_matrix_element(t, tempWeight_1 * &
                                                   tempWeight, 2)

                        call setDeltaB(deltaB + 1, t)

                        tempExcits(:, nExcits) = t
                    end do

                else if ((.not. near_zero(minusWeight)) .and. (.not. near_zero(plusWeight)) &
                         .and. near_zero(zeroWeight)) then
                    ! cont. 0 branches not possible
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (deltaB == -1) then
                            ! do the 0 > 1 switch to the -2 branch
                            set_orb(t, 2 * s - 1)

                            call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order, x1_element=tempWeight_1)

                            call setDeltaB(deltaB - 1, t)

                        else
                            ! do the 0 -> 2 switch to +2 branch
                            set_orb(t, 2 * s)

                            call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order, x1_element=tempWeight_1)

                            call setDeltaB(deltaB + 1, t)
                        end if

                        call encode_matrix_element(t, extract_matrix_element(t, 1) * &
                                                   tempWeight_1, 2)
                        call encode_matrix_element(t, 0.0_dp, 1)

                        tempExcits(:, iEx) = t
                    end do

                else if ((.not. near_zero(minusWeight)) .and. near_zero(plusWeight) &
                         .and. (.not. near_zero(zeroWeight))) then
                    ! cont. + branches not possible
                    ! +2 excitations not possible
                    ! only branching for -1 branch
                    ! but excitation 0->1 everytime possible for both
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        ! make all possible 0->1 excitation first
                        deltaB = getDeltaB(t)

                        set_orb(t, 2 * s - 1)

                        ! todo: see how deltaB access correctly works...
                        ! have to use new deltaB here i Think!!!
                        call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight_0, tempWeight_1)

                        ! need unchanged matrix element for branching
                        tempWeight = extract_matrix_element(t, 1)

                        call encode_matrix_element(t, tempWeight_0 * &
                                                   tempWeight, 1)
                        call encode_matrix_element(t, tempWeight_1 * &
                                                   tempWeight, 2)

                        call setDeltaB(deltaB - 1, t)

                        tempExcits(:, iEx) = t

                        ! and for -1 branch there is also a switch

                        if (deltaB == -1) then
                            ! have to switch from previuosly switched 1 -> 2
                            set_orb(t, 2 * s)
                            clr_orb(t, 2 * s - 1)

                            call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, &
                                                        excitInfo%order, tempWeight_0, tempWeight_1)

                            call encode_matrix_element(t, tempWeight_0 * &
                                                       tempWeight, 1)
                            call encode_matrix_element(t, tempWeight_1 * &
                                                       tempWeight, 2)

                            call setDeltaB(0, t)

                            nExcits = nExcits + 1
                            tempExcits(:, nExcits) = t
                        end if
                    end do

                else if ((.not. near_zero(minusWeight)) .and. near_zero(plusWeight) &
                         .and. near_zero(zeroWeight)) then
                    ! only - branch possible
                    ! so ensure no +1 branch arrives...
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (deltaB == 1) then
                            call print_indices(excitInfo)
                            call write_guga_list(stdout, tempExcits(:, 1:nExcits))
                            ASSERT(deltaB /= 1)
                        end if

                        ! do the 0 > 1 switch to -2 branch
                        set_orb(t, 2 * s - 1)

                        call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, x1_element=tempWeight_1)

                        call setDeltaB(deltaB - 1, t)

                        call encode_matrix_element(t, extract_matrix_element(t, 1) * &
                                                   tempWeight_1, 2)
                        call encode_matrix_element(t, 0.0_dp, 1)

                        tempExcits(:, iEx) = t
                    end do

                else if (near_zero(minusWeight) .and. (.not. near_zero(plusWeight)) &
                         .and. (.not. near_zero(zeroWeight))) then
                    ! cont. - branch not possible
                    ! -2 excitations not possible
                    ! only branching for +1 branch
                    ! 0->2 excitation works for both
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        ! make all possible 0->2 excitation first
                        deltaB = getDeltaB(t)

                        set_orb(t, 2 * s)

                        ! todo: see how deltaB access correctly works...
                        ! have to use new deltaB here i Think!!!
                        call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, &
                                                    excitInfo%order, tempWeight_0, tempWeight_1)

                        ! need unchanged matrix element for branching
                        tempWeight = extract_matrix_element(t, 1)

                        call encode_matrix_element(t, tempWeight_0 * &
                                                   tempWeight, 1)
                        call encode_matrix_element(t, tempWeight_1 * &
                                                   tempWeight, 2)

                        call setDeltaB(deltaB + 1, t)

                        tempExcits(:, iEx) = t

                        ! and for +1 branch there is also a switch

                        if (deltaB == +1) then
                            ! have to switch from previuosly switched 2 -> 1
                            clr_orb(t, 2 * s)
                            set_orb(t, 2 * s - 1)

                            call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, &
                                                        excitInfo%order, tempWeight_0, tempWeight_1)

                            call encode_matrix_element(t, tempWeight_0 * &
                                                       tempWeight, 1)
                            call encode_matrix_element(t, tempWeight_1 * &
                                                       tempWeight, 2)

                            call setDeltaB(0, t)

                            nExcits = nExcits + 1
                            tempExcits(:, nExcits) = t
                        end if
                    end do

                else if (near_zero(minusWeight) .and. (.not. near_zero(plusWeight)) &
                         .and. near_zero(zeroWeight)) then
                    ! only + possible
                    ! so check that no -1 branch arrives
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)
                        deltaB = getDeltaB(t)

                        if (deltaB == -1) then
                            call print_indices(excitInfo)
                            call write_guga_list(stdout, tempExcits(:, 1:nExcits))
                            ASSERT(deltaB /= -1)
                        end if

                        ! do the 0 -> 2 switch to +2 branch
                        set_orb(t, 2 * s)

                        call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                    excitInfo%gen2, bVal, excitInfo%order, x1_element=tempWeight_1)

                        call setDeltaB(deltaB + 1, t)

                        call encode_matrix_element(t, extract_matrix_element(t, 1) * &
                                                   tempWeight_1, 2)
                        call encode_matrix_element(t, 0.0_dp, 1)

                        tempExcits(:, iEx) = t
                    end do

                else if (near_zero(minusWeight) .and. near_zero(plusWeight) &
                         .and. (.not. near_zero(zeroWeight))) then
                    ! only 0 branch possible
                    do iEx = 1, nExcits
                        t = tempExcits(:, iEx)

                        deltaB = getDeltaB(t)

                        if (deltaB == -1) then
                            ! -1 branch 0->2
                            set_orb(t, 2 * s)

                            call getDoubleMatrixElement(2, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order, tempWeight_0, &
                                                        tempWeight_1)

                        else
                            ! +1 branch: 0->1
                            set_orb(t, 2 * s - 1)

                            call getDoubleMatrixElement(1, 0, deltaB, excitInfo%gen1, &
                                                        excitInfo%gen2, bVal, excitInfo%order, tempWeight_0, &
                                                        tempWeight_1)

                        end if

                        call setDeltaB(0, t)

                        call encode_matrix_element(t, tempWeight_1 * &
                                                   extract_matrix_element(t, 1), 2)
                        call update_matrix_element(t, tempWeight_0, 1)

                        tempExcits(:, iEx) = t
                    end do
                else
                    ! something went wrong -> no branches possible
                    call stop_all(this_routine, " no branches possible in lowering semistart!")
                end if

            end select
        end if

    end subroutine calcRaisingSemiStart

    subroutine calcFullstopLowering(ilut, csf_i, excitInfo, excitations, nExcits, &
                                    posSwitches, negSwitches)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcFullstopLowering"

        ! not sure if single or double weight is necessary here...
        type(WeightObj_t) :: weights
        integer(n_int), allocatable :: tempExcits(:, :)
        integer(n_int) :: t(0:nifguga)
        integer :: iOrb, iEx, deltaB, cnt, ierr
        real(dp) :: tempWeight, sig, bVal

        ! create weight object here
        ! i think i only need single excitations weights here, since
        ! the semi stop in this case is like an end step...
        weights = init_singleWeight(csf_i, excitInfo%secondStart)

        ! create start
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                               weights, tempExcits, nExcits)

        excitInfo%currentGen = excitInfo%firstGen
        ! loop until semi-start
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! at semi start with alike generator full stops only the deltaB=0
        ! are compatible, i hope this is correctly accounted by the weights..
        iOrb = excitInfo%secondStart
        bVal = csf_i%B_real(iOrb)

        select case (csf_i%stepvector(iOrb))
        case (1)
            ! only delta -1 branches can lead to deltaB=0 in overlap here
            do iEx = 1, nExcits
                ! correct use of weight gens should exclude this, but check:

                ASSERT(getDeltaB(tempExcits(:, iEx)) == -1)

                t = tempExcits(:, iEx)

                ! change 1 -> 0
                clr_orb(t, 2 * iOrb - 1)

                call getDoubleMatrixElement(0, 1, -1, gen_type%L, gen_type%L, &
                                            bVal, excitInfo%order, tempWeight)

                call update_matrix_element(t, tempWeight, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                call setDeltaB(0, t)

                tempExcits(:, iEx) = t
            end do

        case (2)
            ! only deltaB=+1 branches lead to a 0 branch
            do iEx = 1, nExcits
                ASSERT(getDeltaB(tempExcits(:, iEx)) == +1)

                t = tempExcits(:, iEx)

                ! change 2 -> 0
                clr_orb(t, 2 * iOrb)

                call getDoubleMatrixElement(0, 2, +1, gen_type%L, gen_type%L, &
                                            bVal, 1.0_dp, tempWeight)

                call update_matrix_element(t, tempWeight, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                call setDeltaB(0, t)

                tempExcits(:, iEx) = t
            end do

        case (3)
            ! has to be a 3 in the lowering case.
            ASSERT(isThree(ilut, iOrb))
            ! -1 branches always go to 0 branch
            ! +1 branches might go to 0 branch if bValue high enough
            ! do not have to check weights here, since 0 branch has to have
            ! non-zero weight or it wouldnt even be compatible
            ! no! bullshit! switchin to 0 branch always works!
            ! both branches survive
            do iEx = 1, nExcits
                t = tempExcits(:, iEx)

                deltaB = getDeltaB(t)

                if (deltaB == -1) then

                    ! change 3->2
                    clr_orb(t, 2 * iOrb - 1)

                    call getDoubleMatrixElement(2, 3, deltaB, gen_type%L, gen_type%L, &
                                                bVal, 1.0_dp, tempWeight)

                else
                    ! change 3->1
                    clr_orb(t, 2 * iOrb)

                    call getDoubleMatrixElement(1, 3, deltaB, gen_type%L, gen_type%L, &
                                                bVal, 1.0_dp, tempWeight)

                end if

                call update_matrix_element(t, tempWeight, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                call setDeltaB(0, t)

                tempExcits(:, iEx) = t

            end do
        end select

        ! continue on with double excitation region, only the 0 branch
        ! valid here, where there is no change in stepvector and matrix
        ! element only a sign dependent on the number of open orbitals
        sig = (-1.0_dp)**real(count_open_orbs_ij(csf_i, excitInfo%secondStart + 1, &
                                                 excitInfo%fullEnd - 1), dp)

        ! do the ending
        ASSERT(isZero(ilut, excitInfo%fullEnd))

        iOrb = excitInfo%fullEnd

        cnt = 1

        do iEx = 1, nExcits
            t = tempExcits(:, iEx)

            if (near_zero(extract_matrix_element(t, 1))) cycle
            ! change 0 -> 3
            set_orb(t, 2 * iOrb)
            set_orb(t, 2 * iOrb - 1)

            ! have to cancel 0 weighted excitations

            ! also include the sign of the open orbitals here
            call update_matrix_element(t, sig * Root2, 1)

            tempExcits(:, cnt) = t
            cnt = cnt + 1

        end do

        nExcits = cnt - 1

        ! have to put it in excitations stupid!
        allocate(excitations(0:nifguga, nExcits), stat=ierr)
        excitations = tempExcits(:, 1:nExcits)
        deallocate(tempExcits)

    end subroutine calcFullstopLowering

    subroutine calcFullstopRaising(ilut, csf_i, excitInfo, excitations, nExcits, &
                                   posSwitches, negSwitches)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcFullstopRaising"

        ! not sure if single or double weight is necessary here...
        type(WeightObj_t) :: weights
        integer(n_int), allocatable :: tempExcits(:, :)
        integer(n_int) :: t(0:nifguga)
        integer :: iOrb, iEx, deltaB, cnt, ierr
        real(dp) :: tempWeight, sig, bVal

        ! create weight object here todo
        weights = init_singleWeight(csf_i, excitInfo%secondStart)

        excitInfo%currentGen = excitInfo%firstGen

        ! create start
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                               weights, tempExcits, nExcits)

        ! loop until semi-start
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! at semi start with alike generator full stops only the deltaB=0
        ! are compatible, i hope this is correctly accounted by the weights..
        iOrb = excitInfo%secondStart
        bVal = csf_i%B_real(iOrb)

        select case (csf_i%stepvector(iOrb))
        case (1)
            ! only delta -1 branches can lead to deltaB=0 in overlap here
            do iEx = 1, nExcits
                ! correct use of weight gens should exclude this, but check:
                ASSERT(getDeltaB(tempExcits(:, iEx)) == -1)

                t = tempExcits(:, iEx)

                ! change 1 -> 3
                set_orb(t, 2 * iOrb)

                call getDoubleMatrixElement(3, 1, -1, gen_type%R, gen_type%R, bVal, &
                                            excitInfo%order, tempWeight)

                call update_matrix_element(t, tempWeight, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                call setDeltaB(0, t)

                tempExcits(:, iEx) = t
            end do

        case (2)
            ! only deltaB=+1 branches lead to a 0 branch
            do iEx = 1, nExcits
                ASSERT(getDeltaB(tempExcits(:, iEx)) == +1)

                t = tempExcits(:, iEx)

                ! change 2 -> 3
                set_orb(t, 2 * iOrb - 1)

                call getDoubleMatrixElement(3, 2, +1, gen_type%R, gen_type%R, bVal, &
                                            1.0_dp, tempWeight)

                call update_matrix_element(t, tempWeight, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                call setDeltaB(0, t)

                tempExcits(:, iEx) = t
            end do

        case (0)
            ! has to be a 0 in the raising case.
            ASSERT(isZero(ilut, iOrb))
            ! -1 branches always go to 0 branch
            ! +1 branches might go to 0 branch if bValue high enough
            ! do not have to check weights here, since 0 branch has to have
            ! non-zero weight or it wouldnt even be compatible
            ! no! bullshit! switchin to 0 branch always works!
            ! both branches survive
            do iEx = 1, nExcits
                t = tempExcits(:, iEx)

                deltaB = getDeltaB(t)

                if (deltaB == -1) then

                    ! change 0->2
                    set_orb(t, 2 * iOrb)

                    call getDoubleMatrixElement(2, 0, deltaB, gen_type%R, gen_type%R, bVal, &
                                                1.0_dp, tempWeight)

                else
                    ! change 0->1
                    set_orb(t, 2 * iOrb - 1)

                    call getDoubleMatrixElement(1, 0, deltaB, gen_type%R, gen_type%R, &
                                                bVal, 1.0_dp, tempWeight)

                end if

                call update_matrix_element(t, tempWeight, 1)
                call encode_matrix_element(t, 0.0_dp, 2)

                call setDeltaB(0, t)

                tempExcits(:, iEx) = t

            end do
        end select

        ! continue on with double excitation region, only the 0 branch
        ! valid here, where there is no change in stepvector and matrix
        ! element only a sign dependent on the number of open orbitals
        sig = (-1.0_dp)**real(count_open_orbs_ij(csf_i, excitInfo%secondStart + 1, &
                                                 excitInfo%fullEnd - 1), dp)

        ! do the ending
        ASSERT(isThree(ilut, excitInfo%fullEnd))

        iOrb = excitInfo%fullEnd

        cnt = 1
        do iEx = 1, nExcits
            t = tempExcits(:, iEx)

            if (near_zero(extract_matrix_element(t, 1))) cycle

            ! change 3 -> 0
            clr_orb(t, 2 * iOrb)
            clr_orb(t, 2 * iOrb - 1)

            ! also include the sign of the open orbitals here
            call update_matrix_element(t, sig * Root2, 1)

            tempExcits(:, cnt) = t
            cnt = cnt + 1

        end do

        nExcits = cnt - 1

        allocate(excitations(0:nifguga, nExcits), stat=ierr)
        excitations = tempExcits(:, 1:nExcits)
        deallocate(tempExcits)

    end subroutine calcFullstopRaising

    subroutine calcSingleOverlapLowering(ilut, csf_i, excitInfo, excitations, nExcits, &
                                         posSwitches, negSwitches)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcSingleOverlapLowering"

        integer(n_int), allocatable :: tempExcits(:, :)
        integer(n_int) ::  t(0:nifguga)
        integer :: i, iEx, deltaB, ss
        type(WeightObj_t) :: weights

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(excitInfo%typ == excit_type%single_overlap_lowering)
        ASSERT(excitInfo%firstGen == gen_type%L)
        ASSERT(excitInfo%lastGen == gen_type%L)

        excitInfo%currentGen = excitInfo%firstGen
        ! have to make specific single start to correctly adress the weights
        weights = init_singleOverlapLowering(csf_i, excitInfo%firstEnd, &
                                             excitInfo%fullEnd, negSwitches(excitInfo%firstEnd), posSwitches(excitInfo%firstEnd), &
                                             csf_i%B_real(excitInfo%firstEnd))

        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, weights, &
                               tempExcits, nExcits)

        ss = excitInfo%secondStart
        ! loop until overlap site
        do i = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleUpdate(ilut, csf_i, i, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! do special stuff at lowering site
        ! two lowerings

        ! has only forced switches at switch possibilities
        if (csf_i%stepvector(ss) == 1) then
            ! switch deltaB = -1 branches
            do iEx = 1, nExcits
                deltaB = getDeltaB(tempExcits(:, iEx))
                if (deltaB == -1) then
                    ! switch 1 - > 2
                    t = tempExcits(:, iEx)
                    clr_orb(t, 2 * ss - 1)
                    set_orb(t, 2 * ss)

                    call setDeltaB(1, t)

                    tempExcits(:, iEx) = t

                    ! no change in matrix elements
                end if
            end do

        else if (csf_i%stepvector(ss) == 2) then
            ! switch deltaB = +1
            do iEx = 1, nExcits
                deltaB = getDeltaB(tempExcits(:, iEx))
                if (deltaB == 1) then
                    ! switch 2 -> 1
                    t = tempExcits(:, iEx)
                    clr_orb(t, 2 * ss)
                    set_orb(t, 2 * ss - 1)

                    call setDeltaB(-1, t)

                    tempExcits(:, iEx) = t
                end if
            end do
        end if

        excitInfo%currentGen = excitInfo%lastGen

        ! update weights here?
        weights = weights%ptr

        ! continue with secon region normally
        do i = excitInfo%secondStart + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, csf_i, i, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)

        end do

        ! normal end step
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

    end subroutine calcSingleOverlapLowering

    subroutine calcSingleOverlapRaising(ilut, csf_i, excitInfo, excitations, &
                                        nExcits, posSwitches, negSwitches)
        ! special function to calculate all excitations for a single overlap
        ! double excitations with only raising generators
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcSingleOverlapRaising"

        integer(n_int), allocatable :: tempExcits(:, :)
        integer(n_int) :: t(0:nifguga)
        integer :: i, iEx, deltaB, se, en
        type(WeightObj_t) :: weightObj
        real(dp) :: plusWeight, minusWeight

        ASSERT(isProperCSF_ilut(ilut))
        ASSERT(excitInfo%typ == excit_type%single_overlap_raising)
        ASSERT(excitInfo%firstGen == gen_type%R)
        ASSERT(excitInfo%lastGen == gen_type%R)

        se = excitInfo%secondStart
        en = excitInfo%fullEnd
        ! create weight object here
        weightObj = init_singleOverlapRaising(csf_i, se, en, negSwitches(se), &
                                              posSwitches(se), csf_i%B_real(se))

        excitInfo%currentGen = excitInfo%firstGen
        ! have to make specific single start to correctly adress the weights
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                               weightObj, tempExcits, nExcits)

        ! loop until overlap site
        do i = excitInfo%fullStart + 1, se - 1
            call singleUpdate(ilut, csf_i, i, excitInfo, posSwitches, &
                              negSwitches, weightObj, tempExcits, nExcits)
        end do

        i = excitInfo%secondStart

        ! update weights
        weightObj = init_singleWeight(csf_i, en)
        plusWeight = weightObj%proc%plus(posSwitches(se), &
                                         csf_i%B_real(se), weightObj%dat)
        minusWeight = weightObj%proc%minus(negSwitches(se), &
                                           csf_i%B_real(se), weightObj%dat)

        ASSERT(plusWeight + minusWeight > 0.0_dp)

        ! do special stuff at single overlap raising site
        if (csf_i%stepvector(se) == 1) then
            ! in this case there should come a deltaB=+1 branch, nevertheless
            ! check for now..
            ! have to include probabilistic weights too... which are the normal
            ! single excitations weights at this point
            do iEx = 1, nExcits
                deltaB = getDeltaB(tempExcits(:, iEx))
                if (deltaB == 1) then
                    call stop_all(this_routine, &
                                  "there should not be a DeltaB= +1 at a RR(1) single overlap site")
                end if
                ! only need to check switching weight, since on track has to
                ! be greater than 0, to be even here TODO: think about that!
                ! not quite sure yeah.. since this is a switch possib..
                ! jsut to be save to a check like above
                ! update: to be save change that to check probs

                if (near_zero(minusWeight)) then
                    ! do only switch
                    t = tempExcits(:, iEx)
                    ! change 1 -> 2
                    set_orb(t, 2 * se)
                    clr_orb(t, 2 * se - 1)

                    call setDeltaB(1, t)

                    ! no change in matrix elements
                    tempExcits(:, iEx) = t

                else
                    ! staying on -1 doesnt change anything, so just check if
                    ! a switch is also possible in this case: and add at end

                    if (plusWeight > 0.0_dp) then
                        ! and all deltaB=-1 are branch possibs. so add an excitations
                        nExcits = nExcits + 1
                        t = tempExcits(:, iEx)
                        ! change 1 -> 2
                        set_orb(t, 2 * se)
                        clr_orb(t, 2 * se - 1)

                        call setDeltaB(1, t)

                        ! no change in matrix elements
                        tempExcits(:, nExcits) = t
                    end if
                end if
            end do

        else if (csf_i%stepvector(se) == 2) then
            ! in this case always a switch, and i b allows also a stay
            ! how are the probs here...

            if (near_zero(plusWeight)) then
                ! just switch
                ! only switches in place
                do iEx = 1, nExcits
                    if (getDeltaB(tempExcits(:, iEx)) == -1) then
                        call stop_all(this_routine, &
                                      "there should be no deltaB=-1 at RR(2) single overlap site")
                    end if

                    t = tempExcits(:, iEx)
                    ! change 2 -> 1
                    clr_orb(t, 2 * se)
                    set_orb(t, 2 * se - 1)

                    call setDeltaB(-1, t)

                    tempExcits(:, iEx) = t
                end do

            else
                ! staying doesnt change anything just check if switch is
                ! additionally possible

                if (minusWeight > 0.0_dp) then
                    ! stay and add switch
                    do iEx = 1, nExcits
                        ! still check delta B just to be sure
                        deltaB = getDeltaB(tempExcits(:, iEx))
                        if (deltaB == -1) then
                            call stop_all(this_routine, &
                                          "there should be no deltaB=-1 at RR(2) single overlap site")
                        end if
                        nExcits = nExcits + 1

                        t = tempExcits(:, iEx)
                        ! change 2 -> 1
                        clr_orb(t, 2 * se)
                        set_orb(t, 2 * se - 1)

                        call setDeltaB(-1, t)

                        tempExcits(:, nExcits) = t
                    end do
                end if
            end if
        end if

        ! continue with secon region normally
        ! have to reset weight object here to use "normal" single excitation
        ! also no change in generator type since alike generators
        excitInfo%currentGen = excitInfo%lastGen

        do i = excitInfo%secondStart + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, csf_i, i, excitInfo, posSwitches, negSwitches, &
                              weightObj, tempExcits, nExcits)
        end do

        ! normal end step
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

    end subroutine calcSingleOverlapRaising

    subroutine calcSingleOverlapMixed(ilut, csf_i, excitInfo, excitations, nExcits, &
                                      posSwitches, negSwitches)
        ! control routine to calculate the single overlap excitation with
        ! mixed generators
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(inout) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        character(*), parameter :: this_routine = "calcSingleOverlapMixed"

        type(WeightObj_t) :: weights
        integer :: iOrb, iEx, deltaB
        integer(n_int) :: t(0:nifguga)
        integer(n_int), allocatable :: tempExcits(:, :)
        real(dp) :: tempWeight

        ! todo: create weight objects, here are the normal single excitation
        ! weight operators.
        weights = init_singleWeight(csf_i, excitInfo%fullEnd)

        excitInfo%currentGen = excitInfo%firstGen
        ! create according start
        call createSingleStart(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                               weights, tempExcits, nExcits)

        ! set current gen
        ! do normal updates up to the overlap site
        do iOrb = excitInfo%fullStart + 1, excitInfo%secondStart - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! at single overlap site depending on which generator ends or start
        iOrb = excitInfo%secondStart

        if (excitInfo%firstGen == gen_type%L) then
            ! lowering ends
            ASSERT(isZero(ilut, iOrb))

            do iEx = 1, nExcits
                t = tempExcits(:, iEx)
                deltaB = getDeltaB(t)

                ! have to get double excitation matrix elements in here..
                call getDoubleMatrixElement(3, 0, deltaB, excitInfo%firstGen, &
                                            excitInfo%lastGen, csf_i%B_real(iOrb), 1.0_dp, tempWeight)

                call update_matrix_element(t, tempWeight, 1)

                ! change 0 -> 3
                set_orb(t, 2 * iOrb)
                set_orb(t, 2 * iOrb - 1)

                ! store
                tempExcits(:, iEx) = t

            end do

        else
            ! raising end
            ASSERT(isThree(ilut, iOrb))

            do iEx = 1, nExcits
                t = tempExcits(:, iEx)
                deltaB = getDeltaB(t)

                call getDoubleMatrixElement(0, 3, deltaB, excitInfo%firstGen, &
                                            excitInfo%lastGen, csf_i%B_real(iOrb), 1.0_dp, tempWeight)

                call update_matrix_element(t, tempWeight, 1)

                ! change 3 -> 0
                clr_orb(t, 2 * iOrb)
                clr_orb(t, 2 * iOrb - 1)

                ! store
                tempExcits(:, iEx) = t
            end do
        end if

        ! change excitInfo and continue singleUpdates
        excitInfo%currentGen = excitInfo%lastGen

        do iOrb = excitInfo%secondStart + 1, excitInfo%fullEnd - 1
            call singleUpdate(ilut, csf_i, iOrb, excitInfo, posSwitches, negSwitches, &
                              weights, tempExcits, nExcits)
        end do

        ! lets see if that calling works
        call singleEnd(ilut, csf_i, excitInfo, tempExcits, nExcits, excitations)

    end subroutine calcSingleOverlapMixed

    subroutine calcNonOverlapDouble(ilut, csf_i, excitInfo, excitations, nExcits, &
                                    posSwitches, negSwitches)
        ! specific subroutine to calculate the non overlap double excitations
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer(n_int), allocatable :: tempExcits(:, :), tempExcits2(:, :), tmp_excitations(:, :)
        integer :: tmpNum, iEx, tmpNum2, nOpen1, nOpen2, ierr, iEx2, nMax, i, j
        type(ExcitationInformation_t) :: tmpInfo

        ! plan is to first calculate all lower excitation, and for each state
        ! there calculate the top excitation
        ! have to check that matrix element does not get reset!

        ! modifiy excitInfo so it specifies a single excitation would be better,
        ! or just be lazy for now and call the non excitInfo one.. todo
        tmpInfo = excitInfo
        tmpInfo%fullEnd = excitInfo%firstEnd
        tmpInfo%currentGen = excitInfo%firstGen

        call calcAllExcitations(ilut, csf_i, tmpInfo, posSwitches, negSwitches, &
                                .false., tempExcits, tmpNum)

        ! then change tmpInfo to deal with second excitaion properly
        tmpInfo%fullStart = excitInfo%secondStart
        tmpInfo%fullEnd = excitInfo%fullEnd
        tmpInfo%currentGen = excitInfo%lastGen

        nExcits = 0
        ! have to allocate excitation list to worst casce
        nOpen1 = count_open_orbs_ij(csf_i, excitInfo%fullStart, excitInfo%firstEnd)
        nOpen2 = count_open_orbs_ij(csf_i, excitInfo%secondStart, excitInfo%fullEnd)

        nMax = 4 + 2**(nOpen1 + nOpen2 + 2)
        allocate(tmp_excitations(0:nifguga, nMax), stat=ierr)

        ! and for all created excitations i have to calc. all possible tops
        do iEx = 1, tmpNum
            call calcAllExcitations(tempExcits(:, iEx), csf_i, tmpInfo, posSwitches, &
                                    negSwitches, .false., tempExcits2, tmpNum2)

            ! have to reencode matrix element of tempexcits(:,iEx) as it is
            ! overwritten in the call allexcitation routine
            do iEx2 = 1, tmpNum2
                call update_matrix_element(tempExcits2(:, iEx2), &
                                           extract_matrix_element(tempExcits(:, iEx), 1), 1)

                ! also use the deltaB value of the finished excitations to
                ! handle IC in the NECI code correctly...
                call setDeltaB(2, tempExcits2(:, iEx2))
            end do

            ! and add it up to one list (maybe i can set sorted to true?
            ! since no repetitions in both...
            call add_guga_lists(nExcits, tmpNum2, tmp_excitations, tempExcits2)
        end do

        deallocate(tempExcits)
        deallocate(tempExcits2)

        j = 1
        do i = 1, nExcits
            if (near_zero(extract_matrix_element(tmp_excitations(:, i), 1))) cycle

            tmp_excitations(:, j) = tmp_excitations(:, i)

            j = j + 1
        end do

        nExcits = j - 1

        allocate(excitations(0:nifguga, nExcits), stat=ierr)
        excitations = tmp_excitations(:, 1:nExcits)

        deallocate(tmp_excitations)

    end subroutine calcNonOverlapDouble

    subroutine calcDoubleExcitation_withWeight(ilut, csf_i, excitInfo, excitations, &
                                               nExcits, posSwitches, negSwitches)
        ! subroutine to calculate the double excitations involving a weight
        ! generator, which is really similar so a normal single excitation
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer(n_int), intent(out), allocatable :: excitations(:, :)
        integer, intent(out) :: nExcits
        real(dp), intent(in) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer :: iEx, we
        ! just call excitations first
        ! have to change excitInfo so single excitations are calculated
        ! correctly

        ! in future
        we = excitInfo%weight
        call calcAllExcitations(ilut, csf_i, excitInfo, posSwitches, negSwitches, &
                                .false., excitations, nExcits)

        ! and then modify the matrix element if necessary
        if ((excitInfo%weight /= excitInfo%fullStart) .and. (excitInfo%weight &
                                                             /= excitInfo%fullEnd)) then
            if (isThree(ilut, we)) then
                do iEx = 1, nExcits
                    call update_matrix_element(excitations(:, iEx), 2.0_dp, 1)
                end do
            end if
        end if

        ! also use the deltaB value of the finished excitations to indicate
        ! the IC level for the remaining NECI code
        do iEx = 1, nExcits
            ! is a single excitation in this case!
            call setDeltaB(1, excitations(:, iEx))
        end do

    end subroutine calcDoubleExcitation_withWeight

    subroutine checkCompatibility(csf_i, excitInfo, flag, posSwitches, negSwitches, opt_weight)
        ! depending on the type of excitation determined in the
        ! excitationIdentifier check if the provided ilut and excitation and
        ! the probabilistic weight function allow an excitation
        ! dont do probabilistic weight for now. just check stepvector
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        logical, intent(out) :: flag
        real(dp), intent(out), optional :: posSwitches(nSpatOrbs), &
                                           negSwitches(nSpatOrbs)

        type(WeightObj_t), intent(out), optional :: opt_weight

        real(dp) :: pw, mw, zw
        integer ::  we, st, ss, fe, en, i, j, k, lO
        type(WeightObj_t) :: weights

        ! also include probabilistic weights
        call calcRemainingSwitches_excitInfo_double( &
                csf_i, excitInfo, posSwitches, negSwitches)

        ! todo include probabilistic weights too.
        ! and check all if conditions todo! definetly mistakes there!

        ! i definetly have the occ, b and stepvector info here..

        we = excitInfo%weight
        st = excitInfo%fullStart
        ss = excitInfo%secondStart
        fe = excitInfo%firstEnd
        en = excitInfo%fullEnd
        flag = .true.

        select case (excitInfo%typ)
            ! weight + raising generator:
        case (excit_type%raising)
            if (csf_i%stepvector(we) == 0 .or. csf_i%stepvector(st) == &
                3 .or. csf_i%stepvector(en) == 0 .or. &
                (we == en .and. csf_i%stepvector(en) /= 3)) then
                flag = .false.
                return
            end if

            weights = init_singleWeight(csf_i, en)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)
            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(pw) .and. near_zero(mw) &
                 .or. csf_i%stepvector(st) == 1 .and. near_zero(pw) &
                 .or. csf_i%stepvector(st) == 2 .and. near_zero(mw)) then
                flag = .false.
                return
            end if

            ! weight + lowering generator:
        case (excit_type%lowering)

            if (csf_i%stepvector(we) == 0 &
                .or. csf_i%stepvector(en) == 3 &
                .or. csf_i%stepvector(st) == 0 &
                .or. (we == st .and. csf_i%stepvector(st) /= 3)) then
                flag = .false.
                return
            end if

            weights = init_singleWeight(csf_i, en)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)
            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)

            if ((csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw)) .or. &
                (near_zero(pw + mw))) then
                flag = .false.
                return
            end if

            ! no overlap:
        case (excit_type%non_overlap)

            i = excitInfo%i
            j = excitInfo%j
            k = excitInfo%k
            lO = excitInfo%l

            if (csf_i%stepvector(i) == 3 .or. csf_i%stepvector(j) == 0 .or. &
                csf_i%stepvector(k) == 3 .or. csf_i%stepvector(lO) == 0) then
                flag = .false.
                return
            end if

            weights = init_singleWeight(csf_i, fe)

            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)
            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)

            ! first check lower range
            if (near_zero(pw + mw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! then second
            weights = init_singleWeight(csf_i, en)

            mw = weights%proc%minus(negSwitches(ss), csf_i%B_real(ss), weights%dat)
            pw = weights%proc%plus(posSwitches(ss), csf_i%B_real(ss), weights%dat)

            if (near_zero(pw + mw) .or. &
                (csf_i%stepvector(ss) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(ss) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! single overlap lowering
        case (excit_type%single_overlap_lowering)

            if (csf_i%stepvector(fe) == 0 .or. csf_i%stepvector(en) == &
                3 .or. csf_i%stepvector(st) == 0) then
                flag = .false.
                return
            end if

            ! todo: change single overlap probWeight to correctly include
            ! bvalue restrictions to the calc.!!!
            ! todo

            weights = init_singleOverlapLowering(csf_i, fe, en, negSwitches(fe), &
                                                 posSwitches(fe), csf_i%B_real(fe))

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(pw + mw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! todo: one of the two alike gens single overlap excits has
            ! additional constraints i believe

            ! single overlap raising
        case (excit_type%single_overlap_raising)

            if (csf_i%stepvector(fe) == 0 .or. csf_i%stepvector(st) == &
                3 .or. csf_i%stepvector(en) == 0) then
                flag = .false.
                return
            end if

            weights = init_singleOverlapRaising(csf_i, fe, en, negSwitches(fe), &
                                                posSwitches(fe), csf_i%B_real(fe))

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(pw + mw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! single overlap lowering into raising
        case (excit_type%single_overlap_L_to_R)

            if (csf_i%stepvector(st) == 0 .or. csf_i%stepvector(fe) /= &
                0 .or. csf_i%stepvector(en) == 0) then
                flag = .false.
                return
            end if

            weights = init_singleWeight(csf_i, en)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)
            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(pw + mw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! single overlap raising into lowering
        case (excit_type%single_overlap_R_to_L)

            if (csf_i%stepvector(fe) /= 3 .or. csf_i%stepvector(st) == &
                3 .or. csf_i%stepvector(en) == 3) then
                flag = .false.
                return
            end if

            weights = init_singleWeight(csf_i, en)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)
            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(pw + mw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! normal double two lowering
        case (excit_type%double_lowering)

            if (csf_i%stepvector(st) == 0 .or. csf_i%stepvector(fe) == &
                3 .or. csf_i%stepvector(en) == 3 .or. &
                csf_i%stepvector(ss) == 0) then
                flag = .false.
                return
            end if

            weights = init_fullDoubleWeight(csf_i, ss, fe, en, negSwitches(ss), &
                                            negSwitches(fe), posSwitches(ss), posSwitches(fe), csf_i%B_real(ss), &
                                            csf_i%B_real(fe))

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! normal double two raising
        case (excit_type%double_raising)

            if (csf_i%stepvector(en) == 0 .or. csf_i%stepvector(ss) == &
                3 .or. csf_i%stepvector(st) == 3 .or. &
                csf_i%stepvector(fe) == 0) then
                flag = .false.
                return
            end if

            weights = init_fullDoubleWeight(csf_i, ss, fe, en, negSwitches(ss), &
                                            negSwitches(fe), posSwitches(ss), posSwitches(fe), csf_i%B_real(ss), &
                                            csf_i%B_real(fe))

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! lowering into raising into lowering
        case (excit_type%double_L_to_R_to_L)

            if (csf_i%stepvector(st) == 0 .or. csf_i%stepvector(ss) == &
                3 .or. csf_i%stepvector(en) == 3 .or. &
                csf_i%stepvector(fe) == 0) then
                flag = .false.
                return
            end if

            weights = init_fullDoubleWeight(csf_i, ss, fe, en, negSwitches(ss), &
                                            negSwitches(fe), posSwitches(ss), posSwitches(fe), csf_i%B_real(ss), &
                                            csf_i%B_real(fe))

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! raising into lowering into raising
        case (excit_type%double_R_to_L_to_R)

            if (csf_i%stepvector(en) == 0 .or. csf_i%stepvector(fe) == &
                3 .or. csf_i%stepvector(st) == 3 .or. &
                csf_i%stepvector(ss) == 0) then
                flag = .false.
                return
            end if

            weights = init_fullDoubleWeight(csf_i, ss, fe, en, negSwitches(ss), &
                                            negSwitches(fe), posSwitches(ss), posSwitches(fe), csf_i%B_real(ss), &
                                            csf_i%B_real(fe))

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! lowering into raising double
        case (excit_type%double_L_to_R)

            if (csf_i%stepvector(st) == 0 .or. csf_i%stepvector(ss) == 3 &
                .or. csf_i%stepvector(en) == 0 .or. &
                csf_i%stepvector(fe) == 3) then
                flag = .false.
                return
            end if

            weights = init_fullDoubleWeight(csf_i, ss, fe, en, negSwitches(ss), &
                                            negSwitches(fe), posSwitches(ss), posSwitches(fe), csf_i%B_real(ss), &
                                            csf_i%B_real(fe))

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! raising into lowering double
        case (excit_type%double_R_to_L)

            if (csf_i%stepvector(ss) == 0 .or. csf_i%stepvector(st) == 3 &
                .or. csf_i%stepvector(en) == 3 .or. &
                csf_i%stepvector(fe) == 0) then
                flag = .false.
                return
            end if

            weights = init_fullDoubleWeight(csf_i, ss, fe, en, negSwitches(ss), &
                                            negSwitches(fe), posSwitches(ss), posSwitches(fe), csf_i%B_real(ss), &
                                            csf_i%B_real(fe))

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! full stop two lowering
        case (excit_type%fullstop_lowering)

            if (csf_i%stepvector(st) == 0 .or. csf_i%stepvector(en) /= 0 &
                .or. csf_i%stepvector(ss) == 0) then
                flag = .false.
                return
            end if

            weights = init_singleWeight(csf_i, ss)

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! full stop two raising
        case (excit_type%fullstop_raising)

            if (csf_i%stepvector(st) == 3 .or. csf_i%stepvector(en) /= 3 .or. &
                csf_i%stepvector(ss) == 3) then
                flag = .false.
                return
            end if

            weights = init_singleWeight(csf_i, ss)

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! full stop lowering into raising
        case (excit_type%fullstop_L_to_R)

            if (csf_i%stepvector(st) == 0 .or. csf_i%stepvector(ss) == 3 &
                .or. csf_i%stepvector(en) == 0) then
                flag = .false.
                return
            end if

            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                ! the approximate exchange forces a switch of the
                ! spin-couplings at open-shell orbitals of the
                ! the exchange index. this is esp. problematic for the
                ! full-stop part, where we now have to enforce
                ! a delta-b = +-2 at the end index to ensure a
                ! flip is happening there..
                ! we need new weighting functions for that and also
                ! need to check the compatibility differently as
                ! the flips at the indices are enforced

                weights = init_forced_end_semistart_weight(csf_i, ss, en, &
                                                           negSwitches(ss), posSwitches(ss), csf_i%B_real(ss))

            else

                weights = init_semiStartWeight(csf_i, ss, en, negSwitches(ss), &
                                               posSwitches(ss), csf_i%B_real(ss))

            end if

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! full stop raising into lowering
        case (excit_type%fullstop_R_to_L)

            if (csf_i%stepvector(ss) == 0 .or. csf_i%stepvector(st) == 3 &
                .or. csf_i%stepvector(en) == 0) then
                flag = .false.
                return
            end if

            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                ! todo: the logic

                weights = init_forced_end_semistart_weight(csf_i, ss, en, &
                                                           negSwitches(ss), posSwitches(ss), csf_i%B_real(ss))

            else

                weights = init_semiStartWeight(csf_i, ss, en, negSwitches(ss), &
                                               posSwitches(ss), csf_i%B_real(ss))

            end if

            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! full start two lowering
        case (excit_type%fullstart_lowering)

            if (csf_i%stepvector(st) /= 3 .or. csf_i%stepvector(fe) == 3 &
                .or. csf_i%stepvector(en) == 3) then
                flag = .false.
                return
            end if

            ! in the actual excitation generation i use the the
            ! single weights here.. and this make more sense i must
            ! admit.
            weights = init_singleWeight(csf_i, en)

            ! update! here i shouldnt use the real available switches for
            ! the double overlap region since switches are not allowed in
            ! this kind of excitation! -> just put in 0
            ! doesnt this mean i could just use the singles weight for
            ! the non-overlap region?
            ! and no.. i should check the weights of the single excitation
            ! region..
            pw = weights%proc%plus(posSwitches(fe), csf_i%B_real(fe), weights%dat)
            mw = weights%proc%minus(negSwitches(fe), csf_i%B_real(fe), weights%dat)

            ! only 0 deltab branch valid
            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(fe) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(fe) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! full start two raising
        case (excit_type%fullstart_raising)

            if (csf_i%stepvector(st) /= 0 .or. csf_i%stepvector(fe) == 0 &
                .or. csf_i%stepvector(en) == 0) then
                flag = .false.
                return
            end if

            ! i can actually use just the singles weight..

            weights = init_singleWeight(csf_i, en)

            pw = weights%proc%plus(posSwitches(fe), csf_i%B_real(fe), weights%dat)
            mw = weights%proc%minus(negSwitches(fe), csf_i%B_real(fe), weights%dat)

            ! only 0 deltab branch valid
            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(fe) == 1 .and. near_zero(pw)) .or. &
                (csf_i%stepvector(fe) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

            ! full start lowering into raising
        case (excit_type%fullStart_L_to_R)

            if (csf_i%stepvector(st) == 0 .or. csf_i%stepvector(fe) == 3 &
                .or. csf_i%stepvector(en) == 0) then
                flag = .false.
                return
            end if

            weights = init_fullStartWeight(csf_i, fe, en, negSwitches(fe), posSwitches(fe), &
                                           csf_i%B_real(fe))

            ! then it is actually not a proper double excitation..
            ! and should not be considered here, as it is already
            ! contained in the single excitations
            if (csf_i%stepvector(st) == 3) then
                ! but i need them for the exact excitation
                ! generation
                zw = weights%proc%zero(0.0_dp, 0.0_dp, csf_i%B_real(st), weights%dat)
                pw = 0.0_dp
                mw = 0.0_dp
            else

                zw = weights%proc%zero(negSwitches(st), posSwitches(st), csf_i%B_real(st), &
                                       weights%dat)

                pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
                mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)
            end if

            if (near_zero(mw + pw + zw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(zw + pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(zw + mw)) .or. &
                (csf_i%stepvector(st) == 3 .and. near_zero(zw))) then
                flag = .false.
                return
            end if

            ! full start raising into lowering
        case (excit_type%fullstart_R_to_L)

            if (csf_i%stepvector(st) == 0 .or. csf_i%stepvector(en) == 3 &
                .or. csf_i%stepvector(fe) == 0) then
                flag = .false.
                return
            end if

            weights = init_fullStartWeight(csf_i, fe, en, negSwitches(fe), posSwitches(fe), &
                                           csf_i%B_real(fe))

            ! if its a 3 start no switches in overlap region are possible
            if (csf_i%stepvector(st) == 3) then
                zw = weights%proc%zero(0.0_dp, 0.0_dp, csf_i%B_real(st), weights%dat)
                pw = 0.0_dp
                mw = 0.0_dp
            else
                zw = weights%proc%zero(negSwitches(st), posSwitches(st), csf_i%B_real(st), &
                                       weights%dat)

                pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
                mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)
            end if

            if (near_zero(pw + mw + zw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(zw + pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(zw + mw)) .or. &
                (csf_i%stepvector(st) == 3 .and. near_zero(zw))) then
                flag = .false.
                return
            end if

            ! full start into full stop alike
        case (excit_type%fullstart_stop_alike)
            i = excitInfo%i
            j = excitInfo%j

            if (csf_i%stepvector(j) /= 3 .or. csf_i%stepvector(i) /= 0) then
                flag = .false.
                return
            end if

            ! here i essentially do not need to check the weights..
            ! since no switch is possible anyway and there is only one
            ! connecting CSF..
            weights = init_doubleWeight(csf_i, en)

            zw = weights%proc%zero(negSwitches(st), posSwitches(st), csf_i%B_real(st), &
                                   weights%dat)

            ! again only zero weight counts, as no others allowed.
            if (near_zero(zw)) flag = .false.

            ! full start into full stop mixed
        case (excit_type%fullstart_stop_mixed)

            if (csf_i%stepvector(st) == 0 .or. csf_i%stepvector(en) == 0 &
                .or. csf_i%stepvector(st) == 3 .or. &
                csf_i%stepvector(en) == 3) then
                flag = .false.
                return
            end if

            if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                ! the weights also change for fully-exchange type
                weights = init_forced_end_exchange_weight(csf_i, en)

            else

                weights = init_doubleWeight(csf_i, en)
            end if

            zw = weights%proc%zero(negSwitches(st), posSwitches(st), csf_i%B_real(st), &
                                   weights%dat)
            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)

            ! if only the 0 branch is non-zero, and both + and - branch are
            ! zero, we should abort too, since this means we would produce a
            ! diagonal contribution..
            if (near_zero(mw + pw) .or. &
                (csf_i%stepvector(st) == 1 .and. near_zero(zw + pw)) .or. &
                (csf_i%stepvector(st) == 2 .and. near_zero(zw + mw)) .or. &
                (csf_i%stepvector(st) == 3 .and. near_zero(zw))) then
                flag = .false.
                return
            end if

        end select

        if (present(opt_weight)) opt_weight = weights

    end subroutine checkCompatibility

    function get_excit_level_from_excitInfo(excitInfo) result(ic)
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer :: ic

        character(*), parameter :: this_routine = "get_excit_level_from_excitInfo"

        call stop_all(this_routine, "TODO")
        unused_var(excitInfo)
        ic = 0

    end function get_excit_level_from_excitInfo

    function calc_pgen_mol_guga_single(ilutI, nI, csf_i, ilutJ, nJ, excitInfo) result(pgen)
        integer(n_int), intent(in) :: ilutI(0:niftot), ilutJ(0:niftot)
        integer, intent(in) :: nI(nel), nJ(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        real(dp) :: p_orbs, p_guga
        ! write a pgen calculator for guga-excitations finally..
        ! but only for the molecular for now!

        ! i need the orbital picking part of pgen:
        p_orbs = calc_pgen_mol_guga_single_orbs(ilutI, nI, csf_i, excitInfo)

        p_guga = calc_pgen_mol_guga_single_guga(ilutI, nI, ilutJ, nJ, excitInfo)

        pgen = p_orbs * p_guga

    end function calc_pgen_mol_guga_single

    function calc_pgen_mol_guga_single_guga(ilutI, nI, ilutJ, nJ, excitInfo) result(pgen)
        ! this function calculates the guga branching part of the
        ! generation probability!
        integer(n_int), intent(in) :: ilutI(0:niftot), ilutJ(0:niftot)
        integer, intent(in) :: nI(nel), nJ(nel)
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        character(*), parameter :: this_routine = "calc_pgen_mol_guga_single_guga"

        call stop_all(this_routine, "TODO")

        ! todo
        unused_var(ilutI)
        unused_var(nI)
        unused_var(ilutJ)
        unused_var(nJ)
        unused_var(excitInfo)

        pgen = 0.0_dp

    end function calc_pgen_mol_guga_single_guga

    function calc_pgen_mol_guga_single_orbs(ilut, nI, csf_i, excitInfo) result(pgen)
        integer(n_int), intent(in) :: ilut(0:niftot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        real(dp) :: pgen
        real(dp) :: p_elec, p_orb, cum_sum
        integer :: cc_i, nOrb, ierr, i, a
        real(dp), allocatable :: cum_arr(:)

        call get_orbs_from_excit_info(excitInfo, i, a)

        if (IsDoub(ilut, i)) then
            p_elec = 2.0 / real(nel, dp)
        else if (IsNotOcc(ilut, i)) then
            pgen = 0.0_dp
            return
        else
            p_elec = 1.0 / real(nel, dp)
        end if

        if (IsDoub(ilut, a)) then
            pgen = 0.0_dp
            return
        end if

        if (gtID(i) == gtID(a)) then
            pgen = 0.0_dp
            return
        end if

        cc_i = ClassCountInd(1, SpinOrbSymLabel(a), G1(a)%Ml)

        nOrb = OrbClassCount(cc_i)
        allocate(cum_arr(nOrb), stat=ierr)

        if (IsDoub(ilut, i)) then
            call gen_cum_list_guga_single_3(nI, csf_i, i, cc_i, cum_arr)
        else
            if (is_beta(i)) then
                call gen_cum_list_guga_single_1(nI, csf_i, i, cc_i, cum_arr)
            else
                call gen_cum_list_guga_single_2(nI, csf_i, i, cc_i, cum_arr)
            end if
        end if

        if (near_zero(cum_sum)) then
            pgen = 0.0_dp
            return
        end if

        if (a == 1) then
            p_orb = cum_arr(1) / cum_sum
        else
            p_orb = (cum_arr(a) - cum_arr(a - 1)) / cum_sum
        end if

        pgen = p_orb * p_elec

    end function calc_pgen_mol_guga_single_orbs

    subroutine get_orbs_from_excit_info(excitInfo, a, b, c, d)
        ! routine to extract the orbitals from a excitation information
        ! c and d are optional for double excitations
        type(ExcitationInformation_t), intent(in) :: excitInfo
        integer, intent(out) :: a, b
        integer, intent(out), optional :: c, d

        character(*), parameter :: this_routine = "get_orbs_from_excit_info"

        call stop_all(this_routine, "TODO")
        a = 0
        b = 0
        c = 0
        d = 0
        unused_var(excitInfo)

    end subroutine get_orbs_from_excit_info

    subroutine pickOrbs_sym_uniform_mol_single(ilut, nI, csf_i, excitInfo, pgen)
        ! new implementation to pick single orbitals, more similar to the
        ! other neci implementations
        ! with this new looping over other orbitals it will probably also
        ! be easier to include the neccesarry switch conditions!
        ! this also applies for double excitations!! -> think about that !
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "pickOrbs_sym_uniform_mol_single"

        integer :: elec, cc_i, ierr, nOrb, orb_ind, orb_i, orb_a
        real(dp), allocatable :: cum_arr(:)
        real(dp) :: cum_sum, r, elec_factor

        unused_var(ilut)

        ! first pick completely random from electrons only!
        elec = 1 + floor(genrand_real2_dSFMT() * nEl)
        ! have to adjust pgen if it is a doubly occupied orbital afterwards
        ! -> since twice the chance to pick that orbital then!

        ! pick associated "spin orbital"
        orb_i = nI(elec)

        ! get the symmetry index:
        ! since there is no spin restriction here have to consider both
        ! again
        cc_i = ClassCountInd(1, SpinOrbSymLabel(orb_i), G1(orb_i)%Ml)

        ! get the number of orbitals in this symmetry sector
        nOrb = OrbClassCount(cc_i)
        allocate(cum_arr(nOrb), stat=ierr)

        ! actually only need one of the symmetry list, but then just have to
        ! only check if one of the two spin-orbitals is empty
        ! now have to have specific picking routines depending on the
        ! stepvalue of the already picked orbital i

        ! TODO: hm -> for the singles the matrix element gets calculated
        ! exactly! -> thats NOT EASY in the guga case!
        ! at least not as easy as in the determinant case!!
        ! write an email to simon and ask ali if that makes any sense then
        ! eg. to use the one particle elements as an approximation..
        select case (csf_i%stepvector(gtID(orb_i)))
            ! der stepvalue sagt mir auch, ob es ein alpha oder beta
            ! elektron war..
        case (1)
            elec_factor = 1.0_dp
            call gen_cum_list_guga_single_1(nI, csf_i, orb_i, cc_i, cum_arr)

        case (2)
            ! to do
            elec_factor = 1.0_dp
            call gen_cum_list_guga_single_2(nI, csf_i, orb_i, cc_i, cum_arr)

        case (3)
            ! adjust pgen, the chance to pick a doubly occupied with
            ! spinorbitals is twice as high..
            elec_factor = 2.0_dp
            call gen_cum_list_guga_single_3(nI, csf_i, orb_i, cc_i, cum_arr)

        case default
            call stop_all(this_routine, "should not have picked empty orbital")

        end select

        ! assign the spatial orbital:
        orb_i = gtID(orb_i)

        ! get the orbital
        cum_sum = cum_arr(nOrb)

        if (near_zero(cum_sum)) then
            orb_a = 0
            excitInfo%valid = .false.
            return
        else
            r = genrand_real2_dSFMT() * cum_sum
            orb_ind = binary_search_first_ge(cum_arr, r)
            orb_a = sym_label_list_spat(SymLabelCounts2(1, cc_i) + orb_ind - 1)

            if (orb_ind == 1) then
                pgen = cum_arr(1) / cum_sum
            else
                pgen = (cum_arr(orb_ind) - cum_arr(orb_ind - 1)) / cum_sum
            end if
        end if

        deallocate(cum_arr)

        ASSERT(orb_a /= orb_i)

        ! assign excitInfo and calc. final pgen
        pgen = pgen * elec_factor / real(nEl, dp)

        if (orb_a < orb_i) then
            ! raising generator
            excitInfo = assign_excitInfo_values_single(gen_type%R, orb_a, orb_i, orb_a, orb_i)

        else
            ! lowering generator
            excitInfo = assign_excitInfo_values_single(gen_type%L, orb_a, orb_i, orb_i, orb_a)

        end if

    end subroutine pickOrbs_sym_uniform_mol_single

    subroutine gen_cum_list_real_hub_1(csf_i, orb_i, cum_arr)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orb_i
        real(dp), intent(out) :: cum_arr(nSpatOrbs)

        real(dp) :: cum_sum
        integer :: i, lower, upper
        ! for the case of d(i) = 1 i need to exclude d(j) = 1 when there
        ! is no switch possible in the middle

        cum_sum = 0.0_dp

        call find_switches(csf_i, orb_i, lower, upper)

        ! so from 1 to lower-1 d=1 is possible and allowed!
        ! if lower = 1 (default if no switch is found) no accessing loop anyway

        do i = 1, lower - 1
            ! this means in this regime only doubly occupied orbs get
            ! excluded (i cant be orb_i also, since lower is strictly lower
            ! or 1, which means the loop is not entered!)

            if (csf_i%stepvector(i) == 3) then
                cum_arr(i) = cum_sum
            else
                cum_sum = cum_sum + abs(GetTMatEl(2 * orb_i, 2 * i))
                cum_arr(i) = cum_sum
            end if
        end do

        ! do i loop from lower to upper now? or until orb_i - 1 ?
        ! i can loop until upper already, since d(i) = 1, which means it
        ! get excluded here anyway!

        do i = lower, upper
            ! here d = 1 and d = 3 are excluded!
            if (mod(csf_i%stepvector(i), 2) == 1) then
                cum_arr(i) = cum_sum
            else
                cum_sum = cum_sum + abs(GetTMatEl(2 * orb_i, 2 * i))
                cum_arr(i) = cum_sum
            end if
        end do

        ! from upper+1 until end everything except d = 3 is allowed again
        do i = upper + 1, nSpatOrbs
            if (csf_i%stepvector(i) == 3) then
                cum_arr(i) = cum_sum
            else
                cum_sum = cum_sum + abs(GetTMatEl(2 * orb_i, 2 * i))
                cum_arr(i) = cum_sum
            end if
        end do

    end subroutine gen_cum_list_real_hub_1

    subroutine gen_cum_list_real_hub_2(csf_i, orb_i, cum_arr)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orb_i
        real(dp), intent(out) :: cum_arr(nSpatOrbs)

        real(dp) :: cum_sum
        integer :: i, lower, upper
        ! similar to the case above except the switch restrictrion..

        call find_switches(csf_i, orb_i, lower, upper)

        cum_sum = 0.0_dp

        do i = 1, lower - 1

            if (csf_i%stepvector(i) == 3) then
                cum_arr(i) = cum_sum
            else
                cum_sum = cum_sum + abs(GetTMatEl(2 * orb_i, 2 * i))
                cum_arr(i) = cum_sum
            end if
        end do

        do i = lower, upper
            ! here d = 2 and d = 3 are excluded!
            if (csf_i%stepvector(i) > 1) then
                cum_arr(i) = cum_sum
            else
                cum_sum = cum_sum + abs(GetTMatEl(2 * orb_i, 2 * i))
                cum_arr(i) = cum_sum
            end if
        end do

        do i = upper + 1, nSpatOrbs
            if (csf_i%stepvector(i) == 3) then
                cum_arr(i) = cum_sum
            else
                cum_sum = cum_sum + abs(GetTMatEl(2 * orb_i, 2 * i))
                cum_arr(i) = cum_sum
            end if
        end do

    end subroutine gen_cum_list_real_hub_2

    subroutine gen_cum_list_real_hub_3(csf_i, orb_i, cum_arr)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orb_i
        real(dp), intent(out) :: cum_arr(nSpatOrbs)

        real(dp) :: cum_sum
        integer :: i
        ! no real restrictions here..
        ! notice change to similar molecular routines: here the input is
        ! already in spatial orbitals!

        cum_sum = 0.0_dp

        do i = 1, nSpatOrbs

            ! actually this case is also excluded already since i know it
            ! is d(i) = 3..

            ! actually the matrix element is also always the same ...
            ! i dont actually need this weighting here..

            ! actually i only need to cycle if the stepvector is 3 otherwise
            ! just give it the same pgen increment
            ! no wait a minute.. i know the matrix element.. or?
            ! i have and should check if those are connected through a
            ! hopping possibility.. otherwise it would not make much sense..
            ! but remember: only spin-parallel hops allowed: so access TMAT
            ! always with the same spin-orbital index
            if (csf_i%stepvector(i) == 3) then
                cum_arr(i) = cum_sum
            else
                cum_sum = cum_sum + abs(GetTMatEl(2 * orb_i, 2 * i))
                cum_arr(i) = cum_sum
            end if
        end do

    end subroutine gen_cum_list_real_hub_3

    subroutine pickOrbs_real_hubbard_double(ilut, nI, csf_i, excitInfo, pgen)
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "pickOrbs_real_hubbard_double"

        pgen = 0.0_dp
        unused_var(ilut); unused_var(nI); unused_var(csf_i); unused_var(excitInfo)
        call stop_all(this_routine, &
                      "should not be at double excitations in the real-space hubbard model!")
    end subroutine pickOrbs_real_hubbard_double

    subroutine pickOrbs_real_hubbard_single(ilut, nI, csf_i, excitInfo, pgen)
        ! write a specialized orbital picker for the real-space hubbard
        ! implementation, since we do not need all the symmetry stuff and
        ! we have to take the correct TMAT values for the ml-spin values
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "pickOrbs_real_hubbard_single"

        integer :: elec, orb_i, orb_a
        real(dp) :: cum_arr(nSpatOrbs), elec_factor, cum_sum, r

        unused_var(ilut)

        ! first pick electron randomly:
        elec = 1 + floor(genrand_real2_dSFMT() * nel)

        orb_i = gtID(nI(elec))

        ! still use cum arrays to enable applying guga restrictions!
        ! but all orbitals are possible now of course

        select case (csf_i%stepvector(orb_i))

        case (1)
            ! i need a switch possibility for other d = 1 values
            call gen_cum_list_real_hub_1(csf_i, orb_i, cum_arr)

            elec_factor = 1.0_dp

        case (2)
            ! i need a switch possibility for other d = 2 values
            call gen_cum_list_real_hub_2(csf_i, orb_i, cum_arr)

            elec_factor = 1.0_dp

        case (3)
            ! no restrictions actually
            call gen_cum_list_real_hub_3(csf_i, orb_i, cum_arr)

            ! but twice the chance to have picked this spatial orbital:
            elec_factor = 2.0_dp

        case default
            call stop_all(this_routine, "should not have picked empty orb!")

        end select

        cum_sum = cum_arr(nSpatOrbs)

        if (near_zero(cum_sum)) then
            orb_a = 0
            excitInfo%valid = .false.
            return
        else
            r = genrand_real2_dSFMT() * cum_sum
            orb_a = binary_search_first_ge(cum_arr, r)

            if (orb_a == 1) then
                pgen = cum_arr(1) / cum_sum
            else
                pgen = (cum_arr(orb_a) - cum_arr(orb_a - 1)) / cum_sum

            end if
        end if

        ASSERT(orb_a /= orb_i)

        pgen = pgen * elec_factor / real(nel, dp)

        if (orb_a < orb_i) then
            ! raising generator
            excitInfo = assign_excitInfo_values_single(gen_type%R, orb_a, orb_i, orb_a, orb_i)

        else
            ! lowering generator
            excitInfo = assign_excitInfo_values_single(gen_type%L, orb_a, orb_i, orb_i, orb_a)

        end if

    end subroutine pickOrbs_real_hubbard_single

    subroutine gen_cum_list_guga_single_1(nI, csf_i, orb_i, cc_i, cum_arr)
        ! specific single orbital picker if stepvector of electron (i) is 1
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orb_i, cc_i
        real(dp), intent(out) :: cum_arr(OrbClassCount(cc_i))
        character(*), parameter :: this_routine = "gen_cum_list_guga_single_1"

        integer :: nOrb, i, label_index, j, n_id(nEl), id_i, &
                   lower, upper, s_orb, spin
        real(dp) :: cum_sum, hel

        ! if d(i) = 1 -> i can only pick d(a) = 1 if there is a switch possib
        ! d(j) = 2 inbetween! -> include that in cummulative probabilities!
        ! still not quite sure how to effectively include that...
        nOrb = OrbClassCount(cc_i)
        label_index = SymLabelCounts2(1, cc_i)

        n_id = gtID(nI)
        id_i = gtID(orb_i)

        cum_sum = 0.0_dp

        ! do some precalcing here to determine, which orbitals to exclude due
        ! to GUGA restrictions?..
        call find_switches(csf_i, id_i, lower, upper)

        if (is_beta(orb_i)) then
            spin = 1
        else
            spin = 0
        end if

        do i = 1, nOrb
            s_orb = sym_label_list_spat(label_index + i - 1)

            if (s_orb == id_i) then
                cum_arr(i) = cum_sum
                cycle
            end if

            hel = 0.0_dp

            select case (csf_i%stepvector(s_orb))

                ! include here the guga restrictions...
            case (0)
                ! no restrictions if 0 since both branches allowed

                ! single particle matrix element:
                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb - spin))

                ! do the loop over all the other electrons
                ! (is this always symmetrie allowed?..)

                ! i know both occupation numbers! of start and end!
                ! -> calc. the specific U(i,i,j,i)n(i) and U(i,j,j,j)n(j)
                ! specifically
                ! in case d=0 -> no influence for WR_(s_orb) and WL^(s_orb)
                ! and it can only be one of those
                ! also n(id_i) = 1 which is also 0 -> so exclude this from
                ! loop below

                ! if depending on the type of generator
                ! either the topCont sign (if L)
                ! or botCont sign (if R)
                ! is unknown! -> so for every k < st if R
                ! or k > en if L the sign of the two-particle matrix
                ! elements is unkown! and has to be added in terms of
                ! absolute value!

                ! problem is, if i do not know, all the signs correctly i
                ! have to add up all matrix element contributions with
                ! there absolute value.. otherwise the order of summing
                ! influences the result, and thus cannot be true!
                ! "good" thing is, i do not need to consider, so many
                ! different possibilities

                if (.not. t_mixed_hubbard) then
                    do j = 1, nEl

                        ! todo: finish all contributions later for now only do
                        ! those which are the same for all
                        ! exclude initial orbital, since this case gets
                        ! contributed already outside of loop over electrons!
                        ! but only spin-orbital or spatial??
                        if (n_id(j) == id_i) cycle
                        hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))

                        ! now depending on generator and relation of j to
                        ! st and en -> i know sign or don't

                        hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                    end do
                end if

            case (1)
                ! here i have to somehow find out if there is a
                ! (2) between s_orb and id_i
                ! how to do that?
                ! could also use the loop over nEl to check if there
                ! is a switch between (i) and (a) ->  and set matrix
                ! element to 0 otherwise... -> would make effor O(n) again
                if (s_orb < lower .or. s_orb > upper) then
                    ! only allowed if possible switch

                    hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb - spin))

                    ! also need the weight contribution at start

                    if (.not. t_mixed_hubbard) then
                        hel = hel + abs(get_umat_el(id_i, s_orb, s_orb, s_orb))

                        ! do the loop over all the other electrons
                        ! (is this always symmetrie allowed?..)

                        do j = 1, nEl

                            ! todo: finish all contributions later for now only do
                            ! those which are the same for all
                            if (n_id(j) == id_i .or. n_id(j) == s_orb) cycle
                            hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))

                            hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                        end do
                    end if
                end if

            case (2)
                ! no restrictions for 2 -> 1 excitations
                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb - spin))
                ! do the loop over all the other electrons
                ! (is this always symmetrie allowed?..)

                if (.not. t_mixed_hubbard) then
                    hel = hel + abs(get_umat_el(id_i, s_orb, s_orb, s_orb))

                    do j = 1, nEl
                        ! todo: finish all contributions later for now only do
                        ! those which are the same for all
                        if (n_id(j) == id_i .or. n_id(j) == s_orb) cycle
                        hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))

                        hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                    end do
                end if

            case (3)
                ! do nothing in this case!

            case default
                ! should not be here!
                call stop_all(this_routine, "stepvalue /= {0,1,2,3}! something is wrong!")

            end select
            cum_sum = cum_sum + abs_l1(hel)
            cum_arr(i) = cum_sum

        end do

    end subroutine gen_cum_list_guga_single_1

    subroutine gen_cum_list_guga_single_2(nI, csf_i, orb_i, cc_i, cum_arr)
        ! specific single orbital picker if stepvector of electron (i) is 2
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orb_i, cc_i
        real(dp), intent(out) :: cum_arr(OrbClassCount(cc_i))
        character(*), parameter :: this_routine = "gen_cum_list_guga_single_2"

        integer :: nOrb, i, label_index, j, n_id(nEl), id_i, &
                   lower, upper, s_orb, spin
        real(dp) :: cum_sum, hel

        ! if d(i) = 2 -> i can only pick d(a) = 2 if there is a switch possib
        ! d(j) = 1 inbetween! -> include that in cummulative probabilities!
        ! still not quite sure how to effectively include that...
        nOrb = OrbClassCount(cc_i)
        label_index = SymLabelCounts2(1, cc_i)
        n_id = gtID(nI)
        id_i = gtID(orb_i)

        cum_sum = 0.0_dp

        ! do some precalcing here to determine, which orbitals to exclude due
        ! to GUGA restrictions?..
        call find_switches(csf_i, id_i, lower, upper)

        if (is_beta(orb_i)) then
            spin = 1
        else
            spin = 0
        end if

        do i = 1, nOrb
            s_orb = sym_label_list_spat(label_index + i - 1)

            if (s_orb == id_i) then
                cum_arr(i) = cum_sum
                cycle
            end if

            hel = 0.0_dp

            select case (csf_i%stepvector(s_orb))

                ! include here the guga restrictions...
            case (0)
                ! no restrictions if 0 since both branches allowed

                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb - spin))
                ! do the loop over all the other electrons
                ! (is this always symmetrie allowed?..)

                if (.not. t_mixed_hubbard) then
                    do j = 1, nEl
                        ! todo: finish all contributions later for now only do
                        ! those which are the same for all
                        if (n_id(j) == id_i) cycle
                        hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                        hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                    end do
                end if

            case (1)
                ! no restrictions for 1 -> 2 excitations

                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb - spin))

                ! do the loop over all the other electrons
                ! (is this always symmetrie allowed?..)
                if (.not. t_mixed_hubbard) then
                    hel = hel + abs(get_umat_el(id_i, s_orb, s_orb, s_orb))

                    do j = 1, nEl

                        ! todo: finish all contributions later for now only do
                        ! those which are the same for all
                        if (n_id(j) == id_i .or. n_id(j) == s_orb) cycle
                        hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                        hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                    end do
                end if

            case (2)
                ! here i have to somehow find out if there is a
                ! (1) between s_orb and id_i
                ! how to do that?
                ! could also use the loop over nEl to check if there
                ! is a switch between (i) and (a) ->  and set matrix
                ! element to 0 otherwise... -> would make effor O(n) again
                if (s_orb < lower .or. s_orb > upper) then
                    ! only allowed if possible switch

                    hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb - spin))
                    ! do the loop over all the other electrons
                    ! (is this always symmetrie allowed?..)
                    if (.not. t_mixed_hubbard) then
                        hel = hel + abs(get_umat_el(id_i, s_orb, s_orb, s_orb))

                        do j = 1, nEl

                            ! todo: finish all contributions later for now only do
                            ! those which are the same for all
                            if (n_id(j) == id_i .or. n_id(j) == s_orb) cycle
                            hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                            hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                        end do
                    end if
                end if

            case (3)
                ! do nothing in this case!

            case default
                ! should not be here!
                call stop_all(this_routine, "stepvalue /= {0,1,2,3}! something is wrong!")

            end select
            cum_sum = cum_sum + abs_l1(hel)
            cum_arr(i) = cum_sum

        end do

    end subroutine gen_cum_list_guga_single_2

    subroutine gen_cum_list_guga_single_3(nI, csf_i, orb_i, cc_i, cum_arr)
        ! specific single orbital picker if stepvector of electron (i) is 3
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orb_i, cc_i
        real(dp), intent(out) :: cum_arr(OrbClassCount(cc_i))
        character(*), parameter :: this_routine = "gen_cum_list_guga_single_3"

        integer :: nOrb, i, label_index, j, n_id(nEl), id_i, s_orb, spin
        real(dp) :: cum_sum, hel
        ! in the case of a 3 there are actually no additional, restrictions

        nOrb = OrbClassCount(cc_i)
        label_index = SymLabelCounts2(1, cc_i)

        n_id = gtID(nI)
        id_i = gtID(orb_i)

        cum_sum = 0.0_dp

        if (is_beta(orb_i)) then
            spin = 1
        else
            spin = 0
        end if

        do i = 1, nOrb
            ! ich glaub es ist besser ueber spatial orbitals zu loopen,
            ! sonst ist das so doof mit dem ignorieren der spin symmetrie..
            s_orb = sym_label_list_spat(label_index + i - 1)

            if (s_orb == id_i) then
                cum_arr(i) = cum_sum
                cycle
            end if

            hel = 0.0_dp
            ! here i can exclude it just with check if its not occupied,
            ! since if it a 1 or 2 i will at some point
            ! new spat. orb. impl. : check if not three
            ! or do a select case on stepvector?
            select case (csf_i%stepvector(s_orb))

                ! with the predetermination of the stepvalue at (a) and since
                ! this routine only is called in the case of d(i) = 3 it makes
                ! it much easier to determine the sign of the two-particle
                ! contribution to the single excitation matrix element!
            case (0)

                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb - spin))

                ! here a contribution from orbital id_i

                if (.not. t_mixed_hubbard) then
                    hel = hel + abs(get_umat_el(id_i, id_i, s_orb, id_i))

                    ! do the loop over all the other electrons
                    ! (is this always symmetrie allowed?..)
                    do j = 1, nEl

                        ! todo: finish all contributions later for now only do
                        ! those which are the same for all
                        ! have to exclude both electrons at spatial orb i
                        if (n_id(j) == id_i) cycle
                        hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                        hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                    end do
                end if

            case (1)

                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb - spin))

                if (.not. t_mixed_hubbard) then
                    ! now contribution for both start and end
                    hel = hel + abs(get_umat_el(id_i, id_i, s_orb, id_i))
                    hel = hel + abs(get_umat_el(id_i, s_orb, s_orb, s_orb))

                    ! do the loop over all the other electrons
                    ! (is this always symmetrie allowed?..)
                    do j = 1, nEl

                        ! todo: finish all contributions later for now only do
                        ! those which are the same for all
                        if (n_id(j) == id_i .or. n_id(j) == s_orb) cycle
                        hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                        hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                    end do
                end if

            case (2)

                hel = hel + abs(GetTMatEl(orb_i, 2 * s_orb - spin))
                ! do the loop over all the other electrons
                ! (is this always symmetrie allowed?..)

                if (.not. t_mixed_hubbard) then
                    hel = hel + abs(get_umat_el(id_i, id_i, s_orb, id_i))
                    hel = hel + abs(get_umat_el(id_i, s_orb, s_orb, s_orb))

                    do j = 1, nEl

                        ! todo: finish all contributions later for now only do
                        ! those which are the same for all
                        if (n_id(j) == id_i .or. n_id(j) == s_orb) cycle

                        hel = hel + abs(get_umat_el(id_i, n_id(j), s_orb, n_id(j)))
                        hel = hel + abs(get_umat_el(id_i, n_id(j), n_id(j), s_orb))

                    end do
                end if

            case (3)
                ! do nothing actually

            case default
                ! error happend
                call stop_all(this_routine, "stepvalue /= {0,1,2,3}! something's wrong!")

            end select

            cum_sum = cum_sum + abs_l1(hel)
            cum_arr(i) = cum_sum

        end do

    end subroutine gen_cum_list_guga_single_3

    subroutine pickOrbs_sym_uniform_ueg_double(ilut, nI, csf_i, excitInfo, pgen)
        ! specific orbital picker for hubbard and UEG type models with
        ! full k-point symmetry
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen

        integer :: occ_orbs(2), ind, eleci, elecj, orb, orb_2, orb_arr(nSpatOrbs)
        real(dp) :: cum_sum, cum_arr(nSpatOrbs), pelec, r, cpt1, cpt2
        type(ExcitationInformation_t) :: excit_arr(nBasis)

        ! pick 2 electrons uniformly first:

        ! or since other symmetries are usually not used in the hubbard just:
        ! Pick a pair of electrons (i,j) to generate from.
        ! This uses a triangular mapping to pick them uniformly.
        ind = 1 + int(ElecPairs * genrand_real2_dSFMT())
        elecj = ceiling((1 + sqrt(1 + 8 * real(ind, dp))) / 2)
        eleci = ind - ((elecj - 1) * (elecj - 2)) / 2
        pelec = 1.0_dp / real(ElecPairs, dp)

        ! i pick spatial orbitals out of occupied spin-orbitals ->
        ! so the chance to pick a doubly occupied spatial orbital is twice
        ! as high -> adjust the corresponding probabilities

        ! Obtain the orbitals and their momentum vectors for the given elecs.
        occ_orbs(1) = nI(eleci)
        occ_orbs(2) = nI(elecj)

        call gen_ab_cum_list_ueg(ilut, csf_i, occ_orbs, cum_arr, excit_arr, orb_arr)

        ! then pick a orbital randomly and consider a <> b contribution
        cum_sum = cum_arr(nSpatOrbs)

        if (near_zero(cum_sum)) then
            excitInfo%valid = .false.
            return
        end if

        r = genrand_real2_dSFMT() * cum_sum
        orb = binary_search_first_ge(cum_arr, r)

        ! pick the info
        excitInfo = excit_arr(orb)

        ! and pgens, for that i need second orbital too..
        orb_2 = orb_arr(orb)
        if (orb == 1) then
            cpt1 = cum_arr(1)
        else
            cpt1 = cum_arr(orb) - cum_arr(orb - 1)
        end if

        cpt2 = 0.0_dp
        if (orb /= orb_2) then
            if (orb_2 == 1) then
                cpt2 = cum_arr(1)
            else
                cpt2 = cum_arr(orb_2) - cum_arr(orb_2 - 1)
            end if
        end if
        pgen = pelec * (cpt1 + cpt2) / cum_sum

        if (.not. is_in_pair(occ_orbs(1), occ_orbs(2))) then
            if (csf_i%stepvector(gtID(occ_orbs(1))) == 3) pgen = pgen * 2.0_dp
            if (csf_i%stepvector(gtID(occ_orbs(2))) == 3) pgen = pgen * 2.0_dp
        end if

    end subroutine pickOrbs_sym_uniform_ueg_double

    subroutine pickOrbs_sym_uniform_ueg_single(ilut, nI, csf_i, excitInfo, pgen)
        ! dummy function to abort calculation if single excitation in
        ! hubbard/ueg models gets called incorrectly
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "pickOrbs_sym_uniform_ueg_single"

        pgen = 0.0_dp
        unused_var(ilut); unused_var(nI); unused_var(csf_i); unused_var(excitInfo)

        ! single excitations shouldnt be called in hubbard/ueg simulations
        ! due to k-point symmetry
        call stop_all(this_routine, &
                      "single excitation should not be called in Hubbard/UEG models due to k-point symmetries! abort!")

    end subroutine pickOrbs_sym_uniform_ueg_single

    subroutine gen_ab_cum_list_ueg(ilut, csf_i, occ_orbs, cum_arr, excit_arr, orb_arr)
        ! create the cummulative probability array for (ab) orbital pairs
        ! in the hubbard/UEG case with k-point symmetry
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cum_arr(nSpatOrbs)
        integer, intent(out) :: orb_arr(nSpatOrbs)
        type(ExcitationInformation_t), intent(out) :: excit_arr(nBasis)

        ! determine the GUGA restrictions:
        if (is_in_pair(occ_orbs(1), occ_orbs(2))) then
            call gen_ab_cum_list_3(csf_i, occ_orbs, cum_arr, excit_arr, orb_arr)

        else
            ! determine the different types
            ! actually, due to ki + kj = ka + kb k-point symmetry the
            ! 3 3, 3 1, and 1 3 cases are actually the same!
            ! and no pesky full-stop mixed or full-start mixed are allowed in
            ! the hubbard model!
            ! not sure about this assumption above anymore! have to check that!
            if ((.not. IsDoub(ilut, occ_orbs(1))) .and. (.not. IsDoub(ilut, occ_orbs(2)))) then
                call gen_ab_cum_list_1_1(csf_i, occ_orbs, cum_arr, excit_arr, orb_arr)

            else
                call gen_ab_cum_list_3_3(csf_i, occ_orbs, cum_arr, excit_arr, orb_arr)

            end if
        end if

    end subroutine gen_ab_cum_list_ueg

    subroutine gen_ab_cum_list_1_1(csf_i, occ_orbs, cum_arr, excit_arr, orb_arr)
        ! specific routine when the occupaton of the already picked orbitals
        ! is 2 -> still think if i really want to outpout a cum_arr of lenght
        ! nBasis -> nSpatOrbs would be better and doable... !! todo
        ! only difference to 3_3 case is, that (a) can be (i,j), but which also
        ! forces (b) to be the other of (i,j), but there is the additional
        ! contraint then, that there has to be a possible switch between them!
        ! so, probably a good idea to check if there is a possible switch
        ! between i and j first, and only then allow a = (i,j)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cum_arr(nSpatOrbs)
        type(ExcitationInformation_t), intent(out), optional :: excit_arr(nSpatOrbs)
        integer, intent(out), optional :: orb_arr(nSpatOrbs)
        character(*), parameter :: this_routine = "gen_ab_cum_list_1_1"

        integer :: ki(3), kj(3), n_id(2), orb_a, ka(3), kb(3), orb_b, st, en
        integer :: z_ind
        real(dp) :: cum_sum, contrib
        logical :: tSwitch
        type(ExcitationInformation_t) :: excitInfo

        ! for legacy compatibility
        z_ind = ubound(kPointToBasisFn, 3)

        ki = G1(occ_orbs(1))%k
        kj = G1(occ_orbs(2))%k

        cum_sum = 0.0_dp
        orb_arr = 0
        cum_arr = 0.0_dp

        n_id = gtID(occ_orbs)

        tSwitch = .true.
        ! determine if there is a possible switch between i and j, to check if
        ! a = i/j should be allowed
        ! if the symmetry assumptions below(to be tested!) hold, i do not
        ! need to care about possible switches even here!
        ! meh except that a fullstart into fullstop mixed is possible..
        if (csf_i%stepvector(n_id(1)) == csf_i%stepvector(n_id(2))) then
            if (csf_i%stepvector(n_id(1)) == 1) then
                if (count_alpha_orbs_ij(csf_i, n_id(1), n_id(2)) == 0) tSwitch = .false.
            else if (csf_i%stepvector(n_id(1)) == 2) then
                if (count_beta_orbs_ij(csf_i, n_id(1), n_id(2)) == 0) tSwitch = .false.
            end if
        end if

        ! due to k-point symmetry alot of excitation types:
        ! eg. if i and j are seperate a and b also have to be!
        ! not true! ki + kj = 2*ka can defo be!
        do orb_a = 1, n_id(1) - 1

            contrib = 0.0_dp
            excitInfo%valid = .false.

            ! avoid doubly occupied orbitals
            if (csf_i%stepvector(orb_a) /= 3) then
                ka = G1(2 * orb_a)%k
                kb = ki + kj - ka

                call MomPbcSym(kb, nBasisMax)

                orb_b = gtID(KPointToBasisFn(kb(1), kb(2), z_ind, 1))

                ! just check if k-point restriction works as i thought
                ASSERT(orb_b /= n_id(1) .and. orb_b /= n_id(2))

                if (orb_a == orb_b) then
                    ! (a) must be empty and there must be a switch possible!
                    if (csf_i%stepvector(orb_a) == 0 .and. tSwitch) then
                        ! have to make contrib twice as high here since
                        ! there is no second chance to pick it the other way
                        ! around..
                        contrib = 2.0_dp
                        ! _RR_(ab) > ^RR(i) > ^R(j)
                        excitInfo = assign_excitInfo_values_double( &
                                    excit_type%fullstart_raising, &
                                    gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                    orb_a, n_id(1), orb_a, n_id(2), orb_a, orb_a, n_id(1), n_id(2), &
                                    0, 2, 1.0_dp, 1.0_dp)
                    end if
                else
                    ! check guga restrictions todo
                    if (csf_i%stepvector(orb_b) /= 3) then
                        ! no additional restrictions or?...
                        contrib = 1.0_dp

                        if (orb_b > n_id(2)) then
                            ! _R(a) > _LR(i) > ^RL(j) > ^L(b)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_R_to_L, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                        orb_a, n_id(2), orb_b, n_id(1), orb_a, n_id(1), n_id(2), orb_b, &
                                        0, 4, 1.0_dp, 1.0_dp)

                        else if (orb_b < n_id(1)) then
                            ! check if a == b
                            ! check maxima
                            st = min(orb_a, orb_b)
                            en = max(orb_a, orb_b)
                            ! _R(min) > _RR(max) > ^RR(i) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_raising, &
                                        gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        en, n_id(2), st, n_id(1), st, en, n_id(1), n_id(2), &
                                        0, 4, 1.0_dp, 1.0_dp)
                        else
                            ! b is between i and j
                            ! _R(a) > _LR(i) > ^LR(b) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_R_to_L_to_R, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        orb_a, n_id(2), orb_b, n_id(1), orb_a, n_id(1), orb_b, n_id(2), &
                                        0, 4, 1.0_dp, 1.0_dp)
                        end if
                    end if
                end if
            end if
            cum_sum = cum_sum + contrib
            cum_arr(orb_a) = cum_sum
            excit_arr(orb_a) = excitInfo
            orb_arr(orb_a) = orb_b
        end do

        ! check if orbital (i) is valid to choose from
        if (tSwitch) then
            ! do not have to deal specifically with orbital 1 ...or?
            ! do i need to consider kb restriction... yes i do think so..
            kb = kj
            ! and also need no mapping back here?
            ! hm..

            call MomPbcSym(kb, nBasisMax)
            orb_b = gtID(KPointToBasisFn(kb(1), kb(2), z_ind, 1))

            ! b has to be (j)!
            ASSERT(orb_b == n_id(2))

            contrib = 1.0_dp

            excitInfo = assign_excitInfo_values_double( &
                        excit_type%fullstart_stop_mixed, &
                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                        n_id(1), n_id(2), n_id(2), n_id(1), n_id(1), n_id(1), n_id(2), n_id(2), &
                        0, 2, 1.0_dp, 1.0_dp)

            cum_sum = cum_sum + contrib
            cum_arr(n_id(1)) = cum_sum
            excit_arr(n_id(1)) = excitInfo
            orb_arr(n_id(1)) = n_id(2)

        else
            ! otherwise orb (i) off-limits
            if (n_id(1) == 1) then
                cum_arr(1) = 0.0_dp
            else
                cum_arr(n_id(1)) = cum_arr(n_id(1) - 1)
            end if
        end if

        do orb_a = n_id(1) + 1, n_id(2) - 1

            contrib = 0.0_dp
            excitInfo%valid = .false.

            ! avoid doubly occupied orbitals
            if (csf_i%stepvector(orb_a) /= 3) then
                ka = G1(2 * orb_a)%k
                kb = ki + kj - ka

                call MomPbcSym(kb, nBasisMax)

                orb_b = gtID(KPointToBasisFn(kb(1), kb(2), z_ind, 1))

                ASSERT(orb_b /= n_id(1) .and. orb_b /= n_id(2))

                if (orb_a == orb_b) then
                    if (csf_i%stepvector(orb_a) == 0 .and. tSwitch) then
                        contrib = 2.0_dp
                        ! _L(i) > ^LR_(ab) > ^R(j)
                        excitInfo = assign_excitInfo_values_double( &
                                    excit_type%single_overlap_L_to_R, &
                                    gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                    orb_a, n_id(1), orb_a, n_id(2), n_id(1), orb_a, orb_a, n_id(2), &
                                    0, 2, 1.0_dp, 1.0_dp, 1)
                    end if
                else
                    ! check guga restrictions todo
                    if (csf_i%stepvector(orb_b) /= 3) then
                        ! no additional restrictions or?...
                        contrib = 1.0_dp

                        if (orb_b < n_id(1)) then
                            ! _R(b) > _LR(i) > ^LR(a) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_R_to_L_to_R, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        orb_b, n_id(2), orb_a, n_id(1), orb_b, n_id(1), orb_a, n_id(2), &
                                        0, 4, 1.0_dp, 1.0_dp)

                        else if (orb_b > n_id(2)) then
                            ! _L(i) > _RL(a) > ^RL(j) > ^L(b)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_L_to_R_to_L, &
                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                        orb_b, n_id(1), orb_a, n_id(2), n_id(1), orb_a, n_id(2), orb_b, &
                                        0, 4, 1.0_dp, 1.0_dp)

                        else
                            st = min(orb_a, orb_b)
                            en = max(orb_a, orb_b)
                            ! _L(i) > _RL(min) > ^LR(max) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_L_to_R, &
                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                        en, n_id(1), st, n_id(2), n_id(1), st, en, n_id(2), &
                                        0, 4, 1.0_dp, 1.0_dp)

                        end if
                    end if
                end if
            end if
            cum_sum = cum_sum + contrib
            cum_arr(orb_a) = cum_sum
            excit_arr(orb_a) = excitInfo
            orb_arr(orb_a) = orb_b
        end do

        ! fill in orbital j, depending if a switch is possible
        ! check if orbital (i) is valid to choose from
        if (tSwitch) then
            !todo: think about, that this is actually the same as the
            ! RL > RL considered above, just picked in different order..
            ! do not have to deal specifically with orbital 1 ...or?
            ! do i need to consider kb restriction... yes i do think so..
            kb = ki
            call MomPbcSym(kb, nBasisMax)
            orb_b = gtID(KPointToBasisFn(kb(1), kb(2), z_ind, 1))

            ASSERT(orb_b == n_id(1))
            ! b has to be (j)!
            contrib = 1.0_dp

            excitInfo = assign_excitInfo_values_double( &
                        excit_type%fullstart_stop_mixed, &
                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                        n_id(1), n_id(2), n_id(2), n_id(1), n_id(1), n_id(1), n_id(2), n_id(2), &
                        0, 2, 1.0_dp, 1.0_dp)

            cum_sum = cum_sum + contrib
            cum_arr(n_id(2)) = cum_sum
            excit_arr(n_id(2)) = excitInfo
            orb_arr(n_id(2)) = n_id(1)

        else

            cum_arr(n_id(2)) = cum_arr(n_id(2) - 1)
            excitInfo%valid = .false.
            excit_arr(n_id(2)) = excitInfo

        end if

        do orb_a = n_id(2) + 1, nSpatOrbs

            contrib = 0.0_dp
            excitInfo%valid = .false.
            if (csf_i%stepvector(orb_a) /= 3) then
                ka = G1(2 * orb_a)%k
                kb = ki + kj - ka

                call MomPbcSym(kb, nBasisMax)

                ! there is no "spin" restriction in the guga case
                ! so both possible b orbs have to be checked
                ! its actually NOT possible that ka = ki !! so do not have
                ! to check that case!
                orb_b = gtID(KPointToBasisFn(kb(1), kb(2), z_ind, 1))

                ! this should work as we pick beta orbital first
                ASSERT(orb_b /= n_id(1) .and. orb_b /= n_id(2))

                if (orb_a == orb_b) then
                    if (csf_i%stepvector(orb_a) == 0 .and. tSwitch) then
                        contrib = 2.0_dp
                        ! _L(i) > _LL(j) > ^LL^(ab)
                        excitInfo = assign_excitInfo_values_double( &
                                    excit_type%fullstop_lowering, &
                                    gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                    orb_a, n_id(1), orb_b, n_id(2), n_id(1), n_id(2), orb_a, orb_a, &
                                    0, 2, 1.0_dp, 1.0_dp)
                    end if
                else
                    if (csf_i%stepvector(orb_b) /= 3) then
                        ! only then its a possible excitation
                        contrib = 1.0_dp

                        ! check type of excitation
                        if (orb_b < n_id(1)) then
                            ! _R(b) > _LR(i) > ^RL(j) > ^L(a)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_R_to_L, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                        orb_b, n_id(2), orb_a, n_id(1), orb_b, n_id(1), n_id(2), orb_a, &
                                        0, 4, 1.0_dp, 1.0_dp)

                        else if (orb_b > n_id(2)) then
                            st = min(orb_a, orb_b)
                            en = max(orb_a, orb_b)
                            ! _L(i) > _LL(j) > ^LL(min) > ^L(max)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_lowering, &
                                        gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                        st, n_id(1), en, n_id(2), n_id(1), n_id(2), st, en, &
                                        0, 4, 1.0_dp, 1.0_dp)
                        else
                            ! _L(i) > _RL(b) > ^RL(j) > ^L(a)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_L_to_R_to_L, &
                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                        orb_a, n_id(1), orb_b, n_id(2), n_id(1), orb_b, n_id(2), orb_a, &
                                        0, 4, 1.0_dp, 1.0_dp)
                        end if
                    end if
                end if
            end if
            cum_sum = cum_sum + contrib
            cum_arr(orb_a) = cum_sum
            excit_arr(orb_a) = excitInfo
            orb_arr(orb_a) = orb_b
        end do

    end subroutine gen_ab_cum_list_1_1

    subroutine gen_ab_cum_list_3_3(csf_i, occ_orbs, cum_arr, excit_arr, orb_arr)
        ! specific routine when the occupaton of the already picked orbitals
        ! is 2 -> still think if i really want to outpout a cum_arr of lenght
        ! nBasis -> nSpatOrbs would be better and doable... !! todo
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cum_arr(nSpatOrbs)
        integer, intent(out) :: orb_arr(nSpatOrbs)
        type(ExcitationInformation_t), intent(out) :: excit_arr(nSpatOrbs)
        character(*), parameter :: this_routine = "gen_ab_cum_list_3_3"

        integer :: ki(3), kj(3), n_id(2), orb_a, ka(3), kb(3), orb_b, st, en
        integer :: z_ind

        real(dp) :: cum_sum, contrib
        type(ExcitationInformation_t) :: excitInfo

        ! for legacy compatibility
        z_ind = ubound(kPointToBasisFn, 3)

        ki = G1(occ_orbs(1))%k
        kj = G1(occ_orbs(2))%k

        cum_sum = 0.0_dp
        orb_arr = 0

        n_id = gtID(occ_orbs)

        ! due to k-point symmetry alot of excitation types:
        ! eg. if i and j are seperate a and b also have to be!
        ! not true! ki + kj = 2*ka can defo be!
        do orb_a = 1, n_id(1) - 1

            contrib = 0.0_dp
            excitInfo%valid = .false.

            ! avoid doubly occupied orbitals
            if (csf_i%stepvector(orb_a) /= 3) then
                ka = G1(2 * orb_a)%k
                kb = ki + kj - ka

                call MomPbcSym(kb, nBasisMax)

                ! is b allowed by the size of space? -> TODO: what does that mean?

                orb_b = gtID(KPointToBasisFn(kb(1), kb(2), z_ind, 1))

                ! just check if k-point restriction works as i thought
                ! i think this below still holds..
                ASSERT(orb_b /= n_id(1) .and. orb_b /= n_id(2))

                ! just to be safe also check if (a = b)
                if (orb_a == orb_b) then
                    if (csf_i%stepvector(orb_a) == 0) then
                        contrib = 2.0_dp
                        ! _RR_(ab) > ^RR(i) > ^R(j)
                        excitInfo = assign_excitInfo_values_double( &
                                    excit_type%fullstart_raising, &
                                    gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                    orb_a, n_id(1), orb_a, n_id(2), orb_a, orb_a, n_id(1), n_id(2), &
                                    0, 2, 1.0_dp, 1.0_dp)

                    end if
                else
                    ! check guga restrictions todo
                    if (csf_i%stepvector(orb_b) /= 3) then
                        ! no additional restrictions or?...
                        contrib = 1.0_dp

                        if (orb_b > n_id(2)) then
                            ! _R(a) > _LR(i) > ^RL(j) > ^L(b)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_R_to_L, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                        orb_a, n_id(2), orb_b, n_id(1), orb_a, n_id(1), n_id(2), orb_b, &
                                        0, 4, 1.0_dp, 1.0_dp)

                        else if (orb_b < n_id(1)) then
                            ! check if a == b
                            ! check maxima
                            st = min(orb_a, orb_b)
                            en = max(orb_a, orb_b)
                            ! _R(min) > _RR(max) > ^RR(i) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_raising, &
                                        gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        en, n_id(2), st, n_id(1), st, en, n_id(1), n_id(2), &
                                        0, 4, 1.0_dp, 1.0_dp)
                        else
                            ! b is between i and j
                            ! _R(a) > _LR(i) > ^LR(b) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_R_to_L_to_R, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        orb_a, n_id(2), orb_b, n_id(1), orb_a, n_id(1), orb_b, n_id(2), &
                                        0, 4, 1.0_dp, 1.0_dp)
                        end if
                    end if
                end if
            end if
            cum_sum = cum_sum + contrib
            cum_arr(orb_a) = cum_sum
            excit_arr(orb_a) = excitInfo
            orb_arr(orb_a) = orb_b
        end do
        if (n_id(1) == 1) then
            cum_arr(1) = 0.0_dp
        else
            cum_arr(n_id(1)) = cum_arr(n_id(1) - 1)
        end if

        do orb_a = n_id(1) + 1, n_id(2) - 1

            contrib = 0.0_dp
            excitInfo%valid = .false.

            ! avoid doubly occupied orbitals
            if (csf_i%stepvector(orb_a) /= 3) then
                ka = G1(2 * orb_a)%k
                kb = ki + kj - ka

                call MomPbcSym(kb, nBasisMax)

                ! is b allowed by the size of space? -> TODO: what does that mean?
                orb_b = gtID(KPointToBasisFn(kb(1), kb(2), z_ind, 1))

                ASSERT(orb_b /= n_id(1) .and. orb_b /= n_id(2))

                if (orb_a == orb_b) then
                    if (csf_i%stepvector(orb_a) == 0) then
                        contrib = 2.0_dp
                        ! _L(i) > ^LR_(ab) > ^R(j)
                        excitInfo = assign_excitInfo_values_double( &
                                    excit_type%single_overlap_L_to_R, &
                                    gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                    orb_a, n_id(1), orb_a, n_id(2), n_id(1), orb_a, orb_a, n_id(2), &
                                    0, 2, 1.0_dp, 1.0_dp, 1)

                    end if
                else
                    ! check guga restrictions todo
                    if (csf_i%stepvector(orb_b) /= 3) then
                        ! no additional restrictions or?...
                        contrib = 1.0_dp

                        if (orb_b < n_id(1)) then
                            ! _R(b) > _LR(i) > ^LR(a) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_R_to_L_to_R, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        orb_b, n_id(2), orb_a, n_id(1), orb_b, n_id(1), orb_a, n_id(2), &
                                        0, 4, 1.0_dp, 1.0_dp)

                        else if (orb_b > n_id(2)) then
                            ! _L(i) > _RL(a) > ^RL(j) > ^L(b)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_L_to_R_to_L, &
                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                        orb_b, n_id(1), orb_a, n_id(2), n_id(1), orb_a, n_id(2), orb_b, &
                                        0, 4, 1.0_dp, 1.0_dp)

                        else
                            st = min(orb_a, orb_b)
                            en = max(orb_a, orb_b)
                            ! _L(i) > _RL(min) > ^LR(max) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_L_to_R, &
                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                        en, n_id(1), st, n_id(2), n_id(1), st, en, n_id(2), &
                                        0, 4, 1.0_dp, 1.0_dp)

                        end if
                    end if
                end if
            end if
            cum_sum = cum_sum + contrib
            cum_arr(orb_a) = cum_sum
            excit_arr(orb_a) = excitInfo
            orb_arr(orb_a) = orb_b
        end do

        ! fill in orbital j
        ! not so sure about that anymore.. i think i got that wrong and
        ! this case is possible! damn
        cum_arr(n_id(2)) = cum_arr(n_id(2) - 1)

        do orb_a = n_id(2) + 1, nSpatOrbs

            contrib = 0.0_dp
            excitInfo%valid = .false.
            if (csf_i%stepvector(orb_a) /= 3) then
                ka = G1(2 * orb_a)%k
                kb = ki + kj - ka

                call MomPbcSym(kb, nBasisMax)

                ! there is no "spin" restriction in the guga case
                ! so both possible b orbs have to be checked
                ! its actually NOT possible that ka = ki !! so do not have
                ! to check that case!
                orb_b = gtID(KPointToBasisFn(kb(1), kb(2), z_ind, 1))

                ! this should work as we pick beta orbital first
                ASSERT(orb_b /= n_id(1) .and. orb_b /= n_id(2))

                if (orb_a == orb_b) then
                    if (csf_i%stepvector(orb_a) == 0) then
                        contrib = 2.0_dp
                        ! _L(i) > _LL(j) > ^LL^(ab)
                        excitInfo = assign_excitInfo_values_double( &
                                    excit_type%fullstop_lowering, &
                                    gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                    orb_a, n_id(1), orb_b, n_id(2), n_id(1), n_id(2), orb_a, orb_a, &
                                    0, 2, 1.0_dp, 1.0_dp)
                    end if
                else
                    if (csf_i%stepvector(orb_b) /= 3) then
                        ! only then its a possible excitation
                        contrib = 1.0_dp

                        ! check type of excitation
                        if (orb_b < n_id(1)) then
                            ! _R(b) > _LR(i) > ^RL(j) > ^L(a)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_R_to_L, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                        orb_b, n_id(2), orb_a, n_id(1), orb_b, n_id(1), n_id(2), orb_a, &
                                        0, 4, 1.0_dp, 1.0_dp)

                        else if (orb_b > n_id(2)) then
                            st = min(orb_a, orb_b)
                            en = max(orb_a, orb_b)
                            ! _L(i) > _LL(j) > ^LL(min) > ^L(max)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_lowering, &
                                        gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                        st, n_id(1), en, n_id(2), n_id(1), n_id(2), st, en, &
                                        0, 4, 1.0_dp, 1.0_dp)
                        else
                            ! _L(i) > _RL(b) > ^RL(j) > ^L(a)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_L_to_R_to_L, &
                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                        orb_a, n_id(1), orb_b, n_id(2), n_id(1), orb_b, n_id(2), orb_a, &
                                        0, 4, 1.0_dp, 1.0_dp)
                        end if
                    end if
                end if
            end if
            cum_sum = cum_sum + contrib
            cum_arr(orb_a) = cum_sum
            excit_arr(orb_a) = excitInfo
            orb_arr(orb_a) = orb_b
        end do
        ! make orbital (j) invalid for excitation also... not actually
        ! necessary, as i never pick that orbital anyway,,
        excitInfo%valid = .false.
        excit_arr(n_id(2)) = excitInfo

    end subroutine gen_ab_cum_list_3_3

    subroutine gen_ab_cum_list_3(csf_i, occ_orbs, cum_arr, excit_arr, orb_arr)
        ! specific routine, when 2 already picked orbtitals are from same
        ! spatial orbital
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cum_arr(nSpatOrbs)
        integer, intent(out) :: orb_arr(nSpatOrbs)
        type(ExcitationInformation_t), intent(out) :: excit_arr(nSpatOrbs)

        integer :: ki(3), kj(3), ka(3), kb(3), a, b, en, st, i, z_ind
        real(dp) :: cum_sum, contrib
        type(ExcitationInformation_t) :: excitInfo
        logical :: tSwitch

        z_ind = ubound(kPointToBasisFn, 3)

        ki = G1(occ_orbs(1))%k
        kj = G1(occ_orbs(2))%k

        i = gtID(occ_orbs(1))
        ! redo this!
        cum_sum = 0.0_dp
        orb_arr = 0
        b = 0

        ! GUGA restrictions change depending on n(i), n(j)
        ! I == J -> n(i) = 3
        ! no restrictions on a
        ! todo!! stupid! if i = j => ki = kj => ka = kb! only have to loop
        ! once over orbitals!
        ! NO! 2ki = ka + kb -> kann trotzdem noch mehrere möglichkeiten geben!
        ! change that! (but maybe reuse it for 3 3 case
        ! aber vl. doch spatial orbitals...
        do a = 1, i - 1

            contrib = 0.0_dp
            excitInfo%valid = .false.
            tSwitch = .true.

            if (csf_i%stepvector(a) /= 3) then
                ! determine fiting b by k-point restrictions ->
                ! and only allow b > a
                ka = G1(2 * a)%k
                kb = ki + kj - ka

                call MomPbcSym(kb, nBasisMax)

                ! is b allowed by the size of space? -> TODO: what does that mean?
                ! this is specific for the UEG.. so i do not need it for the
                ! Hubbard model i am implementing right now..

                ! there is no "spin" restriction in the guga case
                ! so both possible b orbs have to be checked
                ! its actually NOT possible that ka = ki !! so do not have
                ! to check that case!
                b = gtID(KPointToBasisFn(kb(1), kb(2), z_ind, 1))

                ! this should work as we pick beta orbital first

                ! orb_a should not be equal orb_b due to k point restrictions
                ! hm... it could happen i guess.. but i have to check it
                ! since i did not include this in the possibilities
                ! below..
                ! ok so i just realised this case can happen...
                ! so i have to fix that and take all the possible
                ! combinations into account..
                if (a == b) then
                    ! then orb a has to be empty!
                    ! actually this type of excitation is not even possible
                    ! with momentum conservation.. todo!
                    if (csf_i%stepvector(a) == 0) then
                        !todo: remove this possibility then!
                        ! _RR_(ab) > ^RR^(ij)
                        excitInfo = assign_excitInfo_values_double( &
                                    excit_type%fullstart_stop_alike, &
                                    gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                    a, i, a, i, a, a, i, i, 0, 2, 1.0_dp, 1.0_dp)

                        ! use uniform, since all the integrals are equal
                        ! anyway
                        ! if (a==b) there is only one entry in the cum-list
                        ! which leads to that excitation, compared to the
                        ! usual 2 for other types of excitations..
                        ! so i could consider multiplying it by 2 just
                        ! to make the pgens more homogenous.. since matrix
                        ! elements usually will be the same i guess..
                        ! no, do it by making it twice as big later..
                        ! NO! i have to do it here or otherwise the actual
                        ! probability of picking this orbital is not twice as
                        ! high..
                        contrib = 2.0_dp

                    end if
                else
                    if (csf_i%stepvector(b) /= 3) then
                        ! only then its a possible excitation
                        ! hm... i guess i should check if there is a possible
                        ! switch in this case too.. if both a and b have the
                        ! same stepvector there needs to be a switch possible
                        ! at some place to lead to a valid excitation..
                        ! i have not yet done that in the molecular excitation
                        ! generator and i do not know why.. was there a reason?

                        st = min(a, b)
                        en = max(a, b)

                        if (csf_i%stepvector(a) == csf_i%stepvector(b)) then
                            if (csf_i%stepvector(a) == 1 .and. &
                                count_alpha_orbs_ij(csf_i, st, en) == 0) tSwitch = .false.
                            if (csf_i%stepvector(a) == 2 .and. &
                                count_beta_orbs_ij(csf_i, st, en) == 0) tSwitch = .false.
                        end if

                        if (tSwitch) then
                            contrib = 1.0_dp

                            ! then have to determine the order of orbitals
                            if (b > i) then
                                ! _R(a) > ^RL_(ij) > ^L(b)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%single_overlap_R_to_L, &
                                            gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                            b, i, a, i, a, i, i, b, &
                                            0, 2, 1.0_dp, 1.0_dp, 1)

                            else
                                ! only the extrema count
                                ! if both have the same stepvalue i need to check
                                ! if there is a switch possible i guess..
                                ! _R(min) > _RR(max) > ^RR^(ij)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%fullstop_raising, &
                                            gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                            st, i, en, i, st, en, i, i, 0, 2, 1.0_dp, 1.0_dp)
                            end if
                        end if
                    end if
                end if
            end if
            cum_sum = cum_sum + contrib
            cum_arr(a) = cum_sum
            excit_arr(a) = excitInfo
            orb_arr(a) = b
        end do
        ! zero out orbital (i)
        if (i == 1) then
            cum_arr(i) = 0.0_dp
        else
            cum_arr(i) = cum_arr(i - 1)
        end if

        do a = i + 1, nSpatOrbs

            contrib = 0.0_dp
            excitInfo%valid = .false.
            tSwitch = .true.

            if (csf_i%stepvector(a) /= 3) then
                ka = G1(2 * a)%k
                kb = ki + kj - ka

                call MomPbcSym(kb, nBasisMax)

                ! there is no "spin" restriction in the guga case
                ! so both possible b orbs have to be checked
                ! its actually NOT possible that ka = ki !! so do not have
                ! to check that case!
                b = gtID(KPointToBasisFn(kb(1), kb(2), z_ind, 1))

                ! this should work as we pick beta orbital first
                ! as i realised this possibility below is possible!
                if (a == b) then
                    !todo: i think with momentum conservation this below is not
                    ! even possible, if really not -> remove it
                    if (csf_i%stepvector(a) == 0) then
                        ! _LL_(ij) > ^LL^(ab)
                        excitInfo = assign_excitInfo_values_double( &
                                    excit_type%fullstart_stop_alike, &
                                    gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                    a, i, a, i, i, i, a, a, 0, 2, 1.0_dp, 1.0_dp)

                        ! same convention as above since only one entry
                        ! for this kind of ab combination in list
                        ! no! just generally add up the 2 combinations later
                        ! on at the end of the orbital picking!
                        contrib = 2.0_dp

                    end if
                else
                    if (csf_i%stepvector(b) /= 3) then
                        ! only then its a possible excitation

                        st = min(a, b)
                        en = max(a, b)

                        if (csf_i%stepvector(a) == csf_i%stepvector(b)) then
                            if (csf_i%stepvector(a) == 1 .and. &
                                count_alpha_orbs_ij(csf_i, st, en) == 0) tSwitch = .false.
                            if (csf_i%stepvector(a) == 2 .and. &
                                count_beta_orbs_ij(csf_i, st, en) == 0) tSwitch = .false.
                        end if

                        if (tSwitch) then
                            contrib = 1.0_dp

                            ! now we know a is bigger than (i)
                            if (b < i) then
                                ! _R(b) > ^RL_(ij) > ^L(a)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%single_overlap_R_to_L, &
                                            gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                            a, i, b, i, b, i, i, a, &
                                            0, 2, 1.0_dp, 1.0_dp, 1)

                            else
                                ! only extrema count
                                ! _LL_(ij) > ^LL(min) > ^L(max)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%fullstart_lowering, &
                                            gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                            st, i, en, i, i, i, st, en, 0, 2, 1.0_dp, 1.0_dp)
                            end if
                        end if
                    end if
                end if
            end if
            ! update the cumsum
            cum_sum = cum_sum + contrib
            cum_arr(a) = cum_sum
            excit_arr(a) = excitInfo
            orb_arr(a) = b
        end do

        excitInfo%valid = .false.
        excit_arr(i) = excitInfo

    end subroutine gen_ab_cum_list_3

    subroutine pickOrbs_sym_uniform_mol_double(ilut, nI, csf_i, excitInfo, pgen)
        ! new orbital picking routine, which is closer to simons already
        ! implemented one for the determinant version
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "pickOrbs_sym_uniform_mol_double"

        integer :: occ_orbs(2), sym_prod, cc_a, cc_b, sum_ml, ind_res, &
                   st, en, a, b, i, j
        real(dp) :: int_contrib(2), cum_sum(2), cum_arr(nSpatOrbs), &
                    int_switch(2), cum_switch(2)
        logical :: range_flag

        ! has to be in the interface for function pointers
        unused_var(ilut)

        ! pick 2 ocupied orbitals randomly:
        call pick_elec_pair_uniform_guga(nI, occ_orbs, sym_prod, sum_ml, &
                                         pgen)
        ! the 2 picked electrons are ALWAYS ordered too! so i already know
        ! that I < J

        ! also store spatial orbitals:
        i = gtID(occ_orbs(1))
        j = gtID(occ_orbs(2))

        ! the occ_orbs are given in "spin-orbital" form.. not sure yet which
        ! implementation to choose!...

        ! then pick orbital a, weighted with FCIDUMP integrals
        call pick_a_orb_guga_mol(csf_i, occ_orbs, int_contrib(1), cum_sum(1), &
                                 cum_arr, a)

        ! changed that orbital a is now a spatial orbital already!!
        ! pick_a_orb now gives me a spatial orbital again!
        ! kill excitation if a == 0
        if (a == 0) then
            excitInfo%valid = .false.
            return
        end if

        ! now the additional GUGA and symmetry restrictions kick in..

        ! TODO: check if symmetries still work with GUGA so easily
        ! here always index it with the beta_spin part by default! and figure
        ! out how to later access it to not restrict any excitations on spin!
        ! change it to 2*a only.. since from now on i am ignoring spin
        ! symmetry anyways
        cc_a = ClassCountInd(2 * a)
        ! cc_a = ClassCountInd(1, SpinOrbSymLabel(2*a), G1(2*a)%Ml)
        ! hm have to think on how to use the symmetries of the spin-orbitals
        ! in the GUGA case..
        ! since original quanities sum_ml and iSpn are not really meaning
        ! anything in the CSF approach

        ! picking i and j in determinental case gives me information if my
        ! electrons have parallel or antiparallel spin ->
        ! depending on that my classcount for b depends on iSpn and the "spin"
        ! of orbitals b -> so i have to somehow exclude this info in the
        ! orbital choosing of b, and only consider spatial symmetries
        ! the spin restriction is not valid in the guga case since i
        ! essentially pick spatial orbitals -> to take into account both
        ! spin-orbital types -> choose cc_a to always take the beta
        ! orbital -> and choose both ispn = 2 and 3 to get all the
        ! possible spin_orbitals -> and then consruct a list of
        ! fitting spatial orbitals in the orbital picker!

        ! for (symorbs) quantity should not matter which spin the picked
        ! electrons have! (still have to check that thoroughly though!)
        ! so it should be enough to pick only one cc_index for b and
        ! convert to spatial orbitals then, and pick randomly from those!
        cc_b = get_paired_cc_ind(cc_a, sym_prod, sum_ml, 2)

        ! determine the GUGA restrictions on the orbitals!
        ! + to determine those restrictions also allows partially determining
        ! the typ of excitation already!
        ! restrctions on b can be:
        ! only exclude orbital (a) -> this is always the case! -> no additional
        ! exclude spin-orbital belonging to I or J
        ! b being strictly lower/higher than a certain orbital

        ! those two are essentially the only restrictions
        ! how to control that? do optional orbital input to the b orbital
        ! picker:
        ! if no additional input -> no additional restriction
        ! if positive integer input: exclude this specific orbital
        ! if integer + logical: if integer positive : excude everything above
        !                       if negative exclude everything below!

        ! do this + the symmetry and integral contributions, then i should
        ! be finished..

        ! and write the corresponding logic out here to determine the cases
        ! and to pre-determine the type of excitation

        ! now first check if the two spin_orb indices belong to the same
        ! spatial orbitals

        ! have to do the logic to determine the type of guga restrictions
        ! and determine the type of excitation here!
        ! and this also influences the pgen calculation!

        ! set default
        ind_res = 0
        range_flag = .false.

        if (is_in_pair(occ_orbs(1), occ_orbs(2))) then
            ! this also implies that there is no x1 matrix element contribution
            ! for the generator matrix elements -> always bias towards
            ! U + U' for orbital b
            ! this implies d(i) == 3 and there is NO additional restriction
            ! on hole a, excpet symmetries and stuff..
            ! in this case:
            ! full-start lowering and
            ! full-stop raising
            ! mixed raising into lowering single overlap
            ! excitations are possible
            ! depending where a and b are compared to I/J

            ! additional there are no restrictions on picking b
            ! todo: i could move the determination of the excitation
            ! information into the pick_b orb function. since there i
            ! already have to check the relation to the other already picked
            ! orbitals.. but not now i guess..
            call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, int_contrib(2), &
                                     cum_sum(2), b)
            ! TODO: still have to decide if i output SPATIAL or SPIN orbital
            ! for picked b... so i might have to convert at some point here
            ! and below!

            if (b == 0) then
                excitInfo%valid = .false.
                return
            end if
            ! changed, that b ouput is already a spatial orb..

            ! short names for stupid reasons
            en = max(a, b)
            st = min(a, b)

            ! also have to consider is a and b are from same spatial orb or?

            if (a == b) then
                if (a > i) then
                    !_LL_(ij) > ^LL^(ab)
                    excitInfo = assign_excitInfo_values_double( &
                                excit_type%fullstart_stop_alike, &
                                gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                b, i, b, i, i, i, b, b, 0, 2, 1.0_dp, 1.0_dp)

                else
                    ! _RR_(ab) > ^RR^(ij)
                    excitInfo = assign_excitInfo_values_double( &
                                excit_type%fullstart_stop_alike, &
                                gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                a, i, a, i, a, a, i, i, 0, 2, 1.0_dp, 1.0_dp)
                end if
            else
                if (a > i .and. b > i) then
                    ! _LL_(ij) > ^LL(min(a,b)) -> ^L(max(a,b))
                    excitInfo = assign_excitInfo_values_double( &
                                excit_type%fullstart_lowering, &
                                gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                en, i, st, i, i, i, st, en, 0, 2, 1.0_dp, 1.0_dp, 2)

                else if (a < i .and. b < i) then
                    ! _R(min(a,b)) > _RR(max(a,b) > ^RR^(ij)
                    excitInfo = assign_excitInfo_values_double( &
                                excit_type%fullstop_raising, &
                                gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                st, i, en, i, st, en, i, i, 0, 2, 1.0_dp, 1.0_dp, 2)

                else
                    ! _R(min(a,b)) > ^RL_(i) > ^L(max(a,b))
                    excitInfo = assign_excitInfo_values_double( &
                                excit_type%single_overlap_R_to_L, &
                                gen_type%R, gen_type%L, gen_type%R, gen_type%R, gen_type%L, &
                                st, i, en, i, st, i, i, en, 0, 2, 1.0_dp, 1.0_dp, 1)

                end if

            end if
            ! could have picked a and b in either way around..
            ! have p(a)*p(b|a) have to determine p(b)*p(a|b)
            ! can resuse cum_arrays form both picking
            ! no.. only first cum_arr -> have to reconstruct second one..
            ! but actually dont need an array
            ! do that below at the end since its always the same
            call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, a, &
                                          int_switch(2), cum_switch(2))

            ! should not happen but assert here that the cummulative
            ! probabilty is not 0
            ASSERT(.not. near_zero(cum_switch(2)))

        else
            ! if they are not in a pair there are more possibilites
            ! i know through the triangular mapping, that the 2 picked
            ! electrons are always ordered!
            if (csf_i%stepvector(i) == 3) then
                ! if the picked stepvector is doubly occupied, since i pick
                ! based on spinorbitals but then only consider spatial orbitals
                ! there are more chances to pick those orbitals then...
                ! but thats probably a not valid bias towards those type
                ! of excitations... -> have to think of a new way to pick
                ! occupied spatial orbitals
                !TODO: yes definetly not pick based on spin-orbitals then.
                ! i unnaturally pick towards doubly occupied sites.. hm..
                pgen = 2.0_dp * pgen
                if (csf_i%stepvector(j) == 3) then

                    ! see above for description
                    pgen = 2.0_dp * pgen

                    ! no additional restrictions in picking b
                    call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                             int_contrib(2), cum_sum(2), b)

                    if (b == 0) then
                        excitInfo%valid = .false.
                        return
                    end if

                    ! the good thing -> independent of the ordering, the pgens
                    ! have the same restrictions independent of the order how
                    ! a and b are picked!
                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, a, &
                                                  int_switch(2), cum_switch(2))

                    ! now determine the type of excitation:
                    if (a == b) then

                        if (a > j) then
                            ! _L(i) -> _LL(j) > ^LL^(ab)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%fullstop_lowering, &
                                        gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                        a, i, a, j, i, j, a, a, 0, 2, 1.0_dp, 1.0_dp, 2)

                        else if (a < i) then
                            ! _RR_(ab) > ^RR(i) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%fullstart_raising, &
                                        gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        a, i, a, j, a, a, i, j, 0, 2, 1.0_dp, 1.0_dp, 2)

                        else
                            ! _L(i) > ^LR_(ab) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%single_overlap_L_to_R, &
                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                        a, i, a, j, i, a, a, j, 0, 2, 1.0_dp, 1.0_dp, 1)

                        end if
                    else
                        ! only the maximums count..
                        en = max(a, b)
                        st = min(a, b)
                        if (en > j) then
                            if (st > j) then
                                ! _L(i) > _LL(j) > ^LL(min) > ^L(max)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%double_lowering, &
                                            gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                            st, i, en, j, i, j, st, en, 0, 4, 1.0_dp, 1.0_dp)
                            else if (st < i) then
                                ! _R(min) > _LR(i) > ^RL(j) > ^L(max)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%double_R_to_L, &
                                            gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                            st, j, en, i, st, i, j, en, 0, 4, 1.0_dp, 1.0_dp)
                            else
                                ! _L(i) > _RL(min) > ^RL(j) > ^L(max)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%double_L_to_R_to_L, &
                                            gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                            en, i, st, j, i, st, j, en, 0, 4, 1.0_dp, 1.0_dp)

                            end if
                        else if (en < i) then
                            ! this implies that st < en < I !
                            ! _R(min) > _RR(max) > ^RR(i) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%double_raising, &
                                        gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        st, j, en, i, st, en, i, j, 0, 4, 1.0_dp, 1.0_dp)

                        else
                            if (st < i) then
                                ! _R(min) > _LR(i) > ^LR(max) > ^R(j)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%double_R_to_L_to_R, &
                                            gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                            st, j, en, i, st, i, en, j, 0, 4, 1.0_dp, 1.0_dp)
                            else
                                ! _L(i) > _RL(min) > ^LR(max) > ^R(j)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%double_L_to_R, &
                                            gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                            en, i, st, j, i, st, en, j, 0, 4, 1.0_dp, 1.0_dp)
                            end if
                        end if
                    end if
                else
                    ! so its a [3 1] configuration to start with
                    ! depending on the position of a there are restrictions
                    ! on b
                    ! have to check if A == J
                    if (a == j) then
                        ! then b has to be strictly lower then j!
                        call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                 int_contrib(2), cum_sum(2), b, -j, .true.)

                        if (b == 0) then
                            excitInfo%valid = .false.
                            return
                        end if

!
                        ! but for picking it the other way around there is no
                        ! restriction on the orbitals
                        ! although for the nasty mixed full-stops i have to
                        ! recalculate the pgens anyway..
                        call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                      a, int_switch(2), cum_switch(2))

                        ! determine excit
                        ! ATTENTION for the 2 below orbital J must not be
                        ! doubly occupied or else its just a single-like
                        ! excitation!! -> check if thats the case!
                        ! since a only can be a empty orbital -> it only can be
                        ! a singly occupied orbital!
                        if (b < i) then
                            ! _R(b) > _LR(i) > ^RL^(ja)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%fullstop_R_to_L, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        b, j, j, i, b, i, j, j, 0, 4, 1.0_dp, 1.0_dp)

                        else
                            ! _L(i) > _RL(b) > ^RL^(ja)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%fullstop_L_to_R, &
                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                        j, i, b, j, i, b, j, j, 0, 4, 1.0_dp, 1.0_dp)

                        end if

                    else
                        ! if its not j
                        if (a > j) then
                            ! b can not be J!
                            call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                     int_contrib(2), cum_sum(2), b, j)

                            if (b == 0) then
                                excitInfo%valid = .false.
                                return
                            end if

                            ! depending on where b is the excitation is defined
                            ! and the pgen influence from picking a <-> b
                            if (b > j) then
                                ! both are on top -> same pgen
                                ! p(b) is always determinable from cum_arr... do
                                ! it outside!
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                              a, int_switch(2), cum_switch(2), j)

                                if (a == b) then
                                    ! _L(i) > _LL(j) > ^LL^(ab)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%fullstop_lowering, &
                                                gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                                a, i, a, j, i, j, a, a, 0, 4, 1.0_dp, 1.0_dp)

                                else
                                    st = min(a, b)
                                    en = max(a, b)
                                    ! _L(i) > _LL(j) > ^LL(min) > ^L(max)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_lowering, &
                                                gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                                st, i, en, j, i, j, st, en, 0, 4, 1.0_dp, 1.0_dp)

                                end if

                            else
                                ! pgen is not restricted
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                              a, int_switch(2), cum_switch(2))
                                if (b < i) then
                                    ! _R(b) > _LR(i) > ^RL(j) > ^L(a)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_R_to_L, &
                                                gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                                b, j, a, i, b, i, j, a, 0, 4, 1.0_dp, 1.0_dp)

                                else
                                    ! _L(i) > _RL(b) > ^RL(j) > ^L(a)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_L_to_R_to_L, &
                                                gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                                a, i, b, j, i, b, j, a, 0, 4, 1.0_dp, 1.0_dp)

                                end if
                            end if
                        else
                            ! there is no restriction in picking b
                            call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                     int_contrib(2), cum_sum(2), b)

                            if (b == 0) then
                                excitInfo%valid = .false.
                                return
                            end if
                            ! check where b is
                            if (b == j) then
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                              b, a, int_switch(2), cum_switch(2), &
                                                              -j, .true.)

                                if (a < i) then
                                    ! _R(a) > _LR(i) > ^RL^(jb)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%fullstop_R_to_L, &
                                                gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                a, j, j, i, a, i, j, j, 0, 2, 1.0_dp, 1.0_dp)

                                else
                                    ! _L(i) > _RL(a) > ^RL^(jb)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%fullstop_L_to_R, &
                                                gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                                j, i, a, j, i, a, j, j, 0, 2, 1.0_dp, 1.0_dp)

                                end if

                            else
                                if (b > j) then
                                    ! a could not have been j
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2), j)
                                    if (a < i) then
                                        ! _R(a) > _LR(i) > ^RL(j) > ^L(b)
                                        excitInfo = assign_excitInfo_values_double( &
                                                    excit_type%double_R_to_L, &
                                                    gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                                    a, j, b, i, a, i, j, b, 0, 4, 1.0_dp, 1.0_dp)

                                    else
                                        ! _L(i) > _RL(a) > ^RL(j) > ^L(b)
                                        excitInfo = assign_excitInfo_values_double( &
                                                    excit_type%double_L_to_R_to_L, &
                                                    gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                                    b, i, a, j, i, a, j, b, 0, 4, 1.0_dp, 1.0_dp)

                                    end if
                                else
                                    ! no restric, on other order
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2))

                                    if (a == b) then
                                        if (a < i) then
                                            ! _RR_(ab) > ^RR(i) > ^R(j)
                                            excitInfo = assign_excitInfo_values_double( &
                                                        excit_type%fullstart_raising, &
                                                        gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                        a, i, a, j, a, a, i, j, 0, 2, 1.0_dp, 1.0_dp)
                                        else
                                            ! _L(i) > ^LR_(ab) > ^R(j)
                                            excitInfo = assign_excitInfo_values_double( &
                                                        excit_type%single_overlap_L_to_R, &
                                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                                        a, i, a, j, i, a, a, j, 0, 2, 1.0_dp, 1.0_dp, 1)
                                        end if
                                    else
                                        ! only min and max importtant
                                        st = min(a, b)
                                        en = max(a, b)

                                        if (en < i) then
                                            ! _R(min) > _RR(max) > ^RR(i) > ^R(j)
                                            excitInfo = assign_excitInfo_values_double( &
                                                        excit_type%double_raising, &
                                                        gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                        st, i, en, j, st, en, i, j, 0, 4, 1.0_dp, 1.0_dp)

                                        else
                                            if (st < i) then
                                                ! _R(min) > _LR(i) > ^LR(max) > ^R(j)
                                                excitInfo = assign_excitInfo_values_double( &
                                                            excit_type%double_R_to_L_to_R, &
                                                            gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                            st, j, en, i, st, i, en, j, 0, 4, 1.0_dp, 1.0_dp)

                                            else
                                                ! _L(i) > _RL(min) > ^LR(max) > ^R(j)
                                                excitInfo = assign_excitInfo_values_double( &
                                                            excit_type%double_L_to_R, &
                                                            gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                                            en, i, st, j, i, st, en, j, 0, 4, 1.0_dp, 1.0_dp)
                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            else
                if (csf_i%stepvector(j) == 3) then
                    pgen = 2.0_dp * pgen
                    ! its a [1 3] configuration -> very similar to [3 1]
                    if (a == i) then
                        ! b has to be strictly higher then I
                        call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                 int_contrib(2), cum_sum(2), b, i, .true.)

                        if (b == 0) then
                            excitInfo%valid = .false.
                            return
                        end if
                        ! no restriction the other way around
                        call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, a, &
                                                      int_switch(2), cum_switch(2))

                        ! ATTENTION: here there have to be the additional
                        ! constraint that spat orb I must not be doubly
                        ! occupied! or otherwise its a single-like excitation!!
                        ! did i consider that?`
                        ! yes! it cant be a doubly occupied orbital since
                        ! orbital a and b can only be non -occupied orbitals!
                        if (b > j) then
                            ! _RL_(ia) > ^RL(j) > ^L(b)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%fullstart_R_to_L, &
                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                        b, i, i, j, i, i, j, b, 0, 2, 1.0_dp, 1.0_dp)

                        else
                            ! _RL_(ia) > ^LR(b) > ^R(j)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%fullstart_L_to_R, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        b, i, i, j, i, i, b, j, 0, 2, 1.0_dp, 1.0_dp)
                        end if

                    else
                        ! if its not i
                        if (a < i) then
                            ! b cannot be I
                            call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                     int_contrib(2), cum_sum(2), b, i)

                            if (b == 0) then
                                excitInfo%valid = .false.
                                return
                            end if

                            if (b < i) then
                                ! same pgen restrictions
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                              b, a, int_switch(2), cum_switch(2), i)

                                if (a == b) then
                                    ! _RR_(ab) > ^RR(i) > ^R(j)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%fullstart_raising, &
                                                1, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                a, i, a, j, a, a, i, j, 0, 2, 1.0_dp, 1.0_dp)
                                else
                                    st = min(a, b)
                                    en = max(a, b)

                                    ! _R(min) > _RR(max) > ^RR(i) > ^R(j)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_raising, &
                                                gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                en, i, st, j, st, en, i, j, 0, 4, 1.0_dp, 1.0_dp)
                                end if
                            else
                                ! pgen is not restricted
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                              b, a, int_switch(2), cum_switch(2))

                                if (b > j) then
                                    ! _R(a) > _LR(i) > ^RL(j) > ^L(b)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_R_to_L, &
                                                gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                                a, j, b, i, a, i, j, b, 0, 4, 1.0_dp, 1.0_dp)

                                else
                                    ! _R(a) > _LR(i)> ^LR(b) > ^R(j)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_R_to_L_to_R, &
                                                gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                a, j, b, i, a, i, b, j, 0, 4, 1.0_dp, 1.0_dp)

                                end if
                            end if
                        else
                            ! there is no restriction on b
                            call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                     int_contrib(2), cum_sum(2), b)

                            if (b == 0) then
                                excitInfo%valid = .false.
                                return
                            end if

                            ! check where b is:
                            if (b == i) then
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                              b, a, int_switch(2), cum_switch(2), i, .true.)

                                if (a > j) then
                                    ! _RL_(ib) > ^RL(j) > ^L(a)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%fullstart_R_to_L, &
                                                gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                                i, j, a, i, i, i, j, a, 0, 2, 1.0_dp, 1.0_dp)

                                else
                                    ! _RL_(ib) > ^LR(a) > ^R(j)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%fullstart_L_to_R, &
                                                gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                i, j, a, i, i, i, a, j, 0, 2, 1.0_dp, 1.0_dp)

                                end if

                            else
                                if (b < i) then
                                    ! a coud not have been i
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2), i)

                                    if (a > j) then
                                        ! _R(b) > _LR(i) > ^RL(j) > ^L(a)
                                        excitInfo = assign_excitInfo_values_double( &
                                                    excit_type%double_R_to_L, &
                                                    gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                                    b, j, a, i, b, i, j, a, 0, 4, 1.0_dp, 1.0_dp)

                                    else
                                        ! _R(b) > _LR(i) > ^LR(a) > ^R(j)
                                        excitInfo = assign_excitInfo_values_double( &
                                                    excit_type%double_R_to_L_to_R, &
                                                    gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                    b, j, a, i, b, i, a, j, 0, 4, 1.0_dp, 1.0_dp)

                                    end if

                                else
                                    ! no restrictions on pgen
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2))

                                    if (a == b) then
                                        if (a > j) then
                                            ! _L(i) > _LL(j) > ^LL^(ab)
                                            excitInfo = assign_excitInfo_values_double( &
                                                        excit_type%fullstop_lowering, &
                                                        gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                                        a, i, a, j, i, j, a, a, 0, 2, 1.0_dp, 1.0_dp)

                                        else
                                            ! _L(i) > ^LR_(ab) > ^R(j)
                                            excitInfo = assign_excitInfo_values_double( &
                                                        excit_type%single_overlap_L_to_R, &
                                                        gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                                        a, i, a, j, i, a, a, j, 0, 2, 1.0_dp, 1.0_dp, 1)

                                        end if

                                    else
                                        ! only extremes coun
                                        st = min(a, b)
                                        en = max(a, b)

                                        if (st > j) then
                                            ! _L(i) > _LL(j) > ^LL(min) > ^L(max)
                                            excitInfo = assign_excitInfo_values_double( &
                                                        excit_type%double_lowering, &
                                                        gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                                        st, i, en, j, i, j, st, en, 0, 4, 1.0_dp, 1.0_dp)

                                        else
                                            if (en > j) then
                                                ! _L(i) > _RL(min) > ^RL(j) > ^L(max)
                                                excitInfo = assign_excitInfo_values_double( &
                                                            excit_type%double_L_to_R_to_L, &
                                                            gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                                            en, i, st, j, i, st, j, en, 0, 4, 1.0_dp, 1.0_dp)

                                            else
                                                ! _L(i) > _RL(min) ^LR(max) > ^R(j)
                                                excitInfo = assign_excitInfo_values_double( &
                                                            excit_type%double_L_to_R, &
                                                            gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                                            en, i, st, j, i, st, en, j, 0, 4, 1.0_dp, 1.0_dp)

                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end if
                else
                    ! [ 1 1 ] config
                    ! in this case, if i dont exclude these cases beforehand
                    ! i have to check if there is a possible switch if
                    ! (a) is (i) or (j)
                    if (a == j) then
                        if (csf_i%stepvector(i) == 1 .and. &
                            csf_i%stepvector(j) == 1) then
                            if (count_alpha_orbs_ij(csf_i, i, j) == 0) then
                                ! no valid excitation
                                excitInfo%valid = .false.
                                return
                            end if
                        else if (csf_i%stepvector(i) == 2 .and. &
                                 csf_i%stepvector(j) == 2) then
                            if (count_beta_orbs_ij(csf_i, i, j) == 0) then
                                excitInfo%valid = .false.
                                return
                            end if
                        end if
                        ! b has to be lower than J
                        call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                 int_contrib(2), cum_sum(2), b, -j, .true.)

                        if (b == 0) then
                            excitInfo%valid = .false.
                            return
                        end if

                        ! have to check where b is to determine switchen
                        ! pgen contribution!
                        if (b == i) then
                            ! a would have to have been > I
                            call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                          a, int_switch(2), cum_switch(2), i, .true.)

                            ! _RL_(ib) > ^RL^(ja)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%fullstart_stop_mixed, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        i, j, j, i, i, i, j, j, 0, 2, 1.0_dp, 1.0_dp)

                        else
                            if (b < i) then
                                ! I is restricted for a
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                              a, int_switch(2), cum_switch(2), i)

                                ! why am i never here??

                                ! _R(b) > _LR(i) > ^RL^(ja)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%fullstop_R_to_L, &
                                            gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                            b, j, j, i, b, i, j, j, 0, 2, 1.0_dp, 1.0_dp)

                            else
                                ! no restrictions
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                              a, int_switch(2), cum_switch(2))

                                ! _L(i) > _RL(b) > ^RL^(ja)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%fullstop_L_to_R, &
                                            gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                            j, i, b, j, i, b, j, j, 0, 2, 1.0_dp, 1.0_dp)

                            end if
                        end if
                    else if (a == i) then
                        ! check if there is a possible switch if both i and j
                        ! have the same stepvalue
                        if (csf_i%stepvector(i) == 1 .and. &
                            csf_i%stepvector(j) == 1) then
                            if (count_alpha_orbs_ij(csf_i, i, j) == 0) then
                                ! no valid excitation
                                excitInfo%valid = .false.
                                return
                            end if
                        else if (csf_i%stepvector(i) == 2 .and. &
                                 csf_i%stepvector(j) == 2) then
                            if (count_beta_orbs_ij(csf_i, i, j) == 0) then
                                excitInfo%valid = .false.
                                return
                            end if
                        end if
                        ! b has to be higher than I
                        call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                 int_contrib(2), cum_sum(2), b, i, .true.)

                        if (b == 0) then
                            excitInfo%valid = .false.
                            return
                        end if

                        ! check where b is

                        if (b == j) then
                            ! a would have to have been < J
                            call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                          a, int_switch(2), cum_switch(2), -j, .true.)

                            ! _RL_(ia) > ^RL^(jb)
                            excitInfo = assign_excitInfo_values_double( &
                                        excit_type%fullstart_stop_mixed, &
                                        gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                        j, i, i, j, i, i, j, j, 0, 2, 1.0_dp, 1.0_dp)
                        else
                            if (b > j) then
                                ! J would be restricted
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                              a, int_switch(2), cum_switch(2), j)

                                ! _RL_(ia) > ^RL(j) > ^L(b)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%fullstart_R_to_L, &
                                            gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                            i, j, b, i, i, i, j, b, 0, 2, 1.0_dp, 1.0_dp)

                            else
                                ! no restrictions
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                              a, int_switch(2), cum_switch(2))

                                ! _RL_(ia) > ^LR(b) > ^R(j)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%fullstart_L_to_R, &
                                            gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                            i, j, b, i, i, i, b, j, 0, 2, 1.0_dp, 1.0_dp)

                            end if
                        end if
                    else
                        ! check were a is
                        if (a > j) then
                            ! b cant be J
                            call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                     int_contrib(2), cum_sum(2), b, j)

                            if (b == 0) then
                                excitInfo%valid = .false.
                                return
                            end if

                            ! check where b is
                            if (b == i) then
                                ! check if there is a possible switch if both i and j
                                ! have the same stepvalue
                                if (csf_i%stepvector(i) == 1 .and. &
                                    csf_i%stepvector(j) == 1) then
                                    if (count_alpha_orbs_ij(csf_i, i, j) == 0) then
                                        ! no valid excitation
                                        excitInfo%valid = .false.
                                        return
                                    end if
                                else if (csf_i%stepvector(i) == 2 .and. &
                                         csf_i%stepvector(j) == 2) then
                                    if (count_beta_orbs_ij(csf_i, i, j) == 0) then
                                        excitInfo%valid = .false.
                                        return
                                    end if
                                end if

                                ! everything below I is restricted
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                              a, int_switch(2), cum_switch(2), i, .true.)

                                ! _RL_(ib) > ^RL(j) > ^L(a)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%fullstart_R_to_L, &
                                            gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                            i, j, a, i, i, i, j, a, 0, 2, 1.0_dp, 1.0_dp)

                            else
                                if (b < i) then
                                    ! I would have been off limits
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2), i)

                                    ! _R(b) > _LR(i) > ^RL(j) > ^L(a)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_R_to_L, &
                                                gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                                b, j, a, i, b, i, j, a, 0, 4, 1.0_dp, 1.0_dp)

                                else if (b > j) then

                                    if (a == b) then
                                        ! J would have been off limits, although in
                                        ! this case i could just copy the probs
                                        ! since a == b .. duh
                                        int_switch(2) = int_contrib(2)
                                        cum_switch(2) = cum_sum(2)

                                        ! and its a:
                                        ! _L(i) -> _LL(j) > ^LL^(a,b)
                                        excitInfo = assign_excitInfo_values_double( &
                                                    excit_type%fullstop_lowering, &
                                                    gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                                    a, i, a, j, i, j, a, a, 0, 2, 1.0_dp, 1.0_dp, 2)

                                    else

                                        ! only extremes count.. and J would have
                                        ! been off-limits
                                        call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                      b, a, int_switch(2), cum_switch(2), j)

                                        st = min(a, b)
                                        en = max(a, b)

                                        ! _L(i) > _LL(j) > ^LL(min) > ^L(max)
                                        excitInfo = assign_excitInfo_values_double( &
                                                    excit_type%double_lowering, &
                                                    gen_type%L, gen_type%L, gen_type%L, gen_type%L, gen_type%L, &
                                                    st, i, en, j, i, j, st, en, 0, 4, 1.0_dp, 1.0_dp)
                                    end if

                                else
                                    ! no restrictions
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2))

                                    ! _L(i) > _RL(b) > ^RL(j) > ^L(a)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_L_to_R_to_L, &
                                                gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                                a, i, b, j, i, b, j, a, 0, 4, 1.0_dp, 1.0_dp)

                                end if
                            end if
                        else if (a < i) then
                            ! b cant be I
                            call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                     int_contrib(2), cum_sum(2), b, i)

                            if (b == 0) then
                                excitInfo%valid = .false.
                                return
                            end if

                            ! check where b is
                            if (b == j) then
                                ! check if there is a possible switch if both i and j
                                ! have the same stepvalue
                                if (csf_i%stepvector(i) == 1 .and. &
                                    csf_i%stepvector(j) == 1) then
                                    if (count_alpha_orbs_ij(csf_i, i, j) == 0) then
                                        ! no valid excitation
                                        excitInfo%valid = .false.
                                        return
                                    end if
                                else if (csf_i%stepvector(i) == 2 .and. &
                                         csf_i%stepvector(j) == 2) then
                                    if (count_beta_orbs_ij(csf_i, i, j) == 0) then
                                        excitInfo%valid = .false.
                                        return
                                    end if
                                end if

                                ! wrong: I would have been off limits
                                ! the above comment is not right i would have
                                ! had to picked something below j, but why
                                ! is it 0? but thats atleast consistent with
                                ! above.. is the umat read in wrong?
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, b, &
                                                              a, int_switch(2), cum_switch(2), -j, .true.)

                                ! _R(a) > _LR(i) > ^RL^(jb)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%fullstop_R_to_L, &
                                            gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                            a, j, j, i, a, i, j, j, 0, 2, 1.0_dp, 1.0_dp)

                            else
                                if (b > j) then
                                    ! J would have been off-limits
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2), j)
                                    ! _R(a) > _LR(i) > ^RL(j) > ^L(b)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_R_to_L, &
                                                gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%L, &
                                                a, j, b, i, a, i, j, b, 0, 4, 1.0_dp, 1.0_dp)

                                else if (b < i) then

                                    ! a can be b! why did i forget that...
                                    if (a == b) then
                                        int_switch(2) = int_contrib(2)
                                        cum_switch(2) = cum_sum(2)

                                        ! _RR_(ab) > ^RR(i) > ^R(j)
                                        excitInfo = assign_excitInfo_values_double( &
                                                    excit_type%fullstart_raising, &
                                                    gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                    a, i, a, j, a, a, i, j, 0, 2, 1.0_dp, 1.0_dp, 2)
                                    else

                                        ! only extremes and I would have been off
                                        ! lmits
                                        call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                      b, a, int_switch(2), cum_switch(2), i)

                                        st = min(a, b)
                                        en = max(a, b)

                                        ! _R(min) > _RR(max) > ^RR(i) > ^R(j)
                                        excitInfo = assign_excitInfo_values_double( &
                                                    excit_type%double_raising, &
                                                    gen_type%R, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                    en, j, st, i, st, en, i, j, 0, 4, 1.0_dp, 1.0_dp)
                                    end if

                                else
                                    ! no restrictions
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2))

                                    ! _R(a) > _LR(i) > ^LR(b) > ^R(j)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_R_to_L_to_R, &
                                                gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                a, j, b, i, a, i, b, j, 0, 4, 1.0_dp, 1.0_dp)

                                end if
                            end if
                        else
                            ! a is between i and j -> no b restrictions
                            call pick_b_orb_guga_mol(csf_i, occ_orbs, a, cc_b, &
                                                     int_contrib(2), cum_sum(2), b)

                            if (b == 0) then
                                excitInfo%valid = .false.
                                return
                            end if
                            if (b < i) then
                                ! I off limits
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                              b, a, int_switch(2), cum_switch(2), i)

                                ! _R(b) > _LR(i) > ^LR(a) > ^R(j)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%double_R_to_L_to_R, &
                                            gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                            b, j, a, i, b, i, a, j, 0, 4, 1.0_dp, 1.0_dp)

                            else if (b > j) then
                                ! J off limits
                                call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                              b, a, int_switch(2), cum_switch(2), j)

                                ! _L(i) > _RL(a) > ^RL(j) > ^L(b)
                                excitInfo = assign_excitInfo_values_double( &
                                            excit_type%double_L_to_R_to_L, &
                                            gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%L, &
                                            b, i, a, j, i, a, j, b, 0, 4, 1.0_dp, 1.0_dp)

                            else

                                ! check where b is
                                if (a == b) then
                                    ! no restrictions on pgen
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2))

                                    ! _L(i) > ^LR_(ab) > ^R(j)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%single_overlap_L_to_R, &
                                                gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                                a, i, a, j, i, a, a, j, 0, 2, 1.0_dp, 1.0_dp, 1)

                                else if (b == i) then
                                    ! a > I
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2), i, .true.)

                                    ! _RL_(ib) > ^LR(a) > ^R(j)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%fullstart_L_to_R, &
                                                gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
                                                i, j, a, i, i, i, a, j, 0, 2, 1.0_dp, 1.0_dp)

                                else if (b == j) then
                                    ! a < J
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2), -j, .true.)

                                    ! _L(i) > _RL(a) > ^RL^(jb)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%fullstop_L_to_R, &
                                                gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                                j, i, a, j, i, a, j, j, 0, 2, 1.0_dp, 1.0_dp)

                                else
                                    ! no restrictions
                                    call pgen_select_orb_guga_mol(csf_i, occ_orbs, &
                                                                  b, a, int_switch(2), cum_switch(2))

                                    ! only extremes count
                                    st = min(a, b)
                                    en = max(a, b)

                                    ! _L(i) > _RL(min) > ^LR(max) > ^R(j)
                                    excitInfo = assign_excitInfo_values_double( &
                                                excit_type%double_L_to_R, &
                                                gen_type%L, gen_type%R, gen_type%L, gen_type%L, gen_type%R, &
                                                en, i, st, j, i, st, en, j, 0, 4, 1.0_dp, 1.0_dp)
                                end if
                            end if
                        end if
                    end if
                end if
            end if
        end if

        ! if the holes are the same spatial orbital it makes no sense of
        ! considering the other order of picking..
        if (a == b) then
            pgen = pgen / 2.0_dp
        end if

        if (b == 1) then
            int_switch(1) = cum_arr(1)
        else
            int_switch(1) = cum_arr(b) - cum_arr(b - 1)
        end if
        cum_switch(1) = cum_sum(1)

        ! have to correctly adress if doubly occupied orbital was chosen for i,j
        ! and think about other stuff too!
        ! but in general it should look like:
        if (any(near_zero(cum_sum)) .or. any(near_zero(cum_switch))) then
            pgen = 0.0_dp
        else
            pgen = pgen * (product(int_contrib) / product(cum_sum) + &
                           product(int_switch) / product(cum_switch))
        end if

    end subroutine pickOrbs_sym_uniform_mol_double

    subroutine pgen_select_orb_guga_mol(csf_i, occ_orbs, orb_b, orb_a, cpt, &
                                        cum_sum, orb_res, range_flag)
        ! routine to recalculate the pgen contribution if orbital (a) and (b)
        ! could have been picked in the opposite order
        ! additional GUGA-restrictions on orbitals are again dealt with
        ! optional input paramters
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_b, orb_a
        real(dp), intent(out) :: cpt, cum_sum
        integer, intent(in), optional :: orb_res
        logical, intent(in), optional :: range_flag
        character(*), parameter :: this_routine = "pgen_select_orb_guga_mol"

        integer :: cc_a, label_index, nOrbs, i, orb
        real(dp) :: tmp
        logical :: tSingle

        ! first get the orbitals from the correct symmetry
        ! again have to settle on a certain spin -> if beta spin was occupied
        ! originally for a -> take alpha spin component
        ! otherwise choose beta
        ! UPDATE: change to only pick spatial orbitals -> only check if orbital
        ! was singly occupied
        if (csf_i%Occ_int(orb_b) == 1) then
            ! then i have to exclude orb_a in the recalculation of p(a|b) prob
            tSingle = .true.
        else
            tSingle = .false.
        end if

        ! need all allowed orbitals independent of actual "spin" of
        ! target orbita
        cc_a = ClassCountInd(2 * orb_a)

        label_index = SymLabelCounts2(1, cc_a)

        nOrbs = OrbClassCount(cc_a)

        cum_sum = 0.0_dp

        if (present(range_flag)) then
            ! if range flag is present, orb_b is already off-limits so no
            ! need to take tSingle into account
            ASSERT(present(orb_res))
            ASSERT(orb_res /= 0)

            if (orb_res < 0) then
                ! only orbitals below restriction
                do i = 1, nOrbs

                    orb = sym_label_list_spat(label_index + i - 1)

                    if (csf_i%stepvector(orb) /= 3 .and. orb < -orb_res) then

                        tmp = get_guga_integral_contrib(occ_orbs, orb_b, orb)

                        cum_sum = cum_sum + tmp
                        if (orb == orb_a) cpt = tmp
                    end if
                end do

            else
                ! only orbitals above restriction
                do i = 1, nOrbs
                    orb = sym_label_list_spat(label_index + i - 1)
                    if (csf_i%stepvector(orb) /= 3 .and. orb > orb_res) then
                        tmp = get_guga_integral_contrib(occ_orbs, orb_b, orb)

                        cum_sum = cum_sum + tmp
                        if (orb == orb_a) cpt = tmp
                    end if
                end do
            end if
        else
            if (present(orb_res)) then
                ASSERT(orb_res /= 0)
                ASSERT(orb_res > 0)
                ! the orbital associated with orb_res is off-limits

                ! also consider if orb b is singly occupied
                if (tSingle) then
                    do i = 1, nOrbs
                        orb = sym_label_list_spat(label_index + i - 1)

                        if (csf_i%stepvector(orb) /= 3 .and. orb /= orb_res &
                            .and. orb /= orb_b) then

                            tmp = get_guga_integral_contrib(occ_orbs, orb_b, orb)

                            cum_sum = cum_sum + tmp
                            if (orb == orb_a) cpt = tmp
                        end if
                    end do
                else
                    ! orb b is also allowed!
                    do i = 1, nOrbs
                        orb = sym_label_list_spat(label_index + i - 1)

                        if (csf_i%stepvector(orb) /= 3 .and. orb /= orb_res) then

                            tmp = get_guga_integral_contrib(occ_orbs, orb_b, orb)

                            cum_sum = cum_sum + tmp
                            if (orb == orb_a) cpt = tmp
                        end if
                    end do
                end if
            else
                ! no restrictions execpt single occupancy maybe
                if (tSingle) then
                    do i = 1, nOrbs
                        orb = sym_label_list_spat(label_index + i - 1)

                        if (csf_i%stepvector(orb) /= 3 .and. orb /= orb_b) then

                            tmp = get_guga_integral_contrib(occ_orbs, orb_b, orb)

                            cum_sum = cum_sum + tmp
                            if (orb == orb_a) cpt = tmp
                        end if
                    end do
                else
                    do i = 1, nOrbs
                        orb = sym_label_list_spat(label_index + i - 1)

                        if (csf_i%stepvector(orb) /= 3) then
                            tmp = get_guga_integral_contrib(occ_orbs, orb_b, orb)

                            cum_sum = cum_sum + tmp

                            if (orb == orb_a) cpt = tmp
                        end if
                    end do
                end if
            end if
        end if

    end subroutine pgen_select_orb_guga_mol

    subroutine pick_b_orb_guga_mol(csf_i, occ_orbs, orb_a, cc_b, int_contrib, &
                                   cum_sum, orb_b, orb_res, range_flag)
        ! restrict the b, if orbital (a) is singly occupied already..
        ! and switch to spatial orbital picking!
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2), orb_a, cc_b
        real(dp), intent(out) :: int_contrib, cum_sum
        integer, intent(out) :: orb_b
        integer, intent(in), optional :: orb_res
        logical, intent(in), optional :: range_flag
        character(*), parameter :: this_routine = "pick_b_orb_guga_mol"

        real(dp) :: cum_arr(OrbClassCount(cc_b)), r
        integer :: nOrbs, label_index, i, orb, orb_index
        logical :: tSingle

        ! rewrite this funcitonality new:
        nOrbs = OrbClassCount(cc_b)

        label_index = SymLabelCounts2(1, cc_b)

        cum_sum = 0.0_dp

        ! have to predetermine is already picked orbital is singly occupied
        ! already: if yes, its not allowd to be picked again in here
        if (csf_i%Occ_int(orb_a) == 1) then
            tSingle = .true.
        else
            tSingle = .false.
        end if

        ! make a spatial orbital list of symmetry allowed orbitals!
        ! but that probably should be done globally in the setup phase, since
        ! this does not change
        ! did that in sym_label_list_spat

        ! depending on input do the specific loops
        if (present(range_flag)) then
            if (range_flag) then
                ASSERT(present(orb_res))
                ASSERT(orb_res /= 0)

                ! the case that a whole range is forbidden due to GUGA restrictions
                ! is only in the case, when orbital (a) is on of the original I,J
                ! orbitals, so in this case a will not be chosen due to the
                ! GUGA restrictions already

                ! the direction of the range is indicated through the sign of the
                ! index restriciton
                if (orb_res < 0) then
                    ! only orbitals below orb_res allowed
                    do i = 1, nOrbs
                        orb = sym_label_list_spat(label_index + i - 1)

                        if (csf_i%stepvector(orb) /= 3 .and. orb < -orb_res) then
                            cum_sum = cum_sum + &
                                      get_guga_integral_contrib(occ_orbs, orb_a, orb)
                        end if

                        cum_arr(i) = cum_sum
                    end do
                else
                    ! only orbitals above orb_res are allowed!
                    do i = 1, nOrbs
                        orb = sym_label_list_spat(label_index + i - 1)

                        if (csf_i%stepvector(orb) /= 3 .and. orb > orb_res) then
                            cum_sum = cum_sum + &
                                      get_guga_integral_contrib(occ_orbs, orb_a, orb)
                        end if

                        cum_arr(i) = cum_sum
                    end do
                end if
            end if
        else
            if (present(orb_res)) then
                ! should i assert here, that orb_res should not be 0?
                ! otherwise it would be stupid to input..
                ASSERT(orb_res /= 0)
                ASSERT(orb_res > 0)
                ! but now i have to include that orb (a) might be off-limits!
                if (tSingle) then
                    ! then orb_a is also off-limits!
                    do i = 1, nOrbs
                        orb = sym_label_list_spat(label_index + i - 1)

                        if (csf_i%stepvector(orb) /= 3 .and. orb /= orb_res &
                            .and. orb /= orb_a) then
                            cum_sum = cum_sum + &
                                      get_guga_integral_contrib(occ_orbs, orb_a, orb)
                        end if

                        cum_arr(i) = cum_sum
                    end do
                else
                    ! orb a is not off-limits
                    do i = 1, nOrbs
                        orb = sym_label_list_spat(label_index + i - 1)

                        if (csf_i%stepvector(orb) /= 3 .and. orb /= orb_res) then
                            cum_sum = cum_sum + &
                                      get_guga_integral_contrib(occ_orbs, orb_a, orb)
                        end if

                        cum_arr(i) = cum_sum

                    end do
                end if
            else
                ! no guga restrictions, only have to check if orb a was single
                if (tSingle) then
                    do i = 1, nOrbs
                        orb = sym_label_list_spat(label_index + i - 1)

                        if (csf_i%stepvector(orb) /= 3 .and. orb /= orb_a) then
                            cum_sum = cum_sum + &
                                      get_guga_integral_contrib(occ_orbs, orb_a, orb)
                        end if

                        cum_arr(i) = cum_sum
                    end do
                else
                    ! no restrictions except double occuations
                    do i = 1, nOrbs
                        orb = sym_label_list_spat(label_index + i - 1)

                        if (csf_i%stepvector(orb) /= 3) then
                            cum_sum = cum_sum + &
                                      get_guga_integral_contrib(occ_orbs, orb_a, orb)
                        end if
                        cum_arr(i) = cum_sum
                    end do
                end if
            end if
        end if

        if (near_zero(cum_sum)) then
            orb_b = 0
            return
        end if

        r = genrand_real2_dSFMT() * cum_sum
        orb_index = binary_search_first_ge(cum_arr, r)

        orb_b = sym_label_list_spat(label_index + orb_index - 1)

        if (orb_index == 1) then
            int_contrib = cum_arr(orb_index)
        else
            int_contrib = cum_arr(orb_index) - cum_arr(orb_index - 1)
        end if

    end subroutine pick_b_orb_guga_mol

    subroutine pick_a_orb_guga_mol(csf_i, occ_orbs, contrib, cum_sum, cum_arr, orb_a)
        ! general routine, which picks orbital a for a  double excitation in
        ! the guga formalism, with symmetry restrictions and weighted
        ! with the FCIDUMP integrals. This is for MOLECULAR calculations,
        ! since in Hubbard and UEG type calculation with existing k-point
        ! restrictions, there is a more efficient and direct way to pick
        ! weighted with the actual matrix elemetn
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: contrib, cum_sum, cum_arr(nSpatOrbs)
        integer, intent(out) :: orb_a

        real(dp) :: r
        ! there are no additional GUGA restrictions on the orbital a yet, only
        ! that it has to be a non-occupied "spin"-orbital
        ! so do it in the same way as already implemented in NECI in
        ! symrandexcit5!

        ! generate the cummulative pgen list:
        if (tGen_guga_weighted) then
            call gen_a_orb_cum_list_guga_mol(csf_i, occ_orbs, cum_arr)
        else
            cum_arr = csf_i%cum_list
        end if

        cum_sum = cum_arr(nSpatOrbs)
        ! check if no excitation is possible
        if (near_zero(cum_sum)) then
            orb_a = 0
            return
        end if

        ! then pick the orbitals according to the list
        r = genrand_real2_dSFMT() * cum_sum
        orb_a = binary_search_first_ge(cum_arr, r)

        if (orb_a == 1) then
            contrib = cum_arr(1)
        else
            contrib = cum_arr(orb_a) - cum_arr(orb_a - 1)
        end if

        ! maybe similar to simon to a DEBUG check if the pgens are correct
        ! TODO

    end subroutine pick_a_orb_guga_mol

    subroutine gen_a_orb_cum_list_guga_mol(csf_i, occ_orbs, cum_arr, tgt_orb, pgen)
        ! subroutine to generate the molecular cumullative probability
        ! distribution. there are no (atleast until now) addiditonal restrictions
        ! or some spin alignement restrictions to generate this list..
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: cum_arr(nSpatOrbs)
        integer, intent(in), optional :: tgt_orb
        real(dp), intent(out), optional :: pgen
        character(*), parameter :: this_routine = "gen_a_orb_cum_list_guga_mol"

        integer :: orb
        real(dp) :: cum_sum
        ! have to think about the ordering of the two electronic indices
        ! in the FCIDUMP integral choosing..
        ! since there is a relative sign, but i have not yet figured out if
        ! i can predetermine this relative sign..

        cum_sum = 0.0_dp

        if (present(tgt_orb)) then
            ASSERT(present(pgen))
            do orb = 1, nSpatOrbs
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)
                end if
                cum_arr(orb) = cum_sum

                if (orb == tgt_orb) then

                end if
            end do

        else
            do orb = 1, nSpatOrbs
                ! only check non-double occupied orbitals
                if (csf_i%stepvector(orb) /= 3) then
                    cum_sum = cum_sum + get_guga_integral_contrib(occ_orbs, orb, -1)
                end if
                cum_arr(orb) = cum_sum
            end do
        end if

    end subroutine gen_a_orb_cum_list_guga_mol

    function get_guga_integral_contrib_spat(occ_orbs, orb_a, orb_b) result(cpt)
        integer, intent(in) :: occ_orbs(2), orb_a, orb_b
        real(dp) :: cpt

        integer :: ind(2)

        ind = occ_orbs
        ! for now, since i dont know how to correctly do it, and for testing
        ! purposes just add a uniformly factor to all of the orbitals.

        ! do a input dependent switch to compare influence on the pgens..
        if (tGen_guga_weighted) then

            if (orb_b < 0) then

                cpt = sqrt(abs_l1(UMat2D(max(ind(1),orb_a), min(ind(1),orb_a)))) &
                    + sqrt(abs_l1(UMat2D(max(ind(2),orb_a), min(ind(2),orb_a))))
            else

                cpt = sqrt(abs(get_umat_el(ind(1),ind(2),orb_a,orb_b)) &
                                 + abs(get_umat_el(ind(1),ind(2),orb_b,orb_a)))

            end if

        else

            if (orb_b < 0) then

                cpt = 1.0_dp

            else

                if (t_guga_pchb) then
                    cpt = abs(get_umat_el(ind(1),ind(2),orb_a,orb_b)) +&
                          abs(get_umat_el(ind(1),ind(2),orb_b,orb_a))

                else
                    cpt = sqrt(abs(get_umat_el(ind(1),ind(2),orb_a,orb_b)) +&
                          abs(get_umat_el(ind(1),ind(2),orb_b,orb_a)))
                  end if
            end if
        end if


    end function get_guga_integral_contrib_spat

    function get_guga_integral_contrib(occ_orbs, orb_a, orb_b) result(cpt)
        ! routine which gets the correct FCIDUMP integral contribution for
        ! orbital a, where electrons i and j are already picked!
        ! have to still figure out how to do that correctly..
        integer, intent(in) :: occ_orbs(2), orb_a, orb_b
        real(dp) :: cpt
        integer :: ind(2)

        ind = gtID(occ_orbs)

        ! ATTENTION! occ_orbs is given in spin orbitals, while orb is a
        ! spatial orbital!

        ! for now, since i dont know how to correctly do it, and for testing
        ! purposes just add a uniformly factor to all of the orbitals.

        ! do a input dependent switch to compare influence on the pgens..
        if (tGen_guga_weighted) then

            if (orb_b < 0) then

                cpt = sqrt(abs_l1(UMat2D(max(ind(1), orb_a), min(ind(1), orb_a)))) &
                      + sqrt(abs_l1(UMat2D(max(ind(2), orb_a), min(ind(2), orb_a))))
            else

                cpt = sqrt(abs(get_umat_el(ind(1), ind(2), orb_a, orb_b)) &
                           + abs(get_umat_el(ind(1), ind(2), orb_b, orb_a)))

            end if

        else

            if (orb_b < 0) then

                cpt = 1.0_dp

            else

                cpt = sqrt(abs(get_umat_el(ind(1), ind(2), orb_a, orb_b)) + &
                           abs(get_umat_el(ind(1), ind(2), orb_b, orb_a)))

            end if
        end if

    end function get_guga_integral_contrib

    subroutine pick_elec_pair_uniform_guga(nI, spin_orbs, sym_prod, sum_ml, &
                                           temp_pgen)
        ! pick two occupied "spin orbitals" uniform-randomly
        integer, intent(in) :: nI(nel)
        integer, intent(out) :: spin_orbs(2), sym_prod, sum_ml

        integer :: i, ind(2)
        real(dp), intent(out) :: temp_pgen

        i = 1 + int(ElecPairs * genrand_real2_dSFMT())

        ind(2) = ceiling((1 + sqrt(1 + 8 * real(i, dp))) / 2)
        ind(1) = i - ((ind(2) - 1) * (ind(2) - 2)) / 2

        spin_orbs = nI(ind)

        sym_prod = RandExcitSymLabelProd(SpinOrbSymLabel(spin_orbs(1)), &
                                         SpinOrbSymLabel(spin_orbs(2)))

        temp_pgen = 1.0_dp / real(ElecPairs, dp)

        sum_ml = sum(G1(spin_orbs)%ml)

    end subroutine pick_elec_pair_uniform_guga

    subroutine pickRandomOrb_restricted(csf_i, start, ende, pgen, orb, occRes)
        ! picks a random orbital from a restricted range start + 1, ende - 1
        ! with optional additional occupation restrictions
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: start, ende
        real(dp), intent(inout) :: pgen
        integer, intent(out) :: orb
        integer, intent(in), optional :: occRes
        character(*), parameter :: this_routine = "pickRandomOrb_restricted"

        integer :: r, nOrbs, ierr
        logical :: mask(nSpatOrbs)
        integer, allocatable :: resOrbs(:)
!         ASSERT(start > 0 .and. start <= nBasis/2)
!         ASSERT(ende > 0 .and. ende <= nBasis/2)

        mask = .true.
        if (present(occRes)) then
            ASSERT(occRes >= 0 .and. occRes <= 2)

            mask = (csf_i%Occ_int /= occRes)

        end if

        mask = (mask .and. (orbitalIndex > start .and. orbitalIndex < ende))

        nOrbs = count(mask)

        if (nOrbs > 0) then
            allocate(resOrbs(nOrbs), stat=ierr)

            resOrbs = pack(orbitalIndex, mask)

            r = 1 + floor(genrand_real2_dSFMT() * real(nOrbs, dp))

            orb = resOrbs(r)

            pgen = pgen / real(nOrbs, dp)

            deallocate(resOrbs)
        else
            orb = 0

            pgen = 0.0_dp
        end if

    end subroutine pickRandomOrb_restricted

    subroutine pickRandomOrb_vector(csf_i, orbRes, pgen, orb, occRes)
        ! this picks a random orb under multiple orbital or occupation
        ! number restrictions
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orbRes(:)
        real(dp), intent(inout) :: pgen
        integer, intent(out) :: orb
        integer, intent(in), optional :: occRes
        character(*), parameter :: this_routine = "pickRandomOrb_vector"

        integer :: r, nOrbs, ierr, num, i
        logical :: mask(nSpatOrbs)
        integer, allocatable :: resOrbs(:)

        ASSERT(all(orbRes >= 0 .and. orbRes <= nSpatOrbs))

        num = size(orbRes)

        ! do not need to check if orbRes is empty, since i call it explicetly
        ! for multiple orbitals
        mask = .true.

        if (present(occRes)) then
            ASSERT(occRes >= 0 .and. occRes <= 2)
            mask = (csf_i%Occ_int /= occRes)
        end if

        do i = 1, num
            mask = (mask .and. orbitalIndex /= orbRes(i))
        end do

        nOrbs = count(mask)

        if (nOrbs > 0) then
            allocate(resOrbs(nOrbs), stat=ierr)

            resOrbs = pack(orbitalIndex, mask)

            r = 1 + floor(genrand_real2_dSFMT() * real(nOrbs, dp))

            orb = resOrbs(r)

            pgen = pgen / real(nOrbs, dp)
            deallocate(resOrbs)
        else
            orb = 0
            pgen = 0.0_dp
        end if

    end subroutine pickRandomOrb_vector

    subroutine pickRandomOrb_forced(csf_i, occRes, pgen, orb, orbRes1)
        ! the version where an orbitals has to have certain occupation.
        ! this never occurs in combination with orbital restrictions!
        ! yes it does!! for fullstart-> fullstop mixed, where i need
        ! n = 1 for both orbitals
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occRes
        real(dp), intent(inout) :: pgen
        integer, intent(out) :: orb
        integer, intent(in), optional :: orbRes1
        character(*), parameter :: this_routine = "pickRandomOrb_forced"

        integer :: r, nOrbs, ierr
        logical :: mask(nSpatOrbs)
        integer, allocatable :: resOrbs(:)

        ASSERT(occRes >= 0 .and. occRes <= 2)

        mask = (csf_i%Occ_int == occRes)

        if (present(orbRes1)) then
            ASSERT(orbRes1 > 0 .and. orbRes1 <= nSpatOrbs)

            mask = (mask .and. orbitalIndex /= orbRes1)
        end if

        nOrbs = count(mask)

        if (nOrbs > 0) then
            allocate(resOrbs(nOrbs), stat=ierr)
            resOrbs = pack(orbitalIndex, mask)

            r = 1 + floor(genrand_real2_dSFMT() * real(nOrbs, dp))

            orb = resOrbs(r)

            pgen = pgen / real(nOrbs, dp)
            deallocate(resOrbs)
        else
            orb = 0
            pgen = 0.0_dp
        end if

    end subroutine pickRandomOrb_forced

    subroutine pickRandomOrb_scalar(csf_i, orbRes, pgen, orb, occRes)
        ! routine to pick a random orbital under certain orbital and/or
        ! occupation number restrictions and gives the probability to pick
        ! this orbital
        ! these orbitals for now are chosen uniformly and not with any
        ! cauchy schwarz like criteria involving the one- and two-particle
        ! integrals
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: orbRes
        real(dp), intent(inout) :: pgen
        integer, intent(out) :: orb
        integer, intent(in), optional :: occRes
        character(*), parameter :: this_routine = "pickRandomOrb_scalar"

        integer :: r, nOrbs, ierr
        logical :: mask(nSpatOrbs)
        integer, allocatable :: resOrbs(:)

        ASSERT(orbRes >= 0 .and. orbRes <= nSpatOrbs)
        ! could make this a global, permanent variable... should even!
!         orbInd = [ (r, r = 1, nBasis/2) ]

        ! create a list of the available orbitals under the provided
        ! restrictions:
        if (orbRes == 0) then
            ! in this case a occRes should be provided or else one could
            ! just create a random number between 1, nBasis/2 without this
            ! function
            ASSERT(present(occRes))
            ASSERT(occRes >= 0 .and. occRes <= 2)

            mask = (csf_i%Occ_int /= occRes)

        else
            ! check i occRes is present
            if (present(occRes)) then
                ASSERT(occRes >= 0 .and. occRes <= 2)
                mask = ((csf_i%Occ_int /= occRes) .and. (orbitalIndex /= orbRes))

            else
                ! only orbital restriction
                mask = (orbitalIndex /= orbRes)

            end if
        end if

        nOrbs = count(mask)
        ! here i could use nOrbs to indicate that no excitation is possible
        ! no -> just check outside if resorbs is size 0
        ! but i get a error further down if i access resOrbs(1)

        if (nOrbs > 0) then
            allocate(resOrbs(nOrbs), stat=ierr)

            resOrbs = pack(orbitalIndex, mask)

            r = 1 + floor(genrand_real2_dSFMT() * real(nOrbs, dp))

            orb = resOrbs(r)

            ! TODO: modify the pgen probabilities!
            ! here it should be 1/nOrbs, since its uniform isnt it?
            pgen = pgen / real(nOrbs, dp)

            deallocate(resOrbs)
        else
            ! indicate that no orbital is fitting with r = 0
            orb = 0

            ! do i have to set pgen to zero? probably
            pgen = 0.0_dp
        end if

    end subroutine pickRandomOrb_scalar

    function excitationIdentifier_double(i, j, k, l) result(excitInfo)
        ! function to identify all necessary information of an excitation
        ! provided with 4 indices of e_{ij,kl}.
        ! determines information on order of indices, involved generator
        ! types, overlap- and non overlap ranges and certain flags, needed for
        ! the correct matrix element calculation. All this information get
        ! stored in a custom type(excitationInformation) defined in the
        ! guga_data module
        integer, intent(in) :: i, j, k, l
        type(ExcitationInformation_t) :: excitInfo
        character(*), parameter :: this_routine = "excitationIdentifier_double"
        integer :: start1, end1, start2, end2

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)
        ASSERT(k <= nSpatOrbs)
        ASSERT(l <= nSpatOrbs)
        ! if accessed with k = l = 0 redirect to single excitation identifier
        ! and exit
        if (k == 0 .or. l == 0) then
            excitInfo = excitationIdentifier(i, j)
            return
        end if
        ! now have to consider all possible i,j,k,l combinations
        excitInfo%i = i
        excitInfo%j = j
        excitInfo%k = k
        excitInfo%l = l

        excitInfo%order = 1.0_dp
        excitInfo%order1 = 1.0_dp

        if (i == j) then
            if (k == l) then
                ! double weight case:
                excitInfo = assign_excitInfo_values_double( &
                            excit_type%weight, &
                            gen_type%W, gen_type%W, gen_type%W, gen_type%W, gen_type%W, i, i, k, k, &
                            0, 0, 0, 0, i, 0, 0.0_dp, 0.0_dp, 0)

            else if (k < l) then
                excitInfo = assign_excitInfo_values_double( &
                            excit_type%raising, &
                            gen_type%W, gen_type%R, gen_type%R, gen_type%R, gen_type%R, i, i, k, l, &
                            k, 0, 0, l, i, 2, 1.0_dp, 0.0_dp, 0)

            else
                excitInfo = assign_excitInfo_values_double( &
                            excit_type%lowering, &
                            gen_type%W, gen_type%L, gen_type%L, gen_type%L, gen_type%L, i, i, k, l, &
                            l, 0, 0, k, i, 2, 1.0_dp, 0.0_dp, 0)

            end if
        else if (k == l) then
            ! other weight combination
            if (i == j) then
                ! double weight case
                excitInfo = assign_excitInfo_values_double( &
                            excit_type%weight, &
                            gen_type%W, gen_type%W, gen_type%W, gen_type%W, gen_type%W, i, i, k, k, &
                            0, 0, 0, 0, k, 0, 0.0_dp, 0.0_dp, 0)

            else if (i < j) then
                excitInfo = assign_excitInfo_values_double( &
                            excit_type%raising, &
                            gen_type%R, gen_type%W, gen_type%R, gen_type%R, gen_type%R, i, j, k, k, &
                            i, 0, 0, j, k, 2, 1.0_dp, 0.0_dp, 0)
            else
                excitInfo = assign_excitInfo_values_double( &
                            excit_type%lowering, &
                            gen_type%L, gen_type%W, gen_type%L, gen_type%L, gen_type%L, i, j, k, l, &
                            j, 0, 0, i, k, 2, 1.0_dp, 0.0_dp, 0)

            end if
        else
            ! no weight generators involved
            start1 = min(i, j)
            end1 = max(i, j)
            start2 = min(k, l)
            end2 = max(k, l)

            excitInfo%fullStart = min(start1, start2)
            excitInfo%fullEnd = max(end1, end2)
            excitInfo%firstEnd = min(end1, end2)
            excitInfo%secondStart = max(start1, start2)

            excitInfo%gen1 = sign(1, j - i)
            excitInfo%gen2 = sign(1, l - k)

            if (excitInfo%firstEnd < excitInfo%secondStart) then
                ! non overlap case

                excitInfo%excitLvl = 4
                excitInfo%typ = excit_type%non_overlap
                excitInfo%overlap = 0
                excitInfo%valid = .true.

                ! maybe need to specify which gen is first and last too
                if (start1 < start2) then
                    excitInfo%firstGen = excitInfo%gen1
                    excitInfo%currentGen = excitInfo%gen1
                    excitInfo%lastGen = excitInfo%gen2

                else
                    excitInfo%firstGen = excitInfo%gen2
                    excitInfo%currentGen = excitInfo%gen2
                    excitInfo%lastGen = excitInfo%gen1
                end if

            else if (excitInfo%firstEnd == excitInfo%secondStart) then
                ! single overlap case:
                excitInfo%overlap = 1

                ! check if that alone is the IC=3 case:
                excitInfo%excitLvl = 3

                excitInfo%valid = .true.

                ! need first and last generators
                if (start1 < start2) then
                    excitInfo%firstGen = excitInfo%gen1
                    excitInfo%currentGen = excitInfo%gen1
                    excitInfo%lastGen = excitInfo%gen2
                else
                    excitInfo%firstGen = excitInfo%gen2
                    excitInfo%currentGen = excitInfo%gen2
                    excitInfo%lastGen = excitInfo%gen1
                end if

                if (excitInfo%firstGen == gen_type%L .and. &
                    excitInfo%lastGen == gen_type%L) then
                    excitInfo%typ = excit_type%single_overlap_lowering

                else if (excitInfo%firstGen == gen_type%R .and. &
                         excitInfo%lastGen == gen_type%R) then
                    excitInfo%typ = excit_type%single_overlap_raising

                else if (excitInfo%firstGen == gen_type%L .and. &
                         excitInfo%lastGen == gen_type%R) then
                    excitInfo%typ = excit_type%single_overlap_L_to_R

                else
                    excitInfo%typ = excit_type%single_overlap_R_to_L

                end if
            else
                ! proper overlap case:
                ! more to determine here...

                ! overlap, and non overlap easiest propably
                ! num overlap entries:
                excitInfo%overlap = excitInfo%firstEnd - excitInfo%secondStart + 1

                excitInfo%valid = .true.

                ! for generator only have to specify which ones are acting in
                ! the non-overlap region, since naturally both of them are
                ! acting in the overlap region simultaniously
                if (start1 < start2) then
                    excitInfo%firstGen = excitInfo%gen1
                    excitInfo%currentGen = excitInfo%gen1

                    if (end1 > end2) then
                        excitInfo%lastGen = excitInfo%gen1
                        if (excitInfo%gen1 == gen_type%L .and. &
                            excitInfo%gen2 == gen_type%L) then
                            excitInfo%typ = excit_type%double_lowering
                            ! here only semi-stop has sign
                            excitInfo%order1 = -1.0_dp

                        else if (excitInfo%gen1 == gen_type%R .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%double_raising
                            ! in this case there are sign changes only at the
                            ! semi-start
                            excitInfo%order = -1.0_dp

                        else if (excitInfo%gen1 == gen_type%L .and. excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%double_L_to_R_to_L

                        else
                            excitInfo%typ = excit_type%double_R_to_L_to_R

                        end if

                    else if (end1 < end2) then
                        excitInfo%lastGen = excitInfo%gen2

                        if (excitInfo%gen1 == gen_type%L .and. &
                            excitInfo%gen2 == gen_type%L) then
                            excitInfo%typ = excit_type%double_lowering
                            ! here no semi has a sign

                        else if (excitInfo%gen1 == gen_type%R .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%double_raising
                            ! here both semi-start and stop have a sign
                            excitInfo%order = -1.0_dp
                            excitInfo%order1 = -1.0_dp

                        else if (excitInfo%gen1 == gen_type%L .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%double_L_to_R

                        else
                            excitInfo%typ = excit_type%double_R_to_L

                        end if

                    else
                        ! set lastGen to gen2 just to make same comparisons
                        excitInfo%lastGen = excitInfo%gen2
                        if (excitInfo%gen1 == gen_type%L .and. &
                            excitInfo%gen2 == gen_type%L) then
                            excitInfo%typ = excit_type%fullstop_lowering

                        else if (excitInfo%gen1 == gen_type%R .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%fullstop_raising

                        else if (excitInfo%gen1 == gen_type%L .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%fullstop_L_to_R

                        else
                            excitInfo%typ = excit_type%fullstop_R_to_L

                        end if

                    end if

                else if (start1 > start2) then
                    excitInfo%firstGen = excitInfo%gen2
                    excitInfo%currentGen = excitInfo%gen2

                    if (end1 > end2) then
                        excitInfo%lastGen = excitInfo%gen1
                        if (excitInfo%gen1 == gen_type%L .and. &
                            excitInfo%gen2 == gen_type%L) then
                            excitInfo%typ = excit_type%double_lowering
                            ! here both have a sign
                            excitInfo%order = -1.0_dp
                            excitInfo%order1 = -1.0_dp

                        else if (excitInfo%gen1 == gen_type%R .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%double_raising
                            ! here both have "normal" sign

                        else if (excitInfo%gen1 == gen_type%L .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%double_R_to_L

                        else
                            excitInfo%typ = excit_type%double_L_to_R

                        end if
                    else if (end1 < end2) then
                        excitInfo%lastGen = excitInfo%gen2
                        if (excitInfo%gen1 == gen_type%L .and. &
                            excitInfo%gen2 == gen_type%L) then
                            excitInfo%typ = excit_type%double_lowering
                            ! here only semi-start has a sign
                            excitInfo%order = -1.0_dp

                        else if (excitInfo%gen1 == gen_type%R .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%double_raising
                            ! here only semi-stop has sign
                            excitInfo%order1 = -1.0_dp

                        else if (excitInfo%gen1 == gen_type%L .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%double_R_to_L_to_R

                        else
                            excitInfo%typ = excit_type%double_L_to_R_to_L

                        end if

                    else
                        ! set lastGen to gen1 just to make same comparisons
                        excitInfo%lastGen = excitInfo%gen1
                        if (excitInfo%gen1 == gen_type%L .and. &
                            excitInfo%gen2 == gen_type%L) then
                            excitInfo%typ = excit_type%fullstop_lowering

                        else if (excitInfo%gen1 == gen_type%R .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%fullstop_raising

                        else if (excitInfo%gen1 == gen_type%L .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%fullstop_R_to_L

                        else
                            excitInfo%typ = excit_type%fullstop_L_to_R

                        end if
                    end if

                else
                    if (end1 > end2) then
                        excitInfo%lastGen = excitInfo%gen1
                        ! set first gen fake to other, to compare it in the
                        ! same way
                        excitInfo%firstGen = excitInfo%gen2

                        if (excitInfo%gen1 == gen_type%L .and. &
                            excitInfo%gen2 == gen_type%L) then
                            excitInfo%typ = excit_type%fullstart_lowering

                        else if (excitInfo%gen1 == gen_type%R .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%fullstart_raising

                        else if (excitInfo%gen1 == gen_type%L .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%fullstart_R_to_L

                        else
                            excitInfo%typ = excit_type%fullStart_L_to_R

                        end if
                    else if (end1 < end2) then
                        excitInfo%lastGen = excitInfo%gen2
                        excitInfo%firstGen = excitInfo%gen1
                        excitInfo%currentGen = excitInfo%gen1
                        if (excitInfo%gen1 == gen_type%L .and. &
                            excitInfo%gen2 == gen_type%L) then
                            excitInfo%typ = excit_type%fullstart_lowering

                        else if (excitInfo%gen1 == gen_type%R .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%fullstart_raising

                        else if (excitInfo%gen1 == gen_type%L .and. &
                                 excitInfo%gen2 == gen_type%R) then
                            excitInfo%typ = excit_type%fullStart_L_to_R

                        else
                            excitInfo%typ = excit_type%fullstart_R_to_L

                        end if
                    else
                        ! check generator types here too.
                        if (excitInfo%gen1 == excitInfo%gen2) then
                            excitInfo%typ = excit_type%fullstart_stop_alike

                        else
                            excitInfo%typ = excit_type%fullstart_stop_mixed
                        end if
                    end if

                end if

                ! TODO: concerning the order flag, there has to be a
                ! decision made. -> todo later

            end if
        end if

    end function excitationIdentifier_double

    function excitationIdentifier_single(i, j) result(excitInfo)
        ! function to identify all necessary information to calculate a
        ! single excitation, provided the two indices i,j. And also to be
        ! able to calculate the matrix elements correctly.
        ! in the single excitation case, there is very little information
        ! necessary, but to keep it coherently with the rest also to it this
        ! way
        ! possible entries for excitation type excitInfo%typ:
        ! single raising
        ! single lowering
        ! single weight
        integer, intent(in) :: i, j
        type(ExcitationInformation_t) :: excitInfo
        character(*), parameter :: this_routine = "excitationIdentifier_single"

        ASSERT(i > 0 .and. i <= nSpatOrbs)
        ASSERT(j > 0 .and. j <= nSpatOrbs)

        ! type of excitation, be more specific here and use the available
        ! information in the calculation of the excitations also.

        ! identify generator, excitation level and start and end indices and
        ! also store nonOverlapRange although maybe not even needed...
        if (i < j) then
            excitInfo = assign_excitInfo_values_single(1, i, j, i, j)

        else if (i > j) then
            excitInfo = assign_excitInfo_values_single(gen_type%L, i, j, j, i)

        else
            excitInfo = assign_excitInfo_values_single(0, i, i, i, i)

        end if

    end function excitationIdentifier_single


    subroutine calcRemainingSwitches_single(csf_i, sOrb, pOrb, &
                                            posSwitches, negSwitches)
        ! subroutine to determine the number of remaining switches for single
        ! excitations between orbitals s, p given in type of excitationInformation.
        ! The switches are given
        ! as a list, to access it for each spatial orbital
        ! stepValue = 1 -> positive delta B switch possibility
        ! stepValue = 2 -> negative delta B switch possibility
        ! assume exitInfo is already calculated
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: sOrb, pOrb
        real(dp), intent(out) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        integer :: iOrb
        real(dp) :: oneCount, twoCount

        ! ignore b > 0 forced switches for now. As they only change the bias
        ! do not make the excitations calculations wrong if ignored.

        oneCount = 0.0_dp
        twoCount = 0.0_dp
        posSwitches = 0.0_dp
        negSwitches = 0.0_dp
        ! have to count from the reversed ilut entries
        do iOrb = max(sOrb, pOrb) - 1, min(sOrb, pOrb), -1
            posSwitches(iOrb) = twoCount
            negSwitches(iOrb) = oneCount

            select case (csf_i%stepvector(iOrb))
            case (1)
                oneCount = oneCount + 1.0_dp
            case (2)
                twoCount = twoCount + 1.0_dp
            end select
        end do

    end subroutine calcRemainingSwitches_single

    function assign_excitInfo_values_double(typ, gen1, gen2, currentGen, firstGen, &
                                            lastGen, i, j, k, l, fullStart, secondStart, firstEnd, fullEnd, &
                                            weight, excitLvl, order, order1, overlap) result(excitInfo)
        integer, intent(in) :: typ, gen1, gen2, currentGen, firstGen, lastGen, &
                               i, j, k, l, fullStart, secondStart, firstEnd, &
                               fullEnd, weight, excitLvl
        integer, intent(in), optional :: overlap
        real(dp), intent(in) :: order, order1
        type(ExcitationInformation_t) :: excitInfo

        ! todo: asserts!
        excitInfo%typ = typ
        excitInfo%gen1 = gen1
        excitInfo%gen2 = gen2
        excitInfo%currentGen = currentGen
        excitInfo%firstGen = firstGen
        excitInfo%lastGen = lastGen
        excitInfo%i = i
        excitInfo%j = j
        excitInfo%k = k
        excitInfo%l = l
        excitInfo%fullStart = fullStart
        excitInfo%secondStart = secondStart
        excitInfo%firstEnd = firstEnd
        excitInfo%fullEnd = fullEnd
        excitInfo%weight = weight
        excitInfo%excitLvl = excitLvl
        excitInfo%order = order
        excitInfo%order1 = order1
        if (present(overlap)) then
            excitInfo%overlap = overlap
        else
            excitInfo%overlap = 2
        end if

        excitInfo%valid = .true.

    end function assign_excitInfo_values_double

    function assign_excitInfo_values_single(gen, i, j, fullStart, fullEnd, typ) &
        result(excitInfo)
        integer, intent(in) :: gen, i, j, fullStart, fullEnd
        integer, intent(in), optional :: typ
        type(ExcitationInformation_t) :: excitInfo

        ! set default values for single excitations: which cause errors if
        ! called in incorrect places
        if (present(typ)) then
            excitInfo%typ = typ
        else
            excitInfo%typ = excit_type%single
        end if

        if (i == j) then
            excitInfo%excitLvl = 0
            excitInfo%weight = i
        else
            excitInfo%excitLvl = 2
            excitInfo%weight = 0
        end if
        excitInfo%k = 0
        excitInfo%l = 0
        excitInfo%secondStart = 0
        excitInfo%firstEnd = 0
        excitInfo%gen2 = -2
        excitInfo%order = 0.0_dp
        excitInfo%order1 = 0.0_dp
        excitInfo%overlap = 0

        ! then set proper values
        excitInfo%i = i
        excitInfo%j = j
        excitInfo%gen1 = gen
        excitInfo%fullStart = fullStart
        excitInfo%fullEnd = fullEnd
        excitInfo%currentGen = gen
        excitInfo%firstGen = gen
        excitInfo%lastGen = gen

        excitInfo%valid = .true.

    end function assign_excitInfo_values_single

    subroutine print_indices(excitInfo)
        type(ExcitationInformation_t), intent(in) :: excitInfo

        print *, "ijkl:", excitInfo%i, excitInfo%j, excitInfo%k, excitInfo%l

    end subroutine print_indices

    subroutine print_excitInfo(excitInfo)
        type(ExcitationInformation_t), intent(in) :: excitInfo

        print *, "Excitation Information: "
        print *, "Typ: ", excitInfo%typ, trim(excit_names(excitInfo%typ)%str)
        print *, "i,j,k,l:"
        print *, excitInfo%i, excitInfo%j,excitInfo%k,excitInfo%l
        print *, "fullStart,secondStart,firstEnd,fullEnd:"
        print *, excitInfo%fullStart, excitInfo%secondStart, &
            excitInfo%firstEnd, excitInfo%fullEnd
        print *, "gen1,gen2,currentGen "
        print *, excitInfo%gen1, excitInfo%gen2, excitInfo%currentGen
        print *, "firstGen,lastGen: "
        print *, excitInfo%firstGen, excitInfo%lastGen
        print *, "valid? ", excitInfo%valid
        print *, "spin_change?", excitInfo%spin_change
        print *, "orders:", excitInfo%order, excitInfo%order1
        print *, "excit-level: ", excitInfo%excitLvl

    end subroutine print_excitInfo

    function checkCompatibility_single(csf_i, excitInfo) result(flag)
        use guga_types, only: WeightObj_t

        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        logical :: flag

        debug_function_name("checkCompatibility_single")

        real(dp) :: pw, mw, posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)
        integer ::  st, en
        type(WeightObj_t) :: weights

        ASSERT(excitInfo%typ == excit_type%single)
        ASSERT(excitInfo%gen1 == gen_type%R .or. excitInfo%gen1 == gen_type%L)

        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitches, negSwitches)

        st = excitInfo%fullStart
        en = excitInfo%fullEnd
        flag = .true.

        if (excitInfo%gen1 == gen_type%R) then
            ! raising
            if (csf_i%stepvector(st) == 3 .or. csf_i%stepvector(en) == 0) then
                flag = .false.
                return
            end if

            weights = init_singleWeight(csf_i, en)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)
            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)

            if ((near_zero(pw) .and. near_zero(mw)) &
                .or. (csf_i%stepvector(st) == 1 .and. near_zero(pw)) &
                .or. (csf_i%stepvector(st) == 2 .and. near_zero(mw))) then
                flag = .false.
                return
            end if

        else if (excitInfo%gen1 == gen_type%L) then
            ! lowering

            if (csf_i%stepvector(en) == 3 .or. csf_i%stepvector(st) == 0) then
                flag = .false.
                return
            end if

            weights = init_singleWeight(csf_i, en)
            mw = weights%proc%minus(negSwitches(st), csf_i%B_real(st), weights%dat)
            pw = weights%proc%plus(posSwitches(st), csf_i%B_real(st), weights%dat)


            if ((csf_i%stepvector(st) == 1 .and. near_zero(pw)) &
                .or. (csf_i%stepvector(st) == 2 .and. near_zero(mw)) &
                .or. (near_zero(pw + mw))) then
                flag = .false.
                return
            end if
        end if
    end function checkCompatibility_single

end module guga_excitations
