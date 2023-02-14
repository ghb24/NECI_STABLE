#include "macros.h"
! GUGA excitations module:
! contains as much excitation related functionality as possible
module guga_main
    use CalcData, only: t_trunc_guga_pgen, &
                        trunc_guga_pgen, t_trunc_guga_pgen_noninits, &
                        t_guga_back_spawn, n_guga_back_spawn_lvl, t_guga_back_spawn_trunc, &
                        t_matele_cutoff, matele_cutoff

    use SystemData, only: nEl, nSpatOrbs, tHub, treal, tgen_guga_crude, &
                          tgen_guga_mixed, t_new_real_space_hubbard, &
                          t_crude_exchange, t_crude_exchange_noninits, &
                          t_approx_exchange, t_approx_exchange_noninits, &
                          is_init_guga, t_heisenberg_model, t_tJ_model, t_mixed_hubbard, &
                          modk_offdiag

    use constants, only: dp, n_int, bn2_, maxExcit, stdout

    use bit_reps, only: decode_bit_det

    use bit_rep_data, only: GugaBits, niftot, nifguga

    use procedure_pointers, only: get_umat_el

    use dSFMT_interface, only: genrand_real2_dSFMT

    use FciMCData, only: excit_gen_store_type, pSingles, pDoubles, &
                         tFillingStochRDMOnFly

    use OneEInts, only: GetTMatEl

    use util_mod, only: near_zero, stop_all

    use GenRandSymExcitNUMod, only: IsMomentumAllowed

    use back_spawn, only: check_electron_location, check_orbital_location, &
                          check_electron_location_spatial, check_orbital_location_spatial

    use util_mod, only: neci_flush

    use guga_procedure_pointers, only: pickOrbitals_single, pickOrbitals_double

    use guga_data, only: &
        ExcitationInformation_t, getDoubleContribution, excit_type

    use guga_bitRepOps, only: isProperCSF_ilut, getDeltaB, &
                              encode_matrix_element, update_matrix_element, &
                              extract_matrix_element, write_det_guga, convert_ilut_toGUGA, &
                              convert_ilut_toNECI, identify_excitation, &
                              extract_h_element, encode_stochastic_rdm_info, &
                              CSF_Info_t, is_compatible, current_csf_i

    use guga_matrixElements, only: &
        calc_guga_matrix_element, calc_integral_contribution_single

    use guga_types, only: WeightObj_t

    use guga_excitations, only: global_excitInfo, print_excitInfo, checkcompatibility, &
        calcsingleoverlapmixedstochastic, calcdoubleloweringstochastic, &
        calcdoubleraisingstochastic, calcdoublel2r2l_stochastic, &
        calcdoubler2l2r_stochastic, calcdoublel2r_stochastic, &
        calcdoubler2l_stochastic, calcfullstoploweringstochastic, &
        calcfullstopl2r_stochastic, calcfullstopraisingstochastic, &
        calcfullstopr2l_stochastic, calcfullstartloweringstochastic, &
        calcfullstartraisingstochastic, calcfullstartl2r_stochastic, &
        calcfullstartr2l_stochastic, calcfullstartfullstopalike, &
        calcfullstartfullstopmixedstochastic, calcremainingswitches_excitinfo_single, &
        init_singleWeight, createstochasticstart_single, &
        singlestochasticupdate, singleStochasticEnd

    use guga_bitRepOps, only: contract_1_rdm_ind

    use guga_crude_approx_mod, only: create_crude_guga_double, create_crude_double, &
        perform_crude_excitation, create_crude_guga_single, create_crude_single


    better_implicit_none

    private
    public :: generate_excitation_guga, &
        createStochasticExcitation_single, createStochasticExcitation_double

contains
    subroutine generate_excitation_guga(nI, ilutI, nJ, ilutJ, exFlag, IC, &
                                        excitMat, tParity, pgen, HElGen, store, part_type)
        !! An API interfacing function for generate_excitation to the rest of NECI:
        !!
        !! Requires guga_bitRepOps::current_csf_i to be set according to the ilutI.
        integer, intent(in) :: nI(nEl), exFlag
        integer(n_int), intent(in) :: ilutI(0:niftot)
        integer, intent(out) :: nJ(nEl), IC, excitMat(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:niftot)
        logical, intent(out) :: tParity
        real(dp), intent(out) :: pgen
        HElement_t(dp), intent(out) :: HElGen
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: part_type
        character(*), parameter :: this_routine = "generate_excitation_guga"

        integer(n_int) :: ilut(0:nifguga), new_ilut(0:nifguga)
        integer :: excit_typ(2)

        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: diff
        HElement_t(dp) :: tmp_mat1
        HElement_t(dp) :: tmp_mat

        unused_var(exFlag); unused_var(part_type); unused_var(store)
        ASSERT(is_compatible(ilutI, current_csf_i))

        ! think about default values and unneeded variables for GUGA, but
        ! which have to be processed anyway to interface to NECI

        ! excitatioin matrix... i could set that up for GUGA too..
        ! but its not needed specifically except for RDM and logging purposes

        ! in new implementation with changing relative probabilites of different
        ! types of excitation, misuse this array to log the type of excitation
        excitMat = 0

        ! the parity flag is also unneccesary in GUGA
        tParity = .true.

        ! the inputted exFlag variable, is also not needed probably..

        ! then choose between single or double excitations..
        ! TODO: still have to think how to set up pSingles and pDoubles in
        ! GUGA...

        ! and before i have to convert to GUGA iluts..
        call convert_ilut_toGUGA(ilutI, ilut)

        ASSERT(isProperCSF_ilut(ilut, .true.))

        ! maybe i need to copy the flags of ilutI onto ilutJ
        ilutJ = ilutI

        if (genrand_real2_dSFMT() < pSingles) then

            IC = 1
            call createStochasticExcitation_single(ilut, nI, current_csf_i, new_ilut, pgen)
            pgen = pgen * pSingles

        else

            IC = 2
            call createStochasticExcitation_double(ilut, nI, current_csf_i, new_ilut, pgen, excit_typ)
            pgen = pgen * pDoubles

        end if

        ! for now add a sanity check to compare the stochastic obtained
        ! matrix elements with the exact calculation..
        ! since something is going obviously wrong..
#ifdef DEBUG_
        call additional_checks()
#endif

        ! check if excitation generation was successful
        if (near_zero(pgen)) then
            ! indicate NullDet to skip spawn step
            nJ(1) = 0
            HElGen = h_cast(0.0_dp)

        else

            ! also store information on type of excitation for the automated
            ! tau-search for the non-weighted guga excitation generator in
            ! the excitMat variable
            excitMat(1, 1:2) = excit_typ

            ! profile tells me this costs alot of time.. so remove it
            ! and only do it in debug mode..
            ! i just have to be sure that no wrong csfs are created..

            ASSERT(isProperCSF_ilut(new_ilut, .true.))
            ! otherwise extract H element and convert to 0

            call convert_ilut_toNECI(new_ilut, ilutJ, HElgen)

            if (t_matele_cutoff .and. abs(HElGen) < matele_cutoff) then
                HElgen = h_cast(0.0_dp)
                nJ(1) = 0
                pgen = 0.0_dp
                return
            end if

            call decode_bit_det(nJ, ilutJ)

            if (tHub .and. .not. treal) then
                if (.not. (IsMomentumAllowed(nJ))) then
                    call write_det_guga(stdout, new_ilut)
                end if
            end if
        end if

        if (modk_offdiag) HElgen = -abs(HElgen)

    contains

#ifdef DEBUG_
        subroutine additional_checks()

            ! for now add a sanity check to compare the stochastic obtained
            ! matrix elements with the exact calculation..
            ! since something is going obviously wrong..
            if (.not. near_zero(pgen)) then
                call convert_ilut_toNECI(new_ilut, ilutJ, HElgen)
                call calc_guga_matrix_element(ilutI, current_csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, tmp_mat, .true.)

                diff = abs(HElGen - tmp_mat)
                if (diff > 1.0e-10_dp) then
                    print *, "WARNING: differing stochastic and exact matrix elements!"
                    call write_det_guga(stdout, ilutI, .true.)
                    call write_det_guga(stdout, ilutJ, .true.)
                    print *, "mat eles and diff:", HElGen, tmp_mat, diff
                    print *, " pgen: ", pgen
                    print *, " deduced excit-info: "
                    call print_excitInfo(excitInfo)
                    print *, " global excit-info: "
                    call print_excitInfo(global_excitInfo)
                    call neci_flush(stdout)
                end if

                ! is the other order also fullfilled?
                call calc_guga_matrix_element(ilutJ, CSF_Info_t(ilutJ), ilutI, current_csf_i, excitInfo, tmp_mat1, .true.)

#ifdef CMPLX_
                diff = abs(tmp_mat1 - conjg(tmp_mat))
#else
                diff = abs(tmp_mat1 - tmp_mat)
#endif
                if (diff > 1.0e-10_dp) then
                    print *, "WARNING: differing sign in matrix elements!"
                    call write_det_guga(stdout, ilutI, .true.)
                    call write_det_guga(stdout, ilutJ, .true.)
                    print *, "mat eles and diff:", tmp_mat, tmp_mat1, diff
                    print *, "<I|H|J> excitInfo:"
                    call print_excitInfo(excitInfo)
                    excitInfo = identify_excitation(ilutI, ilutJ)
                    print *, "<J|H|I> excitInfo:"
                    call print_excitInfo(excitInfo)
                    call neci_flush(stdout)
                end if
            end if
        end subroutine additional_checks
#endif
    end subroutine generate_excitation_guga

    subroutine createStochasticExcitation_single(ilut, nI, csf_i, exc, pgen)
        ! calculate one possible single excitation and the corresponding
        ! probabilistic weight and hamilton matrix element for a given CSF
        ! store matrix element in ilut for now... maybe change that later
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        integer(n_int), intent(out) :: exc(0:nifguga)
        real(dp), intent(out) :: pgen
        character(*), parameter :: this_routine = "createStochasticExcitation_single"

        type(ExcitationInformation_t) :: excitInfo
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs), &
                    branch_pgen, orb_pgen, temp_pgen
        HElement_t(dp) :: integral, mat_ele

        type(WeightObj_t) :: weights
        integer :: iO, st, en, step, i, j, gen, deltaB, step2
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)

        ASSERT(isProperCSF_ilut(ilut))

        ! first pick possible orbitals:
        ! todo: combine that with integrals elements... otherwise it does not
        ! make too much sense
        ! but since there is a sum of contribution, which can lead to the
        ! same excitaiton shouldn't i check if the whole sum is zero, and not
        ! only the one-particle matrix element.

        ! this pickOrbitals picker is the only difference in the symmetric
        ! and non-symmetric excitation generator
        ! so in the initialize function point a general orbital picker to
        ! the specific one, depending on the input..
        call pickOrbitals_single(ilut, nI, csf_i, excitInfo, orb_pgen)

        if (.not. excitInfo%valid) then
            ! if no valid indices were picked, return 0 excitation and return
            exc = 0_n_int
            pgen = 0.0_dp
            return
        end if

        if (t_guga_back_spawn) then
            ! do smth like this:
            ! if i find to increase the excit-lvl with the chosen
            ! orbitals and the current CSF is a non-initiator ->
            ! perform a crude excitation
            if (increase_ex_levl(csf_i, excitInfo) .and. .not. is_init_guga) then

                if (t_guga_back_spawn_trunc) then
                    ! a 2 indicated we want to cancel excit-lvl increasing
                    ! excitations.
                    pgen = 0.0_dp
                    exc = 0_n_int
                    return
                end if

                call create_crude_guga_single(ilut, nI, csf_i, exc, branch_pgen, excitInfo)

                ! there is also this routine I already wrote:
                ! I should combine those two as they do the same job

                pgen = orb_pgen * branch_pgen

                return
            end if
        end if

        ! do the crude approximation here for now..
        if (tgen_guga_crude .and. .not. tgen_guga_mixed) then

            call create_crude_single(ilut, csf_i, exc, branch_pgen, excitInfo)

            if (near_zero(branch_pgen)) then
                exc = 0_n_int
                pgen = 0.0_dp
                return
            end if

            call convert_ilut_toNECI(ilut, ilutI)
            call convert_ilut_toNECI(exc, ilutJ)

            call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, .true.)

            if (near_zero(mat_ele)) then
                exc = 0
                pgen = 0.0_dp

                return
            end if

            call encode_matrix_element(exc, 0.0_dp, 2)
            call encode_matrix_element(exc, mat_ele, 1)

            pgen = orb_pgen * branch_pgen

            return
        end if

        ! have to have these short variable names or otherwise compilation
        ! fails due to line-length restrictions with is Three(), etc. macros
        st = excitInfo%fullStart
        en = excitInfo%fullEnd

        ! reimplement it from scratch
        i = excitInfo%i
        j = excitInfo%j
        ! first the "normal" contribution
        ! not sure if i have to subtract that element or not...
        integral = getTmatEl(2 * i, 2 * j)

        ! ! then calculate the remaing switche given indices
        call calcRemainingSwitches_excitInfo_single(csf_i, excitInfo, posSwitches, negSwitches)
        ! intitialize the weights
        weights = init_singleWeight(csf_i, excitInfo%fullEnd)

        ! create the start randomly(if multiple possibilities)
        ! create the start in such a way to use it for double excitations too
        call createStochasticStart_single(ilut, csf_i, excitInfo, weights, posSwitches, &
                                          negSwitches, exc, branch_pgen)

        ! can it be zero here? maybe due to matrix element issues...
        if (near_zero(branch_pgen)) then
            pgen = 0.0_dp
            exc = 0_n_int
            return
        end if

        ! update at last...

        gen = excitInfo%currentGen

        ! then do the stochastic updates..
        do iO = st + 1, en - 1

            ! need the ongoing deltaB value to access the multFactor table in
            ! the same way as single and double excitations..
            deltaB = getDeltaB(exc)

            call singleStochasticUpdate(ilut, csf_i, iO, excitInfo, weights, posSwitches, &
                                        negSwitches, exc, temp_pgen)

            branch_pgen = branch_pgen * temp_pgen
            ! also get the double contribution during this loop
            ! depends on both stepvalues...

            if (t_trunc_guga_pgen .or. &
                (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
                if (branch_pgen < trunc_guga_pgen) then
                    pgen = 0.0_dp
                    exc = 0_n_int
                    return
                end if
            end if

            ! how to modifiy the values?
            ! one of it is just additional
            if (.not. (treal .or. t_new_real_space_hubbard .or. t_mixed_hubbard &
                       .or. t_tJ_model)) then
                integral = integral + get_umat_el(i, iO, j, iO) * csf_i%Occ_real(iO)

                ! but the r_k part confuses me a bit ...
                step = csf_i%stepvector(iO)
                step2 = getStepvalue(exc, iO)

                integral = integral + get_umat_el(i, iO, iO, j) * &
                           getDoubleContribution(step2, step, deltaB, gen, csf_i%B_real(iO))
            end if
        end do

        ! the end step should be easy in this case. since due to the
        ! correct use of probabilistic weights only valid excitations should
        ! come to this point
        call singleStochasticEnd(csf_i, excitInfo, exc)

        ! maybe but a check here if the matrix element anyhow turned out zero
        if (near_zero(extract_matrix_element(exc, 1))) then
            pgen = 0.0_dp
            exc = 0_n_int
            return
        end if

        pgen = orb_pgen * branch_pgen

        ! do all the integral calulation afterwards..
        ! since it depends on the created excitation.
        ! updates integral variable:
        ! should i intermediately but a if-statement for the real-space
        ! hubbard model here? or should i write a more efficient
        ! single-excitation creator? i guess i should..
        ! and also for the matrix element calculation maybe..
        if (.not. (treal .or. t_new_real_space_hubbard .or. &
                   t_heisenberg_model .or. t_tJ_model .or. t_mixed_hubbard)) then
            ! TODO(@Oskar): Don't calculate here
            call calc_integral_contribution_single(csf_i, CSF_Info_t(exc), i, j, st, en, integral)
        end if

        if (tFillingStochRDMOnFly) then
            ! if we want to do 'fast' GUGA RDMs we need to store the
            ! rdm index and the x0 (for singles here) coupling coefficient
            ! as part of the ilut(0:nifguga). this also necessitates
            ! a 'longer' nifguga (+3 i think for rdm_index and x0 and x1..)
            ! with an accompanying change to niftot.. (this could get a bit
            ! messy in the rest of the code..)
            call encode_stochastic_rdm_info(GugaBits, exc, rdm_ind= &
                                            contract_1_rdm_ind(i, j, excit_lvl=1, excit_typ=excit_type%single), &
                                            x0=extract_matrix_element(exc, 1), x1=0.0_dp)
        end if

        call encode_matrix_element(exc, 0.0_dp, 2)
        call update_matrix_element(exc, integral, 1)

        if (near_zero(extract_matrix_element(exc, 1))) then
            pgen = 0.0_dp
            exc = 0_n_int
        end if

        ! store the most recent excitation information
        global_excitInfo = excitInfo

    end subroutine createStochasticExcitation_single



    subroutine createStochasticExcitation_double(ilut, nI, csf_i, excitation, pgen, excit_typ)
        ! calculate one possible double excitation and the corresponding
        ! probabilistic weight. and hamilton matrix element for a given CSF
        integer(n_int), intent(in) :: ilut(0:nifguga)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: nI(nel)
        integer(n_int), intent(out) :: excitation(0:nifguga)
        real(dp), intent(out) :: pgen
        integer, intent(out) :: excit_typ(2)
        character(*), parameter :: this_routine = "createStochasticExcitation_double"

        type(ExcitationInformation_t) :: excitInfo
        integer(n_int), allocatable :: excitations(:, :)
        integer(n_int) :: ilutI(0:niftot), ilutJ(0:niftot)
        real(dp) :: orb_pgen, branch_pgen
        HElement_t(dp) :: mat_ele
        type(WeightObj_t) :: weights

        logical :: compFlag
        real(dp) :: posSwitches(nSpatOrbs), negSwitches(nSpatOrbs)

        ASSERT(isProperCSF_ilut(ilut))

        ! do that in the calling interface funcition already..
        ! or maybe even in the FciMCPar to use the same b and occvector
        ! for and occupied determinant/CSF

        call pickOrbitals_double(ilut, nI, csf_i, excitInfo, orb_pgen)

        ! check if orbitals were correctly picked
        if ( .not. excitInfo%valid ) then
            excitation = 0_n_int
            pgen = 0.0_dp
            return
        end if

        ! do i need the checkcomp flag here?? how often does it happen that
        ! i create a wrong excitation information? can i avoid to create
        ! a wrong excitation information and thus not use checkComp here?
        ! profiler tells me this takes a lot of time..
        ! and if i have to do it, i could get rid of all the other uncessary
        ! call of calcRemainingSwitches ..
        !todo: since the weights also get initialized within the
        !checkcompatibility function, i should maybe output it also, similar to
        !the posSwitches and negSwitches quantitites..
        ! or just dont do a compatibility check, since it will get aborted
        ! anyway if it does not work in the excitation generation..

        call checkCompatibility(&
                    csf_i, excitInfo, compFlag, posSwitches, negSwitches, weights)

        if (.not. compFlag) then
            excitation = 0
            pgen = 0.0_dp
            return
        end if

        if (t_guga_back_spawn) then
            if (increase_ex_levl(csf_i, excitInfo) .and. .not. is_init_guga) then

                if (t_guga_back_spawn_trunc) then
                    pgen = 0.0_dp
                    excitation = 0_n_int
                    return
                end if

                call create_crude_guga_double(ilut, nI, csf_i, excitation, branch_pgen, excitInfo)

                pgen = orb_pgen * branch_pgen

                return
            end if
        end if

        if (tgen_guga_crude .and. .not. tgen_guga_mixed) then
            ! do the crude approximation here, where i do not switch
            ! in the excitation range, but only at the picked electrons
            ! and orbitals
            ! this includes the change, that I always switch at mixed
            ! start and ends too!

            call create_crude_double(ilut, excitation, branch_pgen, excitInfo)

            if (near_zero(branch_pgen)) then
                excitation = 0_n_int
                pgen = 0.0_dp
                return
            end if

            call convert_ilut_toNECI(ilut, ilutI)
            call convert_ilut_toNECI(excitation, ilutJ)

            call calc_guga_matrix_element(ilutI, csf_i, ilutJ, CSF_Info_t(ilutJ), excitInfo, mat_ele, .true.)

            if (near_zero(mat_ele)) then
                excitation = 0
                pgen = 0.0_dp

                return
            end if

            call encode_matrix_element(excitation, 0.0_dp, 2)
            call encode_matrix_element(excitation, mat_ele, 1)

            pgen = orb_pgen * branch_pgen

            return
        end if

        ! depending on the excitation chosen -> call specific stochastic
        ! excitation calculators. similar to the exact calculation
        ! only certain cases, where the excitation really can not be
        ! described by a single excitation, should arrive here.
        select case (excitInfo%typ)

        case (excit_type%single_overlap_L_to_R)
            ! single overlap lowering into raising
            ! similar to a single excitation except the (predetermined)
            ! single overlap site.
            call calcSingleOverlapMixedStochastic(ilut, csf_i, excitInfo, excitation, &
                                                  branch_pgen, posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%single_overlap_R_to_L) ! single overlap raising into lowering
            ! similar to a single excitation except the (predetermined)
            ! single overlap site.
            !todo: mention on the weight input here: in some routines below,
            ! the weights get reinitialized in the different sectors! so be
            ! careful to just input the weights everywhere, and also check the
            ! checkCompatibility function, if the weights get reinitialized
            ! there correctly!
            call calcSingleOverlapMixedStochastic(ilut, csf_i, excitInfo, excitation, &
                                                  branch_pgen, posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%double_lowering) ! normal double two lowering
            call calcDoubleLoweringStochastic(ilut, csf_i, excitInfo, excitation, branch_pgen, &
                                              posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%double_raising) ! normal double two raising
            call calcDoubleRaisingStochastic(ilut, csf_i, excitInfo, excitation, branch_pgen, &
                                             posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%double_L_to_R_to_L) ! lowergin into raising into lowering
            ! should be able to use the same general function as above to
            ! calculate the excitation, but the matrix element calculation
            ! should be a little bit different... maybe additional input needed
            call calcDoubleL2R2L_stochastic(ilut, csf_i, excitInfo, excitation, branch_pgen, &
                                            posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%double_R_to_L_to_R) ! raising into lowering into raising
            ! dito about matrix elements as above...
            call calcDoubleR2L2R_stochastic(ilut, csf_i, excitInfo, excitation, branch_pgen, &
                                            posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%double_L_to_R) ! lowering into raising double
            call calcDoubleL2R_stochastic(ilut, csf_i, excitInfo, excitation, branch_pgen, &
                                          posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%double_R_to_L) ! raising into lowering double
            call calcDoubleR2L_stochastic(ilut, csf_i, excitInfo, excitation, branch_pgen, &
                                          posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%fullstop_lowering) ! full stop 2 lowering
            ! again the double overlap part is easy to deal with, since its
            ! only the deltaB=0 branch
            call calcFullstopLoweringStochastic(ilut, csf_i, excitInfo, excitation, &
                                                branch_pgen, posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%fullstop_raising) ! full-stop 2 raising
            ! again only deltaB = 0 branch in DE overlap region
            call calcFullstopRaisingStochastic(ilut, csf_i, excitInfo, excitation, &
                                               branch_pgen, posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%fullstop_L_to_R) ! full-stop lowering into raising

            if (t_crude_exchange .or. (t_crude_exchange_noninits .and. (.not. is_init_guga))) then
                ! in case of "crude" exchange excitation perform a
                ! determinant-like excitation without spin-recouplings
                ! but only for non-inits.. so this information has to
                ! be passed in here!
                call perform_crude_excitation(ilut, csf_i, excitInfo, excitation, compFlag)

                ! in this case the pgen is just the orbital pgen, as only
                ! on CSF can be created from it..
                ! but be sure if the excitation is then actually possible
                if (.not. compFlag) then
                    excitation = 0_n_int
                    pgen = 0.0_dp
                    return
                else
                    pgen = orb_pgen
                end if

            else
                call calcFullStopL2R_stochastic(ilut, csf_i, excitInfo, excitation, &
                                                branch_pgen, posSwitches, negSwitches, weights)

                if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                    pgen = branch_pgen * orb_pgen
                else
                    pgen = branch_pgen
                end if

            end if

        case (excit_type%fullstop_R_to_L) ! full-stop raising into lowering

            if (t_crude_exchange .or. (t_crude_exchange_noninits .and. (.not. is_init_guga))) then
                call perform_crude_excitation(ilut, csf_i, excitInfo, excitation, compFlag)

                ! in this case the pgen is just the orbital pgen, as only
                ! on CSF can be created from it..
                ! but be sure if the excitation is then actually possible
                if (.not. compFlag) then
                    excitation = 0_n_int
                    pgen = 0.0_dp
                    return
                else
                    pgen = orb_pgen
                end if

            else

                call calcFullStopR2L_stochastic(ilut, csf_i, excitInfo, excitation, &
                                                branch_pgen, posSwitches, negSwitches, weights)

                if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                    pgen = branch_pgen * orb_pgen
                else
                    pgen = branch_pgen
                end if

            end if

        case (excit_type%fullstart_lowering) ! full-start 2 lowering
            ! again only deltaB = 0 branch in DE overlap region
            call calcFullStartLoweringStochastic(ilut, csf_i, excitInfo, excitation, &
                                                 branch_pgen, posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%fullstart_raising) ! full-start 2 raising
            ! the double overlap part is again really easy here, since only
            ! the deltaB=0 branch is non-zero -> and the second part can be
            ! treated as a single excitation
            call calcFullStartRaisingStochastic(ilut, csf_i, excitInfo, excitation, &
                                                branch_pgen, posSwitches, negSwitches, weights)

            pgen = orb_pgen * branch_pgen

        case (excit_type%fullStart_L_to_R) ! full-start lowering into raising

            if (t_crude_exchange .or. (t_crude_exchange_noninits .and. (.not. is_init_guga))) then
                call perform_crude_excitation(ilut, csf_i, excitInfo, excitation, compFlag)

                ! in this case the pgen is just the orbital pgen, as only
                ! on CSF can be created from it..
                ! but be sure if the excitation is then actually possible
                if (.not. compFlag) then
                    excitation = 0_n_int
                    pgen = 0.0_dp
                    return
                else
                    pgen = orb_pgen
                end if

            else

                call calcFullStartL2R_stochastic(ilut, csf_i, excitInfo, excitation, &
                                                 branch_pgen, posSwitches, negSwitches, weights)

                if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                    pgen = branch_pgen * orb_pgen
                else
                    pgen = branch_pgen
                end if

            end if

        case (excit_type%fullstart_R_to_L) ! full-start raising into lowering

            if (t_crude_exchange .or. (t_crude_exchange_noninits .and. (.not. is_init_guga))) then
                call perform_crude_excitation(ilut, csf_i, excitInfo, excitation, compFlag)

                ! in this case the pgen is just the orbital pgen, as only
                ! on CSF can be created from it..
                ! but be sure if the excitation is then actually possible
                if (.not. compFlag) then
                    excitation = 0_n_int
                    pgen = 0.0_dp
                    return
                else
                    pgen = orb_pgen
                end if

            else

                call calcFullStartR2L_stochastic(ilut, csf_i, excitInfo, excitation, &
                                                 branch_pgen, posSwitches, negSwitches, weights)

                if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                    pgen = branch_pgen * orb_pgen
                else
                    pgen = branch_pgen
                end if

            end if

            ! start, case by case how they appear in the documentary
        case (excit_type%fullstart_stop_alike) ! full-start into full-stop alike
            ! since there is only the deltaB = 0 branch with non-zero weight
            ! only 1 excitation is possible, which is calculated by the
            ! exact version

            ! check matrix element before calculating anything
            if (near_zero(get_umat_el(excitInfo%i, excitInfo%i, excitInfo%j, excitInfo%j))) then
                excitation = 0_n_int
                pgen = 0.0_dp
                return
            end if

            call calcFullStartFullStopAlike(ilut, csf_i, excitInfo, excitations)

            excitation = excitations(:, 1)

            ! have to encode the umat/2 matrix element here to keep this
            ! above function compatible with the exact routines
            call update_matrix_element(excitation, get_umat_el(excitInfo%i, &
                                                               excitInfo%i, excitInfo%j, excitInfo%j) / 2.0_dp, 1)

            ! probWeight is just the index choosing probability
            ! multiply further down
            pgen = orb_pgen
            ! can a full-start-full stop alike excitation have zero
            ! matrix element?
            ! yes since D(b=1,0) and d(b=0,1) can be zero

            ! to use t_trunc_guga_pgen
            branch_pgen = 1.0_dp

            deallocate(excitations)

        case (excit_type%fullstart_stop_mixed) ! full-start into full-stop mixed
            ! here it is garuanteed that its a open orbital at the start and
            ! end to allow for an non-diagonal excitation.
            ! but somehow i have to ensure, that the deltaB=0 path is left
            ! at somepoint... -> chances are slim, but there.. maybe can bias
            ! for that.. since i know how many switche possibilities are
            ! left and could include that somehow... todo
            ! on another note... since i know, that i have to switch at some
            ! point, i could totally ignore the x0-matrix element since it
            ! will be zero in the event of switching...
            ! i also need this behavior in the case of mixed full starts and
            ! full stops... -> so maybe write a specific double update function
            ! for that ...
            ! would have to include some alreadySwitched flag in the excitation
            ! to see and save if a switch already happened to enforce if
            ! necessary..
            ! or use the x0 matrix element as a kind of flag...
            ! since if its 0 it means a switch happended at some point, but
            ! thats seems a bit inefficient.

            if (t_crude_exchange .or. (t_crude_exchange_noninits .and. (.not. is_init_guga))) then
                call perform_crude_excitation(ilut, csf_i, excitInfo, excitation, compFlag)

                ! in this case the pgen is just the orbital pgen, as only
                ! on CSF can be created from it..
                ! but be sure if the excitation is then actually possible
                if (.not. compFlag) then
                    excitation = 0_n_int
                    pgen = 0.0_dp
                    return
                else
                    pgen = orb_pgen
                end if

            else

                call calcFullStartFullStopMixedStochastic(ilut, csf_i, excitInfo, &
                                                          excitation, branch_pgen, posSwitches, negSwitches, weights)

                if (t_approx_exchange .or. (t_approx_exchange_noninits .and. (.not. is_init_guga))) then
                    pgen = branch_pgen * orb_pgen
                else
                    pgen = branch_pgen
                end if

            end if

            ! random notes:

            ! those are the "easy" ones..., which dont actually have something
            ! to do in the DE overlap region.. well except the full-start into
            ! full-stop mixed excitation, where we actually HAVE to do a switch or
            ! else we get a already dealt with single excitation
            ! UPDATE: see the picking of the full-start or stop orbital as picking
            ! the first or last, necessary, stepvector switch. since it has to be
            ! switched somewhere anyway
            ! these are all the possibilities for (ii,jj), (ii,jk) index picking.

            ! a thinking still is to maybe combine all of them into one orbital picker
            ! and allow the second picked orbital to be the same as the first , so
            ! also this would account for the correct probability to pick to similar
            ! ones.. (a little bit of renormalization is needed, since the order
            ! how to pick the orbitals shouldnt matter and thus has to be accounted
            ! for. after here 4 differeing (ij,kl) indices were picked:
            ! that should be about all...

        end select

        ! what if probWeight is 0 for some reason? shouldnt be..
        ! yes it could be since i indicate zero-values excitations in this way

        if (t_trunc_guga_pgen .or. (t_trunc_guga_pgen_noninits .and. .not. is_init_guga)) then
            if (branch_pgen < trunc_guga_pgen) then
                pgen = 0.0_dp
                excitation = 0_n_int
                return
            end if
        end if

        ! check if for some reason the matrix element of the excitation is 0
        if (near_zero(extract_h_element(excitation))) then
            pgen = 0.0_dp
            excitation = 0

        else
            ! also store information of type of excitation for automated tau-search
            ! for the non-weighted guga-excitation-generator
            select case (excitInfo%excitLvl)
            case (0)
                ! (ii,jj) RR/LL excitation
                excit_typ(1) = 2
                excit_typ(2) = 1

            case (1)
                ! (ii,jj) RL excitation
                excit_typ(1) = 2
                excit_typ(2) = 0

            case (2)
                ! (ii,jk) RR/LL excitation
                excit_typ(1) = 3
                excit_typ(2) = 1

            case (3)
                ! (ii,jk) RL excitation
                excit_typ(1) = 3
                excit_typ(2) = 0

            case (4)
                ! (ij,kl) excitation
                excit_typ(1) = 4
                excit_typ(2) = 0

            case default
                call print_excitInfo(excitInfo)
                call stop_all(this_routine, "wrong excit level info!")

            end select
        end if

        select case(excitInfo%typ)
        case (excit_type%single_overlap_L_to_R, &
              excit_type%single_overlap_R_to_L, &
              excit_type%double_lowering, &
              excit_type%double_raising, &
              excit_type%double_L_to_R_to_L, &
              excit_type%double_R_to_L_to_R, &
              excit_type%double_L_to_R, &
              excit_type%double_R_to_L, &
              excit_type%fullstop_raising, &
              excit_type%fullstop_lowering, &
              excit_type%fullstart_lowering, &
              excit_type%fullstart_raising, &
              excit_type%fullstart_stop_alike)

            ! in the other cases global_excitInfo gets assigned before it
            ! gets changed
            global_excitInfo = excitInfo
        end select

    end subroutine createStochasticExcitation_double



    function test_increase_on_loc(loc_elec, loc_orb, ic) result(flag)
        ! test if the excitation increases the excit-lvl based on the
        ! restriction and type of excitation
        integer, intent(in) :: loc_elec, loc_orb, ic
        logical :: flag

        if (ic == 1) then
            ! now the global restriction of n_guga_back_spawn_lvl comes into
            ! play
            select case (n_guga_back_spawn_lvl)
            case (-2)
                ! we want to treat double excitation decreasing the
                ! excit-lvl by 2 fully ..
                ! so single excitations from non-initiators (make this
                ! default!) are always subjected to the approximation
                flag = .true.

            case (-1)
                ! if this excitation decreases the excit-lvl by 1 we
                ! treat it fully
                ! for this to happen the electron must be in the
                ! virtual space of the reference and the orbital must be
                ! in the occupied space of the reference
                if (loc_elec == 0 .and. loc_orb == 0) then
                    flag = .false.
                else
                    flag = .true.
                end if
            case (0)
                ! here we want to only restrict excitation increasing the
                ! excitation lvl with the approximation
                ! this happens if the electron is in the occupied space of
                ! the reference and the orbital in the virtual space
                if (loc_elec == 2 .and. loc_orb == 2) then
                    flag = .true.
                else
                    flag = .false.
                end if
            case (1)
                ! in this case we treat all single excitation fully
                flag = .false.

            end select
        else if (ic == 2) then
            ! maybe i need specific restriction for different types of
            ! GUGA excitations.. figure that out!

            select case (n_guga_back_spawn_lvl)
            case (-2)
                ! only doubles reducing ex-lvl by two get treated fully
                if (loc_elec == 0 .and. loc_orb == 0) then
                    flag = .false.
                else
                    ! everything else gets treated fully
                    flag = .true.
                end if

            case (-1)
                ! also doubles which increase the excit-lvl by 1
                ! get treated fully ..
                ! how does this happen?
                ! at least one electron must hope from the reference
                ! virtuals to the occupied reference space..
                if (loc_elec == 0) then
                    ! both electrons are in the virtual, so atleast
                    ! one orbital must be in the reference
                    if (loc_orb < 2) then
                        flag = .false.
                    else
                        flag = .true.
                    end if
                else if (loc_elec == 1) then
                    ! one electron in occupied and one in virtual
                    ! then both holes must be in the occupied to decrease
                    if (loc_orb == 0) then
                        flag = .false.
                    else
                        flag = .true.
                    end if
                else
                    ! if both electrons are in the occupied space
                    ! it is not possible
                    flag = .true.
                end if

            case (0)
                ! here we also treat excitation leaving the excit-lvl
                ! the same fully..
                if (loc_elec == 0) then
                    ! if both electron are in the virtual space we can not
                    ! increase the excit-lvl
                    flag = .false.

                else if (loc_elec == 1) then
                    ! if one of the electrons is in the occupied space
                    ! atleast one orbital must also be in the virtual space
                    if (loc_orb == 2) then
                        flag = .true.
                    else
                        flag = .false.
                    end if

                else if (loc_elec == 2) then
                    ! if both electrons are in occupied space
                    ! both orbital must also be in the occupied space
                    if (loc_orb == 0) then
                        flag = .false.
                    else
                        flag = .true.
                    end if
                end if

            case (1)
                ! here we also want to treat excitation increasing the
                ! excitation lvl by up to 1 fully

                ! if both electrons are in the virtual space
                ! we do not increase the excit-lvl
                flag = .false.

                ! if only one electron is in the occupied space
                ! we at most increase it by 1, which is fine here
                flag = .false.

                ! if both electrons are in the virtual space
                ! we increase by more than 1 only if both orbs are in
                ! the virtual space

                if (loc_elec == 2 .and. loc_orb == 2) then
                    flag = .true.
                else
                    flag = .false.
                end if
            end select
        end if

    end function test_increase_on_loc


    function increase_ex_levl(csf_i, excitInfo) result(flag)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(in) :: excitInfo
        logical :: flag
        character(*), parameter :: this_routine = "increase_ex_levl"

        integer :: elec_1, elec_2, orb_1, orb_2, orbs(2), elecs(2)
        integer :: d_elec, s_elec, d_orb, s_orb
        integer :: loc_elec, loc_orb

        flag = .true.

        ! first i need to check the location of the picked electrons.
        ! which also has to be adapted for chosen spatial orbitals
        if (excitInfo%typ == excit_type%single) then
            ! single excitation
            elec_1 = excitInfo%j
            orb_1 = excitInfo%i

            ! be general here and maybe a bit too cautious and check if
            ! any spin-orbtial in the reference is occupied
            loc_elec = check_electron_location_spatial([elec_1, 0], 1, 1)

            ! now we also want to check if the orbitals are in the
            ! virtual space of the reference det.
            loc_orb = check_orbital_location_spatial([orb_1, 0], 1, 1)

            flag = test_increase_on_loc(loc_elec, loc_orb, 1)

        else
            ! double excitations
            elec_1 = excitInfo%j
            elec_2 = excitInfo%l

            orb_1 = excitInfo%i
            orb_2 = excitInfo%k

            ! there is now a difference depending on some type of
            ! excitations
            select case (excitInfo%typ)
            case (excit_type%single_overlap_L_to_R, &
                  excit_type%fullstop_lowering, &
                  excit_type%fullstart_raising)

                ! here i know the spatial orbital indices are the same
                ASSERT(orb_1 == orb_2)
                ASSERT(csf_i%stepvector(orb_1) == 0)

                ! here check for spin-orbital as we know the occupation
                orbs = [2 * orb_1, 2 * orb_1 - 1]
                loc_orb = check_orbital_location(orbs, 2, 1)

                ! the electrons need to be specified generally
                ! but are there any spin-restrictions in this case?
                ! due to the possible spin-recouplings in CSFs maybe not..
                ! use a testing function for spatial orbitals
                loc_elec = check_electron_location_spatial([elec_1, elec_2], 2, 1)

                flag = test_increase_on_loc(loc_elec, loc_orb, 2)

            case (excit_type%single_overlap_R_to_L, &
                  excit_type%fullstop_raising, &
                  excit_type%fullstart_lowering)

                ! here i know both spatial electon indices are the same
                ASSERT(elec_1 == elec_2)
                ASSERT(csf_i%stepvector(elec_1) == 3)

                elecs = [2 * elec_1, 2 * elec_1 - 1]

                loc_elec = check_electron_location(elecs, 2, 1)

                loc_orb = check_orbital_location_spatial([orb_1, orb_2], 2, 1)

                flag = test_increase_on_loc(loc_elec, loc_orb, 2)

            case (excit_type%fullstop_L_to_R, &
                  excit_type%fullstop_R_to_L, &
                  excit_type%fullStart_L_to_R, &
                  excit_type%fullstart_R_to_L)

                ! here i know one electron and one hole index are the same
                ASSERT(elec_1 /= elec_2)
                ASSERT(orb_1 /= orb_2)

                ! the occupation in the overlap index does not change..
                ! so we could treat the differing indices as a single
                ! excitation or?
                if (elec_1 == orb_1) then

                    s_elec = elec_1
                    s_orb = orb_1

                    d_elec = elec_2
                    d_orb = orb_2

                else if (elec_1 == orb_2) then

                    s_elec = elec_1
                    s_orb = orb_1

                    d_elec = elec_2
                    d_orb = orb_1

                else if (elec_2 == orb_1) then

                    s_elec = elec_2
                    s_orb = orb_1

                    d_elec = elec_1
                    d_orb = orb_2

                else if (elec_2 == orb_2) then

                    s_elec = elec_2
                    s_orb = orb_2

                    d_elec = elec_1
                    d_orb = orb_1

                end if

                loc_elec = check_electron_location_spatial([d_elec, 0], 1, 1)
                loc_orb = check_orbital_location_spatial([d_orb, 0], 1, 1)

                flag = test_increase_on_loc(loc_elec, loc_orb, 1)

            case (excit_type%fullstart_stop_mixed)
                ! here i do not change the 'orbital excitation level'
                if (n_guga_back_spawn_lvl < 0) then
                    flag = .true.

                else
                    flag = .false.
                end if

            case default
                ! the general 4-index excitations..
                loc_elec = check_electron_location_spatial([elec_1, elec_2], 2, 1)
                loc_orb = check_orbital_location_spatial([orb_1, orb_2], 2, 1)

                flag = test_increase_on_loc(loc_elec, loc_orb, 2)

            end select
        end if

    end function increase_ex_levl




end module guga_main
