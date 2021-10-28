#include "macros.h"

module tJ_model

    use SystemData, only: bhub, nel, nbasis, G1, lattice_type, length_x, &
                          length_y, length_z, nbasis, t_open_bc_x, t_open_bc_y, &
                          t_open_bc_z, ecore, tHPHF, tHub, tReal, t_tJ_model, &
                          t_heisenberg_model, t_new_real_space_hubbard, exchange_j, &
                          t_trans_corr, trans_corr_param, &
                          t_trans_corr_2body, trans_corr_param_2body, &
                          nSpatOrbs, t_bipartite_order

    use constants, only: dp, n_int, EPS, bits_n_int, maxExcit

    use real_space_hubbard, only: lat_tau_factor, t_start_neel_state, &
                                  check_real_space_hubbard_input, init_tmat, &
                                  init_spin_free_tmat

    use procedure_pointers, only: get_umat_el, generate_excitation

    use FciMCData, only: tsearchtau, tsearchtauoption, excit_gen_store_type, &
                         pSingles, pDoubles

    use CalcData, only: t_hist_tau_search_option, t_hist_tau_search, tau

    use bit_rep_data, only: NIfTot, nifguga, nifd, GugaBits

    use umatcache, only: gtid

    use util_mod, only: operator(.div.), near_zero, get_free_unit

    use util_mod_numerical, only: binary_search_first_ge

    use OneEInts, only: GetTMatEl, tmat2d

    use lattice_mod, only: lattice, lat, determine_optimal_time_step, &
                           get_helement_lattice_ex_mat, get_helement_lattice_general

    use DetBitOps, only: FindBitExcitLevel, EncodeBitDet, GetBitExcitation

    use double_occ_mod, only: count_double_orbs

    use FciMCData, only: ilutref

    use bit_reps, only: decode_bit_det

    use lattice_models_utils, only: get_spin_opp_neighbors, make_ilutJ, &
                                    find_elec_in_ni, get_spin_density_neighbors, &
                                    get_occ_neighbors

    use MPI_wrapper, only: iProcIndex

    use get_excit, only: make_single, make_double

    use dsfmt_interface, only: genrand_real2_dsfmt

    use guga_excitations, only: generate_excitation_guga, &
                                assign_excitInfo_values_double, assign_excitInfo_values_single

    use guga_matrixElements, only: calc_guga_matrix_element

    use guga_data, only: ExcitationInformation_t, excit_type, gen_type

    use guga_bitRepOps, only: count_alpha_orbs_ij, count_beta_orbs_ij, &
                              write_det_guga, CSF_Info_t

    implicit none

    real(dp), allocatable :: exchange_matrix(:, :)
    real(dp), allocatable :: spin_free_exchange(:, :)

    interface get_helement_tJ
        module procedure get_helement_tJ_ex_mat
        module procedure get_helement_tJ_general
    end interface get_helement_tJ

    interface get_helement_heisenberg
        module procedure get_helement_heisenberg_ex_mat
        module procedure get_helement_heisenberg_general
    end interface get_helement_heisenberg

contains

    subroutine init_guga_tJ_model
        ! write a specific subroutine to init the spin-free tJ model
        character(*), parameter :: this_routine = "init_guga_tJ_model"

        real(dp) :: tau_opt

        root_print "initializing the spin-free t-J model, with parameter: "
        root_print "t:", bhub
        root_print "J:", exchange_j

        tHub = .false.
        treal = .false.
        t_new_real_space_hubbard = .true.

        call check_real_space_hubbard_input()

        nSpatOrbs = nBasis / 2

        pSingles = 0.1_dp
        pDoubles = 1.0_dp - pSingles

        get_umat_el => get_umat_heisenberg_spin_free

        if (trim(adjustl(lattice_type)) == 'read') then
            ! then i have to construct tmat first
            call stop_all(this_routine, "starting from fcidump not yet implemented!")
            ! and then construct the lattice
            lat => lattice(lattice_type, length_x, length_y, length_z,.not. t_open_bc_x, &
                           .not. t_open_bc_y,.not. t_open_bc_z)
        else
            ! otherwise i have to do it the other way around
            lat => lattice(lattice_type, length_x, length_y, length_z,.not. t_open_bc_x, &
                       .not. t_open_bc_y,.not. t_open_bc_z, t_bipartite_order = t_bipartite_order)
            ! if nbasis was not yet provided:
            if (nbasis <= 0) then
                nbasis = 2 * lat%get_nsites()
            end if

            call init_tmat(lat)
            call init_spin_free_tmat(lat)
            call setup_exchange_matrix(lat)
            call setup_spin_free_exchange(lat)

        end if

        if (nel >= nbasis / 2) then
            call stop_all(this_routine, &
                          " too many electrons for the tJ model! nel >= nbasis/2")
        end if

        ! and also check the double occupancy in the starting det,
        ! no double occupancy allowed!
        if (count_double_orbs(ilutRef) > 0) then
            call stop_all(this_routine, &
                          "incorrect starting state for tJ model: there is an doubly occupied site!")
        end if

        ! i guess i have to setup G1 also.. argh.. i hate this!
        allocate(G1(nbasis))
        G1(1:nbasis - 1:2)%ms = -1
        G1(2:nbasis:2)%ms = 1

        ! Ecore should default to 0, but be sure anyway!
        ecore = 0.0_dp

        tau_opt = determine_optimal_time_step()
        if (tau < EPS) then
            root_print "setting time-step to optimally determined time-step: ", tau_opt
            root_print "times: ", lat_tau_factor
            tau = lat_tau_factor * tau_opt

        else
            root_print "optimal time-step would be: ", tau_opt
            root_print "but tau specified in input!"
        end if

        generate_excitation => generate_excitation_guga

    end subroutine init_guga_tJ_model

    subroutine init_tJ_model
        character(*), parameter :: this_routine = "init_tJ_model"
        real(dp) :: tau_opt

        root_print "initializing tJ-model with parameters: "
        root_print "t: ", bhub
        root_print "J: ", exchange_j

        ! after having used the tHub and treal parameters set them to false
        ! now
        thub = .false.
        treal = .false.

        t_new_real_space_hubbard = .false.

        ! reuse real-space-hubbard stugg
        call check_real_space_hubbard_input()

        get_umat_el => get_umat_el_heisenberg

        if (trim(adjustl(lattice_type)) == 'read') then
            ! then i have to construct tmat first
            call stop_all(this_routine, "starting from fcidump not yet implemented!")
            ! and then construct the lattice
            lat => lattice(lattice_type, length_x, length_y, length_z,.not. t_open_bc_x, &
                           .not. t_open_bc_y,.not. t_open_bc_z)
        else
            ! otherwise i have to do it the other way around
            lat => lattice(lattice_type, length_x, length_y, length_z,.not. t_open_bc_x, &
                       .not. t_open_bc_y,.not. t_open_bc_z, t_bipartite_order = t_bipartite_order)

            ! if nbasis was not yet provided:
            if (nbasis <= 0) then
                nbasis = 2 * lat%get_nsites()
            end if

            call init_tmat(lat)
            call setup_exchange_matrix(lat)

        end if

        if (nel >= nbasis / 2) then
            call stop_all(this_routine, &
                          " too many electrons for the tJ model! nel >= nbasis/2")
        end if

        ! and also check the double occupancy in the starting det,
        ! no double occupancy allowed!
        if (count_double_orbs(ilutRef) > 0) then
            call stop_all(this_routine, &
                          "incorrect starting state for tJ model: there is an doubly occupied site!")
        end if

        ! i guess i have to setup G1 also.. argh.. i hate this!
        allocate(G1(nbasis))
        G1(1:nbasis - 1:2)%ms = -1
        G1(2:nbasis:2)%ms = 1

        ! Ecore should default to 0, but be sure anyway!
        ecore = 0.0_dp

        ! and i have to calculate the optimal time-step for the hubbard models.
        ! where i need the connectivity of the lattice i guess?
        if (.not. tHPHF) then
            generate_excitation => gen_excit_tJ_model
        end if

        tau_opt = determine_optimal_time_step()
        if (tau < EPS) then
            root_print "setting time-step to optimally determined time-step: ", tau_opt
            root_print "times: ", lat_tau_factor
            tau = lat_tau_factor * tau_opt

        else
            root_print "optimal time-step would be: ", tau_opt
            root_print "but tau specified in input!"
        end if

        ! and i have to turn off the time-step search for the hubbard
        tsearchtau = .false.
        ! set tsearchtauoption to true to use the death-tau search option
        tsearchtauoption = .true.

        t_hist_tau_search = .false.
        t_hist_tau_search_option = .false.

        if (t_start_neel_state) then
            root_print "starting from the Neel state: "
            if (nel > nbasis / 2) then
                call stop_all(this_routine, &
                              "more than half-filling! does neel state make sense?")
            end if
        end if

    end subroutine init_tJ_model

    subroutine init_guga_heisenberg_model()
        character(*), parameter :: this_routine = "init_guga_heisenberg_model"
        real(dp) :: tau_opt

        root_print "initializing spin-free Heisenberg model, with "
        root_print "J: ", exchange_j

!         call init_guga()

        tHub = .false.
        treal = .false.
        t_new_real_space_hubbard = .false.

        call check_real_space_hubbard_input()

        pSingles = 0.0_dp
        pDoubles = 1.0_dp

        get_umat_el => get_umat_heisenberg_spin_free

        if (associated(tmat2d)) deallocate(tmat2d)
        allocate(tmat2d(nbasis, nbasis), source=h_cast(0.0_dp))

        if (trim(adjustl(lattice_type)) == 'read') then
            ! then i have to construct tmat first
            ! no need for tmat in the heisenberg model
            call stop_all(this_routine, "starting from fcidump not yet implemented!")
            ! and then construct the lattice
            lat => lattice(lattice_type, length_x, length_y, length_z,.not. t_open_bc_x, &
                           .not. t_open_bc_y,.not. t_open_bc_z)
        else
            lat => lattice(lattice_type, length_x, length_y, length_z,.not. t_open_bc_x, &
                       .not. t_open_bc_y,.not. t_open_bc_z, t_bipartite_order = t_bipartite_order)

            ! if nbasis was not yet provided:
            if (nbasis <= 0) then
                nbasis = 2 * lat%get_nsites()
            end if

            call setup_exchange_matrix(lat)
            call setup_spin_free_exchange(lat)

        end if

        if (nel /= nbasis / 2) then
            call stop_all(this_routine, &
                          "heisenberg model need half filling nel == nbasis/2")
        end if

        if (count_double_orbs(ilutref) > 0) then
            call stop_all(this_routine, &
                          " no double occupancies allowed in the heisenberg model")
        end if

        ! i guess i have to setup G1 also.. argh.. i hate this!
        allocate(G1(nbasis))
        G1(1:nbasis - 1:2)%ms = -1
        G1(2:nbasis:2)%ms = 1

        ! Ecore should default to 0, but be sure anyway!
        ecore = 0.0_dp

        ! TODO: maybe I need the tau-search for the GUGA..
        tau_opt = determine_optimal_time_step()
        if (tau < EPS) then
            root_print "setting time-step to optimally determined time-step: ", tau_opt
            root_print "times: ", lat_tau_factor
            tau = lat_tau_factor * tau_opt

        else
            root_print "optimal time-step would be: ", tau_opt
            root_print "but tau specified in input!"
        end if

        ! and i have to turn off the time-step search for the hubbard
        tsearchtau = .false.
        ! set tsearchtauoption to true to use the death-tau search option
        tsearchtauoption = .true.

        t_hist_tau_search = .false.
        t_hist_tau_search_option = .false.

        generate_excitation => generate_excitation_guga

    end subroutine init_guga_heisenberg_model

    subroutine init_heisenberg_model
        character(*), parameter :: this_routine = "init_heisenberg_model"
        real(dp) :: tau_opt

        root_print "initialising Heisenberg model with "
        root_print "J: ", exchange_j

        thub = .false.
        treal = .false.

        call check_real_space_hubbard_input()

        if (t_trans_corr_2body) then
            if (t_trans_corr) then
                call stop_all(this_routine, &
                              "1-body transcorrelation not allowed in the heisenberg model!")
            end if
        end if

        get_umat_el => get_umat_el_heisenberg

        if (trim(adjustl(lattice_type)) == 'read') then
            ! then i have to construct tmat first
            ! no need for tmat in the heisenberg model
            call stop_all(this_routine, "starting from fcidump not yet implemented!")
            ! and then construct the lattice
            lat => lattice(lattice_type, length_x, length_y, length_z,.not. t_open_bc_x, &
                           .not. t_open_bc_y,.not. t_open_bc_z)
        else
            lat => lattice(lattice_type, length_x, length_y, length_z,.not. t_open_bc_x, &
                       .not. t_open_bc_y,.not. t_open_bc_z, t_bipartite_order = t_bipartite_order)

            ! if nbasis was not yet provided:
            if (nbasis <= 0) then
                nbasis = 2 * lat%get_nsites()
            end if

            call setup_exchange_matrix(lat)

        end if

        if (nel /= nbasis / 2) then
            call stop_all(this_routine, &
                          "heisenberg model need half filling nel == nbasis/2")
        end if

        if (count_double_orbs(ilutref) > 0) then
            call stop_all(this_routine, &
                          " no double occupancies allowed in the heisenberg model")
        end if

        ! i guess i have to setup G1 also.. argh.. i hate this!
        allocate(G1(nbasis))
        G1(1:nbasis - 1:2)%ms = -1
        G1(2:nbasis:2)%ms = 1

        ! Ecore should default to 0, but be sure anyway!
        ecore = 0.0_dp

        ! and i have to calculate the optimal time-step for the hubbard models.
        ! where i need the connectivity of the lattice i guess?
        if (.not. tHPHF) then
            generate_excitation => gen_excit_heisenberg_model
        end if

        tau_opt = determine_optimal_time_step()
        if (tau < EPS) then
            root_print "setting time-step to optimally determined time-step: ", tau_opt
            root_print "times: ", lat_tau_factor
            tau = lat_tau_factor * tau_opt

        else
            root_print "optimal time-step would be: ", tau_opt
            root_print "but tau specified in input!"
        end if

        ! and i have to turn off the time-step search for the hubbard
        tsearchtau = .false.
        ! set tsearchtauoption to true to use the death-tau search option
        tsearchtauoption = .true.

        t_hist_tau_search = .false.
        t_hist_tau_search_option = .false.

        if (t_start_neel_state) then
!             neel_state_ni = create_neel_state(ilut_neel)

            root_print "starting from the Neel state: "
            if (nel > nbasis / 2) then
                call stop_all(this_routine, &
                              "more than half-filling! does neel state make sense?")
            end if

        end if

    end subroutine init_heisenberg_model

    subroutine init_get_helement_tj_guga

        call init_tmat(lat)
        call init_spin_free_tmat(lat)
        call setup_exchange_matrix(lat)
        call setup_spin_free_exchange(lat)

        get_umat_el => get_umat_heisenberg_spin_free

    end subroutine init_get_helement_tj_guga

    subroutine init_get_helement_heisenberg_guga

        if (associated(tmat2d)) deallocate(tmat2d)
        allocate(tmat2d(nbasis, nbasis), source=h_cast(0.0_dp))

        call setup_exchange_matrix(lat)
        call setup_spin_free_exchange(lat)

        get_umat_el => get_umat_heisenberg_spin_free

    end subroutine init_get_helement_heisenberg_guga

    subroutine init_get_helement_tj
        get_helement_lattice_ex_mat => get_helement_tJ_ex_mat
        get_helement_lattice_general => get_helement_tJ_general
        call init_tmat(lat)
        call setup_exchange_matrix(lat)
    end subroutine init_get_helement_tj

    subroutine init_get_helement_heisenberg
        get_helement_lattice_ex_mat => get_helement_heisenberg_ex_mat
        get_helement_lattice_general => get_helement_heisenberg_general
        call setup_exchange_matrix(lat)
    end subroutine init_get_helement_heisenberg

    subroutine gen_excit_heisenberg_model(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                          ex, tParity, pGen, hel, store, run)
        ! the heisenberg excitation generator is only a small modification of
        ! the t-J excitation generator without single excitation hoppings,
        ! due to half-filling

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: run
#ifdef DEBUG_
        character(*), parameter :: this_routine = "gen_excit_heisenberg_model"
#endif
        integer :: elec, src, ind, tgt, src_opp, tgt_opp, elec_opp
        real(dp) :: p_elec, cum_sum, r, p_orb, cpt_opp
        integer, allocatable :: neighbors(:)
        real(dp), allocatable :: cum_arr(:)

        unused_var(exFlag)
        unused_var(store)
        if (present(run)) then
            unused_var(run)
        end if

        hel = h_cast(0.0_dp)

        ASSERT(associated(lat))
        ASSERT(nel == nbasis / 2)

        ic = 2

        ! still pick the first electron at random
        elec = 1 + int(genrand_real2_dsfmt() * nel)

        p_elec = 1.0_dp / real(nel, dp)

        src = nI(elec)

        ! in the heisenberg model i am sure that every site is occupied..
        ! so i could pick the spin-orbit neighbors and check if the spin-orbital
        ! is empty, this would indicate opposite neighboring spins..
        neighbors = lat%get_spinorb_neighbors(src)

        call create_cum_list_heisenberg(ilutI, src, neighbors, cum_arr, cum_sum)

        if (cum_sum < EPS) then
            nJ(1) = 0
            pgen = 0.0_dp
            return
        end if

        r = genrand_real2_dsfmt() * cum_sum

        ind = binary_search_first_ge(cum_arr, r)

        tgt = neighbors(ind)

        if (ind == 1) then
            p_orb = cum_arr(1) / cum_sum
        else
            p_orb = (cum_arr(ind) - cum_arr(ind - 1)) / cum_sum
        end if

        ! and then i have to check for the opposite order generation prob
        if (is_beta(src)) then
            src_opp = neighbors(ind) + 1
            tgt_opp = src + 1
        else
            src_opp = neighbors(ind) - 1
            tgt_opp = src - 1
        end if

        ASSERT(IsOcc(ilutI, src_opp))
        ASSERT(IsNotOcc(ilutI, tgt_opp))

        call create_cum_list_heisenberg(ilutI, src_opp, lat%get_spinorb_neighbors(src_opp), &
                                        cum_arr, cum_sum, tgt_opp, cpt_opp)

        pgen = p_elec * (p_orb + cpt_opp)

        elec_opp = find_elec_in_ni(nI, src_opp)

        call make_double(nI, nJ, elec, elec_opp, tgt, tgt_opp, ex, tParity)

        ilutJ = make_ilutJ(ilutI, ex, 2)

    end subroutine gen_excit_heisenberg_model

    subroutine gen_excit_tJ_model(nI, ilutI, nJ, ilutJ, exFlag, ic, &
                                  ex, tParity, pGen, hel, store, run)

        integer, intent(in) :: nI(nel), exFlag
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(out) :: nJ(nel), ic, ex(2, maxExcit)
        integer(n_int), intent(out) :: ilutJ(0:NifTot)
        real(dp), intent(out) :: pGen
        logical, intent(out) :: tParity
        HElement_t(dp), intent(out) :: hel
        type(excit_gen_store_type), intent(inout), target :: store
        integer, intent(in), optional :: run

        character(*), parameter :: this_routine = "gen_excit_tJ_model"

        integer :: elec, src, id, ind, elec_2, tgt_1, tgt_2
        real(dp) :: p_elec, p_orb, cum_sum, r, cum_sum_opp, cpt_opp
        integer, allocatable :: neighbors(:), ic_list(:), tmp_ic_list(:)
        real(dp), allocatable :: cum_arr(:), cum_arr_opp(:)

        ! the idea for the tJ excitation generator on a lattice is to
        ! still pick the electron at random and then check for its
        ! neighbors, if the neighboring site is empty, do a hop
        ! and if there is a electron of opposite spin do a spin-flip
        ! if it is occupied by an electron of same spin no excitation with
        ! this neighbor is possible

        unused_var(exFlag)
        unused_var(store)
        unused_var(run)

        hel = h_cast(0.0_dp)

        ! use the lattice type like in the real-space hubbard implementation
        ASSERT(associated(lat))
#ifdef WARNING_WORKAROUND_
        hel = 0.0_dp
        if (present(run)) then
            unused_var(run)
        end if
#endif
        unused_var(store)
        unused_var(exflag)

        elec = 1 + int(genrand_real2_dsfmt() * nel)

        p_elec = 1.0_dp / real(nel, dp)

        src = nI(elec)
        id = gtid(src)

        neighbors = lat%get_neighbors(id)

        call create_cum_list_tJ_model(ilutI, src, neighbors, cum_arr, cum_sum, &
                                      ic_list)

        if (cum_sum < EPS) then
            nJ(1) = 0
            pgen = 0.0_dp
            return
        end if

        r = genrand_real2_dsfmt() * cum_sum

        ind = binary_search_first_ge(cum_arr, r)

        ! i just realised that for spin-flip excitation we have to take the
        ! opposite order of orbital picking into account too..
        ! because the same excitation could have been picked with the spin
        ! opposite electron in the neighborhood of the first electron too..

        ic = ic_list(ind)

        ! for spin-flips i have to add a contribution down below!
        if (ind == 1) then
            p_orb = cum_arr(1) / cum_sum
        else
            p_orb = (cum_arr(ind) - cum_arr(ind - 1)) / cum_sum
        end if

        if (ic == 1) then
            if (is_beta(src)) then
                tgt_1 = 2 * neighbors(ind) - 1
            else
                tgt_1 = 2 * neighbors(ind)
            end if

            call make_single(nI, nJ, elec, tgt_1, ex, tParity)

            ilutJ = make_ilutJ(ilutI, ex, 1)

        else if (ic == 2) then
            ! here i have to recalc the contribution if i would have picked
            ! the electron in orbital spin_orb first
            ! but i made some assumptions about the order of the picked
            ! electrons and holes, which is not valid to do..

            if (is_beta(src)) then
                ! need to the the index of electron 2 in nI
                ! the second electron must be alpha
                elec_2 = find_elec_in_ni(nI, 2 * neighbors(ind))
                ! we need the orbital alpha of src
                ! and the beta of the second orbital
                tgt_1 = get_alpha(src)
                tgt_2 = 2 * neighbors(ind) - 1
            else
                ! v.v here
                elec_2 = find_elec_in_ni(nI, 2 * neighbors(ind) - 1)

                tgt_1 = get_beta(src)
                tgt_2 = 2 * neighbors(ind)

            end if

            ! the idea is to target the spin-orbital of the other electron!
            ! (the first in this case!
            call create_cum_list_tJ_model(ilutI, nI(elec_2), &
                lat%get_neighbors(gtid(nI(elec_2))), cum_arr_opp, cum_sum_opp, &
                tmp_ic_list, src, cpt_opp)

            p_orb = p_orb + cpt_opp

            call make_double(nI, nJ, elec, elec_2, tgt_1, tgt_2, ex, tParity)

            ilutJ = make_ilutJ(ilutI, ex, 2)

        else
            ! something went wrong..
            call stop_all(this_routine, &
                          "something went wrong ic > 2!")
        end if

        pgen = p_elec * p_orb

    end subroutine gen_excit_tJ_model

    subroutine create_cum_list_tJ_model(ilutI, src, neighbors, cum_arr, cum_sum, &
                                        ic_list, tgt, cpt)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: src, neighbors(:)
        real(dp), intent(out), allocatable :: cum_arr(:)
        real(dp), intent(out) :: cum_sum
        integer, intent(out), allocatable :: ic_list(:)
        integer, intent(in), optional :: tgt
        real(dp), intent(out), optional :: cpt
#ifdef DEBUG_
        character(*), parameter :: this_routine = "create_cum_list_tJ_model"
#endif
        integer :: i, nI(nel), temp_ex(2, maxExcit)
        integer, allocatable :: single_excits(:)
        integer, allocatable :: spin_flips(:)
        real(dp) :: elem
        logical :: t_single, t_flip, t_single_possible, t_flip_possible

        ASSERT(IsOcc(ilutI, src))
        ASSERT(allocated(exchange_matrix))

        allocate(cum_arr(size(neighbors)))
        allocate(ic_list(size(neighbors)))
        allocate(single_excits(size(neighbors)))
        allocate(spin_flips(size(neighbors)))
        cum_arr = 0
        cum_sum = 0.0_dp
        ic_list = 0

        call decode_bit_det(nI, ilutI)

        temp_ex(1, 1) = src

        if (is_beta(src)) then
            single_excits = 2 * neighbors - 1
            spin_flips = 2 * neighbors
            ! fill in the corresponding alpha orbital
            temp_ex(2, 1) = src + 1
        else
            single_excits = 2 * neighbors
            spin_flips = 2 * neighbors - 1

            ! fill in the beta orbital
            temp_ex(2, 1) = src - 1
        end if

        if (present(tgt)) then
            t_single = .false.
            t_flip = .false.

            ! find the probability of choosing orbital target
            if (is_beta(src) .eqv. is_beta(tgt)) then
                ! then it was definetly a single excitation
                t_single = .true.
            else
                t_flip = .true.
            end if

            ASSERT(present(cpt))
            cpt = 0.0_dp

            do i = 1, ubound(neighbors, 1)
                elem = 0.0_dp
                t_single_possible = .false.
                t_flip_possible = .false.

                if (IsNotOcc(ilutI, single_excits(i)) .and. &
                    IsNotOcc(ilutI, spin_flips(i))) then
                    ! just to be sure use the tmat, so both orbitals are
                    ! definetly connected
                    elem = abs(get_offdiag_helement_tJ(nI, [src, single_excits(i)], .false.))
!                     elem = abs(GetTMatEl(src, single_excits(i)))
                    t_single_possible = .true.

                else if (IsOcc(ilutI, spin_flips(i)) .and. &
                         IsNotOcc(ilutI, single_excits(i))) then

                    temp_ex(1, 2) = spin_flips(i)
                    temp_ex(2, 2) = single_excits(i)
                    elem = abs(get_offdiag_helement_heisenberg(nI, temp_ex, .false.))
!                      elem = abs(get_heisenberg_exchange(src, spin_flips(i)))
                    t_flip_possible = .true.

                end if
                cum_sum = cum_sum + elem

                if (t_single .and. t_single_possible .and. tgt == single_excits(i)) then
                    cpt = elem
                else if (t_flip .and. t_flip_possible .and. tgt == spin_flips(i)) then
                    cpt = elem
                end if
            end do
            if (cum_sum < EPS) then
                cpt = 0.0_dp
            else
                cpt = cpt / cum_sum
            end if
        else
            ! create the list depending on the possible excitations
            do i = 1, ubound(neighbors, 1)
                elem = 0.0_dp
                if (IsNotOcc(ilutI, single_excits(i)) .and. &
                    IsNotOcc(ilutI, spin_flips(i))) then
                    ! then the orbital is empty an we can do a hopping
                    ! reuse the hubbard matrix elements..
                    ! or just the -t element?
                    ! since we only need absolute value of matrix element
                    ! it would be better to just use -t .. anyway.. keep it
                    ! general

                    elem = abs(get_offdiag_helement_tJ(nI, [src, single_excits(i)], .false.))
!                     elem = abs(GetTMatEl(src, single_excits(i)))
                    ic_list(i) = 1

                else if (IsOcc(IlutI, spin_flips(i)) .and. &
                         IsNotOcc(IlutI, single_excits(i))) then
                    ! then we can do a spin flip
                    temp_ex(1, 2) = spin_flips(i)
                    temp_ex(2, 2) = single_excits(i)
                    elem = abs(get_offdiag_helement_heisenberg(nI, temp_ex, .false.))
!                      elem = abs(get_heisenberg_exchange(src, spin_flips(i)))
                    ic_list(i) = 2

                else
                    ! if the spin-parallel is occupied, no exciation
                    ! possible, and also prohibit double occupancies
                    elem = 0.0_dp
                end if

                cum_sum = cum_sum + elem
                cum_arr(i) = cum_sum
            end do
        end if

    end subroutine create_cum_list_tJ_model

    function calc_pgen_tJ_model(ilutI, ex, ic) result(pgen)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, 2), ic
        real(dp) :: pgen
#ifdef DEBUG_
        character(*), parameter :: this_routine = "calc_pgen_tJ_model"
#endif

        integer :: src(2), tgt(2)
        real(dp) :: p_elec, p_orb, cum_sum, cpt_1, cpt_2
        real(dp), allocatable :: cum_arr(:)
        integer, allocatable :: tmp_list(:)

        ASSERT(ic >= 0)

        if (ic == 0 .or. ic > 2) then
            pgen = 0.0_dp
            return
        end if

        src = get_src(ex)
        tgt = get_tgt(ex)

        ASSERT(associated(lat))

        p_elec = 1.0_dp / real(nel, dp)

        if (ic == 1) then
            ! here it is easy..
            ASSERT(is_beta(src(1)) .eqv. is_beta(tgt(1)))
            ASSERT(any(tgt(1) == lat%get_spinorb_neighbors(src(1))))
            ASSERT(IsOcc(ilutI, src(1)))
            ASSERT(IsNotOcc(ilutI, tgt(1)))

            call create_cum_list_tJ_model(ilutI, src(1), lat%get_neighbors(gtid(src(1))), &
                                          cum_arr, cum_sum, tmp_list, tgt(1), p_orb)

        else if (ic == 2) then
            ! here we have to find the correct orbital..
            ! and i just realised that we maybe have to take into account
            ! of having picked the orbitals in a different order..
            ! for HPHF reasons i have to return here if the orbitals are
            ! not "correct"
            if (same_spin(src(1), src(2)) .or. same_spin(tgt(1), tgt(2))) then
                pgen = 0.0_dp
                return
            end if
            if (.not. (is_in_pair(src(1), tgt(1)) .or. is_in_pair(src(1), tgt(2)))) then
                pgen = 0.0_dp
                return
            end if
            if (.not. (is_in_pair(src(2), tgt(1)) .or. is_in_pair(src(2), tgt(2)))) then
                pgen = 0.0_dp
                return
            end if

            call create_cum_list_tJ_model(ilutI, src(1), lat%get_neighbors(gtid(src(1))), &
                                          cum_arr, cum_sum, tmp_list, src(2), cpt_1)
            call create_cum_list_tJ_model(ilutI, src(2), lat%get_neighbors(gtid(src(2))), &
                                          cum_arr, cum_sum, tmp_list, src(1), cpt_2)

            p_orb = cpt_1 + cpt_2

#ifdef DEBUG_
            if (is_beta(src(1)) .eqv. is_beta(tgt(1))) then
                ! then those to orbitls were chosen
                ASSERT(is_beta(src(2)) .eqv. is_beta(tgt(2)))
            else if (is_beta(src(1)) .eqv. is_beta(tgt(2))) then
                ASSERT(is_beta(src(2)) .eqv. is_beta(tgt(1)))
            end if

        else
            ! something went wrong
            call stop_all(this_routine, "something went wrong!")
#endif

        end if

        pgen = p_elec * p_orb

    end function calc_pgen_tJ_model

    subroutine pick_orbitals_guga_tJ(ilut, nI, csf_i, excitInfo, orb_pgen)
        ! orbital picking routine for the GUGA t-J model.
        ! the most effective way would be to pick a hole first, instead
        ! of an random electron, so the chance of a succesful hop
        ! increases.. but I might be too lazy for now to do that
        ! actually.
        integer(n_int), intent(in) :: ilut(0:nifguga)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: orb_pgen
        character(*), parameter :: this_routine = "pick_orbitals_guga_tJ"

        integer :: elec, id, ind, tgt
        real(dp) :: p_elec, cum_sum, r, p_orb
        integer, allocatable :: neighbors(:)
        real(dp), allocatable :: cum_arr(:)

        unused_var(ilut)

        elec = 1 + int(genrand_real2_dsfmt() * nel)
        p_elec = 1.0_dp / real(nel, dp)

        id = gtID(nI(elec))

        neighbors = lat%get_neighbors(lat%get_site_index(id))

        call gen_guga_tJ_cum_list(csf_i, id, cum_arr)

        cum_sum = cum_arr(size(neighbors))

        if (cum_sum < EPS) then
            excitInfo%valid = .false.
            orb_pgen = 0.0_dp
            return
        end if

        r = genrand_real2_dsfmt() * cum_sum
        ind = binary_search_first_ge(cum_arr, r)

        tgt = neighbors(ind)

        if (ind == 1) then
            p_orb = cum_arr(1) / cum_sum
        else
            p_orb = (cum_arr(ind) - cum_arr(ind - 1)) / cum_sum
        end if

        if (tgt < id) then
            excitInfo = assign_excitInfo_values_single(gen_type%R, tgt, id, tgt, id)
        else
            excitInfo = assign_excitInfo_values_single(gen_type%L, tgt, id, id, tgt)
        end if

        orb_pgen = p_elec * p_orb

    end subroutine pick_orbitals_guga_tJ

    subroutine gen_guga_tJ_cum_list(csf_i, id, cum_arr, tgt, tgt_pgen)
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: id
        real(dp), allocatable, intent(out) :: cum_arr(:)
        integer, intent(in), optional :: tgt
        real(dp), intent(out), optional :: tgt_pgen
        character(*), parameter :: this_routine = "gen_guga_tJ_cum_list"

        integer, allocatable :: neighbors(:)
        integer :: i, n
        real(dp) :: cum_sum, tmp

        neighbors = lat%get_neighbors(id)

        allocate(cum_arr(size(neighbors)), source=0.0_dp)

        cum_sum = 0.0_dp
        tmp = 0.0_dp

        if (present(tgt)) then
            tgt_pgen = 0.0_dp

            do i = 1, size(neighbors)

                n = neighbors(i)

                if (csf_i%stepvector(n) == 0) then
                    tmp = 1.0_dp
                    cum_sum = cum_sum + tmp
                    if (tgt == n) tgt_pgen = tmp
                end if
            end do
        else
            do i = 1, size(neighbors)
                n = neighbors(i)

                if (csf_i%stepvector(n) == 0) then
                    cum_sum = cum_sum + 1.0_dp
                end if

                cum_arr(i) = cum_sum
            end do
        end if

    end subroutine gen_guga_tJ_cum_list

    subroutine pick_orbitals_guga_heisenberg(ilut, nI, csf_i, excitInfo, orb_pgen)
        ! i "just" need to implement a custom orbital picker for the
        ! spin-free Heisenberg exchange
        integer(n_int), intent(in) :: ilut(0:GugaBits%len_tot)
        integer, intent(in) :: nI(nel)
        type(CSF_Info_t), intent(in) :: csf_i
        type(ExcitationInformation_t), intent(out) :: excitInfo
        real(dp), intent(out) :: orb_pgen
        character(*), parameter :: this_routine = "pick_orbitals_guga_heisenberg"

        integer :: elec, src, id, ind, tgt, start, ende
        real(dp) :: p_elec, cum_sum, r, p_orb
        integer, allocatable :: neighbors(:)
        real(dp), allocatable :: cum_arr(:)

        ! Exists for the function pointer interface
        unused_var(ilut)

        ! here i want to only pick nearest neighbor electrons, where a
        ! spin recoupling is possible

        ! i need to pick one orbital at random: since it should be
        ! general for the t-J model too, pick an occupied orbital
        elec = 1 + int(genrand_real2_dsfmt() * nel)
        p_elec = 1.0_dp / real(nel, dp)

        src = nI(elec)
        ! spatial orbital:
        id = gtID(src)

        neighbors = lat%get_neighbors(lat%get_site_index(id))

        call gen_guga_heisenberg_cum_list(csf_i, id, cum_arr)

        cum_sum = cum_arr(size(neighbors))

        if (cum_sum < EPS) then
            excitInfo%valid = .false.
            orb_pgen = 0.0_dp
            return
        end if

        r = genrand_real2_dsfmt() * cum_sum
        ind = binary_search_first_ge(cum_arr, r)

        tgt = neighbors(ind)

        if (ind == 1) then
            p_orb = cum_arr(1) / cum_sum
        else
            p_orb = (cum_arr(ind) - cum_arr(ind - 1)) / cum_sum
        end if

        ! and now we have to to the GUGA excitation:
        ! were I just realised, everything will get recalculated anyway..
        ! so i have to adapt the pgen recalculation.. to fit this
        ! scheme above..
        start = min(id, tgt)
        ende = max(id, tgt)

        excitInfo = assign_excitInfo_values_double(excit_type%fullstart_stop_mixed, &
               gen_type%L, gen_type%R, gen_type%R, gen_type%R, gen_type%R, &
               start, ende, ende, start, start, start, ende, ende, 0, 2, 1.0_dp, 1.0_dp)

        orb_pgen = p_elec * p_orb

    end subroutine pick_orbitals_guga_heisenberg

    subroutine gen_guga_heisenberg_cum_list(csf_i, id, cum_arr, tgt, tgt_pgen)
        ! make a routine for this, for easy pgen recalculation
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: id
        real(dp), allocatable, intent(out) :: cum_arr(:)
        integer, intent(in), optional :: tgt
        real(dp), intent(out), optional :: tgt_pgen
        character(*), parameter :: this_routine = "gen_guga_heisenberg_cum_list"

        integer :: step, i, n
        integer, allocatable :: neighbors(:)
        real(dp) :: cum_sum, tmp

        neighbors = lat%get_neighbors(lat%get_site_index(id))

        ! then check if the neighbors are available for exchange
        allocate(cum_arr(size(neighbors)), source=0.0_dp)

        cum_sum = 0.0_dp
        tmp = 0.0_dp

        step = csf_i%stepvector(id)

        if (present(tgt)) then

            tgt_pgen = 0.0_dp
            ! if target is not in the neighbors list we can exit
            if (.not. any(tgt == neighbors)) return

            do i = 1, size(neighbors)
                n = neighbors(i)

                if (csf_i%stepvector(n) == 0) then
                    cycle

                else if (csf_i%stepvector(n) == step) then
                    if (abs(id - n) == 1) then
                        cycle
                    else

                        if (step == 1 .and. count_alpha_orbs_ij(csf_i, min(id, n), max(id, n)) == 0) cycle

                        if (step == 2 .and. count_beta_orbs_ij(csf_i, min(id, n), max(id, n)) == 0) cycle

                        tmp = 1.0_dp

                        cum_sum = cum_sum + tmp

                        if (tgt == n) tgt_pgen = tmp
                    end if

                else if (csf_i%stepvector(n) /= 0 &
                         .and. csf_i%stepvector(n) /= step) then

                    if (id - n == -1 &
                            .and. csf_i%stepvector(id) == 1 &
                            .and.  csf_i%B_int(id) == 1) then
                        cycle
                    end if

                    if (id - n == 1 &
                            .and. csf_i%stepvector(n) == 1 &
                            .and. csf_i%B_int(n) == 1) then
                        cycle
                    end if

                    tmp = 1.0_dp

                    cum_sum = cum_sum + tmp

                    if (tgt == n) tgt_pgen = tmp
                end if
            end do

            if (.not. near_zero(cum_sum)) then
                tgt_pgen = tgt_pgen / cum_sum
            end if

        else
            do i = 1, size(neighbors)
                n = neighbors(i)
                if (csf_i%stepvector(n) == 0) then
                    ! t-J case
                    cum_arr(i) = cum_sum
                    cycle

                else if (csf_i%stepvector(n) == step) then
                    ! this can only be if we have space
                    ! between the step-vectors
                    if (abs(id - n) == 1) then
                        cum_arr(i) = cum_sum
                        cycle
                    else
                        if (step == 1 &
                                .and. count_alpha_orbs_ij(csf_i, min(id, n), max(id, n)) == 0) then
                            cum_arr(i) = cum_sum
                            cycle
                        end if

                        if (step == 2 &
                                .and. count_beta_orbs_ij(csf_i, min(id, n), max(id, n)) == 0) then
                            cum_arr(i) = cum_sum
                            cycle
                        end if

                        cum_arr(i) = cum_sum + 1.0_dp
                        cum_sum = cum_sum + 1.0_dp
                    end if

                else if (csf_i%stepvector(n) /= 0 .and. csf_i%stepvector(n) /= step) then

                    if (id - n == -1 &
                            .and. csf_i%stepvector(id) == 1 &
                            .and.  csf_i%B_int(id) == 1) then
                        cum_arr(i) = cum_sum
                        cycle
                    end if

                    if (id - n == 1 &
                            .and. csf_i%stepvector(n) == 1 &
                            .and. csf_i%B_int(n) == 1) then
                        cum_arr(i) = cum_sum
                        cycle
                    end if

                    cum_arr(i) = cum_sum + 1.0_dp
                    cum_sum = cum_sum + 1.0_dp
                end if
            end do
        end if

    end subroutine gen_guga_heisenberg_cum_list

    subroutine calc_orbital_pgen_contr_heisenberg(csf_i, occ_orbs, above_cpt, below_cpt)
        ! and I also need an orbital pgen recalculator for the
        ! exchange type excitations
        type(CSF_Info_t), intent(in) :: csf_i
        integer, intent(in) :: occ_orbs(2)
        real(dp), intent(out) :: above_cpt, below_cpt
        character(*), parameter :: this_routine = "calc_orbital_pgen_contr_heisenberg"

        real(dp) :: p_elec
        real(dp), allocatable :: cum_arr(:)
        integer :: sp_orbs(2)

        ! here I have to somehow recalculate the probability of
        ! picking a pair (i,j)
        p_elec = 1.0_dp / real(nel, dp)

        ! i could have taken both of the orbitals in any order.. so I have
        ! to take this into account

        ASSERT(all(occ_orbs > 0))
        ASSERT(all(occ_orbs <= nBasis))
        ! the occ_orbs are spin-orbitals!  so convert!
        sp_orbs = gtID(occ_orbs)

        call gen_guga_heisenberg_cum_list(csf_i, minval(sp_orbs), cum_arr, &
                                          maxval(sp_orbs), below_cpt)

        call gen_guga_heisenberg_cum_list(csf_i, maxval(sp_orbs), cum_arr, &
                                          minval(sp_orbs), above_cpt)

        above_cpt = above_cpt * p_elec
        below_cpt = below_cpt * p_elec

    end subroutine calc_orbital_pgen_contr_heisenberg

    subroutine create_cum_list_heisenberg(ilutI, src, neighbors, cum_arr, cum_sum, &
                                          tgt, cpt)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: src, neighbors(:)
        real(dp), intent(out), allocatable :: cum_arr(:)
        real(dp), intent(out) :: cum_sum
        integer, intent(in), optional :: tgt
        real(dp), intent(out), optional :: cpt
#ifdef DEBUG_
        character(*), parameter :: this_routine = "create_cum_list_heisenberg"
#endif
        integer :: flip, i, temp_ex(2, maxExcit), nI(nel)
        real(dp) :: elem

        ASSERT(IsOcc(ilutI, src))

        allocate(cum_arr(size(neighbors)))
        cum_arr = 0.0_dp
        cum_sum = 0.0_dp

        ! for the heisenberg model, where we know that every site is singly
        ! occupied, we search for empty spin-orbital neighbors, to see if
        ! a spin-flip is possible! so we do not need to check for the type
        ! of spin of orbital src here!
        ! although for the matrix element i need the opposite spin!
        ! add the according flip to get the other spin!

        ! also use an excitation matrix array to effectively calculate the
        ! transcorrelation factor everywhere needed!
        temp_ex(1, 1) = src

        call decode_bit_det(nI, ilutI)

        if (is_beta(src)) then
            flip = +1
        else
            flip = -1
        end if

        ! add the spinflipped
        temp_ex(2, 1) = src + flip

        if (present(tgt)) then
            ASSERT(present(cpt))
            cpt = 0.0_dp

            do i = 1, ubound(neighbors, 1)
                elem = 0.0_dp
                if (IsNotOcc(ilutI, neighbors(i))) then
                    temp_ex(1, 2) = neighbors(i) + flip
                    temp_ex(2, 2) = neighbors(i)
                    elem = abs(get_offdiag_helement_heisenberg(nI, temp_ex, .false.))
!                     elem = abs(get_heisenberg_exchange(src, neighbors(i)+flip))
                end if
                if (neighbors(i) == tgt) then
                    cpt = elem
                end if
                cum_sum = cum_sum + elem
            end do
            if (cum_sum < EPS) then
                cpt = 0.0_dp
            else
                cpt = cpt / cum_sum
            end if

        else
            do i = 1, ubound(neighbors, 1)
                elem = 0.0_dp
                if (IsNotOcc(ilutI, neighbors(i))) then
                    ! this is a valid orbital to choose from
                    ! but for the matrix element calculation, we need to
                    ! have the opposite spin of the neighboring orbital!
                    temp_ex(1, 2) = neighbors(i) + flip
                    temp_ex(2, 2) = neighbors(i)
                    elem = abs(get_offdiag_helement_heisenberg(nI, temp_ex, .false.))
!                     elem = abs(get_heisenberg_exchange(src, neighbors(i)+flip))
                end if
                cum_sum = cum_sum + elem
                cum_arr(i) = cum_sum
            end do
        end if

    end subroutine create_cum_list_heisenberg

    function calc_pgen_heisenberg_model(ilutI, ex, ic) result(pgen)
        integer(n_int), intent(in) :: ilutI(0:NIfTot)
        integer, intent(in) :: ex(2, ic), ic
        real(dp) :: pgen
#ifdef DEBUG_
        character(*), parameter :: this_routine = "calc_pgen_heisenberg_model"
#endif
        integer :: src(2), tgt(2)
        real(dp) :: p_elec, cum_sum, cpt_1, cpt_2
        real(dp), allocatable :: cum_arr(:)

        ASSERT(ic >= 0)
        ASSERT(associated(lat))
        if (ic /= 2) then
            pgen = 0.0_dp
            return
        end if

        src = get_src(ex)
        tgt = get_tgt(ex)

        ASSERT(.not. same_spin(src(1), src(2)))
        ASSERT(.not. same_spin(tgt(1), tgt(2)))

        p_elec = 1.0_dp / real(nel, dp)

        if (is_beta(src(1)) .eqv. is_beta(tgt(1))) then
            ASSERT(is_in_pair(src(1), tgt(2)))
            ASSERT(is_in_pair(src(2), tgt(1)))
            ASSERT(same_spin(src(2), tgt(2)))

            call create_cum_list_heisenberg(ilutI, src(1), lat%get_spinorb_neighbors(src(1)), &
                                            cum_arr, cum_sum, tgt(1), cpt_1)
            call create_cum_list_heisenberg(ilutI, src(2), lat%get_spinorb_neighbors(src(2)), &
                                            cum_arr, cum_sum, tgt(2), cpt_2)

        else if (is_beta(src(1)) .eqv. is_beta(tgt(2))) then
            ASSERT(is_in_pair(src(1), tgt(1)))
            ASSERT(is_in_pair(src(2), tgt(2)))
            ASSERT(same_spin(src(2), tgt(1)))
            call create_cum_list_heisenberg(ilutI, src(1), lat%get_spinorb_neighbors(src(1)), &
                                            cum_arr, cum_sum, tgt(2), cpt_1)
            call create_cum_list_heisenberg(ilutI, src(2), lat%get_spinorb_neighbors(src(2)), &
                                            cum_arr, cum_sum, tgt(1), cpt_2)
#ifdef DEBUG_
        else
            call stop_all(this_routine, "something went wrong!")
#endif
        end if

        pgen = p_elec * (cpt_1 + cpt_2)

    end function calc_pgen_heisenberg_model

    function get_heisenberg_exchange(src, tgt) result(hel)
        ! this is the wrapper function to get the heisenberg exchange
        ! contribution. which will substitute get_umat in the
        ! necessary places..
        ! NOTE: only in the debug mode the neighboring condition is tested
        ! also that that both orbitals src and tgt are occupied by opposite
        ! spins have to be checked before this function call!
        integer, intent(in) :: src, tgt
        HElement_t(dp) :: hel
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_heisenberg_exchange"
#endif

        ASSERT(src > 0)
        ASSERT(src <= nbasis)
        ASSERT(tgt > 0)
        ASSERT(tgt <= nbasis)

        ASSERT(associated(lat))
        ! i also want to get 0 matrix elements ofc.. so leave the possibility
        ASSERT(allocated(exchange_matrix))

        hel = h_cast(exchange_matrix(src, tgt))

    end function get_heisenberg_exchange

    subroutine setup_spin_free_exchange(in_lat)
        ! spin-free version of the exchange matrix
        class(lattice), intent(in), pointer :: in_lat
        character(*), parameter :: this_routine = "setup_spin_free_exchange"

        integer :: i, ind, iunit

        ASSERT(associated(in_lat))

        if (allocated(spin_free_exchange)) deallocate(spin_free_exchange)
        allocate(spin_free_exchange(nBasis / 2, nBasis / 2), source=0.0_dp)

        ASSERT(in_lat%get_nsites() == nBasis / 2)

        iunit = get_free_unit()
        open(iunit, file = 'spatial-exchange', status = 'replace')

        do i = 1, in_lat%get_nsites()
            ind = in_lat%get_site_index(i)

            ASSERT(ind > 0)
            ASSERT(ind <= nBasis / 2)

            associate(next => in_lat%get_neighbors(ind))
                ASSERT(all(next > 0))
                ASSERT(all(next <= nBasis / 2))

                ! in the spin-free form I have to divide by one 1/2 more!
                spin_free_exchange(ind, next) = -exchange_j / 2.0_dp
                write(iunit, *) ind, next, -exchange_j / 2.0_dp

            end associate

        end do

        close(iunit)

    end subroutine setup_spin_free_exchange

    subroutine setup_exchange_matrix(in_lat)
        ! by convention encode the exchange matrix element by the two
        ! involved electrons with opposite spins! so we can encode this
        ! as a 2D matrix
        class(lattice), intent(in), pointer :: in_lat
        character(*), parameter :: this_routine = "setup_exchange_matrix"
        integer :: i, ind

        ASSERT(associated(in_lat))
        ! create the exchange matrix from the given lattice
        ! connections
        if (allocated(exchange_matrix)) deallocate(exchange_matrix)
        allocate(exchange_matrix(nbasis, nbasis))
        exchange_matrix = 0.0_dp

        ASSERT(in_lat%get_nsites() == nbasis / 2)
        do i = 1, in_lat%get_nsites()
            ind = in_lat%get_site_index(i)
            ! print *, "ind, next:", ind, in_lat%get_neighbors(ind)
            associate(next => in_lat%get_neighbors(ind))
                exchange_matrix(2 * ind - 1, 2 * next) = exchange_j / 2.0_dp
                exchange_matrix(2 * ind, 2 * next - 1) = exchange_j / 2.0_dp

                ASSERT(all(next > 0))
                ASSERT(all(next <= nbasis / 2))
            end associate
            ASSERT(ind > 0)
            ASSERT(ind <= nbasis / 2)
        end do

    end subroutine setup_exchange_matrix

    function get_helement_tJ_ex_mat(nI, ic, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel), ic, ex(2, ic)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel

        if (ic == 0) then
            ! is the diagonal element the same as in the pure
            ! heisenberg model? yes or?
            hel = get_diag_helement_heisenberg(nI)

        else if (ic == 1) then
            ! the single excitation is the same as in the hubbard model
            hel = get_offdiag_helement_tJ(nI, ex(:, 1), tpar)

        else if (ic == 2) then
            hel = get_offdiag_helement_heisenberg(nI, ex, tpar)

        else
            hel = h_cast(0.0_dp)

        end if
    end function get_helement_tJ_ex_mat

    function get_helement_tJ_general(nI, nJ, ic_ret) result(hel)
        integer, intent(in) :: nI(nel), nJ(nel)
        integer, intent(inout), optional :: ic_ret
        HElement_t(dp) :: hel

        integer :: ic, ex(2, maxExcit)
        logical :: tpar
        integer(n_int) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)

        if (present(ic_ret)) then
            if (ic_ret == 0) then
                hel = get_diag_helement_heisenberg(nI)

            else if (ic_ret == 1) then
                ex(1, 1) = 1
                call GetExcitation(nI, nJ, nel, ex, tpar)
                hel = get_offdiag_helement_tJ(nI, ex(:, 1), tpar)

            else if (ic_ret == 2) then
                ex(1, 1) = 2
                call GetExcitation(nI, nJ, nel, ex, tpar)
                hel = get_offdiag_helement_heisenberg(nI, ex, tpar)

            else if (ic_ret == -1) then
                call EncodeBitDet(nI, ilutI)
                call EncodeBitDet(nJ, ilutJ)

                ic_ret = FindBitExcitLevel(ilutI, ilutJ)

                if (ic_ret == 0) then
                    hel = get_diag_helement_heisenberg(nI)

                else if (ic_ret == 1) then
                    ex(1, 1) = 1
                    call GetBitExcitation(ilutI, ilutJ, ex, tpar)
                    hel = get_offdiag_helement_tJ(nI, ex(:, 1), tpar)

                else if (ic_ret == 2) then
                    ex(1, 1) = 2
                    hel = get_offdiag_helement_heisenberg(nI, ex, tpar)

                else
                    hel = h_cast(0.0_dp)
                end if
            else
                hel = h_cast(0.0_dp)
            end if
        else
            call EncodeBitDet(nI, ilutI)
            call EncodeBitDet(nJ, ilutJ)

            ic = FindBitExcitLevel(ilutI, ilutJ)

            if (ic == 0) then
                hel = get_diag_helement_heisenberg(nI)
            else if (ic == 1) then
                ex(1, 1) = 1
                call GetBitExcitation(ilutI, ilutJ, ex, tpar)
                hel = get_offdiag_helement_tJ(nI, ex(:, 1), tPar)
            else if (ic == 2) then
                ex(1, 1) = 2
                call GetBitExcitation(ilutI, ilutJ, ex, tpar)
                hel = get_offdiag_helement_heisenberg(nI, ex, tpar)
            else
                hel = h_cast(0.0_dp)
            end if
        end if

    end function get_helement_tJ_general

    function get_helement_heisenberg_ex_mat(nI, ic, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel), ic, ex(2, ic)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel

        if (ic == 0) then
            hel = get_diag_helement_heisenberg(nI)

        else if (ic == 2) then
            hel = get_offdiag_helement_heisenberg(nI, ex, tpar)

        else
            hel = h_cast(0.0_dp)

        end if

    end function get_helement_heisenberg_ex_mat

    function get_helement_heisenberg_general(nI, nJ, ic_ret) result(hel)
        integer, intent(in) :: nI(nel), nJ(nel)
        integer, intent(inout), optional :: ic_ret
        HElement_t(dp) :: hel

        integer :: ic, ex(2, maxExcit)
        logical :: tpar
        integer(n_int) :: ilutI(0:NIfTot), ilutJ(0:NIfTot)

        if (present(ic_ret)) then
            if (ic_ret == 0) then
                hel = get_diag_helement_heisenberg(nI)

            else if (ic_ret == 2) then
                ex(1, 1) = 2
                call GetExcitation(nI, nJ, nel, ex, tpar)
                hel = get_offdiag_helement_heisenberg(nI, ex, tpar)

            else if (ic_ret == -1) then
                call EncodeBitDet(nI, ilutI)
                call EncodeBitDet(nJ, ilutJ)

                ic_ret = FindBitExcitLevel(ilutI, ilutJ)

                if (ic_ret == 0) then
                    hel = get_diag_helement_heisenberg(nI)

                else if (ic_ret == 2) then
                    ex(1, 1) = 2
                    call GetBitExcitation(ilutI, ilutJ, ex, tpar)

                    hel = get_offdiag_helement_heisenberg(nI, ex, tpar)

                else
                    hel = h_cast(0.0_dp)
                end if
            else
                hel = h_cast(0.0_dp)
            end if
        else
            call EncodeBitDet(nI, ilutI)
            call EncodeBitDet(nJ, ilutJ)

            ic = FindBitExcitLevel(ilutI, ilutJ)

            if (ic == 0) then
                hel = get_diag_helement_heisenberg(nI)

            else if (ic == 2) then
                ex(1, 1) = 2
                call GetBitExcitation(ilutI, ilutJ, ex, tpar)
                hel = get_offdiag_helement_heisenberg(nI, ex, tpar)

            else
                hel = h_cast(0.0_dp)
            end if
        end if

    end function get_helement_heisenberg_general

    function get_diag_helement_heisenberg(nI) result(hel)
        integer, intent(in) :: nI(nel)
        HElement_t(dp) :: hel
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_diag_helement_heisenberg"
#endif
        integer :: i, j, src, flip
        integer, allocatable :: spin_neighbors(:)
        integer(n_int) :: ilut(0:NIfTot)

        ! here i have to sum over all the nearest neigbhor pairs and
        ! find the spin alignment, if it is parallel or not..
        ! in the heisenberg case the contribution form n_i n_j form
        ! nearest neighbors is always the same for every determinant nI
        ! but for the tJ model with less than half-filling this quantitiy is
        ! not a constant.. so i should differentiate in here between tJ
        ! and heisenberg models ..
        ! and what about double counting? should i loop over all the
        ! orbitals and their neighbors or must i differentiate between already
        ! counted contributions?

        ! it is easier to check occupancy in the ilut format

        ASSERT(associated(lat))

        call EncodeBitDet(nI, ilut)
        hel = h_cast(0.0_dp)

        do i = 1, nel
            src = ni(i)
            spin_neighbors = lat%get_spinorb_neighbors(src)
            if (is_beta(src)) then
                flip = +1
            else
                flip = -1
            end if
            do j = 1, size(spin_neighbors)
                if (IsOcc(ilut, spin_neighbors(j))) then
                    ! then it is same spin
                    ! but i really think that we have to take the
                    ! occupancy into account
                    ! and in the tJ model this cancels for parallel spin
                    if (t_heisenberg_model) then
                        hel = hel + h_cast(exchange_j / 4.0_dp)
                    end if
                else if (IsOcc(ilut, spin_neighbors(j) + flip)) then
                    if (t_heisenberg_model) then
                        hel = hel - h_cast(exchange_j / 4.0_dp)
                    else
                        hel = hel - h_cast(exchange_j / 2.0_dp)
                    end if
                    ! it can be empty too, then there is no contribution
                end if
            end do

        end do

        ! i think i would double count if i do not divide by 2..
        ! i hope i am right..
        hel = hel / 2.0_dp

    end function get_diag_helement_heisenberg

    function get_offdiag_helement_heisenberg(nI, ex, tpar) result(hel)
        integer, intent(in) :: nI(nel), ex(2, 2)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_offdiag_helement_heisenberg"
#endif
        integer(n_int) :: ilutI(0:NIfTot)
        integer :: src(2), tgt(2)
        real(dp) :: ni_z, nj_z, spin_fac

        src = get_src(ex)
        tgt = get_tgt(ex)

        ! oh and i have to check the occupancy of nI here or? otherwise this
        ! does not make any sense and is just based on the ex-matrix
        ! it is done in the same way in the slater condon routines, to not
        ! check occupancy.. do it in the debug mode atleast. .
        call EncodeBitDet(nI, ilutI)

        ASSERT(IsOcc(ilutI, src(1)))
        ASSERT(IsOcc(ilutI, src(2)))
        ASSERT(IsNotOcc(ilutI, tgt(1)))
        ASSERT(IsNotOcc(ilutI, tgt(2)))

        ! should i assert the same "same-spinness" here or just return
        ! zero is spins and orbitals do not fit?? i guess that would be
        ! better
        if (same_spin(src(1), src(2))) then
            hel = h_cast(0.0_dp)
        else
            ! i have to check if the orbitals fit.. ex is sorted here or?
            ! can i be sure about that?? check that!

            if (.not. (is_in_pair(src(1), tgt(1)) .and. is_in_pair(src(2), tgt(2)))) then
                hel = h_cast(0.0_dp)

            else
                ! this matrix access checks if the orbitals are connected
                hel = get_heisenberg_exchange(src(1), src(2))

            end if
        end if

        if (t_trans_corr_2body) then
            ! i just realise the sign of the spin densities around the
            ! involved orbitals depend on the spin to be flipped!
            ni_z = get_spin_density_neighbors(ilutI, gtid(src(1)))
            nj_z = get_spin_density_neighbors(ilutI, gtid(src(2)))

            ! beta is defined as negative ms!
            if (is_beta(src(1))) then
                spin_fac = ni_z - nj_z
            else
                spin_fac = nj_z - ni_z
            end if

            hel = hel * exp(2.0_dp * trans_corr_param_2body * (spin_fac - 1.0_dp))

        end if

        if (tpar) hel = -hel

    end function get_offdiag_helement_heisenberg

    function determine_optimal_time_step_heisenberg() result(time_step)
        real(dp) :: time_step
#ifdef DEBUG_
        character(*), parameter :: this_routine = "determine_optimal_time_step_heisenberg"
#endif
        real(dp) :: p_elec, p_hole, mat_ele, max_diag

        p_elec = 1.0_dp / real(nel, dp)

        p_hole = 1.0_dp / real(lat%get_nconnect_max(), dp)

        mat_ele = real(abs(exchange_j), dp)

        time_step = p_elec * p_hole / mat_ele

    end function determine_optimal_time_step_heisenberg

    function get_umat_heisenberg_spin_free(i, j, k, l) result(hel)
        ! for the spin-free form, I do not need information about
        ! the spin-orbitals
        integer, intent(in) :: i, j, k, l
        HElement_t(dp) :: hel
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_umat_heisenberg_spin_free"
#endif

        ASSERT(allocated(spin_free_exchange))

        if (i == j) then
            hel = h_cast(0.0_dp)
        else
            if (i == k .and. j == l) then
                hel = h_cast(spin_free_exchange(i, j))
            else if (i == l .and. j == k) then
                hel = h_cast(spin_free_exchange(i, j))
            else
                hel = h_cast(0.0_dp)
            end if
        end if

    end function get_umat_heisenberg_spin_free

    function get_umat_el_heisenberg(i, j, k, l) result(hel)
        integer, intent(in) :: i, j, k, l
        HElement_t(dp) :: hel
#ifdef DEBUG_
        character(*), parameter :: this_routine = "get_umat_el_heisenberg"
#endif
        ! how exactly do i do that? i access umat with spatial orbitals
        ! so i lose the info if the spins fit.. al i know is that
        ! <ij|kl> i and j are the electrons which must be neighors
        ! and k and l must be in the same orbital as one of the electrons
        ! and what about <ij|kl> = -<ij|lk> symmetry? hm..
        ! i have to think about that and the influence on the matrix element
        ASSERT(allocated(exchange_matrix))

        if (i == j) then
            hel = h_cast(0.0_dp)
        else
            if (i == k .and. j == l) then
                hel = h_cast(exchange_matrix(2 * i, 2 * j - 1))
            else if (i == l .and. j == k) then
                hel = h_cast(exchange_matrix(2 * i, 2 * j - 1))
            else
                hel = h_cast(0.0_dp)
            end if
        end if

    end function get_umat_el_heisenberg

    function get_offdiag_helement_tJ(nI, ex, tpar) result(hel)
        ! if we want to use a transcorrelated Hamiltonian in the tJ case
        ! it is different than in the hubbard model.. so write a one
        ! routine to calculate the one-body matrix element!
        integer, intent(in) :: nI(nel), ex(2)
        logical, intent(in) :: tpar
        HElement_t(dp) :: hel

        integer(n_int) :: ilut(0:niftot)

        ! the first part is exactly the same..
        hel = GetTMatEl(ex(1), ex(2))

        if (tpar) hel = -hel

        if (t_trans_corr) then
            call EncodeBitDet(nI, ilut)

            ! i am not yet sure about the order here.. to have it
            ! "hermitian"
            hel = hel * exp(trans_corr_param) * exp(trans_corr_param * &
                                                    (get_occ_neighbors(ilut, gtid(ex(1))) - get_occ_neighbors(ilut, gtid(ex(2)))))

        end if

        if (t_trans_corr_2body) then
            call EncodeBitDet(nI, ilut)
            hel = hel * exp(trans_corr_param_2body * ( &
                            get_spin_opp_neighbors(ilut, ex(1)) - get_spin_opp_neighbors(ilut, ex(2))))
        end if

    end function get_offdiag_helement_tJ

end module tJ_model
