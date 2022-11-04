#include "macros.h"

module tau_search_conventional

    use SystemData, only: AB_elec_pairs, par_elec_pairs, tGen_4ind_weighted, &
                          AB_hole_pairs, par_hole_pairs, tGen_4ind_reverse, &
                          nOccAlpha, nOccBeta, tUEG, tGen_4ind_2, tReltvy, &
                          t_k_space_hubbard, t_trans_corr_2body, &
                          t_new_real_space_hubbard, &
                          tGUGA, t_mixed_hubbard, t_olle_hubbard, &
                          t_trans_corr, tHub, t_trans_corr_hop, &
                          t_exclude_3_body_excits, t_ueg_3_body, &
                          t_pchb_excitgen, tGAS

    use FciMCData, only: tRestart, pSingles, pDoubles, pParallel, &
                         ProjEDet, ilutRef, pExcit2, pExcit4, &
                         pExcit2_same, pExcit3_same, &
                         pSing_spindiff1, pDoub_spindiff1, pDoub_spindiff2

    use util_mod, only: clamp

    use tau_main, only: min_tau, max_tau, possible_tau_search_methods, &
            tau_search_method, tau_start_val, possible_tau_start, &
            tau, assign_value_to_tau, max_death_cpt, MaxWalkerBloom

    use tc_three_body_data, only: pTriples

    use bit_reps, only: getExcitationType

    use Parallel_neci, only: MPIAllReduce, MPIAllLORLogical, MPI_MAX, MPI_SUM, iProcIndex

    use constants, only: dp, EPS, stdout, n_int, maxExcit

    use CalcData, only: InitiatorWalkNo, t_consider_par_bias, tTruncInitiator

    use util_mod, only: near_zero, operator(.isclose.), stop_all

    use lattice_mod, only: get_helement_lattice

    use lattice_models_utils, only: gen_all_excits_k_space_hubbard

    use pchb_excitgen, only: FCI_PCHB_particle_selection

    use gasci, only: possible_GAS_exc_gen, GAS_exc_gen

    use gasci_pchb, only: GAS_PCHB_particle_selection

    use gasci_pc_select_particles, only: PCHB_particle_selections

    implicit none
    private

    public :: log_spawn_magnitude, init_tau_search_conventional, &
        finalize_tau_search_conventional, tau_search_stats, update_tau


    type :: TauSearchConventionalStats_t
        real(dp) :: gamma_sing = 0._dp, gamma_doub = 0._dp, gamma_trip = 0._dp, &
                    gamma_opp = 0._dp, gamma_par = 0._dp, &
                    gamma_sing_spindiff1 = 0._dp, gamma_doub_spindiff1 = 0._dp, &
                    gamma_doub_spindiff2 = 0._dp
        integer :: cnt_sing = 0, cnt_doub = 0, cnt_opp = 0, cnt_par = 0, cnt_trip = 0
        logical :: enough_sing = .false., enough_doub = .false., &
                   enough_trip = .false., enough_opp = .false., &
                   enough_par = .false.
    end type

    type(TauSearchConventionalStats_t) :: tau_search_stats

    logical :: consider_par_bias = .false.

    ! this is to keep probabilities of generating excitations of allowed classes above zero
    real(dp), parameter :: prob_min_thresh = 1e-8_dp

contains

    subroutine init_tau_search_conventional()
        ! Don't overwrite stats that were read from popsfile
        if (tau_start_val /= possible_tau_start%from_popsfile) then
            tau_search_stats = TauSearchConventionalStats_t()
        end if

        ! Are we considering parallel-spin bias in the doubles?
        ! Do this logic here, so that if we add opposite spin bias to more
        ! excitation generators, then there is only one place that this logic
        ! needs to be updated!
        if (tGen_4ind_weighted .or. tGen_4ind_2) then
            consider_par_bias = .true.
        else if (tGen_4ind_reverse) then
            consider_par_bias = .true.
        else if (t_k_space_hubbard .and. t_trans_corr_2body) then
            ! for the 2-body transcorrelated k-space hubbard we also have
            ! possible parallel excitations now. and to make the tau-search
            ! working we need to set this to true ofc:
            consider_par_bias = .true.
        else if ((t_pchb_excitgen .and. .not. tGUGA) &
                .and. (FCI_PCHB_particle_selection == PCHB_particle_selections%UNIFORM)) then
            ! The default pchb excitgen also uses parallel biases
            consider_par_bias = .true.
        else if (tGAS &
            .and. (GAS_exc_gen /= possible_GAS_exc_gen%PCHB &
                    .or. GAS_PCHB_particle_selection == PCHB_particle_selections%UNIFORM)) then
            consider_par_bias = .true.
        else
            consider_par_bias = .false.
        end if

        ! If there are only a few electrons in the system, then this has
        ! impacts for the choices that can be made.
        if (nOccAlpha == 0 .or. nOccBeta == 0) then
            consider_par_bias = .false.
            pParallel = 1.0_dp
            tau_search_stats%enough_opp = .true.
        end if
        if (nOccAlpha == 1 .and. nOccBeta == 1) then
            consider_par_bias = .false.
            pParallel = 0.0_dp
            tau_search_stats%enough_par = .true.
        end if

        if (t_mixed_hubbard .or. t_olle_hubbard) then
            pParallel = 0.0_dp
        end if

        t_consider_par_bias = consider_par_bias

    end subroutine init_tau_search_conventional

    subroutine finalize_tau_search_conventional()
        tau_search_stats = TauSearchConventionalStats_t()
    end subroutine

    subroutine log_spawn_magnitude(ic, ex, matel, prob)

        integer, intent(in) :: ic, ex(2, ic)
        real(dp), intent(in) :: prob, matel
        real(dp) :: tmp_gamma, tmp_prob
        integer, parameter :: cnt_threshold = 50

        associate(t_s => tau_search_stats)
        select case (getExcitationType(ex, ic))
        case (1)
            ! no spin changing
            ! Log the details if necessary!
            tmp_prob = prob / pSingles
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > t_s%gamma_sing) t_s%gamma_sing = tmp_gamma
            ! And keep count!
            if (.not. t_s%enough_sing .and. t_s%gamma_sing > 0) then
                t_s%cnt_sing = t_s%cnt_sing + 1
                if (t_s%cnt_sing > cnt_threshold) t_s%enough_sing = .true.
            end if

        case (3)
            ! single spin changing
            ! Log the details if necessary!
            tmp_prob = prob / pSing_spindiff1
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > t_s%gamma_sing_spindiff1) &
                t_s%gamma_sing_spindiff1 = tmp_gamma

            ! And keep count!
            if (.not. t_s%enough_sing .and. tmp_gamma > 0) then
                t_s%cnt_sing = t_s%cnt_sing + 1
                if (t_s%cnt_sing > cnt_threshold) t_s%enough_sing = .true.
            end if

        case (2)
            ! We need to unbias the probability for pDoubles
            tmp_prob = prob / pDoubles

            ! We are not playing around with the same/opposite spin bias
            ! then we should just treat doubles like the singles

            if (consider_par_bias) then
                if (same_spin(ex(1, 1), ex(1, 2))) then
                    tmp_prob = tmp_prob / pParallel
                    tmp_gamma = abs(matel) / tmp_prob
                    if (tmp_gamma > t_s%gamma_par) then
                        t_s%gamma_par = tmp_gamma
                    end if

                    ! And keep count
                    if (.not. t_s%enough_par) then
                        t_s%cnt_par = t_s%cnt_par + 1
                        if (t_s%cnt_par > cnt_threshold) t_s%enough_par = .true.
                        if (t_s%enough_opp .and. t_s%enough_par) t_s%enough_doub = .true.
                    end if
                else
                    tmp_prob = tmp_prob / (1.0_dp - pParallel)
                    tmp_gamma = abs(matel) / tmp_prob
                    if (tmp_gamma > t_s%gamma_opp) then
                        t_s%gamma_opp = tmp_gamma
                    end if

                    ! And keep count
                    if (.not. t_s%enough_opp) then
                        t_s%cnt_opp = t_s%cnt_opp + 1
                        if (t_s%cnt_opp > cnt_threshold) t_s%enough_opp = .true.
                        if (t_s%enough_opp .and. t_s%enough_par) t_s%enough_doub = .true.
                    end if
                end if
            else
                ! We are not playing around with the same/opposite spin bias
                ! then we should just treat doubles like the singles
                tmp_gamma = abs(matel) / tmp_prob
                if (tmp_gamma > t_s%gamma_doub) t_s%gamma_doub = tmp_gamma

                ! And keep count
                if (.not. t_s%enough_doub .and. tmp_gamma > 0) then
                    t_s%cnt_doub = t_s%cnt_doub + 1
                    if (t_s%cnt_doub > cnt_threshold) t_s%enough_doub = .true.
                end if
            end if

        case (4)
            ! We need to unbias the probability for pDoubles
            tmp_prob = prob / pDoub_spindiff1
            ! We are not playing around with the same/opposite spin bias
            ! then we should just treat doubles like the singles
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > t_s%gamma_doub_spindiff1) &
                t_s%gamma_doub_spindiff1 = tmp_gamma
            ! And keep count
            if (.not. t_s%enough_doub .and. tmp_gamma > 0) then
                t_s%cnt_doub = t_s%cnt_doub + 1
                if (t_s%cnt_doub > cnt_threshold) t_s%enough_doub = .true.
            end if

        case (5)
            ! We need to unbias the probability for pDoubles
            tmp_prob = prob / pDoub_spindiff2

            ! We are not playing around with the same/opposite spin bias
            ! then we should just treat doubles like the singles
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > t_s%gamma_doub_spindiff2) &
                t_s%gamma_doub_spindiff2 = tmp_gamma
            ! And keep count
            if (.not. t_s%enough_doub .and. tmp_gamma > 0) then
                t_s%cnt_doub = t_s%cnt_doub + 1
                if (t_s%cnt_doub > cnt_threshold) t_s%enough_doub = .true.
            end if

        case (6)
            ! also treat triple excitations now.
            ! NOTE: but for now this is only done in the transcorrelated
            ! k-space hubbard model, where there are still no single
            ! excitations -> so reuse the quantities for the the singles
            ! instead of introducing yet more variables
            if (.not. t_exclude_3_body_excits) then
                tmp_prob = prob / pTriples
                tmp_gamma = abs(matel) / tmp_prob
            else
                tmp_gamma = 0.0_dp
            end if

            if (tmp_gamma > t_s%gamma_trip) t_s%gamma_trip = tmp_gamma
            ! And keep count!
            if (.not. t_s%enough_trip) then
                t_s%cnt_trip = t_s%cnt_trip + 1
                if (t_s%cnt_trip > cnt_threshold) t_s%enough_trip = .true.
            end if

        end select
        end associate
    end subroutine

    subroutine update_tau()

        real(dp) :: psingles_new, tau_new, mpi_tmp, tau_death, pParallel_new, pTriples_new
        real(dp) :: pSing_spindiff1_new, pDoub_spindiff1_new, pDoub_spindiff2_new
        logical :: mpi_ltmp
        character(*), parameter :: this_routine = "update_tau"
        real(dp) :: gamma_sum


        ASSERT(tau_search_method == possible_tau_search_methods%CONVENTIONAL)

        ! default value for pTriples_new
        pTriples_new = pTriples

        ! What needs doing depends on the number of parameters that are being
        ! updated.

        associate(t_s => tau_search_stats)

        call MPIAllLORLogical(t_s%enough_sing, mpi_ltmp)
        t_s%enough_sing = mpi_ltmp
        call MPIAllLORLogical(t_s%enough_doub, mpi_ltmp)
        t_s%enough_doub = mpi_ltmp
        call MPIAllLORLogical(t_s%enough_trip, mpi_ltmp)
        t_s%enough_trip = mpi_ltmp

        ! Only considering a direct singles/doubles/triples bias
        call MPIAllReduce(t_s%gamma_sing, MPI_MAX, mpi_tmp)
        t_s%gamma_sing = mpi_tmp
        call MPIAllReduce(t_s%gamma_doub, MPI_MAX, mpi_tmp)
        t_s%gamma_doub = mpi_tmp
        call MPIAllReduce(t_s%gamma_trip, MPI_MAX, mpi_tmp)
        t_s%gamma_trip = mpi_tmp
        if (tReltvy) then
            call MPIAllReduce(t_s%gamma_sing_spindiff1, MPI_MAX, mpi_tmp)
            t_s%gamma_sing_spindiff1 = mpi_tmp
            call MPIAllReduce(t_s%gamma_doub_spindiff1, MPI_MAX, mpi_tmp)
            t_s%gamma_doub_spindiff1 = mpi_tmp
            call MPIAllReduce(t_s%gamma_doub_spindiff2, MPI_MAX, mpi_tmp)
            t_s%gamma_doub_spindiff2 = mpi_tmp
            gamma_sum = t_s%gamma_sing + t_s%gamma_sing_spindiff1 + t_s%gamma_doub + t_s%gamma_doub_spindiff1 + t_s%gamma_doub_spindiff2
        else
            gamma_sum = t_s%gamma_sing + t_s%gamma_doub
        end if

        if (consider_par_bias) then
            if (.not. tReltvy) then
                call MPIAllReduce(t_s%gamma_sing, MPI_MAX, mpi_tmp)
                t_s%gamma_sing = mpi_tmp
                call MPIAllReduce(t_s%gamma_opp, MPI_MAX, mpi_tmp)
                t_s%gamma_opp = mpi_tmp
                call MPIAllReduce(t_s%gamma_par, MPI_MAX, mpi_tmp)
                t_s%gamma_par = mpi_tmp
                call MPIAllLORLogical(t_s%enough_opp, mpi_ltmp)
                t_s%enough_opp = mpi_ltmp
                call MPIAllLORLogical(t_s%enough_par, mpi_ltmp)
                t_s%enough_par = mpi_ltmp

                if (t_s%enough_sing .and. t_s%enough_doub) then
                    pparallel_new = t_s%gamma_par / (t_s%gamma_opp + t_s%gamma_par)
                    psingles_new = t_s%gamma_sing * pparallel_new &
                                   / (t_s%gamma_par + t_s%gamma_sing * pparallel_new)
                    tau_new = psingles_new * MaxWalkerBloom &
                              / t_s%gamma_sing
                else
                    pparallel_new = pParallel
                    psingles_new = pSingles
                    if (t_s%gamma_sing > EPS .and. t_s%gamma_par > EPS .and. t_s%gamma_opp > EPS) then
                        tau_new = MaxWalkerBloom * &
                                  min(pSingles / t_s%gamma_sing, &
                                      min(pDoubles * pParallel / t_s%gamma_par, &
                                          pDoubles * (1.0 - pParallel) / t_s%gamma_opp))
                    else
                        ! if no spawns happened, do nothing
                        tau_new = tau
                    end if
                end if

!               checking for triples
                if (t_s%enough_trip) then
                    pTriples_new = t_s%gamma_trip / (t_s%gamma_par + t_s%gamma_sing * pparallel_new + t_s%gamma_trip)
                end if
                ! We only want to update the opposite spins bias here, as we only
                ! consider it here!
                if (t_s%enough_opp .and. t_s%enough_par) then
                    if (abs(pParallel_new - pParallel) / pParallel > 0.0001_dp) then
                        root_print "Updating parallel-spin bias; new pParallel = ", &
                            pParallel_new
                    end if
                    pParallel = pParallel_new
                end if
            else
                call stop_all(this_routine, "Parallel bias is incompatible with magnetic excitation classes")
            end if
        else

            ! Get the probabilities and tau that correspond to the stored
            ! values
            if ((tUEG .or. t_s%enough_sing) .and. t_s%enough_doub) then
                psingles_new = max(t_s%gamma_sing / gamma_sum, prob_min_thresh)
                if (t_s%enough_trip) pTriples_new = max(t_s%gamma_trip / &
                                                    (gamma_sum + t_s%gamma_trip), prob_min_thresh)
                if (tReltvy) then
                    pSing_spindiff1_new = t_s%gamma_sing_spindiff1 / gamma_sum
                    pDoub_spindiff1_new = t_s%gamma_doub_spindiff1 / gamma_sum
                    pDoub_spindiff2_new = t_s%gamma_doub_spindiff2 / gamma_sum
                end if
                tau_new = MaxWalkerBloom / gamma_sum
            else if (t_new_real_space_hubbard .and. t_s%enough_sing .and. &
                     (t_trans_corr_2body .or. t_trans_corr)) then
                ! for the transcorrelated real-space hubbard we could
                ! actually also adapt the time-step!!
                ! but psingles stays 1
                psingles_new = pSingles
                tau_new = MaxWalkerBloom / gamma_sum
            else
                psingles_new = pSingles

                if (tReltvy) then
                    pSing_spindiff1_new = pSing_spindiff1
                    pDoub_spindiff1_new = pDoub_spindiff1
                    pDoub_spindiff2_new = pDoub_spindiff2
                end if
                ! If no single/double spawns occurred, they are also not taken into account
                ! (else would be undefined)
                if (abs(t_s%gamma_doub) > EPS .and. abs(t_s%gamma_sing) > EPS) then
                    tau_new = MaxWalkerBloom * &
                              min(pSingles / t_s%gamma_sing, pDoubles / t_s%gamma_doub)
                else if (abs(t_s%gamma_doub) > EPS) then
                    ! If only doubles were counted, take them
                    tau_new = MaxWalkerBloom * pDoubles / t_s%gamma_doub
                else if (abs(t_s%gamma_sing) > eps) then
                    ! else, we had to have some singles
                    tau_new = MaxWalkerBloom * pSingles / t_s%gamma_sing
                else if (abs(t_s%gamma_trip) > eps) then
                    tau_new = MaxWalkerBloom * PTriples / t_s%gamma_trip
                else
                    ! no spawns
                    tau_new = tau
                end if
            end if

        end if

        ! The range of tau is restricted by particle death. It MUST be <=
        ! the value obtained to restrict the maximum death-factor to 1.0.
        call MPIAllReduce(max_death_cpt, MPI_MAX, mpi_tmp)
        max_death_cpt = mpi_tmp
        ! If there is no death logged, dont do anything
        if (abs(max_death_cpt) > EPS) then
            tau_death = 1.0_dp / max_death_cpt
            if (tau_death < tau_new) then
                if (tau_death < min_tau) then
                    root_print "min-tau reduced, due to death events! reset min_tau to:", tau_death
                    min_tau = tau_death
                end if
                tau_new = tau_death
            end if
        end if

        ! If the calculated tau is less than the current tau, we should ALWAYS
        ! update it. Once we have a reasonable sample of excitations, then we
        ! can permit tau to increase if we have started too low.
        ! make the right if-statements here!
        ! remember t_s%enough_sing is (mis)used for triples in the
        ! 2-body transcorrelated k-space hubbard
        tau_new = clamp(tau_new, min_tau, max_tau)
        if (tau_new < tau .or. (t_s%enough_sing .and. t_s%enough_doub) .or. &
            ((tUEG .and. .not. t_ueg_3_body) .or. tHub .or. t_s%enough_sing .or. &
             (t_k_space_hubbard .and. .not. t_trans_corr_2body) .and. t_s%enough_doub) .or. &
            (t_new_real_space_hubbard .and. t_s%enough_sing .and. &
             (t_trans_corr_2body .or. t_trans_corr)) .or. &
            (t_new_real_space_hubbard .and. t_trans_corr_hop .and. t_s%enough_doub)) then

            ! Make the final tau smaller than tau_new by a small amount
            ! so that we don't get spawns exactly equal to the
            ! initiator threshold, but slightly below it instead.
            ! does this make sense in the new implmentation? this way
            ! i will always decrease the time-step even if its not necessary..
            call assign_value_to_tau(tau_new * 0.99999_dp, 'Conventional tau search')
        end if

        ! Make sure that we have at least some of both singles and doubles
        ! before we allow ourselves to change the probabilities too much...
        if (t_s%enough_sing .and. t_s%enough_doub .and. psingles_new > 1e-5_dp &
            .and. psingles_new < (1.0_dp - 1e-5_dp)) then

            if (abs(psingles - psingles_new) / psingles > 0.0001_dp) then
                if (tReltvy) then
                    root_print "Updating spin-excitation class biases. pSingles(s->s) = ", &
                        psingles_new, ", pSingles(s->s') = ", psing_spindiff1_new, &
                        ", pDoubles(st->st) = ", 1.0_dp - pSingles - pSing_spindiff1_new - pDoub_spindiff1_new - pDoub_spindiff2, &
                        ", pDoubles(st->s't) = ", pDoub_spindiff1_new, &
                        ", pDoubles(st->s't') = ", pDoub_spindiff2_new
                else
                    root_print "Updating singles/doubles bias. pSingles = ", psingles_new
                    root_print " pDoubles = ",(1.0_dp - pSingles_new) * (1.0 - pTriples_new)
                end if
            end if

            pSingles = pSingles_new
            if (tReltvy) then
                pSing_spindiff1 = max(pSing_spindiff1_new, prob_min_thresh)
                pDoub_spindiff1 = max(pDoub_spindiff1_new, prob_min_thresh)
                pDoub_spindiff2 = max(pDoub_spindiff2_new, prob_min_thresh)
                pDoubles = max(1.0_dp - pSingles - pSing_spindiff1_new - pDoub_spindiff1_new - pDoub_spindiff2_new, prob_min_thresh)
                ASSERT(pDoubles - t_s%gamma_doub / gamma_sum < prob_min_thresh)
            else
                pDoubles = 1.0_dp - pSingles
            end if
        end if

        !checking whether we have enouigh triples
        if (t_s%enough_trip) then
            if (abs(pTriples_new - pTriples) / pTriples > 0.0001_dp) then
                root_print "Updating triple-excitation bias. pTriples =", pTriples_new
                pTriples = pTriples_new
            end if
        end if
        end associate

    end subroutine update_tau
end module
