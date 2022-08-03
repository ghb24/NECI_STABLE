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

    use tau_search, only: min_tau, max_tau, possible_tau_search_methods, &
                          input_tau_search_method, tau_search_method, &
                          t_scale_tau_to_death, scale_tau_to_death_triggered, max_death_cpt, &
                          tau_start_val, possible_tau_start, tau, assign_value_to_tau

    use tc_three_body_data, only: pTriples

    use bit_reps, only: getExcitationType

    use Parallel_neci, only: MPIAllReduce, MPIAllLORLogical, MPI_MAX, MPI_SUM, iProcIndex

    use constants, only: dp, EPS, stdout, n_int, maxExcit

    use CalcData, only: InitiatorWalkNo, MaxWalkerBloom, t_consider_par_bias, tReadPops, &
                        tTruncInitiator, tWalkContGrow

    use util_mod, only: near_zero, operator(.isclose.)

    use lattice_mod, only: get_helement_lattice

    use lattice_models_utils, only: gen_all_excits_k_space_hubbard

    implicit none
    private

    public :: find_tau_from_refdet_conn, log_spawn_magnitude, init_tau_search_conventional

    public :: gamma_sing, gamma_doub, gamma_trip, gamma_opp, gamma_par, &
              cnt_doub, cnt_opp, cnt_par, cnt_sing, cnt_trip, &
              enough_sing, enough_doub, enough_trip, enough_opp, enough_par, &
              update_tau, gamma_sing_spindiff1, gamma_doub_spindiff1, gamma_doub_spindiff2

! make variables for automated tau determination, globally available
! 4ind-weighted variables:
    real(dp) :: gamma_sing = 0._dp, gamma_doub = 0._dp, gamma_trip = 0._dp, &
                gamma_opp = 0._dp, gamma_par = 0._dp, &
                gamma_sing_spindiff1 = 0._dp, gamma_doub_spindiff1 = 0._dp, &
                gamma_doub_spindiff2 = 0._dp
    integer :: cnt_sing = 0, cnt_doub = 0, cnt_opp = 0, cnt_par = 0, cnt_trip = 0
    logical :: enough_sing = .false., enough_doub = .false., &
               enough_trip = .false., enough_opp = .false., &
               enough_par = .false., consider_par_bias = .false.
    real(dp) :: gamma_sum, max_permitted_spawn

! guga-specific:
    integer :: n_opp, n_par

    ! this is to keep probabilities of generating excitations of allowed classes above zero
    real(dp), parameter :: prob_min_thresh = 1e-8_dp

    interface
        ! This is implemented in a submodule
        module subroutine find_tau_from_refdet_conn()
        end subroutine
    end interface

contains

    subroutine init_tau_search_conventional()
        ! N.B. This must be called BEFORE a popsfile is read in, otherwise
        !      we screw up the gamma values that have been carefully read in.

        ! We want to start off with zero-values
        gamma_sing = 0.0_dp
        gamma_doub = 0.0_dp
        gamma_trip = 0.0_dp
        gamma_opp = 0.0_dp
        gamma_par = 0.0_dp
        if (tReltvy) then
            gamma_sing_spindiff1 = 0
            gamma_doub_spindiff1 = 0
            gamma_doub_spindiff2 = 0
        end if

        ! And what is the maximum death-component found
        max_death_cpt = 0

        ! And the counts are used to make sure we don't update anything too
        ! early
        cnt_sing = 0
        cnt_doub = 0
        cnt_trip = 0

        cnt_opp = 0
        cnt_par = 0
        enough_sing = .false.
        enough_doub = .false.
        enough_opp = .false.
        enough_par = .false.
        enough_trip = .false.

        write(stdout, *) 'Using initial time-step: ', tau

        ! Set the maximum spawn size
        if (MaxWalkerBloom.isclose.-1._dp) then
            ! No maximum manually specified, so we set the limit of spawn
            ! size to either the initiator criterion, or to 5 otherwise
            if (tTruncInitiator) then
                max_permitted_spawn = InitiatorWalkNo
            else
                max_permitted_spawn = 5.0_dp
            end if
        else
            ! This is specified manually
            max_permitted_spawn = real(MaxWalkerBloom, dp)
        end if

        if (.not.(tReadPops .and. .not. tWalkContGrow)) then
            write(stdout, "(a,f10.5)") "Will dynamically update timestep to &
                         &limit spawning probability to", max_permitted_spawn
        end if

        ! Are we considering parallel-spin bias in the doubles?
        ! Do this logic here, so that if we add opposite spin bias to more
        ! excitation generators, then there is only one place that this logic
        ! needs to be updated!
        if (tGen_4ind_weighted .or. tGen_4ind_2) then
            consider_par_bias = .true.
            !n_opp = AB_elec_pairs
            !n_par = par_elec_pairs
        else if (tGen_4ind_reverse) then
            consider_par_bias = .true.
            n_opp = AB_hole_pairs
            n_par = par_hole_pairs
        else if (t_k_space_hubbard .and. t_trans_corr_2body) then
            ! for the 2-body transcorrelated k-space hubbard we also have
            ! possible parallel excitations now. and to make the tau-search
            ! working we need to set this to true ofc:
            consider_par_bias = .true.
        else if (t_pchb_excitgen .and. .not. tGUGA) then
            ! The default pchb excitgen also uses parallel biases
            consider_par_bias = .true.
        else if (tGAS) then
            consider_par_bias = .true.
        else
            consider_par_bias = .false.
        end if

        ! If there are only a few electrons in the system, then this has
        ! impacts for the choices that can be made.
        if (nOccAlpha == 0 .or. nOccBeta == 0) then
            consider_par_bias = .false.
            pParallel = 1.0_dp
            enough_opp = .true.
        end if
        if (nOccAlpha == 1 .and. nOccBeta == 1) then
            consider_par_bias = .false.
            pParallel = 0.0_dp
            enough_par = .true.
        end if

        if (t_mixed_hubbard .or. t_olle_hubbard) then
            pParallel = 0.0_dp
        end if

        t_consider_par_bias = consider_par_bias

    end subroutine init_tau_search_conventional

    subroutine log_spawn_magnitude(ic, ex, matel, prob)

        integer, intent(in) :: ic, ex(2, ic)
        real(dp), intent(in) :: prob, matel
        real(dp) :: tmp_gamma, tmp_prob
        integer, parameter :: cnt_threshold = 50

        select case (getExcitationType(ex, ic))
        case (1)
            ! no spin changing
            ! Log the details if necessary!
            tmp_prob = prob / pSingles
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > gamma_sing) gamma_sing = tmp_gamma
            ! And keep count!
            if (.not. enough_sing .and. gamma_sing > 0) then
                cnt_sing = cnt_sing + 1
                if (cnt_sing > cnt_threshold) enough_sing = .true.
            end if

        case (3)
            ! single spin changing
            ! Log the details if necessary!
            tmp_prob = prob / pSing_spindiff1
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > gamma_sing_spindiff1) &
                gamma_sing_spindiff1 = tmp_gamma

            ! And keep count!
            if (.not. enough_sing .and. tmp_gamma > 0) then
                cnt_sing = cnt_sing + 1
                if (cnt_sing > cnt_threshold) enough_sing = .true.
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
                    if (tmp_gamma > gamma_par) then
                        gamma_par = tmp_gamma
                    end if

                    ! And keep count
                    if (.not. enough_par) then
                        cnt_par = cnt_par + 1
                        if (cnt_par > cnt_threshold) enough_par = .true.
                        if (enough_opp .and. enough_par) enough_doub = .true.
                    end if
                else
                    tmp_prob = tmp_prob / (1.0_dp - pParallel)
                    tmp_gamma = abs(matel) / tmp_prob
                    if (tmp_gamma > gamma_opp) then
                        gamma_opp = tmp_gamma
                    end if

                    ! And keep count
                    if (.not. enough_opp) then
                        cnt_opp = cnt_opp + 1
                        if (cnt_opp > cnt_threshold) enough_opp = .true.
                        if (enough_opp .and. enough_par) enough_doub = .true.
                    end if
                end if
            else
                ! We are not playing around with the same/opposite spin bias
                ! then we should just treat doubles like the singles
                tmp_gamma = abs(matel) / tmp_prob
                if (tmp_gamma > gamma_doub) gamma_doub = tmp_gamma

                ! And keep count
                if (.not. enough_doub .and. tmp_gamma > 0) then
                    cnt_doub = cnt_doub + 1
                    if (cnt_doub > cnt_threshold) enough_doub = .true.
                end if
            end if

        case (4)
            ! We need to unbias the probability for pDoubles
            tmp_prob = prob / pDoub_spindiff1
            ! We are not playing around with the same/opposite spin bias
            ! then we should just treat doubles like the singles
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > gamma_doub_spindiff1) &
                gamma_doub_spindiff1 = tmp_gamma
            ! And keep count
            if (.not. enough_doub .and. tmp_gamma > 0) then
                cnt_doub = cnt_doub + 1
                if (cnt_doub > cnt_threshold) enough_doub = .true.
            end if

        case (5)
            ! We need to unbias the probability for pDoubles
            tmp_prob = prob / pDoub_spindiff2

            ! We are not playing around with the same/opposite spin bias
            ! then we should just treat doubles like the singles
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > gamma_doub_spindiff2) &
                gamma_doub_spindiff2 = tmp_gamma
            ! And keep count
            if (.not. enough_doub .and. tmp_gamma > 0) then
                cnt_doub = cnt_doub + 1
                if (cnt_doub > cnt_threshold) enough_doub = .true.
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

            if (tmp_gamma > gamma_trip) gamma_trip = tmp_gamma
            ! And keep count!
            if (.not. enough_trip) then
                cnt_trip = cnt_trip + 1
                if (cnt_trip > cnt_threshold) enough_trip = .true.
            end if

        end select
    end subroutine

    subroutine update_tau()

        real(dp) :: psingles_new, tau_new, mpi_tmp, tau_death, pParallel_new, pTriples_new
        real(dp) :: pSing_spindiff1_new, pDoub_spindiff1_new, pDoub_spindiff2_new
        logical :: mpi_ltmp
        character(*), parameter :: this_routine = "update_tau"


        ASSERT(tau_search_method == possible_tau_search_methods%CONVENTIONAL)
        ASSERT(.not. scale_tau_to_death_triggered)

        ! default value for pTriples_new
        pTriples_new = pTriples

        ! What needs doing depends on the number of parameters that are being
        ! updated.

        call MPIAllLORLogical(enough_sing, mpi_ltmp)
        enough_sing = mpi_ltmp
        call MPIAllLORLogical(enough_doub, mpi_ltmp)
        enough_doub = mpi_ltmp
        call MPIAllLORLogical(enough_trip, mpi_ltmp)
        enough_trip = mpi_ltmp

        ! Only considering a direct singles/doubles/triples bias
        call MPIAllReduce(gamma_sing, MPI_MAX, mpi_tmp)
        gamma_sing = mpi_tmp
        call MPIAllReduce(gamma_doub, MPI_MAX, mpi_tmp)
        gamma_doub = mpi_tmp
        call MPIAllReduce(gamma_trip, MPI_MAX, mpi_tmp)
        gamma_trip = mpi_tmp
        if (tReltvy) then
            call MPIAllReduce(gamma_sing_spindiff1, MPI_MAX, mpi_tmp)
            gamma_sing_spindiff1 = mpi_tmp
            call MPIAllReduce(gamma_doub_spindiff1, MPI_MAX, mpi_tmp)
            gamma_doub_spindiff1 = mpi_tmp
            call MPIAllReduce(gamma_doub_spindiff2, MPI_MAX, mpi_tmp)
            gamma_doub_spindiff2 = mpi_tmp
            gamma_sum = gamma_sing + gamma_sing_spindiff1 + gamma_doub + gamma_doub_spindiff1 + gamma_doub_spindiff2
        else
            gamma_sum = gamma_sing + gamma_doub
        end if

        if (consider_par_bias) then
            if (.not. tReltvy) then
                call MPIAllReduce(gamma_sing, MPI_MAX, mpi_tmp)
                gamma_sing = mpi_tmp
                call MPIAllReduce(gamma_opp, MPI_MAX, mpi_tmp)
                gamma_opp = mpi_tmp
                call MPIAllReduce(gamma_par, MPI_MAX, mpi_tmp)
                gamma_par = mpi_tmp
                call MPIAllLORLogical(enough_opp, mpi_ltmp)
                enough_opp = mpi_ltmp
                call MPIAllLORLogical(enough_par, mpi_ltmp)
                enough_par = mpi_ltmp

                if (enough_sing .and. enough_doub) then
                    pparallel_new = gamma_par / (gamma_opp + gamma_par)
                    psingles_new = gamma_sing * pparallel_new &
                                   / (gamma_par + gamma_sing * pparallel_new)
                    tau_new = psingles_new * max_permitted_spawn &
                              / gamma_sing
                else
                    pparallel_new = pParallel
                    psingles_new = pSingles
                    if (gamma_sing > EPS .and. gamma_par > EPS .and. gamma_opp > EPS) then
                        tau_new = max_permitted_spawn * &
                                  min(pSingles / gamma_sing, &
                                      min(pDoubles * pParallel / gamma_par, &
                                          pDoubles * (1.0 - pParallel) / gamma_opp))
                    else
                        ! if no spawns happened, do nothing
                        tau_new = tau
                    end if
                end if

!               checking for triples
                if (enough_trip) then
                    pTriples_new = gamma_trip / (gamma_par + gamma_sing * pparallel_new + gamma_trip)
                end if
                ! We only want to update the opposite spins bias here, as we only
                ! consider it here!
                if (enough_opp .and. enough_par) then
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
            if ((tUEG .or. enough_sing) .and. enough_doub) then
                psingles_new = max(gamma_sing / gamma_sum, prob_min_thresh)
                if (enough_trip) pTriples_new = max(gamma_trip / &
                                                    (gamma_sum + gamma_trip), prob_min_thresh)
                if (tReltvy) then
                    pSing_spindiff1_new = gamma_sing_spindiff1 / gamma_sum
                    pDoub_spindiff1_new = gamma_doub_spindiff1 / gamma_sum
                    pDoub_spindiff2_new = gamma_doub_spindiff2 / gamma_sum
                end if
                tau_new = max_permitted_spawn / gamma_sum
            else if (t_new_real_space_hubbard .and. enough_sing .and. &
                     (t_trans_corr_2body .or. t_trans_corr)) then
                ! for the transcorrelated real-space hubbard we could
                ! actually also adapt the time-step!!
                ! but psingles stays 1
                psingles_new = pSingles
                tau_new = max_permitted_spawn / gamma_sum
            else
                psingles_new = pSingles

                if (tReltvy) then
                    pSing_spindiff1_new = pSing_spindiff1
                    pDoub_spindiff1_new = pDoub_spindiff1
                    pDoub_spindiff2_new = pDoub_spindiff2
                end if
                ! If no single/double spawns occurred, they are also not taken into account
                ! (else would be undefined)
                if (abs(gamma_doub) > EPS .and. abs(gamma_sing) > EPS) then
                    tau_new = max_permitted_spawn * &
                              min(pSingles / gamma_sing, pDoubles / gamma_doub)
                else if (abs(gamma_doub) > EPS) then
                    ! If only doubles were counted, take them
                    tau_new = max_permitted_spawn * pDoubles / gamma_doub
                else if (abs(gamma_sing) > eps) then
                    ! else, we had to have some singles
                    tau_new = max_permitted_spawn * pSingles / gamma_sing
                else if (abs(gamma_trip) > eps) then
                    tau_new = max_permitted_spawn * PTriples / gamma_trip
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
            if (tau_death < tau_new .and. tau_death < min_tau) then
                min_tau = tau_death
                tau_new = tau_death
                root_print "time-step reduced, due to death events! reset min_tau to:", tau_new
            end if
        end if

        ! And a last sanity check/hard limit
        tau_new = min(tau_new, max_tau)

        ! If the calculated tau is less than the current tau, we should ALWAYS
        ! update it. Once we have a reasonable sample of excitations, then we
        ! can permit tau to increase if we have started too low.
        ! make the right if-statements here!
        ! remember enough_sing is (mis)used for triples in the
        ! 2-body transcorrelated k-space hubbard
        if (tau_new < tau .or. (enough_sing .and. enough_doub) .or. &
            ((tUEG .and. .not. t_ueg_3_body) .or. tHub .or. enough_sing .or. &
             (t_k_space_hubbard .and. .not. t_trans_corr_2body) .and. enough_doub) .or. &
            (t_new_real_space_hubbard .and. enough_sing .and. &
             (t_trans_corr_2body .or. t_trans_corr)) .or. &
            (t_new_real_space_hubbard .and. t_trans_corr_hop .and. enough_doub)) then

            ! Make the final tau smaller than tau_new by a small amount
            ! so that we don't get spawns exactly equal to the
            ! initiator threshold, but slightly below it instead.
            tau_new = clamp(tau_new, min_tau, max_tau) * 0.99999_dp

            if (abs(tau - tau_new) / tau > 0.001_dp) then
                call assign_value_to_tau(tau_new)
                root_print "Updating time-step. New time-step = ", tau_new, "in: ", this_routine
            end if
        end if

        ! Make sure that we have at least some of both singles and doubles
        ! before we allow ourselves to change the probabilities too much...
        if (enough_sing .and. enough_doub .and. psingles_new > 1e-5_dp &
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
                ASSERT(pDoubles - gamma_doub / gamma_sum < prob_min_thresh)
            else
                pDoubles = 1.0_dp - pSingles
            end if
        end if

        !checking whether we have enouigh triples
        if (enough_trip) then
            if (abs(pTriples_new - pTriples) / pTriples > 0.0001_dp) then
                root_print "Updating triple-excitation bias. pTriples =", pTriples_new
                pTriples = pTriples_new
            end if
        end if

    end subroutine update_tau
end module
