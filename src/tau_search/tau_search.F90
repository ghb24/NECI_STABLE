#include "macros.h"

module tau_search

    use SystemData, only: AB_elec_pairs, par_elec_pairs, tGen_4ind_weighted, &
                          tHPHF, tKpntSym, nel, G1, nbasis, &
                          AB_hole_pairs, par_hole_pairs, tGen_4ind_reverse, &
                          nOccAlpha, nOccBeta, tUEG, tGen_4ind_2, tReltvy, &
                          t_3_body_excits, t_k_space_hubbard, t_trans_corr_2body, &
                          t_uniform_excits, t_new_real_space_hubbard, &
                          tGUGA, t_mixed_hubbard, t_olle_hubbard, max_ex_level, &
                          t_trans_corr, tHub, t_trans_corr_hop, tNoSinglesPossible, &
                          t_exclude_3_body_excits, t_mol_3_body, t_ueg_3_body, &
                          t_pchb_excitgen, tGAS

    use CalcData

    use FciMCData, only: tRestart, pSingles, pDoubles, pParallel, &
                         ProjEDet, ilutRef, pExcit2, pExcit4, &
                         pExcit2_same, pExcit3_same, &
                         pSing_spindiff1, pDoub_spindiff1, pDoub_spindiff2

    use GenRandSymExcitNUMod, only: construct_class_counts, &
                                    init_excit_gen_store, clean_excit_gen_store

    use tc_three_body_data, only: pTriples

    use SymExcit3, only: GenExcitations3

    use Determinants, only: get_helement

    use SymExcitDataMod, only: excit_gen_store_type

    use bit_rep_data, only: NIfTot

    use bit_reps, only: getExcitationType, decode_bit_det

    use DetBitOps, only: FindBitExcitLevel, TestClosedShellDet, &
                         EncodeBitDet, GetBitExcitation

    use sym_general_mod, only: SymAllowedExcit

    use Parallel_neci

    use constants

    use util_mod, only: near_zero, operator(.isclose.)

    use lattice_mod, only: get_helement_lattice

    use lattice_models_utils, only: gen_all_excits_k_space_hubbard

    implicit none
    ! private

    public :: tSearchTau, tSearchTauOption, MaxTau

    public :: FindMaxTauDoubs, log_spawn_magnitude, init_tau_search, &
        max_death_cpt, log_death_magnitude
    public :: integrate_frequency_histogram_spec

    public :: gamma_sing, gamma_doub, gamma_trip, gamma_opp, gamma_par, &
        cnt_doub, cnt_opp, cnt_par, cnt_sing, cnt_trip, &
        enough_sing, enough_doub, enough_trip, enough_opp, enough_par, &
        update_tau, gamma_sing_spindiff1, gamma_doub_spindiff1, gamma_doub_spindiff2

    ! Tau searching variables
    ! tSearchTau specifies if we are searching tau
    ! tSearchTauOption specifies if we have ever searched for tau
    ! tSearchTauDeath is an override - if we need to adjust tau due to
    !     particle death, when tSearchTau is disabled, but tSearchTauOption
    !     is enabled.
    logical :: tSearchTau, tSearchTauOption
    logical :: tSearchTauDeath
    real(dp) :: MaxTau

! introduce a min_tau value to set a minimum of tau for the automated tau
! search
    logical :: t_min_tau = .false.
    real(dp) :: min_tau_global = 1.0e-7_dp

! alis suggestion: have an option after restarting to keep the time-step
! fixed to the values obtained from the POPSFILE
    logical :: t_keep_tau_fixed = .false.

    logical :: t_test_hist_tau = .false.

    logical :: t_hist_tau_search = .false.
    logical :: t_fill_frequency_hists = .false.


! make variables for automated tau determination, globally available
! 4ind-weighted variables:
    real(dp) :: gamma_sing, gamma_doub, gamma_opp, gamma_par, max_death_cpt, &
                max_permitted_spawn
    real(dp) :: gamma_trip
    logical :: enough_sing, enough_doub, enough_opp, enough_par, consider_par_bias
    logical :: enough_trip
    real(dp) :: gamma_sum

    real(dp) :: gamma_sing_spindiff1, gamma_doub_spindiff1, gamma_doub_spindiff2
    integer :: cnt_sing, cnt_doub, cnt_opp, cnt_par, cnt_trip
! guga-specific:
    integer :: cnt_four, cnt_three_same, cnt_three_mixed, cnt_two_same, cnt_two_mixed
    integer :: n_opp, n_par
    integer :: cnt_sing_hist, cnt_doub_hist, cnt_opp_hist, cnt_par_hist

! guga non-weighted excitation generator tau-update variables
    real(dp) :: gamma_two_same, gamma_two_mixed, gamma_three_same, gamma_three_mixed, &
                gamma_four
    logical :: enough_two, enough_two_same, enough_two_mixed, enough_three, &
               enough_three_same, enough_three_mixed, enough_four

    ! this is to keep probabilities of generating excitations of allowed classes above zero
    real(dp) :: prob_min_thresh

    interface
        ! This is implemented in a submodule
        module subroutine FindMaxTauDoubs()
        end subroutine
    end interface

contains

    subroutine init_tau_search()
        ! N.B. This must be called BEFORE a popsfile is read in, otherwise
        !      we screw up the gamma values that have been carefully read in.
        character(*), parameter :: this_routine = "init_tau_search"

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

        ! Unless it is already specified, set an initial value for tau
        if (.not. tRestart .and. .not. tReadPops .and. near_zero(tau)) then
            if (tGUGA) then
                if (near_zero(MaxTau)) then
                    call stop_all(this_routine, &
                        "please specify a sensible 'max-tau' in input for GUGA calculations!")
                else
                    tau = MaxTau
                end if
            else
                call FindMaxTauDoubs()
            end if
        end if

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

        if (.not. (tReadPops .and. .not. tWalkContGrow)) then
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
        else if(t_pchb_excitgen .and.  .not. tGUGA) then
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

        prob_min_thresh = 1e-8_dp

        t_consider_par_bias = consider_par_bias

    end subroutine init_tau_search

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

    subroutine log_death_magnitude(mult)

        ! The same as above, but for particle death

        real(dp) :: mult

        if (mult > max_death_cpt) then
            tSearchTauDeath = .true.
            max_death_cpt = mult
        end if

    end subroutine

    subroutine update_tau()

        real(dp) :: psingles_new, tau_new, mpi_tmp, tau_death, pParallel_new, pTriples_new
        real(dp) :: pSing_spindiff1_new, pDoub_spindiff1_new, pDoub_spindiff2_new
        logical :: mpi_ltmp
        character(*), parameter :: this_routine = "update_tau"

        ! This is an override. In case we need to adjust tau due to particle
        ! death rates, when it otherwise wouldn't be adjusted
        if (.not. tSearchTau) then

            ! Check that the override has actually occurred.
            ASSERT(tSearchTauOption)
            ASSERT(tSearchTauDeath)

            ! The range of tau is restricted by particle death. It MUST be <=
            ! the value obtained to restrict the maximum death-factor to 1.0.
            call MPIAllReduce(max_death_cpt, MPI_MAX, mpi_tmp)
            max_death_cpt = mpi_tmp
            if (abs(max_death_cpt) > EPS) then
                tau_death = 1.0_dp / max_death_cpt

                ! If this actually constrains tau, then adjust it!
                if (tau_death < tau) then
                    tau = tau_death

                    root_print "******"
                    root_print "WARNING: Updating time step due to particle death &
                         &magnitude"
                    root_print "This occurs despite variable shift mode"
                    root_print "Updating time-step. New time-step = ", tau
                    root_print "******"
                end if
            end if

            ! Condition met --> no need to do this again next iteration
            tSearchTauDeath = .false.

            if (.not. t_hist_tau_search) then
                return
            end if

        end if

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
            if (tau_death < tau_new) then
                if (t_min_tau) then
                    root_print "time-step reduced, due to death events! reset min_tau to:", tau_death
                    min_tau_global = tau_death
                end if
                tau_new = tau_death
            end if
        end if

        ! And a last sanity check/hard limit
        tau_new = min(tau_new, MaxTau)

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
            tau_new = tau_new * 0.99999_dp

            if (abs(tau - tau_new) / tau > 0.001_dp) then
                if (t_min_tau) then
                    if (tau_new < min_tau_global) then
                        root_print "new time-step less than min_tau! set to min_tau:", min_tau_global

                        tau_new = min_tau_global

                    else
                        root_print "Updating time-step. New time-step = ", tau_new, "in: ", this_routine
                    end if
                else
                    root_print "Updating time-step. New time-step = ", tau_new, "in: ", this_routine

                end if
            end if
            tau = tau_new
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
                    root_print " pDoubles = ", (1.0_dp - pSingles_new) * (1.0 - pTriples_new)
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

    subroutine fill_frequency_histogram_nosym_diff(mat_ele, pgen, ic, typ, diff)
        ! specific frequency fill routine for the nosym guga implementation
        ! where no differentiating between mixed or same generators is done
        ! type of excitation is stored in the excitation matrix in the
        ! GUGA implementation!
        ! this ist the implememtation considering diff bias
        real(dp), intent(in) :: mat_ele, pgen
        integer, intent(in) :: ic, typ, diff
        character(*), parameter :: this_routine = "fill_frequency_histogram_nosym_diff"

        real(dp) :: ratio
        integer :: ind
        integer, parameter :: cnt_threshold = 50

        real(dp) :: pBranch2, pBranch3

        ASSERT(pgen > EPS)
        ASSERT(ic == 1 .or. ic == 2)
        ASSERT(typ == 2 .or. typ == 3 .or. typ == 4)
        ASSERT(diff == 0 .or. diff == 1)

        if (mat_ele < EPS) return

        ratio = mat_ele / pgen

        if (ic == 1) then
            ! single excitation
            ratio = ratio * pSingles

            if (ratio < max_frequency_bound) then

                if (.not. enough_sing_hist) then
                    cnt_sing_hist = cnt_sing_hist + 1
                    if (cnt_sing_hist > cnt_threshold) enough_sing_hist = .true.
                end if

                ind = int(ratio / frq_step_size) + 1

                frequency_bins_singles(ind) = frequency_bins_singles(ind) + 1

            end if

        else
            ! double excitation -> check type
            if (typ == 2) then
                pBranch2 = pDoubles * (1.0_dp - pExcit4) * pExcit2
                if (diff == 1) then

                    ratio = ratio * pBranch2 * (1.0_dp - pExcit2_same)

                    if (ratio < max_frequency_bound) then
                        if (.not. enough_two_same) then
                            cnt_type2_same = cnt_type2_same + 1
                            if (cnt_type2_same > cnt_threshold) enough_two_same = .true.
                        end if

                        ind = int(ratio / frq_step_size) + 1

                        frequency_bins_type2(ind) = frequency_bins_type2(ind) + 1

                    end if
                else if (diff == 0) then

                    ratio = ratio * pBranch2 * pExcit2_same

                    if (ratio < max_frequency_bound) then
                        if (.not. enough_two_mixed) then
                            cnt_type2_same = cnt_type2_same + 1
                            if (cnt_type2_same > cnt_threshold) enough_two_mixed = .true.
                        end if

                        ind = int(ratio / frq_step_size) + 1

                        frequency_bins_type2_diff(ind) = frequency_bins_type2_diff(ind) + 1

                    end if
                end if

                if (enough_two_same .and. enough_two_mixed) enough_two = .true.

            else if (typ == 3) then

                pBranch3 = pDoubles * (1.0_dp - pExcit4) * (1.0_dp - pExcit2)
                if (diff == 1) then
                    ratio = ratio * pBranch3 * (1.0_dp - pExcit3_same)

                    if (ratio < max_frequency_bound) then
                        if (.not. enough_three_same) then
                            cnt_type3_same = cnt_type3_same + 1
                            if (cnt_type3_same > cnt_threshold) enough_three_same = .true.
                        end if
                        ind = int(ratio / frq_step_size) + 1

                        frequency_bins_type3(ind) = frequency_bins_type3(ind) + 1

                    end if
                else if (diff == 0) then
                    ratio = ratio * pBranch3 * pExcit3_same

                    if (ratio < max_frequency_bound) then
                        if (.not. enough_three_mixed) then
                            cnt_type3_diff = cnt_type3_diff + 1
                            if (cnt_type3_diff > cnt_threshold) enough_three_mixed = .true.
                        end if
                        ind = int(ratio / frq_step_size) + 1

                        frequency_bins_type3_diff(ind) = frequency_bins_type3_diff(ind) + 1

                    end if
                end if

                if (enough_three_same .and. enough_three_mixed) enough_three = .true.

            else if (typ == 4) then
                ratio = ratio * pDoubles * pExcit4

                if (ratio < max_frequency_bound) then
                    if (.not. enough_four) then
                        cnt_type4 = cnt_type4 + 1
                        if (cnt_type4 > cnt_threshold) enough_four = .true.
                    end if

                    ind = int(ratio / frq_step_size) + 1

                    frequency_bins_type4(ind) = frequency_bins_type4(ind) + 1

                end if
            end if
            if (enough_two .and. enough_three .and. enough_four) enough_doub_hist = .true.

        end if

    end subroutine fill_frequency_histogram_nosym_diff

    subroutine fill_frequency_histogram_nosym_nodiff(mat_ele, pgen, ic, typ)
        ! specific frequency fill routine for the nosym guga implementation
        ! where no differentiating between mixed or same generators is done
        ! type of excitation is stored in the excitation matrix in the
        ! GUGA implementation!
        real(dp), intent(in) :: mat_ele, pgen
        integer, intent(in) :: ic, typ
        character(*), parameter :: this_routine = "fill_frequency_histogram_nosym_nodiff"

        real(dp) :: ratio
        integer :: ind
        integer, parameter :: cnt_threshold = 50

        ASSERT(pgen > EPS)
        ASSERT(ic == 1 .or. ic == 2)
#ifdef DEBUG_
        if (ic == 2) then
            ASSERT(typ == 2 .or. typ == 3 .or. typ == 4)
        end if
#endif

        if (mat_ele < EPS) return

        ratio = mat_ele / pgen

        if (ic == 1) then
            ! single excitation
            ratio = ratio * pSingles

            if (ratio < max_frequency_bound) then
                if (.not. enough_sing_hist) then
                    cnt_sing_hist = cnt_sing_hist + 1
                    if (cnt_sing_hist > cnt_threshold) enough_sing_hist = .true.
                end if
                ind = int(ratio / frq_step_size) + 1

                frequency_bins_singles(ind) = frequency_bins_singles(ind) + 1

            end if

        else
            ! double excitation -> check type
            if (typ == 2) then
                ratio = ratio * pDoubles * (1.0_dp - pExcit4) * pExcit2

                if (ratio < max_frequency_bound) then
                    if (.not. enough_two) then
                        cnt_type2_same = cnt_type2_same + 1
                        if (cnt_type2_same > cnt_threshold) enough_two = .true.
                    end if
                    ind = int(ratio / frq_step_size) + 1

                    frequency_bins_type2(ind) = frequency_bins_type2(ind) + 1

                end if

            else if (typ == 3) then
                ratio = ratio * pDoubles * (1.0_dp - pExcit4) * (1.0_dp - pExcit2)

                if (ratio < max_frequency_bound) then
                    if (.not. enough_three) then
                        cnt_type3_same = cnt_type3_same + 1
                        if (cnt_type3_same > cnt_threshold) enough_three = .true.
                    end if
                    ind = int(ratio / frq_step_size) + 1

                    frequency_bins_type3(ind) = frequency_bins_type3(ind) + 1

                end if

            else if (typ == 4) then
                ratio = ratio * pDoubles * pExcit4

                if (ratio < max_frequency_bound) then
                    if (.not. enough_four) then
                        cnt_type4 = cnt_type4 + 1
                        if (cnt_type4 > cnt_threshold) enough_four = .true.
                    end if
                    ind = int(ratio / frq_step_size) + 1

                    frequency_bins_type4(ind) = frequency_bins_type4(ind) + 1

                end if
            end if
            if (enough_two .and. enough_three .and. enough_four) enough_doub_hist = .true.
        end if

    end subroutine fill_frequency_histogram_nosym_nodiff

    subroutine integrate_frequency_histogram_spec(spec_frequency_bins, ratio)
        ! specific histogram integration routine which sums up the inputted
        ! frequency_bins
        integer, intent(in) :: spec_frequency_bins(n_frequency_bins)
        real(dp), intent(out) :: ratio
        character(*), parameter :: this_routine = "integrate_frequency_histogram_spec"

        integer :: all_frequency_bins(n_frequency_bins)
        integer :: i, threshold
        integer :: n_elements, cnt
        real(dp) :: test_ratio, all_test_ratio
        logical :: mpi_ltmp

        ! test a change to the tau-search by now integrating on each
        ! processor seperately and communicate the maximas
        if (t_test_hist_tau) then
            test_ratio = 0.0_dp
            n_elements = sum(spec_frequency_bins)
            if (n_elements == 0) then
                test_ratio = 0.0_dp

            else if (n_elements < 0) then
                test_ratio = -1.0_dp
                ! if any of the frequency_ratios is full i guess i should
                ! also end the histogramming tau-search or?
                ! yes i have to communicate that.. or else it gets
                ! fucked up..

                t_fill_frequency_hists = .false.

            else

                threshold = int(frq_ratio_cutoff * real(n_elements, dp))
                cnt = 0
                i = 0
                do while (cnt < threshold)
                    i = i + 1
                    cnt = cnt + spec_frequency_bins(i)
                end do

                test_ratio = i * frq_step_size

            end if

            ! how do i best deal with the mpi communication.
            ! i could use a mpialllor on (.not. t_fill_frequency_hists) to
            ! check if one of them is false on any processor..
            call MPIAllLORLogical(.not. t_fill_frequency_hists, mpi_ltmp)
            if (mpi_ltmp) then
                ! then i know one of the frequency histograms is full.. so
                ! stop on all nodes!
                t_fill_frequency_hists = .false.
                ratio = -1.0_dp
                return
            else
                all_test_ratio = 0.0_dp
                call MPIAllReduce(test_ratio, MPI_MAX, all_test_ratio)

                ratio = all_test_ratio
            end if

            return
        end if

        ! MPI communicate
        all_frequency_bins = 0
        call MPIAllReduce(spec_frequency_bins, MPI_SUM, all_frequency_bins)

        n_elements = sum(all_frequency_bins)

        ! have to check if no elements are yet stored into the histogram!
        if (n_elements == 0) then
            ratio = 0.0_dp
            return

        else if (n_elements < 0) then
            ! i reached an integer overflow.. and should stop histogramming
            ! this also means i should make an additional flag for only
            ! the histogramming option without the tau-search to it
            ! so i can also stop just the histogramming after an int
            ! overflow in the histograms
            ! TODO: in this case i also have to decide if i want to print
            ! it at this moment.. or maybe still at the end of the
            ! calculation.. but yes, maybe i want to, by default, always
            ! print them to be able to continue from a certain setting
            call stop_all(this_routine, "Overflow reached")
            ratio = -1.0_dp
            t_fill_frequency_hists = .false.
            return
        end if

        threshold = int(frq_ratio_cutoff * real(n_elements, dp))

        cnt = 0
        i = 0
        do while (cnt < threshold)
            i = i + 1
            cnt = cnt + all_frequency_bins(i)
        end do

        ratio = i * frq_step_size

    end subroutine integrate_frequency_histogram_spec

end module
