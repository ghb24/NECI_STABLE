#include "macros.h"

module tau_search

    use SystemData, only: AB_elec_pairs, par_elec_pairs, tGen_4ind_weighted, &
                          tHPHF, tCSF, tKpntSym, nel, G1, nbasis, &
                          AB_hole_pairs, par_hole_pairs, tGen_4ind_reverse, &
                          nOccAlpha, nOccBeta, tUEG, tGen_4ind_2, tReltvy, & 
                          t_3_body_excits, t_k_space_hubbard, t_trans_corr_2body, &
                          t_uniform_excits, t_new_real_space_hubbard, & 
                          t_trans_corr, tHub, t_trans_corr_hop, umateps

    use CalcData, only: tTruncInitiator, tReadPops, MaxWalkerBloom, tau, &
                        InitiatorWalkNo, tWalkContGrow, t_min_tau, min_tau_global, &
                        t_consider_par_bias, &
                    enough_two_same, enough_two_mixed, enough_three_same, &
                    enough_three_mixed, enough_four, enough_two, enough_three, &
                    frequency_bins, & 
                    max_frequency_bound, n_frequency_bins, &
                    frq_step_size, frequency_bins_singles, &
                    frequency_bins_para, frequency_bins_anti,  &
                    frequency_bins_doubles, &
                    frequency_bins_type2, frequency_bins_type3, &
                    frequency_bins_type4, &
                    frequency_bins_type2_diff, frequency_bins_type3_diff, & 
                    frq_ratio_cutoff, t_hist_tau_search,  enough_sing_hist, &
                    enough_doub_hist, enough_par_hist, enough_opp_hist, & 
                    t_fill_frequency_hists, cnt_type3_same, &
                    cnt_type2_same, cnt_type3_diff, &
                    cnt_type4

    use FciMCData, only: tRestart, pSingles, pDoubles, pParallel, &
                         ProjEDet, ilutRef, MaxTau, tSearchTau, &
                         tSearchTauOption, tSearchTauDeath, pExcit2, pExcit4, &
                         pExcit2_same, pExcit3_same, &
                         pSing_spindiff1, pDoub_spindiff1, pDoub_spindiff2

    use GenRandSymExcitNUMod, only: construct_class_counts, &
                                    init_excit_gen_store, clean_excit_gen_store

    use SymExcit3, only: GenExcitations3

    use Determinants, only: get_helement

    use HPHFRandExcitMod, only: ReturnAlphaOpenDet, CalcPGenHPHF, &
                                CalcNonUniPGen

    use HPHF_integrals, only: hphf_off_diag_helement_norm

    use SymExcitDataMod, only: excit_gen_store_type

    use bit_rep_data, only: NIfTot

    use bit_reps, only: getExcitationType, decode_bit_det

    use DetBitOps, only: FindBitExcitLevel, TestClosedShellDet, &
                         EncodeBitDet

    use sym_general_mod, only: SymAllowedExcit

    use Parallel_neci

    use constants

    use k_space_hubbard, only: calc_pgen_k_space_hubbard_uniform_transcorr, &
                               calc_pgen_k_space_hubbard_transcorr, &
                               calc_pgen_k_space_hubbard

    use lattice_mod, only: get_helement_lattice

    use lattice_models_utils, only: gen_all_excits_k_space_hubbard

    implicit none

    real(dp) :: gamma_sing, gamma_doub, gamma_opp, gamma_par, max_death_cpt
    real(dp) :: gamma_sing_spindiff1, gamma_doub_spindiff1, gamma_doub_spindiff2
    real(dp) :: gamma_sum
    real(dp) :: max_permitted_spawn
    integer :: cnt_sing, cnt_doub, cnt_opp, cnt_par
    ! guga-specific:
    integer :: cnt_four, cnt_three_same, cnt_three_mixed, cnt_two_same, cnt_two_mixed
    integer :: n_opp, n_par
    integer :: cnt_sing_hist, cnt_doub_hist, cnt_opp_hist, cnt_par_hist
    logical :: enough_sing, enough_doub, enough_opp, enough_par
    logical :: consider_par_bias

    ! this is to keep probabilities of generating excitations of allowed classes above zero
    real(dp) :: prob_min_thresh


contains

    subroutine init_tau_search ()
        ! N.B. This must be called BEFORE a popsfile is read in, otherwise
        !      we screw up the gamma values that have been carefully read in.
        character(*), parameter :: this_routine = "init_tau_search"

        ! We want to start off with zero-values
        gamma_sing = 0
        gamma_doub = 0

        gamma_opp = 0
        gamma_par = 0
        if (tReltvy) then
            gamma_sing_spindiff1 = 0
            gamma_doub_spindiff1 = 0
            gamma_doub_spindiff2 = 0
        endif

        cnt_opp = 0
        cnt_par = 0

        enough_opp = .false.
        enough_par = .false.

        ! And what is the maximum death-component found
        max_death_cpt = 0

        ! And the counts are used to make sure we don't update anything too
        ! early
        cnt_sing = 0
        cnt_doub = 0

        enough_sing = .false.
        enough_doub = .false.
        ! Unless it is already specified, set an initial value for tau
        if (.not. tRestart .and. .not. tReadPops .and. tau == 0) then
            call FindMaxTauDoubs()
        end if

        write(6,*) 'Using initial time-step: ', tau
        
        ! Set the maximum spawn size
        if (MaxWalkerBloom == -1) then
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
            write(iout, "(a,f10.5)") "Will dynamically update timestep to &
                         &limit spawning probability to", max_permitted_spawn
        end if

        ! Are we considering parallel-spin bias in the doubles?
        ! Do this logic here, so that if we add opposite spin bias to more
        ! excitation generators, then there is only one place that this logic
        ! needs to be updated!
        if (tGen_4ind_weighted .or. tGen_4ind_2 ) then
            !consider_par_bias = .false.
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
        else
            consider_par_bias = .false.
        end if

        t_consider_par_bias = consider_par_bias

        ! If there are only a few electrons in the system, then this has
        ! impacts for the choices that can be made.
        if (nOccAlpha == 0 .or. nOccBeta == 0) then
            consider_par_bias = .false.
            pParallel = 1.0_dp
            enough_opp = .true.
            call stop_all(this_routine, "no electrons in the system?")
        end if
        if (nOccAlpha == 1 .and. nOccBeta == 1) then
            consider_par_bias = .false.
            pParallel = 0.0_dp
            enough_par = .true.
            call stop_all(this_routine, & 
                "do we really need a tau-search for 2 electrons?")
        end if

        prob_min_thresh = 1e-8_dp

        ! since only this routine is called if both tau-search options 
        ! are turned on call it from here too
!         if (t_hist_tau_search) call init_hist_tau_search()

    end subroutine init_tau_search


    subroutine log_spawn_magnitude(ic, ex, matel, prob)

        integer, intent(in) :: ic, ex(2,2)
        real(dp), intent(in) :: prob, matel
        real(dp) :: tmp_gamma, tmp_prob
        integer, parameter :: cnt_threshold = 50
#ifdef __DEBUG
        character(*), parameter :: this_routine = "log_spawn_magnitude"
#endif
        ! i need some changes for 3 body excitations for dynamic tau-search!
!         ASSERT(.not. t_3_body_excits)

        select case(getExcitationType(ex, ic))
        case(1)
            ! no spin changing
            ! Log the details if necessary!
            tmp_prob = prob / pSingles
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > gamma_sing) &
                gamma_sing = tmp_gamma
            
            ! And keep count!
            if (.not. enough_sing) then
                cnt_sing = cnt_sing + 1
                if (cnt_sing > cnt_threshold) enough_sing = .true.
            endif

        case(3)
            ! single spin changing
            ! Log the details if necessary!
            tmp_prob = prob / pSing_spindiff1
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > gamma_sing_spindiff1) &
                gamma_sing_spindiff1 = tmp_gamma
            
            ! And keep count!
            if (.not. enough_sing) then
                cnt_sing = cnt_sing + 1
                if (cnt_sing > cnt_threshold) enough_sing = .true.
            endif
           
        case(2) 
            ! We need to unbias the probability for pDoubles
            tmp_prob = prob / pDoubles

            if (consider_par_bias) then 
                if (same_spin(ex(1,1),ex(1,2))) then 
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
                if (.not. enough_doub) then
                    cnt_doub = cnt_doub + 1
                    if (cnt_doub > cnt_threshold) enough_doub = .true.
                end if
            end if
     
        case(4) 
            ! We need to unbias the probability for pDoubles
            tmp_prob = prob / pDoub_spindiff1
            ! We are not playing around with the same/opposite spin bias
            ! then we should just treat doubles like the singles
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > gamma_doub_spindiff1) &
                gamma_doub_spindiff1 = tmp_gamma
            ! And keep count
            if (.not. enough_doub) then
                cnt_doub = cnt_doub + 1
                if (cnt_doub > cnt_threshold) enough_doub = .true.
            endif

        case(5) 
            ! We need to unbias the probability for pDoubles
            tmp_prob = prob / pDoub_spindiff2

            ! We are not playing around with the same/opposite spin bias
            ! then we should just treat doubles like the singles
            tmp_gamma = abs(matel) / tmp_prob
            if (tmp_gamma > gamma_doub_spindiff2) &
                gamma_doub_spindiff2 = tmp_gamma
            ! And keep count
            if (.not. enough_doub) then
                cnt_doub = cnt_doub + 1
                if (cnt_doub > cnt_threshold) enough_doub = .true.
            endif

        case(6) 
            ! also treat triple excitations now. 
            ! NOTE: but for now this is only done in the transcorrelated 
            ! k-space hubbard model, where there are still no single 
            ! excitations -> so reuse the quantities for the the singles 
            ! instead of introducing yet more variables
            tmp_prob = prob / (1.0_dp - pDoubles)
            tmp_gamma = abs(matel) / tmp_prob

            if (tmp_gamma > gamma_sing) gamma_sing = tmp_gamma
            ! And keep count!
            if (.not. enough_sing) then
                cnt_sing = cnt_sing + 1
                if (cnt_sing > cnt_threshold) enough_sing = .true.
            endif

        end select

     end subroutine

    subroutine log_death_magnitude (mult)

        ! The same as above, but for particle death

        real(dp) :: mult

        if (mult > max_death_cpt) then
            tSearchTauDeath = .true.
            max_death_cpt = mult
        end if

    end subroutine

    subroutine update_tau()

        use FcimCData, only: iter

        real(dp) :: psingles_new, tau_new, mpi_tmp, tau_death, pParallel_new
        real(dp) :: pSing_spindiff1_new, pDoub_spindiff1_new, pDoub_spindiff2_new
        logical :: mpi_ltmp
        character(*), parameter :: this_routine = "update_tau"

        integer :: itmp, itmp2

        ! This is an override. In case we need to adjust tau due to particle
        ! death rates, when it otherwise wouldn't be adjusted
        if (.not. tSearchTau) then
            
            ! Check that the override has actually occurred.
            ASSERT(tSearchTauOption)
            ASSERT(tSearchTauDeath)

            ! The range of tau is restricted by particle death. It MUST be <=
            ! the value obtained to restrict the maximum death-factor to 1.0.
            call MPIAllReduce (max_death_cpt, MPI_MAX, mpi_tmp)
            max_death_cpt = mpi_tmp
            if(abs(max_death_cpt) > EPS) then
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

        ! What needs doing depends on the number of parameters that are being
        ! updated.

        call MPIAllLORLogical(enough_sing, mpi_ltmp)
        enough_sing = mpi_ltmp
        call MPIAllLORLogical(enough_doub, mpi_ltmp)
        enough_doub = mpi_ltmp

        ! Only considering a direct singles/doubles bias
        call MPIAllReduce (gamma_sing, MPI_MAX, mpi_tmp)
        gamma_sing = mpi_tmp
        call MPIAllReduce (gamma_doub, MPI_MAX, mpi_tmp)
        gamma_doub = mpi_tmp
        if (tReltvy) then
            call MPIAllReduce (gamma_sing_spindiff1, MPI_MAX, mpi_tmp)
            gamma_sing_spindiff1 = mpi_tmp
            call MPIAllReduce (gamma_doub_spindiff1, MPI_MAX, mpi_tmp)
            gamma_doub_spindiff1 = mpi_tmp
            call MPIAllReduce (gamma_doub_spindiff2, MPI_MAX, mpi_tmp)
            gamma_doub_spindiff2 = mpi_tmp
            gamma_sum = gamma_sing + gamma_sing_spindiff1 + gamma_doub + gamma_doub_spindiff1 + gamma_doub_spindiff2
        else
            gamma_sum = gamma_sing + gamma_doub
        endif

        if (consider_par_bias) then
            if (.not. tReltvy) then
                call MPIAllReduce (gamma_sing, MPI_MAX, mpi_tmp)
                gamma_sing = mpi_tmp
                call MPIAllReduce (gamma_opp, MPI_MAX, mpi_tmp)
                gamma_opp = mpi_tmp
                call MPIAllReduce (gamma_par, MPI_MAX, mpi_tmp)
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
                    if(gamma_sing > EPS .and. gamma_par > EPS .and. gamma_opp > EPS) then 
                       tau_new = max_permitted_spawn * &
                            min(pSingles / gamma_sing, &
                            min(pDoubles * pParallel / gamma_par, &
                                pDoubles * (1.0_dp - pParallel) / gamma_opp))
                    else                       
                       ! if no spawns happened, do nothing
                       tau_new = tau
                    endif
                end if

                ! We only want to update the opposite spins bias here, as we only
                ! consider it here!
                if (enough_opp .and. enough_par) then
                    if (abs(pParallel_new-pParallel) / pParallel > 0.0001_dp) then
                        root_print "Updating parallel-spin bias; new pParallel = ", &
                            pParallel_new
                    end if
                    pParallel = pParallel_new
                end if
            else
                call stop_all(this_routine, "Parallel bias is incompatible with magnetic excitation classes")
            endif
        else

            ! Get the probabilities and tau that correspond to the stored
            ! values
            if ((tUEG .or. tHub .or. t_k_space_hubbard .or. enough_sing) .and. enough_doub) then
                psingles_new = gamma_sing / gamma_sum
                if (tReltvy) then
                    pSing_spindiff1_new = gamma_sing_spindiff1/gamma_sum
                    pDoub_spindiff1_new = gamma_doub_spindiff1/gamma_sum
                    pDoub_spindiff2_new = gamma_doub_spindiff2/gamma_sum
                endif
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
                endif
           ! If no single/double spawns occurred, they are also not taken into account
           ! (else would be undefined)
                if(abs(gamma_doub) > EPS .and. abs(gamma_sing) > EPS) then
                   tau_new = max_permitted_spawn * &
                        min(pSingles / gamma_sing, pDoubles / gamma_doub)
                else if(abs(gamma_doub) > EPS) then
                   ! If only doubles were counted, take them
                   tau_new = max_permitted_spawn * pDoubles / gamma_doub
                else
                   ! else, we had to have some singles
                   tau_new = max_permitted_spawn * pSingles / gamma_sing
                endif
            end if

        end if

        ! The range of tau is restricted by particle death. It MUST be <=
        ! the value obtained to restrict the maximum death-factor to 1.0.
        call MPIAllReduce (max_death_cpt, MPI_MAX, mpi_tmp)
        max_death_cpt = mpi_tmp
        ! If there is no death logged, dont do anything
        if(abs(max_death_cpt) > EPS) then
           tau_death = 1.0_dp / max_death_cpt
           if (tau_death < tau_new) then
              if (t_min_tau) then
                 root_print "time-step reduced, due to death events! reset min_tau to:", tau_death
                 min_tau_global = tau_death
              end if
              tau_new = tau_death
           end if
        endif

        ! And a last sanity check/hard limit
        tau_new = min(tau_new, MaxTau)

        ! If the calculated tau is less than the current tau, we should ALWAYS
        ! update it. Once we have a reasonable sample of excitations, then we
        ! can permit tau to increase if we have started too low.
        ! make the right if-statements here! 
        ! remember enough_sing is (mis)used for triples in the 
        ! 2-body transcorrelated k-space hubbard 
        if (tau_new < tau .or. & 
            (tUEG .or. tHub .or. enough_sing .or. & 
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
                        root_print "new time-step less then min_tau! set to min_tau!", min_tau_global

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
                    root_print " pDoubles = ", 1.0_dp - psingles_new
                endif
            end if
            pSingles = max(psingles_new, prob_min_thresh)
            if (tReltvy) then
                pSing_spindiff1 = max(pSing_spindiff1_new, prob_min_thresh)
                pDoub_spindiff1 = max(pDoub_spindiff1_new, prob_min_thresh)
                pDoub_spindiff2 = max(pDoub_spindiff2_new, prob_min_thresh)
                pDoubles = max(1.0_dp - pSingles - pSing_spindiff1_new - pDoub_spindiff1_new - pDoub_spindiff2_new, prob_min_thresh)
                ASSERT(pDoubles-gamma_doub/gamma_sum < prob_min_thresh)
            else
                pDoubles = 1.0_dp - pSingles
            endif
        end if

    end subroutine update_tau

    subroutine FindMaxTauDoubs()

        ! Routine to find an upper bound to tau, by consideration of the
        ! singles and doubles connected to the reference determinant
        ! 
        ! Obviously, this make assumptions about the possible range of pgen,
        ! so may actually give a tau that is too SMALL for the latest
        ! excitation generators, which is exciting!

        use neci_intfce
        use SymExcit4, only : GenExcitations4, ExcitGenSessionType
        type(excit_gen_store_type) :: store, store2
        logical :: tAllExcitFound,tParity,tSameFunc,tSwapped,tSign
        character(len=*), parameter :: t_r="FindMaxTauDoubs"
        character(len=*), parameter :: this_routine ="FindMaxTauDoubs"
        integer :: ex(2,2),ex2(2,2),exflag,iMaxExcit,nStore(6),nExcitMemLen(1)
        integer, allocatable :: Excitgen(:)
        real(dp) :: nAddFac,MagHel,pGen,pGenFac
        HElement_t(dp) :: hel
        integer :: ic,nJ(nel),nJ2(nel),ierr,iExcit,ex_saved(2,2)
        integer(kind=n_int) :: iLutnJ(0:niftot),iLutnJ2(0:niftot)

        type(ExcitGenSessionType) :: session

        integer(n_int), allocatable :: det_list(:,:) 
        integer :: n_excits, i, ex_3(2,3)

        if(tCSF) call stop_all(t_r,"TauSearching needs fixing to work with CSFs or MI funcs")

        if(MaxWalkerBloom.eq.-1) then
            !No MaxWalkerBloom specified
            !Therefore, assume that we do not want blooms larger than n_add if initiator,
            !or 5 if non-initiator calculation.
            if(tTruncInitiator) then
                nAddFac = InitiatorWalkNo
            else
                nAddFac = 5.0_dp    !Won't allow more than 5 particles at a time
            endif
        else
            nAddFac = real(MaxWalkerBloom,dp) !Won't allow more than MaxWalkerBloom particles to spawn in one event. 
        endif

        Tau = 1000.0_dp

        ! NOTE: test if the new real-space implementation works with this 
        ! function! maybe i also have to use a specific routine for this ! 
        ! since it might be necessary in the transcorrelated approach to 
        ! the real-space hubbard
!         if (t_new_real_space_hubbard) then 
!             call Stop_All(this_routine, "does this routine work correctly? test it!")
!         end if

        ! bypass everything below for the new k-space hubbard implementation
        if (t_k_space_hubbard) then 
            if (tHPHF) then 
                call Stop_All(this_routine, &
                    "not yet implemented with HPHF, since gen_all_excits not atapted to it!")
            end if

            call gen_all_excits_k_space_hubbard(ProjEDet(:,1), n_excits, det_list)

            ! now loop over all of them and determine the worst case H_ij/pgen ratio
            do i = 1, n_excits 
                call decode_bit_det(nJ, det_list(:,i))
                ! i have to take the right direction in the case of the 
                ! transcorrelated, due to non-hermiticity..
                ic = FindBitExcitlevel(det_list(:,i), ilutRef(:,1))
                ASSERT(ic == 2 .or. ic == 3)
                if (ic == 2) then 
                    call GetBitExcitation(ilutRef(:,1), det_list(:,i), ex, tParity)
                else if (ic == 3) then 
                    call GetBitExcitation(ilutRef(:,1), det_list(:,i), ex_3, tParity)
                end if

                MagHel = abs(get_helement_lattice(nJ, ProjEDet(:,1)))
                ! and also get the generation probability 
                if (t_trans_corr_2body) then 
                    if (t_uniform_excits) then 
                        ! i have to setup pDoubles and the other quantities 
                        ! before i call this functionality! 
                        pgen = calc_pgen_k_space_hubbard_uniform_transcorr(& 
                            ProjEDet(:,1), ilutRef(:,1), ex_3, ic)
                    else 
                        pgen = calc_pgen_k_space_hubbard_transcorr(&
                            ProjEDet(:,1), ilutRef(:,1), ex_3, ic)
                    end if
                else
                    pgen = calc_pgen_k_space_hubbard(&
                            ProjEDet(:,1), ilutRef(:,1), ex, ic)
                end if

                if (MagHel > EPS) then 
                    pGenFac = pgen * nAddFac / MagHel

                    if (tau > pGenFac .and. pGenFac > EPS) then 
                        tau = pGenFac
                    end if
                end if
            end do

            if(tau.gt.0.075_dp) then
                tau=0.075_dp
                write(iout,"(A,F8.5,A)") "Small system. Setting initial timestep to be ",Tau," although this &
                                                &may be inappropriate. Care needed"
            else
                write(iout,"(A,F18.10)") "From analysis of reference determinant and connections, &
                                         &an upper bound for the timestep is: ",Tau
            endif

            return
        end if

        tAllExcitFound=.false.
        Ex_saved(:,:)=0
        exflag=3
        tSameFunc=.false.
        call init_excit_gen_store(store)
        call init_excit_gen_store(store2)
        store%tFilled = .false.
        store2%tFilled = .false.
        CALL construct_class_counts(ProjEDet(:,1), store%ClassCountOcc, &
                                    store%ClassCountUnocc)
        store%tFilled = .true.
        if(tKPntSym) then
            !TODO: It REALLY needs to be fixed so that we don't need to do this!!
            !Setting up excitation generators that will work with kpoint sampling
            iMaxExcit=0
            nStore(:)=0
            CALL GenSymExcitIt2(ProjEDet(:,1),NEl,G1,nBasis,.TRUE.,nExcitMemLen,nJ,iMaxExcit,nStore,exFlag)
            ALLOCATE(EXCITGEN(nExcitMemLen(1)),stat=ierr)
            IF(ierr.ne.0) CALL Stop_All(t_r,"Problem allocating excitation generator")
            EXCITGEN(:)=0
            CALL GenSymExcitIt2(ProjEDet(:,1),NEl,G1,nBasis,.TRUE.,EXCITGEN,nJ,iMaxExcit,nStore,exFlag)
        !    CALL GetSymExcitCount(EXCITGEN,DetConn)
        endif

        do while (.not.tAllExcitFound)
            if(tKPntSym) then
                call GenSymExcitIt2(ProjEDet(:,1),nel,G1,nBasis,.false.,EXCITGEN,nJ,iExcit,nStore,exFlag)
                if(nJ(1).eq.0) exit
                !Calculate ic, tParity and Ex
                call EncodeBitDet (nJ, iLutnJ)
                Ex(:,:)=0
                ic = FindBitExcitlevel(iLutnJ,iLutRef(:,1),2)
                ex(1,1) = ic
                call GetExcitation(ProjEDet(:,1),nJ,Nel,ex,tParity)
            else
                if (tReltvy) then
                    call GenExcitations4(session, ProjEDet(:,1), nJ, exflag, ex_saved, tParity, tAllExcitFound, .false.)
                else
                    CALL GenExcitations3(ProjEDet(:,1),iLutRef(:,1),nJ,exflag,Ex_saved,tParity,tAllExcitFound,.false.)
                endif

                IF(tAllExcitFound) EXIT
                Ex(:,:) = Ex_saved(:,:)
                if(Ex(2,2).eq.0) then
                    ic=1
                else
                    ic=2
                endif
                call EncodeBitDet (nJ, iLutnJ)
            endif

            ! Exclude an excitation if it isn't symmetry allowed.
            ! Note that GenExcitations3 is not perfect, especially if there
            ! additional restrictions, such as LzSymmetry.
            if (.not. SymAllowedExcit(ProjEDet(:,1), nJ, ic, ex)) &
                cycle

            !Find Hij
            if (tGUGA) then 
                call stop_all(t_r, "modify get_helement for GUGA")
            end if
            if(tHPHF) then
                if(.not.TestClosedShellDet(iLutnJ)) then
                    CALL ReturnAlphaOpenDet(nJ,nJ2,iLutnJ,iLutnJ2,.true.,.true.,tSwapped)
                    if(tSwapped) then
                        !Have to recalculate the excitation matrix.
                        ic = FindBitExcitLevel(iLutnJ, iLutRef(:,1), 2)
                        ex(:,:) = 0
                        ASSERT(.not. t_3_body_excits)
                        if(ic.le.2) then
                            ex(1,1) = ic
                            call GetBitExcitation(iLutRef(:,1),iLutnJ,Ex,tParity)
                        endif
                    endif
                endif
                hel = hphf_off_diag_helement_norm(ProjEDet(:,1),nJ,iLutRef(:,1),iLutnJ)
            else
                hel = get_helement(ProjEDet(:,1),nJ,ic,ex,tParity)
            endif

            MagHel = abs(hel)

            !Find pGen (nI -> nJ)
            if(tHPHF) then
                call CalcPGenHPHF(ProjEDet(:,1),iLutRef(:,1),nJ,iLutnJ,ex,store%ClassCountOcc,    &
                            store%ClassCountUnocc,pDoubles,pGen,tSameFunc)
            else
                call CalcNonUnipGen(ProjEDet(:,1),ilutRef(:,1),ex,ic,store%ClassCountOcc,store%ClassCountUnocc,pDoubles,pGen)
            endif
            if(tSameFunc) cycle
            if(MagHel.gt.0.0_dp) then
                pGenFac = pGen*nAddFac/MagHel
                if(Tau.gt.pGenFac .and. pGenFac > EPS) then
                    Tau = pGenFac
                endif
            endif

            !Find pGen(nJ -> nI)
            CALL construct_class_counts(nJ, store2%ClassCountOcc, &
                                        store2%ClassCountUnocc)
            store2%tFilled = .true.
            if(tHPHF) then
                ic = FindBitExcitLevel(iLutnJ, iLutRef(:,1), 2)
                ex2(:,:) = 0
                ASSERT(.not. t_3_body_excits)
                if(ic.le.2) then
                    ex2(1,1) = ic

                    call GetBitExcitation(iLutnJ,iLutRef(:,1),Ex2,tSign)
                endif
                call CalcPGenHPHF(nJ,iLutnJ,ProjEDet(:,1),iLutRef(:,1),ex2,store2%ClassCountOcc,    &
                            store2%ClassCountUnocc,pDoubles,pGen,tSameFunc)
            else
                ex2(1,:) = ex(2,:)
                ex2(2,:) = ex(1,:)
                call CalcNonUnipGen(nJ,ilutnJ,ex2,ic,store2%ClassCountOcc,store2%ClassCountUnocc,pDoubles,pGen)
            endif
            if(tSameFunc) cycle
            if(MagHel.gt.0.0_dp) then
                pGenFac = pGen*nAddFac/MagHel
                if(Tau.gt.pGenFac .and. pGenFac > EPS) then
                    Tau = pGenFac
                endif
            endif

        enddo
                
        call clean_excit_gen_store (store)
        call clean_excit_gen_store (store2)
        if(tKPntSym) deallocate(EXCITGEN)

        if(tau.gt.0.075_dp) then
            tau=0.075_dp
            write(iout,"(A,F8.5,A)") "Small system. Setting initial timestep to be ",Tau," although this &
                                            &may be inappropriate. Care needed"
        else
            write(iout,"(A,F18.10)") "From analysis of reference determinant and connections, &
                                     &an upper bound for the timestep is: ",Tau
        endif

    end subroutine FindMaxTauDoubs

    subroutine fill_frequency_histogram_4ind(mat_ele, pgen, ic, t_parallel)
        ! this is the specific routine to fill up the frequency histograms 
        ! for the 4ind-weighted excitation generators, which use pParallel too
        ! (do they always by default?? check that! todo! 
        ! i just realised, due to the linear and ordered bins, i actually dont
        ! need to binary search in them but can determine the index just 
        ! through the ratio and frequency step size
        ! i want to stop filling the histograms after a certain number of
        ! excitations is stored. 1.: since i dont want integer overflows and 
        ! i think it is unnecessary to store it as an 64bit integer array 
        ! and 2.: the distribution shouldn't change too much after a 
        ! certain point.. 
        ! but this would also mean, that the tau-update would not give 
        ! any new values after i dont change the histograms anymore 
        ! so i could stop the tau-adaptation after that point too.. 
        ! so maybe first experiment a bit with 64bit array to see if 
        ! the time does change more after the 32 bit limit is reached 
        ! and then change it back to 32 bits..
        real(dp), intent(in) :: mat_ele, pgen 
        integer, intent(in) :: ic 
        logical, intent(in) :: t_parallel
        character(*), parameter :: this_routine = "fill_frequency_histogram_4ind"
        integer, parameter :: cnt_threshold = 50

        real(dp) :: ratio
        integer :: new_n_bins, old_n_bins, i, ind 
        integer, allocatable :: save_bins(:)


        ASSERT(pgen > EPS) 
        ASSERT(ic == 1 .or. ic == 2)

        if (mat_ele < EPS) return 

        ratio = mat_ele / pgen

        ! then i have to decide which histogram to fill 

        if (ic == 1) then 

            ! IMPORTANT change: unbias the stored ratios!! 
            ratio = ratio * pSingles

            ! if i ignore the ratios above the upper limit i also can 
            ! only count these excitations if they do not get ignored..

            if (ratio < max_frequency_bound) then

                ! also keep track of the number of done excitations
                if (.not. enough_sing_hist) then 
                    cnt_sing_hist = cnt_sing_hist + 1
                    if (cnt_sing_hist > cnt_threshold) enough_sing_hist = .true.
                end if

                ! find where to put the ratio
                ! since ordered and linear bin bounds i can just divide.. 
                ! i just hope there are no numerical errors creeping in ..
                ! if the ratio happens to be exactly the upper bound, it 
                ! still has to be put in the last bin or? yes i guess.. 
                ! so take the minimum of the ind and the max bin index
                ind = int(ratio / frq_step_size) + 1
                frequency_bins_singles(ind) = frequency_bins_singles(ind) + 1

            end if

        else
            ! check if parallel or anti-parallel 
            if (t_parallel) then 

                ! unbias:
                ratio = ratio * pDoubles * pParallel

                if (ratio < max_frequency_bound) then

                    if (.not. enough_par_hist) then 
                        cnt_par_hist = cnt_par_hist + 1
                        if (cnt_par_hist > 2*cnt_threshold) enough_par_hist = .true.
                    end if

                    ind = int(ratio / frq_step_size) + 1

                    frequency_bins_para(ind) = frequency_bins_para(ind) + 1

                end if
            else 

                ! unbias:
                ratio = ratio * pDoubles * (1.0_dp - pParallel)

                if (ratio < max_frequency_bound) then

                    if (.not. enough_opp_hist) then 
                        cnt_opp_hist = cnt_opp_hist + 1
                        if (cnt_opp_hist > cnt_threshold) enough_opp_hist = .true. 
                    end if

                    ind = int(ratio / frq_step_size) + 1

                    frequency_bins_anti(ind) = frequency_bins_anti(ind) + 1

                end if
            end if
            
            ! combine par and opp into enough_doub
            if (enough_par_hist .and. enough_opp_hist) enough_doub_hist = .true.

        end if

    end subroutine fill_frequency_histogram_4ind

    subroutine fill_frequency_histogram_sd(mat_ele, pgen, ic)
        ! this is the routine, which for now i can use in the symmetry 
        ! adapted GUGA code, where i only distinguish between single and 
        ! double excitation for now.. 
        ! if i adapt that to the already used method in nosym_guga to 
        ! distinguish between more excitaiton i have to write an additional 
        ! subroutine which deals with that 
        ! and this can also be used in the case where we dont differentiate 
        ! between parallel or antiparallel double excitations in the old 
        ! excitation generators
        real(dp), intent(in) :: mat_ele, pgen 
        integer, intent(in) :: ic 
        character(*), parameter :: this_routine = "fill_frequency_histogram_sd"
        integer, parameter :: cnt_threshold = 50

        real(dp) :: ratio
        integer :: ind, new_n_bins, i, old_n_bins
        integer, allocatable :: save_bins(:)

        ASSERT(pgen > EPS) 
        ASSERT( ic == 1 .or. ic == 2)

        if (mat_ele < umateps) return

        ratio = mat_ele / pgen 

        if (ic == 1) then 

            ! unbias: 
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

            ! unbias: 
            ratio = ratio * pDoubles

            if (ratio < max_frequency_bound) then

                if (.not. enough_doub_hist) then 
                    cnt_doub_hist = cnt_doub_hist + 1
                    if (cnt_doub_hist > cnt_threshold) enough_doub_hist = .true.
                end if

                ind = int(ratio / frq_step_size) + 1

                frequency_bins_doubles(ind) = frequency_bins_doubles(ind) + 1

            end if
        end if

    end subroutine fill_frequency_histogram_sd

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
        integer :: ind, new_n_bins, i, old_n_bins 
        integer, allocatable :: save_bins(:) 
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
        integer :: ind, new_n_bins, i, old_n_bins 
        integer, allocatable :: save_bins(:) 
        integer, parameter :: cnt_threshold = 50

        ASSERT(pgen > EPS)
        ASSERT(ic == 1 .or. ic == 2) 
#ifdef __DEBUG 
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
                if (.not.enough_sing_hist) then
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

    subroutine fill_frequency_histogram(mat_ele, pgen)
        ! routine to accumulate the H_ij/pgen ration into the frequency bins
        ! keep parallelism in mind, especially if we have to adjust the 
        ! bin list and boundary list on the fly.. this has to happen on 
        ! all the processors then or? or can i just do it in the end of a loop
        ! adapt this function now, to discriminate between single and double 
        ! excitations, and it pParallel is used, which can be determined 
        ! through the excitation matrix and the 2 involved determinants
        ! (or even the starting determinant) and fill up to 3 frequency 
        ! histograms .. this also means that we have to reset them maybe 
        ! after pSingles or pParallel have been changed since the pgens 
        ! and so the H_ij/pgen ratios change (or? can i modify that somehow 
        ! on the fly?) think about that later and with ali.. 
        ! and for the guga stuff, if i finally implement similar pgens like 
        ! pParallel, as it is already done in the nosym_guga case, i have 
        ! to use even more histograms, but that should be fine.
        real(dp), intent(in) :: mat_ele, pgen 
        character(*), parameter :: this_routine = "fill_frequency_histogram"
        integer, parameter :: cnt_threshold = 50

        real(dp) :: ratio
        integer :: ind, new_n_bins, i, old_n_bins
        integer, allocatable :: save_bins(:)
        ! first have to take correct matrix element, dependent if we use 
        ! complex code or not, or no: for now just input the absolute value 
        ! of H_ij, so its always a real then..
        ! nah.. it puts in 0 mat_eles too.. so just return if 0 mat_ele
        ASSERT(pgen > EPS)

        ! if the matrix element is 0, no excitation will or would be done 
        ! and if the pgen is 0 i also shouldnt be here i guess.. so assert that
        if (mat_ele < EPS) return

        ! then i have to first check if i have to make the histogram bigger...
        ratio = mat_ele / pgen
       
        if (ratio < max_frequency_bound) then
            ! the ratio fits in to the bins now i just have to find the 
            ! correct one
            if (.not. enough_doub_hist) then
                cnt_doub_hist = cnt_doub_hist + 1
                if (cnt_doub_hist > cnt_threshold) enough_doub_hist = .true.
            end if

            ind = int(ratio / frq_step_size) + 1
            
            ! increase counter 
            frequency_bins(ind) = frequency_bins(ind) + 1

        end if

    end subroutine fill_frequency_histogram

    ! write a general histogram communication routine, which takes 
    ! specific histogram as input 
    subroutine comm_frequency_histogram_spec(spec_size, spec_frequency_bins, &
                     all_spec_frq_bins)
        integer, intent(in) :: spec_size
        integer, intent(in) :: spec_frequency_bins(spec_size)
!         integer, allocatable, intent(out) :: all_spec_frq_bins(:)
        integer, intent(out) :: all_spec_frq_bins(n_frequency_bins)
        character(*), parameter :: this_routine = "comm_frequency_histogram_spec"

        integer :: max_size
!         integer, allocatable :: temp_bins(:) 

        ! with the new scheme i do not need to adjust the lengths since they 
        ! are all the same all the time
        all_spec_frq_bins = 0

        call MPIAllReduce(spec_frequency_bins, MPI_SUM, all_spec_frq_bins)
!         
    end subroutine comm_frequency_histogram_spec

    subroutine comm_frequency_histogram(all_frequency_bins)
        ! routine to communicate the frequency histogram data across all 
        ! processors
        integer, intent(out) :: all_frequency_bins(n_frequency_bins)
        character(*), parameter :: this_routine = "comm_frequency_histogram"

        ! do the new implementation! with fixed sizes of the frequency bins
        all_frequency_bins = 0

        call MPIAllReduce(frequency_bins, MPI_SUM, all_frequency_bins)

    end subroutine comm_frequency_histogram

    subroutine integrate_frequency_histogram_spec(spec_size, spec_frequency_bins, &
            ratio)
        ! specific histogram integration routine which sums up the inputted 
        ! frequency_bins
        integer, intent(in) :: spec_size
        integer, intent(in) :: spec_frequency_bins(spec_size) 
        real(dp), intent(out) :: ratio
        character(*), parameter :: this_routine = "integrate_frequency_histogram_spec"

!         integer, allocatable :: all_frequency_bins(:)
        integer :: all_frequency_bins(n_frequency_bins)
        integer :: i, threshold, n_bins
        integer :: n_elements, cnt

        call comm_frequency_histogram_spec(spec_size, spec_frequency_bins, &
            all_frequency_bins)

        n_bins = size(all_frequency_bins)
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

    subroutine integrate_frequency_histogram(ratio)
        ! routine to integrate the entries of the frequency histogram to 
        ! determine the time-step through this. 
        ! use a predefined or inputted threshold, which determines how much 
        ! of the histogram has to be covered to determine the new time-step
        real(dp), intent(out) :: ratio 
        character(*), parameter :: this_routine = "integrate_frequency_histogram"

        integer :: all_frequency_bins(n_frequency_bins)
        integer :: i, threshold, n_bins, n_elements, cnt

        ! have to communicate all the histograms across all cores 

        call comm_frequency_histogram(all_frequency_bins) 

        ! then loop over the histogram and check when the threshold is reached
        n_bins = size(all_frequency_bins) 
        n_elements = sum(all_frequency_bins) 

        if (n_elements == 0) then
            ratio = 0.0_dp
            return

        else if (n_elements < 0) then 
            ! integer overflow -> stop the hist search!
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

        ! (i) determines the boundary which is the ratio then 
        ! use a fixed step-size from the start, so no numericall error 
        ! creeps in.. 
        ratio = i * frq_step_size

    end subroutine integrate_frequency_histogram

end module
