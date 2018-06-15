#include "macros.h"

module guga_tausearch 
    use CalcData, only: gamma_sing, gamma_doub, gamma_two_same, gamma_two_mixed, &
                        gamma_three_same, gamma_three_mixed, gamma_four, &
                        max_death_cpt, enough_sing, enough_doub, enough_two_same, &
                        enough_two_mixed, enough_three_same, enough_three_mixed, &
                        enough_four, tReadPops, tau, MaxWalkerBloom, tTruncInitiator, &
                        InitiatorWalkNo, tWalkContGrow, max_permitted_spawn, &
                        enough_three, enough_two, &
                    frequency_bins_type2, frequency_bounds_type2, frequency_bins_type3, &
                    frequency_bounds_type3, frequency_bins_type4, frequency_bounds_type4, &
                    frequency_bins_type2_diff, frequency_bins_type3_diff, & 
                    frequency_bounds_type2_diff, frequency_bounds_type3_diff, &
                    frequency_bins_singles, frequency_bounds_singles, &
                    t_frequency_analysis, frq_step_size, max_frequency_bound, &
                    n_frequency_bins, t_min_tau, min_tau_global, t_hist_tau_search, &
                    cnt_type2_same, cnt_type2_diff, cnt_type3_same, cnt_type3_diff, &
                    cnt_type4, cnt_sing_hist, cnt_doub_hist, enough_sing_hist, &
                    enough_doub_hist, t_hist_tau_search_option
    use SystemData, only: tUEG, t_consider_diff_bias
    use FciMCData, only: tRestart, pSingles, pDoubles, pExcit2, pExcit4, &
                         pExcit2_same, pExcit3_same, MaxTau, tSearchTau, &
                         tSearchTauOption, tSearchTauDeath

    use constants, only: dp, EPS

    use tau_search, only: FindMaxTauDoubs, integrate_frequency_histogram_spec

    use Parallel_neci

    implicit none
    integer :: cnt_sing, cnt_four, cnt_two_same, cnt_two_mixed, cnt_three_same, &
        cnt_three_mixed

contains
    ! put the previous specifically defined variables in tau_search 
    ! in some general data module 

    subroutine init_tau_search_guga_nosym
        integer :: i

        ! guga version of the tau search routine for the old "nosym" and 
        ! non-weighted excitation generator

        ! the first two are the same as the usual implementation
        gamma_sing = 0.0_dp
        gamma_doub = 0.0_dp

        ! then i have to consider the different general types of 
        ! excitations
        ! for (ii,jj) excitations i could differentiate between the easy
        ! _RR_ -> ^RR^
        ! _LL_ -> ^LL^
        ! and the notorious:
        ! _LR_ -> ^LR^ 
        ! excitation
        gamma_two_same = 0.0_dp
        gamma_two_mixed = 0.0_dp
        ! for the (ii,jk):
        ! here i can differentiate between full-stop and full-start mixed
        ! x -> ^RL^
        ! _RL_ -> x 
        ! which are problematic 
        ! and the "other"
        gamma_three_same = 0.0_dp
        gamma_three_mixed = 0.0_dp
        ! for (ijkl) there is not really reason to differentiate
        gamma_four = 0.0_dp

        ! have to consider deaths too
        max_death_cpt = 0.0_dp

        ! the counts to avoid preemptive updates:
        cnt_sing = 0
      
        enough_sing = .false.
        enough_doub = .false.
        ! check if tau was already set.. -> have to change that routine 
        ! for guga also! does not give correct value yet, is only 
        ! implemented for determinal 4ind-weighted excitation generator
        ! but atleast still runs.. and since tau is updated on the fly then
        ! this shouldn't be too big of a problem..
        if (.not. tRestart .and. (.not. tReadPops) .and. tau < EPS) then
            call FindMaxTauDoubs()
        end if
        write(6,*) "Using initial time-step: ", tau
        write(6,*) "NOTE: this is not yet correctly adapted for the GUGA implementation"
        write(6,*) " -> so use this with caution and check for erroneous values!"

        ! check maximum spawn size: 
        if (MaxWalkerBloom == -1) then
            ! just copy that from tau_search.F90
            if (tTruncInitiator) then
                max_permitted_spawn = InitiatorWalkNo
            else
                max_permitted_spawn = 5.0_dp
            end if
        else
            ! otherwise take inputted value 
            max_permitted_spawn = real(MaxWalkerBloom, dp) 
        end if

        if (.not. (tReadPops .and. .not. tWalkContGrow)) then
            write(iout, "(a,f10.5)") "Will dynamically update timestep to &
                         &limit spawning probability to", max_permitted_spawn
        end if

        if (t_hist_tau_search) call init_hist_tau_search_guga_nosym

    end subroutine init_tau_search_guga_nosym

    subroutine init_hist_tau_search_guga_nosym()


        ! not yet 100% sure about implementation: 
        ! do i really want to distinguish between case(3) e
        ! determine the global and fixed step-size quantitiy! 
        frq_step_size = max_frequency_bound / real(n_frequency_bins, dp)
        ! here i can differentiate between the different types of 
        ! excitations.. 

        cnt_sing_hist = 0
        enough_sing_hist = .false. 

        cnt_doub_hist = 0
        enough_doub_hist = .false.

        cnt_type2_same = 0

        cnt_type3_same = 0 

        cnt_type4 = 0 
        enough_four = .false. 

        enough_two = .false. 
        enough_three = .false.

        allocate(frequency_bins_singles(n_frequency_bins))
        frequency_bins_singles = 0

        ! use the "normal" type2 etc. bins always and just allocate 
        ! additional ones if t_consider_diff_bias 
        allocate(frequency_bins_type2(n_frequency_bins))
        frequency_bins_type2 = 0

        allocate(frequency_bins_type3(n_frequency_bins)) 
        frequency_bins_type3 = 0 

        allocate(frequency_bins_type4(n_frequency_bins))
        frequency_bins_type4 = 0

        if (t_consider_diff_bias) then 

            enough_two_same = .false.
            enough_three_same = .false. 

            cnt_type2_diff = 0 
            enough_two_mixed = .false. 

            cnt_type3_diff = 0
            enough_three_mixed = .false. 

            ! allocate the additional bins and bounds 
            allocate(frequency_bins_type2_diff(n_frequency_bins))
            frequency_bins_type2_diff = 0

            allocate(frequency_bins_type3_diff(n_frequency_bins))
            frequency_bins_type3_diff = 0

        end if

        ! also need to setup all the other quantities necessary for the 
        ! "normal" tau-search if they have not yet been setup if we only 
        ! use the new tau-search 
        if (.not. tSearchTauOption) then 

            ! And what is the maximum death-component found
            max_death_cpt = 0

            ! And the counts are used to make sure we don't update anything too
            ! early
            ! should we use the same variables in both tau-searches?? 

            ! Unless it is already specified, set an initial value for tau
            if (.not. tRestart .and. .not. tReadPops .and. tau == 0) &
                call FindMaxTauDoubs()
            write(6,*) 'Using initial time-step: ', tau
     
            ! Set the maximum spawn size
            if (MaxWalkerBloom == -1) then
                ! No maximum manually specified, so we set the limit of spawn
                ! size to either the initiator criterion, or to 5 otherwise
                if (tTruncInitiator) then
                    max_permitted_spawn = InitiatorWalkNo
                else
                    ! change here to the "old" algorithm, since the time-step
                    ! will be orders of magnitude larger we should limit 
                    ! the max_permitted_spawn to 1. (or 2 maybe.. lets see!)
                    max_permitted_spawn = 1.0_dp
                end if
            else
                ! This is specified manually
                max_permitted_spawn = real(MaxWalkerBloom, dp)
            end if

            if (.not. (tReadPops .and. .not. tWalkContGrow)) then
                write(iout, "(a,f10.5)") "Will dynamically update timestep to &
                             &limit spawning probability to", max_permitted_spawn
            end if

        end if

    end subroutine init_hist_tau_search_guga_nosym

    subroutine log_spawn_magnitude_guga_nosym(ic, ex, matel, pgen)
        ! to allow to use these sort of routines as function pointers, 
        ! dependend on the type of calculation run, also use the 2x2 ex 
        ! matrix, but store the type of guga non-weighted excitation! 
        integer, intent(in) :: ic, ex(2,2)
        real(dp), intent(in) :: pgen, matel
        
        real(dp) :: tmp_gamma, tmp_prob
        integer :: guga_type, same_ind
        integer, parameter :: cnt_threshold = 50
        ! ask simon, how he came to the threshold of 50!?

        ! if i use the hist_tau_search i probably shouldn't do anything in 
        ! here.. 

        if (t_hist_tau_search) return

        guga_type = ex(1,1)
        same_ind = ex(1,2)

        if (ic == 1) then
            ! single excitation
            tmp_prob = pgen / pSingles 
            tmp_gamma = abs(matel) / tmp_prob
            ! this condition picks out, the "worst case" H_ij/pgen relation
            ! -> so, essentially this is why simon means, that the worst case
            ! excitation determines the global FCIQMC time-step..
            ! is there a other way to update the time-step, so this worst
            ! case doesn't influence the time-step so directly?! 
            ! -> talk with ali about that! 
            if (tmp_gamma > gamma_sing) gamma_sing = tmp_gamma

            ! and keep count of excitation (does this have to be reset to 0
            ! at some point?)
            if (.not. enough_sing) then 
                cnt_sing = cnt_sing + 1
                if (cnt_sing > cnt_threshold) enough_sing = .true.
            end if

        else
            ! here i have to determine the type guga-double excitation
            tmp_prob = pgen / pDoubles

            ! for now, consider all type of differentiations! 
            if (guga_type == 2) then
                ! excit_level 2 excitation (ii,jj)
                tmp_prob = tmp_prob / ((1.0_dp - pExcit4) * pExcit2 )

                if (t_consider_diff_bias) then
                    if (same_ind == 1) then
                        ! this means its an excitation without _RL_ or ^RL^
                        tmp_prob = tmp_prob / pExcit2_same 

                        tmp_gamma = abs(matel) / tmp_prob

                        if (tmp_gamma > gamma_two_same) gamma_two_same = tmp_gamma

                        if (.not. enough_two_same) then
                            cnt_two_same = cnt_two_same + 1
                            if (cnt_two_same > cnt_threshold) enough_two_same = .true.

                            ! differentiate here too between the specifics? 
                            ! yes! for now
                            if (enough_two_same .and. enough_two_mixed) enough_two = .true.

                            if (enough_two .and. enough_three .and. enough_four) enough_doub = .true.

    !                         if (all([enough_two_same,enough_two_mixed,enough_three_same,&
    !                             enough_three_mixed,enough_four])) enough_doub = .true.
                        end if
        
                    else
                        ! this means its a _RL_ -> ^RL^ excitation
                        tmp_prob = tmp_prob / (1.0_dp - pExcit2_same) 

                        tmp_gamma = abs(matel) / tmp_prob 

                        if (tmp_gamma > gamma_two_mixed) gamma_two_mixed = tmp_gamma

                        if (.not. enough_two_mixed) then
                            cnt_two_mixed = cnt_two_mixed + 1
                            
                            if (cnt_two_mixed > cnt_threshold) enough_two_mixed = .true.

                            if (enough_two_same .and. enough_two_mixed) enough_two = .true.

                            if (enough_two .and. enough_three .and. enough_four) enough_doub = .true.

                        end if
                    end if

                else
                    ! use the _same variable to store it in this case 
                    tmp_gamma = abs(matel) / tmp_prob

                    if (tmp_gamma > gamma_two_same) gamma_two_same = tmp_gamma

                    if (.not. enough_two) then 
                        cnt_two_same = cnt_two_same + 1 

                        if (cnt_two_same > cnt_threshold) enough_two = .true. 

                        if (enough_two .and. enough_three .and. enough_four) enough_doub = .true. 
                    end if
                end if

            else if (guga_type == 3) then
                ! (ii,jk) type excitations 
                tmp_prob = tmp_prob / ((1.0_dp - pExcit4) * (1.0_dp - pExcit2))

                if (t_consider_diff_bias) then
                    if (same_ind == 1) then
                        ! excitation do NOT involve _RL_ or ^RL^ types
                        tmp_prob = tmp_prob / pExcit3_same

                        tmp_gamma = abs(matel) / tmp_prob 

                        if (tmp_gamma > gamma_three_same) gamma_three_same = tmp_gamma

                        if (.not. enough_three_same) then
                            cnt_three_same = cnt_three_same + 1

                            if (cnt_three_same > cnt_threshold) enough_three_same = .true.

                            if (enough_three_same .and. enough_three_mixed) enough_three = .true.

                            if (enough_two .and. enough_three .and. enough_four) enough_doub = .true.

                        end if
                    else
                        ! its a _RL_ or ^RL^ excitation! 
                        tmp_prob = tmp_prob / (1.0_dp - pExcit3_same) 

                        tmp_gamma = abs(matel) / tmp_prob 

                        if (tmp_gamma > gamma_three_mixed) gamma_three_mixed = tmp_gamma 

                        if (.not. enough_three_mixed) then 
                            cnt_three_mixed = cnt_three_mixed + 1

                            if (cnt_three_mixed > cnt_threshold) enough_three_mixed = .true.

                            if (enough_three_same .and. enough_three_mixed) enough_three = .true. 

                            if (enough_two .and. enough_three .and. enough_four) enough_doub = .true.
                        end if
                    end if
                else 
                    tmp_gamma = abs(matel) / tmp_prob 

                    if (tmp_gamma > gamma_three_same) gamma_three_same = tmp_gamma

                    if (.not. enough_three) then 
                        cnt_three_same = cnt_three_same + 1

                        if (cnt_three_same > cnt_threshold) enough_three = .true.

                        if (enough_two .and. enough_three .and. enough_four) enough_doub = .true.
                    end if
                end if
            else if (guga_type == 4) then 
                ! (ij,kl) excitation! 
                ! no distinction between excitation in here 
                tmp_prob = tmp_prob / pExcit4

                tmp_gamma = abs(matel) / tmp_prob 

                if (tmp_gamma > gamma_four) gamma_four = tmp_gamma 

                if (.not. enough_four) then 
                    cnt_four = cnt_four + 1

                    if (cnt_four > cnt_threshold) enough_four = .true. 

                    if (enough_two .and. enough_three .and. enough_four) enough_doub = .true. 

                end if
            end if
        end if

!         root_print "# of diff. excitations: ", cnt_two_same, cnt_two_mixed, &
!             cnt_three_same, cnt_three_mixed, cnt_four
!         print *, "2 same: ", cnt_two_same
!         print *, "2 mixed:", cnt_two_mixed, enough_two_mixed
!         print *, "3 same: ", cnt_three_same
!         print *, "3 mixed: ", cnt_three_mixed
!         print *, "4 : ", cnt_four

    end subroutine log_spawn_magnitude_guga_nosym

    subroutine update_tau_guga_nosym () 
        ! specialised tau update routine for the guga non-weighted 
        ! excitation generator, which uses no symmetry 
        ! but i have to update that depending if consider diff bias is used 
        ! or not...
        real(dp) :: pSingles_new, tau_new, mpi_tmp, tau_death, pParallel_new
        logical :: mpi_ltmp
        character(*), parameter :: this_routine = "update_tau_guga_nosym"

        real(dp) :: pExcit4_new, pExcit3_same_new, &
                    pExcit2_new, pExcit2_same_new, pBranch2, pBranch3
        real(dp) :: ratio_singles, ratio_type2, ratio_type2_diff, &
                    ratio_type3, ratio_type3_diff, ratio_type4
        real(dp) :: tau_test, psingles_test, pExcit2_same_test, pexcit3_same_test, &
                    pExcit2_test, pexcit4_test, ratio_doubles
        ! from the "old" routine 
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

            ! Condition met --> no need to do this again next iteration
            tSearchTauDeath = .false.
            return

        end if

        ! What needs doing depends on the number of parametrs that are being
        ! updated.
        call MPIAllLORLogical(enough_sing, mpi_ltmp)
        enough_sing = mpi_ltmp
        call MPIAllLORLogical(enough_doub, mpi_ltmp)
        enough_doub = mpi_ltmp

        ! if enough_doub is .true. it also implicates enough_two_same etc. 
        ! are also .true.

        ! Considering two types of double exctitaion...
        ! change that implementation to depent on the input if 
        ! t_consider_diff_bias is true
        call MPIAllReduce (gamma_sing, MPI_MAX, mpi_tmp)
        gamma_sing = mpi_tmp

        call MPIAllReduce (gamma_two_same, MPI_MAX, mpi_tmp)
        gamma_two_same = mpi_tmp

        call MPIAllReduce (gamma_three_same, MPI_MAX, mpi_tmp)
        gamma_three_same = mpi_tmp

        call MPIAllReduce (gamma_four, MPI_MAX, mpi_tmp)
        gamma_four = mpi_tmp

        call MPIAllLORLogical(enough_two_same, mpi_ltmp)
        enough_two_same = mpi_ltmp

        call MPIAllLORLogical(enough_three_same, mpi_ltmp)
        enough_three_same = mpi_ltmp

        call MPIAllLORLogical(enough_four, mpi_ltmp)
        enough_four = mpi_ltmp

        gamma_doub = gamma_two_same + gamma_three_same + gamma_four

        if (t_consider_diff_bias) then
            call MPIAllReduce (gamma_two_mixed, MPI_MAX, mpi_tmp)
            gamma_two_mixed = mpi_tmp

            call MPIAllReduce (gamma_three_mixed, MPI_MAX, mpi_tmp)
            gamma_three_mixed = mpi_tmp

            call MPIAllLORLogical(enough_two_mixed, mpi_ltmp)
            enough_two_mixed = mpi_ltmp

            call MPIAllLORLogical(enough_three_mixed, mpi_ltmp)
            enough_three_mixed = mpi_ltmp

            gamma_doub = gamma_two_same + gamma_two_mixed + gamma_three_same &
                    + gamma_three_mixed + gamma_four
        end if

        ! change: already unbiased in the histogram:
        pBranch2 = pDoubles * (1.0_dp - pExcit4) * pExcit2
        pBranch3 = pDoubles * (1.0_dp - pExcit4) * (1.0_dp - pExcit2)
        
        ! do i want to finetune the different sorts of double excitations
        ! even if there are not enough of the other types of excitations?...

        if (t_consider_diff_bias) then
            if (enough_sing .and. enough_doub) then
                ! this implies all sort of excitations have enough counts..
                pExcit3_same_new = gamma_three_same / (gamma_three_same + gamma_three_mixed)

                pExcit2_same_new = gamma_two_same / (gamma_two_same + gamma_two_mixed)

                pExcit4_new = gamma_four / gamma_doub

                pExcit2_new = (gamma_two_same + gamma_two_mixed) / &
                    (gamma_two_same + gamma_two_mixed + gamma_three_same + gamma_three_mixed) 

                pSingles_new = gamma_sing / (gamma_sing + gamma_doub) 

    !             tau_new = pSingles_new * max_permitted_spawn / gamma_sing
                
                ! this is the same, but gives more insight on what it is 
                tau_new = max_permitted_spawn / ( gamma_sing + gamma_doub )

            else
                ! do the default worst case, according to simons implementation
                pSingles_new = pSingles
                pExcit2_new = pExcit2
                pExcit2_same_new = pExcit2_same
                pExcit3_same_new = pExcit3_same
                pExcit4_new = pExcit4

                tau_new = max_permitted_spawn * min( pSingles / gamma_sing, &
                    pDoubles * pExcit4 / gamma_four, &
                    pBranch2 * pExcit2_same / gamma_two_same, & 
                    pBranch2 * (1.0_dp - pExcit2_same) / gamma_two_mixed, & 
                    pBranch3 * pExcit3_same / gamma_three_same, & 
                    pBranch3 * (1.0_dp - pExcit3_same) / gamma_three_mixed)
            
            end if
        else
            ! do also an implementation without the diff bias 
            if (enough_sing .and. enough_doub) then 
                pExcit4_new = gamma_four / gamma_doub 

                pExcit2_new = gamma_two_same / (gamma_two_same + gamma_three_same) 

                pSingles_new = gamma_sing / (gamma_sing + gamma_doub) 
                
                tau_new = max_permitted_spawn / (gamma_sing + gamma_doub) 

            else 
                pSingles_new = pSingles
                pExcit2_new = pExcit2
                pExcit4_new = pExcit4

                tau_new = max_permitted_spawn * min( & 
                    pSingles / gamma_sing, & 
                    pDoubles * pExcit4 / gamma_four, & 
                    pBranch2 / gamma_two_same) 
            end if
        end if

        ! i only change something if i have enough of all excitations.. 
        ! so i think i do not need to check if there are enough of them 
        ! here 
        if (t_consider_diff_bias) then
            if (abs(pExcit3_same_new - pExcit3_same) / pExcit3_same > 0.0001_dp) then
                root_print "Updating pExcit3_same! new pExcit3_same = ", pExcit3_same_new
            end if

            if (abs(pExcit4_new - pExcit4) / pExcit4 > 0.0001_dp) then 
                root_print "Updating pExcit4! new pExcit4 = ", pExcit4_new
            end if

            if (abs(pExcit2_same_new - pExcit2_same) / pExcit2_same > 0.0001_dp) then
                root_print "Updating pExcit2_same! new pExcit2_same = ", pExcit2_same_new
            end if

            if (abs(pExcit2_new - pExcit2) / pExcit2 > 0.0001_dp) then
                root_print "Updating pExcit2! new pExcit2 = ", pExcit2_new 
            end if
        else 
            if (abs(pExcit4_new - pExcit4) / pExcit4 > 0.0001_dp) then
                root_print "Updating pExcit4! new pExcit4 = ", pExcit4_new
            end if
            if (abs(pExcit2_new - pExcit2) / pExcit2 > 0.0001_dp) then
                root_print "Updating pExcit2! new pExcit2 = ", pExcit2_new 
            end if
        end if

        ! only print out the information if its a major change... 
        ! but doesnt that mean that gradual changes get ignored in the 
        ! output atleast...
        pExcit4 = pExcit4_new
        pExcit2 = pExcit2_new
        if (t_consider_diff_bias) then
            pExcit3_same = pExcit3_same_new
            pExcit2_same = pExcit2_same_new
        end if

        ! simon implemented some hard-limits on single excitation probability
        ! maybe i should think about that too here! 
        
        ! Make sure that we have at least some of both singles and doubles
        ! before we allow ourselves to change the probabilities too much...
        if (enough_sing .and. enough_doub .and. psingles_new > 1e-5_dp &
            .and. psingles_new < (1.0_dp - 1e-5_dp)) then

            if (abs(psingles - psingles_new) / psingles > 0.0001_dp) then
                root_print "Updating singles/doubles bias. pSingles = ", &
                    psingles_new, ", pDoubles = ", 1.0_dp - psingles_new
            end if
            pSingles = psingles_new
            pDoubles = 1.0_dp - pSingles
        end if

        ! in the end the death restriction finally restricts tau.

        ! The range of tau is restricted by particle death. It MUST be <=
        ! the value obtained to restrict the maximum death-factor to 1.0.
        call MPIAllReduce (max_death_cpt, MPI_MAX, mpi_tmp)
        max_death_cpt = mpi_tmp
        tau_death = 1.0_dp / max_death_cpt
        if (tau_death < tau_new) then
            if (t_min_tau) then
                root_print "time-step reduced, due to death events! reset min_tau to:", tau_death
                min_tau_global = tau_death
            end if
            tau_new = tau_death
        end if

        ! And a last sanity check/hard limit
        tau_new = min(tau_new, MaxTau)

        ! If the calculated tau is less than the current tau, we should ALWAYS
        ! update it. Once we have a reasonable sample of excitations, then we
        ! can permit tau to increase if we have started too low.
        if (tau_new < tau .or. ((tUEG .or. enough_sing) .and. enough_doub))then

            ! Make the final tau smaller than tau_new by a small amount
            ! so that we don't get spawns exactly equal to the
            ! initiator threshold, but slightly below it instead.
            tau_new = tau_new * 0.99999_dp

            if (abs(tau - tau_new) / tau > 0.001_dp) then
                if (t_min_tau) then 
                    if (tau_new < min_tau_global) then
                        root_print "new time-step less then min_tau! set to min_tau", min_tau_global

                        tau_new = min_tau_global
                    else
                        root_print "Updating time-step. New time-step = ", tau_new
                    end if
                else
                    root_print "Updating time-step. New time-step = ", tau_new
                end if
            end if
            tau = tau_new
        end if

    end subroutine update_tau_guga_nosym 

    subroutine update_hist_tau_guga_nosym
 
        character(*), parameter :: this_routine = "update_hist_tau_guga_nosym"

        real(dp) :: mpi_tmp, tau_death, ratio_singles, ratio_type2, &
                    ratio_type2_diff, ratio_type3, ratio_type3_diff, ratio_type4, &
                    pExcit3_same_new, pExcit2_same_new, pExcit4_new, pExcit2_new, &
                    pSingles_new, tau_new, ratio_doubles, pBranch2, pBranch3


        logical :: mpi_ltmp

        if (.not. t_hist_tau_search) then 
            ! this means the option was turned on but got turned off due to 
            ! entering var. shift mode, or because the histogramms are full
            ! but the death events should still be considered 
 
            ! Check that the override has actually occurred.
            ASSERT(t_hist_tau_search_option)
            ASSERT(tSearchTauDeath)

            ! The range of tau is restricted by particle death. It MUST be <=
            ! the value obtained to restrict the maximum death-factor to 1.0.
            call MPIAllReduce (max_death_cpt, MPI_MAX, mpi_tmp)
            max_death_cpt = mpi_tmp
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

            ! Condition met --> no need to do this again next iteration
            tSearchTauDeath = .false.

            return
        end if

        ! What needs doing depends on the number of parametrs that are being
        ! updated.
        call MPIAllLORLogical(enough_sing_hist, mpi_ltmp)
        enough_sing_hist = mpi_ltmp
        call MPIAllLORLogical(enough_doub_hist, mpi_ltmp)
        enough_doub_hist = mpi_ltmp

        call integrate_frequency_histogram_spec(size(frequency_bins_singles), &
            frequency_bins_singles, ratio_singles)

        ! if i have a integer overflow i should probably deal here with it..
        if (ratio_singles < 0.0_dp) then 
            ! this means i had an int overflow and should stop the tau-searching
            root_print "The single excitation histogram is full!"
            root_print "stop the hist_tau_search with last time-step: ", tau
            root_print "and pSingles and pDoubles:", pSingles, pDoubles
            t_hist_tau_search = .false. 

            return 
        end if

        ! change: already unbiased in histogram
!         ratio_singles = ratio_singles * pSingles

        pBranch2 = pDoubles * (1.0_dp - pExcit4) * pExcit2
        pBranch3 = pDoubles * (1.0_dp - pExcit4) * (1.0_dp - pExcit2)

        if (t_consider_diff_bias) then 
            ! i have to do that above too.. only check mixed excitaitons 
            ! if i actually use consider_diff_bias 
            call integrate_frequency_histogram_spec(size(frequency_bins_type2), &
                frequency_bins_type2, ratio_type2) 

            ! if i have a integer overflow i should probably deal here with it..
            if (ratio_type2 < 0.0_dp) then 
                ! this means i had an int overflow and should stop the tau-searching
                root_print "The type2_same excitation histogram is full!"
                root_print "stop the hist_tau_search with last time-step: ", tau
                root_print "and pSingles and pDoubles:", pSingles, pDoubles
                t_hist_tau_search = .false. 

                return 
            end if

            call integrate_frequency_histogram_spec(size(frequency_bins_type2_diff), &
                frequency_bins_type2_diff, ratio_type2_diff) 

            ! if i have a integer overflow i should probably deal here with it..
            if (ratio_type2_diff < 0.0_dp) then 
                ! this means i had an int overflow and should stop the tau-searching
                root_print "The type2_diff excitation histogram is full!"
                root_print "stop the hist_tau_search with last time-step: ", tau
                root_print "and pSingles and pDoubles:", pSingles, pDoubles
                t_hist_tau_search = .false. 

                return 
            end if

            call integrate_frequency_histogram_spec(size(frequency_bins_type3), &
                frequency_bins_type3, ratio_type3) 

            ! if i have a integer overflow i should probably deal here with it..
            if (ratio_type3 < 0.0_dp) then 
                ! this means i had an int overflow and should stop the tau-searching
                root_print "The type3_same excitation histogram is full!"
                root_print "stop the hist_tau_search with last time-step: ", tau
                root_print "and pSingles and pDoubles:", pSingles, pDoubles
                t_hist_tau_search = .false. 

                return 
            end if

            call integrate_frequency_histogram_spec(size(frequency_bins_type3_diff), &
                frequency_bins_type3_diff, ratio_type3_diff)

            ! if i have a integer overflow i should probably deal here with it..
            if (ratio_type3_diff < 0.0_dp) then 
                ! this means i had an int overflow and should stop the tau-searching
                root_print "The type_3_diff excitation histogram is full!"
                root_print "stop the hist_tau_search with last time-step: ", tau
                root_print "and pSingles and pDoubles:", pSingles, pDoubles
                t_hist_tau_search = .false. 

                return 
            end if

            call integrate_frequency_histogram_spec(size(frequency_bins_type4), &
                frequency_bins_type4, ratio_type4)

            ! if i have a integer overflow i should probably deal here with it..
            if (ratio_type4 < 0.0_dp) then 
                ! this means i had an int overflow and should stop the tau-searching
                root_print "The type4 excitation histogram is full!"
                root_print "stop the hist_tau_search with last time-step: ", tau
                root_print "and pSingles and pDoubles:", pSingles, pDoubles
                t_hist_tau_search = .false. 

                return 
            end if

            ! unbias already in histograms!
!             ratio_type2 = ratio_type2 * pDoubles * (1.0_dp - pExcit4) * pExcit2 * &
!                 pExcit2_same
! 
!             ratio_type2_diff = ratio_type2_diff * pDoubles * (1.0_dp - pExcit4) * pExcit2 * &
!                 (1.0_dp - pExcit2_same)
! 
!             ratio_type3 = ratio_type3 * pDoubles * (1.0_dp - pExcit4) * &
!                 (1.0_dp - pExcit2) * pExcit3_same
! 
!             ratio_type3_diff = ratio_type3_diff * pDoubles * (1.0_dp - pExcit4) * &
!                 (1.0_dp - pExcit2) * (1.0_dp - pExcit3_same)
! 
!             ratio_type4 = ratio_type4 * pDoubles * pExcit4
! 
            ratio_doubles = ratio_type2 + ratio_type2_diff + ratio_type3 + &
                ratio_type3_diff + ratio_type4

            if (enough_sing_hist .and. enough_doub_hist) then 
                pExcit3_same_new = ratio_type3 / (ratio_type3 + ratio_type3_diff) 
                pExcit2_same_new = ratio_type2 / (ratio_type2 + ratio_type2_diff) 
                pExcit4_new = ratio_type4 / ratio_doubles 
                pExcit2_new = (ratio_type2 + ratio_type2_diff) / & 
                    (ratio_type2 + ratio_type2_diff + ratio_type3 + ratio_type3_diff) 
                pSingles_new = ratio_singles / (ratio_singles + ratio_doubles) 
                
                tau_new = max_permitted_spawn / (ratio_singles + ratio_doubles) 
                if (pSingles_new > 1e-5_dp .and. &
                    pSingles_new < (1.0_dp - 1e-5_dp)) then 

                    root_print "Updating singles/doubles bias. pSingles = ", &
                        psingles_new, ", pDoubles = ", 1.0_dp - psingles_new

                    pSingles = pSingles_new
                    pDoubles = 1.0_dp - pSingles

                end if

                ! i should also set some boundaries so these quantitities 
                ! do not get incredibly low or very close to 1...
                ! 10% maybe?? less more?! check in runs!
                if (pExcit3_same_new > 1e-1_dp .and. pExcit3_same_new < (1.0_dp - 1e-1_dp)) then
                    if (abs(pExcit3_same_new - pExcit3_same) / pExcit3_same > 0.0001_dp) then
                        root_print "new pExcit3_same_new: ", pExcit3_same_new
                    end if
                    pExcit3_same = pExcit3_same_new
                end if

                if (pExcit2_same_new > 1e-1_dp .and. pExcit2_same_new < (1.0_dp - 1e-1_dp)) then
                    if (abs(pExcit2_same_new - pExcit2_same) / pExcit2_same > 0.0001_dp) then
                        root_print "new pExcit2_same_new: ", pExcit2_same_new
                    end if
                    pExcit2_same = pExcit2_same_new
                end if

                if (pExcit4_new > 1e-1_dp .and. pExcit4_new < (1.0_dp - 1e-1_dp)) then
                    if (abs(pExcit4_new - pExcit4) / pExcit4 > 0.0001_dp) then
                        root_print "new pExcit4_new: ", pExcit4_new
                    end if
                    pExcit4 = pExcit4_new
                end if

                if (pExcit2_new > 1e-1_dp .and. pExcit2_new < (1.0_dp - 1e-1_dp)) then 
                    if (abs(pExcit2_new - pExcit2) / pExcit2 > 0.0001_dp) then
                        root_print "new pExcit2_new: ", pExcit2_new
                    end if
                    pExcit2 = pExcit2_new
                end if

            else 
                tau_new = max_permitted_spawn * min(&
                    pSingles / ratio_singles, &
                    pDoubles * pExcit4 / ratio_type4, & 
                    pBranch2 * pExcit2_same / ratio_type2, & 
                    pBranch2 * (1.0_dp - pExcit2_same) / ratio_type2_diff, &
                    pBranch3 * pExcit3_same / ratio_type3, &
                    pBranch3 * (1.0_dp - pExcit3_same) / ratio_type3_diff)

            end if


        else 
            ! no differentiating between mixed and alike type 2 and 3 
            ! excitations 
            call integrate_frequency_histogram_spec(size(frequency_bins_type2), &
                frequency_bins_type2, ratio_type2) 

            ! if i have a integer overflow i should probably deal here with it..
            if (ratio_type2 < 0.0_dp) then 
                ! this means i had an int overflow and should stop the tau-searching
                root_print "The type2_same excitation histogram is full!"
                root_print "stop the hist_tau_search with last time-step: ", tau
                root_print "and pSingles and pDoubles:", pSingles, pDoubles
                t_hist_tau_search = .false. 

                return 
            end if

            call integrate_frequency_histogram_spec(size(frequency_bins_type3), &
                frequency_bins_type3, ratio_type3) 

            ! if i have a integer overflow i should probably deal here with it..
            if (ratio_type3 < 0.0_dp) then 
                ! this means i had an int overflow and should stop the tau-searching
                root_print "The type3_same excitation histogram is full!"
                root_print "stop the hist_tau_search with last time-step: ", tau
                root_print "and pSingles and pDoubles:", pSingles, pDoubles
                t_hist_tau_search = .false. 

                return 
            end if

            call integrate_frequency_histogram_spec(size(frequency_bins_type4), &
                frequency_bins_type4, ratio_type4)

            ! if i have a integer overflow i should probably deal here with it..
            if (ratio_type4 < 0.0_dp) then 
                ! this means i had an int overflow and should stop the tau-searching
                root_print "The type4 excitation histogram is full!"
                root_print "stop the hist_tau_search with last time-step: ", tau
                root_print "and pSingles and pDoubles:", pSingles, pDoubles
                t_hist_tau_search = .false. 

                return 
            end if

            ! unbias already in the histograms!
!             ratio_type2 = ratio_type2 * pDoubles * (1.0_dp - pExcit4) * pExcit2 
! 
!             ratio_type3 = ratio_type3 * pDoubles * (1.0_dp - pExcit4) * &
!                 (1.0_dp - pExcit2)
! 
!             ratio_type4 = ratio_type4 * pDoubles * pExcit4

            ratio_doubles = ratio_type2 + ratio_type3 + ratio_type4

            ! hm.. i have to change the pgen calc. again in the GUGA 
            ! approach.. i should not include pSingles, pDoubles etc. 
            ! within the excitation generation routines, but only apply 
            ! it afterwards, in the outer loops, so i can change 
            ! pDoubles etc. outside without changing the ratio 
            ! mat_ele / pgen ... so i dont have to change anything 
            ! in the histograms.. todo! now! 
            ! na des passt schon.. ich dividiers halt hier wieder raus 
            ! um die "richtige" wahrscheinlichkeit zu kriegen.. 

            if (enough_sing_hist .and. enough_doub_hist) then 
                pExcit4_new = ratio_type4 / ratio_doubles 
                pExcit2_new = ratio_type2 / (ratio_type2 + ratio_type3) 
                pSingles_new = ratio_singles / (ratio_singles + ratio_doubles) 
                
                tau_new = max_permitted_spawn / (ratio_singles + ratio_doubles) 

                if (pSingles_new > 1e-5_dp .and. &
                    pSingles_new < (1.0_dp - 1e-5_dp)) then 
                    root_print "new psingles test: ", pSingles_new

                    pSingles = pSingles_new
                    pDoubles = 1.0_dp - pSingles

                end if

                if (pExcit4_new > 1e-1_dp .and. pExcit4_new < (1.0_dp - 1e-1_dp)) then
                    if (abs(pExcit4_new - pExcit4) / pExcit4 > 0.0001_dp) then
                        root_print "new pExcit4_new: ", pExcit4_new
                    end if
                    pExcit4 = pExcit4_new
                end if

                if (pExcit2_new > 1e-1_dp .and. pExcit2_new < (1.0_dp - 1e-1_dp)) then
                    if (abs(pExcit2_new - pExcit2) / pExcit2 > 0.0001_dp) then
                        root_print "new pExcit2_new: ", pExcit2_new
                    end if
                    pExcit2 = pExcit2_new
                end if

            else 
                tau_new = max_permitted_spawn * min(&
                    pSingles / ratio_singles, &
                    pDoubles * pExcit4 / ratio_type4, & 
                    pBranch2 / ratio_type2)

            end if
        end if

        ! to deatch check again and finally update time-step
        ! The range of tau is restricted by particle death. It MUST be <=
        ! the value obtained to restrict the maximum death-factor to 1.0.
        call MPIAllReduce (max_death_cpt, MPI_MAX, mpi_tmp)
        max_death_cpt = mpi_tmp
        tau_death = 1.0_dp / max_death_cpt

        ! and have to carefully check if i have enough excitations of all sorts
        ! before i individually update the probabilities
        if (tau_death < tau_new) then
            if (t_min_tau) then
                root_print "time-step reduced, due to death events! reset min_tau to:", tau_death
                min_tau_global = tau_death
            end if
            tau_new = tau_death
        end if

        ! And a last sanity check/hard limit
        tau_new = min(tau_new, MaxTau)

        ! If the calculated tau is less than the current tau, we should ALWAYS
        ! update it. Once we have a reasonable sample of excitations, then we
        ! can permit tau to increase if we have started too low.
        if (tau_new < tau .or. ((tUEG .or. enough_sing) .and. enough_doub))then

            ! Make the final tau smaller than tau_new by a small amount
            ! so that we don't get spawns exactly equal to the
            ! initiator threshold, but slightly below it instead.
!             tau_new = tau_new * 0.99999_dp

            if (abs(tau - tau_new) / tau > 0.0001_dp) then
                if (t_min_tau) then
                    if (tau_new < min_tau_global) then
                        root_print "new time-step less then min_tau! set to min_tau!", min_tau_global

                        tau_new = min_tau_global

                    else 
                        root_print "Updating time-step. New time-step = ", tau_new
                    end if
                else
                    root_print "Updating time-step. New time-step = ", tau_new

                end if
            end if
            tau = tau_new
        end if

    end subroutine update_hist_tau_guga_nosym

end module
