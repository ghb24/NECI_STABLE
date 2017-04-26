#include "macros.h"

module tau_search_hist

    use SystemData, only: tGen_4ind_weighted, AB_hole_pairs, par_hole_pairs,tHub, & 
                          tGen_4ind_reverse, nOccAlpha, nOccBeta, tUEG, tGen_4ind_2, &
                          UMatEps
    use CalcData, only: tTruncInitiator, tReadPops, MaxWalkerBloom, tau, &
                        InitiatorWalkNo, tWalkContGrow, &                
                        t_min_tau, min_tau_global, & 
                        max_frequency_bound, n_frequency_bins, &
                        frq_ratio_cutoff, t_hist_tau_search,  &
                        t_fill_frequency_hists, t_hist_tau_search_option, &
                        t_truncate_spawns, t_mix_ratios, mix_ratio
    use FciMCData, only: tRestart, pSingles, pDoubles, pParallel, &
                         MaxTau, tSearchTau, tSearchTauOption, tSearchTauDeath
    use Parallel_neci, only: MPIAllReduce, MPI_MAX, MPI_SUM, MPIAllLORLogical
    use ParallelHelper, only: iprocindex
    use constants, only: dp, EPS, iout
    use tau_search, only: FindMaxTauDoubs
    use MemoryManager, only: LogMemAlloc, LogMemDealloc, TagIntType

    use procedure_pointers, only: get_umat_el
    use UMatCache, only: gtid, UMat2d
    use util_mod, only: abs_l1
    implicit none
    ! variables which i might have to define differently:
    logical :: consider_par_bias
    real(dp) :: max_permitted_spawn, max_death_cpt
    integer :: n_opp, n_par
    ! do i have to define this here or in the CalcData:??
    integer :: cnt_sing_hist, cnt_doub_hist, cnt_opp_hist, cnt_par_hist
    integer :: above_max_singles, above_max_para, above_max_anti, above_max_doubles
    logical :: enough_sing_hist, enough_doub_hist, enough_par_hist, enough_opp_hist

    ! store the necessary quantities here and not in CalcData! 
    real(dp) :: frq_step_size = 0.1_dp

    ! i need bin arrays for all types of possible spawns: 
    integer, allocatable :: frequency_bins_singles(:), frequency_bins_para(:), &
                            frequency_bins_anti(:), frequency_bins_doubles(:), &
                            frequency_bins(:)

    integer(TagIntType) :: mem_tag_histograms = 0

contains

    subroutine init_hist_tau_search 
        ! split the new and old tau-search routines up, to no mix up 
        ! too much stuff 
        character(*), parameter :: this_routine = "init_hist_tau_search"

        integer :: ierr

        ! at the beginning to some inut checking: 
        if (tSearchTauOption .or. tSearchTau) then 
            ! is it already too late here? maybe so better stop! but for now 
            ! try to turn it off! 
            write(iout, '("WARNING: standard and histogramming tau-search are&
                & chosen! TURNING OFF STANDARD TAU-SEARCH!")')
            tSearchTauOption = .false.
            tSearchTau = .false. 
        end if

        ! if no truncating spawns is chosen warn here againg that that might 
        ! cause problems! 
        if (.not. t_truncate_spawns) then 
            write(iout, '("WARNING: NO spawn truncation chosen with keyword: &
                &truncate-spawns [float] in input! this might cause &
                &bloom problems with histogramming tau-search! BE CAUTIOUS!")')
        end if

        ! i have to check if i am restarting or i do a fresh calculation:
        if (tReadPops) then 
            ! TODO i have to figure out what exactly to do in this case! 
            ! maybe use Pablos idea of storing histograms always and checking
            ! here if a histogram is present and then calculating the 
            ! time-step and psingles etc. from it! 
            ! and also output the read-in or calculated quantities here! 
        end if

        ! print out the standard quantities: 
        print *, "Setup of the Histogramming tau-search: "
        print *, "  Integration cut-off: ", frq_ratio_cutoff
        print *, "  Number of bins: ", n_frequency_bins
        print *, "  Max. ratio: ", max_frequency_bound

        ! do the initialization of the frequency analysis here.. 
        ! i think otherwise it is not done on all the nodes.. 
        ! don't need to check if t_frequency_analysis here anymore
        ! determine the global and fixed step-size quantitiy! 
        frq_step_size = max_frequency_bound / real(n_frequency_bins, dp)

        print *, "  Bin-width: ", frq_step_size

        ! and do the rest of the initialisation:

        ! Are we considering parallel-spin bias in the doubles?
        ! Do this logic here, so that if we add opposite spin bias to more
        ! excitation generators, then there is only one place that this logic
        ! needs to be updated!
        if (tGen_4ind_weighted .or. tGen_4ind_2) then
            consider_par_bias = .true.
        else if (tGen_4ind_reverse) then
            consider_par_bias = .true.
            n_opp = AB_hole_pairs
            n_par = par_hole_pairs
        else
            consider_par_bias = .false.
        end if

        ! If there are only a few electrons in the system, then this has
        ! impacts for the choices that can be made.
        if (nOccAlpha == 0 .or. nOccBeta == 0) then
            ! do i really need a histogramming tau-search in this case?? 
            ! come on..
            call stop_all(this_routine, &
                "Do you really need a tau-search for such a small system?")
            consider_par_bias = .false.
            pParallel = 1.0_dp
        end if
        if (nOccAlpha == 1 .and. nOccBeta == 1) then
            consider_par_bias = .false.
            call stop_all(this_routine, &
                "Do you really need a tau-search for such a small system?")
            pParallel = 0.0_dp
        end if

        cnt_sing_hist = 0
        enough_sing_hist = .false.
        above_max_singles = 0

        cnt_doub_hist = 0 
        enough_doub_hist = .false.
        above_max_doubles = 0
        ! dependent if we use pParallel or not init specific hists
        if (consider_par_bias) then
            
            cnt_par_hist = 0 
            cnt_opp_hist = 0

            enough_par_hist = .false. 
            enough_opp_hist = .false. 

            above_max_para = 0
            above_max_anti = 0

            ! i always use the singles histogram dont I? i think so.. 
            ! Log the memory here! TODO
            allocate(frequency_bins_singles(n_frequency_bins), stat = ierr)
            frequency_bins_singles = 0

            allocate(frequency_bins_para(n_frequency_bins), stat = ierr)
            frequency_bins_para = 0

            allocate(frequency_bins_anti(n_frequency_bins), stat = ierr)
            frequency_bins_anti = 0

            call LogMemAlloc('frequency_bins', n_frequency_bins * 3, 4, &
                this_routine, mem_tag_histograms, ierr)

        else 
            if (tHub .or. tUEG) then
                ! only one histogram is used! 
                allocate(frequency_bins(n_frequency_bins), stat = ierr)

                call LogMemAlloc('frequency_bins', n_frequency_bins, 4, &
                    this_routine, mem_tag_histograms, ierr)

            else

                ! i always use the singles histogram dont I? i think so.. 
                allocate(frequency_bins_singles(n_frequency_bins), stat = ierr)
                frequency_bins_singles = 0

                ! for now use only pSingles and pDoubles for GUGA implo
                allocate(frequency_bins_doubles(n_frequency_bins), stat = ierr)
                frequency_bins_doubles = 0

                call LogMemAlloc('frequency_bins', n_frequency_bins * 2, 4, &
                    this_routine, mem_tag_histograms, ierr)
            end if
        end if

        ! also need to setup all the other quantities necessary for the 
        ! "normal" tau-search if they have not yet been setup if we only 
        ! use the new tau-search 

        ! And what is the maximum death-component found
        max_death_cpt = 0

        ! And the counts are used to make sure we don't update anything too
        ! early
        ! should we use the same variables in both tau-searches?? 

        ! Unless it is already specified, set an initial value for tau
        if (.not. tRestart .and. .not. tReadPops .and. tau < EPS) then
            call FindMaxTauDoubs()
        end if
        if (tReadPops) then
            write(iout,*) "Using time-step from POPSFILE!"
        else
            write(6,*) 'Using initial time-step: ', tau
        end if
 
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

    end subroutine init_hist_tau_search


    subroutine update_tau_hist() 
        ! split up the new tau-update routine from the old one! 
        character(*), parameter :: this_routine = "update_tau_hist"

        real(dp) :: psingles_new, tau_new, mpi_tmp, tau_death, pParallel_new
        real(dp) :: ratio_singles, ratio_anti, ratio_para, ratio_doubles, ratio
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

        ! singles is always used.. 
        ! thats not quite right.. for the hubbard/UEG case it is not.. 
        if (tUEG .or. tHub) then
            ! here i could implement the summing in the case of Hubbard 
            ! and UEG models.. although I could "just" implement the 
            ! optimal time-step in the case of Hubbard models! 
            if (tHub) then
                call stop_all(this_routine, &
                    "in the case of the Hubbard model there is a optimal time-step &
                    &analytically calculatable! So do this!")
            end if
            ! for UEG not i guess.. 
            call integrate_frequency_histogram_spec(frequency_bins, ratio)

            if (ratio < 0.0_dp) then 
                ! TODO: and also print out in this case! 
                ! this means i had an int overflow and should stop the tau-searching
                root_print "The single excitation histogram is full!"
                root_print "stop the hist_tau_search with last time-step: ", tau
                t_hist_tau_search = .false. 

                return
            end if

            if (enough_doub_hist) then 
                ! in this case the enough_doub_hist flag is (mis)used to 
                ! indicate enough spawning events! 
                tau_new = max_permitted_spawn / ratio

                ! and use the mixing now: 
                if (t_mix_ratios) then
                    tau_new = (1.0_dp - mix_ratio) * tau + mix_ratio * tau_new
                end if

            else 
                ! what do i do in the case if not enough spawns.. 
                ! nothing i guess.. 
                tau_new = tau
            end if

        else

            call integrate_frequency_histogram_spec(frequency_bins_singles, &
                ratio_singles) 

            ! for analyis print this also out: 
#ifdef __DEBUG
            print *, "ratio_singles: ", ratio_singles
#endif
            ! if i have a integer overflow i should probably deal here with it..
            if (ratio_singles < 0.0_dp) then 
                ! TODO: and also print out in this case! 
                ! this means i had an int overflow and should stop the tau-searching
                root_print "The single excitation histogram is full!"
                root_print "stop the hist_tau_search with last time-step: ", tau
                root_print "and pSingles and pDoubles:", pSingles, pDoubles
                t_hist_tau_search = .false. 

                return 
            end if

            ! change the storage of the ratios so that they are already 
            ! unbiased!
!             ratio_singles = ratio_singles * pSingles

            if (tGen_4ind_weighted .or. tGen_4ind_2 .or. tGen_4ind_reverse) then 

                ASSERT(consider_par_bias)

                ! sum up all 3 ratios: single, parallel and anti-parallel
                call integrate_frequency_histogram_spec(frequency_bins_para, &
                    ratio_para)
#ifdef __DEBUG
                print *, " ratio_para: ", ratio_para
#endif

                if (ratio_para < 0.0_dp) then 
                    root_print "The parallel excitation histogram is full!" 
                    root_print "stop the hist_tau_search with last time-step: ", tau
                    root_print "and pSingles and pParallel: ", pSingles, pParallel

                    t_hist_tau_search = .false. 

                    return
                end if

                call integrate_frequency_histogram_spec(frequency_bins_anti, &
                    ratio_anti)
#ifdef __DEBUG
                print *, "ratio_anti: ", ratio_anti
#endif

                if (ratio_anti < 0.0_dp) then 
                    root_print "The anti-parallel excitation histogram is full!" 
                    root_print "stop the hist_tau_search with last time-step: ", tau
                    root_print "and pSingles and pParallel: ", pSingles, pParallel

                    t_hist_tau_search = .false. 

                    return
                end if

                ! to compare the influences on the time-step:
                ! change that the ratios are already stored in an unbiased 
                ! way, so no feedback happens!
!                 ratio_para = ratio_para * pDoubles * pParallel
!                 ratio_anti = ratio_anti * pDoubles * (1.0_dp - pParallel)

                ! also calculate new time-step through this method and 
                ! check the difference to the old method 
                if (enough_sing_hist .and. enough_doub_hist) then 
                    ! puh.. for the mixing i am not sure how to exactly do 
                    ! that.. since all of the ratios depend on each other..
                    ! test that!
                    pparallel_new = ratio_para / (ratio_anti + ratio_para)
                    if (t_mix_ratios) then
                        pparallel_new = (1.0_dp - mix_ratio) * pParallel + &
                                        mix_ratio * pparallel_new
                    end if

                    psingles_new = ratio_singles * pparallel_new / &
                        (ratio_para + ratio_singles * pparallel_new) 

                    if (t_mix_ratios) then
                        psingles_new = (1.0_dp - mix_ratio) * psingles + & 
                                        mix_ratio * psingles_new
                    end if 

                    tau_new = psingles_new * max_permitted_spawn / ratio_singles

                    if (t_mix_ratios) then
                        tau_new = (1.0_dp - mix_ratio) * tau + mix_ratio * tau_new
                    end if

                    if (psingles_new > 1e-5_dp .and. &
                        psingles_new < (1.0_dp - 1e-5_dp)) then

                        root_print "Updating singles/doubles bias. pSingles = ", &
                            psingles_new, ", pDoubles = ", 1.0_dp - psingles_new, "in: ", this_routine

                        pSingles = psingles_new
                        pDoubles = 1.0_dp - pSingles
                    end if

                    ! although checking for enough doubles is probably more efficient 
                    ! than always checking the reals below..
                    if (pParallel_new  > 1e-1_dp .and. pParallel_new < (1.0_dp - 1e-1_dp)) then
                        ! enough_doub implies that both enough_opp and enough_par are 
                        ! true.. so this if statement makes no sense 
                        ! and otherwise pParallel_new is the same as before
                        if (abs(pparallel_new-pParallel) / pParallel > 0.0001_dp) then
                            root_print "Updating parallel-spin bias; new pParallel = ", &
                                pParallel_new, "in: ", this_routine
                        end if
                        ! in this new implementation the weighting make pParallel 
                        ! smaller and smaller.. so also limit it to some lower bound
                            pParallel = pParallel_new
                    end if

                else
                    pparallel_new = pParallel
                    psingles_new = pSingles
                    ! im not sure if this possible division by 0 is cool...
                    ! and if this is somehow cool in the histogramming 
                    ! approach.. this again is some kind of worst case 
                    ! adaptation.. hm.. todo
                    tau_new = max_permitted_spawn * min(&
                        pSingles / ratio_singles, &
                        pDoubles * pParallel / ratio_para, &
                        pDoubles * (1.0_dp - pParallel) / ratio_anti)

                end if

            else 
                ! here i only use doubles for now
                call integrate_frequency_histogram_spec(frequency_bins_doubles, &
                    ratio_doubles)

                if (ratio_doubles < 0.0_dp) then 
                    root_print "The double excitation histogram is full!" 
                    root_print "stop the hist_tau_search with last time-step: ", tau
                    root_print "and pSingles and pDoubles: ", pSingles, pDoubles

                    t_hist_tau_search = .false. 

                    return
                end if

#ifdef __DEBUG
                print *, "ratio_doubles: ", ratio_doubles
#endif
                ! to compare the influences on the time-step:
                ! change that so it does store the ratio unbiased
!                 ratio_doubles = ratio_doubles * pDoubles

                if (enough_sing_hist .and. enough_doub_hist) then 
                    psingles_new = ratio_singles / (ratio_doubles + ratio_singles)

                    if (t_mix_ratios) then
                        psingles_new = (1.0_dp - mix_ratio) * pSingles + &
                                        mix_ratio * psingles_new
                    end if

                    tau_new = max_permitted_spawn / (ratio_doubles + ratio_singles)

                    if (t_mix_ratios) then
                        tau_new = (1.0_dp - mix_ratio) * tau + mix_ratio * tau_new
                    end if

                    if (psingles_new > 1e-5_dp .and. &
                        psingles_new < (1.0_dp - 1e-5_dp)) then

                        if (abs(psingles - psingles_new) / psingles > 0.0001_dp) then
                            root_print "Updating singles/doubles bias. pSingles = ", &
                                psingles_new, ", pDoubles = ", 1.0_dp - psingles_new, "in: ", this_routine
                        end if

                        pSingles = psingles_new
                        pDoubles = 1.0_dp - pSingles
                    end if

                else 
                    psingles_new = pSingles
                    tau_new = max_permitted_spawn * & 
                        min(pSingles / ratio_singles, pDoubles / ratio_doubles)

                end if
            end if
        end if

        ! to deatch check again and finally update time-step
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
        if (tau_new < tau .or. (((tUEG .or. tHub) .or. enough_sing_hist) .and. enough_doub_hist))then

            ! Make the final tau smaller than tau_new by a small amount
            ! so that we don't get spawns exactly equal to the
            ! initiator threshold, but slightly below it instead.
            ! does this make sense in the new implmentation? this way 
            ! i will always decrease the time-step even if its not necessary.. 
!             tau_new = tau_new * 0.99999_dp

            ! also does the restriction on the output make sense? since i am 
            ! always changing it anyway... atleast make it smaller..
            if (abs(tau - tau_new) / tau > 0.0001_dp) then
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

    end subroutine update_tau_hist

   
    subroutine fill_frequency_histogram_4ind(mat_ele, pgen, ic, t_parallel, ex)
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
        integer, intent(in), optional :: ex(2,2)
        character(*), parameter :: this_routine = "fill_frequency_histogram_4ind"
        ! i think in my histogramming tau-search i have to increase this 
        ! threshold by a LOT to avoid fluctuating behavior at the 
        ! beginning of the calculation
        integer, parameter :: cnt_threshold = 50

        real(dp) :: ratio
        integer :: ind
        integer :: indi,indj,inda,indb

        ASSERT(pgen > EPS) 
        ASSERT(ic == 1 .or. ic == 2)

        ! i can't make this exception here without changing alot in the 
        ! other parts of the code.. and i have to talk to ali first about 
        ! that
        if (mat_ele < EPS) then
#ifdef __DEBUG
            print *, "zero matele should not be here!"
            print *, "mat_ele: ", mat_ele
            print *, "pgen: ", pgen 
            print *, "ic: ", ic
            print *, "parallel: ", t_parallel
            print *, "ex-matrix: ", ex
#endif
            return 
        end if

#ifdef __DEBUG
        if (pgen < EPS) then
            print *, "zero pgen! should not be here!" 
            print *, "mat_ele: ", mat_ele
            print *, "pgen: ", pgen 
            print *, "ic: ", ic
            print *, "parallel? ", t_parallel
            print *, "ex-matrix: ", ex
            
            ! i think i have to stop all.. since walkers will explode.. 
            call stop_all(this_routine, "zero pgen! how did this happen?")

        end if

        if (isnan(pgen)) then
            print *, "pgen is nan for some reason.." 
            print *, "mat_ele: ", mat_ele
            print *, "pgen: ", pgen 
            print *, "ic: ", ic
            print *, "parallel? ", t_parallel
            print *, "ex-matrix: ", ex

            call stop_all(this_routine, "pgen is nan! what happened?")
        end if
#endif

        ratio = mat_ele / pgen

        ! then i have to decide which histogram to fill 

        if (ic == 1) then 

            ! if i ignore the ratios above the upper limit i also can 
            ! only count these excitations if they do not get ignored..

            ! i have to change the way how we store those excitations in the 
            ! histograms! i have to unbias it against the psingles, etc. 
            ! quantities, so the histograms do not have feedback! 

            ratio = ratio * pSingles

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

            else
                ! store the number of excitation which exceed the upper limit!
                above_max_singles = above_max_singles + 1
                print *, "Warning: single excitation H_ij/pgen above max_frequency_bound!" 
                print *, " H_ij: ", mat_ele, ", pgen: ", pgen, ", pSingles: ", pSingles
                print *, " excitation-matrix: ", ex
                print *, " H_ij/pgen: ", ratio, " ; bound: ", max_frequency_bound
                print *, " Consider increasing the bound!"

            end if

        else
            ! check if parallel or anti-parallel 
            if (t_parallel) then 

                ratio = ratio * (pDoubles * pParallel)
#ifdef __DEBUG
                ! analyse the really low and really high ratios: 
                if (ratio < 0.001_dp) then 
                    print *, "******************"
                    print *, "parallel excitation:"
                    print *, "ratio: ", ratio
                    print *, "mat_ele: ", mat_ele
                    print *, "pgen: ", pgen
                    print *, "ex-maxtrix: ", gtid(ex)
                    indi = gtid(ex(1,1))
                    indj = gtid(ex(1,2))
                    inda = gtid(ex(2,1))
                    indb = gtid(ex(2,2))
                    print *, "umat (ij|ab) ", get_umat_el(indi,indj,inda,indb)
                    print *, "umat (ij|ba) ", get_umat_el(indi,indj,indb,inda)
                    print *, "diff: ", abs(get_umat_el(indi,indj,inda,indb) - &
                        get_umat_el(indi,indj,indb,inda))
                    print *, "(ii|aa):", abs_l1(UMat2d(max(indi,inda),min(indi,inda)))
                    print *, "(jj|aa):", abs_l1(UMat2d(max(indj,inda),min(indj,inda)))
                    print *, "(ii|bb): ",abs_l1(UMat2d(max(indi,indb),min(indi,indb)))
                    print *, "(jj|bb): ", abs_l1(UMat2d(max(indj,indb),min(indj,indb)))
                    print *, "******************"

                end if
#endif
               if (ratio < max_frequency_bound) then
                    if (.not. enough_par_hist) then 
                        cnt_par_hist = cnt_par_hist + 1
                        if (cnt_par_hist > 2*cnt_threshold) enough_par_hist = .true.
                    end if

                    ind = int(ratio / frq_step_size) + 1

                    frequency_bins_para(ind) = frequency_bins_para(ind) + 1
                else
                    above_max_para = above_max_para + 1
                    print *, "Warning: parallel excitation H_ij/pgen above max_frequency_bound!" 
                    print *, " H_ij/pgen: ", ratio, " ; bound: ", max_frequency_bound
                    print *, " Consider increasing the bound!"
                end if
            else 

                ratio = ratio * (pDoubles * (1.0_dp - pParallel))
                ! analyse the really low and really high ratios: 
#ifdef __DEBUG
                if (ratio < 0.001_dp) then 
                    print *, "******************"
                    print *, "anti-parallel excitation:"
                    print *, "ratio: ", ratio
                    print *, "mat_ele: ", mat_ele
                    print *, "pgen: ", pgen
                    print *, "ex-maxtrix: ", ex
                    indi = gtid(ex(1,1))
                    indj = gtid(ex(1,2))
                    inda = gtid(ex(2,1))
                    indb = gtid(ex(2,2))
                    print *, "umat (ij|ab) ", get_umat_el(indi,indj,inda,indb)
                    print *, "umat (ij|ba) ", get_umat_el(indi,indj,indb,inda)
                    print *, "diff: ", abs(get_umat_el(indi,indj,inda,indb) - &
                        get_umat_el(indi,indj,indb,inda))
                    print *, "(ii|aa):", abs_l1(UMat2d(max(indi,inda),min(indi,inda)))
                    print *, "(jj|aa):", abs_l1(UMat2d(max(indj,inda),min(indj,inda)))
                    print *, "(ii|bb): ",abs_l1(UMat2d(max(indi,indb),min(indi,indb)))
                    print *, "(jj|bb): ", abs_l1(UMat2d(max(indj,indb),min(indj,indb)))
                    print *, "******************"
                end if
#endif

                if (ratio < max_frequency_bound) then
                    if (.not. enough_opp_hist) then 
                        cnt_opp_hist = cnt_opp_hist + 1
                        if (cnt_opp_hist > cnt_threshold) enough_opp_hist = .true. 
                    end if

                    ind = int(ratio / frq_step_size) + 1

                    frequency_bins_anti(ind) = frequency_bins_anti(ind) + 1

                else 
                    above_max_anti = above_max_anti + 1
                    print *, "Warning: anti-parallel excitation H_ij/pgen above max_frequency_bound!" 
                    print *, " H_ij/pgen: ", ratio, " ; bound: ", max_frequency_bound
                    print *, " Consider increasing the bound!"
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
        integer :: ind

        ASSERT(pgen > EPS) 
        ASSERT( ic == 1 .or. ic == 2)

        if (mat_ele < EPS) return

        if (pgen < EPS) then
            ! something went wrong then or?? but what exactly? 
            call stop_all(this_routine, "pgen is zero! something went wrong!")

        end if

        ratio = mat_ele / pgen 

        if (ic == 1) then 

            ! change to unbias the stored ratios: 
            ratio = ratio * pSingles

            if (ratio < max_frequency_bound) then

                if (.not. enough_sing_hist) then 
                    cnt_sing_hist = cnt_sing_hist + 1
                    if (cnt_sing_hist > cnt_threshold) enough_sing_hist = .true.
                end if

                ind = int(ratio / frq_step_size) + 1

                frequency_bins_singles(ind) = frequency_bins_singles(ind) + 1

            else
                above_max_singles = above_max_singles + 1
                print *, "Warning: single excitation H_ij/pgen above max_frequency_bound!" 
                print *, " H_ij/pgen: ", ratio, " ; bound: ", max_frequency_bound
                print *, " Consider increasing the bound!"
            end if
        else 

            ratio = ratio * pDoubles

            if (ratio < max_frequency_bound) then

                if (.not. enough_doub_hist) then 
                    cnt_doub_hist = cnt_doub_hist + 1
                    if (cnt_doub_hist > cnt_threshold) enough_doub_hist = .true.
                end if

                ind = int(ratio / frq_step_size) + 1

                frequency_bins_doubles(ind) = frequency_bins_doubles(ind) + 1

            else 
                above_max_doubles = above_max_doubles + 1
                print *, "Warning: double excitation H_ij/pgen above max_frequency_bound!" 
                print *, " H_ij/pgen: ", ratio, " ; bound: ", max_frequency_bound
                print *, " Consider increasing the bound!"
            end if
        end if

    end subroutine fill_frequency_histogram_sd


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
        integer :: ind
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

        else
            above_max_doubles = above_max_doubles + 1
            print *, "Warning: excitation H_ij/pgen above max_frequency_bound!" 
            print *, " H_ij/pgen: ", ratio, " ; bound: ", max_frequency_bound
            print *, " Consider increasing the bound!"
        end if

    end subroutine fill_frequency_histogram


    subroutine integrate_frequency_histogram_spec(spec_frequency_bins, ratio)
        ! specific histogram integration routine which sums up the inputted 
        ! frequency_bins
        integer, intent(in) :: spec_frequency_bins(n_frequency_bins) 
        real(dp), intent(out) :: ratio
        character(*), parameter :: this_routine = "integrate_frequency_histogram_spec"

        integer :: all_frequency_bins(n_frequency_bins)
        integer :: i, threshold
        integer :: n_elements, cnt

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

    ! also provide the printing routines here: 
    subroutine print_frequency_histograms
        ! try to write one general one and not multiple as in my GUGA branch
        use constants, only: int64
        use util_mod, only: get_free_unit, get_unique_filename
        use ParallelHelper, only: root
        
        character(*), parameter :: this_routine = "print_frequency_histograms"
        
        character(255) :: filename 
        integer :: all_frequency_bins_spec(n_frequency_bins)
        integer :: all_frequency_bins(n_frequency_bins)
        integer :: iunit, i, max_size
        real(dp) :: step_size, norm
        integer(int64) :: sum_all

        all_frequency_bins = 0

        ! maybe first check if we have only singles or only doubles like in 
        ! the real-space or momentum space hubbard: 
        if (tHub .or. tUEG) then 
            ! we only need to print one frequency_histogram:

            ! i need to change the rest of the code to use this 
            ! frequency_bins histogram in this case! TODO
            call MPIAllReduce(frequency_bins, MPI_SUM, all_frequency_bins)

            if (iProcIndex == root) then
                iunit = get_free_unit()
                call get_unique_filename('frequency_histogram', .true., &
                    .true., 1, filename)
                open(iunit, file = filename, status = 'unknown')

                do i = 1, n_frequency_bins
                    write(iunit, "(f16.7)", advance = 'no') frq_step_size * i
                    write(iunit, "(i12)") all_frequency_bins(i)
                end do

                close(iunit)
            end if

        else 
            ! in the other cases we definetly have singles and more
            ! first the singles:
            all_frequency_bins_spec = 0
            call MPIAllReduce(frequency_bins_singles, MPI_SUM, all_frequency_bins_spec) 

            if (iProcIndex == root) then
                iunit = get_free_unit()
                call get_unique_filename('frequency_histogram_singles', .true., &
                    .true., 1, filename)
                open(iunit, file = filename, status = 'unknown')

                do i = 1, n_frequency_bins
                    write(iunit, "(f16.7)", advance = 'no') frq_step_size * i
                    write(iunit, "(i12)") all_frequency_bins_spec(i)
                end do

                close(iunit)

                ! and add them up for the final normed one
                all_frequency_bins = all_frequency_bins_spec 

            end if

            all_frequency_bins_spec = 0

            ! do the cases where there is antiparallel or parallel
            if (tGen_4ind_weighted .or. tGen_4ind_2 .or. tGen_4ind_reverse) then 

                ! then the parallel:
                call MPIAllReduce(frequency_bins_para, MPI_SUM, all_frequency_bins_spec)

                if (iProcIndex == root) then
                    iunit = get_free_unit()
                    call get_unique_filename('frequency_histogram_para', .true., &
                        .true., 1, filename)
                    open(iunit, file = filename, status = 'unknown')

                    do i = 1, n_frequency_bins
                        write(iunit, "(f16.7)", advance = 'no') frq_step_size * i
                        write(iunit, "(i12)") all_frequency_bins_spec(i)
                    end do

                    close(iunit)

                    ! and add them up for the final normed one
                    all_frequency_bins = all_frequency_bins + all_frequency_bins_spec 

                end if

                all_frequency_bins_spec = 0
                ! then anti: 
                call MPIAllReduce(frequency_bins_anti, MPI_SUM, all_frequency_bins_spec)

                if (iProcIndex == root) then
                    iunit = get_free_unit()
                    call get_unique_filename('frequency_histogram_anti', .true., &
                        .true., 1, filename)
                    open(iunit, file = filename, status = 'unknown')

                    do i = 1, n_frequency_bins
                        write(iunit, "(f16.7)", advance = 'no') frq_step_size * i
                        write(iunit, "(i12)") all_frequency_bins_spec(i)
                    end do

                    close(iunit)

                    ! and add them up for the final normed one
                    all_frequency_bins = all_frequency_bins + all_frequency_bins_spec 

                end if

                all_frequency_bins_spec = 0

            else 
                ! we only have additional doubles:
                call MPIAllReduce(frequency_bins_doubles, MPI_SUM, all_frequency_bins_spec)

                if (iProcIndex == root) then
                    iunit = get_free_unit()
                    call get_unique_filename('frequency_histogram_doubles', .true., &
                        .true., 1, filename)
                    open(iunit, file = filename, status = 'unknown')

                    do i = 1, n_frequency_bins
                        write(iunit, "(f16.7)", advance = 'no') frq_step_size * i
                        write(iunit, "(i12)") all_frequency_bins_spec(i)
                    end do

                    close(iunit)

                    ! and add them up for the final normed one
                    all_frequency_bins = all_frequency_bins + all_frequency_bins_spec 

                end if
            end if
        end if

        ! and the norm then is the same for all cases:

        ! also print out a normed version for comparison
        sum_all = sum(all_frequency_bins)

        ! check if integer overflow:
        if (.not. sum_all < 0) then
            norm = real(sum_all, dp)

            iunit = get_free_unit()
            call get_unique_filename('frequency_histogram_normed', .true., &
                .true., 1, filename)
            open(iunit, file = filename, status = 'unknown')

            do i = 1, n_frequency_bins
                write(iunit, "(f16.7)", advance = 'no') frq_step_size * i
                write(iunit, "(f16.7)") real(all_frequency_bins(i),dp) / norm
            end do

            close(iunit)
        else
            write(iout,*) "Integer overflow in normed frequency histogram!" 
            write(iout,*) "DO NOT PRINT IT!"
        end if

    end subroutine print_frequency_histograms

    subroutine deallocate_histograms
        ! TODO: also do mem-logging here and in the allocation!
        character(*), parameter :: this_routine = "deallocate_histograms"

        call LogMemDealloc(this_routine, mem_tag_histograms)
        if (allocated(frequency_bins)) deallocate(frequency_bins)
        if (allocated(frequency_bins_singles)) deallocate(frequency_bins_singles)
        if (allocated(frequency_bins_doubles)) deallocate(frequency_bins_doubles)
        if (allocated(frequency_bins_para)) deallocate(frequency_bins_para)
        if (allocated(frequency_bins_anti)) deallocate(frequency_bins_anti)

    end subroutine deallocate_histograms

    subroutine read_frequency_histograms
        ! also write a routine, which in case of a restart reads in 
        ! saved frequency histograms and calculates time-step and 
        ! psingles etc. from it! although i am not sure if this is 
        ! possible without prior knowledge of the psingles etc..
        character(*), parameter :: this_routine = "read_frequency_histograms"

    end subroutine read_frequency_histograms

end module
