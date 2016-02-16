#include "macros.h"

module real_time

    use real_time_init, only: init_real_time_calc_single
    use real_time_procs, only: save_current_dets, reset_spawned_list, &
                               reload_current_dets
    use CalcData, only: pops_norm
    use real_time_data, only: gf_type
    use FciMCData, only: pops_pert, walker_time, iter, ValidSpawnedList, &
                         spawn_ht, FreeSlot, iStartFreeSlot, iEndFreeSlot, &
                         fcimc_iter_data, InitialSpawnedSlots, iter_data_fciqmc, &
                         TotWalkers
    use kp_fciqmc_data_mod, only: overlap_pert
    use timing_neci, only: set_timer, halt_timer
    use FciMCParMod, only: rezero_iter_stats_each_iter
    use hash, only: clear_hash_table
    use constants, only: int64, sizeof_int
    use AnnihilationMod, only: DirectAnnihilation

    implicit none

! main module file for the real-time implementation of the FCIQMC algorithm
! created on 04.02.2016 by Werner Dobrautz

! first i have to do a brainstorm on how to implement Ole's idea of the 
! real time version and also discuss with Ali, Simon and George to avoid 
! unnecessary effort

! implementation idea: 
! sample the real time Schroedinger equation by integrating:
! i d/dt |Psi(t)> = (H - E0 - ie)|Psi(t)> 

! with the 2nd order runge-kutta method: 
! d/dt y(t) = f(t,y)

! -> y(n+1) = y(n) + h f(n dt, y(n))
! t = n dt
! y(n+1) = y(n) + k2
! k1 = dt f(n dt, y(n)) 
! k2 = dt f(n + dt/2, y(n) + k1/2) 

! what does that mean in the dynamics? 

! from brainstorming on the 2nd order runge kutta implementation of the 
! original imaginary time FCIQMC algorithm, a few clarifications surfaced: 

! the RK2 needs the application of the square of the Hamiltonian essentially
! this implies huge changes to the underlying machinery in the current 
! FCIQMC implementation, which often relies on the type of exciations
! (single, doubles) 
! since the intermediatly created list y(n) + k1 / 2 already contains 
! single and double excitations, spawning from this list, increases the 
! possible excitations in the final y(n) + k2 -> y(n+1) to 
! singles, double, triples and quadrupels

! but the essential goal in the real-time code will be to obtain the 
! overlap: 
! <y(0)|y(t)> 

! to obtain the spectral information of the system up to a certain time 
! t_max

! i should devise an action plan to optimally implement this method 

! first of all i definetly need to use kneci to be able to handle the 
! imaginary and real parts of the equation
! dy(t)/dt = -i(H - E0 -ie)y(t) 

! then, since Ole showed already, and what i already know, back from the 
! work in Graz, a averaging over multiple runs, starting from different 
! stochastic representations of the converged ground state wave function is 
! necessary! 

! so probably a multiple kneci version, as developed by Dongxia is needed

! so whats the ideas?

! we have to start from a converged imaginary-time ground state wave-function! 
! a way could be to use a printed result from a previous calculationm through
! a popsfile. 
! but since we want average over multiple calculations, this would need 
! unnecessary amount of storage, I/O etc.
! but we should definetly be able to have this option! in the case the 
! imaginary-time calculation is very involved and already takes a lot of time 
! so maybe be able to read in #n popsfiles, printed at different times of 
! the imag-time calculation and start that amount of mkneci processes, which 
! do the real-time calculation then. 

! or run a "normal" imag-time run and then after convergence, start the 
! #n real-time calculations from diffferent starting points, as indicated by 
! Ole in this unterlagen
! or start with running multiple #n mneci runs with different seeds for the 
! random number generator, and then at the same time in equilibration, start 
! the real-time calculation

! this would imply a completly seperated code-basis in which the real-time 
! implementation is performed. which probably is a good idea anyway to keep 
! the code clean. since we need totally different information and stats 
! anyway..

! at each time step, we could then average <y(0)|y(t)> from the different 
! mkneci processes and just store one quantity and additional statistical info

! back in Graz i also came to the conclusion, that calculating the overlap 
! between different runs was of great help -> we could do that here too

! on the implementation of the new dynamics: 
! the underlying differential equation changes to (for the greater GF)
! dy(t)/dt = -i(H - E0 - ie)y(t) 
! this means 
! 1) we need both real and imaginary walkers 
! 2) when using the 2nd order RK method: 
!       y(n+1) = y(n) + k2(n)
!       k1(n) = -i dt (H - E0 - ie)y(n)
!       k2(n) = -i dt (H - E0 - ie)(y(n) + k1 / 2) 
!       
!       we could implement that by creating an intermediate determinant 
!       list obtained by the spawn(+annihilation) and death/cloning step
!       from y(n) + k1 / 2 
!       and based upon the spawns from that list we can combine that to the 
!       new y(n+1) = y(n) + k2 for the list at the next time-step
! 
!   or we can work out the final recursion formula:
!       y(n+1) = [(1 - e dt)(1 - i dt H') - dt^2 / 2 H'^2] y(n) (see doc.)
!
!   and do the spawning/cloning/death etc. directly based on this equation
!   lets see and talk to ali/ole about that 

! so the workflow of real-time propagations will be: 

! 1) run a normal imaginary time propagation until convergence and print out 
!       #n popsfiles as the groundstate information

! 2) terminate the imaginary run and start the real-time propagation with #n
!       threads mneci style

! 3) specify a set of specific orbitals j, on which either a^+_j or a_j acts
!       with possible symmetry constraints also, or specify k-space calculation
!       like for the Hubbard or UEG model. 
!       for each j or k also the coresponding <0|a_i/a^+_i or multiple if 
!       possible act, to calculate the overlaps. applying multiple <0|a_i is 
!       cheap, so all of them should be done. (for k-space k'=k due to symmetry)

! 4) do the actual real time propagation! 

contains

    subroutine perform_real_time_fciqmc
        ! main real-time calculation routine
        ! do all the setup, read-in and calling of the "new" real-time MC loop
        use real_time_procs, only: update_gf_overlap
        use real_time_data, only: gf_overlap
        implicit none

        character(*), parameter :: this_routine = "perform_real_time_fciqmc"
        integer :: n_determ_states

        print *, " ========================================================== "
        print *, " ------------------ Real-time FCIQMC ---------------------- "
        print *, " ========================================================== "

        ! call the real-time setup routine and all the initialization
        call init_real_time_calc_single()

        print *, " Real-time FCIQMC initialized! "
        ! rewrite the major original neci core loop here and adapt it to 
        ! the new necessary real-time stuff
        ! check nicks kp code, to have a guideline in how to go into that! 

        ! as a first test if everything is set up correctly calculate the 
        ! initial overlap for the same indiced of creation and annihilation
        ! operators <y(0)|a^+_i a_i |y(0)> = n_i 
        ! and       <y(0)|a_i a^+_i |y(0)> = 1 - n_i 
        ! if normed correctly
        ! in the quantity perturbed_ground the left hand side and in 
        ! CurrentDets the right hand side should be stored ..
        ! should write a routine which calulated the overlap of CurrentDets
        ! with the perturbed groundstate and stores it in a pre-allocated list
        ! shouldnt geourge have some of those routines..

        call update_gf_overlap()

        print *, "test on overlap at t = 0: "
        if (gf_type == -1) then
            print *, " for lesser GF  <y(0)| a^+_i a_j |y(0) >; i,j: ", &
                overlap_pert(1)%ann_orbs(1), pops_pert(1)%ann_orbs(1)
        else if (gf_type == 1) then
            print *, " for greater GF <y(0)| a_i a^+_j |y(0)> ; i,j: ", &
                overlap_pert(1)%crtn_orbs(1), pops_pert(1)%crtn_orbs(1)
        end if

        print *, " overlap: ",  sqrt(gf_overlap(0)) / sqrt(pops_norm) 

        ! main part of the initialization, concerning the walker lists 
        ! and operator application works now.. 
        ! so next step is to finish up the whole initialization to be able 
        ! to run a (mk)neci with the new walker dynamics provided by the 
        ! real-time Schroedinger equation

        ! enter the main real-time fciqmc loop here
        fciqmc_loop: do while (.true.)

            ! the timing stuff has to be done a bit differently in the 
            ! real-time fciqmc, since there are 2 spawing and annihilation 
            ! steps involved...
            call set_timer(walker_time)

            iter = iter + 1

            
            ! perform the actual iteration(excitation generation etc.) 
            call perform_real_time_iteration() 

            ! update various variables.. time-step, initiator criteria etc.
            call update_real_time_iteration()

            ! update, print and log all the global variables and interesting 
            ! quantities
            call log_real_time_iteration() 

            ! check if somthing happpened to stop the iteration or something
            call check_real_time_iteration()

            call stop_all(this_routine, "stop for now to avoid endless loop")

        end do fciqmc_loop
        

    end subroutine perform_real_time_fciqmc

    subroutine update_real_time_iteration()
        ! routine to update certain global variables each loop iteration in 
        ! the real-time fciqmc 
        character(*), parameter :: this_routine = "update_real_time_iteration"

    end subroutine update_real_time_iteration

    subroutine log_real_time_iteration()
        ! routine to log all the interesting quantities in the real-time 
        ! fciqmc 
        character(*), parameter :: this_routine = "log_real_time_iteration"

    end subroutine log_real_time_iteration

    subroutine check_real_time_iteration()
        ! routine to check if somthing wrong happened during the main 
        ! real-time fciqmc loop or the external CHANGEVARS utility does smth
        character(*), parameter :: this_routine = "check_real_time_iteration"

    end subroutine check_real_time_iteration

    subroutine init_real_time_iteration(iter_data, n_determ_states)
        ! routine to reinitialize all the necessary variables and pointers 
        ! for a sucessful real-time fciqmc iteration
        type(fcimc_iter_data), intent(inout) :: iter_data
        integer, intent(out) :: n_determ_states
        character(*), parameter :: this_routine = "init_real_time_iteration"

        ! reuse parts of nicks routine and add additional real-time fciqmc
        ! specific quantities
        
        ! Reset positions to spawn into in the spawning array.
        ValidSpawnedList = InitialSpawnedSlots

        ! Reset the array which holds empty slots in CurrentDets.
        ! hm, where is the iEndFreeSlot var. set? 
        ! well i guess it uses the value from the last iteration and then 
        ! gets reset..
        FreeSlot(1:iEndFreeSlot) = 0
        iStartFreeSlot = 1
        iEndFreeSlot = 0


        ! Index for counting deterministic states.
        n_determ_states = 1

        ! Clear the hash table for the spawning array.
        call clear_hash_table(spawn_ht)

        call rezero_iter_stats_each_iter(iter_data)

        ! additionaly i have to copy the CurrentDets array and all the 
        ! associated pointers and hashtable related stuff to the 
        ! temporary 2nd list 
        call save_current_dets() 

    end subroutine init_real_time_iteration
    
    subroutine perform_real_time_iteration()
        ! routine which performs one real-time fciqmc iteration
        character(*), parameter :: this_routine = "perform_real_time_iteration"

        integer :: n_determ_states, idet 
        integer(int64) :: temp_totWalkers, TotWalkersNew

        ! 0)
        ! do all the necessary preperation(resetting pointers etc.)
        call init_real_time_iteration(iter_data_fciqmc, n_determ_states)
        ! 1)
        ! do a "normal" spawning step and combination to y(n) + k1/2
        ! into CurrentDets: 
        !   
        do idet = 1, int(TotWalkers, sizeof_int) 
            ! todo: includes walker_death too..
            ! so think on how to exclude that
        end do

        TotWalkersNew = int(TotWalkers, sizeof_int)

        ! Annihilation is done after loop over walkers
        call DirectAnnihilation (TotWalkersNew, iter_data_fciqmc, .false.)

        TotWalkers = int(TotWalkersNew, sizeof_int)

        ! 2)
        ! reset the spawned list and do a second spawning step to create 
        ! the spawend list k2 
        ! but DO NOT yet recombine with stored walker list 
        call reset_spawned_list(n_determ_states) 

        ! create a second spawned list from y(n) + k1/2
        do idet = 1, int(TotWalkers, sizeof_int)

            ! have to think to exclude death_step here and store this 
            ! information into the spawned k2 list..
            ! quick solution would be to loop again over reloaded y(n)
            ! and do a death step for wach walker

        end do

        ! 3) 
        ! reload stored temp_det_list y(n) into CurrentDets 
        call reload_current_dets()

        ! 4)
        ! for the death_step for now: loop once again over the walker list 
        ! and do a death_step for each walker..
        ! meh.. that seems really inefficient, better do it more cleverly 
        ! in the creation of the k2 spawned list + annihilation!

!         do idet = 1, int(TotWalkers, sizeof_int)
! 
!             ! do death_step only.. maybe..
!         end do

        ! combine y(n+1) = y(n) + k2 into CurrentDets to finish time-step
        ! this should be done with a single Annihilation step between
        ! y(n) = CurrentDets and k2 = SpawnedWalkers

        TotWalkersNew = int(TotWalkers, sizeof_int)

        ! Annihilation is done after loop over walkers
        call DirectAnnihilation (TotWalkersNew, iter_data_fciqmc, .false.)

        TotWalkers = int(TotWalkersNew, sizeof_int)

    end subroutine perform_real_time_iteration

end module real_time

! wrapper (dont know why this is necessary quite..)
subroutine perform_real_time_fciqcm_wrap
    use real_time, only: perform_real_time_fciqmc
    implicit none

    call perform_real_time_fciqmc

end subroutine perform_real_time_fciqcm_wrap
