#include "macros.h"

module real_time

    use real_time_init, only: init_real_time_calc_single, dealloc_real_time_memory, &
         rotate_time
    use real_time_procs, only: save_current_dets, reset_spawned_list, merge_spawn, &
                               reload_current_dets, walker_death_realtime, &
                               walker_death_spawn, attempt_die_realtime, trunc_shift, &
                               create_diagonal_as_spawn, count_holes_in_currentDets, &
                               DirectAnnihilation_diag, check_update_growth, &
                               update_gf_overlap, calc_norm, adjust_decay_channels, &
                               update_shift_damping, real_time_determ_projection, &
                               refresh_semistochastic_space, update_peak_walker_number, &
                               makePopSnapshot, update_elapsed_time, logTimeCurve, &
                               get_current_alpha_from_cache, closeTauContourFile, &
                               get_corespace_from_buf
    use real_time_data, only: gf_type,  tVerletSweep, &
                              pert_norm, second_spawn_iter_data, runge_kutta_step,&
                              current_overlap, SumWalkersCyc_1, DiagParts, stepsAlpha, &
                              elapsedRealTime, elapsedImagTime, TotPartsPeak, tVerletScheme, &
                              tau_real, tau_imag, t_rotated_time, temp_iendfreeslot, &
                              temp_freeslot, overlap_real, overlap_imag, dyn_norm_psi, &
                              NoatHF_1, shift_damping, tDynamicCoreSpace, dyn_norm_red, &
                              normsize, gf_count, tRealTimePopsfile, tStabilizerShift, &
                              tLimitShift, tDynamicAlpha, dpsi_cache, dpsi_size, tGZero, &
                              tDynamicDamping, stabilizerThresh, popSnapshot, spawnBufSize, &
                              tLogTrajectory, tReadTrajectory, tGenerateCoreSpace, &
                              numSnapShotOrbs, core_space_buf, csbuf_size, corespace_log_interval, &
                              real_time_info
    use verlet_aux, only: init_verlet_iteration, obtain_h2_psi, update_delta_psi, &
         init_verlet_sweep, check_verlet_sweep, end_verlet_sweep
    use CalcData, only: pops_norm, tTruncInitiator, tPairedReplicas, ss_space_in, &
                        tDetermHFSpawning, AvMCExcits, tSemiStochastic, StepsSft, &
                        tChangeProjEDet, DiagSft, nmcyc, tau, InitWalkers, &
                        s_global_start, StepsSft, semistoch_shift_iter
    use FciMCData, only: pops_pert, walker_time, iter, ValidSpawnedList, spawnedParts, &
                         spawn_ht, FreeSlot, iStartFreeSlot, iEndFreeSlot,  & 
                         fcimc_iter_data, InitialSpawnedSlots, iter_data_fciqmc, &
                         TotWalkers, fcimc_excit_gen_store, ilutRef, max_calc_ex_level, &
                         iLutHF_true, indices_of_determ_states, partial_determ_vecs, &
                         exFlag, CurrentDets, TotParts, ilutHF, SumWalkersCyc, IterTime, &
                         HFCyc, norm_psi, NoatHF, NoDied, tTimeExit, maxTimeExit, &
                         Annihilated, NoBorn, tSinglePartPhase, AllSumNoatHF, AllTotParts
    use kp_fciqmc_data_mod, only: overlap_pert
    use timing_neci, only: set_timer, halt_timer
    use FciMCParMod, only: rezero_iter_stats_each_iter
    use hash, only: clear_hash_table
    use constants, only: int64, sizeof_int, n_int, lenof_sign, dp, EPS, inum_runs, bits_n_int, &
         iout
    use AnnihilationMod, only: DirectAnnihilation, AnnihilateSpawnedParts, &
         deterministic_annihilation
    use bit_reps, only: extract_bit_rep, decode_bit_det
    use SystemData, only: nel, tRef_Not_HF, tAllSymSectors, nOccAlpha, nOccBeta, &
                          nbasis
    use DetBitOps, only: FindBitExcitLevel, return_ms
    use semi_stoch_procs, only: check_determ_flag, determ_projection, write_core_space
    use semi_stoch_gen, only: init_semi_stochastic, write_most_pop_core_at_end
    use global_det_data, only: det_diagH
    use fcimc_helper, only: CalcParentFlag, decide_num_to_spawn, &
                            create_particle_with_hash_table, walker_death, &
                            SumEContrib, end_iter_stats, check_semistoch_flags
    use procedure_pointers, only: generate_excitation, encode_child, &
                                  attempt_create, new_child_stats
    use bit_rep_data, only: tUseFlags, nOffFlag, niftot, extract_sign
    use bit_reps, only: set_flag, flag_deterministic, flag_determ_parent, test_flag
    use fcimc_iter_utils, only: update_iter_data, collate_iter_data, iter_diagnostics, &
                                population_check, update_shift, calculate_new_shift_wrapper
    use soft_exit, only: ChangeVars, tSoftExitFound
    use fcimc_initialisation, only: CalcApproxpDoubles
    use LoggingData, only: tPopsFile, write_end_core_size, tWriteCoreEnd
    use PopsFileMod, only: WriteToPopsfileParOneArr
    use load_balance, only: test_hash_table, tLoadBalanceBlocks, adjust_load_balance, &
         CalcHashTableStats
    use Parallel_neci
    use util_mod, only : neci_etime
    use ParallelHelper, only: MPI_SUM, iProcIndex, root
    use adi_references, only: setup_reference_space
    use adi_data, only: allDoubsInitsDelay, nRefs, tDelayGetRefs

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

  subroutine check_walker_number(iter_data, message)
    implicit none
    type(fcimc_iter_data), intent(in) :: iter_data
    character(len=*), intent(in) :: message
    real(dp) :: growth(lenof_sign)
    
    growth = iter_data%nborn &
         - iter_data%ndied - iter_data%nannihil &
         - iter_data%naborted - iter_data%nremoved
    
    print *, message, "Iter data growth:", growth
    
  end subroutine check_walker_number


    subroutine perform_real_time_fciqmc
        ! main real-time calculation routine
        ! do all the setup, read-in and calling of the "new" real-time MC loop
        use real_time_data, only: gf_overlap
        use FciMCData, only : TotImagTime
        use real_time_procs, only: expand_corespace_buf
        use fcimc_output, only: PrintHighPops
        implicit none

        character(*), parameter :: this_routine = "perform_real_time_fciqmc"
        integer :: j, i, iterRK
        real(dp) :: s_start, s_end, tstart(2), tend(2)
        real(dp) :: totalTime
        complex(dp), allocatable :: overlap_buf(:)
        complex(dp), allocatable :: norm_buf(:)
        character (255) :: rtPOPSFILE_name

        rtPOPSFILE_name = 'TIME_EVOLVED_POP'

        write(iout,*) " ========================================================== "
        write(iout,*) " ------------------ Real-time FCIQMC ---------------------- "
        write(iout,*) " ========================================================== "

        ! call the real-time setup routine and all the initialization
       
        call init_real_time_calc_single()

        ! counts the number of iterations in Runge-Kutta when creating 
        ! initial states for verlet
        iterRK = 0
        
        write(iout,*)  " Real-time FCIQMC initialized! "
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

        ! normsize and gf_count are initialized in init_real_time_calc_single
        ! -> cannot be used before
        allocate(overlap_buf(gf_count),stat = i)
        allocate(norm_buf(normsize), stat = i)

        call update_gf_overlap()
        if(iProcIndex == root) then
           print *, "test on overlap at t = 0: "
           if(.not. tRealTimePopsfile) then
              if (gf_type == -1) then
                 print *, " for lesser GF  <y(0)| a^+_i a_j |y(0) >; i,j: ", &
                      overlap_pert(1)%ann_orbs(1), pops_pert(1)%ann_orbs(1)
              else if (gf_type == 1) then
                 print *, " for greater GF <y(0)| a_i a^+_j |y(0)> ; i,j: ", &
                      overlap_pert(1)%crtn_orbs(1), pops_pert(1)%crtn_orbs(1)
              end if
           endif
              print *, "Current GF:", gf_overlap(1,1)/pert_norm(1,1), pert_norm(1,1), normsize
              print *, "Normalization", pert_norm(1,1), dyn_norm_red(1,1)
        endif

        ! enter the main real-time fciqmc loop here
        fciqmc_loop: do while (.true.)

            ! the timing stuff has to be done a bit differently in the 
            ! real-time fciqmc, since there are 2 spawing and annihilation 
            ! steps involved...
            
            call set_timer(walker_time)

            if(iProcIndex == root) s_start = neci_etime(tstart)

            ! the iter data is used in updating the shift, so it has to be reset before 
            ! output
            ! this is a bad implementation : iter should be local

            ! if the trajectory is logged, print alpha and tau here
            if(tLogTrajectory) call logTimeCurve()

            ! update the overlap each time
            ! rmneci_setup: computation of instantaneous projected norm is shifted to here
            if(mod(iter, StepsSft) == 0) then 
               ! current overlap is now the one after iteration
               ! update the normalization due to the shift

               norm_buf = calc_norm(CurrentDets,int(TotWalkers))
               call MPIReduce(norm_buf,MPI_SUM,dyn_norm_psi)

               call update_gf_overlap()

               do j = 1, gf_count
                  do i = 1, normsize
                     current_overlap(i,j) = gf_overlap(i,j)/dyn_norm_red(i,j) * &
                          exp(shift_damping(((i-1)/inum_runs+1)))
                  end do

                  !normalize the greens function
                  ! averaging should probably be done over the normalized
                  ! GFs
                  overlap_buf(j) = sum(gf_overlap(:,j))/sum(dyn_norm_red(:,j)) * &
                       sum(exp(shift_damping))/inum_runs

                  overlap_real(j) = real(overlap_buf(j))
                  overlap_imag(j) = aimag(overlap_buf(j))
               end do
               call update_real_time_iteration()
            endif

            if(mod(iter,stepsAlpha)==0 .and. (tDynamicAlpha .or. tDynamicDamping)) then
               call adjust_decay_channels()
               ! the verlet scheme only works for constant time step, hence, upon adjusting
               ! the time step, we need to do another iteration of RK2
               if(tVerletScheme) call end_verlet_sweep()
            endif

            if(tVerletScheme) call check_verlet_sweep(iterRK)

            ! if a threshold value is set, check it
            if(tLimitShift) call trunc_shift()

            if(tGenerateCoreSpace .and. (mod(iter,corespace_log_interval) .eq. 0)) then
               call expand_corespace_buf(core_space_buf, csbuf_size)
               write(6,*) "New corespace buffer size", csbuf_size
            endif
            
            ! perform the actual iteration(excitation generation etc.) 
            if(iterRK .eq. 0) iter = iter + 1
            if(tVerletScheme .and. .not. tVerletSweep) iterRK = iterRK + 1            
            if(tVerletSweep) then
               call perform_verlet_iteration()
            else
               call perform_real_time_iteration() 
            endif

            ! check if somthing happpened to stop the iteration or something
            call check_real_time_iteration()

            ! If required, set up the references
            if(Iter == allDoubsInitsDelay + 1 .and. nRefs > 1 .and. tDelayGetRefs) &
                 call setup_reference_space(.true.)

            if(iProcIndex.eq.root) then 
               s_end = neci_etime(tend)
               totalTime = real(s_end - s_global_start, dp)
            endif
            call MPIBcast(totalTime)

            if(tTimeExit .and. totalTime > MaxTimeExit) then
               ! reset the maximum number of iterations to the index of the next one
               nmcyc = Iter + StepsSft
               ! to prevent nmcyc from being updated
               tTimeExit = .false.
            endif

            ! load balancing
            if(tLoadBalanceBlocks .and. mod(iter,1000) == 0 .and. (.not. tSemiStochastic)) then
               call adjust_load_balance(iter_data_fciqmc)
            endif

            if(mod(iter,400) == 0 .and. tSemiStochastic .and. tDynamicCoreSpace) &
                 call refresh_semistochastic_space()
            if(iProcIndex == root) then
               s_end = neci_etime(tend)
               IterTime = IterTime + (s_end - s_start)
            endif

            ! update the normalization of the greensfunction according to damping (dynamic)
            call update_shift_damping()

            if (tSoftExitFound) exit fciqmc_loop
            ! we have to stop here since the gf_overlap array is now full
            if (iter == nmcyc) exit fciqmc_loop

        end do fciqmc_loop

        call PrintHighPops()
        
        if (tPopsFile) then 
           ! as both the elapsed real-time and the elapsed imaginary time are required
           ! for dynamic alpha, both are stored. That is, the total elapsed real time
           ! is now stored instead of tau. If no tau is supplied upon read-in of the
           ! generated popsfile, it is estimated using the number of cycles, the current
           ! angle of rotation and the total elapsed real time
           TotImagTime = elapsedImagTime
           tau = elapsedRealTime
           ! THIS IS A HACK: We dont want to alter the POPSFILE functions themselves
           ! so we sneak in the shift_damping into some slot unimportant to rneci
           AllSumNoatHF(1:inum_runs) = shift_damping
           ! Another hack: we also sneak the value of alpha somewhere here, the shift
           ! is not meaningful in real-time anyway
           DiagSft = real_time_info%time_angle
           call WriteToPopsfileParOneArr(CurrentDets,TotWalkers,rtPOPSFILE_name)
        endif

        ! We finish writing any extra output files like the tauContour and
        ! the CORESPACE file
        if(tLogTrajectory) call closeTauContourFile()

        if(tWriteCoreEnd) call write_most_pop_core_at_end(write_end_core_size)
        
        ! GENERATE-CORESPACE has precendence over WRITE-CORE-END
        if(tGenerateCoreSpace) call get_corespace_from_buf(core_space_buf, csbuf_size)
        
        deallocate(norm_buf, stat = i)
        deallocate(overlap_buf, stat = i)

        ! rt_iter_adapt : avoid memleaks
        call dealloc_real_time_memory()

    end subroutine perform_real_time_fciqmc

    subroutine update_real_time_iteration()
        ! routine to update certain global variables each loop iteration in 
        ! the real-time fciqmc 
        ! from the 2 distince spawn/death/cloning info stored in the 
        ! two iter_data vars, i have to combine the general updated 
        ! statistics for the actual time step
      implicit none
        character(*), parameter :: this_routine = "update_real_time_iteration"
        integer :: run

        ! how to combine those 2? 
        ! in iter_data_fciqmc the info on the born, died, aborted, removed and 
        ! annihilated particles of the first spawn and y(n) + k1/2 step is stored

        ! in the second_spawn_iter_data, the number of born particles in the 
        ! spawing step is stored first.. 
        ! note: in the create_particle routine the global variable 
        ! acceptances gets updated, and i essentially update that 
        ! quantity twice when calling it in the first and second spawn 
        ! loop -> i think i have to store 2 of those variables and 
        ! update and reset both of them seperately to keep track of the 
        ! statistics correctly
        ! and the NoBorn, NoDied and similar variables also get used in the 
        ! statistics about the simulation.. maybe i need to adjust them too
        ! so take the 2 RK loops into account
        ! do i want to use the new_shift_wrapper here? .. 
        ! no i think i just want to combine the important infos from both the 
        ! iter_datas so to get the correct and valid info for the full 
        ! time-step ... hm.. 

        ! do a correct combination of the essential parts of the new_shift_wrapper
        ! in the end i could just use the calc_new_shift_wrapper.. 

        ! still have to do more combination of necessary data..

        ! combine log_real_time into this routine too! 
        ! get the norm of the state

        if(tReadTrajectory) call get_current_alpha_from_cache

        call calculate_new_shift_wrapper(second_spawn_iter_data, totParts, &
             tPairedReplicas)

        if(tStabilizerShift) then
           if(iProcIndex == Root) then
              ! check if the walker number started to decay uncontrolled
              call update_peak_walker_number()
              do run = 1,inum_runs
                 if((AllTotParts(min_part_type(run))+ AllTotParts(max_part_type(run)))&
                      < stabilizerThresh*TotPartsPeak(run) .and. tSinglePartPhase(run)) then
                    ! if it is, enable dynamic shift to enforce a sufficiently high walker number
                    tSinglePartPhase(run) = .false.
                    write(6,*) "Walker number dropped below threshold, enabling dynamic shift"
                 end if
              end do
           end if
           ! and do not forget to communicate the decision
           call MPIBcast(tSinglePartPhase)
        end if     

        call rotate_time()
    end subroutine update_real_time_iteration

    subroutine log_real_time_iteration()
        ! routine to log all the interesting quantities in the real-time 
        ! fciqmc 
        ! have to figure out what i want to have an output of.. and then 
        ! print that to a FCIMCstats file!
        character(*), parameter :: this_routine = "log_real_time_iteration"

    end subroutine log_real_time_iteration

    subroutine check_real_time_iteration()
        ! routine to check if somthing wrong happened during the main 
        ! real-time fciqmc loop or the external CHANGEVARS utility does smth
        character(*), parameter :: this_routine = "check_real_time_iteration"
        logical :: tSingBiasChange, tWritePopsFound

        if (mod(iter, StepsSft) == 0) then
            call ChangeVars(tSingBiasChange, tWritePopsFound)
            if (tWritePopsFound) call WriteToPopsfileParOneArr(CurrentDets, TotWalkers)
            if (tSingBiasChange) call CalcApproxpDoubles()
        end if

        if(semistoch_shift_iter/=0) then
           if(Iter == semistoch_shift_iter + 1) then
              tSemiStochastic = .true.
              call init_semi_stochastic(ss_space_in)
           endif
        endif

    end subroutine check_real_time_iteration

    subroutine init_real_time_iteration(iter_data, iter_data2)
        ! routine to reinitialize all the necessary variables and pointers 
        ! for a sucessful real-time fciqmc iteration
      ! RT_M_Merge: Added dummy argument for rezero_iter_stats_each_iter
      use rdm_data, only: rdm_definitions_t
        type(fcimc_iter_data), intent(inout) :: iter_data
        type(fcimc_iter_data), intent(inout), optional :: iter_data2
        type(rdm_definitions_t) :: dummy
!         integer, intent(out) :: n_determ_states
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
        temp_freeslot = 0
        temp_iendfreeslot = 0

        ! also reset the temporary variables
        ! this probably has not to be done, since at the end of the first 
        ! spawn i set it to the freeslot array anyway..
        ! with changed temp_ var. usage i have to reset them! 

        ! Index for counting deterministic states.
!         n_determ_states = 1

        ! reset population snapshot
        popSnapshot = 0

        ! Clear the hash table for the spawning array.
        call clear_hash_table(spawn_ht)

        call rezero_iter_stats_each_iter(iter_data,dummy)
        if (present(iter_data2)) then
            call rezero_iter_stats_each_iter(iter_data2, dummy)
        end if

        ! additionaly i have to copy the CurrentDets array and all the 
        ! associated pointers and hashtable related stuff to the 
        ! temporary 2nd list 
        call save_current_dets() 
        call update_elapsed_time()
    end subroutine init_real_time_iteration


    subroutine second_real_time_spawn()
        ! routine for the second spawning step in the 2nd order RK method
        ! to create the k2 spawning list. An important change: 
        ! this "spawning" list also has to contain the diagonal death/cloning
        ! step influence to combine it with the original y(n) 
!         integer, intent(inout) :: n_determ_states
        character(*), parameter :: this_routine = "second_real_time_spawn"

        ! mimic the most of this routine to the already written first
        ! spawning step, but with a different death step and without the 
        ! annihilation step at the end! 
        integer :: idet, parent_flags, nI_parent(nel), unused_flags, ex_level_to_ref, &
                   ireplica, nspawn, ispawn, nI_child(nel), ic, ex(2,2), &
                   ex_level_to_hf 
        integer(n_int) :: ilut_child(0:niftot)
        real(dp) :: parent_sign(lenof_sign), parent_hdiag, prob, child_sign(lenof_sign), &
                    unused_sign(lenof_sign), unused_rdm_real, diag_sign(lenof_sign)
        logical :: tParentIsDeterm, tParentUnoccupied, tParity, break
        HElement_t(dp) :: HelGen
        integer :: determ_index
        real(dp) :: prefactor

        ! declare this is the second runge kutta step
        runge_kutta_step = 2

        ! index for counting deterministic states
        determ_index = 1
        ! prefactor for unbiasing if the number of spawns is cut off
        prefactor = 1.0_dp
        ! use part of nicks code, and remove the parts, that dont matter 
        do idet = 1, int(TotWalkers, sizeof_int)

            parent_flags = 0_n_int

            ! Indicate that the scratch storage used for excitation generation from the
            ! same walker has not been filled (it is filled when we excite from the first
            ! particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.

            call extract_bit_rep(CurrentDets(:,idet), nI_parent, parent_sign, unused_flags, &
                                  fcimc_excit_gen_store)

            ex_level_to_ref = FindBitExcitLevel(iLutRef, CurrentDets(:,idet), max_calc_ex_level)
            if(tRef_Not_HF) then
                ex_level_to_hf = FindBitExcitLevel (iLutHF_true, CurrentDets(:,idet), &
                     max_calc_ex_level)
            else
                ex_level_to_hf = ex_level_to_ref
            endif

            tParentIsDeterm = check_determ_flag(CurrentDets(:,idet))
            tParentUnoccupied = IsUnoccDet(parent_sign)            
            
            if(tParentIsDeterm) then
               indices_of_determ_states(determ_index) = idet
               partial_determ_vecs(:,determ_index) = parent_sign
               determ_index = determ_index + 1
               if(IsUnoccDet(parent_sign)) cycle
            endif

            if (tTruncInitiator) call CalcParentFlag(idet, parent_flags)

            ! If this slot is unoccupied (and also not a core determinant) then add it to
            ! the list of free slots and cycle.
            ! actually dont need to update this here since nothing gets merged
            ! at the end.. only k2 gets created
            if (tParentUnoccupied) then
               ! this is unnecessary in the end, the merge is done 
               ! with the original ensemble, with the original FreeSlots
               ! iEndFreeSlot = iEndFreeSlot + 1
               ! FreeSlot(iEndFreeSlot) = idet
               cycle
            end if
            ! The current diagonal matrix element is stored persistently.
            parent_hdiag = det_diagH(idet)


            ! UPDATE: call this routine anyway to update info on noathf 
            ! and noatdoubs, for the intermediate step
            call SumEContrib (nI_parent, ex_level_to_ref, parent_sign, CurrentDets(:,idet), &
                              parent_hdiag, 1.0_dp, tPairedReplicas, idet)

            ! If we're on the Hartree-Fock, and all singles and
            ! doubles are in the core space, then there will be
            ! no stochastic spawning from this determinant, so
            ! we can the rest of this loop.
            if (ss_space_in%tDoubles .and. ex_level_to_hf == 0 .and. tDetermHFSpawning) cycle

            ! i dont thin i need this below..
!             if (tAllSymSectors) then
!                 ms_parent = return_ms(ilut_parent)
!                 nOccAlpha = (nel+ms_parent)/2
!                 nOccBeta = (nel-ms_parent)/2
!             end if
! 
            ! If this condition is not met (if all electrons have spin up or all have spin down)
            ! then there will be no determinants to spawn to, so don't attempt spawning.
            ! thats a really specific condition.. shouldnt be checked each 
            ! cycle.. since this is only input dependent..
            if(.not. tGZero) then ! skip this if we only want the corespace-evolution
               do ireplica = 1, lenof_sign

                  call decide_num_to_spawn(parent_sign(ireplica), parent_hdiag, AvMCExcits, nspawn)
                  !call merge_spawn(nspawn,prefactor)
                  do ispawn = 1, nspawn

                     ! Zero the bit representation, to ensure no extraneous data gets through.
                     ilut_child = 0_n_int
                     call generate_excitation (nI_parent, CurrentDets(:,idet), nI_child, &
                          ilut_child, exFlag, ic, ex, tParity, prob, &
                          HElGen, fcimc_excit_gen_store)
                     ! If a valid excitation.
                     if (.not. IsNullDet(nI_child)) then

                        call encode_child (CurrentDets(:,idet), ilut_child, ic, ex)
                        if (tUseFlags) ilut_child(nOffFlag) = 0_n_int

                        if (tSemiStochastic) then
                           break = check_semistoch_flags(ilut_child, nI_child, tParentIsDeterm)
                           if(break) cycle
                        endif

                        ! unbias if the number of spawns was truncated
                        child_sign = attempt_create (nI_parent, CurrentDets(:,idet), parent_sign, &
                             nI_child, ilut_child, prob, HElGen, ic, ex, tParity, &
                             ex_level_to_ref, ireplica, unused_sign, unused_rdm_real)
                        child_sign = prefactor*child_sign
                     else
                        child_sign = 0.0_dp
                     end if

                     ! If any (valid) children have been spawned.
                     if ((any(abs(child_sign) > EPS)) .and. (ic /= 0) .and. (ic <= 2)) then

                        ! not quite sure about how to collect the child
                        ! stats in the new rt-fciqmc ..
                        call new_child_stats (second_spawn_iter_data, CurrentDets(:,idet), &
                             nI_child, ilut_child, ic, ex_level_to_ref, &
                             child_sign, parent_flags, ireplica)
                        call create_particle_with_hash_table (nI_child, ilut_child, child_sign, &
                             ireplica, CurrentDets(:,idet), &
                             second_spawn_iter_data)
                     end if ! If a child was spawned.

                  end do ! Over mulitple particles on same determinant.

               end do ! Over the replicas on the same determinant.
            endif
            ! If this is a core-space determinant then the death step is done in
            ! determ_projection.
            if (.not. tParentIsDeterm) then
                ! essentially have to treat the diagonal clone/death step like 
                ! a spawning step and store the result into the spawned array..
!                 call walker_death_spawn ()!iter_data_fciqmc, nI_parent,  &
!                     ilut_parent, parent_hdiag, parent_sign, idet, ex_level_to_ref)

                ! also not quite sure how to do the child_stats here ...
                ! or in general for now in the rt-fciqmc
                ! have to write all new book-keeping routines i guess.. 
                ! i also have to change the sign convention here, since in 
                ! the annihilation(and in the spawing routine above) the 
                ! particles get merged with a + instead of the - in the 
                ! original death routine for a "normal" death
                diag_sign = -attempt_die_realtime(parent_hdiag, &
                    parent_sign, ex_level_to_ref)
                if (any(abs(diag_sign) > EPS)) then
                   call create_diagonal_as_spawn(CurrentDets(:,idet), &
                        diag_sign, second_spawn_iter_data)
                end if

            end if

        end do ! Over all determinants. 

        ! the deterministic time evoulution is performed here according to the 
        ! RK scheme (annihilation is done separately as deterministic_annihilation)

        if (tSemiStochastic) call real_time_determ_projection() 

    end subroutine second_real_time_spawn

    subroutine first_real_time_spawn()
        ! routine which first loops over the CurrentDets array and creates the 
        ! first spawning list k1 and combines it to y(n) + k1/2
!         integer, intent(inout) :: n_determ_states
      implicit none
        character(*), parameter :: this_routine = "first_real_time_spawn"
        integer :: idet, parent_flags, nI_parent(nel), unused_flags, ex_level_to_ref, &
                   ireplica, nspawn, ispawn, nI_child(nel), ic, ex(2,2), &
                   ex_level_to_hf
        integer(n_int) :: ilut_child(0:niftot)
        real(dp) :: parent_sign(lenof_sign), parent_hdiag, prob, child_sign(lenof_sign), &
                    unused_sign(lenof_sign), prefactor, unused_rdm_real
        logical :: tParentIsDeterm, tParentUnoccupied, tParity, break
        HElement_t(dp) :: HelGen
        integer :: TotWalkersNew, run, determ_index

        ! declare this is the first runge kutta step
        runge_kutta_step = 1
        ! use part of nicks code, and remove the parts, that dont matter 
        prefactor = 1.0_dp
        determ_index = 1

        do idet = 1, int(TotWalkers, sizeof_int)

            ! The 'parent' determinant from which spawning is to be attempted.
            parent_flags = 0_n_int

            ! Indicate that the scratch storage used for excitation generation from the
            ! same walker has not been filled (it is filled when we excite from the first
            ! particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.

            call extract_bit_rep(CurrentDets(:,idet), nI_parent, parent_sign, unused_flags, &
                                  fcimc_excit_gen_store)

            ex_level_to_ref = FindBitExcitLevel(iLutRef, CurrentDets(:,idet), max_calc_ex_level)
            if(tRef_Not_HF) then
                ex_level_to_hf = FindBitExcitLevel (iLutHF_true, CurrentDets(:,idet), &
                     max_calc_ex_level)
            else
                ex_level_to_hf = ex_level_to_ref
            endif

            tParentIsDeterm = check_determ_flag(CurrentDets(:,idet))
            tParentUnoccupied = IsUnoccDet(parent_sign)

            if(tParentIsDeterm) then
               indices_of_determ_states(determ_index) = idet
               partial_determ_vecs(:,determ_index) = parent_sign
               determ_index = determ_index + 1
               if(tParentUnoccupied) cycle
            endif

            ! If this slot is unoccupied (and also not a core determinant) then add it to
            ! the list of free slots and cycle.
            if (tParentUnoccupied) then
                iEndFreeSlot = iEndFreeSlot + 1
                FreeSlot(iEndFreeSlot) = idet
                temp_iendfreeslot = temp_iendfreeslot + 1
                temp_freeslot(temp_iendfreeslot) = idet
                cycle
            end if

            ! get new population of observed orbitals
            if(numSnapShotOrbs .gt. 0) call makePopSnapshot(idet)

            ! The current diagonal matrix element is stored persistently.
            parent_hdiag = det_diagH(idet)

            if (tTruncInitiator) call CalcParentFlag(idet, parent_flags)

            ! do i need to calc. the energy contributions in the rt-fciqmc?
            ! leave it for now.. and figure out later..
            call SumEContrib (nI_parent, ex_level_to_ref, parent_sign, CurrentDets(:,idet), &
                               parent_hdiag, 1.0_dp, tPairedReplicas, idet)
            
            ! If we're on the Hartree-Fock, and all singles and
            ! doubles are in the core space, then there will be
            ! no stochastic spawning from this determinant, so
            ! we can the rest of this loop.
            if (ss_space_in%tDoubles .and. ex_level_to_hf == 0 .and. tDetermHFSpawning) cycle

            ! i dont thin i need this below..
!             if (tAllSymSectors) then
!                 ms_parent = return_ms(CurrentDets(:,idet))
!                 nOccAlpha = (nel+ms_parent)/2
!                 nOccBeta = (nel-ms_parent)/2
!             end if
! 
            ! If this condition is not met (if all electrons have spin up or all have spin down)
            ! then there will be no determinants to spawn to, so don't attempt spawning.
            ! thats a really specific condition.. shouldnt be checked each 
            ! cycle.. since this is only input dependent..
            do ireplica = 1, lenof_sign

                call decide_num_to_spawn(parent_sign(ireplica), parent_hdiag, AvMCExcits, nspawn)
                !call merge_spawn(nspawn,prefactor)
                do ispawn = 1, nspawn

                    ! Zero the bit representation, to ensure no extraneous data gets through.
                    ilut_child = 0_n_int

                    call generate_excitation (nI_parent, CurrentDets(:,idet), nI_child, &
                                        ilut_child, exFlag, ic, ex, tParity, prob, &
                                        HElGen, fcimc_excit_gen_store)

                    ! If a valid excitation.
                    if (.not. IsNullDet(nI_child)) then

                        call encode_child (CurrentDets(:,idet), ilut_child, ic, ex)
                        if (tUseFlags) ilut_child(nOffFlag) = 0_n_int

                        if (tSemiStochastic) then
                           break = check_semistoch_flags(ilut_child, nI_child, tParentIsDeterm)
                           if(break) cycle
                        endif

                        child_sign = attempt_create (nI_parent, CurrentDets(:,idet), parent_sign, &
                                            nI_child, ilut_child, prob, HElGen, ic, ex, tParity, &
                                            ex_level_to_ref, ireplica, unused_sign, unused_rdm_real)
                        child_sign = child_sign*prefactor
                    else
                        child_sign = 0.0_dp
                    end if

                    ! If any (valid) children have been spawned.
                    if ((any(child_sign /= 0)) .and. (ic /= 0) .and. (ic <= 2)) then

                        ! not quite sure about how to collect the child
                        ! stats in the new rt-fciqmc ..
                        call new_child_stats (iter_data_fciqmc, CurrentDets(:,idet), &
                                              nI_child, ilut_child, ic, ex_level_to_ref, &
                                              child_sign, parent_flags, ireplica)

                        call create_particle_with_hash_table (nI_child, ilut_child, child_sign, &
                                                               ireplica, CurrentDets(:,idet), &
                                                               iter_data_fciqmc)

                    end if ! If a child was spawned.

                end do ! Over mulitple particles on same determinant.

            end do ! Over the replicas on the same determinant.

            ! If this is a core-space determinant then the death step is done in
            ! determ_projection.

            ! in the 2nd RK loop of the real-time fciqmc the death step 
            ! should not act on the currently looped over walker list y(n)+k1/2
            ! but should be added to the spawned list k2, to then apply it to 
            ! the original y(n) + k2 
            ! the iEndFreeSlot variable also gets influenced in the 
            ! death-step (duh) -> so keep count of the temp iEndFreeSlot above
            if (.not. tParentIsDeterm) then
                call walker_death_realtime (iter_data_fciqmc, nI_parent,  &
                    CurrentDets(:,idet), parent_hdiag, parent_sign, idet, ex_level_to_ref)
            end if

        end do ! Over all determinants.

        if (tSemiStochastic) call real_time_determ_projection()

        ! also update the temp. variables to reuse in y(n) + k2 comb.
        ! this should be done before the annihilaiton step, as there these 
        ! values get changed! 
        ! have to do that above in the loop as it gets influenced in the 
        ! death-step too, 


        ! this is the original number of dets.
        TotWalkersNew = int(TotWalkers, sizeof_int)

        ! have to call end_iter_stats to get correct acceptance rate
!         call end_iter_stats(TotWalkersNew)
        ! but end iter stats for me is only uses to get SumWalkersCyc .. 
        do run = 1, inum_runs
           SumWalkersCyc_1(run) = SumWalkersCyc_1(run) + &
                sum(TotParts(min_part_type(run):max_part_type(run)))
        enddo

        ! the number TotWalkersNew changes below in annihilation routine
        ! Annihilation is done after loop over walkers

        call DirectAnnihilation (TotWalkersNew, iter_data_fciqmc, .false.)


        TotWalkers = int(TotWalkersNew, sizeof_int)

    end subroutine first_real_time_spawn

    subroutine perform_real_time_iteration()
        ! routine which performs one real-time fciqmc iteration
      implicit none
        character(*), parameter :: this_routine = "perform_real_time_iteration"
        
        integer :: TotWalkersNew, run
        real(dp) :: tmp_sign(lenof_sign), tau_real_tmp, tau_imag_tmp
        logical :: both, rkone, rktwo
        rkone = .true.
        rktwo = .true.
        both = .false.
        if(rkone .and. rktwo) both = .true.
        ! 0)
        ! do all the necessary preperation(resetting pointers etc.)
        ! concerning the statistics: i could use the "normal" iter_data
        ! for the first spawn, except change it for the new death-step, as 
        ! there the particles also change from Re <-> Im 

        call init_real_time_iteration(iter_data_fciqmc, second_spawn_iter_data)
        ! 1)
        ! do a "normal" spawning step and combination to y(n) + k1/2
        ! into CurrentDets: 
if(rkone) then
   if(iProcIndex == root .and. .false.) then
      print *, "TotParts and totDets before first spawn: ", TotParts, TotWalkers
   endif

if(both) then
        tau_real_tmp = tau_real
        tau_imag_tmp = tau_imag
        tau_real = tau_real/2.0
        tau_imag = tau_imag/2.0
endif
        call first_real_time_spawn()
if(both) then
        tau_real = tau_real_tmp
        tau_imag = tau_imag_tmp
endif

if(iProcIndex == root .and. .false.) then
        print *, "ValidSpawnedList", ValidSpawnedList
        print *, "TotParts and totDets after first spawn: ", TotParts, TotWalkers
        print *, "=========================="
     endif

#ifdef __DEBUG
        call check_update_growth(iter_data_fciqmc,"Error in first RK step")
#endif

        ! for now update the iter data here, although in the final 
        ! implementation i only should do that after the 2nd RK step
        ! or keep track of two different iter_datas for first and second 
        ! spawning..

        call update_iter_data(iter_data_fciqmc)
endif
if(rktwo) then
        ! 2)
        ! reset the spawned list and do a second spawning step to create 
        ! the spawend list k2 
        ! but DO NOT yet recombine with stored walker list 
        ! if i want to keep track of the two distinct spawns in 2 different 
        ! iter_datas i probably have to reset some values before the 
        ! second spawn.. to then keep the new values in the 2nd list

        call reset_spawned_list() 

        NoBorn = 0
        Annihilated = 0
        NoDied = 0

        ! create a second spawned list from y(n) + k1/2
        ! have to think to exclude death_step here and store this 
        ! information into the spawned k2 list..
        ! quick solution would be to loop again over reloaded y(n)
        ! and do a death step for wach walker
         call second_real_time_spawn()

        ! 3) 
        ! reload stored temp_det_list y(n) into CurrentDets 
        ! have to figure out how to effectively save the previous hash_table
        ! or maybe just use two with different types of update functions..

if(both) then
        call reload_current_dets()
endif

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

        ! UPDATE! have changed the 2nd diagonal event, so these particles get 
        ! stored in a seperate DiagParts array -> so i have to do two 
        ! annihilation events, first with the diagonal list and then with 
        ! the actual spawned particles, to best mimick the old algorithm and 
        ! also to correctly keep the stats of the events! 

        TotWalkersNew = int(TotWalkers, sizeof_int)

        ! also have to set the SumWalkersCyc before the "proper" annihilaiton 
        do run = 1, inum_runs
           SumWalkersCyc(run) = SumWalkersCyc(run) + &
                sum(TotParts(min_part_type(run):max_part_type(run)))
        enddo

        call DirectAnnihilation_diag(TotWalkersNew, second_spawn_iter_data)
         TotWalkersNew = int(TotWalkersNew, sizeof_int)
        
        ! and then do the "normal" annihilation with the SpawnedParts array!
        ! Annihilation is done after loop over walkers
        call DirectAnnihilation (TotWalkersNew, second_spawn_iter_data, .false.)

#ifdef __DEBUG
        call check_update_growth(second_spawn_iter_data,"Error in second RK step")
#endif


        ! for debugging comfort: if the second step is to be used on its own
        if(.not. both) then
           NoatHF = NoatHF_1
        endif

        TotWalkers = int(TotWalkersNew, sizeof_int)

        ! also do the update on the second_spawn_iter_data to combine both of 
        ! them outside this function 

        call update_iter_data(second_spawn_iter_data)
else
   SumWalkersCyc = SumWalkersCyc_1
endif

    end subroutine perform_real_time_iteration

    subroutine perform_verlet_iteration
      implicit none
      integer :: TotWalkersNew

      call init_verlet_iteration()
      call update_elapsed_time()

      ! load H^2 psi into spawnedParts

      call obtain_h2_psi()
      
      ! merge delta_psi and spawnedParts into the new delta_psi (which is stored in
      ! spawnedParts)

      call update_delta_psi()
      
      ! merge delta_psi (now spawnedParts) into CurrentDets 
      ! We need to cast TotWalkers to a regular int to pass it to the annihilation
      ! as it is modified, we need to pass an lvalue and cannot just pass int(TotWalkers)
      TotWalkersNew = int(TotWalkers,sizeof_int)
      call end_iter_stats(TotWalkersNew)
      ! for semistochastic method, we add in the core -> core spawns
      ! if(tSemiStochastic) call deterministic_annihilation(iter_data_fciqmc)
      call AnnihilateSpawnedParts(spawnBufSize,TotWalkersNew,iter_data_fciqmc)
      ! Updating the statistics is usually done in the annihilation, but since we
      ! explicitly carry out the annihilation, this has to be included explicitly
      ! (We can not use DirectAnnihilation because we need the communicated spawnedParts
      !  in between to update delta_psi)
      call CalcHashTableStats(TotWalkersNew,iter_data_fciqmc)
      TotWalkers = TotWalkersNew
      
    end subroutine perform_verlet_iteration

end module real_time

! wrapper (dont know why this is necessary quite..)
subroutine perform_real_time_fciqcm_wrap
    use real_time, only: perform_real_time_fciqmc
    implicit none

    call perform_real_time_fciqmc

end subroutine perform_real_time_fciqcm_wrap
