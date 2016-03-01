#include "macros.h"

module real_time

    use real_time_init, only: init_real_time_calc_single
    use real_time_procs, only: save_current_dets, reset_spawned_list, &
                               reload_current_dets, walker_death_realtime, &
                               walker_death_spawn, attempt_die_realtime, &
                               create_diagonal_as_spawn, count_holes_in_currentDets, &
                               DirectAnnihilation_diag
    use real_time_data, only: gf_type, temp_freeslot, temp_iendfreeslot, wf_norm, &
                              pert_norm, second_spawn_iter_data, runge_kutta_step
    use CalcData, only: pops_norm, tTruncInitiator, tPairedReplicas, ss_space_in, &
                        tDetermHFSpawning, AvMCExcits, tSemiStochastic, StepsSft, &
                        tChangeProjEDet, tInstGrowthRate
    use FciMCData, only: pops_pert, walker_time, iter, ValidSpawnedList, &
                         spawn_ht, FreeSlot, iStartFreeSlot, iEndFreeSlot, &
                         fcimc_iter_data, InitialSpawnedSlots, iter_data_fciqmc, &
                         TotWalkers, fcimc_excit_gen_store, ilutRef, max_calc_ex_level, &
                         iLutHF_true, indices_of_determ_states, partial_determ_vecs, &
                         exFlag, CurrentDets, TotParts, ilutHF
    use kp_fciqmc_data_mod, only: overlap_pert
    use timing_neci, only: set_timer, halt_timer
    use FciMCParMod, only: rezero_iter_stats_each_iter
    use hash, only: clear_hash_table
    use constants, only: int64, sizeof_int, n_int, lenof_sign, dp, EPS, inum_runs
    use AnnihilationMod, only: DirectAnnihilation
    use bit_reps, only: extract_bit_rep
    use SystemData, only: nel, tRef_Not_HF, tAllSymSectors, nOccAlpha, nOccBeta, &
                          nbasis
    use DetBitOps, only: FindBitExcitLevel, return_ms
    use semi_stoch_procs, only: check_determ_flag, is_core_state, determ_projection
    use global_det_data, only: det_diagH
    use fcimc_helper, only: CalcParentFlag, decide_num_to_spawn, &
                            create_particle_with_hash_table, walker_death, &
                            SumEContrib
    use procedure_pointers, only: generate_excitation, encode_child, &
                                  attempt_create, new_child_stats
    use bit_rep_data, only: tUseFlags, nOffFlag, niftot, extract_sign
    use bit_reps, only: set_flag, flag_deterministic, flag_determ_parent
    use fcimc_iter_utils, only: update_iter_data
    use soft_exit, only: ChangeVars, tSoftExitFound
    use fcimc_initialisation, only: CalcApproxpDoubles
    use LoggingData, only: tPopsFile
    use PopsFileMod, only: WriteToPopsfileParOneArr

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

!         print *, " overlap: ",  sqrt(gf_overlap(0)) / sqrt(pops_norm) 
!         print *, " overlap2:", sqrt(gf_overlap(0)) / sqrt(wf_norm(0))
!         print *, pops_norm, pert_norm, wf_norm(0)

        ! main part of the initialization, concerning the walker lists 
        ! and operator application works now.. 
        ! so next step is to finish up the whole initialization to be able 
        ! to run a (mk)neci with the new walker dynamics provided by the 
        ! real-time Schroedinger equation

        print *, "toto: changeref: ", tChangeProjEDet
        print *, "toto: instgrowrate: ", tInstGrowthRate
        ! enter the main real-time fciqmc loop here
        fciqmc_loop: do while (.true.)


            ! update the overlap each time step.. could do that in some of 
            ! the loops down there too.. but for now do it here seperately 
            ! to make it more organized

            call update_gf_overlap()

            ! do also need to update the norm of the wavefunction, but this 
            ! i should definetly do in the first RK loop over the CurrentDets

            ! the timing stuff has to be done a bit differently in the 
            ! real-time fciqmc, since there are 2 spawing and annihilation 
            ! steps involved...
            call set_timer(walker_time)

            iter = iter + 1

            ! perform the actual iteration(excitation generation etc.) 
            call perform_real_time_iteration() 

            print *, "overlap: ", gf_overlap(1,iter-1) / wf_norm(iter-1), &
                gf_overlap(2,iter-1) / wf_norm(iter-1), gf_overlap(2,iter-1)

            ! update various variables.. time-step, initiator criteria etc.
            call update_real_time_iteration()

            ! update, print and log all the global variables and interesting 
            ! quantities
            call log_real_time_iteration() 

            ! check if somthing happpened to stop the iteration or something
            call check_real_time_iteration()

            if (tSoftExitFound) exit fciqmc_loop

!             if (iter == 1000) call stop_all(this_routine, "stop for now to avoid endless loop")

        end do fciqmc_loop
        
        if (tPopsFile) call WriteToPopsfileParOneArr(CurrentDets,TotWalkers)

    end subroutine perform_real_time_fciqmc

    subroutine update_real_time_iteration()
        ! routine to update certain global variables each loop iteration in 
        ! the real-time fciqmc 
        ! from the 2 distince spawn/death/cloning info stored in the 
        ! two iter_data vars, i have to combine the general updated 
        ! statistics for the actual time step
        character(*), parameter :: this_routine = "update_real_time_iteration"

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
        ! 
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

    end subroutine check_real_time_iteration

    subroutine init_real_time_iteration(iter_data, iter_data2)
        ! routine to reinitialize all the necessary variables and pointers 
        ! for a sucessful real-time fciqmc iteration
        type(fcimc_iter_data), intent(inout) :: iter_data
        type(fcimc_iter_data), intent(inout), optional :: iter_data2
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

        ! also reset the temporary variables
        ! this probably has not to be done, since at the end of the first 
        ! spawn i set it to the freeslot array anyway..
        ! with changed temp_ var. usage i have to reset them! 
        temp_freeslot(1:temp_iendfreeslot) = 0
        temp_iendfreeslot = 0

        ! Index for counting deterministic states.
!         n_determ_states = 1

        ! Clear the hash table for the spawning array.
        call clear_hash_table(spawn_ht)

        call rezero_iter_stats_each_iter(iter_data)

        if (present(iter_data2)) then
            call rezero_iter_stats_each_iter(iter_data2)
        end if

        ! additionaly i have to copy the CurrentDets array and all the 
        ! associated pointers and hashtable related stuff to the 
        ! temporary 2nd list 
        call save_current_dets() 

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
                   ms_parent, ireplica, nspawn, ispawn, nI_child(nel), ic, ex(2,2), &
                   ex_level_to_hf
        integer(n_int), pointer :: ilut_parent(:) 
        integer(n_int) :: ilut_child(0:niftot)
        real(dp) :: parent_sign(lenof_sign), parent_hdiag, prob, child_sign(lenof_sign), &
                    unused_sign(lenof_sign), unused_rdm_real, diag_sign(lenof_sign)
        logical :: tParentIsDeterm, tParentUnoccupied, tParity, tChildIsDeterm
        HElement_t(dp) :: HelGen
        integer(int64) :: TotWalkersNew

        ! declare this is the second runge kutta step
        runge_kutta_step = 2

        ! use part of nicks code, and remove the parts, that dont matter 
        do idet = 1, int(TotWalkers, sizeof_int)

            ! The 'parent' determinant from which spawning is to be attempted.
            ilut_parent => CurrentDets(:,idet)
            parent_flags = 0_n_int

            ! Indicate that the scratch storage used for excitation generation from the
            ! same walker has not been filled (it is filled when we excite from the first
            ! particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.

            call extract_bit_rep(ilut_parent, nI_parent, parent_sign, unused_flags, &
                                  fcimc_excit_gen_store)

            ex_level_to_ref = FindBitExcitLevel(iLutRef, ilut_parent, max_calc_ex_level)
            if(tRef_Not_HF) then
                ex_level_to_hf = FindBitExcitLevel (iLutHF_true, ilut_parent, max_calc_ex_level)
            else
                ex_level_to_hf = ex_level_to_ref
            endif

            tParentIsDeterm = check_determ_flag(ilut_parent)
            tParentUnoccupied = IsUnoccDet(parent_sign)

            ! If this determinant is in the deterministic space then store the relevant
            ! data in arrays for later use.
            ! do not use the deterministic option yet in the rt-fciqmc, as i
            ! am not yet sure how to do the deterministic projection with the 
            ! new -i H walker dynamics... todo!
!             if (tParentIsDeterm) then
!                 ! Store the index of this state, for use in annihilation later.
!                 indices_of_determ_states(n_determ_states) = idet
!                 ! Add the amplitude to the deterministic vector.
!                 partial_determ_vecs(:,n_determ_states) = parent_sign
!                 n_determ_states = n_determ_states + 1
! 
!                 ! The deterministic states are always kept in CurrentDets, even when
!                 ! the amplitude is zero. Hence we must check if the amplitude is zero
!                 ! and, if so, skip the state.
!                 if (tParentUnoccupied) cycle
!             end if

            ! If this slot is unoccupied (and also not a core determinant) then add it to
            ! the list of free slots and cycle.
            ! actually dont need to update this here since nothing gets merged
            ! at the end.. only k2 gets created
            if (tParentUnoccupied) then
                iEndFreeSlot = iEndFreeSlot + 1
                FreeSlot(iEndFreeSlot) = idet
                cycle
            end if

            ! The current diagonal matrix element is stored persistently.
            parent_hdiag = det_diagH(idet)

            if (tTruncInitiator) call CalcParentFlag(idet, parent_flags, parent_hdiag)

            ! do i need to calc. the energy contributions in the rt-fciqmc?
            ! leavi it for now.. and figure out later..
            ! definetly do not need it in the second spawn, since it is only 
            ! an intermediate list, and not an actual representation of the 
            ! wavefunction
!             call SumEContrib (nI_parent, ex_level_to_ref, parent_sign, ilut_parent, &
!                                parent_hdiag, 1.0_dp, tPairedReplicas, idet)

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
            do ireplica = 1, lenof_sign

                call decide_num_to_spawn(parent_sign(ireplica), AvMCExcits, nspawn)
                
                do ispawn = 1, nspawn

                    ! Zero the bit representation, to ensure no extraneous data gets through.
                    ilut_child = 0_n_int

                    call generate_excitation (nI_parent, ilut_parent, nI_child, &
                                        ilut_child, exFlag, ic, ex, tParity, prob, &
                                        HElGen, fcimc_excit_gen_store)

                    ! If a valid excitation.
                    if (.not. IsNullDet(nI_child)) then

                        call encode_child (ilut_parent, ilut_child, ic, ex)
                        if (tUseFlags) ilut_child(nOffFlag) = 0_n_int

                        if (tSemiStochastic) then
                            tChildIsDeterm = is_core_state(ilut_child, nI_child)

                            ! Is the parent state in the core space?
                            if (tParentIsDeterm) then
                                ! If spawning is both from and to the core space, cancel it.
                                if (tChildIsDeterm) cycle
                                call set_flag(ilut_child, flag_determ_parent)
                            else
                                if (tChildIsDeterm) call set_flag(ilut_child, flag_deterministic)
                            end if
                        end if

                        child_sign = attempt_create (nI_parent, ilut_parent, parent_sign, &
                                            nI_child, ilut_child, prob, HElGen, ic, ex, tParity, &
                                            ex_level_to_ref, ireplica, unused_sign, unused_rdm_real)

                    else
                        child_sign = 0.0_dp
                    end if

                    ! If any (valid) children have been spawned.
                    if ((any(child_sign /= 0)) .and. (ic /= 0) .and. (ic <= 2)) then

                        ! not quite sure about how to collect the child
                        ! stats in the new rt-fciqmc ..
                        call new_child_stats (second_spawn_iter_data, ilut_parent, &
                                              nI_child, ilut_child, ic, ex_level_to_ref, &
                                              child_sign, parent_flags, ireplica)

                        call create_particle_with_hash_table (nI_child, ilut_child, child_sign, &
                                                               ireplica, ilut_parent, second_spawn_iter_data)

                    end if ! If a child was spawned.

                end do ! Over mulitple particles on same determinant.

            end do ! Over the replicas on the same determinant.

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
                diag_sign = -attempt_die_realtime(nI_parent, parent_hdiag, &
                    parent_sign, ex_level_to_ref)
                
                if (any(abs(diag_sign) > EPS)) then

                    call create_diagonal_as_spawn(nI_parent, ilut_parent, &
                        diag_sign, second_spawn_iter_data)
                end if

            end if

        end do ! Over all determinants.

        if (tSemiStochastic) call determ_projection()

    end subroutine second_real_time_spawn

    subroutine first_real_time_spawn()
        ! routine which first loops over the CurrentDets array and creates the 
        ! first spawning list k1 and combines it to y(n) + k1/2
!         integer, intent(inout) :: n_determ_states
        character(*), parameter :: this_routine = "first_real_time_spawn"
        integer :: idet, parent_flags, nI_parent(nel), unused_flags, ex_level_to_ref, &
                   ms_parent, ireplica, nspawn, ispawn, nI_child(nel), ic, ex(2,2), &
                   ex_level_to_hf, i
        integer(n_int), pointer :: ilut_parent(:) 
        integer(n_int) :: ilut_child(0:niftot)
        real(dp) :: parent_sign(lenof_sign), parent_hdiag, prob, child_sign(lenof_sign), &
                    unused_sign(lenof_sign), unused_rdm_real, tmp_norm
        logical :: tParentIsDeterm, tParentUnoccupied, tParity, tChildIsDeterm
        HElement_t(dp) :: HelGen
        integer(int64) :: TotWalkersNew

        ! declare this is the first runge kutta step
        runge_kutta_step = 1

        ! also calculate the norm of the current y(n) 
        tmp_norm = 0.0_dp

        ! use part of nicks code, and remove the parts, that dont matter 
        do idet = 1, int(TotWalkers, sizeof_int)

            ! The 'parent' determinant from which spawning is to be attempted.
            ilut_parent => CurrentDets(:,idet)
            parent_flags = 0_n_int

            ! Indicate that the scratch storage used for excitation generation from the
            ! same walker has not been filled (it is filled when we excite from the first
            ! particle on a determinant).
            fcimc_excit_gen_store%tFilled = .false.

            call extract_bit_rep(ilut_parent, nI_parent, parent_sign, unused_flags, &
                                  fcimc_excit_gen_store)

            ex_level_to_ref = FindBitExcitLevel(iLutRef, ilut_parent, max_calc_ex_level)
            if(tRef_Not_HF) then
                ex_level_to_hf = FindBitExcitLevel (iLutHF_true, ilut_parent, max_calc_ex_level)
            else
                ex_level_to_hf = ex_level_to_ref
            endif

            tParentIsDeterm = check_determ_flag(ilut_parent)
            tParentUnoccupied = IsUnoccDet(parent_sign)

            ! If this determinant is in the deterministic space then store the relevant
            ! data in arrays for later use.
            ! for now, do not use deterministic option in the real-time fciqmc
            ! add that later.. todo!
!             if (tParentIsDeterm) then
!                 ! Store the index of this state, for use in annihilation later.
!                 indices_of_determ_states(n_determ_states) = idet
!                 ! Add the amplitude to the deterministic vector.
!                 partial_determ_vecs(:,n_determ_states) = parent_sign
!                 n_determ_states = n_determ_states + 1
! 
!                 ! The deterministic states are always kept in CurrentDets, even when
!                 ! the amplitude is zero. Hence we must check if the amplitude is zero
!                 ! and, if so, skip the state.
!                 if (tParentUnoccupied) cycle
!             end if

            ! If this slot is unoccupied (and also not a core determinant) then add it to
            ! the list of free slots and cycle.
            if (tParentUnoccupied) then
                iEndFreeSlot = iEndFreeSlot + 1
                FreeSlot(iEndFreeSlot) = idet

                ! also update the temporary variables here, as the also get 
                ! influenced in the death-step below
                temp_iendfreeslot = temp_iendfreeslot + 1
                temp_freeslot(temp_iendfreeslot) = idet

                cycle
            end if

            ! add up norm here: 
!             do i = 1, lenof_sign
!                 tmp_norm = tmp_norm + parent_sign(i)**2
!             end do
            tmp_norm = tmp_norm + sum(parent_sign**2)

            ! The current diagonal matrix element is stored persistently.
            parent_hdiag = det_diagH(idet)

            if (tTruncInitiator) call CalcParentFlag(idet, parent_flags, parent_hdiag)

            ! do i need to calc. the energy contributions in the rt-fciqmc?
            ! leave it for now.. and figure out later..
            call SumEContrib (nI_parent, ex_level_to_ref, parent_sign, ilut_parent, &
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
            do ireplica = 1, lenof_sign

                call decide_num_to_spawn(parent_sign(ireplica), AvMCExcits, nspawn)
                
                do ispawn = 1, nspawn

                    ! Zero the bit representation, to ensure no extraneous data gets through.
                    ilut_child = 0_n_int

                    call generate_excitation (nI_parent, ilut_parent, nI_child, &
                                        ilut_child, exFlag, ic, ex, tParity, prob, &
                                        HElGen, fcimc_excit_gen_store)

                    ! If a valid excitation.
                    if (.not. IsNullDet(nI_child)) then

                        call encode_child (ilut_parent, ilut_child, ic, ex)
                        if (tUseFlags) ilut_child(nOffFlag) = 0_n_int

                        if (tSemiStochastic) then
                            tChildIsDeterm = is_core_state(ilut_child, nI_child)

                            ! Is the parent state in the core space?
                            if (tParentIsDeterm) then
                                ! If spawning is both from and to the core space, cancel it.
                                if (tChildIsDeterm) cycle
                                call set_flag(ilut_child, flag_determ_parent)
                            else
                                if (tChildIsDeterm) call set_flag(ilut_child, flag_deterministic)
                            end if
                        end if

                        child_sign = attempt_create (nI_parent, ilut_parent, parent_sign, &
                                            nI_child, ilut_child, prob, HElGen, ic, ex, tParity, &
                                            ex_level_to_ref, ireplica, unused_sign, unused_rdm_real)

                    else
                        child_sign = 0.0_dp
                    end if

                    ! If any (valid) children have been spawned.
                    if ((any(child_sign /= 0)) .and. (ic /= 0) .and. (ic <= 2)) then

                        ! not quite sure about how to collect the child
                        ! stats in the new rt-fciqmc ..
                        call new_child_stats (iter_data_fciqmc, ilut_parent, &
                                              nI_child, ilut_child, ic, ex_level_to_ref, &
                                              child_sign, parent_flags, ireplica)

                        call create_particle_with_hash_table (nI_child, ilut_child, child_sign, &
                                                               ireplica, ilut_parent, iter_data_fciqmc)

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
                    ilut_parent, parent_hdiag, parent_sign, idet, ex_level_to_ref)
            end if

        end do ! Over all determinants.

        if (tSemiStochastic) call determ_projection()

        ! also update the temp. variables to reuse in y(n) + k2 comb.
        ! this should be done before the annihilaiton step, as there these 
        ! values get changed! 
        ! have to do that above in the loop as it gets influenced in the 
        ! death-step too, 

        ! this is the original number of dets.
        TotWalkersNew = int(TotWalkers, sizeof_int)

        ! the number TotWalkersNew changes below in annihilation routine
        ! Annihilation is done after loop over walkers
        call DirectAnnihilation (TotWalkersNew, iter_data_fciqmc, .false.)

        TotWalkers = int(TotWalkersNew, sizeof_int)

        ! also store the new norm information 
        wf_norm(iter-1) = sqrt(pert_norm * tmp_norm)

    end subroutine first_real_time_spawn

    subroutine perform_real_time_iteration()
        ! routine which performs one real-time fciqmc iteration
        character(*), parameter :: this_routine = "perform_real_time_iteration"

        integer :: idet !, n_determ_states 
        integer(int64) :: TotWalkersNew
        real(dp) :: tmp_sign(lenof_sign)

        ! 0)
        ! do all the necessary preperation(resetting pointers etc.)
        ! concerning the statistics: i could use the "normal" iter_data
        ! for the first spawn, except change it for the new death-step, as 
        ! there the particles also change from Re <-> Im 
        call init_real_time_iteration(iter_data_fciqmc, second_spawn_iter_data)
        ! 1)
        ! do a "normal" spawning step and combination to y(n) + k1/2
        ! into CurrentDets: 
        
        print *, "nruns: ", inum_runs

        print *, "TotParts and totDets before first spawn: ", TotParts, TotWalkers
        call extract_sign(CurrentDets(:,1), tmp_sign)
        print *, "hf occ before first spawn:", tmp_sign
        call first_real_time_spawn()
        call extract_sign(CurrentDets(:,1), tmp_sign)
        print *, "hf occ after first spawn:", tmp_sign
        print *, "TotParts and totDets after first spawn: ", TotParts, TotWalkers
        print *, "=========================="

        

        ! for now update the iter data here, although in the final 
        ! implementation i only should do that after the 2nd RK step
        ! or keep track of two different iter_datas for first and second 
        ! spawning..
        call update_iter_data(iter_data_fciqmc)
        ! 2)
        ! reset the spawned list and do a second spawning step to create 
        ! the spawend list k2 
        ! but DO NOT yet recombine with stored walker list 
        ! if i want to keep track of the two distinct spawns in 2 different 
        ! iter_datas i probably have to reset some values before the 
        ! second spawn.. to then keep the new values in the 2nd list
        call reset_spawned_list() 

        ! create a second spawned list from y(n) + k1/2
        ! have to think to exclude death_step here and store this 
        ! information into the spawned k2 list..
        ! quick solution would be to loop again over reloaded y(n)
        ! and do a death step for wach walker

!         iStartFreeSlot = 1
!         iEndFreeSlot = 0
!         FreeSlot = 0

        call second_real_time_spawn()
        call extract_sign(CurrentDets(:,1), tmp_sign)
        print *, "hf occ after second spawn:", tmp_sign

        ! 3) 
        ! reload stored temp_det_list y(n) into CurrentDets 
        ! have to figure out how to effectively save the previous hash_table
        ! or maybe just use two with different types of update functions..
        call reload_current_dets()

        call extract_sign(CurrentDets(:,1), tmp_sign)
        print *, "hf occ after reload:", tmp_sign

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

        call DirectAnnihilation_diag(TotWalkersNew, second_spawn_iter_data, .false.)

        call extract_sign(CurrentDets(:,1), tmp_sign)
        print *, "hf occ after second death:", tmp_sign
!         TotWalkersNew = int(TotWalkersNew, sizeof_int)

        ! and then do the "normal" annihilation with the SpawnedParts array!
        ! Annihilation is done after loop over walkers
        call DirectAnnihilation (TotWalkersNew, second_spawn_iter_data, .false.)

        call extract_sign(CurrentDets(:,1), tmp_sign)
        print *, "hf occ after second annihil:", tmp_sign
        TotWalkers = int(TotWalkersNew, sizeof_int)

    end subroutine perform_real_time_iteration

end module real_time

! wrapper (dont know why this is necessary quite..)
subroutine perform_real_time_fciqcm_wrap
    use real_time, only: perform_real_time_fciqmc
    implicit none

    call perform_real_time_fciqmc

end subroutine perform_real_time_fciqcm_wrap
