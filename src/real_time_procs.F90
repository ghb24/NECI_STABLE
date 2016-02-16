#include "macros.h"

! module containing useful functions and subroutines needed in the real-time
! implementation of the FCIQMC algotrithm

module real_time_procs
    use hash, only: hash_table_lookup, init_hash_table, clear_hash_table
    use SystemData, only: nel
    use real_time_data, only: gf_overlap, TotWalkers_orig, TotWalkers_pert
    use kp_fciqmc_data_mod, only: perturbed_ground, overlap_pert
    use constants, only: dp, lenof_sign, int64
    use bit_reps, only: decode_bit_det
    use bit_rep_data, only: extract_sign, nifdbo, niftot
    use FciMCData, only: CurrentDets, HashIndex, popsfile_dets, MaxWalkersPart, &
                         WalkVecDets, freeslot, spawn_ht, nhashes_spawn, MaxSpawned, &
                         iStartFreeSlot, iEndFreeSlot, ValidSpawnedList, &
                         InitialSpawnedSlots
    use perturbations, only: apply_perturbation
    use real_time_data, only: temp_det_list, temp_det_pointer, temp_det_hash, &
                              temp_freeslot
    use util_mod, only: int_fmt

    implicit none

contains

    subroutine save_current_dets() 
        ! routine to copy the currentDets array and all the associated 
        ! pointers an hashtable related quantities to the 2nd temporary 
        ! list, from which the first spawn and y(n) + k1/2 addition is done 
        ! and the k2 spawing list is created to then use CurrentDets to go to
        ! the next time-step y(n+1) = y(n) + k2
        ! new idea to reduce the amount of rewriting of old routines:
        ! save the CurrentDets and all associated quantities in the temporary
        ! variables, then normally work on CurrentDets to create the k1 
        ! spawning and also the CurrentDets = y(n) + k1 / 2 combination
        ! then spawn again from this list to create k2 spawning array
        ! but then before combing reload the saved CurrentDets from the 
        ! temporary variable
        ! probably have to think about parallelism issues here..
        character(*), parameter :: this_routine = "save_current_dets"

        ! save the WalkVecDets variable, i think thats the only necessary 
        ! variable, the pointers don't count 
        temp_det_list = WalkVecDets
        
        ! for now also store the pointer, but thats not needed i guess
        temp_det_pointer => temp_det_list

        ! also need the hash table and the freeslot i guess
        temp_det_hash = HashIndex

        ! and the freeslot.. although this one gets reinitialized to 0 
        ! every iteration or not? yeah it is.. so i only have to reset it 
        ! twice in the rt-fciqmc before the y(n) + k2 combination
        ! do that in the reload_current_dets routine!
        ! same with n_determ_states var.
        
        ! test 
        print *, WalkVecDets(:,1:10)
        print *, temp_det_list(:,1:10)

    end subroutine save_current_dets

    subroutine reload_current_dets()
        ! routine to reload the saved y(n) CurrentDets array for the final
        ! y(n) + k/2 combination to move to the next time step y(n+1)
        ! have also to think about the death-step and annihilation step 
        ! influences.. maybe have to write new death/born routines to split 
        ! that from the spawned list creation..
        character(*), parameter :: this_routine = "reload_current_dets" 

        ! copy the list
        WalkVecDets = temp_det_list

        ! and point to it
        CurrentDets => WalkVecDets 

        ! and the hash
        HashIndex = temp_det_hash

        ! for correct load Balancin i also have to reset the the freeslot var.
        FreeSlot(1:iEndFreeSlot) = 0
        iStartFreeSlot = 1
        iEndFreeSlot = 0

        ! not sure if i have to reset this here, or in the routine to reset 
        ! the spawned list
        ! i think that has to be done in the reset_spawned_list
!         n_determ_states = 1

        ! test 
        print *, WalkVecDets(:,1:10)
        print *, temp_det_list(:,1:10)

    end subroutine reload_current_dets

    subroutine reset_spawned_list(n_determ_states)
        ! also need a routine to reset the spawned lists before the second 
        ! spawning step for the 2nd order RK method in the rt-fciqmc
        integer, intent(out) :: n_determ_states
        character(*), parameter :: this_routine = "reset_spawned_list"

        ! Reset positions to spawn into in the spawning array.
        ValidSpawnedList = InitialSpawnedSlots

        ! Clear the hash table for the spawning array.
        call clear_hash_table(spawn_ht)

        ! think i have to reset the deterministic counter here
        n_determ_states = 1

    end subroutine reset_spawned_list



    subroutine setup_temp_det_list()
        ! setup the second list to temporaly store the list of determinants
        ! necessary in the real-time fciqmc list
        ! determine the necessary size from the already setup CurrentDets
        character(*), parameter :: this_routine = "setup_temp_det_list"
        integer :: ierr, tmp_siz1, tmp_siz2, i, spawn_ht_mem

        tmp_siz1 = size(WalkVecDets,dim=1)
        tmp_siz2 = size(WalkVecDets,dim=2)

        ! allocate the array
        allocate(temp_det_list(0:tmp_siz1-1,tmp_siz2), stat = ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error in allocation")

        ! and init it
        temp_det_list(0:tmp_siz1,1:tmp_siz2) = 0
        
        ! and point to it 
        temp_det_pointer => temp_det_list

        ! and also allocate the hash-table
        tmp_siz1 = size(HashIndex)
        allocate(temp_det_hash(tmp_siz1), stat = ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error in allocation")

        ! and initialize it to 0
        do i = 1, tmp_siz1
            temp_det_hash(i)%ind = 0
        end do

        ! also need a corresponding freeslot array
        tmp_siz1 = size(freeslot)
        allocate(temp_freeslot(tmp_siz1), stat = ierr)
        if(ierr.ne.0) call stop_all(this_routine,"Error in allocation")

        ! and init it
        temp_freeslot = 0

        ! also use the spawn_ht hash table, so also allocate it here! 

        ! Allocate the hash table to the spawning array.
        ! The number of MB of memory required to allocate spawn_ht.
        ! Each node requires 16 bytes.
        nhashes_spawn = int(0.8_dp * real(MaxSpawned,dp))
        spawn_ht_mem = nhashes_spawn*16/1000000
        write(6,'(a78,'//int_fmt(spawn_ht_mem,1)//')') "About to allocate hash table to the spawning array. &
                                       &Memory required (MB):", spawn_ht_mem
        write(6,'(a13)',advance='no') "Allocating..."; call neci_flush(6)
        allocate(spawn_ht(nhashes_spawn), stat=ierr)
        if (ierr /= 0) then
            write(6,'(1x,a11,1x,i5)') "Error code:", ierr
            call stop_all(this_routine, "Error allocating spawn_ht array.")
        else
            write(6,'(1x,a5)') "Done."
            write(6,'(a106)') "Note that the hash table uses linked lists, and the memory usage will &
                              &increase as further nodes are added."
        end if

        call init_hash_table(spawn_ht)

    end subroutine setup_temp_det_list

    ! subroutine to calculate the overlap of the current y(t) = a_j(a^+_j)(t)y(0)>
    ! time evolved wavefunction to the saved <y(0)|a^+_i(a_i) 
    subroutine update_gf_overlap() 
        ! this routine only deals with globally defined variables
        integer :: idet, nI(nel), det_ind, hash_val
        real(dp) :: overlap, real_sign_1(lenof_sign), real_sign_2(lenof_sign)
        logical :: tDetFound

        overlap = 0.0_dp

        do idet = 1, size(perturbed_ground, dim = 2)

            call extract_sign(perturbed_ground(:,idet), real_sign_1) 

            ! why should i store zero valued dets in the perturbed_gs? 
            ! and i have to expand that to complex code later on..
            ! although.. the perturbed ground state wavefunction is always
            ! real-valued -> so i only need the actual real-part of the 
            ! time evolved wavefunction! 
            if (IsUnoccDet(real_sign_1)) cycle

            call decode_bit_det(nI, perturbed_ground(:,idet))

            ! search for the hash table associated with the time evolved 
            ! wavefunction -> is this already initialized correctly? 
            call hash_table_lookup(nI, perturbed_ground(:,idet), nifdbo, &
                HashIndex, CurrentDets, det_ind, hash_val, tDetFound)

            if (tDetFound) then
                call extract_sign(CurrentDets(:,det_ind), real_sign_2)

                overlap = overlap + real_sign_1(1) * real_sign_2(1)
                
            end if
        end do

        ! need the timestep here... or the cycle of the current real-time loop
        gf_overlap(0) = overlap 

    end subroutine update_gf_overlap

    subroutine create_perturbed_ground()
        ! routine to create from the already read in popsfile info in 
        ! popsfile_dets the left hand <y(0)| by applying the corresponding 
        ! creation or annihilation operator
        character(*), parameter :: this_routine = "create_perturbed_ground"
        integer(int64) :: tmp_totwalkers
        integer :: ierr

        tmp_totwalkers = TotWalkers_orig

        print *, "Creating the wavefunction to projected on!"
        print *, "Initial number of walkers: ", tmp_totwalkers


        allocate(perturbed_ground(0:niftot,TotWalkers_orig), stat = ierr)

        call apply_perturbation(overlap_pert(1), tmp_totwalkers, popsfile_dets,&
            perturbed_ground)

        TotWalkers_pert = int(tmp_totwalkers, int64)

        print *, "Walkers remaining in perturbed ground state:" , TotWalkers_pert

        ! also need to create and associated hash table to this list 
!         call clear_hash_table(perturbed_ground_ht)
        ! or maybe not... lets see later on..


    end subroutine create_perturbed_ground

end module real_time_procs
