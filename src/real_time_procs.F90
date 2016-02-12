#include "macros.h"

! module containing useful functions and subroutines needed in the real-time
! implementation of the FCIQMC algotrithm

module real_time_procs
    use hash, only: hash_table_lookup
    use SystemData, only: nel
    use real_time_data, only: gf_overlap, TotWalkers_orig, TotWalkers_pert
    use kp_fciqmc_data_mod, only: perturbed_ground, overlap_pert
    use constants, only: dp, lenof_sign, int64
    use bit_reps, only: decode_bit_det
    use bit_rep_data, only: extract_sign, nifdbo, niftot
    use FciMCData, only: CurrentDets, HashIndex, popsfile_dets
    use perturbations, only: apply_perturbation

    implicit none

contains

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
        gf_overlap(1) = overlap 

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
