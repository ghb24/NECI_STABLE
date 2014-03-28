#include "macros.h"

! This module contains the routine which is used when NECI is run with the option
! determ-proj. This takes the HF state as an initial state and applies the
! deterministic part of the full projection operator *only*. So, it is different
! to a normal FCIQMC run where the rest of the space is projected out stochastically.
! Only the iteration number and energy is output.

! This is mainly useful for debugging.

module determ_proj

    use bit_rep_data, only: flag_deterministic, NIfD
    use bit_reps, only: test_flag
    use CalcData, only: NMCyc, tSemiStochastic
    use constants
    use DetBitOps, only: DetBitLT
    use FciMCData, only: HFDet, ilutHF, iHFProc, CurrentDets, core_hamiltonian, &
                         determ_proc_sizes, determ_space_size, partial_determ_vector, &
                         full_determ_vector, determ_proc_indices, TotWalkers
    use Parallel_neci, only: iProcIndex, MPIAllGatherV, MPISum
    use semi_stoch_procs, only: deterministic_projection

    implicit none

contains

    subroutine perform_determ_proj()

        integer :: counter, iter, comp, hf_index, ierr
        integer(int64) :: i
        real(dp), allocatable, dimension(:) :: wavefunction
        real(dp), allocatable, dimension(:) :: ham_times_hf
        real(dp) :: energy_num, energy_denom, tot_e_num, tot_e_denom

        if ((.not. tSemiStochastic) .or. (.not. allocated(core_hamiltonian))) &
            call stop_all("determ_projection_only", "You must use the semi-stochastic &
                &option and define a core space to use the determ-proj option.") 

        allocate(wavefunction(determ_proc_sizes(iProcIndex)))
        allocate(ham_times_hf(determ_proc_sizes(iProcIndex)))

        write(6,'()')
        write(6,'(a83)') "Performing a deterministic projection using the defined &
                         &semi-stochastic core space."
        write(6,'()')

        iter = 1
        energy_denom = 0.0_dp

        ! Find the index of the HF state in the vectors of the HF processor.
        if (iProcIndex == iHFProc) then
            counter = 0
            do i = 1, TotWalkers
                if (test_flag(CurrentDets(:,i),flag_deterministic)) then
                    counter = counter + 1
                    comp = DetBitLT(CurrentDets(:,i), ilutHF, NIfD, .false.)
                    if (comp == 0) hf_index = counter
                end if
            end do
        end if

        wavefunction = 0.0_dp
        if (iProcIndex == iHFProc) wavefunction(hf_index) = 1.0_dp

        call MPIAllGatherV(wavefunction, full_determ_vector, determ_proc_sizes, &
                            determ_proc_indices)

        call dgemv('N', &
                   determ_proc_sizes(iProcIndex), &
                   determ_space_size, &
                   1.0_dp, &
                   core_hamiltonian, &
                   determ_proc_sizes(iProcIndex), &
                   full_determ_vector, &
                   1, &
                   0.0_dp, &
                   ham_times_hf, &
                   1)

        write(6,'(a9,7X,a6)') "Iteration", "Energy"
        call neci_flush(6)

        do while(iter <= NMCyc .or. NMCyc == -1)

            partial_determ_vector = wavefunction

            call deterministic_projection()

            wavefunction = wavefunction + partial_determ_vector

            energy_num = dot_product(ham_times_hf, wavefunction)
            if (iProcIndex == iHFProc) energy_denom = wavefunction(hf_index)

            call MPISum(energy_num, tot_e_num)
            call MPISum(energy_denom, tot_e_denom)

            write(6,'(i9,7X,f13.10)') iter, tot_e_num/tot_e_denom
            call neci_flush(6)

            iter = iter + 1

        end do

        deallocate(wavefunction)
        deallocate(ham_times_hf)

    end subroutine perform_determ_proj

end module determ_proj
