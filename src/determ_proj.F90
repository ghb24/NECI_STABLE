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
    use CalcData, only: NMCyc, tSemiStochastic, tOrthogonaliseReplicas
    use CalcData, only: tau, DiagSft
    use constants
    use DetBitOps, only: DetBitLT
    use FciMCData, only: HFDet, ilutHF, iRefProc, CurrentDets, &
                         TotWalkers, core_run
    use Parallel_neci, only: iProcIndex, MPIAllGatherV, MPISum
    use semi_stoch_procs, only: determ_projection, determ_proj_approx
    use core_space_util, only: cs_replicas
    implicit none

contains

    subroutine perform_determ_proj()

        integer :: counter, iter, comp, hf_index, ierr
        integer(int64) :: i, j
        real(dp), allocatable, dimension(:) :: wavefunction
        real(dp), allocatable, dimension(:) :: wavefunction_old
        real(dp), allocatable, dimension(:) :: ham_times_hf
        real(dp), allocatable, dimension(:) :: ham_times_wf
        real(dp) :: var_e_num, var_e_denom, energy_num, energy_denom
        real(dp) :: tot_var_e_num, tot_var_e_denom, tot_e_num, tot_e_denom
        character(*), parameter :: this_routine = 'perform_determ_proj'

        associate(rep => cs_replicas(core_run))

            if ((.not. tSemiStochastic) .or. (.not. allocated(rep%sparse_core_ham))) &
                call stop_all(this_routine, "You must use the semi-stochastic &
                    &option and define a core space to use the determ-proj option.")

            allocate(wavefunction(rep%determ_sizes(iProcIndex)))
            allocate(wavefunction_old(rep%determ_sizes(iProcIndex)))
            allocate(ham_times_hf(rep%determ_sizes(iProcIndex)))
            allocate(ham_times_wf(rep%determ_sizes(iProcIndex)))

            write(6, '()')
            write(6, '(a83)') "Performing a deterministic projection using the defined &
                &semi-stochastic core space."
            write(6, '()')

            iter = 1
            energy_denom = 0.0_dp

            ! Find the index of the HF state in the vectors of the HF processor.
            ASSERT(.not. tOrthogonaliseReplicas)
            if (iProcIndex == iRefProc(1)) then
                counter = 0
                do i = 1, TotWalkers
                    if (test_flag(CurrentDets(:, i), flag_deterministic(core_run))) then
                        counter = counter + 1
                        comp = DetBitLT(CurrentDets(:, i), ilutHF, NIfD)
                        if (comp == 0) hf_index = counter
                    end if
                end do
            end if

            wavefunction = 0.0_dp
            ASSERT(.not. tOrthogonaliseReplicas)
            if (iProcIndex == iRefProc(1)) wavefunction(hf_index) = 1.0_dp

            call MPIAllGatherV(wavefunction, rep%full_determ_vecs(1, :), rep%determ_sizes, &
                               rep%determ_displs)

            rep%partial_determ_vecs = 0.0_dp
            ham_times_hf = 0.0_dp

            do i = 1, rep%determ_sizes(iProcIndex)
                do j = 1, rep%sparse_core_ham(i)%num_elements
                    ham_times_hf(i) = ham_times_hf(i) + &
                                      rep%sparse_core_ham(i)%elements(j) * rep%full_determ_vecs(1, rep%sparse_core_ham(i)%positions(j))
                end do
            end do

            write(6, '(a11,7X,a12,7X,a11)') "# Iteration", "Proj. Energy", "Var. Energy"
            call neci_flush(6)

            do while (iter <= NMCyc .or. NMCyc == -1)

                wavefunction_old = wavefunction

                rep%partial_determ_vecs(1, :) = wavefunction

                call determ_projection()

                ham_times_wf = -rep%partial_determ_vecs(1, :) / tau + (DiagSft(1) * wavefunction)

                wavefunction = wavefunction + rep%partial_determ_vecs(1, :)

                energy_num = dot_product(ham_times_hf, wavefunction)
                ASSERT(.not. tOrthogonaliseReplicas)
                if (iProcIndex == iRefProc(1)) energy_denom = wavefunction(hf_index)

                var_e_num = dot_product(ham_times_wf, wavefunction_old)
                var_e_denom = dot_product(wavefunction_old, wavefunction_old)

                call MPISum(var_e_num, tot_var_e_num)
                call MPISum(var_e_denom, tot_var_e_denom)
                call MPISum(energy_num, tot_e_num)
                call MPISum(energy_denom, tot_e_denom)

                write(6, '(i9,7X,f13.10,7X,f13.10)') iter, tot_e_num / tot_e_denom, tot_var_e_num / tot_var_e_denom
                call neci_flush(6)

                iter = iter + 1

            end do

            deallocate(wavefunction)
            deallocate(ham_times_hf)
            deallocate(ham_times_wf)
        end associate

    end subroutine perform_determ_proj

    subroutine perform_determ_proj_approx_ham()

        integer :: counter, iter, comp, hf_index, ierr
        integer(int64) :: i, j
        real(dp), allocatable, dimension(:) :: wavefunction
        real(dp), allocatable, dimension(:) :: wavefunction_old
        real(dp), allocatable, dimension(:) :: ham_times_hf
        real(dp), allocatable, dimension(:) :: ham_times_wf
        real(dp) :: var_e_num, var_e_denom, energy_num, energy_denom
        real(dp) :: tot_var_e_num, tot_var_e_denom, tot_e_num, tot_e_denom
        character(*), parameter :: this_routine = 'perform_determ_proj'

        associate(rep => cs_replicas(core_run))

            if ((.not. tSemiStochastic) .or. (.not. allocated(rep%sparse_core_ham))) &
                call stop_all(this_routine, "You must use the semi-stochastic &
                    &option and define a core space to use the determ-proj option.")

            allocate(wavefunction(rep%determ_sizes(iProcIndex)))
            allocate(wavefunction_old(rep%determ_sizes(iProcIndex)))
            allocate(ham_times_hf(rep%determ_sizes(iProcIndex)))
            allocate(ham_times_wf(rep%determ_sizes(iProcIndex)))

            write(6, '()')
            write(6, '(a83)') "Performing a deterministic projection using the defined &
                             &semi-stochastic core space."
            write(6, '()')

            iter = 1
            energy_denom = 0.0_dp

            ! Find the index of the HF state in the vectors of the HF processor.
            ASSERT(.not. tOrthogonaliseReplicas)
            if (iProcIndex == iRefProc(1)) then
                counter = 0
                do i = 1, TotWalkers
                    if (test_flag(CurrentDets(:, i), flag_deterministic(core_run))) then
                        counter = counter + 1
                        comp = DetBitLT(CurrentDets(:, i), ilutHF, NIfD)
                        if (comp == 0) hf_index = counter
                    end if
                end do
            end if

            wavefunction = 0.0_dp
            ASSERT(.not. tOrthogonaliseReplicas)
            if (iProcIndex == iRefProc(1)) wavefunction(hf_index) = 1.0_dp

            call MPIAllGatherV(wavefunction, rep%full_determ_vecs(1, :), rep%determ_sizes, &
                               rep%determ_displs)

            rep%partial_determ_vecs = 0.0_dp
            ham_times_hf = 0.0_dp

            do i = 1, rep%determ_sizes(iProcIndex)
                do j = 1, rep%sparse_core_ham(i)%num_elements
                    ham_times_hf(i) = ham_times_hf(i) + &
                                      rep%sparse_core_ham(i)%elements(j) * rep%full_determ_vecs(1, rep%sparse_core_ham(i)%positions(j))
                end do
            end do

            write(6, '(a11,7X,a12,7X,a11)') "# Iteration", "Proj. Energy", "Var. Energy"
            call neci_flush(6)

            do while (iter <= NMCyc .or. NMCyc == -1)

                ! First, use the full Hamiltonian to get energy estimates

                wavefunction_old = wavefunction
                rep%partial_determ_vecs(1, :) = wavefunction

                call determ_projection()
                ham_times_wf = -rep%partial_determ_vecs(1, :) / tau + (DiagSft(1) * wavefunction)

                wavefunction = wavefunction + rep%partial_determ_vecs(1, :)

                energy_num = dot_product(ham_times_hf, wavefunction)
                ASSERT(.not. tOrthogonaliseReplicas)
                if (iProcIndex == iRefProc(1)) energy_denom = wavefunction(hf_index)

                var_e_num = dot_product(ham_times_wf, wavefunction_old)
                var_e_denom = dot_product(wavefunction_old, wavefunction_old)

                call MPISum(var_e_num, tot_var_e_num)
                call MPISum(var_e_denom, tot_var_e_denom)
                call MPISum(energy_num, tot_e_num)
                call MPISum(energy_denom, tot_e_denom)

                write(6, '(i9,7X,f13.10,7X,f13.10)') iter, tot_e_num / tot_e_denom, tot_var_e_num / tot_var_e_denom
                call neci_flush(6)

                ! Perform the actual projection used, with the approximate Hamiltonian

                wavefunction = wavefunction_old
                rep%partial_determ_vecs(1, :) = wavefunction

                call determ_proj_approx()

                wavefunction = wavefunction + rep%partial_determ_vecs(1, :)

                iter = iter + 1

            end do

            deallocate(wavefunction)
            deallocate(ham_times_hf)
            deallocate(ham_times_wf)
        end associate

    end subroutine perform_determ_proj_approx_ham

end module determ_proj
